"""Research-informed proof of concept using the actual No004/003 source adapters."""
import numpy as np
from numba import njit,prange
from .render_v2 import Renderer as ReusedRenderer,sample
from .physics_v3 import Nebula
from .evolution_v3 import Track
from .evolution import PC,RSUN,smooth

@njit(cache=True,parallel=True)
def project_material(field,loglo,dlog,w,h,halfwidth,starR,star_image,basis,rnodes,ir_gain):
    image=np.zeros((h,w,3));pixel=2*halfwidth/w;cy=.575*h
    ircol=np.array([1.,.42,.16]);ircol/=np.dot(ircol,np.array([.2126729,.7151522,.0721750]))
    for iy in prange(h):
        y=(cy-iy-.5)*pixel
        for ix in range(w):
            x=(ix+.5-w*.5)*pixel;b2=x*x+y*y
            I=np.zeros(3);Iir=0.;zstar=-1e100
            if b2<starR*starR:
                zstar=np.sqrt(starR*starR-b2)
                for k in range(3):I[k]=star_image[iy,ix,k]
            # Actual radial-cell/opacity/emission importance nodes. Both sides of
            # the star are integrated; source occultation is retained.
            for side in range(2):
                for j in range(len(rnodes)-1):
                    k=len(rnodes)-2-j if side==0 else j
                    lo=rnodes[k];hi=rnodes[k+1]
                    if hi*hi<=b2:continue
                    zl=np.sqrt(max(0,lo*lo-b2));zh=np.sqrt(max(0,hi*hi-b2))
                    z0=-zh if side==0 else zl;z1=-zl if side==0 else zh
                    z0=max(z0,zstar)
                    if z1<=z0:continue
                    # Two-point quadrature resolves angular variations as the
                    # ray crosses a radial material layer; no density contours.
                    for q in range(2):
                        zz=z0+(q+.5)*.5*(z1-z0);ds=.5*(z1-z0)
                        pos=x*basis[0]+y*basis[1]+zz*basis[2]
                        f=sample(field,pos,loglo,dlog)
                        tau=max(0,f[3]*ds);att=np.exp(-min(700,tau));fac=ds if tau<1e-6 else ds*(-np.expm1(-tau))/tau
                        for c in range(3):I[c]=I[c]*att+f[c]*PC*fac
                        ti=max(0,f[5]*ds);ai=np.exp(-min(700,ti));fi=ds if ti<1e-6 else ds*(-np.expm1(-ti))/ti
                        Iir=Iir*ai+f[4]*PC*fi
            for c in range(3):image[iy,ix,c]=I[c]+ir_gain*Iir*ircol[c]
    return image

class Renderer(ReusedRenderer):
    def __init__(self,w=288,h=512,depth=112,dynamics='radial-v3t-12x24-144.npz'):
        track=Track();model=Nebula(nr=384,dynamics=dynamics)
        super().__init__(w,h,depth,model=model,track=track,compact_preview_floor=1.2)
        self.compact_optics.cfg['scatter_fraction']=.025
        self.last_layers={}
    def radial_nodes(self,state):
        r=state['r'];e=state['edges'];fields=state['fields'];dv=self.model.omega*np.diff(e**3)/3
        power=np.maximum(fields[...,:3].sum(axis=-1),0).sum(axis=(0,1))*dv
        op=fields[...,3].sum(axis=(0,1))*np.diff(e)
        # Importance changes integration sampling only; not the physical field.
        weights=.75*power/max(power.sum(),1e-99)+.25*op/max(op.sum(),1e-99)
        cdf=np.r_[0,np.cumsum(weights)];cdf/=max(cdf[-1],1e-99)
        imp=np.interp(np.linspace(0,1,self.depth),cdf,e)
        return np.unique(np.r_[0.,imp,np.geomspace(e[0],e[-1],32),e[-1]])
    def geometry_v3(self,t,star,state):
        giant_scale=float(smooth(np.log10(max(star['R']*PC/RSUN,1.))/2.2))
        near=star['R']/(.50+.27*giant_scale)
        power=state['escaped_ir_power'].sum(axis=(0,1));cum=np.cumsum(power)
        ir95=float(np.interp(.95*cum[-1],cum,state['r'])) if cum[-1]>0 else star['R']
        env=max(near,ir95/.65)
        thick=float(smooth((state['report']['dust_IR_tau_median']-.3)/.5))
        w=max(float(smooth((t-8)/6)),thick)
        half=np.exp((1-w)*np.log(near)+w*np.log(env))
        if t>=13:
            e=state['edges']
            if 'line_fields' in state:
                p=state['line_fields'].sum(axis=-1).sum(axis=(0,1))*np.diff(e**3)
                c=np.cumsum(np.maximum(0,p));rr=float(np.interp(.96*c[-1],c,state['r'])) if c[-1]>0 else 0.
                target=max(.003,rr/.78)
            else:target=.003
            w=float(smooth((t-13)/5))
            half=np.exp((1-w)*np.log(half)+w*np.log(target))
        return float(half)
    def frame(self,t,fps=10,mode='full'):
        tr=float(self.track.clock(t));state=self.model.state(tr);star=state['star'].copy();basis,az=self.basis(t)
        # Small representative AGB radial pulsation; fixed L, no extra radiative
        # energy. Its display clock is separate from the secular mass-loss clock.
        agb=float(smooth((star['R']*PC/RSUN-40)/80))*float(smooth((4500-star['T'])/1000))
        pulse=1+.025*agb*np.sin(2*np.pi*(t-8)/1.8)
        star['R']*=pulse;star['T']/=np.sqrt(pulse)
        half=self.geometry_v3(t,state['star'],state);pix=2*half/self.w;rpx=star['R']/pix
        sim,sm=self.stellar_image(t,star,basis,rpx)
        center=self.compact.surface(np.array([1.]),star)[0]@np.array([.2126729,.7151522,.0721750])
        source_exposure=.95/max(center,1e-20)
        fade=float(smooth((t-14)/5));nebula_exposure=3e12*(max(100.,tr)/16000)**3
        base_exposure=float(np.exp((1-fade)*np.log(source_exposure)+fade*np.log(nebula_exposure)))
        flux=np.zeros(3);point_reference=0.
        if rpx<.65:
            mu,wg=np.polynomial.legendre.leggauss(24);mu=(mu+1)/2;wg/=2
            flux=2*np.pi*rpx*rpx*np.sum(self.compact.surface(mu,star)*(mu*wg)[:,None],axis=0)
            point_gate=float(smooth((t-15)/2))*float(smooth((star['T']-6500)/3500))
            # Exposure uses the intrinsic flux reference; physical extinction is
            # still applied afterward. It cannot reveal an opaque hidden source.
            point_reference=point_gate*2/max(flux@np.array([.2126729,.7151522,.0721750]),1e-30)
        exposure=max(base_exposure,point_reference)
        ir_window=float(smooth((t-9)/1))*(1-float(smooth((t-16)/2)))
        ir_detector_exposure=8000.
        irgain=(1-ir_window)*.12+ir_window*ir_detector_exposure/exposure
        nodes=self.radial_nodes(state)
        fields=state['fields'].copy()
        line_sensitivity=100. # explicit long narrowband / short stellar-continuum HDR composite
        if 'line_fields' in state:fields[...,:3]+=(line_sensitivity-1)*state['line_fields']
        unresolved=state['r']<pix
        irpoint=float(state['escaped_ir_power'][...,unresolved].sum())/(4*np.pi*(pix*PC)**2*self.track.Yref)
        fields[...,unresolved,4]=0
        hdr=project_material(fields,np.log(state['r'][0]),np.log(state['r'][1]/state['r'][0]),self.w,self.h,half,star['R'] if rpx>=.65 else 0.,sim,basis,nodes,irgain)
        ircolor=np.array([1.,.42,.16]);ircolor/=ircolor@np.array([.2126729,.7151522,.0721750])
        hdr[round(self.h*.575),self.w//2]+=irgain*irpoint*ircolor
        star_tau=None;transmitted_point=0.
        if rpx<.65:
            obs=basis[2];u=int(np.clip((obs[2]+1)*.5*self.model.nmu,0,self.model.nmu-1));ph=np.arctan2(obs[1],obs[0])%(2*np.pi);p=int(ph/(2*np.pi)*self.model.nphi)%self.model.nphi
            tau=float(np.sum(state['fields'][u,p,:,3]*np.diff(state['edges'])))
            hdr[round(self.h*.575),self.w//2]+=flux*np.exp(-min(700,tau))
            star_tau=tau;transmitted_point=float(flux@np.array([.2126729,.7151522,.0721750]))*np.exp(-min(700,tau))
        # Original No004 ultraviolet projector is retained only at solar stage.
        solar=1-float(smooth((star['R']*PC/RSUN-1.05)/.95))
        if solar>1e-4 and tr < -1e6 and star['T']<7000:
            cam=dict(motion='uniform-orbit',inclination_deg=67.,azimuth_start_deg=-15.,revolutions=1.1/(2*np.pi))
            uv=self.solar.atmosphere_image(t,self.w,self.h,[.5,.575],rpx/self.w,cam)
            rho=np.hypot(self.xx+.5-self.w/2,self.yy+.5-self.h*.575)/max(rpx,1e-12)
            channel=.8*np.arcsinh(uv*np.minimum(20.,np.exp(np.minimum(30,np.maximum(0,rho-1)/.04)))/1e4)/np.arcsinh(10.)
            col=np.array([1.,.23,.035]);col/=col@np.array([.2126729,.7151522,.0721750]);hdr+=solar*channel[...,None]*col
        ow=float(smooth((t-3)/13))
        if ow<1:
            self.solar.calls['src/optics.py:Optics.process']=self.solar.calls.get('src/optics.py:Optics.process',0)+1
            a=self.solar_optics.process(hdr*exposure)
        if ow>0:
            linear=self.compact_optics.linear_process(hdr*exposure)
            # Fixed hue-preserving maximum-channel Reinhard display map. It
            # retains the calculated line-color ordering instead of whitening
            # every bright nebular region independently in R,G and B.
            b=linear/(1+np.max(linear,axis=-1)[...,None])
        im=a if ow==0 else (b if ow==1 else (1-ow)*a+ow*b)
        self.last_layers=dict(hdr=hdr,star_unattenuated=sim,state=state,halfwidth=half,exposure=exposure)
        im=im*(1-self.alpha[...,None])+self.alpha[...,None]*.9
        srgb=np.where(im<=.0031308,12.92*im,1.055*np.maximum(0,im)**(1/2.4)-.055)
        info=dict(playback_s=float(t),physical_age_yr=float(tr+self.track.t0),post_AGB_reference_yr=tr,M=star['M'],R_Rsun=star['R']*PC/RSUN,Teff_K=star['T'],L_Lsun=state['report']['L_Lsun'],stellar_radius_pixels=rpx,camera_halfwidth_pc=half,AGB_radius_pulsation_factor=float(pulse),representative_pulsation_display_period_s=1.8,representative_pulsation_period_yr=1.388*(state['star']['R']*PC/RSUN/372)**1.5/np.sqrt(star['M']),IR_display_gain=irgain,line_to_continuum_display_sensitivity=line_sensitivity,exposure=exposure,ray_radial_nodes=len(nodes),accounting=state['report'],**sm)
        info.update(IR_separate_exposure_window=ir_window,IR_detector_exposure=ir_detector_exposure,base_visible_exposure=base_exposure,intrinsic_point_reference_exposure=point_reference)
        info['central_components']=dict(resolved_stellar_surface=bool(rpx>=.65),stellar_radial_tau_V=star_tau,transmitted_stellar_point_Y=transmitted_point,unresolved_dust_IR_point_Y=float(irgain*irpoint),scope='detector-normalized contributions beforePSF/tone; extended gas not included; not bolometric luminosities')
        self.records.append(info)
        return np.uint16(np.clip(srgb,0,1)*65535+.5)
    def evidence(self):
        report=super().evidence()
        report.update(interface=self.model.interface.report,photon_history=str(self.model.photo_path),
                      detector_sampling=self.compact_optics.sampling_check,
                      display='100x line-to-continuum sensitivity; separate8000 IR exposure in9–18s window; intrinsic-point exposure reference retains extinction; line-power camera framing; no intrinsic luminosity increase')
        return report
