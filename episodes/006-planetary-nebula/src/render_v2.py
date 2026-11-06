"""Single low-resolution narrative preview; no expensive production render."""
from pathlib import Path
import sys,argparse,time,json
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.fft import fftfreq,rfftfreq
from numba import njit,prange
from PIL import Image,ImageDraw,ImageFont
from .evolution import PC,RSUN,LSUN,smooth
from .physics_v2 import Nebula
from .imaging import Imaging
from .video_export import encoder,metadata
ROOT=Path(__file__).resolve().parents[1]

@njit(cache=True,inline='always')
def sample(field,pos,loglo,dlog):
    r=np.sqrt((pos*pos).sum())
    if r<=0:return np.zeros(6)
    fr=(np.log(r)-loglo)/dlog
    if fr<0 or fr>field.shape[2]-1:return np.zeros(6)
    mu=pos[2]/r;phi=np.arctan2(pos[1],pos[0])%(2*np.pi)
    fu=(mu+1)*.5*field.shape[0]-.5;fp=phi/(2*np.pi)*field.shape[1]-.5
    iu=int(np.floor(fu));ip=int(np.floor(fp));ir=min(int(fr),field.shape[2]-2)
    du=fu-iu;dp=fp-ip;dr=fr-ir
    iu0=max(0,min(field.shape[0]-1,iu));iu1=max(0,min(field.shape[0]-1,iu+1))
    ip0=ip%field.shape[1];ip1=(ip+1)%field.shape[1]
    out=np.zeros(6)
    for a in range(2):
        u=iu0 if a==0 else iu1;wa=1-du if a==0 else du
        for b in range(2):
            p=ip0 if b==0 else ip1;wb=1-dp if b==0 else dp
            for c in range(2):
                wc=1-dr if c==0 else dr
                for k in range(6):out[k]+=wa*wb*wc*field[u,p,ir+c,k]
    return out

@njit(cache=True,parallel=True)
def project(field,loglo,dlog,w,h,halfwidth,starR,star_image,basis,depth_steps,far_radius):
    image=np.zeros((h,w,3),np.float64);pixel=2*halfwidth/w;cy=.575*h
    # Orthographic external camera. World coordinates never rotate the gas.
    for iy in prange(h):
        yy=(cy-iy-.5)*pixel
        for ix in range(w):
            xx=(ix+.5-w*.5)*pixel;b2=xx*xx+yy*yy
            b=np.sqrt(b2)
            if b>=far_radius:continue
            zmax=np.sqrt(max(0,far_radius**2-b2))
            I=np.zeros(3);Iir=0.;zstart=-zmax
            if b<starR:
                zstart=np.sqrt(starR*starR-b2)
                for k in range(3):I[k]=star_image[iy,ix,k]
            # Log-spaced ray quadrature around small star, near-uniform at nebula scale.
            soft=max(b*.5,starR*.5,1e-12)
            umax=np.arcsinh(zmax/soft);umin=np.arcsinh(zstart/soft)
            du=(umax-umin)/depth_steps
            for iz in range(depth_steps):
                u=umin+(iz+.5)*du
                zz=soft*np.sinh(u);ds=soft*np.cosh(u)*du
                pos=xx*basis[0]+yy*basis[1]+zz*basis[2]
                v=sample(field,pos,loglo,dlog)
                tau=max(0,v[3]*ds);att=np.exp(-min(700,tau))
                fac=ds if tau<1e-6 else ds*(-np.expm1(-tau))/tau
                for k in range(3):I[k]=I[k]*att+v[k]*PC*fac
                tauir=max(0,v[5]*ds);attir=np.exp(-min(700,tauir))
                facir=ds if tauir<1e-6 else ds*(-np.expm1(-tauir))/tauir
                Iir=Iir*attir+v[4]*PC*facir
            ircol=np.array([1.,.18,.012]);ircol/=np.dot(ircol,np.array([.2126729,.7151522,.0721750]))
            for k in range(3):image[iy,ix,k]=I[k]+Iir*ircol[k]
    return image


from .evolution_v2 import Track
from .reuse_no004 import No004Adapter
from .reuse_no003 import CompactAdapter,CompactOptics,reuse_evidence,parameters
from copy import deepcopy

class Renderer:
    def __init__(self,w=360,h=640,depth=128,model=None,track=None,compact_preview_floor=2.5):
        self.w=w;self.h=h;self.depth=depth
        self.model=Nebula(nr=256) if model is None else model
        self.track=Track() if track is None else track;self.model.track=self.track
        self.solar=No004Adapter(spatial_average=2)
        self.compact=CompactAdapter(y_reference=self.track.Yref)
        self.compact_optics=CompactOptics((h,w),100000.,exposure=1.,scatter_fraction=.06,preview_floor_px=compact_preview_floor)
        cfg=deepcopy(self.solar.cfg['optics']);cfg.update(extra_halo_outside_disk_only=False,halo_display_gain=2.0,exposure=1.05)
        # Original No004 optical functions, with detector sampling floor for preview.
        cfg['airy_first_zero_550nm_at1080']=max(cfg['airy_first_zero_550nm_at1080'],2.5*1080/w)
        cfg['jitter_sigma_at1080']=max(cfg['jitter_sigma_at1080'],.4*1080/w)
        self.solar_optics=self.solar.optics((h,w),cfg)
        self.alpha=self.title();self.records=[]
        self.yy,self.xx=np.mgrid[:h,:w]
        self.stellar_metrics=[]
    def title(self):
        ss=3;s=self.w/1080;im=Image.new('L',(self.w*ss,self.h*ss));d=ImageDraw.Draw(im)
        for x,y,txt,size,font,anchor in [(84,250,'@sideratlas',34,'lmroman10-regular.otf','lt'),(996,253,'NO.006',25,'lmroman10-regular.otf','rt'),(540,320,'THE PLANETARY NEBULA',61,'lmroman12-bold.otf','mt')]:
            d.text((x*s*ss,y*s*ss),txt,font=ImageFont.truetype(str(ROOT/'fonts'/font),round(size*s*ss)),fill=255,anchor=anchor)
        return np.asarray(im.resize((self.w,self.h),Image.Resampling.LANCZOS),'f4')/255
    def basis(self,t):
        # Use No004's actual orientation convention so its real magnetic complex
        # is initially visible, rather than accidentally placed on the far side.
        az=-15+np.rad2deg(1.1)*t/32
        basis=self.solar.geometry.orientation(np.eye(3),67.,az)
        return basis,az
    def stellar_image(self,t,star,basis,pixel_radius,source_only=False):
        if pixel_radius<.65:
            # Do not render an enlarged core; sum the original spectral surface
            # over quadrature mu and deposit its flux only after attenuation.
            return np.zeros((self.h,self.w,3)),dict(resolved=False,source='No003 spectral flux + original PSF')
        evolved=float(smooth((star['R']*PC/RSUN-1.05)/1.6))
        compact=float(smooth((star['T']-9000)/14000))
        if compact>.999:
            image=self.compact.render_disk(self.w,self.h,pixel_radius,star,center=(self.w/2,self.h*.575))
            return image,dict(resolved=True,source='Original No003 Samples/Tile.project + colour_radiance',compact_weight=1.)
        x=(self.xx+.5-self.w/2)/pixel_radius;y=(self.h*.575-self.yy-.5)/pixel_radius
        q=x*x+y*y;mask=q<1;mu=np.sqrt(np.maximum(0,1-q[mask]))
        # Original No004's -20degree detector roll is retained in source coordinates.
        roll=np.deg2rad(-20);nx=np.cos(roll)*x[mask]+np.sin(roll)*y[mask];ny=-np.sin(roll)*x[mask]+np.cos(roll)*y[mask]
        normals=np.stack([nx,ny,mu],axis=-1)@basis
        observer=basis[2]
        solar=self.solar.source_solar(normals,mu,t,observer,grade=True,active=True) if evolved<1 else None
        if solar is not None:solar*=self.track.rgb(star['T'])/self.track.rgb(5772.)
        ratio,info=self.solar.morphology(normals,mu,t,observer,star['R']*PC/RSUN,star['T'],star['M'],contrast=2.0)
        # Original solar radiation fluctuations are a declared morphology surrogate.
        # Keep their evolving measured spatial structure; replace their thermodynamics.
        ratio=np.maximum(.05,ratio)
        # Visible-disk normalization is for the source proxy only; not gas luminosity.
        sphere_ratio,_=self.solar.morphology(self.model.patch_normals,np.ones(self.model.mu.size),t,observer,star['R']*PC/RSUN,star['T'],star['M'],contrast=2.0)
        ratio/=np.mean(sphere_ratio)
        base=self.compact.surface(mu,star)
        local_temp=star['T']*ratio**.25
        chrom=self.track.rgb(local_temp)/np.maximum(self.track.rgb(star['T']),1e-30)
        giant=base*chrom
        if solar is not None:values=(1-evolved)*solar+evolved*giant
        else:values=giant
        if compact>0:values=(1-compact)*values+compact*base
        image=np.zeros((self.h,self.w,3));image[mask]=values
        # Exact original No003 circle/pixel area for limb antialiasing, avoiding
        # a new generic hard disk even while using No004's stellar source.
        from .reuse_no003 import circle_rectangle_area
        rho=np.sqrt(q);edge=np.abs(rho-1)*pixel_radius<.8
        ey,ex=np.nonzero(edge)
        for iy,ix in zip(ey,ex):
            area=float(circle_rectangle_area(ix-self.w/2,ix+1-self.w/2,self.h*.575-iy-1,self.h*.575-iy,pixel_radius))
            if area<=0:continue
            xx=(ix+.5-self.w/2)/pixel_radius;yy=(self.h*.575-iy-.5)/pixel_radius
            rad=max(np.hypot(xx,yy),1e-10);xx*=min(1, (1-1e-8)/rad);yy*=min(1,(1-1e-8)/rad)
            mm=np.sqrt(max(0,1-xx*xx-yy*yy));nn=np.array([[np.cos(roll)*xx+np.sin(roll)*yy,-np.sin(roll)*xx+np.cos(roll)*yy,mm]])@basis
            # Nearest interior source pixel keeps actual source patterns at the edge;
            # exact overlap area supplies coverage rather than invented emission.
            iy0=int(np.clip(self.h*.575-yy*pixel_radius,0,self.h-1));ix0=int(np.clip(self.w/2+xx*pixel_radius,0,self.w-1))
            source=image[iy0,ix0] if image[iy0,ix0].sum()>0 else self.compact.surface(np.array([mm]),star)[0]
            image[iy,ix]=source*area
        return image,dict(info,resolved=True,source='Original No004 sample_patches/active_rectangle/grade + No003 spectral integration',evolved_weight=evolved,compact_weight=compact)
    def geometry(self,t,star,mode):
        if mode=='A' or mode=='B':return star['R']/.76
        if t<13:
            return star['R']/(.73+.05*float(smooth(t/6)))
        if t<18:
            # A continuous camera pullback; no physics radius or trajectory edits.
            start=self.track.state(float(self.track.clock(13)))['R']/.78
            return np.exp((1-float(smooth((t-13)/5)))*np.log(start)+float(smooth((t-13)/5))*np.log(.075))
        tr=float(self.track.clock(t));sh=self.model.shell_solution.sol(min(tr,45000))[:self.model.mu.size]
        return max(.075,float(np.quantile(sh,.95))/.72)
    def frame(self,t,fps=12,mode='full'):
        tr=float(self.track.clock(t));state=self.model.state(tr);star=state['star'];basis,az=self.basis(t)
        half=self.geometry(t,star,mode)
        wind_view=0.
        if mode=='full':
            wind_view=float(smooth((t-8)/1.0))*(1-float(smooth((t-16)/2.0)))*float(smooth((state['report']['dust_IR_tau_median']-.15)/1.3))
            weights=state['escaped_ir_power'].sum(axis=(0,1));cum=np.cumsum(weights)
            focus=float(np.interp(.8*cum[-1],cum,state['r'])) if cum[-1]>0 else 0.
            half=np.exp((1-wind_view)*np.log(half)+wind_view*np.log(max(half,focus/.65)))
        pix=2*half/self.w;rpx=star['R']/pix
        starim,sm=self.stellar_image(t,star,basis,rpx,source_only=mode in ['A','B'])
        info=dict(playback_s=t,physical_age_yr=tr+self.track.t0,post_AGB_reference_yr=tr,M=star['M'],R_Rsun=star['R']*PC/RSUN,Teff_K=star['T'],L_Lsun=star['L']/LSUN,stellar_radius_pixels=rpx,camera_halfwidth_pc=half,**sm)
        source_only=mode in ['A','B']
        compact=float(smooth((star['T']-9000)/14000))
        if source_only:
            hdr=starim.copy()
        else:
            fields=state['fields'].copy()
            # No 400x orange wash: independent IR display is bounded/compressed,
            # without touching physical opacity or IR-power records.
            ir_display_gain=4+44*wind_view
            fields[...,4]*=ir_display_gain
            unresolved=state['r']<pix
            irpoint=float(state['escaped_ir_power'][...,unresolved].sum())/(4*np.pi*(pix*PC)**2*self.track.Yref)
            fields[...,unresolved,4]=0
            hdr=project(fields,np.log(state['r'][0]),np.log(state['r'][1]/state['r'][0]),self.w,self.h,half,star['R'] if rpx>=.65 else 0.,starim,basis,self.depth,float(state['r'][-1]))
            ircolor=np.array([1,.22,.025]);ircolor/=ircolor@np.array([.2126729,.7151522,.0721750])
            hdr[round(self.h*.575),self.w//2]+=ir_display_gain*irpoint*ircolor
            if rpx<.65:
                # Original No003 spectral/limb integrator, full apparent disk flux.
                mu,wgt=np.polynomial.legendre.leggauss(24);mu=(mu+1)/2;wgt=wgt/2
                radiance=self.compact.surface(mu,star)
                flux=2*np.pi*rpx*rpx*np.sum(radiance*(mu*wgt)[:,None],axis=0)
                observer=basis[2];u=int(np.clip((observer[2]+1)*.5*self.model.nmu,0,self.model.nmu-1));ph=np.arctan2(observer[1],observer[0])%(2*np.pi);p=int(ph/(2*np.pi)*self.model.nphi)%self.model.nphi
                tau=float(np.sum(state['fields'][u,p,:,3]*np.diff(state['edges'])))
                hdr[round(self.h*.575),self.w//2]+=flux*np.exp(-min(tau,700))
                info['central_source']='Flux-integrated No003 atmosphere method, MIST R/T/M, physical dust extinction'
        # Actual solar C IV source only while the source remains close to solar.
        solar_activity=1-float(smooth((star['R']*PC/RSUN-1.05)/.95))
        if solar_activity>1e-4 and mode!='B' and tr < -1e6 and star['T']<7000:
            camera=dict(motion='uniform-orbit',inclination_deg=67.,azimuth_start_deg=-15.,revolutions=1.1/(2*np.pi))
            uv=self.solar.atmosphere_image(t,self.w,self.h,[.5,.575],rpx/self.w,camera)
            rho=np.hypot(self.xx+.5-self.w/2,self.yy+.5-self.h*.575)/max(rpx,1e-12)
            visibility=np.minimum(20.,np.exp(np.minimum(30,np.maximum(0,rho-1)/.04)))
            channel=.8*np.arcsinh(uv*visibility/1e4)/np.arcsinh(10.)
            col=np.array([1.,.23,.035]);col/=col@np.array([.2126729,.7151522,.0721750])
            hdr+=solar_activity*channel[...,None]*col
            info['CIV']='Actual No004 Atmosphere.image, faded to zero before giant; original independent UV exposure'
        if source_only:
            # Normalize to disk-center spectral radiance to inspect the source,
            # explicitly a camera/exposure test, not a bolometric luminosity plot.
            center=self.compact.surface(np.array([1.]),star)[0]@np.array([.2126729,.7151522,.0721750])
            exposure=.95/max(center,1e-20)
        else:
            fade=float(smooth((t-15)/4));ne=5e12*(max(100.,tr)/16000.)**3
            center=self.compact.surface(np.array([1.]),star)[0]@np.array([.2126729,.7151522,.0721750])
            photos_exposure=.95/max(center,1e-20)
            exposure=float(np.exp((1-fade)*np.log(photos_exposure)+fade*np.log(ne)))
        # Execute actual No004 optics at the start, actual No003 optics for the
        # evolved/compact source; the transition is a declared imaging blend.
        optic_weight=float(smooth((t-3)/13)) if mode!='B' else 1.
        if optic_weight<=0:
            self.solar.calls['src/optics.py:Optics.process']=self.solar.calls.get('src/optics.py:Optics.process',0)+1
            im=self.solar_optics.process(hdr*exposure)
        elif optic_weight>=1:im=self.compact_optics.process(hdr*exposure)
        else:
            self.solar.calls['src/optics.py:Optics.process']=self.solar.calls.get('src/optics.py:Optics.process',0)+1
            im=(1-optic_weight)*self.solar_optics.process(hdr*exposure)+optic_weight*self.compact_optics.process(hdr*exposure)
        im=im*(1-self.alpha[...,None])+self.alpha[...,None]*.9
        srgb=np.where(im<=.0031308,12.92*im,1.055*np.maximum(0,im)**(1/2.4)-.055)
        info['wind_view_weight']=wind_view;info['IR_display_gain']=4+44*wind_view
        info['mode']=mode;info['exposure']=exposure;info['accounting']=state['report'];self.records.append(info)
        return np.uint16(np.clip(srgb,0,1)*65535+.5)
    def evidence(self):
        return dict(no004=self.solar.provenance(),no003=reuse_evidence(),scope='Actual source invocation inventory; not inherited validation of changed physics')
