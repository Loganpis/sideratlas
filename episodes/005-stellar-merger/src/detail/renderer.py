"""Native EWA cooling structures and energy-conserving luminous gas photosphere."""
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from src.motion_renderer import MotionRenderer
from src.detail.fields import coarse_sigma,fine_covariances,splat_emission,optical_depth,warm_chromaticity,gas_display_grade,Y
from src.detail.reprocessing import render_photosphere,RM,MM,CV,SIGMA

class DetailRenderer(MotionRenderer):
    def camera_geometry(self,playback,t=None):
        from src.dynamics.prescribed import smooth
        t=float(self.clock(playback)) if t is None else t;L=float(self.camera_width(min(t,390)));blend=float(smooth((t-240)/60));L=(1-blend)*L+blend*self.late_width;inc=.82-.12*np.clip((min(t,300)-215)/175,0,1);az=-np.pi/2+.35*min(playback/32,.75)
        if getattr(self,'fixed_diagnostic_camera',None) is not None:L,inc,az=self.fixed_diagnostic_camera
        eye=np.array([np.sin(inc)*np.cos(az),np.sin(inc)*np.sin(az),np.cos(inc)]);right=np.array([-np.sin(az),np.cos(az),0]);up=np.cross(eye,right)
        return L,inc,az,np.column_stack([right,up,eye])
    def __init__(self,source,width=1080):
        super().__init__(source,width,128)
        folder=Path(__file__).resolve().parents[2]/'cache/radiation';a=np.load(folder/'radiation-layer.npz');self.radiation={k:a[k] for k in a.files};self.angular_tree=cKDTree(self.radiation['directions']);self.layers={}
        # The field must include the actual emitting photosphere; keep late
        # camera scale fixed, not a zoom that cancels expansion.
        d=self.radiation['directions'];R=self.radiation['radius'][-1];inc=.82-.12*(300-215)/175;az=-np.pi/2+.35*.75;eye=np.array([np.sin(inc)*np.cos(az),np.sin(inc)*np.sin(az),np.cos(inc)]);right=np.array([-np.sin(az),np.cos(az),0]);up=np.cross(eye,right);x=d*R[:,None];needed=np.maximum(abs(x@right),np.maximum((x@up)/(1.28*16/9),-(x@up)/(.72*16/9)));self.late_width=max(float(self.camera_width(390)),1.15*float(np.quantile(needed,.98)))
    def stellar_optical_depth(self,t,q,ids,sigma,p,depth,camera,pixel):
        return optical_depth(p,sigma,np.full(len(ids),self.motion.mass),depth,self.W,self.H,pixel,.64*self.H,.01*MM/RM**2) if len(ids) else np.zeros(depth.shape)
    def gas_chromaticity(self,T,ids,t,q):
        return warm_chromaticity(T)
    def gas_display_mapping(self,rgb):
        return gas_display_grade(rgb)
    def radiation_at(self,t):
        a=self.radiation;i=max(0,min(len(a['t'])-2,np.searchsorted(a['t'],t)-1));f=np.clip((t-a['t'][i])/(a['t'][i+1]-a['t'][i]),0,1)
        return (1-f)*a['radius'][i]+f*a['radius'][i+1],(1-f)*a['luminosity'][i]+f*a['luminosity'][i+1]
    def emitting_layer(self,t,q,right,up,eye,L,power):
        a=self.radiation;i=max(0,min(len(a['t'])-2,np.searchsorted(a['t'],t)-1));f=np.clip((t-a['t'][i])/(a['t'][i+1]-a['t'][i]),0,1);tau=(1-f)*a['outward_tau'][i]+f*a['outward_tau'][i+1] if 'outward_tau' in a else None;rr=(a['radial_edges'][:-1]+a['radial_edges'][1:])/2
        levels,weights=np.polynomial.laguerre.laggauss(3);r=np.empty((len(power),3))
        if 'layer_radii' in a:r=(1-f)*a['layer_radii'][i]+f*a['layer_radii'][i+1]
        else:
            for j in range(len(power)):
                for k,level in enumerate(levels):r[j,k]=max(q['radius'],np.interp(level,tau[j,::-1],rr[::-1]))
        dirs=a['directions'];omega=a['omega'];mu=np.maximum(0,dirs@eye);visible=4*mu*(.4+.6*mu)/.8;camera=np.column_stack([right,up,eye]);pp=(dirs[:,None,:]*r[:,:,None]+q['donor'])@camera;Teff=(power[:,None]/np.maximum(omega[:,None]*(r*RM)**2*SIGMA,1e-100))**.25
        raw=np.stack([np.interp(np.log(np.maximum(Teff,300)),np.log(self.gas_lut[0]),self.gas_lut[1][:,j]) for j in range(3)],-1);color=warm_chromaticity(Teff)*(raw@Y)[...,None]/np.maximum(SIGMA*np.maximum(Teff,1)**4,1e-100)[...,None]
        tang=r*np.sqrt(omega[:,None]/(2*np.pi))*1.1;radial=np.minimum(.15*r,np.maximum(.015,(r[:,0]-r[:,2])[:,None]/4));cov=tang[:,:,None,None]**2*np.eye(3)+(radial**2-tang**2)[:,:,None,None]*dirs[:,None,:,None]*dirs[:,None,None,:];cov=np.einsum('ij,nmjk,kl->nmil',camera.T,cov,camera)
        pp=pp.reshape(-1,3);cov=cov.reshape(-1,3,3);base_power=power[:,None]*weights[None,:]
        modulation=self.emissivity_modulation((dirs[:,None,:]*r[:,:,None]+q['donor']).reshape(-1,3),t,q) if hasattr(self,'emissivity_modulation') else None
        if modulation is not None:
            modulation=modulation.reshape(base_power.shape);weighted=base_power*modulation;weighted*=base_power.sum()/max(weighted.sum(),1e-100);base_power=weighted
        P=(base_power*visible[:,None]).ravel();color=color.reshape(-1,3);pixel=2*L/self.W
        img,used=splat_emission(pp,cov,P,color,np.full((self.H,self.W),-np.inf),self.W,self.H,pixel,.64*self.H,False);img/=4*(pixel*RM)**2
        # Fixed photospheric display shoulder, not a new luminosity source.
        y=img@Y;scale=.02*np.log1p(y/.02)/np.maximum(y,1e-100);img*=scale[...,None]
        self.emitting_metrics=dict(allocated_physical_luminosity_W=float(power.sum()),view_isotropic_equivalent_source_W=float(P.sum()),view_power_in_native_footprints_W=float(used),last_interaction_weight_sum=float(weights.sum()))
        return img,(y>1e-5).astype('f4'),Teff.ravel()
    def curved_cooling(self,t,q,ids,p,cov,source,colors,camera,depth,pixel):
        # Five-point material-kernel quadrature. Curvature is the measured
        # trajectory acceleration normal to velocity, not a drawn spline.
        old=self.motion.at(max(0,t-.15));new=self.motion.at(min(390,t+.15));acc=(new['v'][ids]-old['v'][ids])/.3;v=q['v'][ids];speed=np.linalg.norm(v,axis=1);axis=v/np.maximum(speed[:,None],1e-12);normal=acc-axis*np.sum(acc*axis,axis=1)[:,None];curvature=normal/np.maximum(speed[:,None]**2,.01)
        length=np.sqrt(np.einsum('ni,nij,nj->n',axis,cov,axis));bend=np.linalg.norm(curvature,axis=1)*length;curvature*=np.minimum(1,.35/np.maximum(bend,1e-12))[:,None]
        nodes,w=np.polynomial.hermite.hermgauss(5);nodes*=np.sqrt(2);w/=np.sqrt(np.pi);sub=cov-.85*length[:,None,None]**2*axis[:,:,None]*axis[:,None,:];val,vec=np.linalg.eigh(sub);val=np.maximum(val,.025**2);sub=np.einsum('nij,nkj->nik',vec*np.sqrt(val)[:,None,:],vec*np.sqrt(val)[:,None,:]);sub=np.einsum('ij,njk,kl->nil',camera.T,sub,camera)
        result=np.zeros((self.H,self.W,3),np.float32);used=0.
        for node,weight in zip(nodes,w):
            shift=np.sqrt(.85)*length[:,None]*node*axis+.5*.85*length[:,None]**2*(node*node-1)*curvature
            image,pow_in=splat_emission(p+shift@camera,sub,.85*source*weight,colors,depth,self.W,self.H,pixel,.64*self.H,False);result+=image;used+=pow_in
        return result,used
    def render_linear(self,playback):
        from src.dynamics.prescribed import smooth
        t=float(self.clock(playback));q=self.motion.at(t);L,inc,az,camera=self.camera_geometry(playback,t);right,up,eye=camera.T;pixel=2*L/self.W
        lum=float(np.interp(t,self.times,self.data['gas_luminosity_code']+self.data['remnant_luminosity_code']))*self.code_lum;star,depth,mask,Lgas=self.stars_native(t,q['phase'],right,up,eye,L,np.zeros(3),lum)
        ids=np.flatnonzero(q['live']);n=len(ids);sigma=coarse_sigma(q['x'][ids],q['radius'],q['donor']) if n else np.zeros(0);p=q['x'][ids]@camera;tau=self.stellar_optical_depth(t,q,ids,sigma,p,depth,camera,pixel);direct=star*np.exp(-np.minimum(tau,100))[...,None]
        rph,power=self.radiation_at(t);photo,photo_mask,photo_T=self.emitting_layer(t,q,right,up,eye,L,power)
        photo_depth=np.full(depth.shape,-np.inf)
        fine=np.zeros_like(star);coarse=np.zeros_like(star);fit=dict(parcels=0,fitted=0);sourcepower=0.;in_power=0.
        if n:
            cov,fit=fine_covariances(self.motion,t,q,ids,sigma);cc=np.einsum('ij,njk,kl->nil',camera.T,cov,camera);broad=np.eye(3)[None]*sigma[:,None,None]**2
            source=q['T'][ids]*self.motion.mass*MM*CV/(self.config['gas_cooling_time_code']*self.report['physical_time_unit_seconds']);sourcepower=float(source.sum());source*=min(1,Lgas/max(sourcepower,1e-100));colors=self.gas_chromaticity(q['T'][ids],ids,t,q)
            # Diagnostic bolometric cooling overlay: NOT a second physical
            # escaped luminosity. Fixed coarse/fine split preserves its scalar
            # source integral; physical photons are in the reservoir ledger.
            fine,fi=self.curved_cooling(t,q,ids,p,cov,source,colors,camera,depth,pixel);coarse_fraction=1-float(getattr(self,'cooling_fine_fraction',.85))
            if coarse_fraction>0:coarse,co=splat_emission(p,broad,coarse_fraction*source,colors,depth,self.W,self.H,pixel,.64*self.H,False)
            else:coarse=np.zeros_like(fine);co=0.
            in_power=fi+co
            factor=.01/(4*(pixel*RM)**2*SIGMA*self.report['donor']['Teff']**4);total=(fine+coarse)*factor;brightness=total@Y;scale=.06*np.arcsinh(brightness/.06)/np.maximum(brightness,1e-100);fine*=factor*scale[...,None];coarse*=factor*scale[...,None]
        gas_rgb=self.gas_display_mapping(photo+fine+coarse);composite=direct+gas_rgb
        visible=np.divide(direct.sum(axis=-1),star.sum(axis=-1),out=np.zeros_like(mask),where=star.sum(axis=-1)>1e-20);final_mask=np.maximum(mask*np.clip(visible,0,1),photo_mask)
        self.layers=dict(stellar_emission=star,optical_depth=tau,attenuated_stars=direct,reprocessed_photosphere=photo,fine_cooling=fine,coarse_cooling=coarse,composite_linear=composite,photosphere_depth=photo_depth)
        gas=n*self.motion.mass;m2=q['m2'];M1=1.68-gas-m2
        self.records.append(dict(playback_s=float(playback),t_code=t,width=self.W,gas_mass_msun=gas,mass_identity_error=abs(M1+m2+gas-1.68),source_cooling_power_W=sourcepower,assigned_gas_power_W=Lgas,rendered_photosphere_luminosity_W=float(power.sum()),photosphere_radii_Rsun=np.quantile(rph,[.1,.5,.9]).tolist(),photosphere_Teff_K=np.quantile(photo_T[photo_T>0],[.1,.5,.9]).tolist() if np.any(photo_T>0) else [0,0,0],stellar_flux_survival=float(direct.sum()/max(star.sum(),1e-100)),camera_halfwidth_Rsun=L,late_width_Rsun=self.late_width,camera_azimuth=az,camera_inclination=inc,source_overlay_cooling_fraction=.85,fine_footprints=fit,frame_overlay_power_in_W=float(in_power),emitting_layer=self.emitting_metrics,warm_grade='Fixed monotone physical-T to display color-temperature 1700..5200 K, luminance preserved; display colors not measured temperatures',overlay_is_not_extra_escaped_radiation=True,cooling_projection_is_unoccluded_energy_field=True,source_cache_unchanged=True))
        self.records[-1]['source_overlay_cooling_fraction']=float(getattr(self,'cooling_fine_fraction',.85))
        return composite,final_mask
