"""100% cooling-display source on persistent synchronized compact emitters."""
import numpy as np
from numba import njit
from scipy.ndimage import gaussian_filter
from scipy.optimize import brentq
from src.detail.synchronized_microflow import SynchronizedRenderer,SynchronizedFlow
from src.detail.microflow import MicroFlow,balanced_nodes
from src.detail.fields import Y,warm_chromaticity
from src.detail.reprocessing import RM,SIGMA

@njit(cache=True)
def b2(x):
    x=abs(x)
    if x<.5:return .75-x*x
    if x<1.5:return .5*(1.5-x)**2
    return 0.

@njit(cache=True)
def compact_emission(p,power,colors,W,H,pixel,cy):
    """Positive quadratic B-spline: compact 3-pixel support, unit integral."""
    image=np.zeros((H,W,3),np.float32);used=0.
    for j in range(len(p)):
        if power[j]<=0:continue
        x=W/2+p[j,0]/pixel-.5;y=cy-p[j,1]/pixel-.5;ix=int(np.floor(x+.5));iy=int(np.floor(y+.5))
        if x< -1.5 or x>W+.5 or y< -1.5 or y>H+.5:continue
        for b in range(iy-1,iy+2):
            if b<0 or b>=H:continue
            wy=b2(b-y)
            for a in range(ix-1,ix+2):
                if a<0 or a>=W:continue
                w=wy*b2(a-x);value=power[j]*w;used+=value
                for c in range(3):image[b,a,c]+=value*colors[j,c]
    return image,used

@njit(cache=True)
def compact_path_emission(p0,p1,power,colors,W,H,pixel,cy,step):
    """Bounded model-path integral, at most step pixels per spatial sample."""
    image=np.zeros((H,W,3),np.float32);used=0.
    for j in range(len(p0)):
        if power[j]<=0:continue
        x0=W/2+p0[j,0]/pixel-.5;y0=cy-p0[j,1]/pixel-.5;x1=W/2+p1[j,0]/pixel-.5;y1=cy-p1[j,1]/pixel-.5
        if max(x0,x1)<-1.5 or min(x0,x1)>W+.5 or max(y0,y1)<-1.5 or min(y0,y1)>H+.5:continue
        count=max(1,int(np.ceil(np.hypot(x1-x0,y1-y0)/step)));value0=power[j]/count
        for k in range(count):
            f=(k+.5)/count;x=x0+f*(x1-x0);y=y0+f*(y1-y0);ix=int(np.floor(x+.5));iy=int(np.floor(y+.5))
            for b in range(iy-1,iy+2):
                if b<0 or b>=H:continue
                wy=b2(b-y)
                for a in range(ix-1,ix+2):
                    if a<0 or a>=W:continue
                    v=value0*wy*b2(a-x);used+=v
                    for c in range(3):image[b,a,c]+=v*colors[j,c]
    return image,used

class ParticleFlow(SynchronizedFlow):
    def cloud(self,t,ids):
        # Exact same material positions; only raster footprints change.
        return MicroFlow.cloud(self,t,ids,return_footprints=False)

class ParticleRenderer(SynchronizedRenderer):
    def __init__(self,source,width=1080,settings_path=None):
        super().__init__(source,width,settings_path=settings_path);self.micro=ParticleFlow(self.control_micro,self.motion);self.cooling_fine_fraction=1.;self.particle_bloom_fraction=.005;self.photosphere_samples_per_cell_node=256;self.path_segments=2;self.path_pixel_step=.35
        if self.micro.settings.get('force_quintic_accelerations'):
            from src.detail.microflow import ROOT
            self.motion.enable_force_quintic(ROOT/self.micro.settings['force_quintic_accelerations'])
        dirs=self.radiation['directions'];east=np.cross(dirs,np.array([0.,0.,1.]));bad=np.linalg.norm(east,axis=1)<1e-6;east[bad]=np.cross(dirs[bad],np.array([1.,0.,0.]));east/=np.linalg.norm(east,axis=1,keepdims=True);self.ph_east=east;self.ph_north=np.cross(dirs,east);z=balanced_nodes(self.photosphere_samples_per_cell_node)[:,:3];labels=[]
        for j in range(len(dirs)):
            Q,_=np.linalg.qr(np.random.default_rng(5005999+j*53).normal(size=(3,3)));labels.append(z@Q)
        self.ph_labels=np.asarray(labels)
    def photo_radii(self,t):
        a=self.radiation;i=max(0,min(len(a['t'])-2,np.searchsorted(a['t'],t)-1));f=np.clip((t-a['t'][i])/(a['t'][i+1]-a['t'][i]),0,1)
        return (1-f)*a['layer_radii'][i]+f*a['layer_radii'][i+1]
    def photo_world(self,t,cell_ids=None):
        ii=slice(None) if cell_ids is None else cell_ids;r=self.photo_radii(t)[ii];dirs=self.radiation['directions'][ii];omega=self.radiation['omega'][ii];tang=r*np.sqrt(omega[:,None]/(2*np.pi))*1.1;radial=np.minimum(.15*r,np.maximum(.015,(r[:,0]-r[:,2])[:,None]/4));z=self.ph_labels[ii];center=self.motion.at(t)['donor']
        return center+r[:,:,None,None]*dirs[:,None,None,:]+tang[:,:,None,None]*(z[:,None,:,0,None]*self.ph_east[ii,None,None,:]+z[:,None,:,1,None]*self.ph_north[ii,None,None,:])+radial[:,:,None,None]*z[:,None,:,2,None]*dirs[:,None,None,:]
    def emitting_layer(self,t,q,right,up,eye,L,power):
        # Continuous field retained for physical geometry and the halo mask;
        # none of its Gaussian emission enters the displayed gas image.
        reference,mask,_=super().emitting_layer(t,q,right,up,eye,L,power);self.reprocessed_reference=reference;r=self.photo_radii(t);dirs=self.radiation['directions'];omega=self.radiation['omega'];_,weights=np.polynomial.laguerre.laggauss(3);mu=np.maximum(0,dirs@eye);visible=4*mu*(.4+.6*mu)/.8;Teff=(power[:,None]/np.maximum(omega[:,None]*(r*RM)**2*SIGMA,1e-100))**.25
        if not np.any(power>0):
            self.emitting_metrics=dict(allocated_physical_luminosity_W=0.,view_isotropic_equivalent_source_W=0.,view_power_in_native_footprints_W=0.,last_interaction_weight_sum=float(weights.sum()),reprocessing_source_samples=len(dirs)*3*self.photosphere_samples_per_cell_node,samples_per_cell_node=self.photosphere_samples_per_cell_node,continuous_emission_displayed=False,surface_samples_are_radiation_not_material=True,persistent_source_labels=True,bounded_model_paths=True,path_segments=0,path_pixel_step=self.path_pixel_step)
            return np.zeros_like(reference),mask,Teff.ravel()
        raw=np.stack([np.interp(np.log(np.maximum(Teff,300)),np.log(self.gas_lut[0]),self.gas_lut[1][:,j]) for j in range(3)],-1);color=warm_chromaticity(Teff)*(raw@Y)[...,None]/np.maximum(SIGMA*np.maximum(Teff,1)**4,1e-100)[...,None];eligible=np.flatnonzero((power>0)&(mu>0));P=power[eligible,None]*weights[None,:]*visible[eligible,None];P=np.repeat((P/self.photosphere_samples_per_cell_node)[...,None],self.photosphere_samples_per_cell_node,axis=2).ravel();colors=np.repeat(color[eligible,:,None,:],self.photosphere_samples_per_cell_node,axis=2).reshape(-1,3);pixel=2*L/self.W;image=np.zeros((self.H,self.W,3),np.float32);used=0.;timeline=self.bin_timeline(self.radiation['t']);duration=timeline[-1]-timeline[0];previous=None
        for k,playback in enumerate(timeline):
            tt=float(self.clock(playback));camera,inner_pixel=self.view(playback);world=self.photo_world(tt,eligible);current=world.reshape(-1,3)@(camera*(pixel/inner_pixel))
            if previous is not None:
                part,u=compact_path_emission(previous,current,P*(playback-timeline[k-1])/duration,colors,self.W,self.H,pixel,.64*self.H,self.path_pixel_step);image+=part;used+=u
            previous=current
        image/=4*(pixel*RM)**2;y=image@Y;image*=(.02*np.log1p(y/.02)/np.maximum(y,1e-100))[...,None]
        self.emitting_metrics=dict(allocated_physical_luminosity_W=float(power.sum()),view_isotropic_equivalent_source_W=float(P.sum()),view_power_in_native_footprints_W=float(used),last_interaction_weight_sum=float(weights.sum()),reprocessing_source_samples=len(dirs)*3*self.photosphere_samples_per_cell_node,active_positive_view_samples=len(P),samples_per_cell_node=self.photosphere_samples_per_cell_node,continuous_emission_displayed=False,surface_samples_are_radiation_not_material=True,persistent_source_labels=True,bounded_model_paths=True,path_segments=len(timeline)-1,path_pixel_step=self.path_pixel_step)
        return image,mask,Teff.ravel()
    def fine_at(self,playback,ids,source,colors,depth,outer_pixel,allocation):
        t=float(self.clock(playback));q=self.motion.at(t);camera,pixel=self.view(playback);camera=camera*(outer_pixel/pixel);offsets,_=self.micro.cloud(t,ids);world=q['x'][ids,None,:]+offsets;weights=self.micro.weights[self.micro.lookup[ids]];P=(source[:,None]*weights/weights.sum(axis=1)[:,None]).ravel()
        return compact_emission((world@camera).reshape(-1,3),P,np.repeat(colors,len(self.micro.nodes),axis=0),self.W,self.H,outer_pixel,.64*self.H)
    def bin_timeline(self,events):
        width=.5/self.micro_fps/self.micro_sample_count;lo=self.micro_playback-width/2;hi=self.micro_playback+width/2;tl=float(self.clock(lo));th=float(self.clock(hi));inside=np.asarray(events);inside=np.unique(inside[np.isfinite(inside)&(inside>tl)&(inside<th)]);extra=[brentq(lambda p:float(self.clock(p))-t,lo,hi,xtol=1e-14) for t in inside];return np.unique(np.r_[np.linspace(lo,hi,self.path_segments+1),extra])
    def cooling_positions(self,playback,ids,pixel):
        t=float(self.clock(playback));q=self.motion.at(t);camera,inner_pixel=self.view(playback);offsets,_=self.micro.cloud(t,ids);return (q['x'][ids,None,:]+offsets).reshape(-1,3)@(camera*(pixel/inner_pixel))
    def curved_cooling(self,t,q,ids,p,cov,source,colors,camera,depth,pixel):
        events=np.r_[self.micro.times,self.motion.birth[ids],self.motion.death[ids],self.motion.event_time];timeline=self.bin_timeline(events);dur=np.diff(timeline);states=[self.motion.at(float(self.clock((a+b)/2)))['live'][ids] for a,b in zip(timeline[:-1],timeline[1:])];active=np.asarray(states);normal=np.sum(dur[:,None]*active,axis=0);result=np.zeros((self.H,self.W,3),np.float32);used=0.;allocated=np.zeros(len(ids));weights=self.micro.weights[self.micro.lookup[ids]];weights=weights/weights.sum(axis=1)[:,None]
        if np.any(normal<=0):raise RuntimeError('Midpoint material must have a nonzero valid shutter interval')
        for k,(a,b) in enumerate(zip(timeline[:-1],timeline[1:])):
            good=active[k]
            if not np.any(good):continue
            power=source[good]*dur[k]/normal[good];P=(power[:,None]*weights[good]).ravel();start=self.cooling_positions(float(a),ids[good],pixel);end=self.cooling_positions(float(b),ids[good],pixel);image,u=compact_path_emission(start,end,P,np.repeat(colors[good],len(self.micro.nodes),axis=0),self.W,self.H,pixel,.64*self.H,self.path_pixel_step);result+=image;used+=u;allocated[good]+=power
        offsets,_=self.micro.cloud(t,ids);center=np.sum(offsets*weights[:,:,None],axis=1);self.micro_metrics=dict(parent_tracers=len(ids),micro_emitters=len(ids)*len(self.micro.nodes),emitters_per_parent=len(self.micro.nodes),source_power_W=float(allocated.sum()),expected_fine_power_W=float(source.sum()),actual_micro_fraction_of_fine=1.,particle_cooling_fraction=1.,smooth_fine_cooling_fraction=0.,broad_cooling_fraction=0.,max_parent_centroid_error_Rsun=float(np.max(abs(offsets.mean(axis=1)))),max_weighted_parent_centroid_error_Rsun=float(np.max(abs(center))),extra_mass=0.,per_frame_resampling=False,contours_used=False,physical_opacity_unchanged=True,micro_is_visual_quadrature_not_mass=True,gaussian_temporal_blur=False,shutter_fraction_frame=.5,common_exact_time_support=True,independent_orientation_transport=False,bounded_model_paths=True,path_segments=len(timeline)-1,path_pixel_step=self.path_pixel_step,max_parent_temporal_power_error=float(np.max(abs(allocated/np.maximum(source,1e-100)-1))))
        return result,used
    def render_linear(self,playback):
        if hasattr(self,'micro_metrics'):del self.micro_metrics
        image,mask=super().render_linear(playback);photo=self.layers['reprocessed_photosphere'];cool=self.layers['fine_cooling'];cY=np.maximum(0,cool@Y);pY=np.maximum(0,photo@Y);gas=np.maximum(0,image-self.layers['attenuated_stars']);self.particle_display=gas;self.base_display=self.layers['attenuated_stars']
        if 'microflow' not in self.records[-1]:
            self.records[-1]['microflow']=dict(parent_tracers=0,micro_emitters=0,emitters_per_parent=len(self.micro.nodes),source_power_W=0.,expected_fine_power_W=0.,particle_cooling_fraction=1.,smooth_fine_cooling_fraction=0.,broad_cooling_fraction=0.,max_parent_centroid_error_Rsun=0.,max_weighted_parent_centroid_error_Rsun=0.,extra_mass=0.,per_frame_resampling=False,contours_used=False,physical_opacity_unchanged=True,micro_is_visual_quadrature_not_mass=True,gaussian_temporal_blur=False,shutter_fraction_frame=.5,common_exact_time_support=True,independent_orientation_transport=False,bounded_model_paths=True,path_segments=0,path_pixel_step=self.path_pixel_step,max_parent_temporal_power_error=0.)
        self.layers.update(particle_display=self.particle_display,stellar_display=self.base_display,reprocessed_reference=self.reprocessed_reference)
        self.records[-1]['particle_only_cooling']=True;self.records[-1]['all_gas_emission_particle_sampled']=True;self.records[-1]['smooth_cooling_layer_max']=float(np.max(abs(self.layers['coarse_cooling'])));self.records[-1]['cooling_pregrade_luminance']=float(cY.sum());self.records[-1]['reprocessed_pregrade_luminance']=float(pY.sum());return image,mask
    def process_components(self,base,particles,mask):
        opt=self.optics;core=np.stack([opt.convolve(base[...,c],opt.core) for c in range(3)],-1);halo=np.stack([opt.convolve(base[...,c],opt.wing) for c in range(3)],-1);pcore=np.stack([opt.convolve(particles[...,c],opt.core) for c in range(3)],-1);phalo=np.stack([opt.convolve(particles[...,c],opt.wing) for c in range(3)],-1);eta=.035;ep=self.particle_bloom_fraction;physical=(1-eta)*core+eta*halo+(1-ep)*pcore+ep*phalo;total=base.sum(dtype='f8')+particles.sum(dtype='f8');opt.flux_errors.append(float(abs(physical.sum(dtype='f8')/max(total,1e-100)-1)));p=opt.pad;opt.crop_fractions.append(float(physical[p:p+self.H,p:p+self.W].sum()/max(physical.sum(),1e-100)));outside=1-gaussian_filter(np.pad(mask.astype('f4'),p),max(1.,self.W/180));hdr=np.maximum((physical+eta*2*halo*outside[...,None])[p:p+self.H,p:p+self.W],0);a=1.6*hdr;return a/(1+a.max(axis=-1))[...,None]
    def frame(self,playback,fps=30,temporal_samples=4):
        self.micro_sample_count=temporal_samples;self.micro_fps=fps;base=np.zeros((self.H,self.W,3),np.float32);particle=np.zeros_like(base);mask=np.zeros((self.H,self.W),np.float32)
        for k in range(temporal_samples):
            offset=((k+.5)/temporal_samples-.5)*.5/fps;_,m=self.render_linear(float(np.clip(playback+offset,0,32)));base+=self.base_display/temporal_samples;particle+=self.particle_display/temporal_samples;mask+=m/temporal_samples
        mapped=self.process_components(base,particle,mask);mapped*=1-self.title_outline[...,None]*.85;mapped=mapped*(1-self.title_core[...,None])+self.title_core[...,None]*.91;srgb=np.where(mapped<=.0031308,12.92*mapped,1.055*np.maximum(mapped,0)**(1/2.4)-.055);return np.uint16(np.clip(srgb,0,1)*65535+.5)
