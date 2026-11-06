"""Conservative, persistent emissivity quadrature inside CURRENT flow kernels.

Samples are visualization labels, not newly calculated fluid particles or mass.
The smooth physical opacity and radiation fields remain unchanged.
"""
from pathlib import Path
import json,numpy as np
from scipy.special import ndtri
from scipy.stats import qmc
from scipy.spatial.transform import Rotation
from scipy.optimize import brentq
from src.detail.renderer import DetailRenderer
from src.detail.fields import splat_emission,warm_chromaticity

ROOT=Path(__file__).resolve().parents[2]

def balanced_nodes(count=2048,seed=5005004):
    u=qmc.Sobol(4,scramble=True,seed=seed).random_base2(int(np.log2(count//2)))
    z=ndtri(np.clip(u,1e-7,1-1e-7));z=np.r_[z,-z]
    # Exact finite-sample mean and covariance; antithetic odd moments vanish.
    C=z.T@z/len(z);val,vec=np.linalg.eigh(C);return z@((vec/np.sqrt(val))@vec.T)

def sqrt_spd(C,floor=1e-12):
    val,vec=np.linalg.eigh(C);return np.einsum('nij,nkj->nik',vec*np.sqrt(np.maximum(val,floor))[:,None,:],vec)

def curved_basis(cov,v,acc):
    speed=np.linalg.norm(v,axis=1);axis=v/np.maximum(speed[:,None],1e-12)
    normal=acc-axis*np.sum(acc*axis,axis=1)[:,None];curvature=normal/np.maximum(speed[:,None]**2,.01)
    length=np.sqrt(np.einsum('ni,nij,nj->n',axis,cov,axis));bend=np.linalg.norm(curvature,axis=1)*length;curvature*=np.minimum(1,.35/np.maximum(bend,1e-12))[:,None]
    sub=cov-.85*length[:,None,None]**2*axis[:,:,None]*axis[:,None,:]
    val,vec=np.linalg.eigh(sub);val=np.maximum(val,.025**2);sub=np.einsum('nij,nkj->nik',vec*np.sqrt(val)[:,None,:],vec*np.sqrt(val)[:,None,:])
    return sub,axis,length,curvature

class MicroFlow:
    def __init__(self,path=None,settings_path=None):
        self.settings=json.loads(Path(settings_path or ROOT/'calculations/visualization-v1/microflow-display.json').read_text());support=Path(path or ROOT/self.settings['material_support_source'])
        if support.is_dir():
            manifest=json.loads((support/'manifest.json').read_text());self.data={k:np.load(support/f,mmap_mode='r') for k,f in manifest['arrays'].items()};self.identity_orientation=manifest.get('identity_orientation',False)
        else:
            a=np.load(support);self.data={k:a[k] for k in a.files};self.identity_orientation=False
        self.times=self.data['t'];self.ids=self.data['ids'];self.lookup=np.full(8192,-1,int);self.lookup[self.ids]=np.arange(len(self.ids));count=self.settings['emitters_per_parent'];self.nodes=balanced_nodes(count)
        if self.settings.get('nested_512_labels'):
            if count!=1024:raise ValueError('Nested production sampling requires 1024 nodes')
            self.nodes=np.concatenate([balanced_nodes(512),balanced_nodes(512,seed=5005105)])
        self.labels=[]
        for j in self.ids:
            # Immutable material-label transform, never reseeded by frame.
            rng=np.random.default_rng(5005004+int(j)*97);Q,_=np.linalg.qr(rng.normal(size=(4,4)));self.labels.append(self.nodes@Q)
        self.labels=np.asarray(self.labels)
        self.weights=[]
        for z in self.labels:
            f=1/(1+np.exp(np.clip(1.5*(np.sum(z[:,:3]**2,axis=1)-5),-50,50)));d=z[:,3]**2-1
            # Positive material taper with exact longitudinal second moment:
            # suppress isolated tails without shifting curved source centers.
            lam=brentq(lambda l:np.sum(f*np.exp(l*d)*d),-1.,1.);w=f*np.exp(lam*d);w/=max(1.,float(w.max()));self.weights.append(w)
        self.weights=np.asarray(self.weights)
    def at(self,t,ids):
        if not self.times[0]<=t<=self.times[-1]:raise RuntimeError('Microflow cache covers only requested comparison interval')
        k=np.clip(np.searchsorted(self.times,t)-1,0,len(self.times)-2);f=(t-self.times[k])/(self.times[k+1]-self.times[k]);ii=self.lookup[ids]
        if np.any(ii<0):raise RuntimeError('No material support for a live parent')
        C=(1-f)*self.data['cov'][k,ii]+f*self.data['cov'][k+1,ii]
        if self.identity_orientation:R=np.tile(np.eye(3),(len(ids),1,1))
        else:
            qq0=self.data['quaternion'][k,ii];qq1=self.data['quaternion'][k+1,ii];qq1=qq1*np.where(np.sum(qq0*qq1,axis=1)<0,-1,1)[:,None];qq=(1-f)*qq0+f*qq1;qq/=np.linalg.norm(qq,axis=1,keepdims=True);R=Rotation.from_quat(qq).as_matrix()
        acc=(1-f)*self.data['acc'][k,ii]+f*self.data['acc'][k+1,ii];v=(1-f)*self.data['velocity'][k,ii]+f*self.data['velocity'][k+1,ii]
        return C,R,acc,v,ii
    def cloud(self,t,ids,return_footprints=True):
        C,R,acc,v,ii=self.at(t,ids);sub,axis,length,curvature=curved_basis(C,v,acc)
        # Fixed physical display size at 15 s, then local material expansion.
        ratio=(np.linalg.det(C)/np.maximum(self.data['reference_det'][ii],1e-100))**(1/6);ratio=np.clip(ratio,.5,3.)
        minor=self.settings['reference_pixel_Rsun']*self.settings['micro_sigma_transverse_px']*ratio;major=self.settings['reference_pixel_Rsun']*self.settings['micro_sigma_longitudinal_px']*ratio
        K=minor[:,None,None]**2*np.eye(3)+(major*major-minor*minor)[:,None,None]*axis[:,:,None]*axis[:,None,:]
        cap=np.minimum(1.,.4*np.linalg.eigvalsh(sub)[:,0]/np.maximum(np.linalg.eigvalsh(K)[:,-1],1e-100));K*=cap[:,None,None]
        z=self.labels[ii];r2=np.sum(z*z,axis=-1);density_width=1+2*np.clip((r2-3)/8,0,1)**2
        # Sparse Gaussian tails get broader, dimmer footprints, avoiding dust
        # beads while dense material retains small-pixel granularity. Labels
        # and width factors stay fixed in material coordinates.
        average_width2=np.mean(density_width[:,:512]**2,axis=1) if self.settings.get('nested_512_labels') else np.mean(density_width**2,axis=1);cap=np.minimum(1.,.4*np.linalg.eigvalsh(sub)[:,0]/np.maximum(np.linalg.eigvalsh(K)[:,-1]*average_width2,1e-100));K*=cap[:,None,None]
        averageK=K*average_width2[:,None,None]
        # Remove kernel variance from center distribution, preserving envelope
        # moments rather than broadening or adding mass to CURRENT support.
        B=sqrt_spd(sub-averageK);offset=np.einsum('nij,nmj->nmi',B@R,z[:,:,:3]);s=np.sqrt(.85)*length[:,None]*z[:,:,3];offset+=s[:,:,None]*axis[:,None,:]+.5*.85*length[:,None,None]**2*(z[:,:,3]**2-1)[:,:,None]*curvature[:,None,:]
        offset-=offset.mean(axis=1,keepdims=True)
        return offset,K[:,None]*density_width[:,:,None,None]**2 if return_footprints else None

class MicroRenderer(DetailRenderer):
    def __init__(self,source,width=1080,settings_path=None):
        super().__init__(source,width);self.micro=MicroFlow(settings_path=settings_path);self.micro_sample_count=4;self.micro_fps=30
    def frame(self,playback,fps=30,temporal_samples=4):
        self.micro_sample_count=temporal_samples;self.micro_fps=fps
        return super().frame(playback,fps,temporal_samples)
    def render_linear(self,playback):
        self.micro_playback=playback;image,mask=super().render_linear(playback)
        if hasattr(self,'micro_metrics'):self.records[-1]['microflow']=self.micro_metrics
        return image,mask
    def curved_cooling(self,t,q,ids,p,cov,source,colors,camera,depth,pixel):
        settings=self.micro.settings;fraction=settings['fine_component_micro_fraction'];weights=self.micro.weights[self.micro.lookup[ids]];allocated=fraction*weights.mean(axis=1);smooth,used=super().curved_cooling(t,q,ids,p,cov,source*(1-allocated),colors,camera,depth,pixel)
        offsets,K=self.micro.cloud(t,ids);world=q['x'][ids,None,:]+offsets
        eps=.01;lo=max(self.micro.times[0],t-eps);hi=min(self.micro.times[-1],t+eps);om,_=self.micro.cloud(lo,ids);op,_=self.micro.cloud(hi,ids);velocity=q['v'][ids,None,:]+(op-om)/(hi-lo)
        # Linear motion integration within each of the unchanged four shutter
        # bins avoids blinking/dashed sampling of subpixel moving emitters.
        dt=float(self.clock.derivative()(self.micro_playback))*.5/self.micro_fps/self.micro_sample_count
        vc=velocity@camera;CC=np.einsum('ij,nmjk,kl->nmil',camera.T,K,camera)+vc[:,:,:,None]*vc[:,:,None,:]*(dt*dt/12)
        count=len(self.micro.nodes);P=(.85*fraction*source[:,None]/count*weights).ravel();col=np.repeat(colors,count,axis=0);xyz=(world@camera).reshape(-1,3);CC=CC.reshape(-1,3,3)
        detail,detail_used=splat_emission(xyz,CC,P,col,depth,self.W,self.H,pixel,.64*self.H,False);image=smooth+detail
        self.micro_metrics=dict(parent_tracers=len(ids),micro_emitters=len(xyz),emitters_per_parent=count,micro_fraction_of_fine_upper_bound=fraction,actual_micro_fraction_of_fine=float(np.sum(allocated*source)/max(source.sum(),1e-100)),source_power_W=float(P.sum()+.85*np.sum((1-allocated)*source)),expected_fine_power_W=float(.85*source.sum()),max_parent_centroid_error_Rsun=float(np.max(abs(offsets.mean(axis=1)))),max_weighted_parent_centroid_error_Rsun=float(np.max(abs(np.sum(offsets*weights[:,:,None],axis=1)/weights.sum(axis=1)[:,None]))),max_single_emitter_fraction=float(1/count),static_seed_per_parent=True,per_frame_resampling=False,contours_used=False,physical_opacity_unchanged=True,linear_subshutter_motion_blur=True,subshutter_duration_code=dt,micro_is_visual_quadrature_not_mass=True,density_adaptive_material_footprint=True,unused_tail_power_returned_to_smooth_field=True)
        return image,used+detail_used
