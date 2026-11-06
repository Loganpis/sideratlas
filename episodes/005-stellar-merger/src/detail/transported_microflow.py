"""Common material transport and bounded, synchronized shutter quadrature."""
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from src.detail.microflow import MicroRenderer,ROOT
from src.detail.fields import splat_emission

class TransportedFlow:
    def __init__(self,reference,path=None):
        self.reference=reference;self.settings=reference.settings;self.weights=reference.weights;self.nodes=reference.nodes;self.lookup=reference.lookup
        a=np.load(path or ROOT/'output/transported-microflow-limited/transport.npz');self.data={k:a[k] for k in a.files};self.times=self.data['t'];self.ids=self.data['ids']
    def matrix(self,t,ids):
        if not self.times[0]<=t<=self.times[-1]:raise RuntimeError('Transport outside qualified interval')
        k=int(np.clip(np.searchsorted(self.times,t,side='right')-1,0,len(self.times)-2));ii=self.lookup[ids];A=self.data['A'][k,ii];G=expm(A*(t-self.times[k]))@self.data['G'][k,ii]
        return G,A,ii
    def cloud(self,t,ids):
        G,A,ii=self.matrix(t,ids);offset=np.einsum('nij,nmj->nmi',G,self.data['offset0'][ii]);K=np.einsum('nij,nmjk,nlk->nmil',G,self.data['K0'][ii],G)
        return offset,K

class TransportedRenderer(MicroRenderer):
    def __init__(self,source,width=1080):
        super().__init__(source,width);self.control_micro=self.micro;self.micro=TransportedFlow(self.control_micro);self.shutter_order=3
    def view(self,playback):
        L,inc,az,camera=self.camera_geometry(playback)
        return camera,2*L/self.W
    def fine_at(self,playback,ids,source,colors,depth,outer_pixel,allocation):
        t=float(self.clock(playback));q=self.motion.at(t);camera,pixel=self.view(playback);scale=outer_pixel/pixel;camera=camera*scale;G,A,ii=self.micro.matrix(t,ids);d=self.micro.data;sub=G@d['sub0'][ii]@G.transpose(0,2,1);sub=np.einsum('ij,njk,kl->nil',camera.T,sub,camera);axis=np.einsum('nij,nj->ni',G,d['axis0'][ii]);curv=np.einsum('nij,nj->ni',G,d['curvature0'][ii]);L=d['length0'][ii];p=q['x'][ids]@camera;result=np.zeros((self.H,self.W,3),np.float32);used=0.;nodes,w=np.polynomial.hermite.hermgauss(5);nodes*=np.sqrt(2);w/=np.sqrt(np.pi)
        for z,weight in zip(nodes,w):
            shift=np.sqrt(.85)*L[:,None]*z*axis+.5*.85*L[:,None]**2*(z*z-1)*curv;im,u=splat_emission(p+shift@camera,sub,.85*(1-allocation)*source*weight,colors,depth,self.W,self.H,outer_pixel,.64*self.H,False);result+=im;used+=u
        offsets,K=self.micro.cloud(t,ids);world=q['x'][ids,None,:]+offsets;CC=np.einsum('ij,nmjk,kl->nmil',camera.T,K,camera);weights=self.micro.weights[ii];P=(.85*self.micro.settings['fine_component_micro_fraction']*source[:,None]/len(self.micro.nodes)*weights).ravel();col=np.repeat(colors,len(self.micro.nodes),axis=0);im,u=splat_emission((world@camera).reshape(-1,3),CC.reshape(-1,3,3),P,col,depth,self.W,self.H,outer_pixel,.64*self.H,False);result+=im;used+=u
        return result,used
    def curved_cooling(self,t,q,ids,p,cov,source,colors,camera,depth,pixel):
        fraction=self.micro.settings['fine_component_micro_fraction'];weights=self.micro.weights[self.micro.lookup[ids]];allocation=fraction*weights.mean(axis=1);nodes,w=np.polynomial.legendre.leggauss(self.shutter_order);w=w/2;bin_width=.5/self.micro_fps/self.micro_sample_count;sample_playback=self.micro_playback+nodes*bin_width/2
        states=[self.motion.at(float(self.clock(playback))) for playback in sample_playback];live=np.stack([s['live'][ids] for s in states]);normal=np.sum(w[:,None]*live,axis=0);result=np.zeros((self.H,self.W,3),np.float32);used=0.
        for k,playback in enumerate(sample_playback):
            good=live[k]&(normal>0)
            if not np.any(good):continue
            power=source[good]*w[k]/normal[good];im,u=self.fine_at(float(playback),ids[good],power,colors[good],depth,pixel,allocation[good]);result+=im;used+=u
        missing=normal==0
        if np.any(missing):
            im,u=self.fine_at(self.micro_playback,ids[missing],source[missing],colors[missing],depth,pixel,allocation[missing]);result+=im;used+=u
        offsets,K=self.micro.cloud(t,ids);centroid=np.sum(offsets*weights[:,:,None],axis=1)/weights.sum(axis=1)[:,None];microP=.85*np.sum(allocation*source)
        self.micro_metrics=dict(parent_tracers=len(ids),micro_emitters=len(ids)*len(self.micro.nodes),emitters_per_parent=len(self.micro.nodes),source_power_W=float(microP+.85*np.sum((1-allocation)*source)),expected_fine_power_W=float(.85*source.sum()),max_parent_centroid_error_Rsun=float(np.max(abs(offsets.mean(axis=1)))),max_weighted_parent_centroid_error_Rsun=float(np.max(abs(centroid))),actual_micro_fraction_of_fine=float(np.sum(allocation*source)/max(source.sum(),1e-100)),per_frame_resampling=False,contours_used=False,physical_opacity_unchanged=True,micro_is_visual_quadrature_not_mass=True,common_smooth_and_micro_transport=True,gaussian_temporal_blur=False,bounded_shutter_nodes=self.shutter_order,shutter_fraction_frame=.5,short_lifetime_midpoint_fallback_power_fraction=float(source[missing].sum()/max(source.sum(),1e-100)),parent_trajectory_unchanged=True)
        return result,used
