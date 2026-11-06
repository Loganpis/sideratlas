"""Minimal correction: one continuous support for smooth and fine emissivity.

There is no independently integrated orientation or derivative-based blur.
The existing resolved trajectories and finite-deformation estimator are reused.
"""
import numpy as np
from src.detail.microflow import MicroFlow,MicroRenderer,curved_basis
from src.detail.transported_microflow import TransportedRenderer
from src.detail.fields import coarse_sigma,fine_covariances,splat_emission

class SynchronizedFlow(MicroFlow):
    def __init__(self,reference,motion):
        self.__dict__.update(reference.__dict__);self.motion=motion;self.reference=reference
    def at(self,t,ids):
        C,R,acc,v,ii=self.reference.at(t,ids)
        return C,np.tile(np.eye(3),(len(ids),1,1)),acc,v,ii

class SynchronizedRenderer(TransportedRenderer):
    def __init__(self,source,width=1080,settings_path=None):
        MicroRenderer.__init__(self,source,width,settings_path=settings_path);self.control_micro=self.micro;self.micro=SynchronizedFlow(self.control_micro,self.motion);self.shutter_order=3
    def fine_at(self,playback,ids,source,colors,depth,outer_pixel,allocation):
        t=float(self.clock(playback));C,R,acc,v,ii=self.micro.at(t,ids);q=self.motion.at(t);sub,axis,L,curv=curved_basis(C,v,acc);camera,pixel=self.view(playback);camera=camera*(outer_pixel/pixel);sub=np.einsum('ij,njk,kl->nil',camera.T,sub,camera);p=q['x'][ids]@camera;result=np.zeros((self.H,self.W,3),np.float32);used=0.;nodes,w=np.polynomial.hermite.hermgauss(5);nodes*=np.sqrt(2);w/=np.sqrt(np.pi)
        for z,weight in zip(nodes,w):
            shift=np.sqrt(.85)*L[:,None]*z*axis+.5*.85*L[:,None]**2*(z*z-1)*curv;im,u=splat_emission(p+shift@camera,sub,.85*(1-allocation)*source*weight,colors,depth,self.W,self.H,outer_pixel,.64*self.H,False);result+=im;used+=u
        offsets,K=self.micro.cloud(t,ids);world=q['x'][ids,None,:]+offsets;CC=np.einsum('ij,nmjk,kl->nmil',camera.T,K,camera);weights=self.micro.weights[ii];P=(.85*self.micro.settings['fine_component_micro_fraction']*source[:,None]/len(self.micro.nodes)*weights).ravel();im,u=splat_emission((world@camera).reshape(-1,3),CC.reshape(-1,3,3),P,np.repeat(colors,len(self.micro.nodes),axis=0),depth,self.W,self.H,outer_pixel,.64*self.H,False);return result+im,used+u
    def curved_cooling(self,*args):
        image,used=super().curved_cooling(*args);self.micro_metrics.update(common_smooth_and_micro_transport=False,common_exact_time_support=True,independent_orientation_transport=False,cached_parent_velocity_used_for_translation=False,cached_support_basis_used=True)
        return image,used
