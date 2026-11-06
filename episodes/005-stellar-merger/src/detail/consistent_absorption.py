"""Density-consistency experiment; visible material and kappa are unchanged."""
import numpy as np
from src.detail.material_only import MaterialOnlyRenderer
from src.detail.density_consistency import material_density_covariance,gaussian_columns
from src.detail.reprocessing import MM,RM

def validate_density_scale(value):
 value=float(value)
 if not np.isfinite(value) or value<=0:raise ValueError('Density width must be finite and strictly positive: width zero gives a singular zero-volume mass kernel, not an absorption-off setting')
 return value

class ConsistentAbsorptionRenderer(MaterialOnlyRenderer):
 density_neighbors=12
 legacy_isotropic=False
 co_spatial=False
 density_scale=1.
 def __init__(self,*args,**kwargs):
  super().__init__(*args,**kwargs);self.co_spatial=bool(self.micro.settings.get('cospatial_absorption',False));self.density_scale=validate_density_scale(self.micro.settings.get('density_scale',1.));self.density_profile=self.micro.settings.get('density_profile','gaussian')
 def stellar_optical_depth(self,t,q,ids,sigma,p,depth,camera,pixel):
  if self.co_spatial:validate_density_scale(self.density_scale)
  if not len(ids):return np.zeros(depth.shape)
  if self.co_spatial:
   k=np.clip(np.searchsorted(self.micro.times,t)-1,0,len(self.micro.times)-2);f=np.clip((t-self.micro.times[k])/(self.micro.times[k+1]-self.micro.times[k]),0,1);ii=self.micro.lookup[ids];C=((1-f)*self.micro.geometry[k,ii]+f*self.micro.geometry[k+1,ii])*self.density_scale**2;metrics={'co_spatial_density_emissivity_envelope_assumption':True,'density_scale':self.density_scale,'mass_per_kernel_unchanged':True,'kappa_unchanged':True}
  elif self.legacy_isotropic:C=np.eye(3)[None]*sigma[:,None,None]**2;metrics={'isotropic':True,'full_optical_tails':True}
  else:C,metrics=material_density_covariance(self.motion,q,ids,self.density_neighbors)
  self.last_density_world=C;cc=np.einsum('ij,njk,kl->nil',camera.T,C,camera)
  if np.min(np.linalg.eigvalsh(C))<1e-24:raise ValueError('Density width collapses a kernel below the reliable numerical scale; no silent variance floor is permitted')
  if self.density_profile=='compact':
   from src.detail.compact_density import compact_columns
   _,front,_=compact_columns(p,cc,np.full(len(ids),self.motion.mass),q['T'][ids],depth,self.W,self.H,pixel,.64*self.H,foreground_only=True)
  else:_,front,_=gaussian_columns(p,cc,np.full(len(ids),self.motion.mass),q['T'][ids],depth,self.W,self.H,pixel,.64*self.H,.01*MM/RM**2,foreground_only=True)
  metrics['density_profile']=self.density_profile;self.density_metrics=metrics;return (front*.01*MM/RM**2).astype('f4')
 def render_linear(self,p):
  im,mask=super().render_linear(p);self.records[-1]['density_reconstruction']=getattr(self,'density_metrics',{});self.records[-1]['microflow']['physical_opacity_unchanged']=False;self.records[-1].update(absorption_coefficient_unchanged=True,density_reconstruction_changed=True,prior_radiation_cache_is_archival_reference=True);return im,mask
