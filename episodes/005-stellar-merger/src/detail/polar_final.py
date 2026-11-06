"""Final more-polar camera and one controlled compact-density alternative."""
import numpy as np
from src.detail.remnant_focus import RemnantFocusRenderer
from src.detail.reprocessing import MM,RM

def start_ease(u):return 6*u**3-8*u**4+3*u**5

def stop_ease(u):return u+4*u**3-7*u**4+3*u**5

class PolarFinalRenderer(RemnantFocusRenderer):
 def __init__(self,*a,**k):
  super().__init__(*a,**k);self.camera_offset=float(self.micro.settings.get('camera_inclination_offset_degrees',-10.));self.kernel_profile=self.micro.settings.get('density_profile','gaussian');self.kernel_order=int(self.micro.settings.get('wendland_quadrature_order',16));self.kernel_rays=int(self.micro.settings.get('wendland_pixel_ray_factor',2));self.density_metrics=dict(density_scale=self.density_scale,density_profile=self.kernel_profile,mass_per_kernel_unchanged=True,kappa_unchanged=True)
 def camera_geometry(self,p,t=None):
  L,inc,az,cam=super().camera_geometry(p,t);t=float(self.clock(p)) if t is None else t;z=min(t,300.);k=.12/175
  if z<=215:d=0.
  elif z<225:d=10*start_ease((z-215)/10)
  elif z<=290:d=z-215
  else:d=75+10*stop_ease((z-290)/10)
  inc=.82-k*d+np.deg2rad(getattr(self,'camera_offset',-10.))
  if p<22:az=-np.pi/2+.35*p/32
  elif p<24:az=-np.pi/2+.35*22/32+.35*2/32*stop_ease((p-22)/2)
  else:az=-np.pi/2+.35*.75
  eye=np.array([np.sin(inc)*np.cos(az),np.sin(inc)*np.sin(az),np.cos(inc)]);right=np.array([-np.sin(az),np.cos(az),0]);up=np.cross(eye,right);return L,inc,az,np.column_stack([right,up,eye])
 def stars_native(self,*a,**k):
  q=super().stars_native(*a,**k);self.true_depth=q[1];self.true_mask=q[2];return q
 def compact_field(self,pos,C,depth,pixel):
  from src.detail.wendland_c2 import foreground_columns
  n=self.kernel_rays;trans=np.zeros((self.H,self.W),np.float64)
  for iy in range(n):
   for ix in range(n):
    shift=pos.copy();shift[:,0]-=((ix+.5)/n-.5)*pixel;shift[:,1]+=((iy+.5)/n-.5)*pixel
    tau=foreground_columns(shift,C,np.full(len(pos),self.motion.mass),depth,self.W,self.H,pixel,.64*self.H,self.kernel_order)*.01*MM/RM**2;trans+=np.exp(-np.minimum(tau,100))/(n*n)
  return (-np.log(np.maximum(trans,np.exp(-100)))).astype('f4')
 def compact_tau(self,p):
  t=float(self.clock(p));q=self.motion.at(t);ids=np.flatnonzero(q['live']);L,_,_,cam=self.camera_geometry(p,t)
  if not len(ids):return np.zeros((self.H,self.W),np.float32)
  C=np.einsum('ij,njk,kl->nil',cam.T,self.last_density_world,cam)
  return self.compact_field(q['x'][ids]@cam,C,self.true_depth,2*L/self.W)
 def stellar_optical_depth(self,t,q,ids,sigma,pos,depth,cam,pixel):
  # Establish the same physical covariance by the existing conservative path.
  base=super().stellar_optical_depth(t,q,ids,sigma,pos,depth,cam,pixel)
  if self.kernel_profile!='wendland_c2' or not len(ids):return base
  C=np.einsum('ij,njk,kl->nil',cam.T,self.last_density_world,cam);tau=self.compact_field(pos,C,depth,pixel);self.density_metrics.update(density_profile='wendland_c2',mass_normalized=True,covariance_matched=True,compact_support_sigma=float(np.sqrt(15)),quadrature_order=self.kernel_order,pixel_ray_factor=self.kernel_rays,effective_pixel_optical_depth=True)
  return tau.astype('f4')
 def render_linear(self,p):
  q=super().render_linear(p);self.records[-1].update(camera_inclination_offset_degrees=self.camera_offset,camera_angular_easing=True,opacity_profile=self.kernel_profile);return q
