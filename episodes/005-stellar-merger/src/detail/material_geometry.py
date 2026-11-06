"""Persistent emission sampling using material-neighborhood support geometry."""
import numpy as np
from src.detail.particle_cooling import ParticleFlow
from src.detail.microflow import sqrt_spd,ROOT

class MaterialGeometryFlow(ParticleFlow):
 def __init__(self,reference,motion):
  super().__init__(reference,motion);self.geometry=np.load(ROOT/self.settings['material_geometry_support_path'],mmap_mode='r')
  # Weighted, immutable 3-D quadrature. Fixed positive weights and antithetic
  # pairs preserve centroid; whitening restores the support covariance.
  z=self.labels[:,:,:3];w=self.weights/self.weights.sum(axis=1)[:,None];mean=np.einsum('nm,nmj->nj',w,z);z=z-mean[:,None];C=np.einsum('nm,nmi,nmj->nij',w,z,z);val,vec=np.linalg.eigh(C);inv=np.einsum('nij,nkj->nik',vec/np.sqrt(val)[:,None,:],vec);self.material_labels=np.einsum('nij,nmj->nmi',inv,z)
 def cloud(self,t,ids):
  k=np.clip(np.searchsorted(self.times,t)-1,0,len(self.times)-2);f=np.clip((t-self.times[k])/(self.times[k+1]-self.times[k]),0,1);ii=self.lookup[ids]
  C=(1-f)*self.geometry[k,ii]+f*self.geometry[k+1,ii];off=np.einsum('nij,nmj->nmi',sqrt_spd(C),self.material_labels[ii]);return off,None
