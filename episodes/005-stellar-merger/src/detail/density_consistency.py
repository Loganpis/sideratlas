"""Conservative local-geometry density estimates and Gaussian foreground columns.

Post-processing only: not fluid pressure, new mass, or new trajectories. Kappa
is unchanged. Bandwidth is measured from parent neighborhoods, not fitted to
make the image transparent or copied from the emissivity contrast controls.
"""
import numpy as np
from scipy.spatial import cKDTree
from numba import njit
from math import erfc

def material_density_covariance(motion,q,ids,neighbors=12):
 x=q['x'][ids];v=q['v'][ids];n=len(ids)
 if n<3:
  from src.detail.fields import coarse_sigma
  s=coarse_sigma(x,q['radius'],q['donor']);return np.eye(3)[None]*s[:,None,None]**2,dict(neighbors=neighbors,fallback=n)
 distance,near=cKDTree(x).query(x,k=min(neighbors+1,n),workers=2);dx=x[near]-x[:,None];dv=v[near]-v[:,None];birth=motion.birth[ids];same=abs(birth[near]-birth[:,None])<8.;scale=np.maximum(.04,.45*np.linalg.norm(v,axis=1));h=np.maximum(distance[:,-1],1e-8)
 w=np.exp(-.5*(distance/h[:,None])**2)*same*np.exp(-.5*np.sum(dv*dv,axis=-1)/scale[:,None]**2);w[:,0]=0;sw=w.sum(axis=1);raw=np.einsum('nk,nki,nkj->nij',w,dx,dx)/np.maximum(sw,1e-30)[:,None,None];neff=sw**2/np.maximum(np.sum(w*w,axis=1),1e-30);raw*=np.maximum(neff,2)[:,None,None]**(-2/7)
 eig,vec=np.linalg.eigh(raw);resolution=(np.maximum(distance[:,1],1e-6)*.25)**2;floor=np.maximum(resolution,eig[:,-1]/16);eig=np.maximum(eig,floor[:,None]);bad=np.sum(w>.05,axis=1)<3
 from src.detail.fields import coarse_sigma
 old=coarse_sigma(x,q['radius'],q['donor']);eig[bad]=old[bad,None]**2
 C=np.einsum('nij,nkj->nik',vec*eig[:,None,:],vec)
 return C,dict(neighbors=neighbors,fallback=int(bad.sum()),median_effective_neighbors=float(np.median(neff)),axis_ratio_quantiles=np.quantile(np.sqrt(eig[:,-1]/eig[:,0]),[.1,.5,.9]).tolist(),sigma_quantiles_Rsun=np.quantile(np.sqrt(eig),[.1,.5,.9]).tolist(),kappa_unchanged=True,mass_per_kernel_unchanged=True)

@njit(cache=True)
def gaussian_columns(p,cov,mass,T,depth,W,H,pixel,cy,kappa,tau_tolerance=1e-8,legacy_cut=False,foreground_only=False):
 total=np.zeros((H,W),np.float64);front=total.copy();heat=total.copy();loX=0;hiX=W-1;loY=0;hiY=H-1
 if foreground_only:
  loX=W;hiX=-1;loY=H;hiY=-1
  for y in range(H):
   for x in range(W):
    if np.isfinite(depth[y,x]):loX=min(loX,x);hiX=max(hiX,x);loY=min(loY,y);hiY=max(hiY,y)
 for j in range(len(p)):
  C=cov[j];a=C[0,0]+pixel*pixel/12;b=C[0,1];c=C[1,1]+pixel*pixel/12;det=max(a*c-b*b,1e-30);ia=c/det;ib=-b/det;ic=a/det;cx=W/2+p[j,0]/pixel-.5;yy=cy-p[j,1]/pixel-.5;peak=mass[j]/(2*np.pi*np.sqrt(det));limit=16. if legacy_cut else max(16.,2*np.log(max(peak*kappa/tau_tolerance,1.)));sx=np.sqrt(a*limit)/pixel;sy=np.sqrt(c*limit)/pixel;zx=C[2,0]*ia+C[2,1]*ib;zy=C[2,0]*ib+C[2,1]*ic;sz=np.sqrt(max(C[2,2]-zx*C[0,2]-zy*C[1,2],1e-30))
  for y in range(max(loY,int(np.floor(yy-sy))),min(hiY,int(np.ceil(yy+sy)))+1):
   for x in range(max(loX,int(np.floor(cx-sx))),min(hiX,int(np.ceil(cx+sx)))+1):
    if foreground_only and not np.isfinite(depth[y,x]):continue
    dx=(x-cx)*pixel;dy=-(y-yy)*pixel;e=ia*dx*dx+2*ib*dx*dy+ic*dy*dy
    if e>limit:continue
    value=peak*np.exp(-.5*e);total[y,x]+=value
    if np.isfinite(depth[y,x]):
     mean=p[j,2]+zx*dx+zy*dy;fraction=.5*erfc((depth[y,x]-mean)/(np.sqrt(2)*sz));v=value*fraction;front[y,x]+=v;heat[y,x]+=v*T[j]
 return total,front,heat

@njit(cache=True)
def radial_stellar_columns(x,cov,mass,dirs,R,kappa):
 """Exact Gaussian columns along outward radial rays from a spherical source.

 This is an irradiation/extinction bound, not local thermalization, scattering,
 or an observer-directed reprocessing calculation.
 """
 cf=cov.astype(np.float64);inv=np.empty_like(cf);det=np.empty(len(x));v=np.empty_like(x);c=np.empty(len(x))
 for j in range(len(x)):
  inv[j]=np.linalg.inv(cf[j]);det[j]=np.linalg.det(cf[j]);v[j]=inv[j]@x[j];c[j]=np.dot(x[j],v[j])
 result=np.zeros(len(dirs))
 for k in range(len(dirs)):
  d=dirs[k];value=0.
  for j in range(len(x)):
   a=np.dot(d,inv[j]@d);b=np.dot(d,v[j]);e=max(0.,c[j]-b*b/a)
   if e>140:continue
   value+=mass[j]*kappa/(4*np.pi*np.sqrt(det[j]*a))*np.exp(-.5*e)*erfc(np.sqrt(a/2)*(R-b/a))
  result[k]=value
 return result
