"""Positive covariance-matched compact density using No.004's Wendland profile.

A density-estimation closure, not new hydrodynamics. Kernel mass and covariance
are normalized analytically; the absorption coefficient is not altered.
"""
import numpy as np
from numba import njit
NODES,WEIGHTS=np.polynomial.legendre.leggauss(16)
I2=8/495
SUPPORT=np.sqrt(39/2)

@njit(cache=True)
def _segment(lo,hi,center,a,qmin):
 if hi<=lo:return 0.
 s=0.;half=(hi-lo)/2;mid=(hi+lo)/2
 for k in range(len(NODES)):
  z=mid+half*NODES[k];q=np.sqrt(max(0.,qmin+a*(z-center)**2))
  if q<1:s+=WEIGHTS[k]*(1-q)**6*(1+6*q+(35/3)*q*q)
 return half*s

@njit(cache=True)
def _integral(lo,hi,center,a,qmin):
 if lo<center<hi:return _segment(lo,center,center,a,qmin)+_segment(center,hi,center,a,qmin)
 return _segment(lo,hi,center,a,qmin)

@njit(cache=True)
def _columns(p,C,inv,det,mass,T,depth,W,H,pixel,cy,foreground_only):
 total=np.zeros((H,W),np.float64);front=total.copy();heat=total.copy();lx=0;hx=W-1;ly=0;hy=H-1
 if foreground_only:
  lx=W;hx=-1;ly=H;hy=-1
  for y in range(H):
   for x in range(W):
    if np.isfinite(depth[y,x]):lx=min(lx,x);hx=max(hx,x);ly=min(ly,y);hy=max(hy,y)
 for j in range(len(p)):
  Q=inv[j]/SUPPORT**2;cx=W/2+p[j,0]/pixel-.5;yy=cy-p[j,1]/pixel-.5;sx=SUPPORT*np.sqrt(C[j,0,0])/pixel;sy=SUPPORT*np.sqrt(C[j,1,1])/pixel;norm=mass[j]/(4*np.pi*I2*SUPPORT**3*np.sqrt(det[j]));a=Q[2,2]
  for y in range(max(ly,int(np.floor(yy-sy))),min(hy,int(np.ceil(yy+sy)))+1):
   for x in range(max(lx,int(np.floor(cx-sx))),min(hx,int(np.ceil(cx+sx)))+1):
    if foreground_only and not np.isfinite(depth[y,x]):continue
    dx=(x-cx)*pixel;dy=-(y-yy)*pixel;b=Q[2,0]*dx+Q[2,1]*dy;c=Q[0,0]*dx*dx+2*Q[0,1]*dx*dy+Q[1,1]*dy*dy;qmin=max(0.,c-b*b/a)
    if qmin>=1:continue
    center=-b/a;half=np.sqrt((1-qmin)/a);lo=center-half;hi=center+half
    if np.isfinite(depth[y,x]):
     split=min(hi,max(lo,depth[y,x]-p[j,2]));v=norm*_integral(split,hi,center,a,qmin);back=norm*_integral(lo,split,center,a,qmin);total[y,x]+=v+back;front[y,x]+=v;heat[y,x]+=v*T[j]
    else:total[y,x]+=norm*_integral(lo,hi,center,a,qmin)
 return total,front,heat

def compact_columns(p,cov,mass,T,depth,W,H,pixel,cy,foreground_only=True):
 C=cov.astype('f8').copy();C[:,0,0]+=pixel*pixel/12;C[:,1,1]+=pixel*pixel/12
 return _columns(p,C,np.linalg.inv(C),np.linalg.det(C),mass,T,depth,W,H,pixel,cy,foreground_only)

@njit(cache=True)
def _radial(x,Q,norm,dirs,R):
 tau=np.zeros(len(dirs))
 for k in range(len(dirs)):
  d=dirs[k];v=0.
  for j in range(len(x)):
   a=np.dot(d,Q[j]@d);b=-np.dot(d,Q[j]@x[j]);c=np.dot(x[j],Q[j]@x[j]);qmin=max(0.,c-b*b/a)
   if qmin>=1:continue
   center=-b/a;half=np.sqrt((1-qmin)/a);v+=norm[j]*_integral(max(R,center-half),center+half,center,a,qmin)
  tau[k]=v
 return tau

def compact_radial_columns(x,C,mass,dirs,R,kappa):
 C=C.astype('f8');Q=np.linalg.inv(C)/SUPPORT**2;norm=mass*kappa/(4*np.pi*I2*SUPPORT**3*np.sqrt(np.linalg.det(C)));return _radial(x,Q,norm,dirs,R)
