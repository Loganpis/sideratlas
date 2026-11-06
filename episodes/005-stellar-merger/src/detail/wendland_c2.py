"""Mass- and covariance-matched Wendland C2 opacity closure (not new hydro).

3D profile (1-q)^4(1+4q), q<1. I2=1/42, I4=1/210;
H=sqrt(15) gives unit per-axis covariance in whitened coordinates.
"""
import numpy as np
from numba import njit
SUPPORT=np.sqrt(15.)

@njit(cache=True)
def segment(lo,hi,center,a,qmin,nodes,weights):
 if hi<=lo:return 0.
 half=.5*(hi-lo);mid=.5*(hi+lo);value=0.
 for k in range(len(nodes)):
  z=mid+half*nodes[k];q=np.sqrt(max(0.,qmin+a*(z-center)**2))
  if q<1:value+=weights[k]*(1-q)**4*(1+4*q)
 return half*value

@njit(cache=True)
def integral(lo,hi,center,a,qmin,nodes,weights):
 if lo<center<hi:return segment(lo,center,center,a,qmin,nodes,weights)+segment(center,hi,center,a,qmin,nodes,weights)
 return segment(lo,hi,center,a,qmin,nodes,weights)

@njit(cache=True)
def columns_core(p,C,inv,det,mass,depth,W,H,pixel,cy,nodes,weights):
 front=np.zeros((H,W),np.float64);yy,xx=np.nonzero(np.isfinite(depth));lx=xx.min();hx=xx.max();ly=yy.min();hy=yy.max()
 for j in range(len(p)):
  Q=inv[j]/15.;cx=W/2+p[j,0]/pixel-.5;yc=cy-p[j,1]/pixel-.5;sx=SUPPORT*np.sqrt(C[j,0,0])/pixel;sy=SUPPORT*np.sqrt(C[j,1,1])/pixel;norm=mass[j]*21/(2*np.pi*SUPPORT**3*np.sqrt(det[j]));a=Q[2,2]
  for y in range(max(ly,int(np.floor(yc-sy))),min(hy,int(np.ceil(yc+sy)))+1):
   for x in range(max(lx,int(np.floor(cx-sx))),min(hx,int(np.ceil(cx+sx)))+1):
    if not np.isfinite(depth[y,x]):continue
    dx=(x-cx)*pixel;dy=-(y-yc)*pixel;b=Q[2,0]*dx+Q[2,1]*dy;c=Q[0,0]*dx*dx+2*Q[0,1]*dx*dy+Q[1,1]*dy*dy;qmin=max(0.,c-b*b/a)
    if qmin>=1:continue
    center=-b/a;half=np.sqrt((1-qmin)/a);hi=center+half;lo=max(center-half,depth[y,x]-p[j,2])
    if lo>=hi:continue
    front[y,x]+=norm*integral(lo,hi,center,a,qmin,nodes,weights)
 return front

def foreground_columns(p,cov,mass,depth,W,H,pixel,cy,order=16):
 C=cov.astype('f8').copy();C[:,0,0]+=pixel*pixel/12;C[:,1,1]+=pixel*pixel/12;nodes,weights=np.polynomial.legendre.leggauss(order)
 return columns_core(p,C,np.linalg.inv(C),np.linalg.det(C),mass,depth,W,H,pixel,cy,nodes,weights)
