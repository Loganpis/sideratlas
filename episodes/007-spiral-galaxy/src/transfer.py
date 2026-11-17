"""Finite-height extinction, line footprints and conservative light projection."""
import numpy as np
from numba import njit,prange
GAS_COLUMN_TO_G_CM2=1.98847e33/(3.085677581491367e21)**2
NODES,WEIGHTS=np.polynomial.legendre.leggauss(24)

def ccm_ratio(wavelength_nm):
    """CCM89 optical polynomial, RV=3.1; extinction, not absorption alone."""
    y=1000/np.asarray(wavelength_nm)-1.82
    a=1+.17699*y-.50447*y**2-.02427*y**3+.72085*y**4+.01979*y**5-.77530*y**6+.32999*y**7
    b=1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-.62251*y**5+5.30260*y**6-2.09002*y**7
    return a+b/3.1

EXT=np.array(ccm_ratio(np.array([650.,550.,450.])),dtype=float)
LINE_EXT=np.array(ccm_ratio(np.array([656.3,486.1])),dtype=float)

@njit(cache=True)
def sample_grid(field,x,y,dx):
    n=field.shape[0]; u=(x+.5*n*dx)/dx-.5; v=(y+.5*n*dx)/dx-.5
    i=int(np.floor(u)); j=int(np.floor(v)); a=u-i; b=v-j
    result=0.
    for jj in range(2):
        for ii in range(2):
            if 0<=j+jj<n and 0<=i+ii<n:
                result+=field[j+jj,i+ii]*(a if ii else 1-a)*(b if jj else 1-b)
    return result

@njit(cache=True)
def dust_column(points,sigma,dx,nview):
    """Finite-height foreground column; no transparent-star override."""
    out=np.zeros(len(points)); h=.08; zmax=.6
    for i in range(len(points)):
        z0=max(-zmax,points[i,2]); length=zmax-z0
        if length<=0: continue
        for k in range(len(NODES)):
            z=z0+.5*length*(NODES[k]+1)
            s=(z-points[i,2])/nview[2]
            x=points[i,0]+s*nview[0]; y=points[i,1]+s*nview[1]
            rho=sample_grid(sigma,x,y,dx)*np.exp(-.5*(z/h)**2)/(np.sqrt(2*np.pi)*h)
            out[i]+=.5*length*WEIGHTS[k]*rho/nview[2]
    return out

@njit(cache=True)
def cic_rgb(xy,lum,h,w):
    out=np.zeros((h,w,3),np.float64)
    for k in range(len(xy)):
        x=xy[k,0]-.5; y=xy[k,1]-.5
        i=int(np.floor(x)); j=int(np.floor(y)); u=x-i; v=y-j
        for jj in range(2):
            for ii in range(2):
                if 0<=j+jj<h and 0<=i+ii<w:
                    weight=(u if ii else 1-u)*(v if jj else 1-v)
                    for c in range(3): out[j+jj,i+ii,c]+=weight*lum[k,c]
    return out

@njit(cache=True)
def compact_sources(xy,lum,sigma,h,w):
    """Flux-normalized subpixel Gaussian approximation to unresolved HII
    spheres (same second moment), before a separate optical PSF.
    """
    out=np.zeros((h,w,3),np.float64)
    for k in range(len(xy)):
        if lum[k,0]+lum[k,1]+lum[k,2]<1e-30: continue
        s=max(.42,sigma[k]); r=int(np.ceil(3.5*s))
        x=xy[k,0]-.5; y=xy[k,1]-.5; i=int(np.floor(x)); j=int(np.floor(y))
        norm=0.
        for jj in range(j-r,j+r+2):
            for ii in range(i-r,i+r+2):
                norm+=np.exp(-.5*((ii-x)**2+(jj-y)**2)/s**2)
        for jj in range(max(0,j-r),min(h,j+r+2)):
            for ii in range(max(0,i-r),min(w,i+r+2)):
                a=np.exp(-.5*((ii-x)**2+(jj-y)**2)/s**2)/norm
                for c in range(3): out[jj,ii,c]+=a*lum[k,c]
    return out

def bulge_column(R):
    x=np.maximum(R/.55,1e-5); b=np.empty_like(x)
    lo=x<.997; hi=x>1.003; mid=~(lo|hi)
    z=x[lo]; X=np.arccosh(1/z)/np.sqrt(1-z*z)
    b[lo]=((2+z*z)*X-3)/(1-z*z)**2
    z=x[hi]; X=np.arccos(1/z)/np.sqrt(z*z-1)
    b[hi]=((2+z*z)*X-3)/(1-z*z)**2
    d=x[mid]-1
    b[mid]=4/15-16/35*d+8/15*d*d
    return (1.5e9/3)/(2*np.pi*.55**2)*b

@njit(cache=True)
def columns_at_points(points,sigma,dx):
    out=np.empty(len(points))
    for i in range(len(points)):out[i]=sample_grid(sigma,points[i,0],points[i,1],dx)
    return out

@njit(cache=True, parallel=True)
def slab_transfer_parallel(old,mid,aged,gas,x0,y0,view,dx,vertical,color0,color1,B,bulge_color,zmax):
    """Formal homogeneous-slab transfer per cell, behind -> observer.

    Vertical emissivity PDFs are measured from current orbital populations.
    The aged association continuum is in its calculated birth/orbit midplane.
    Bulge rear/front halves use the explicitly declared thin-screen closure.
    """
    h,w=x0.shape;nz=vertical.shape[1];dz=2*zmax/nz;ds=dz/view[2]
    out=np.zeros((4,h,w,3));tauV=np.zeros((h,w))
    for j in prange(h):
        for i in range(w):
            intensity=np.zeros((4,3))
            for c in range(3):intensity[3,c]=.5*B[j,i]*bulge_color[c]
            column=0.
            for k in range(nz):
                z=-zmax+(k+.5)*dz
                x=x0[j,i]+z*view[0]/view[2];y=y0[j,i]+z*view[1]/view[2]
                r=np.hypot(x,y);u=min(r/.1,vertical.shape[2]-1.000001);ir=int(u);f=u-ir
                p0=(1-f)*vertical[0,k,ir]+f*vertical[0,k,ir+1]
                p1=(1-f)*vertical[1,k,ir]+f*vertical[1,k,ir+1]
                m0=sample_grid(old,x,y,dx)*p0*ds;m1=sample_grid(mid,x,y,dx)*p1*ds
                rho=sample_grid(gas,x,y,dx)*np.exp(-.5*(z/.08)**2)/(np.sqrt(2*np.pi)*.08)
                dtau=rho*ds*GAS_COLUMN_TO_G_CM2*220.;column+=dtau
                for c in range(3):
                    a=dtau*EXT[c]
                    if a<1e-5:
                        tr=1-a+.5*a*a;escape=1-.5*a+a*a/6
                    else:
                        tr=np.exp(-a);escape=-np.expm1(-a)/a
                    intensity[0,c]=intensity[0,c]*tr+m0*color0[c]*escape
                    intensity[1,c]=intensity[1,c]*tr+m1*color1[c]*escape
                    intensity[2,c]*=tr;intensity[3,c]*=tr
                if k==nz//2-1:
                    for c in range(3):intensity[2,c]+=sample_grid(aged[c],x0[j,i],y0[j,i],dx)/view[2]
            for c in range(3):intensity[3,c]+=.5*B[j,i]*bulge_color[c]
            for q in range(4):
                for c in range(3):out[q,j,i,c]=intensity[q,c]
            tauV[j,i]=column
    return out,tauV
