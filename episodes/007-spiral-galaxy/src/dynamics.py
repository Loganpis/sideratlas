"""Prescribed background gravity, KDK orbits, and original initial sampling."""
import numpy as np
from numba import njit
from . import CONFIG
G=4.498502151575286e-12
KM_S=0.001022712165045695
VREF2=(190*KM_S)**2
BASE=np.array(CONFIG['model']['base_parameters'])
EMPTY=np.empty((0,9))
WAVES=EMPTY

@njit(cache=True)
def envelope(t, w):
    u = (t-w[5])/(w[6]-w[5])
    if u <= 0 or u >= 1:
        return 0., 0.
    return (w[4]*VREF2*np.sin(np.pi*u)**2,
            w[4]*VREF2*np.pi/(w[6]-w[5])*np.sin(2*np.pi*u))

@njit(cache=True)
def radial_basis(r, scale):
    u = (r/1.5)**4
    inner = -np.expm1(-u)
    outer = np.exp((4.-r)/scale-(r/17.)**6)
    f = inner*outer
    df = (4*r**3/1.5**4*np.exp(-u)
          +inner*(-1/scale-6*r**5/17.**6))*outer
    return f, df

@njit(cache=True)
def field(x, y, z, t, base=BASE, waves=WAVES):
    r2 = x*x+y*y+z*z
    r = max(np.sqrt(r2), 1e-8)
    R = max(np.sqrt(x*x+y*y), 1e-8)
    ax=0.; ay=0.; az=0.; phi=0.
    for off in (0,3):
        M,a,b = base[off],base[off+1],base[off+2]
        zz=np.sqrt(z*z+b*b); B=a+zz
        d=np.sqrt(R*R+B*B); fac=-G*M/d**3
        ax+=fac*x; ay+=fac*y; az+=fac*B*z/zz
        phi-=G*M/d
    Mb,ab,Mh,rs,fc = base[6],base[7],base[8],base[9],base[10]
    fac=-G*Mb/(r*(r+ab)**2)
    ax+=fac*x; ay+=fac*y; az+=fac*z; phi-=G*Mb/(r+ab)
    u=r/rs; fm=np.log1p(u)-u/(1+u)
    if u<1e-4:
        fm=u*u*(.5-2*u/3+3*u*u/4)
    fac=-G*Mh*fm/(fc*r**3)
    ax+=fac*x; ay+=fac*y; az+=fac*z
    phi-=G*Mh*np.log1p(u)/(fc*r)
    power=0.; torque=0.
    if R>.02:
        theta=np.arctan2(y,x)
        for j in range(len(waves)):
            w=waves[j]; A,Adot=envelope(t,w)
            if A==0 and Adot==0:
                continue
            f,df=radial_basis(R,w[8]); vertical=np.exp(-.5*(z/w[7])**2)
            phase=w[0]*(theta-w[2]*t)-w[1]*np.log(R/4.)+w[3]
            c=np.cos(phase); s=np.sin(phase)
            ar=A*vertical*(df*c+f*w[1]/R*s)
            at=-A*vertical*f*w[0]/R*s
            ax+=ar*x/R-at*y/R; ay+=ar*y/R+at*x/R
            az-=A*f*vertical*z/w[7]**2*c
            phi-=A*f*vertical*c
            tj=-A*f*vertical*w[0]*s
            torque+=tj; power+=w[2]*tj-Adot*f*vertical*c
    return ax,ay,az,phi,power,torque

@njit(cache=True)
def evaluate(x, t, base=BASE, waves=WAVES, scalar=field):
    out=np.empty((len(x),6))
    for i in range(len(x)):
        out[i]=scalar(x[i,0],x[i,1],x[i,2],t,base,waves)
    return out

@njit(cache=True)
def orbit_step(x,v,a,dt,t,base=BASE,waves=WAVES,scalar=field,weights=None):
    """KDK + trapezoidal explicit-potential work and external torque."""
    work=0.; impulse=0.
    for i in range(len(x)):
        for k in range(3):
            v[i,k]+=.5*dt*a[i,k]
            x[i,k]+=dt*v[i,k]
        b=scalar(x[i,0],x[i,1],x[i,2],t+dt,base,waves)
        weight=1. if weights is None else weights[i]
        work+=weight*.5*dt*(a[i,4]+b[4])
        impulse+=weight*.5*dt*(a[i,5]+b[5])
        for k in range(3):
            v[i,k]+=.5*dt*b[k]
        for k in range(6):
            a[i,k]=b[k]
    return work,impulse

def circular_speed(R):
    R=np.atleast_1d(R).astype(float)
    a=evaluate(np.c_[R,R*0,R*0],0.,BASE,WAVES[:0])
    return np.sqrt(np.maximum(0.,-R*a[:,0]))

def initial_stars(n, seed=7007,scale_length=3.,height0=.13,height_slope=.006,sr_base=10.,sr_extra=14.):
    """Exponential LIGHT tracers with approximate Jeans/epicyclic support.

    Not an exact DF for the MN mass source; the independent control measures
    the resulting adjustment rather than asserting exact equilibrium.
    """
    rng=np.random.default_rng(seed)
    R=rng.gamma(2.,scale_length,n)
    bad=(R<.45)|(R>17)
    while np.any(bad):
        R[bad]=rng.gamma(2.,scale_length,bad.sum()); bad=(R<.45)|(R>17)
    theta=rng.uniform(0,2*np.pi,n)
    h=height0+height_slope*R
    z=rng.normal(0,h)
    vc=circular_speed(R)
    eps=.001
    vp=circular_speed(R*(1+eps)); vm=circular_speed(R*(1-eps))
    dlnv=np.log(vp/vm)/(2*eps)
    ratio=np.clip(.5*(1+dlnv),.25,1.)
    sr=KM_S*(sr_base+sr_extra*np.exp(-R/8))
    deriv=-R/scale_length-(2*sr_extra*R/8*np.exp(-R/8))/(sr_base+sr_extra*np.exp(-R/8))
    vmean=np.sqrt(np.maximum(vc*vc+sr*sr*(1-ratio+deriv),.1*vc*vc))
    vr=rng.normal(0,sr); vt=vmean+rng.normal(0,sr*np.sqrt(ratio))
    probe=np.c_[R,R*0,np.full(n,.001)]
    nu2=-evaluate(probe,0.,BASE,WAVES[:0])[:,2]/.001
    vz=rng.normal(0,h*np.sqrt(np.maximum(nu2,0)))
    x=np.c_[R*np.cos(theta),R*np.sin(theta),z]
    v=np.c_[vr*np.cos(theta)-vt*np.sin(theta),
            vr*np.sin(theta)+vt*np.cos(theta),vz]
    return x,v,rng

def initial_cohorts(n=200000):
    n0=round(.6*n);n1=n-n0
    x0,v0,_=initial_stars(n0,seed=7007,scale_length=2.7,height0=.20,height_slope=.004,sr_base=16,sr_extra=20)
    x1,v1,_=initial_stars(n1,seed=9007,scale_length=3.8,height0=.10,height_slope=.004,sr_base=9,sr_extra=11)
    x=np.r_[x0,x1];v=np.r_[v0,v1]
    cohort=np.r_[np.zeros(n0,dtype='i1'),np.ones(n1,dtype='i1')]
    mass=np.r_[np.full(n0,BASE[0]*.6/n0),np.full(n1,BASE[0]*.4/n1)]
    return x,v,cohort,mass
