"""Published MESA/MIST 1 Msun track. No new stellar evolution is solved."""
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1]
MSUN=1.98847e33; RSUN=6.957e10; LSUN=3.828e33; PC=3.0856775814913673e18
YR=31557600.; H=6.62607015e-27; C=2.99792458e10; KB=1.380649e-16; MH=1.6735575e-24
SIGMA=5.670374419e-5; KM_S=1e5*YR/PC

def smooth(x):
    x=np.clip(x,0,1);return x*x*(3-2*x)

class Track:
    def __init__(self):
        p=ROOT/'data/evolution/00100M.track.eep'
        names=next(l[1:].split() for l in p.read_text().splitlines() if l.startswith('#') and 'star_age' in l)
        a=np.loadtxt(p);self.d={k:a[:,i] for i,k in enumerate(names)}
        self.age=self.d['star_age']; self.t0=self.age[np.flatnonzero((self.d['phase']==6)&(self.d['log_Teff']>3.85))[0]]
        self.rel=self.age-self.t0
        self.clock_s=np.array([0,4,6,7.2,8,10,11,13,16,20,24,32.])
        self.clock_age=np.array([4.57e9,10.9e9,self.age[604],11.44e9,self.age[807],self.t0-4000,self.t0-3000,self.t0-2997,self.t0,self.t0+2500,self.t0+8000,self.t0+16000])
        self.clock=PchipInterpolator(self.clock_s,self.clock_age-self.t0)
        self.massfn=PchipInterpolator(self.rel,self.d['star_mass'])
        self.M0=float(self.massfn(0)); self.mdotfn=self.massfn.derivative()
        self.X=float(self.d['surface_h1'][1454]); self.O_H=float(self.d['surface_o16'][1454]/16/self.X);self.N_H=float(self.d['surface_n14'][1454]/14/self.X)
        self.cmf=np.loadtxt(ROOT/'data/spectral/CIE_xyz_1931_2deg.csv',delimiter=',')
        self.xyzrgb=np.array([[3.2404542,-1.5371385,-.4985314],[-.969266,1.8760108,.041556],[.0556434,-.2040259,1.0572252]])
        self.Ts=np.geomspace(2000,200000,400);w=self.cmf[:,0]*1e-7
        B=2*H*C*C/w[None,:]**5/np.expm1(H*C/(w[None,:]*KB*self.Ts[:,None]))
        xyz=np.trapz(B[:,:,None]*self.cmf[None,:,1:4],w,axis=1)
        self.rgbs=np.maximum(0,xyz@self.xyzrgb.T)
        self.Yref=float(np.interp(5772,self.Ts,xyz[:,1]))
    def state(self,t):
        vals={k:float(np.interp(t,self.rel,self.d[k])) for k in ['star_mass','log_L','log_Teff','log_R','phase']}
        vals.update(M=vals['star_mass'],L=10**vals['log_L']*LSUN,T=10**vals['log_Teff'],R=10**vals['log_R']*RSUN/PC)
        return vals
    def rgb(self,T):
        return np.stack([np.interp(np.asarray(T),self.Ts,self.rgbs[:,i]) for i in range(3)],axis=-1)/self.Yref
    def ionizing(self,t,threshold=13.6):
        s=self.state(t);x=threshold*1.602176634e-12/(KB*s['T'])
        f=quad(lambda z:z*z*np.exp(-z)/(-np.expm1(-z)),x,max(100,x+70),epsabs=1e-13)[0]
        q=15*s['L']/np.pi**4/(KB*s['T'])*f
        en=quad(lambda z:z**3*np.exp(-z)/(-np.expm1(-z)),x,max(100,x+70),epsabs=1e-13)[0]
        return q,15*s['L']/np.pi**4*en
