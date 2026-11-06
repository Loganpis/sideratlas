"""No.004-style Airy + Moffat imaging with positive pixel-integrated kernels.

The inherited frequency-cutoff PSF at subpixel widths can have negative sampled
lobes around an unresolved source. Direct pixel integration avoids that artifact.
No halo is treated as gas or additional physical stellar luminosity.
"""
import numpy as np
from scipy.special import j1
from scipy.signal import fftconvolve
class Imaging:
    def __init__(self,shape):
        h,w=shape;scale=w/1080;zero=1.3*scale
        y,x=np.mgrid[-8:9,-8:9];kernel=np.zeros(x.shape,float)
        for dx in (np.arange(12)+.5)/12-.5:
            for dy in (np.arange(12)+.5)/12-.5:
                q=3.83170597*np.hypot(x+dx,y+dy)/zero
                a=np.ones_like(q);np.divide(2*j1(q),q,out=a,where=q>0)
                kernel+=a*a
        self.core=kernel/kernel.sum()
        radius=round(230*scale);y,x=np.mgrid[-radius:radius+1,-radius:radius+1]
        self.halo=(1+(x*x+y*y)/(75*scale)**2)**-2.5;self.halo/=self.halo.sum()
        self.flux_errors=[];self.crop_fractions=[]
    def process(self,rgb):
        core=np.stack([fftconvolve(rgb[...,c],self.core,mode='same') for c in range(3)],-1)
        halo=np.stack([fftconvolve(rgb[...,c],self.halo,mode='same') for c in range(3)],-1)
        hdr=np.maximum(.96*core+.04*1.8*halo,0)
        y=hdr@np.array([.2126729,.7151522,.0721750]);out=hdr*np.divide(-np.expm1(-y),y,out=np.ones_like(y),where=y>1e-30)[...,None]
        return out/np.maximum(1.,out.max(axis=-1))[...,None]
