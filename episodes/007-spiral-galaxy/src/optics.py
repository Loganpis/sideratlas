"""Accepted Bradford/D65 display adaptation and bounded radial contrast."""
import numpy as np
from scipy.ndimage import gaussian_filter1d
from .populations import YROW
from .spectrum import XYZ_RGB
def white_balance_matrix(solar_rgb):
    # ICC Bradford cone-space transform. Choosing a5772K photographic white
    # reference is a declared display convention, not altered source spectra.
    B=np.array([[.8951,.2664,-.1614],[-.7502,1.7135,.0367],[.0389,-.0685,1.0296]])
    RGB_XYZ=np.linalg.inv(XYZ_RGB);source=RGB_XYZ@solar_rgb;source/=source[1]
    target=np.array([.95047,1.,1.08883]);M=np.linalg.inv(B)@np.diag((B@target)/(B@source))@B
    return XYZ_RGB@M@RGB_XYZ

def radial_tone(linear,pixel_kpc,matrix,mode,exposure=1.3e8,projected_radius=None):
    if mode not in ['radial','radial-balance']:raise ValueError('Only accepted radial branch is qualified')
    physical_Y=np.maximum(linear@YROW,1e-90);rgb=linear
    if mode in ['balance','combined','radial-balance']:
        rgb=np.maximum(linear@matrix.T,0)
        # Keep the computed photopic luminance reference. Only photographic
        # chromatic appearance changes, with no assigned young-star gain.
        rgb*=np.divide(physical_Y,np.maximum(rgb@YROW,1e-90))[:,:,None]
    if mode in ['radial','radial-balance']:
        if projected_radius is None:raise ValueError('Radial display requires the actual frozen camera coordinates')
        # A local display reference from observed linear light, not an arm
        # template. Preserve faint/dark light; boost only measured positive
        # departures from the smoothed annular mean. No source is reweighted.
        bins=np.minimum((projected_radius/.1).astype(int),499)
        weight=np.bincount(bins.ravel(),minlength=500)
        total=np.bincount(bins.ravel(),weights=physical_Y.ravel(),minlength=500)
        mean=np.divide(gaussian_filter1d(total,4.,mode='nearest'),np.maximum(gaussian_filter1d(weight.astype(float),4.,mode='nearest'),1e-90))
        reference=np.interp(projected_radius,(np.arange(500)+.5)*.1,mean)
        excess=np.maximum(np.log((physical_Y+1e5)/(reference+1e5)),0)
        gain=1+1.0*np.tanh(excess/.5)*physical_Y/(physical_Y+3e6)
    elif mode in ['positive','combined']:
        gain=1+.65*np.tanh(np.maximum(ratio,0)/.8)*physical_Y/(physical_Y+3e6)
    else:gain=np.exp(.25*np.tanh(ratio/.8)*physical_Y/(physical_Y+3e6))
    Y=physical_Y*gain;mapped=1-np.exp(-.8*np.arcsinh((Y/exposure)**1.12))
    rgb=rgb*(mapped/physical_Y)[:,:,None]
    chroma=1/(1+.08*np.maximum(Y/exposure-3,0))
    rgb=mapped[:,:,None]+chroma[:,:,None]*(rgb-mapped[:,:,None])
    rgb/=np.maximum(1,np.max(rgb,axis=2))[:,:,None];rgb=np.clip(rgb,0,1)
    return np.where(rgb<=.0031308,12.92*rgb,1.055*rgb**(1/2.4)-.055)
