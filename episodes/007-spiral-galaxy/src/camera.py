"""Accepted orthographic inspection camera, independent of physical time."""
import numpy as np
from scipy.interpolate import PchipInterpolator

class Camera3:
    """Right/up/view vectors rotate the 3D rays, not a rendered image."""
    def __init__(self,composition='A'):
        s=[0,8,16,24,32]
        if composition!='A':raise ValueError('Only accepted camera A is distributed')
        di,dw=0.,0.
        self.width=PchipInterpolator(s,np.array([18.8,18.4,18.0,17.2,18.8])+dw)
        self.inc=PchipInterpolator(s,np.deg2rad(np.array([53,50,48,50,54])+di))
        self.az=PchipInterpolator(s,[1.34,1.37,1.40,1.43,1.46])
        self.roll=PchipInterpolator(s,np.deg2rad([-59,-60,-61,-60.5,-59.5]))
        self.composition=composition

    def at(self,s):
        inc=float(self.inc(s));az=float(self.az(s));roll=float(self.roll(s))
        n=np.array([np.sin(inc)*np.cos(az),np.sin(inc)*np.sin(az),np.cos(inc)])
        basic=np.array([-np.sin(az),np.cos(az),0.]);up=np.cross(n,basic)
        return float(self.width(s)),np.cos(roll)*basic+np.sin(roll)*up,-np.sin(roll)*basic+np.cos(roll)*up,n
