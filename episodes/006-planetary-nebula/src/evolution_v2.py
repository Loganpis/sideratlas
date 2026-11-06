"""One MIST progenitor; explicit small calibration to the No004 initial reference."""
import numpy as np
from scipy.interpolate import PchipInterpolator
from .evolution import Track as OriginalTrack,LSUN,PC,SIGMA,smooth
class Track(OriginalTrack):
    def __init__(self):
        super().__init__()
        initial=super().state(4.57e9-self.t0)
        self.initial_L=initial['L'];self.initial_T=initial['T']
        self.clock_s=np.array([0,4,6,7.2,8,10,11,13,16,20,24,28,32.])
        self.clock_age=np.array([4.57e9,10.9e9,self.age[604],11.44e9,self.age[807],self.t0-4000,self.t0-3000,self.t0-2997,self.t0,self.t0+2500,self.t0+8000,self.t0+16000,self.t0+45000])
        self.clock=PchipInterpolator(self.clock_s,self.clock_age-self.t0)
    def state(self,t):
        s=super().state(t)
        # A declared reference-match, NOT recomputed solar-calibrated evolution.
        # The published 1Msun mass/composition history is unchanged. The tiny
        # L/T calibration smoothly vanishes within the first evolutionary Gyr.
        if hasattr(self,'initial_L'):
            w=1-float(smooth((t+self.t0-4.57e9)/1e9))
            s['L']*= (LSUN/self.initial_L)**w
            s['T']*= (5772/self.initial_T)**w
            s['R']=np.sqrt(s['L']/(4*np.pi*SIGMA*s['T']**4))/PC
            s['reference_calibration_weight']=w
        return s
