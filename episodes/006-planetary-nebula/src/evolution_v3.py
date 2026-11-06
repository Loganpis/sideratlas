"""Same MIST star; continuous secular clock without the old3yr camera hold."""
import numpy as np
from scipy.interpolate import PchipInterpolator
from .evolution_v2 import Track as PriorTrack
class Track(PriorTrack):
    def __init__(self):
        super().__init__()
        self.clock_s=np.array([0,4,6,7.2,8,10,12,14,16,18,20,24,28,32.])
        # End on the mature nebula, before the earlier45kyr fading-halo endpoint.
        # The object is a hot post-AGB core, not yet No003's cold ZZCeti star.
        self.clock_age=np.array([4.57e9,10.9e9,self.age[604],11.44e9,self.age[807],self.t0-22000,self.t0-4000,self.t0-1000,self.t0+100,self.t0+600,self.t0+1200,self.t0+2500,self.t0+8000,self.t0+16000])
        self.clock=PchipInterpolator(self.clock_s,self.clock_age-self.t0)
