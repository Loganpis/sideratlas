"""Selected narrow repair; the tested AGB extension is explicitly excluded."""
from .production_temporal_fix import TemporalSamplingRenderer,PlaybackCoordinate,ResolvedEvolutionClock
from .production_camera_fix import CameraContinuityRenderer

class PublicationRenderer(TemporalSamplingRenderer):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.sampling_plan=dict(self.sampling_plan,windows=[w for w in self.sampling_plan['windows'] if w['end_s']<8])
        self.coordinate=PlaybackCoordinate(self.sampling_plan)
        self.track.clock=ResolvedEvolutionClock(self.original_clock,self.coordinate)

    def geometry_v3(self,t,star,state):
        return CameraContinuityRenderer.geometry_v3(self,t,star,state)
