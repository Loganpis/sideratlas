"""Production timing repair; no stellar data, gas physics or image blending."""
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicHermiteSpline,PchipInterpolator
from .production_camera_fix import CameraContinuityRenderer
from .evolution import smooth
ROOT=Path(__file__).resolve().parents[1]

class PlaybackCoordinate:
    def __init__(self,plan):
        self.windows=[(w['start_s'],w['end_s'],CubicHermiteSpline(w['display_knots_s'],w['original_playback_knots_s'],w['derivatives'])) for w in plan['windows']]
    def __call__(self,t):
        t=np.asarray(t,dtype=float);result=t.copy()
        for lo,hi,curve in self.windows:
            result=np.where((t>lo)&(t<hi),curve(t),result)
        return result

class ResolvedEvolutionClock:
    def __init__(self,original,coordinate):self.original=original;self.coordinate=coordinate
    def __call__(self,t):return self.original(self.coordinate(t))

class TemporalSamplingRenderer(CameraContinuityRenderer):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.original_clock=self.track.clock
        self.sampling_plan=json.loads((ROOT/'production/TEMPORAL-SAMPLING.json').read_text())
        self.coordinate=PlaybackCoordinate(self.sampling_plan)
        self.camera_reference=PchipInterpolator(self.sampling_plan['camera_reference_times_s'],self.sampling_plan['camera_reference_log_halfwidth_pc'])
        self.track.clock=ResolvedEvolutionClock(self.original_clock,self.coordinate)
        self.actual_surface_phases=[];original_quiet=self.solar._quiet
        def monitored_quiet(*args,**kwargs):
            phase=float(args[2] if len(args)>2 else kwargs['phase']);self.actual_surface_phases.append(phase)
            return original_quiet(*args,**kwargs)
        self.solar._quiet=monitored_quiet

    def geometry_v3(self,t,star,state):
        original=super().geometry_v3(t,star,state)
        if hasattr(self,'coordinate') and 8<t<12:
            target=float(self.camera_reference(float(self.coordinate(t))))
            edge=self.sampling_plan['camera_endpoint_blend_s']
            blend=float(smooth((t-8)/edge))*float(smooth((12-t)/edge))
            return float(np.exp((1-blend)*np.log(original)+blend*target))
        return original

    def production_frame(self,index):
        self.actual_surface_phases=[]
        image,ledger=super().production_frame(index)
        times=np.array([s['playback_s'] for s in ledger['samples']])
        for phase in self.actual_surface_phases:assert np.min(np.abs(times-phase))<1e-12
        for sample in ledger['samples']:
            t=sample['playback_s'];sample['original_playback_coordinate_s']=float(self.coordinate(t))
            sample['pulsation_phase_radians']=float(2*np.pi*(t-8)/1.8)
            sample['camera_azimuth_deg']=float(-15+np.rad2deg(1.1)*t/32)
        ledger['temporal_sampling_revision']='production/TEMPORAL-SAMPLING.json'
        ledger['original_playback_coordinate_s']=float(self.coordinate(ledger['playback_s']))
        ledger['surface_clock_observed_s']=self.actual_surface_phases.copy()
        ledger['source_clock_rule']='original No004 source_start_s+20*global playback_s; no reset or retiming'
        return image,ledger
