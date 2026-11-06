"""Local camera-continuity repair; accepted physics and imaging stay inherited."""
import json
from pathlib import Path
import numpy as np
from .production_render import ProductionRenderer
from .evolution import smooth

ROOT = Path(__file__).resolve().parents[1]

def camera_halfwidth(t, original, plan):
    for w in plan['windows']:
        lo, hi = w['start_s'], w['end_s']
        if lo < t < hi:
            f = (t - lo) / (hi - lo)
            target = (1-f)*np.log(w['halfwidth_start_pc']) + f*np.log(w['halfwidth_end_pc'])
            blend = float(smooth((t-lo)/w['edge_blend_s'])) * float(smooth((hi-t)/w['edge_blend_s']))
            return float(np.exp((1-blend)*np.log(original) + blend*target))
    return original

class CameraContinuityRenderer(ProductionRenderer):
    def __init__(self, *args, **kwargs):
        self.camera_plan = json.loads((ROOT/'production/CAMERA-CONTINUITY.json').read_text())
        super().__init__(*args, **kwargs)

    def geometry_v3(self, t, star, state):
        original = super().geometry_v3(t, star, state)
        return camera_halfwidth(t, original, self.camera_plan)

    def production_frame(self, index):
        image, ledger = super().production_frame(index)
        ledger['camera_continuity_revision'] = 'production/CAMERA-CONTINUITY.json'
        return image, ledger
