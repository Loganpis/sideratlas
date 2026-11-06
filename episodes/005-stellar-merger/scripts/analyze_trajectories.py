"""Inspect public cached dynamics without the two external surface maps."""
from pathlib import Path
import sys,json,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.motion_view import MotionData
def main():
    m=MotionData(ROOT/'cache/motion');m.enable_force_quintic(ROOT/'cache/motion/playback-accelerations.npz');rows=[]
    for playback in [0,8,16,24,32]:
        q=m.at(float(m.clock(playback)));live=q['live'];counts={name:int(np.sum(live&(q['classes']==code))) for name,code in [('outward_energy_candidates',0),('other_instantaneously_bound',1),('returning',2),('other',3)]}
        rows.append(dict(playback_s=playback,model_time=q['t'],live_tracers=int(live.sum()),gas_mass_msun=float(live.sum()*m.mass),stellar_radius_Rsun=float(q['radius']),classes=counts))
    print(json.dumps(dict(scope='Cached gravity-guided dynamics, not hydro; instantaneous energy is not guaranteed escape',surface_maps_required=False,states=rows),indent=2))
if __name__=='__main__':main()
