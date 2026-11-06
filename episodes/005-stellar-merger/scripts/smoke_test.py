"""Small rendering/accounting test; never launches SPH or a new gas integration."""
from pathlib import Path
import sys,json,numpy as np,argparse,hashlib
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.verify_inputs import verify
from src.rendering.final_scene import create
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'output/smoke-test.json');a=p.parse_args();verify();r=create(144);rows=[]
    for t in [10.5,18.5,28.5]:
        im=r.frame(t,30,4);q=r.records[-1];m=q['microflow']
        tests=dict(finite_native_computation=im.shape==(256,144,3) and np.isfinite(im).all(),mass=q['mass_identity_error']<1e-10,cooling=abs(m['source_power_W']-m['expected_fine_power_W'])<=1e-12*max(1,m['expected_fine_power_W']),centroids=m['max_weighted_parent_centroid_error_Rsun']<1e-12,material_only=q['all_gas_emission_particle_sampled'] and q['smooth_cooling_layer_max']==0,no_new_light=q['displayed_reprocessed_luminosity_W']==0,opacity=q['opacity_profile']=='wendland_c2' and q['density_reconstruction']['density_scale']==.4,selected_camera=abs(q['camera_halfwidth_Rsun']-55)<1e-10 if t>24 else True)
        rows.append(dict(playback_s=t,tests={k:bool(v) for k,v in tests.items()},rgb16_sha256=hashlib.sha256(im.tobytes()).hexdigest()))
    result=dict(passed=all(all(x['tests'].values()) for x in rows),width=144,height=256,frames=rows,full_reel_rerendered=False,not_hydrodynamic_validation=True)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');assert result['passed'];print(json.dumps(result),flush=True)
if __name__=='__main__':main()
