"""Single supported release scene: the selected 55-Rsun production camera."""
from pathlib import Path
import json
from src.detail.polar_final import PolarFinalRenderer
ROOT=Path(__file__).resolve().parents[2]
def create(width=1080):
    p=ROOT/'config/final-release.json';c=json.loads(p.read_text())
    assert c['material_only_emission'] and c['emitters_per_parent']==512
    assert c['density_profile']=='wendland_c2' and c['density_scale']==.4
    assert c['wendland_pixel_ray_factor']==4 and c['camera_inclination_offset_degrees']==-10
    return PolarFinalRenderer(ROOT/'cache/motion',width,p)
