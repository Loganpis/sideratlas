"""Dimensionless No.004 morphology surrogate; no solar thermodynamics copied."""
from pathlib import Path
import json,sys,time,hashlib
import numpy as np
from scipy.ndimage import gaussian_filter
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
SOURCE=None
Y=np.array([.2126729,.7151522,.0721750])
def main():
    global SOURCE
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);a=p.parse_args();SOURCE=a.source
    start=time.monotonic();out=ROOT/'data/quality-surface';out.mkdir(parents=True,exist_ok=True)
    target=out/'morphology.npz'
    if target.exists():raise RuntimeError('Reuse existing derived source rather than overwrite')
    files=sorted(p for p in (SOURCE/'data/continuum').glob('*-17-2-solar.npz') if p.with_suffix('.json').exists())
    maps=[];times=[];manifest=[]
    for p in files:
        d=np.load(p);intensity=np.maximum(d['rgb'][[0,2,4]]@Y,1e-20)
        z=np.log(intensity);z-=z.mean(axis=(-2,-1),keepdims=True);z/=np.maximum(z.std(axis=(-2,-1),keepdims=True),1e-8)
        # Dimensionless topology/contrast only. Mean spectrum, temperature,
        # gravity, fields and absolute solar intensity are discarded.
        maps.append(np.int16(np.clip(z,-7.9,7.9)*4096));times.append(float(d['time_s']))
        manifest.append(dict(path=str(p.relative_to(SOURCE)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    np.savez_compressed(target,z=np.stack(maps),time=np.array(times),mu=np.array([.1,.55,1.]),scale=np.array(4096.))
    active=np.load(SOURCE/'data/active-region/intensity-sequence.npz');t=active['time_s'];a=active['intensity_cgs'][::10]
    quiet=float(np.quantile(a[0],.70));ratio=np.clip(a/quiet,.05,1.5)
    # Dark-region morphology remains a separately labelled illustrative solar
    # template, not computed giant magnetic geometry or calibrated spot Teff.
    np.savez_compressed(out/'illustrative-active-region.npz',ratio=ratio.astype('f4'),time=t[::10],quiet=quiet)
    # Large-scale surface modulation is derived by coarse filtering the same
    # radiation sequence, rather than independent random image texture.
    center=np.stack([m[-1,0].astype('f4')/4096 for m in maps]);large=gaussian_filter(center,(0,10,10),mode=('nearest','wrap','wrap'))
    large/=np.maximum(large.std(axis=(1,2),keepdims=True),1e-6)
    np.savez_compressed(out/'mesoscales.npz',z=large.astype('f4'),time=np.array(times))
    meta=dict(scope=__doc__,morphology_source='SiderAtlas No.004 calculated Bifrost continuum maps',sources=manifest,source_time_range=[times[0],times[-1]],source_mu=[.1,.55,1.],solar_physical_state_reused=False,statistical_giant_morphology_surrogate=True,illustrative_active_region='Cheung et al. 2018/2019 MURaM outgoing intensity, CC BY 3.0, rethermalized morphology only',magnetic_field_calculated_for_donor=False,source_solar_reference_turnover_s=300,donor_granule_diameter_Hp=15,donor_display_magnification=6,donor_assumed_bolometric_RMS=.18,granulation_display_contrast_exponent=3.25,convection_clock_slowdown=32,donor_turnover_s=2548.2776557705824*15/8,elapsed_s=time.monotonic()-start)
    (out/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps({k:v for k,v in meta.items() if k!='sources'},indent=2))
if __name__=='__main__':main()
