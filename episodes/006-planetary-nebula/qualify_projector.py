"""Bounded CPU/Metal comparison on the selected 27-second nebular state.

This qualifies projection arithmetic, not the physical model or all-frame GPU
convergence. It reuses saved dynamics/ionization and never creates a history.
"""
from pathlib import Path
import hashlib, json, os, sys
os.environ.setdefault('NUMBA_NUM_THREADS','3')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import numpy as np

def main():
    # Refuse missing histories rather than letting the historical model create one.
    inputs=json.loads((ROOT/'inputs.json').read_text())
    required=inputs['archives'][0]['files']+inputs['archives'][1]['files']+inputs['files']
    for row in required:
        path=ROOT/row['path']
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:
            raise RuntimeError('Restore exact numerical inputs first: '+row['path'])
    import src.physics_v3 as physical
    def refuse_generation(*args,**kwargs):
        raise RuntimeError('Selected ionization cache does not match; no automatic calculation permitted')
    physical.calculate_ions=refuse_generation
    Nebula=physical.Nebula
    from src.evolution_v3 import Track
    from src.render_v3 import project_material
    from src.metal_projector import MetalProjector
    track=Track();model=Nebula(nr=384,dynamics='radial-v3t-12x24-144.npz')
    state=model.state(float(track.clock(27.)))
    field=state['fields'];edges=state['edges'];nodes=np.unique(np.r_[0.,edges])
    width,height=144,256;halfwidth=float(edges[-1]/.82)
    basis=np.eye(3);star=np.zeros((height,width,3))
    args=[field,float(np.log(state['r'][0])),float(np.log(state['r'][1]/state['r'][0])),width,height,halfwidth,0.,star,basis,nodes,0.]
    reference=project_material(*args)
    gpu=MetalProjector(1)
    try:
        result=gpu(*args);again=gpu(*args)
    finally:
        gpu.close()
    scale=max(float(np.abs(reference).sum()),1e-200)
    error=float(np.abs(result-reference).sum()/scale)
    bias=float(abs(result.sum()-reference.sum())/scale)
    peak=float(np.abs(result-reference).max()/max(float(np.abs(reference).max()),1e-200))
    record=dict(passed=error<.001 and bias<.001 and peak<.005 and np.array_equal(result,again),
                raw_HDR_relative_L1=error,raw_HDR_flux_relative_bias=bias,
                max_absolute_error_over_peak=peak,repeat_bit_exact=bool(np.array_equal(result,again)),
                width=width,height=height,selected_state_playback_s=27.,
                thresholds=dict(HDR_L1=.001,flux_bias=.001,peak_normalized=.005),
                binary_sha256=hashlib.sha256((ROOT/'native/build/projector').read_bytes()).hexdigest(),
                scope='Bounded arithmetic check; original native-resolution qualification remains archived in validation')
    (ROOT/'cache').mkdir(exist_ok=True)
    (ROOT/'cache/projector-qualification.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))
    if not record['passed']:raise RuntimeError('CPU/Metal discrepancy; do not render with this build')

if __name__=='__main__':main()
