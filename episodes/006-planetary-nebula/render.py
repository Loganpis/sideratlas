"""Rerender selected global frame indices from frozen inputs; no physical solves."""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys
ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parents[1]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--stop', type=int, default=960)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 0 <= args.start < args.stop <= 960:
        raise ValueError('Use valid global frame indices, with no retiming')
    subprocess.run([sys.executable, str(REPOSITORY/'tools/fetch_inputs.py'), '006', '--check'], check=True)
    projector = ROOT/'native/build/projector'
    qualification = ROOT/'cache/projector-qualification.json'
    if not projector.is_file() or not qualification.is_file():
        raise RuntimeError('Build and qualify the projector first; see README')
    record = json.loads(qualification.read_text())
    if not record['passed'] or record['binary_sha256'] != hashlib.sha256(projector.read_bytes()).hexdigest():
        raise RuntimeError('This projector build is not qualified against the CPU reference')
    if args.output.exists():
        raise FileExistsError('Preserve existing output; choose a fresh filename')
    for name in ['cache', 'validation']:
        (ROOT/name).mkdir(exist_ok=True)
    os.environ.setdefault('NUMBA_CACHE_DIR', str(ROOT/'cache/numba'))
    os.environ.setdefault('NUMBA_NUM_THREADS', '3')
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
    sys.path.insert(0, str(ROOT))
    import src.physics_v3 as physical
    def refuse_generation(*args, **kwargs):
        raise RuntimeError('Selected ionization cache differs; no automatic calculation permitted')
    physical.calculate_ions = refuse_generation
    from src.production_giant_cells import GiantCellRenderer
    from src.video_export import encoder
    renderer = None
    process = None
    records = []
    try:
        renderer = GiantCellRenderer(spatial_samples=4, temporal_samples=2)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        process = encoder(args.output, 1080, 1920, 30)
        for index in range(args.start, args.stop):
            image, ledger = renderer.production_frame(index)
            process.stdin.write(image.tobytes())
            records.append(ledger)
            if index % 30 == 0:
                print('Rendered global frame', index, flush=True)
        process.stdin.close()
        if process.wait():
            raise RuntimeError('Encoder failed')
        args.output.with_suffix('.json').write_text(json.dumps(dict(start=args.start, stop=args.stop, frames=records), indent=2)+'\n')
    finally:
        if process is not None and not process.stdin.closed:
            process.stdin.close()
            process.wait()
        if renderer is not None:
            renderer.close()

if __name__ == '__main__':
    main()
