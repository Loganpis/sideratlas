"""Render the selected scene, using the original frame clock and sample counts."""
from pathlib import Path
import sys,argparse,json,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.verify_inputs import verify
from src.rendering.final_scene import create
from src.video_export import encoder
def main():
    p=argparse.ArgumentParser();p.add_argument('--width',type=int,default=1080);p.add_argument('--start',type=float,default=0);p.add_argument('--end',type=float,default=32);p.add_argument('--output',type=Path,default=ROOT/'output/reel.mp4');a=p.parse_args()
    assert 0<=a.start<a.end<=32 and a.width>=144
    if a.output.exists():raise RuntimeError('Do not overwrite an existing render')
    verify();a.output.parent.mkdir(parents=True,exist_ok=True);r=create(a.width);counts=json.loads((ROOT/'config/exposure-samples.json').read_text());lo=round(a.start*30);hi=round(a.end*30);enc=encoder(a.output,a.width,round(a.width*16/9),30)
    try:
        for i in range(lo,hi):
            enc.stdin.write(r.frame((i+.5)/30,30,counts[i]).astype('<u2').tobytes())
            if i%30==0:print(f'{i/30:.1f} / {a.end:.1f} seconds',flush=True)
    finally:enc.stdin.close()
    assert enc.wait()==0
    a.output.with_suffix('.json').write_text(json.dumps(dict(source_interval_s=[lo/30,hi/30],frames=hi-lo,width=a.width,height=r.H,fps=30,records=r.records,scope='Cached physics-informed visualization; no first-principles merger prediction'),indent=2)+'\n')
if __name__=='__main__':main()
