"""Identify the released, unchanged full-resolution film."""
from pathlib import Path
import argparse,json,subprocess,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.verify_inputs import sha
def main():
    p=argparse.ArgumentParser();p.add_argument('video',type=Path);a=p.parse_args();expected=json.loads((ROOT/'validation/final-checks/encoded-video.json').read_text())
    result=json.loads(subprocess.check_output(['ffprobe','-v','error','-count_frames','-show_streams','-show_format','-of','json',str(a.video)],text=True));s=result['streams'][0]
    assert sha(a.video)==expected['movie_sha256'];assert (s['width'],s['height'])==(1080,1920) and s['r_frame_rate']=='30/1' and int(s['nb_read_frames'])==960 and abs(float(result['format']['duration'])-32)<.001
    print('Exact selected movie verified: 32 s, 1080 × 1920, 30 fps, 960 frames')
if __name__=='__main__':main()
