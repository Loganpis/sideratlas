"""Download and safely restore checksum-identified permitted release assets."""
from pathlib import Path
import sys,json,urllib.request,zipfile,tempfile,argparse,shutil
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.verify_inputs import sha,verify
def restore_archive(path,expected):
    if sha(path)!=expected:raise RuntimeError('Release archive checksum mismatch')
    with zipfile.ZipFile(path) as z:
        for item in z.infolist():
            dest=(ROOT/item.filename).resolve()
            if not dest.is_relative_to(ROOT) or item.is_dir() or item.external_attr>>16&0o170000==0o120000:raise RuntimeError('Unexpected archive member')
            dest.parent.mkdir(parents=True,exist_ok=True)
            with z.open(item) as source,dest.open('wb') as target:shutil.copyfileobj(source,target)
def main():
    p=argparse.ArgumentParser();p.add_argument('--archive',type=Path,help='Use the downloaded render-inputs.zip instead of network');p.add_argument('--local-surface-source',type=Path,help='Original project root containing the two undistributed, checksum-matching maps');a=p.parse_args()
    manifest=json.loads((ROOT/'data/release-assets.json').read_text())
    if a.archive:restore_archive(a.archive,manifest['sha256'])
    else:
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'render-inputs.zip';urllib.request.urlretrieve(manifest['url'],path);restore_archive(path,manifest['sha256'])
    if a.local_surface_source:
        for row in json.loads((ROOT/'data/input-manifest.json').read_text())['files']:
            if row['distribution']!='external-surface-input':continue
            source=a.local_surface_source/row['path'];assert sha(source)==row['sha256'],row['path'];dest=ROOT/row['path'];dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest)
    missing=verify(required=False)
    if missing:print('Permitted inputs restored. Externally supplied surface maps still required:\n'+'\n'.join(missing)+'\nSee METHODS.md')
    else:print('All rendering inputs restored and verified')
if __name__=='__main__':main()
