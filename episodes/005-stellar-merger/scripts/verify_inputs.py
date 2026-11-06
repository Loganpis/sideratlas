"""Verify exact cached inputs without loading NumPy or a renderer."""
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1]
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def verify(required=True):
    rows=json.loads((ROOT/'data/input-manifest.json').read_text())['files'];missing=[]
    for row in rows:
        p=ROOT/row['path']
        if not p.is_file():missing.append(row['path']);continue
        if sha(p)!=row['sha256']:raise RuntimeError('Input checksum mismatch: '+row['path'])
    if missing and required:raise RuntimeError('Missing inputs:\n'+'\n'.join(missing)+'\nSee METHODS.md; two surface inputs are intentionally not redistributed')
    return missing
if __name__=='__main__':
    missing=verify();print('All exact rendering inputs verified' if not missing else missing)
