"""Read-only source, link, provenance and release-registry audit; no physics runs."""
from pathlib import Path
import ast,hashlib,json,re,subprocess,urllib.parse
ROOT=Path(__file__).resolve().parents[1]

def main():
    raw=subprocess.check_output(['git','ls-files','-z'],cwd=ROOT)
    files=[ROOT/p.decode() for p in raw.split(b'\0') if p]
    problems=[];links=0;size=0
    for p in files:
        rel=str(p.relative_to(ROOT));size+=len(p.read_bytes())
        if p.is_symlink() and (not p.exists() or not p.resolve().is_relative_to(ROOT)):
            problems.append('Invalid symlink: '+rel)
        if len(p.read_bytes())>10_000_000 or p.suffix in ['.mp4','.mov','.zip','.dylib','.so','.o','.log']:
            problems.append('Unwanted large/generated file: '+rel)
        if any(s in p.parts for s in ['__pycache__','.venv','native-work','build','cache']):
            problems.append('Generated directory: '+rel)
        if p.suffix=='.py':ast.parse(p.read_text())
        if p.suffix=='.json':json.loads(p.read_text())
        if p.suffix not in ['.py','.md','.json','.text','.txt','.toml','.cpp','.swift','.metal','']:
            continue
        try:text=p.read_text()
        except UnicodeDecodeError:continue
        if p.name not in ['check_release.py','check_repository.py']:
            if any(s in text for s in ['/Users/','/private/tmp/','/var/folders/','-----BEGIN PRIVATE KEY-----']):
                problems.append('Private/local marker: '+rel)
            if re.search(r'gh[pousr]_[A-Za-z0-9]{25,}|github_pat_[A-Za-z0-9_]{25,}|AKIA[A-Z0-9]{16}',text):
                problems.append('Credential-shaped text: '+rel)
        if p.suffix=='.md':
            targets=re.findall(r'\]\(([^)]+)\)',text)+re.findall(r'(?:src|href)="([^"]+)"',text)
            for t in targets:
                t=t.strip('<>').split('#')[0]
                if not t or re.match(r'^[a-zA-Z]+:',t):continue
                links+=1
                if not (p.parent/urllib.parse.unquote(t)).exists():problems.append('Missing link: '+rel+' -> '+t)
    episodes=json.loads((ROOT/'releases.json').read_text())
    assert [e['number'] for e in episodes]==['001','002','003','004','005','006']
    for e in episodes:
        p=ROOT/e['path'];j=json.loads((p/'provenance.json').read_text())
        if e['sha256']!=j['sha256'] or e['bytes']!=j['bytes']:problems.append('Video registry mismatch: '+e['number'])
        for r in j['source_files']:
            path=p/r['path']
            if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=r['sha256']:
                problems.append('Source fingerprint mismatch: '+str(path.relative_to(ROOT)))
        if e['url'] not in (p/'README.md').read_text():problems.append('Episode video link mismatch: '+e['number'])
    result=dict(passed=not problems,tracked_files=len(files),tracked_worktree_bytes=size,relative_links_checked=links,
                released_video_bytes=sum(e['bytes'] for e in episodes),problems=problems,
                scope='Read-only static collection checks; not a new astrophysical validation')
    print(json.dumps(result,indent=2))
    if problems:raise SystemExit(1)

if __name__=='__main__':main()
