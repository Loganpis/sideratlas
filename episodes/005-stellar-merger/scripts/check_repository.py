"""Static release hygiene and Markdown navigation checks, not a physics test."""
from pathlib import Path
import ast,json,re,subprocess,argparse,urllib.parse
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'validation/final-checks/repository-audit.json');a=p.parse_args()
    files=[ROOT/x for x in subprocess.check_output(['git','ls-files','--cached'],cwd=ROOT,text=True).splitlines()];problems=[];links=0;total=0;maxbytes=0
    for path in files:
        size=path.stat().st_size;total+=size;maxbytes=max(maxbytes,size)
        if size>20_000_000:problems.append('Unexpected large Git file: '+str(path.relative_to(ROOT)))
        if any(x in path.parts for x in ['.venv','__pycache__','cache']):problems.append('Unwanted directory tracked: '+str(path.relative_to(ROOT)))
        if path.suffix not in ['.py','.md','.json','.txt','.text','.data','']:continue
        try:text=path.read_text()
        except UnicodeDecodeError:continue
        for marker in ['/Users/','/private/tmp/','/var/folders/','github_pat_','-----BEGIN PRIVATE KEY-----']:
            # This checker necessarily declares the patterns it checks for.
            if path.name!='check_repository.py' and marker in text:problems.append('Sensitive/local marker in '+str(path.relative_to(ROOT)))
        if re.search(r'gh[pousr]_[A-Za-z0-9]{25,}|AKIA[A-Z0-9]{16}',text):problems.append('Credential pattern in '+str(path.relative_to(ROOT)))
        if path.suffix=='.py':ast.parse(text)
        if path.suffix=='.md':
            for target in re.findall(r'\]\(([^)]+)\)',text)+re.findall(r'(?:src|href)="([^"]+)"',text):
                target=target.strip('<>').split('#')[0]
                if not target or re.match(r'^[a-zA-Z]+:',target):continue
                links+=1
                if (path.parent/urllib.parse.unquote(target)).resolve()!=a.output.resolve() and not (path.parent/urllib.parse.unquote(target)).exists():problems.append('Broken relative link: '+str(path.relative_to(ROOT))+' -> '+target)
    result=dict(passed=not problems,tracked_files=len(files),tracked_bytes=total,largest_tracked_file_bytes=maxbytes,relative_links_checked=links,problems=problems,scope='Static source/link/size/sensitive-pattern checks; no credentials printed, no astrophysical validation')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True);assert result['passed']
if __name__=='__main__':main()
