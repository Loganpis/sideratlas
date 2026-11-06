"""Restore only listed, SHA-256 verified inputs; never run a simulation.

Run from any directory. External provider terms still apply. --check is read-only;
--archive-dir can reuse already downloaded release assets without network access.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import tarfile
import tempfile
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
EPISODES = {'002': '002-pulsar', '003': '003-white-dwarf', '006': '006-planetary-nebula'}


def sha(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def valid(path, record):
    return path.is_file() and path.stat().st_size == record['bytes'] and sha(path) == record['sha256']


def target(base, relative):
    path = (base / relative).resolve()
    if not path.is_relative_to(ROOT):
        raise ValueError('Input destination leaves this repository')
    return path


def install(stream, path, record):
    if path.exists():
        raise FileExistsError(f'Existing input differs; preserve it for review: {path.relative_to(ROOT)}')
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as output:
        partial = Path(output.name)
        shutil.copyfileobj(stream, output)
    try:
        if not valid(partial, record):
            raise ValueError('Input checksum/size differs: ' + record['path'])
        partial.replace(path)
    finally:
        partial.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('episode', choices=EPISODES)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--archive-dir', type=Path)
    args = parser.parse_args()
    base = ROOT / 'episodes' / EPISODES[args.episode]
    manifest = json.loads((base / 'inputs.json').read_text())
    missing = []
    for archive in manifest.get('archives', []):
        needed = [r for r in archive['files'] if not valid(target(base, r['path']), r)]
        if args.check:
            missing.extend(r['path'] for r in needed)
            continue
        if not needed:
            continue
        with tempfile.TemporaryDirectory() as directory:
            source = (args.archive_dir / archive['name']) if args.archive_dir else Path(directory) / archive['name']
            if not source.exists():
                if args.archive_dir:
                    raise FileNotFoundError(source)
                print('Downloading', archive['name'], flush=True)
                urllib.request.urlretrieve(archive['url'], source)
            if not valid(source, archive):
                raise ValueError('Archive checksum/size differs: ' + archive['name'])
            if archive['format'] == 'zip':
                with zipfile.ZipFile(source) as contents:
                    for record in needed:
                        with contents.open(record['member']) as stream:
                            install(stream, target(base, record['path']), record)
            else:
                with tarfile.open(source, 'r:*') as contents:
                    for record in needed:
                        member = contents.getmember(record['member'])
                        if not member.isfile():
                            raise ValueError('Expected an ordinary input file')
                        with contents.extractfile(member) as stream:
                            install(stream, target(base, record['path']), record)
    for record in manifest.get('files', []):
        destination = target(base, record['path'])
        if valid(destination, record):
            continue
        if args.check:
            missing.append(record['path'])
        else:
            with urllib.request.urlopen(record['url'], timeout=120) as stream:
                install(stream, destination, record)
    for record in manifest.get('external_files', []):
        if not valid(target(base, record['path']), record):
            missing.append(record['path'])
    if missing:
        print('Additional inputs required:\n' + '\n'.join(missing))
        if args.check:
            raise SystemExit(1)
        print('No simulations were started. See the episode README for external-data reconstruction.')
    else:
        print('Every listed input matches its saved checksum')


if __name__ == '__main__':
    main()
