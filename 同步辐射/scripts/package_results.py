#!/usr/bin/env python3
"""Package verified calculated outputs; exclude bulky reproducible scratch files."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile
from workflow import read_json, job_dir, job_status, stamp


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--config',required=True)
    p.add_argument('--profile',required=True)
    a=p.parse_args();config_path=Path(a.config).resolve();root=config_path.parent;c=read_json(config_path)
    manifest=root/c.get('paths',{}).get('manifests_dir','manifests')/f'{a.profile}_jobs.json'
    jobs=read_json(manifest)
    if isinstance(jobs,dict):jobs=jobs['jobs']
    incomplete=[j['id'] for j in jobs if job_status(j,root)!='complete']
    if incomplete:raise RuntimeError('Refusing final package with incomplete jobs: '+', '.join(incomplete[:5]))
    results=root/c.get('postprocess',{}).get('output_dir','results')/a.profile
    complete=sorted(p for p in results.iterdir() if p.is_dir() and p.name.startswith('COMPLETE'))
    if not complete:raise RuntimeError('No complete postprocessed result directory exists')
    latest=complete[-1]
    files=set()
    for directory in ['structures','scripts','docs','validation']:
        files.update(p for p in (root/directory).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    files.update(p for p in latest.rglob('*') if p.is_file())
    for job in jobs:
        files.update(p for p in job_dir(job,root).iterdir() if p.is_file() and p.suffix in {'.inp','.dat','.json','.log'})
    for name in ['README.md','AGENTS.md','CODEX_TASK.txt','config.json','production_config.json','requirements.txt',
                 'structure_checksums.json','run_all.sh','software/python_environment.lock.txt',
                 'software/feff10/install_complete.json','software/feff10/source/LICENSE']:
        path=root/name
        if path.is_file():files.add(path)
    files.add(manifest);files.add(config_path)
    listing=io.StringIO();writer=csv.writer(listing);writer.writerow(['relative_path','bytes','sha256'])
    for path in sorted(files):
        with path.open('rb') as f:checksum=hashlib.file_digest(f,'sha256').hexdigest()
        writer.writerow([path.relative_to(root).as_posix(),path.stat().st_size,checksum])
    folder=root/'deliverables';folder.mkdir(exist_ok=True)
    archive=folder/f'CALCULATED_XAS_{a.profile}_{stamp()}.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        z.comment=b'Calculated static-geometry XAS; not measured operando spectra'
        for path in sorted(files):z.write(path,arcname=path.relative_to(root).as_posix())
        z.writestr('SOURCE_DATA_SHA256.csv',listing.getvalue())
        z.writestr('PACKAGE_SCOPE.txt',
                   'Calculated FEFF spectra with source CIFs, inputs, raw textual outputs, logs, processed data and plots.\n'
                   'Potential assignments and original VASP electronic-state transfer are not established.\n'
                   'Reproducible large scratch binaries are omitted; path tables and software provenance are retained.\n')
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:raise RuntimeError('ZIP integrity failed')
    print(json.dumps({'archive':str(archive),'files':len(files)+2,'bytes':archive.stat().st_size}))


if __name__=='__main__':main()
