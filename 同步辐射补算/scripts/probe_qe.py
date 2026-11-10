#!/usr/bin/env python3
"""Reject serial QE before launching multiple ranks in a shared scratch directory."""
from pathlib import Path
import argparse,json,re,subprocess,tempfile,hashlib

def probe(executable,ranks=16):
    executable=Path(executable).resolve()
    if not executable.is_file():raise RuntimeError('QE executable not found')
    with tempfile.TemporaryDirectory(prefix='qe_banner_probe_') as directory:
        p=subprocess.run([str(executable),'-h'],cwd=directory,input='',capture_output=True,text=True,timeout=30)
    text=p.stdout+'\n'+p.stderr
    mpi=bool(re.search(r'Parallel\s+version[^\n]*MPI',text,re.I))
    if int(ranks)>1 and not mpi:
        raise RuntimeError('Could not verify an MPI-enabled QE banner. Inspect pw.x -h and load the matching MPI build before multi-rank execution; no blind mpirun launch.')
    return {'executable':str(executable),'sha256':hashlib.sha256(executable.read_bytes()).hexdigest(),
            'MPI_banner_verified':mpi,'requested_ranks':int(ranks),'banner_excerpt':text[:6000]}

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('executable');a.add_argument('--ranks',type=int,default=16);a.add_argument('--output');p=a.parse_args();r=probe(p.executable,p.ranks)
    if p.output:Path(p.output).write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8')
    print('QE MPI executable verified')
