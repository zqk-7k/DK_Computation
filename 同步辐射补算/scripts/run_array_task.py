#!/usr/bin/env python3
"""One scheduler-array task; no shell interpolation of input paths."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--task-list',required=True)
    p.add_argument('--config',required=True)
    p.add_argument('--index',type=int)
    args=p.parse_args()
    tasks=json.loads(Path(args.task_list).read_text(encoding='utf-8'))
    index=args.index
    if index is None:
        for key in ['SLURM_ARRAY_TASK_ID','PBS_ARRAY_INDEX','PBS_ARRAYID']:
            if key in os.environ:
                index=int(os.environ[key]); break
    if index is None or not 0<=index<len(tasks):
        p.error('No valid zero-based scheduler array index was provided')
    root=Path(__file__).resolve().parent.parent
    config=Path(args.config).resolve()
    job=Path(tasks[index]['path'])
    if not job.is_absolute(): job=config.parent/job
    job=job.resolve()
    if not job.is_relative_to(config.parent.resolve()):
        raise ValueError('Job path must remain within the workflow directory')
    expected=tasks[index].get('feff_input_sha256')
    if expected and expected!=hashlib.sha256((job/'feff.inp').read_bytes()).hexdigest():
        raise RuntimeError('FEFF input changed after submission. Refusing mixed-input execution.')
    result=job/'run_result.json'
    if result.exists():
        from workflow import job_status
        if job_status(tasks[index],config.parent)=='complete':
            print('Already complete:',job); return
    completed=subprocess.run([sys.executable,str(root/'scripts'/'run_feff_job.py'),
                              '--job',str(job),'--config',str(config)])
    raise SystemExit(completed.returncode)


if __name__=='__main__': main()
