"""Explicit finite Slurm parent chain; dry-run unless --execute is supplied.

No hard-coded 7/30 day entitlement. The operator supplies a CURRENT verified
queue/QOS maximum; the actual EndTime is checked again by every worker.
"""
import argparse
import getpass
import json
import os
import re
from pathlib import Path
import shlex
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from clusterflow.common import writej
from clusterflow.runtime import site_config


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nodes',type=int,required=True)
    p.add_argument('--hours',type=int,required=True)
    p.add_argument('--verified-max-hours',type=int,required=True)
    p.add_argument('--waves',type=int,default=3)
    p.add_argument('--execute',action='store_true')
    a=p.parse_args()
    if not 1<=a.nodes<=64 or not 1<=a.waves<=10 or not 1<=a.hours<=a.verified_max_hours:
        raise ValueError('Invalid nodes/waves or requested walltime exceeds the verified site limit')
    site=site_config();slurm=site['slurm']
    base=['sbatch','--parsable','--account='+slurm['account'],'--partition='+slurm['partition'],
          '--qos='+slurm['qos'],'--nodes='+str(a.nodes),'--time='+str(a.hours)+':00:00']
    excluded=[path.stem for path in (ROOT/'operations/packed/quarantine').glob('*.json')]
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]+',node) for node in excluded):raise ValueError('Invalid quarantined node name')
    if excluded:base+=['--exclude='+','.join(sorted(excluded))]
    report={'status':'DRY_RUN_NO_SUBMISSION','nodes_per_parent':a.nodes,'parent_waves':a.waves,
            'allocated_cpus_per_parent':a.nodes*32,'walltime_hours':a.hours,
            'notice':'afterany is finite continuation, not guaranteed checkpoint resume; shared persistent state and audits survive early parent exit.',
            'example':shlex.join(base+[str(ROOT/'scripts/scnet_packed.sbatch')])}
    if not a.execute:
        print(json.dumps(report,ensure_ascii=False,indent=2));return
    if not os.environ.get('WORK2_RUNTIME_EVIDENCE'):raise ValueError('Runtime evidence environment variable required')
    queue=subprocess.run(['squeue','-h','-u',getpass.getuser(),'-o','%i'],capture_output=True,text=True,timeout=30,check=True)
    if len(queue.stdout.splitlines())+a.waves>190:
        raise ValueError('Leave headroom below historical 200-job QOS; no submission')
    env=os.environ.copy();env['WORK2_PACKAGE_ROOT']=str(ROOT)
    journal=ROOT/'operations'/('submission_'+str(time.time_ns())+'.json')
    submitted=[];writej(journal,{**report,'status':'SUBMISSION_INTENT','jobs':submitted})
    for _ in range(a.waves):
        command=base+(['--dependency=afterany:'+submitted[-1]] if submitted else [])+[str(ROOT/'scripts/scnet_packed.sbatch')]
        try:
            result=subprocess.run(command,env=env,capture_output=True,text=True,timeout=60,check=True)
        except (OSError,subprocess.SubprocessError) as exc:
            writej(journal,{**report,'status':'SUBMISSION_UNCERTAIN_QUERY_SQUEUE_BEFORE_RETRY',
                           'jobs':submitted,'reason':str(exc)})
            raise
        job=result.stdout.strip().split(';')[0]
        if not job.isdigit():raise ValueError('Unrecognized sbatch response; inspect journal and queue, never blindly resubmit')
        submitted.append(job)
        writej(journal,{**report,'status':'SUBMITTED_NOT_RUNNING','jobs':submitted})
    print(json.dumps({'journal':str(journal),'jobs':submitted,'status':'SUBMITTED_NOT_RUNNING'},ensure_ascii=False))


if __name__=='__main__':main()
