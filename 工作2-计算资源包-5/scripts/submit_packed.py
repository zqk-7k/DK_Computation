"""Explicit finite Slurm parent chain; dry-run unless --execute is supplied.

Default request is 30 days, validated by live sbatch --test-only on execute.
This is not an entitlement; every worker checks the actual allocation EndTime.
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
from clusterflow.common import readj, writej
from clusterflow.runtime import site_config
from clusterflow.operations import scheduler_terminal
from scripts.packed_worker import queue_lock


def verify_previous_submissions(root):
    """Never stack another parent chain on active or uncertain submissions."""
    for path in sorted((root/'operations').glob('submission_*.json')):
        record=readj(path)
        if record.get('status')!='SUBMITTED_NOT_RUNNING' or not record.get('jobs'):
            raise ValueError('Unresolved submission journal; reconcile queue before retry: '+str(path))
        if any(not scheduler_terminal({'slurm_job_id':job}) for job in record['jobs']):
            raise ValueError('Previous parent chain active or not confirmed terminal: '+str(path))


def submit_chain(base,report,waves,env,launcher=None):
    # Separate submission mutex: long scheduler calls must not block workers'
    # queue mutex, and concurrent launchers must not create duplicate chains.
    with queue_lock(ROOT/'operations/submission_lock'):
        verify_previous_submissions(ROOT)
        journal=ROOT/'operations'/('submission_'+str(time.time_ns())+'.json')
        submitted=[];writej(journal,{**report,'status':'SUBMISSION_INTENT','jobs':submitted})
        for _ in range(waves):
            command=base+(['--dependency=afterany:'+submitted[-1]] if submitted else [])+[str(launcher or ROOT/'scripts/scnet_packed.sbatch')]
            try:
                result=subprocess.run(command,env=env,capture_output=True,text=True,timeout=60,check=True)
                job=result.stdout.strip().split(';')[0]
                if not job.isdigit():raise ValueError('Unrecognized sbatch response: '+result.stdout)
            except (OSError,subprocess.SubprocessError,ValueError) as exc:
                writej(journal,{**report,'status':'SUBMISSION_UNCERTAIN_QUERY_SQUEUE_BEFORE_RETRY',
                               'jobs':submitted,'reason':str(exc)})
                raise
            submitted.append(job)
            # Intermediate success is not a fully accounted-for chain. If
            # this process dies, require reconciliation instead of resubmission.
            writej(journal,{**report,'status':'SUBMISSION_IN_PROGRESS','jobs':submitted})
        writej(journal,{**report,'status':'SUBMITTED_NOT_RUNNING','jobs':submitted})
        return {'journal':str(journal),'jobs':submitted,'status':'SUBMITTED_NOT_RUNNING'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nodes',type=int,required=True)
    p.add_argument('--hours',type=int,help='Default 720 hours (30 days); live sbatch test-only required on execute')
    p.add_argument('--verified-max-hours',type=int,help='Optional additional operator ceiling; never replaces live validation')
    p.add_argument('--profile',choices=['balanced','large_memory','branch_parallel'],default='balanced')
    p.add_argument('--waves',type=int,default=3)
    p.add_argument('--execute',action='store_true')
    p.add_argument('--pilot',action='store_true',help='Bounded real-input pilot; bootstrap runtime evidence on compute node')
    a=p.parse_args()
    if a.pilot and (a.nodes>2 or a.waves!=1):
        raise ValueError('Pilot requires at most two nodes and --waves 1; audit before continuation or expansion')
    launcher=ROOT/'scripts'/('scnet_pilot.sbatch' if a.pilot else 'scnet_packed.sbatch')
    os.environ['WORK2_PACKED_PROFILE']=a.profile
    site=site_config();slurm=site['slurm']
    a.hours=a.hours if a.hours is not None else slurm.get('default_hours',720)
    if not 1<=a.nodes<=64 or not 1<=a.waves<=10 or a.hours<1 or (a.verified_max_hours is not None and a.hours>a.verified_max_hours):
        raise ValueError('Invalid nodes/waves or requested walltime exceeds the supplied ceiling')
    if a.nodes*site['packed']['slots_per_node']>66:
        raise ValueError('More branch workers than product branches; choose fewer nodes or branch_parallel profile')
    base=['sbatch','--parsable','--account='+slurm['account'],'--partition='+slurm['partition'],
          '--qos='+slurm['qos'],'--nodes='+str(a.nodes),'--time='+str(a.hours)+':00:00',
          '--cpus-per-task='+str(site['packed']['slots_per_node']*site['packed']['slot_allocated_cpus']),
          '--mem='+str(site['node']['requested_memory_mb'])+'M']
    excluded=[path.stem for path in (ROOT/'operations/packed/quarantine').glob('*.json')]
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]+',node) for node in excluded):raise ValueError('Invalid quarantined node name')
    if excluded:base+=['--exclude='+','.join(sorted(excluded))]
    report={'status':'DRY_RUN_NO_SUBMISSION','nodes_per_parent':a.nodes,'parent_waves':a.waves,
            'allocated_cpus_per_parent':a.nodes*32,'walltime_hours':a.hours,
            'profile':a.profile,'slots_per_node':site['packed']['slots_per_node'],
            'notice':'afterany is finite continuation, not guaranteed checkpoint resume; shared persistent state and audits survive early parent exit.',
            'pilot':a.pilot,
            'example':shlex.join(base+[str(launcher)])}
    if not a.execute:
        print(json.dumps(report,ensure_ascii=False,indent=2));return
    if not a.pilot and not os.environ.get('WORK2_RUNTIME_EVIDENCE'):raise ValueError('Runtime evidence environment variable required')
    # Scheduler checks the actual account/partition/QOS combination without
    # submitting a job. Refuse rather than silently reduce a long request.
    validation=subprocess.run(base+['--test-only',str(launcher)],
                              capture_output=True,text=True,timeout=60,check=True)
    report['scheduler_test_only']={'stdout':validation.stdout,'stderr':validation.stderr,'returncode':validation.returncode}
    queue=subprocess.run(['squeue','-h','-u',getpass.getuser(),'-o','%i'],capture_output=True,text=True,timeout=30,check=True)
    if len(queue.stdout.splitlines())+a.waves>190:
        raise ValueError('Leave headroom below historical 200-job QOS; no submission')
    env=os.environ.copy();env['WORK2_PACKAGE_ROOT']=str(ROOT)
    env['WORK2_PYTHON']=site['python']
    print(json.dumps(submit_chain(base,report,a.waves,env,launcher),ensure_ascii=False))


if __name__=='__main__':main()
