"""Add a disjoint staged cohort without interrupting an existing parent.

The normal submit_packed guard is unchanged. This entry point requires explicit
scopes, one finite parent, shared canonical runtime identity and a project-wide
node cap. Unknown scheduler/submission state blocks submission.
"""
import argparse
import getpass
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from clusterflow.common import readj, writej, sha, inside
from clusterflow.runtime import site_config, check_storage
from scripts.packed_worker import queue_lock, selected_species

TERMINAL = {'COMPLETED', 'CANCELLED', 'FAILED', 'TIMEOUT', 'OUT_OF_MEMORY', 'NODE_FAIL', 'PREEMPTED'}
ACTIVE = {'PENDING', 'RUNNING', 'CONFIGURING', 'COMPLETING', 'SUSPENDED', 'RESIZING', 'REQUEUED'}


def scheduler_state(job):
    if not re.fullmatch(r'\d+', str(job)): raise ValueError('Invalid parent ID')
    result = subprocess.run(['sacct', '-n', '-P', '-j', str(job), '--format=JobIDRaw,State'],
                            capture_output=True, text=True, timeout=30, check=True)
    values = [row[1].split()[0].split('+')[0] for row in
              (line.split('|') for line in result.stdout.splitlines()) if len(row)>1 and row[0]==str(job)]
    if len(values)!=1 or values[0] not in TERMINAL|ACTIVE:
        raise ValueError('Unknown parent state: '+str(job))
    return values[0]


def active_inventory(root):
    active = []
    for path in sorted((root/'operations').glob('submission_*.json')):
        record = readj(path)
        if record.get('status') != 'SUBMITTED_NOT_RUNNING' or not record.get('jobs'):
            raise ValueError('Reconcile unresolved submission before expansion: '+str(path))
        for job in record['jobs']:
            state = scheduler_state(job)
            if state in TERMINAL: continue
            if record.get('scope_file'):
                scope = inside(root, record['scope_file'])
                if sha(scope) != record.get('scope_sha256'): raise ValueError('Active scope changed')
                ids = selected_species(root, scope)
            elif record.get('pilot') is True:
                scope = root/'config/pilot_scope_20261014.json'
                ids = selected_species(root, scope)
            else:
                raise ValueError('Active parent has no provably disjoint scope: '+str(job))
            nodes = record.get('nodes_per_parent')
            if type(nodes) is not int or nodes<1: raise ValueError('Active node budget unknown')
            active.append({'job_id':str(job), 'nodes':nodes, 'species_ids':ids, 'state':state})
    return active


def validate_expansion(active, species, nodes, max_nodes):
    if not 1<=nodes<=max_nodes<=12: raise ValueError('Invalid bounded node request')
    if sum(row['nodes'] for row in active)+nodes>max_nodes:
        raise ValueError('Project node cap would be exceeded')
    for row in active:
        overlap=set(species)&set(row['species_ids'])
        if overlap: raise ValueError('Scope overlaps active parent '+row['job_id']+': '+','.join(sorted(overlap)))


def submit(root, name, scope, nodes, *, execute=False, require_goat=False, max_nodes=12):
    root=Path(root).resolve(); scope=Path(scope).resolve()
    if not re.fullmatch(r'[a-z0-9_]{1,48}',name): raise ValueError('Invalid cohort name')
    if not scope.is_relative_to(root): raise ValueError('Scope must be inside package')
    species=selected_species(root,scope)
    if len(species)<nodes*4: raise ValueError('Scope too small for requested balanced slots')
    if os.environ.get('WORK2_PACKED_PROFILE','balanced')!='balanced': raise ValueError('Staged cohort requires balanced profile')
    site=site_config(); slurm=site['slurm']
    runtime=Path(os.environ.get('WORK2_RUNTIME_EVIDENCE',''))
    if not runtime.is_file() or not runtime.resolve().is_relative_to(root):
        raise ValueError('Canonical runtime evidence inside this package required')
    evidence=readj(runtime)
    if evidence.get('version')!='6.1.0' or evidence.get('site_sha256')!=sha(root/'config/site_scnet.json'):
        raise ValueError('Canonical runtime does not match site configuration')
    journal=root/'operations'/('submission_cohort_'+name+'.json')
    with queue_lock(root/'operations/submission_lock'):
        if journal.exists():
            old=readj(journal)
            if old.get('status')!='SUBMITTED_NOT_RUNNING' or not old.get('jobs'):
                raise ValueError('Uncertain existing cohort submission; reconcile rather than retry')
            if old.get('scope_sha256')!=sha(scope) or old.get('nodes_per_parent')!=nodes:
                raise ValueError('Existing cohort request differs')
            return {'status':'ALREADY_SUBMITTED','jobs':old['jobs'],'states':[scheduler_state(j) for j in old['jobs']]}
        active=active_inventory(root)
        validate_expansion(active,species,nodes,max_nodes)
        if require_goat:
            from clusterflow.stage_release import audit
            checked=audit(root)
            if not checked['validated_goat']:
                return {'status':'WAITING_VALID_GOAT','quantum_jobs_launched':0,'nodes_submitted':0}
        excluded=sorted(p.stem for p in (root/'operations/packed/quarantine').glob('*.json'))
        if any(not re.fullmatch(r'[A-Za-z0-9_.-]+',n) for n in excluded):raise ValueError('Invalid node exclusion')
        base=['sbatch','--parsable','--account='+slurm['account'],'--partition='+slurm['partition'],
              '--qos='+slurm['qos'],'--nodes='+str(nodes),'--time=720:00:00','--cpus-per-task=32',
              '--mem='+str(site['node']['requested_memory_mb'])+'M','--job-name=work2-64-'+name,
              '--chdir='+str(root)]
        if excluded:base.append('--exclude='+','.join(excluded))
        report={'status':'DRY_RUN_NO_SUBMISSION','name':name,'nodes_per_parent':nodes,'parent_waves':1,
                'walltime_hours':720,'scope_file':str(scope.relative_to(root)),'scope_sha256':sha(scope),
                'species_ids':species,'pilot':False,'staged_release':True,'priority_policy':'remaining_work',
                'project_max_nodes':max_nodes,'coexisting_parents':active,'excluded_nodes':excluded,
                'canonical_runtime_evidence':str(runtime.resolve()),'runtime_sha256':sha(runtime),
                'jobs':[],'created_epoch':time.time()}
        if not execute:return report
        total_slots=(sum(row['nodes'] for row in active)+nodes)*4
        storage=check_storage(site['shared_scratch'],site['disk']['admit_free_gib']+
                              total_slots*site['disk']['unmeasured_task_growth_gib'])
        report['storage']=storage
        launcher=root/'scripts/scnet_cohort.sbatch'
        validation=subprocess.run(base+['--test-only',str(launcher)],capture_output=True,text=True,timeout=60,check=True)
        report['scheduler_test_only']={'stdout':validation.stdout,'stderr':validation.stderr}
        queue=subprocess.run(['squeue','-h','-u',getpass.getuser(),'-o','%i'],capture_output=True,text=True,timeout=30,check=True)
        if len(queue.stdout.splitlines())>=190: raise ValueError('Leave queue headroom below 200')
        env=os.environ.copy()
        env.update({'WORK2_PACKAGE_ROOT':str(root),'WORK2_PYTHON':site['python'],
                    'WORK2_COHORT_SCOPE':str(scope),'WORK2_COHORT_SCOPE_SHA256':sha(scope),
                    'WORK2_STAGED_RELEASE':'1','WORK2_PRIORITY_POLICY':'remaining_work',
                    'WORK2_RUNTIME_EVIDENCE':str(runtime.resolve())})
        writej(journal,{**report,'status':'SUBMISSION_INTENT'})
        try:
            result=subprocess.run(base+[str(launcher)],env=env,capture_output=True,text=True,timeout=60,check=True)
            job=result.stdout.strip().split(';')[0]
            if not job.isdigit():raise ValueError('Unrecognized sbatch output')
        except (OSError,subprocess.SubprocessError,ValueError) as exc:
            writej(journal,{**report,'status':'SUBMISSION_UNCERTAIN_QUERY_SQUEUE_BEFORE_RETRY','reason':str(exc)})
            raise
        writej(journal,{**report,'status':'SUBMITTED_NOT_RUNNING','jobs':[job]})
        return {'status':'SUBMITTED_NOT_RUNNING','jobs':[job],'journal':str(journal)}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--name',required=True);p.add_argument('--scope',type=Path,required=True)
    p.add_argument('--nodes',type=int,required=True);p.add_argument('--execute',action='store_true')
    p.add_argument('--require-goat',action='store_true')
    a=p.parse_args()
    print(json.dumps(submit(ROOT,a.name,a.scope,a.nodes,execute=a.execute,require_goat=a.require_goat),ensure_ascii=False,indent=2))


if __name__=='__main__':main()
