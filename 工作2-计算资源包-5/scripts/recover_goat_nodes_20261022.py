"""Explicit, one-wave NODE_FAIL recovery; retain old scratch and scientific inputs.

The 2026-10-22 operator-authorized recovery uses existing source workspaces,
new attempt tokens, and separate immutable priority scopes. No infinite retry.
"""
import argparse
import contextlib
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from clusterflow.common import readj, writej, sha, inside
from clusterflow.priority import estimate
from clusterflow.scheduling import single_execution
from clusterflow.goat_retry import retry_token
from scripts.goat_job import validate_job, INPUT_HASH_FIELDS
from scripts.goat_portable import IMMUTABLE_FILES
from scripts.packed_worker import queue_lock
from scripts.submit_goat_bulk import IDENTITIES

BASE = 'operations/retry_20261022'
AUDIT = 'operations/recovery_audit_20261022'
ROOT_NAMES = ('', 'operations/expansion_20261020', 'operations/recovery_cohort_20261020')


def exact_failed(accounting, jid):
    matches = [r for r in accounting if r[0] == str(jid)]
    if len(matches) != 1 or matches[0][2] != 'NODE_FAIL':
        raise ValueError('Exact failed allocation not established: ' + str(jid))
    if any(r[0].startswith(str(jid)+'.') and r[2] in ('RUNNING','PENDING','COMPLETING') for r in accounting):
        raise ValueError('Prior allocation still has active steps')
    return matches[0]


def scheduler(ids):
    p = subprocess.run(['sacct','-j',','.join(sorted(set(ids))),'-S','2026-10-14',
                        '-n','-P','--format=JobIDRaw,JobID,State,End,NodeList'],
                       capture_output=True,text=True,check=True,timeout=90)
    return [l.split('|') for l in p.stdout.splitlines() if l.strip()]


def ranking(root, jobs):
    rows = {r['job_id']:r for r in readj(root/'registry/jobs.json')}
    branches = {r['species_id']:r for r in readj(root/'registry/batch_plan_v5.json')['species']}
    scores = {rows[j]['species_id']:None for j in jobs}
    for sid in scores:
        scores[sid] = estimate(root, branches[sid], rows)
    ordered = sorted(jobs, key=lambda j: (-scores[rows[j]['species_id']]['remaining_serial_hours_estimate'],
                                         rows[j]['species_id'],j))
    return ordered, scores


def prepare(root, probe_path):
    root = Path(root).resolve(); audit = root/AUDIT; journal = audit/'PREPARED.json'
    if journal.exists():
        return readj(journal)
    if (audit/'PREPARATION_INTENT.json').exists():
        raise ValueError('Incomplete recovery preparation; reconcile journal before continuing')
    probe = readj(probe_path); proofs = probe['proofs']
    if len(proofs) != 137 or len({p['job_id'] for p in proofs}) != 137:
        raise ValueError('Expected the 137 explicitly reviewed interrupted starts')
    allowed = {str((root/n).resolve()) for n in ROOT_NAMES}
    if any(p['root'] not in allowed for p in proofs):
        raise ValueError('Unexpected source workspace')
    accounting = scheduler([str(p['owner']['slurm_job_id']) for p in proofs])
    roots = [root/n for n in ROOT_NAMES]
    indexed = {str(r):{j['job_id']:j for j in readj(r/'registry/jobs.json')} for r in roots}
    old_bases = [r/'operations/goat_bulk_20261015' for r in roots]
    old_bases.append(roots[1]/'operations/pilot_20261020')
    with queue_lock(root/'operations/submission_lock'), contextlib.ExitStack() as locks:
        for b in sorted(old_bases):locks.enter_context(queue_lock(b))
        queues = {str(b):readj(b/'queue.json') for b in old_bases}
        for p in proofs:
            locks.enter_context(single_execution(root,'goat_'+hashlib.sha256(p['job_id'].encode()).hexdigest()))
        verified = []
        for p in proofs:
            jid=p['job_id'];source=Path(p['root']);row=indexed[str(source)][jid];d=inside(source,row['job_directory'])
            if (d/'delegated_goat.json').exists() or (d/'goat_retry_authorization.json').exists():
                raise ValueError('Job already delegated/retried: '+jid)
            if sha(d/'.goat-run.lock')!=p['owner_sha256']:raise ValueError('Owner changed: '+jid)
            owner=readj(d/'.goat-run.lock');allocation=exact_failed(accounting,str(owner['slurm_job_id']))
            stage=Path(owner['portable_stage_root'])/'job'
            if str(stage)!=p['stage'] or sha(stage/'receipt.json')!=p['receipt_sha256']:
                raise ValueError('Prior attempt changed: '+jid)
            if readj(stage/'receipt.json')['status']!='RUNNING':raise ValueError('Recover terminal result first')
            if (stage/'goat.finalensemble.xyz').exists():raise ValueError('Review existing ensemble before retry')
            with (stage/'goat.out').open('rb') as f:
                f.seek(max(0,(stage/'goat.out').stat().st_size-32768));tail=f.read()
            if b'ORCA TERMINATED NORMALLY' in tail:raise ValueError('Review normal termination before retry')
            h=stage/'runtime_heartbeat.json'
            if h.exists():
                try:heartbeat=readj(h)
                except json.JSONDecodeError:heartbeat={}
                if time.time()-heartbeat.get('updated_epoch',0)<3600:raise ValueError('Fresh producer heartbeat')
            fields=dict(INPUT_HASH_FIELDS,**{'provenance.json':'provenance_sha256'})
            inputs={n:sha(d/n) for n in IMMUTABLE_FILES}
            if any(inputs[n]!=row[field] for n,field in fields.items()):raise ValueError('Input identity changed')
            containing=[b for b in old_bases if jid in queues[str(b)]['active']]
            if len(containing)!=1:raise ValueError('Missing/duplicate old claim: '+jid)
            b=containing[0];claim=queues[str(b)]['active'][jid]
            if str(claim['slurm_job_id'])!=str(owner['slurm_job_id']):raise ValueError('Claim allocation changed')
            if (source/'operations/packed/blocked'/(row['species_id']+'.json')).exists():
                raise ValueError('Species blocked; independent review required: '+jid)
            verified.append(dict(p,scheduler=allocation,scheduler_state='NODE_FAIL',old_base=str(b),
                                 inputs=inputs,manifest_sha256=sha(source/'registry/jobs.json')))
        # Transfer only unclaimed starts under the original queue lock.
        exp=roots[1];old=exp/'operations/goat_bulk_20261015';old_scope=readj(old/'scope.json');q=queues[str(old)]
        fresh=[]
        for jid in old_scope['job_ids'][q['cursor']:]:
            d=exp/indexed[str(exp)][jid]['job_directory']
            if any((d/n).exists() for n in ('.goat-run.lock','receipt.json','portable_stageout.json','delegated_goat.json')) or any((d/'stage_intents').glob('*.json')):
                raise ValueError('Pending start already owned: '+jid)
            validate_job(exp,jid);fresh.append(jid)
        all_ids=[p['job_id'] for p in verified]+fresh
        if len(set(all_ids))!=len(all_ids):raise ValueError('Duplicate recovery scope')
        ordered,scores=ranking(root,all_ids)
        token=uuid.uuid4().hex
        writej(audit/'PREPARATION_INTENT.json',dict(token=token,epoch=time.time(),proofs=verified,
               fresh=fresh,ranking=scores,ordered=ordered,authorization='User explicitly requested recovering all 137 interrupted starts, restoring monitoring and complex-first GOAT on 2026-10-22.'))
        groups=[]
        for source in roots:
            ids=[j for j in ordered if any(p['job_id']==j and p['root']==str(source) for p in verified) or (source==exp and j in fresh)]
            if not ids:continue
            cpus=8 if source==root else 16;base=source/BASE
            if (base/'scope.json').exists():raise ValueError('Existing retry scope; reconcile')
            scope=dict(schema=1,stage='goat',job_ids=ids,operations_directory=BASE,cpus_per_job=cpus,
                       slots_per_node=32//cpus,growth_gib_per_job=8,free_floor_gib=128,legacy_active_jobs=469,
                       identity_sha256={n:sha(source/n) for n in IDENTITIES},created_epoch=time.time(),
                       canonical_root=str(root),recovery_token=token,
                       notice='One authorized retry wave plus transferred initial starts; descending estimated remaining work. No acceptance or checkpoint claim.')
            writej(base/'scope.json',scope)
            groups.append(dict(root=str(source),base=str(base),jobs=ids,cpus_per_job=cpus,nodes=math.ceil(len(ids)/(32//cpus)),
                               priority_score=max(scores[indexed[str(source)][j]['species_id']]['remaining_serial_hours_estimate'] for j in ids)))
        for p in verified:
            jid=p['job_id'];source=Path(p['root']);d=source/p['job_directory'];b=Path(p['old_base'])
            archive=d/'recovery_history'/('node_fail_20261022_'+token);archive.mkdir(parents=True)
            writej(archive/'proof.json',p)
            (d/'.goat-run.lock').rename(archive/'original_owner.json')
            (d/'stage_intents').rename(archive/'stage_intents')
            writej(archive/'old_active_claim.json',queues[str(b)]['active'].pop(jid))
            cp=b/'claims'/(jid+'.json')
            if cp.exists():cp.rename(archive/'old_pool_claim.json')
            outcome=b/'results'/(jid+'.json')
            if outcome.exists():outcome.rename(archive/'old_pool_result.json')
            writej(outcome,dict(job_id=jid,finished_epoch=time.time(),result=dict(status='NODE_FAIL_RETRY_SCHEDULED',new_base=str(source/BASE),prior_stage=p['stage'])))
            authorization=dict(status='PREPARED_VERIFIED_NODE_FAIL',job_id=jid,attempt_token=uuid.uuid4().hex,
                               manifest_sha256=p['manifest_sha256'],inputs=p['inputs'],
                               proof_path=str((archive/'proof.json').relative_to(source)),proof_sha256=sha(archive/'proof.json'),created_epoch=time.time())
            writej(d/'goat_retry_authorization.json',authorization)
            validate_job(source,jid);retry_token(source,d,jid,p['manifest_sha256'],p['inputs'])
        for jid in fresh:q['skipped'].append(dict(job_id=jid,reason='TRANSFERRED_TO_COMPLEX_FIRST_SCOPE',new_base=str(exp/BASE)))
        q['cursor']=len(old_scope['job_ids'])
        for b in old_bases:writej(b/'queue.json',queues[str(b)])
        # Move only records belonging to conclusively failed allocations into history.
        failed_ids={str(p['owner']['slurm_job_id']) for p in verified}
        for b in old_bases:
            for folder in ('nodes','parent_guards'):
                for path in (b/folder).glob('*.json'):
                    if path.stem.split('_')[0] not in failed_ids:continue
                    destination=audit/'old_worker_records'/str(old_bases.index(b))/folder/path.name
                    destination.parent.mkdir(parents=True,exist_ok=True);path.rename(destination)
        report=dict(status='PREPARED_NOT_SUBMITTED',retry_tasks=len(verified),fresh_tasks=len(fresh),groups=groups,
                    total_nodes=sum(g['nodes'] for g in groups),cpus=sum(g['nodes']*32 for g in groups),
                    ranking=scores,ordered=ordered,created_epoch=time.time(),old_scratch_retained=True)
        writej(journal,report);return report


def submit(root):
    root=Path(root).resolve();report=readj(root/AUDIT/'PREPARED.json');site=readj(root/'config/site_scnet.json')
    with queue_lock(root/'operations/submission_lock'):
        jobs=[]
        for i,g in enumerate(sorted(report['groups'],key=lambda g:-g['priority_score'])):
            source=Path(g['root']);base=Path(g['base']);ledger=base/'submission.json'
            if ledger.exists():
                previous=readj(ledger)
                if previous['status'] not in ('SUBMITTED_NOT_YET_VERIFIED_RUNNING','RUNNING_REAL_GOAT_VERIFIED'):raise ValueError('Reconcile submission intent')
                jobs.append(previous);continue
            scope=base/'scope.json';n=g['nodes']
            live=subprocess.check_output(['squeue','-A',site['slurm']['account'],'-r','-h','-o','%T|%C'],text=True,timeout=60).splitlines()
            allocated=sum(int(l.split('|')[1]) for l in live if l.split('|')[0] in ('RUNNING','COMPLETING','CONFIGURING'))
            if allocated+n*32>10000 or len(live)+n>195:raise ValueError('Insufficient account/queue headroom')
            launcher=source/'scripts/scnet_goat_expansion.sbatch'
            cmd=['sbatch','--parsable','--account='+site['slurm']['account'],'--partition='+site['slurm']['partition'],
                 '--qos='+site['slurm']['qos'],'--time=720:00:00','--time-min=720:00:00','--nodes=1',
                 '--array=0-'+str(n-1)+'%'+str(n),'--ntasks-per-node=1','--cpus-per-task=32','--mem=114000M',
                 '--job-name=work2-64-r22-'+str(i),'--chdir='+str(source),'--exclude=a10r2n13,k18r2n[10-17]',
                 '--output='+str(base/'slurm-%A_%a.out'),'--error='+str(base/'slurm-%A_%a.err')]
            env=dict(os.environ,WORK2_PACKAGE_ROOT=str(source),WORK2_PYTHON=site['python'],
                     WORK2_PACKED_PROFILE='balanced' if g['cpus_per_job']==8 else 'large_memory',
                     WORK2_RUNTIME_EVIDENCE=str(root/'operations/runtime_124440267.json'),
                     WORK2_GOAT_BULK_BASE=str(base),WORK2_GOAT_BULK_SCOPE=str(scope),WORK2_GOAT_BULK_SCOPE_SHA256=sha(scope))
            test=subprocess.run(cmd+['--test-only',str(launcher)],env=env,capture_output=True,text=True,check=True,timeout=60)
            record=dict(status='SUBMISSION_INTENT',command=cmd,created_epoch=time.time(),nodes=n,cpus=n*32,scope_sha256=sha(scope),
                        tasks=len(g['jobs']),allocated_cpus_before=allocated,priority_score=g['priority_score'],preflight=test.stderr)
            writej(ledger,record)
            try:
                p=subprocess.run(cmd+[str(launcher)],env=env,capture_output=True,text=True,check=True,timeout=60)
                jid=p.stdout.strip().split(';')[0]
                if not jid.isdigit():raise ValueError('Unexpected submission response')
            except Exception as exc:
                writej(ledger,dict(record,status='SUBMISSION_UNCERTAIN',error=repr(exc)));raise
            record.update(status='SUBMITTED_NOT_YET_VERIFIED_RUNNING',job_id=jid);writej(ledger,record);jobs.append(record)
        return jobs

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','submit']);parser.add_argument('--probe',type=Path)
    args=parser.parse_args();print(json.dumps(prepare(ROOT,args.probe) if args.action=='prepare' else submit(ROOT),ensure_ascii=False,indent=2))
