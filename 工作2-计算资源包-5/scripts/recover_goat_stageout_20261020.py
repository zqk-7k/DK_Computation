"""Recover hash-bound terminal GOAT output after a full-disk stage-out failure.

The per-job execution flock must be free. A terminal receipt is mandatory;
RUNNING receipts and orphan quantum executions are never called successful.
Conflicting partial copies are retained in the audit directory.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import uuid

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from clusterflow.common import readj, writej, sha, inside
from clusterflow.scheduling import single_execution
from clusterflow.runner import _safe_file, _copy_new, _commit_receipt, _lock
from scripts.goat_portable import IMMUTABLE_FILES
from scripts.packed_worker import queue_lock


def validate_terminal(root, job, intent):
    indexed={r['job_id']:r for r in readj(root/'registry/jobs.json')}
    original=inside(root,indexed[job]['job_directory'])
    stage=Path(intent['stage']).resolve()
    site=readj(root/'config/site_scnet.json')
    if str(original)!=intent['original'] or intent['job_identity']!=job:
        raise ValueError('Stage intent identity mismatch')
    if not stage.is_relative_to(Path(site['shared_scratch']).resolve()):
        raise ValueError('Stage escapes the dedicated scratch root')
    for name in IMMUTABLE_FILES:
        if sha(_safe_file(original,name))!=intent['inputs'][name] or sha(_safe_file(stage,name))!=intent['inputs'][name]:
            raise ValueError('Stage or canonical input changed')
    receipt=readj(_safe_file(stage,'receipt.json'))
    if receipt.get('job_id')!=job or receipt.get('status') not in ('FAILED','COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW'):
        raise ValueError('No actual terminal GOAT receipt')
    for name in set(IMMUTABLE_FILES)-{'isotope_map.json'}:
        if receipt.get('input_sha256',{}).get(name)!=intent['inputs'][name]:
            raise ValueError('Terminal input binding missing or different')
    outputs=receipt.get('output_sha256',{})
    required={'goat.out'}
    if receipt['status']=='COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW':
        required.add('goat.finalensemble.xyz')
    if not required<=set(outputs):raise ValueError('Missing terminal output bindings')
    for name,digest in outputs.items():
        if name in IMMUTABLE_FILES or not name.startswith('goat.') or name.startswith('.'):
            raise ValueError('Unsafe terminal output inventory')
        if sha(_safe_file(stage,name))!=digest:
            raise ValueError('Scratch output changed: '+name)
    return original,stage,receipt


def recover(root,job):
    root=Path(root).resolve()
    indexed={r['job_id']:r for r in readj(root/'registry/jobs.json')}
    original=inside(root,indexed[job]['job_directory'])
    intents=list((original/'stage_intents').glob('*.json'))
    if len(intents)!=1:raise ValueError('Exactly one unrecovered stage intent required')
    intent=readj(intents[0])
    key='goat_'+hashlib.sha256(job.encode()).hexdigest()
    with single_execution(root,key), _lock(original,'.stageout-recovery.lock'):
        original,stage,receipt=validate_terminal(root,job,intent)
        stamp=uuid.uuid4().hex
        archive=original/'recovery_history'/('stageout_20261020_'+stamp)
        archive.mkdir(parents=True)
        writej(archive/'INTENT.json',dict(intent=intent,receipt_sha256=sha(stage/'receipt.json'),
                                         reason='Terminal output verified and producer execution flock released'))
        # A producer commits its terminal receipt after wait() has returned.
        # Together with the execution flock, this is stronger than the state of
        # an unrelated still-running parent allocation.
        for directory in (original,stage):
            lock=directory/'.goat-run.lock'
            if lock.exists():
                owner=readj(lock)
                if owner.get('kind')=='RESOURCE_COHORT_DELEGATION':
                    raise ValueError('Job delegated to another execution')
                lock.rename(archive/('original_lock.json' if directory==original else 'stage_lock.json'))
        copies=[];conflicts=[]
        for name,digest in receipt['output_sha256'].items():
            src=_safe_file(stage,name);dst=_safe_file(original,name,existing=False)
            if dst.exists() and sha(dst)!=digest:
                old_hash=sha(dst)
                saved=archive/'partial_outputs'/name;saved.parent.mkdir(parents=True,exist_ok=True)
                dst.rename(saved)
                conflicts.append(dict(name=name,retained_path=str(saved),old_sha256=old_hash))
            if not dst.exists():
                _copy_new(src,dst,digest);copies.append(name)
        digest=sha(stage/'receipt.json')
        if (original/'receipt.json').exists():
            if sha(original/'receipt.json')!=digest:
                raise ValueError('A different original receipt exists; explicit review required')
        else:
            _commit_receipt(stage/'receipt.json',original/'receipt.json',digest,
                            lambda: validate_terminal(root,job,intent))
        result=dict(status='RECOVERED_EXISTING_TERMINAL_RECEIPT',job_id=job,
                    receipt_status=receipt['status'],receipt_sha256=digest,quantum_jobs_launched=0,
                    copied_files=len(copies),retained_conflicting_partials=conflicts,
                    audit=str(archive),at_epoch=time.time())
        writej(original/'stageout_recovery.json',result)
        writej(archive/'COMPLETE.json',result)
    base=root/'operations/goat_bulk_20261015'
    with queue_lock(base):
        queue=readj(base/'queue.json')
        old=queue['active'].pop(job,None)
        outcome_path=base/'results'/(job+'.json')
        if outcome_path.exists():
            prior=readj(outcome_path);writej(archive/'previous_pool_result.json',prior)
        else:prior=old or {'job_id':job,'species_id':indexed[job]['species_id']}
        if old:writej(archive/'previous_active_claim.json',old)
        writej(outcome_path,dict(prior,result=dict(status=receipt['status'],recovered_stageout=result),
                                finished_epoch=time.time()))
        claim=base/'claims'/(job+'.json')
        if claim.exists():claim.rename(archive/'pool_claim.json')
        writej(base/'queue.json',queue)
        block=base/'blocked'/(indexed[job]['species_id']+'.json')
        if receipt['status']=='COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW' and block.exists() and readj(block).get('job_id')==job:
            block.rename(archive/'retired_branch_block.json')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--job',required=True)
    args=p.parse_args()
    print(json.dumps(recover(ROOT,args.job),ensure_ascii=False,indent=2))
