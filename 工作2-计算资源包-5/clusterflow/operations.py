"""Durable coverage audit and bounded recovery. Never delete scientific data."""
from pathlib import Path
import os
import re
import socket
import subprocess
import time
import uuid
from .common import readj, writej, sha, inside
from .runtime import classify_failure, lock_state, owner_record


def scheduler_terminal(record):
    job = str(record.get('slurm_job_id') or '')
    if record.get('slurm_array_task_id') is not None: job += '_' + str(record['slurm_array_task_id'])
    if not re.fullmatch(r'\d+(?:_\d+)?', job): return False
    try:
        result = subprocess.run(['sacct','-n','-P','-j',job,'--format=JobIDRaw,State'],
                                capture_output=True, text=True, timeout=20, check=False)
    except (OSError, subprocess.TimeoutExpired): return False
    if result.returncode: return False
    states = [p[1].split()[0].split('+')[0] for p in (v.split('|') for v in result.stdout.splitlines())
              if len(p)>1 and p[0] == job]
    return len(states)==1 and states[0] in {'COMPLETED','CANCELLED','FAILED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED'}


def retire_terminal_lock(path):
    """Reconcile ONLY scheduler-confirmed terminal owners, retaining evidence.

    A dead local PID alone is insufficient: an orphan child could still run.
    Caller must hold the outer recovery/queue mutex, never use age as evidence.
    """
    path = Path(path)
    record = readj(path)
    digest = sha(path)
    if not scheduler_terminal(record): return False
    if record.get('hostname') == socket.gethostname() and lock_state(path) == 'RUNNING': return False
    if sha(path) != digest: raise ValueError('Lock owner changed during reconciliation')
    archive = path.parent/'lock_history'
    archive.mkdir(exist_ok=True)
    path.rename(archive/(path.name.lstrip('.')+'_'+uuid.uuid4().hex+'.json'))
    return True


def tail(path, limit=1024*1024):
    path = Path(path)
    if not path.is_file(): return ''
    with path.open('rb') as stream:
        stream.seek(max(0, path.stat().st_size-limit))
        return stream.read().decode('utf-8', errors='replace')


def diagnosis(directory):
    directory = Path(directory)
    evidence = {}
    text = ''
    for name in ('job.out','job.stderr','goat.out','goat.stderr','receipt.json'):
        path = directory/name
        if path.is_file():
            evidence[name] = sha(path)
            text += '\n' + tail(path)
    kind = classify_failure(text)
    if (directory/'receipt.json').is_file():
        receipt = readj(directory/'receipt.json')
        if receipt.get('timed_out') or receipt.get('failure_reason') in ('timeout','keyboard_interrupt'):
            kind = 'INTERRUPTED'
    action = {
        'OOM':'new resource-scoped branch: lower actual MPI/increase per-rank budget within node total; never repeat unchanged',
        'NODE_READ_ONLY':'quarantine node; same scientific input may retry on healthy node, at most two infrastructure attempts',
        'DISK_OR_QUOTA':'stop NEW admissions; verify quota/inodes/free space, recover stage-out before considering QM retry',
        'INTERRUPTED':'preserve last coordinates/GBW/Hessian; inspect restartability; not a promise of exact GOAT resume',
        'CONTAINER_IO':'preserve logs; private sandbox plus healthy node; validate runtime hashes',
        'MPI_OR_SEGFAULT_REVIEW':'distinguish shared-memory/transport/heap/module fault from OOM before retry',
        'MATRIX_IO_REVIEW':'inspect disk/quota/lost matrices and container; do not assume OOM',
        'GOAT_NO_VALID_ENSEMBLE':'do not promote seed or normal termination; inspect identity/geometry then release independent starts',
        'MODULE_FAILURE_REVIEW':'inspect Startup/PROPINT inputs, basis, environment and memory; module name is not root cause',
        'UNKNOWN_REVIEW':'retain evidence and diagnose; no blind repeat'
    }[kind]
    return {'directory': str(directory), 'category': kind, 'action': action,
            'evidence_sha256': evidence, 'automatic_scientific_approval': False}


def coverage(root, write=False):
    root = Path(root).resolve()
    products = {}
    for model in readj(root/'registry/batch_plan_v5.json')['species']:
        if model['compound_id'] == 'TMS': continue
        sid = model['species_id']; base = root/'automation'/sid
        state = readj(base/'state.json') if (base/'state.json').is_file() else {}
        latest = readj(base/'latest_status.json') if (base/'latest_status.json').is_file() else {}
        if (base/'.autoflow.lock').exists():
            status = lock_state(base/'.autoflow.lock')
        elif not state: status = 'NOT_STARTED_READY_FOR_ADMISSION'
        elif state.get('phase') == 'REVIEW': status = 'SCIENTIFIC_REVIEW_REQUIRED'
        elif str(latest.get('status','')).startswith('EXECUTED_'):
            status = 'READY_FOR_NEXT_CONTROLLER_STEP'  # Previous execution has returned, not RUNNING.
        else: status = latest.get('status', 'CONTROLLER_RECONCILIATION_REQUIRED')
        products.setdefault(model['compound_id'], []).append({
            'species_id': sid, 'status': status, 'phase': state.get('phase'),
            'completed_exploratory_waves': sum('prediction' in w for w in state.get('waves', [])),
            'last_controller_update': state.get('updated_utc'),
            'needs_followup': status not in ('RUNNING',),
            'final_delivery_approved': False})
    report = {'timestamp_epoch': time.time(), 'products': products,
              'product_count': len(products), 'model_count': sum(map(len, products.values())),
              'notice': 'Not a spectrum-delivery certificate; no parent COMPLETED or old EXECUTED flag is called active chemistry.'}
    if write: writej(root/'operations/coverage.json', report)
    return report


def retry_dft_infrastructure(root, relative, orca, scratch_root, execute=False):
    """Only explicit, diagnosed node failure with terminal receipt; max two.

    Scientific/unknown/OOM/interruption cases produce repair plans instead of
    an identical expensive restart. Old output is archived by the runner.
    """
    root = Path(root).resolve(); directory = inside(root, relative)
    report = diagnosis(directory)
    if any((directory/name).exists() for name in ('.orca-run.lock','.clusterflow-run.lock')):
        raise ValueError('Owner/stage-out reconciliation required before any retry')
    receipt = readj(directory/'receipt.json')
    if receipt.get('status') != 'failed': raise ValueError('Retry requires a terminal failed DFT receipt')
    if report['category'] != 'NODE_READ_ONLY':
        return {**report, 'status':'TARGETED_REPAIR_REQUIRED_NO_BLIND_RETRY'}
    if receipt.get('hostname') == socket.gethostname():
        return {**report, 'status':'WAIT_HEALTHY_DIFFERENT_NODE'}
    ledger = directory/'recovery'
    attempts = list(ledger.glob('attempt_*.json')) if ledger.is_dir() else []
    if len(attempts) >= 2: return {**report, 'status':'RETRY_LIMIT_REVIEW'}
    if not execute: return {**report, 'status':'DIAGNOSED_RETRY_PLAN_ONLY'}
    ledger.mkdir(exist_ok=True)
    path = ledger/('attempt_'+uuid.uuid4().hex+'.json')
    writej(path, {**report, **owner_record('recovery'), 'status':'INTENT_PERSISTED_BEFORE_LAUNCH'})
    from .runner import run_one
    result = run_one(root, relative, orca, acknowledge=True, retry=True, scratch_root=scratch_root)
    writej(path, {**report, 'status':'ATTEMPT_RETURNED', 'result':result})
    return result


def recover_stageout(root, intent_path, execute=False):
    """Recover a terminal scratch receipt, not restart QM or invent success.

    Requires terminal scheduler owner and matching original/staged inputs.
    Equal partial copies resume; conflicts fail. Receipt last; retain scratch.
    """
    from .runner import _safe_file, _copy_new, _commit_receipt, _lock
    from .runtime import site_config
    root=Path(root).resolve();intent_path=Path(intent_path).resolve()
    if not intent_path.is_relative_to(root) or intent_path.parent.name!='stage_intents':
        raise ValueError('Use a recorded stage intent inside this run copy')
    intent=readj(intent_path)
    original=Path(intent['original']).resolve();stage=Path(intent['stage']).resolve()
    if intent_path.parent.parent!=original or not original.is_relative_to(root):
        raise ValueError('Intent original-directory mismatch')
    scratch=Path(site_config()['shared_scratch']).resolve()
    if not stage.is_relative_to(scratch) or stage==scratch or not scheduler_terminal(intent):
        raise ValueError('Shared scratch/terminal scheduler owner not verified')
    locks=[d/name for d in (original,stage) for name in ('.goat-run.lock','.orca-run.lock','.clusterflow-run.lock') if (d/name).exists()]
    for lock in locks:
        if not scheduler_terminal(readj(lock)):
            raise ValueError('A stage/original owner is still active or unknown')
    def validate():
        for name,digest in intent['inputs'].items():
            if sha(_safe_file(original,name))!=digest or sha(_safe_file(stage,name))!=digest:
                raise ValueError('Stage/original inputs changed')
        receipt=readj(_safe_file(stage,'receipt.json'))
        if receipt.get('job_id',receipt.get('id'))!=intent['job_identity']:
            raise ValueError('Wrong scratch job identity')
        if receipt.get('status') not in ('completed','failed','FAILED','COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW'):
            raise ValueError('No terminal receipt; preserve partial output for restart diagnosis')
        required_inputs=set(intent['inputs'])-({'isotope_map.json'} if intent['kind']=='GOAT' else set())
        if not required_inputs<=set(receipt.get('input_sha256',{})):
            raise ValueError('Incomplete terminal input bindings')
        for name,digest in receipt['input_sha256'].items():
            if intent['inputs'].get(name)!=digest:raise ValueError('Receipt input mismatch')
        outputs=receipt.get('output_sha256')
        required={'goat.out'} if intent['kind']=='GOAT' else {'job.out'}
        if intent['kind']=='GOAT' and receipt['status'] in ('completed','COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW'):
            required.add('goat.finalensemble.xyz')
        if not isinstance(outputs,dict) or not required<=outputs.keys():raise ValueError('Incomplete terminal output bindings')
        for name,digest in outputs.items():
            if name in intent['inputs'] or name=='receipt.json' or name.startswith('.'):
                raise ValueError('Unsafe output inventory')
            if sha(_safe_file(stage,name))!=digest:raise ValueError('Scratch output mismatch')
            target=_safe_file(original,name,existing=False)
            if target.exists() and sha(target)!=digest:raise ValueError('Conflicting partial stage-out; no overwrite')
        return receipt
    receipt=validate()
    result={'status':'STAGEOUT_RECOVERY_PLAN','intent':str(intent_path),'receipt_status':receipt['status'],
            'quantum_jobs_launched':0,'scratch_retained':True}
    if not execute:return result
    with _lock(original,'.stageout-recovery.lock'):
        for lock in locks:
            if lock.exists() and not retire_terminal_lock(lock):raise ValueError('Lock reconciliation failed')
        receipt=validate();digest=sha(stage/'receipt.json')
        for name,expected in receipt['output_sha256'].items():
            target=_safe_file(original,name,existing=False)
            if not target.exists():_copy_new(_safe_file(stage,name),target,expected)
        if (original/'receipt.json').exists():
            if sha(original/'receipt.json')!=digest:raise ValueError('Original receipt differs; no overwrite')
        else:_commit_receipt(stage/'receipt.json',original/'receipt.json',digest,validate)
        writej(original/'stageout_recovery.json',{**result,'status':'RECOVERED_EXISTING_TERMINAL_RECEIPT',
               'receipt_sha256':digest,'intent_sha256':sha(intent_path)})
    return {**result,'status':'RECOVERED_EXISTING_TERMINAL_RECEIPT'}


def cleanup_plan(root, write=False):
    """Only proven redundant SUCCESSFUL scratch; generates plan, deletes nothing.

    Execution sources/cache/GBW/ensemble in the original run tree stay pinned.
    Failed, live, unknown-owner or incompletely returned attempts are excluded.
    """
    from .runtime import site_config
    from .gates import inspect_job
    from .pipeline import _goat_source
    root=Path(root).resolve();scratch=Path(site_config()['shared_scratch']).resolve()
    candidates=[];excluded=[];seen=set()
    for intent_path in root.glob('**/stage_intents/*.json'):
        intent=readj(intent_path);stage=Path(intent['stage']).resolve();original=Path(intent['original']).resolve()
        if str(stage) in seen:continue
        seen.add(str(stage))
        try:
            if not stage.is_relative_to(scratch) or stage==scratch or not original.is_relative_to(root):raise ValueError('scope')
            if not scheduler_terminal(intent):raise ValueError('live/unknown scheduler owner')
            if any((d/n).exists() for d in (stage,original) for n in ('.goat-run.lock','.orca-run.lock','.clusterflow-run.lock','.stageout-recovery.lock')):raise ValueError('lock present')
            receipt=readj(original/'receipt.json')
            if receipt.get('status') not in ('completed','COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW'):raise ValueError('not successful')
            if sha(stage/'receipt.json')!=sha(original/'receipt.json'):raise ValueError('receipt differs')
            for name,digest in intent['inputs'].items():
                if sha(original/name)!=digest or sha(stage/name)!=digest:raise ValueError('input changed')
            for name,digest in receipt['output_sha256'].items():
                if sha(original/name)!=digest or sha(stage/name)!=digest:raise ValueError('incomplete stage-out')
            workspace=next(p for p in original.parents if (p/'config/protocol.json').is_file())
            if intent['kind']=='DFT':
                if not inspect_job(original,readj(workspace/'config/protocol.json'))['accepted']:raise ValueError('scientific gate incomplete')
            else:_goat_source(workspace,intent['job_identity'])
            files=list(stage.rglob('*'))
            if any(p.is_symlink() for p in files):raise ValueError('linked scratch tree')
            candidates.append({'redundant_scratch':str(stage),'retain_original':str(original),
                'intent_sha256':sha(intent_path),'receipt_sha256':sha(original/'receipt.json'),
                'bytes':sum(p.stat().st_size for p in files if p.is_file()),
                'action':'RECHECK_THEN_SEPARATELY_AUTHORIZE_CLEANUP_NO_DELETION_PERFORMED'})
        except (ValueError,OSError,KeyError,StopIteration) as exc:
            excluded.append({'stage':str(stage),'reason':str(exc)})
    report={'status':'AUDIT_ONLY_NO_DELETION','candidates':candidates,'excluded':excluded,
            'candidate_bytes':sum(r['bytes'] for r in candidates)}
    if write:writej(root/'operations/cleanup_plan.json',report)
    return report
