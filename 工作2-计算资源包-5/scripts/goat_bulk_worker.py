"""Finite GOAT-only start-level pool; immutable inputs and original receipts.

Unlike the branch controller, independent released starts may run concurrently.
The branch controllers must be draining. Unknown ownership is never reclaimed.
No DFT, retry, scientific acceptance or scratch deletion occurs here.
"""
import argparse
import fcntl
import hashlib
import json
import multiprocessing
from datetime import datetime
import os
from pathlib import Path
import socket
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from clusterflow.common import readj, writej, sha, inside
from clusterflow.runtime import site_config, owner_record, check_storage, StorageAdmissionError
from clusterflow.scheduling import single_execution, admit_quantum, record_timing, DeferredWork
from scripts.packed_worker import queue_lock


def allocation_deadline(base):
    """One scheduler query per parent, shared by all its compute-node workers."""
    job = os.environ['SLURM_JOB_ID']
    if not job.isdigit():
        raise ValueError('Numeric Slurm allocation required')
    directory = base/'allocations'/job
    with queue_lock(directory):
        path = directory/'deadline.json'
        if path.exists():
            record = readj(path)
            if record.get('job_id') != job or record['deadline_epoch'] <= time.time():
                raise ValueError('Expired or mismatched allocation deadline')
            return record['deadline_epoch']
        last = None
        for _ in range(3):
            try:
                result = subprocess.run(['scontrol', 'show', 'job', job, '-o'],
                                        capture_output=True, text=True, timeout=60, check=True)
                fields = dict(v.split('=', 1) for v in result.stdout.split() if '=' in v)
                if fields.get('JobId') != job or fields.get('JobState') != 'RUNNING':
                    raise ValueError('Allocation is not the requested running parent')
                deadline = datetime.fromisoformat(fields['EndTime']).timestamp()
                if deadline <= time.time():
                    raise ValueError('Allocation deadline already passed')
                writej(path, {'job_id': job, 'deadline_epoch': deadline,
                              'array_job_id': fields.get('ArrayJobId'), 'array_task_id': fields.get('ArrayTaskId'),
                              'scheduler_fields': fields, 'checked_epoch': time.time()})
                return deadline
            except (subprocess.SubprocessError, OSError) as exc:
                last = exc
                time.sleep(2)
        raise RuntimeError('Unable to verify parent deadline: '+str(last))


def verify_scope(root, scope, expected):
    if sha(scope) != expected:
        raise ValueError('Frozen GOAT scope changed')
    policy = readj(scope)
    for name, expected_hash in policy['identity_sha256'].items():
        if sha(inside(root, name)) != expected_hash:
            raise ValueError('GOAT source identity changed: ' + name)
    indexed = {row['job_id']: row for row in readj(root/'registry/jobs.json')}
    plans = readj(root/'registry/batch_plan_v5.json')['species']
    allowed = {j for row in plans for j in row['goat_start_ids'][:row['initial_goat_starts']]}
    ids = policy['job_ids']
    if not ids or len(set(ids)) != len(ids):
        raise ValueError('Empty or duplicated GOAT scope')
    for job in ids:
        if job not in allowed or indexed[job].get('initial_release') is not True:
            raise ValueError('Scope contains an unreleased start')
        inside(root, indexed[job]['job_directory'])
    if policy['cpus_per_job'] != 8 or policy['slots_per_node'] != 4:
        raise ValueError('Only verified 4 x 8 CPU layout is supported')
    if policy['growth_gib_per_job'] < 8 or policy['free_floor_gib'] < 128:
        raise ValueError('GOAT storage reservation below reviewed floor')
    return policy, indexed


def claim(root, base, policy, indexed, storage_check=check_storage):
    with queue_lock(base):
        if (base/'PAUSE_NEW_CLAIMS').exists():
            return None
        state_path = base/'queue.json'
        state = readj(state_path) if state_path.exists() else {'cursor': 0, 'active': {}, 'skipped': []}
        # Existing branch controllers finish only their current GOAT starts.
        # Reserve their full growth allowance even after an individual finishes.
        reserve = policy['free_floor_gib'] + policy['growth_gib_per_job'] * (
            len(state['active']) + policy['legacy_active_jobs'] + 1)
        storage_check(site_config()['shared_scratch'], reserve)
        while state['cursor'] < len(policy['job_ids']):
            job = policy['job_ids'][state['cursor']]
            state['cursor'] += 1
            row = indexed[job]; sid = row['species_id']
            directory = inside(root, row['job_directory'])
            markers = ['receipt.json', '.goat-run.lock', 'portable_stageout.json']
            occupied = any((directory/name).exists() for name in markers)
            occupied = occupied or any((directory/'stage_intents').glob('*.json'))
            occupied = occupied or (base/'results'/(job+'.json')).exists()
            blocked = (root/'operations/packed/blocked'/(sid+'.json')).exists()
            blocked = blocked or (base/'blocked'/(sid+'.json')).exists()
            if occupied or blocked:
                state['skipped'].append({'job_id': job, 'reason': 'EXISTING_ATTEMPT_OR_BLOCKED_SPECIES'})
                continue
            record = {**owner_record(uuid.uuid4().hex), 'job_id': job, 'species_id': sid,
                      'growth_gib': policy['growth_gib_per_job']}
            state['active'][job] = record
            writej(state_path, state)
            writej(base/'claims'/(job+'.json'), record)
            return record
        writej(state_path, state)
        return None


def finish(root, base, record, result):
    with queue_lock(base):
        state_path = base/'queue.json'
        state = readj(state_path)
        if state['active'].get(record['job_id'], {}).get('token') != record['token']:
            raise ValueError('GOAT claim ownership changed')
        outcome = {**record, 'finished_epoch': time.time(), 'result': result}
        writej(base/'results'/(record['job_id']+'.json'), outcome)
        if result.get('status') not in {'COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW', 'completed'}:
            writej(base/'blocked'/(record['species_id']+'.json'), outcome)
        del state['active'][record['job_id']]
        writej(state_path, state)
        path = base/'claims'/(record['job_id']+'.json')
        if path.exists() and readj(path).get('token') == record['token']:
            history = base/'claim_history'; history.mkdir(exist_ok=True)
            path.rename(history/(record['job_id']+'_'+record['token']+'.json'))


def worker(slot, root, base, policy, indexed, deadline):
    site = site_config()
    available = sorted(os.sched_getaffinity(0))
    if len(available) < 32:
        raise ValueError('Full 32 CPU node allocation required')
    os.sched_setaffinity(0, available[slot*8:(slot+1)*8])
    os.environ.update({'SLURM_NTASKS': '1', 'SLURM_CPUS_PER_TASK': '8',
                       'WORK2_SLOT_MEMORY_MB': '28500', 'WORK2_DEADLINE_EPOCH': str(deadline)})
    from scripts.goat_portable import run_one
    while time.time() < deadline - 3600:
        if (base/'PAUSE_NEW_CLAIMS').exists() or (base/('DRAIN_'+os.environ['SLURM_JOB_ID'])).exists():
            break
        try:
            check_storage(site['local_runtime_root'], site['disk']['local_admit_gib'])
            check_storage('/dev/shm', site['disk']['shm_min_gib'], minimum_inodes=100)
            record = claim(root, base, policy, indexed)
        except StorageAdmissionError as exc:
            writej(base/'admission'/(socket.gethostname()+'_'+str(slot)+'.json'),
                   {'status': 'WAITING_STORAGE', 'reason': str(exc), 'time': time.time()})
            # Release this finite worker; never leave a whole empty node polling forever.
            break
        if record is None:
            break
        job = record['job_id']; row = indexed[job]
        result = None
        started = time.monotonic()
        try:
            digest = 'goat_' + hashlib.sha256(job.encode()).hexdigest()
            with single_execution(root, digest):
                admit_quantum(root, 'goat', row['species_id'])
                result = run_one(root, job, Path(site['orca_wrapper']), Path(site['xtb']),
                                 Path(site['shared_scratch']), acknowledge_models=True)
            if result['status'] == 'COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW':
                record_timing(root, 'goat', row['species_id'], time.monotonic()-started)
        except Exception as exc:
            result = {'status': 'REVIEW_REQUIRED', 'reason': type(exc).__name__+': '+str(exc),
                      'notice': 'Inspect retained attempt; never retry or fabricate completion.'}
        finally:
            if result is not None:
                finish(root, base, record, result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scope', type=Path, required=True)
    parser.add_argument('--scope-sha256', required=True)
    parser.add_argument('--recover-startup', action='store_true',
                        help='Only launch on a node with no prior worker record and no surviving pool process')
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('Compute allocation required')
    policy, indexed = verify_scope(ROOT, args.scope, args.scope_sha256)
    base = inside(ROOT, policy['operations_directory'])
    site = site_config()
    node_record = base/'nodes'/(os.environ['SLURM_JOB_ID']+'_'+socket.gethostname()+'.json')
    guard_path = base/'node_guards'/(socket.gethostname()+'.lock')
    guard_path.parent.mkdir(parents=True, exist_ok=True)
    guard = guard_path.open('a')
    fcntl.flock(guard, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if args.recover_startup:
        if node_record.exists():
            raise ValueError('Existing node worker record; startup recovery is not authorized to replace it')
        for entry in Path('/proc').iterdir():
            if not entry.name.isdigit() or int(entry.name) == os.getpid():
                continue
            try:
                command = (entry/'cmdline').read_bytes().split(b'\0')
            except (FileNotFoundError, PermissionError, ProcessLookupError):
                continue
            if any(arg.endswith(b'/goat_bulk_worker.py') or arg == b'goat_bulk_worker.py' for arg in command):
                raise ValueError('A pool worker process already exists on this node')
    try:
        if int(os.environ.get('SLURM_MEM_PER_NODE', '0')) < 114000:
            raise ValueError('Node memory allocation below four verified slots')
        if (ROOT/'operations/packed/quarantine'/(socket.gethostname()+'.json')).exists():
            raise ValueError('Node quarantined')
        check_storage(site['local_runtime_root'], site['disk']['local_admit_gib'])
        check_storage('/dev/shm', site['disk']['shm_min_gib'], minimum_inodes=100)
        from scripts.check_packed_mpi import main as mpi_probe
        mpi_probe()
    except Exception as exc:
        writej(node_record, {'status': 'NODE_REJECTED_BEFORE_QM', 'reason': str(exc),
                            'hostname': socket.gethostname(), 'updated_epoch': time.time()})
        return
    deadline = allocation_deadline(base)
    children = [multiprocessing.Process(target=worker, args=(i, ROOT, base, policy, indexed, deadline))
                for i in range(4)]
    for child in children:
        child.start()
    while True:
        states = [{'pid': p.pid, 'alive': p.is_alive(), 'exitcode': p.exitcode} for p in children]
        live = sum(s['alive'] for s in states)
        writej(node_record, {'status': 'RUNNING' if live else 'EXITED', 'workers': states,
                            'live_slots': live, 'hostname': socket.gethostname(),
                            'slurm_job_id': os.environ['SLURM_JOB_ID'], 'updated_epoch': time.time()})
        if not live:
            break
        time.sleep(15)
    for child in children:
        child.join()
    if any(p.exitcode != 0 for p in children):
        raise RuntimeError('Worker exited unexpectedly; retained claims require review')


if __name__ == '__main__':
    main()
