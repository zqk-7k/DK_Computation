"""Read-only GOAT pool audit. Running processes are not accepted chemistry."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from collections import Counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from clusterflow.common import readj as _readj, writej
from clusterflow.runtime import site_config


def readj(path):
    # ParaStor can briefly lose an old file handle during a remote atomic
    # replacement. Retry observation only; never infer completion from absence.
    for attempt in range(10):
        try:
            return _readj(path)
        except (FileNotFoundError, json.JSONDecodeError):
            if attempt == 9:
                raise
            time.sleep(.2)


def audit(root):
    base = root/'operations/goat_bulk_20261015'
    scope = readj(base/'scope.json')
    jobs = {r['job_id']: r for r in readj(root/'registry/jobs.json')}
    queue = readj(base/'queue.json') if (base/'queue.json').exists() else {'active': {}, 'cursor': 0, 'skipped': []}
    live = []; stale = []; results = Counter(); errors = []
    for job, record in queue['active'].items():
        lock = root/jobs[job]['job_directory']/'.goat-run.lock'
        current = False
        if lock.exists():
            owner = readj(lock)
            stage = Path(owner.get('portable_stage_root', '/nonexistent'))/'job'
            heartbeat = stage/'runtime_heartbeat.json'
            if heartbeat.exists() and time.time()-readj(heartbeat).get('updated_epoch', 0) < 180:
                current = True
        (live if current else stale).append(record)
    for path in (base/'results').glob('*.json'):
        row = readj(path)
        status = row['result']['status']; results[status] += 1
        if status not in ('COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW', 'completed',
                          'TECHNICALLY_RECOVERED_AWAITING_SCIENTIFIC_REVIEW',
                          'NODE_FAIL_DELEGATED_TO_RESOURCE_COHORT', 'BOUNDED_RECOVERY_DELEGATED'):
            errors.append({'job_id': row['job_id'], 'result': row['result']})
    nodes = [readj(p) for p in (base/'nodes').glob('*.json')]
    healthy = [r for r in nodes if r.get('status') == 'RUNNING' and time.time()-r['updated_epoch'] < 90]
    rejected = [r for r in nodes if r.get('status') == 'NODE_REJECTED_BEFORE_QM']
    guards = [readj(p) for p in (base/'parent_guards').glob('*.json')]
    command = subprocess.run(['squeue', '-h', '-r', '-u', 'qikai', '-o', '%i|%j|%T|%C|%D|%R'],
                             capture_output=True, text=True, timeout=30, check=True)
    slurm = [line for line in command.stdout.splitlines() if '|work2-64-' in line]
    submission = readj(base/'submission.json') if (base/'submission.json').exists() else {}
    report = {'updated_epoch': time.time(), 'status': 'GOAT_ONLY_POOL_AUDIT',
              'requested_new_cpus': submission.get('cpus'), 'queue': slurm, 'frozen_pending_starts': len(scope['job_ids']),
              'cursor': queue['cursor'], 'unclaimed_starts': len(scope['job_ids'])-queue['cursor'],
              'active_claims': len(queue['active']), 'fresh_runtime_heartbeats': len(live),
              'fresh_runtime_configured_cpus': len(live)*scope['cpus_per_job'], 'claims_awaiting_start_or_review': len(stale),
              'active_species': len({r['species_id'] for r in live}),
              'nodes_with_fresh_quantum_heartbeats': len({r['hostname'] for r in live}),
              'healthy_worker_nodes': len(healthy), 'live_worker_slots': sum(r['live_slots'] for r in healthy),
              'parent_guards': [{'job_id': r.get('slurm_job_id'), 'status': r['status'],
                                 'age_seconds': round(time.time()-r.get('updated_epoch', 0))} for r in guards],
              'rejected_nodes': rejected, 'results': dict(results), 'failures': errors[:20],
              'shared_free_gib': round(shutil.disk_usage(site_config()['shared_scratch']).free/1024**3, 2),
              'old_branch_controllers_draining': (root/'operations/packed/STOP_NEW_CLAIMS').exists(),
              'notice': 'Counts exclude original parent computations; heartbeats do not imply 100% CPU utilization or chemical acceptance.'}
    writej(base/'latest.json', report)
    return report


if __name__ == '__main__':
    print(json.dumps(audit(ROOT), ensure_ascii=False, indent=2))
