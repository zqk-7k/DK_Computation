"""Submit one explicitly authorized 80-node GOAT-only pool, idempotently."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from clusterflow.common import readj, writej, sha
from clusterflow.runtime import site_config, check_storage
from clusterflow.priority import estimate
from scripts.packed_worker import queue_lock

IDENTITIES = ('registry/jobs.json', 'registry/batch_plan_v5.json', 'registry/species.json',
              'config/protocol.json', 'config/automation.json', 'config/site_scnet.json')


def build_scope(root):
    plans = readj(root/'registry/batch_plan_v5.json')['species']
    jobs = {r['job_id']: r for r in readj(root/'registry/jobs.json')}
    # Long remaining workloads start first; one shared reference is not duplicated.
    ranked = sorted(plans, key=lambda r: (r['compound_id'] != 'TMS',
                    -estimate(root, r, jobs)['remaining_serial_hours_estimate'], r['species_id']))
    selected = []; excluded = []; legacy = []
    for row in ranked:
        sid = row['species_id']
        blocked = (root/'operations/packed/blocked'/(sid+'.json')).exists()
        for job in row['goat_start_ids'][:row['initial_goat_starts']]:
            entry = jobs[job]
            if entry.get('initial_release') is not True:
                raise ValueError('Initial start was not released in the immutable registry')
            directory = root/entry['job_directory']
            if (directory/'.goat-run.lock').exists():
                legacy.append({'job_id': job, **readj(directory/'.goat-run.lock')})
                excluded.append({'job_id': job, 'reason': 'EXISTING_GOAT_OWNER'})
            elif blocked or any((directory/name).exists() for name in ('receipt.json', 'portable_stageout.json')) or any((directory/'stage_intents').glob('*.json')):
                excluded.append({'job_id': job, 'reason': 'BLOCKED_OR_EXISTING_ATTEMPT'})
            else:
                selected.append(job)
    return {'schema': 1, 'created_epoch': time.time(), 'stage': 'goat',
            'operations_directory': 'operations/goat_bulk_20261015',
            'job_ids': selected, 'excluded': excluded, 'legacy_owners': legacy,
            'legacy_active_jobs': len(legacy), 'cpus_per_job': 8, 'slots_per_node': 4,
            'growth_gib_per_job': 8, 'free_floor_gib': 128,
            'storage_basis': 'GOAT-only reservation: current 12 scratch samples 0.87-2.46 GiB; reserve both scratch and stage-out plus growth. DFT 64-GiB reservation is unchanged. Recheck free space under every claim lock.',
            'identity_sha256': {name: sha(root/name) for name in IDENTITIES},
            'notice': 'Initial released starts only; no retry, held starts, DFT or automatic scientific acceptance.'}


def submit(root, execute=False):
    root = Path(root).resolve(); site = site_config(); slurm = site['slurm']
    base = root/'operations/goat_bulk_20261015'; base.mkdir(parents=True, exist_ok=True)
    with queue_lock(root/'operations/submission_lock'):
        journal = base/'submission.json'
        if journal.exists():
            old = readj(journal)
            if old.get('status') != 'SUBMITTED_NOT_RUNNING' or not old.get('job_id'):
                raise ValueError('Reconcile existing submission intent; do not resubmit')
            live = subprocess.run(['squeue', '-h', '-r', '-j', old['job_id'], '-o', '%i|%T|%C|%D|%R'],
                                  capture_output=True, text=True, timeout=30, check=True)
            return {**old, 'scheduler_queue': live.stdout}
        if not all((root/name).exists() for name in ('operations/ROLLOUT_PAUSED', 'operations/packed/STOP_NEW_CLAIMS')):
            raise ValueError('Existing controllers must drain before start-level pooling')
        runtime = root/'operations/runtime_124440267.json'
        evidence = readj(runtime)
        if evidence.get('version') != '6.1.0' or evidence.get('site_sha256') != sha(root/'config/site_scnet.json'):
            raise ValueError('Canonical runtime evidence mismatch')
        scope = build_scope(root)
        if len(scope['job_ids']) < 320:
            raise ValueError('Insufficient unclaimed starts for 80 nodes')
        reserve = 128 + (320 + scope['legacy_active_jobs']) * scope['growth_gib_per_job']
        storage = check_storage(site['shared_scratch'], reserve)
        excluded = sorted(p.stem for p in (root/'operations/packed/quarantine').glob('*.json'))
        if any(not re.fullmatch(r'[A-Za-z0-9_.-]+', name) for name in excluded):
            raise ValueError('Invalid node exclusion')
        command = ['sbatch', '--parsable', '--account='+slurm['account'], '--partition='+slurm['partition'],
                   '--qos='+slurm['qos'], '--nodes=8', '--array=0-9%10', '--time=720:00:00', '--ntasks-per-node=1',
                   '--cpus-per-task=32', '--mem=114000M', '--job-name=work2-64-goat2560',
                   '--chdir='+str(root), '--output='+str(base/'slurm-%A_%a.out'), '--error='+str(base/'slurm-%A_%a.err')]
        if excluded:
            command.append('--exclude='+','.join(excluded))
        launcher = root/'scripts/scnet_goat_bulk.sbatch'
        tested = subprocess.run(command+['--test-only', str(launcher)], capture_output=True, text=True, timeout=60, check=True)
        queued = subprocess.run(['squeue', '-h', '-u', 'qikai', '-o', '%i'], capture_output=True, text=True, timeout=30, check=True)
        if len(queued.stdout.splitlines()) >= 190:
            raise ValueError('Insufficient queue headroom')
        scope_path = base/'scope.json'
        writej(scope_path, scope)
        report = {'status': 'DRY_RUN_NO_SUBMISSION', 'nodes': 80, 'cpus': 2560, 'slots': 320,
                  'array_tasks': 10, 'nodes_per_task': 8,
                  'walltime_hours': 720, 'pending_goat_starts': len(scope['job_ids']),
                  'legacy_active_jobs': scope['legacy_active_jobs'], 'excluded_nodes': excluded,
                  'scope_sha256': sha(scope_path), 'storage_reserve_gib': reserve, 'storage': storage,
                  'test_only': {'stdout': tested.stdout, 'stderr': tested.stderr},
                  'created_epoch': time.time(), 'scientific_inputs_changed': False}
        writej(base/'preflight.json', report)
        if not execute:
            return report
        env = os.environ.copy()
        env.update({'WORK2_PACKAGE_ROOT': str(root), 'WORK2_PYTHON': site['python'],
                    'WORK2_PACKED_PROFILE': 'balanced', 'WORK2_RUNTIME_EVIDENCE': str(runtime),
                    'WORK2_GOAT_BULK_BASE': str(base), 'WORK2_GOAT_BULK_SCOPE': str(scope_path),
                    'WORK2_GOAT_BULK_SCOPE_SHA256': sha(scope_path)})
        writej(journal, {**report, 'status': 'SUBMISSION_INTENT'})
        try:
            result = subprocess.run(command+[str(launcher)], env=env, capture_output=True, text=True, timeout=60, check=True)
            job = result.stdout.strip().split(';')[0]
            if not job.isdigit():
                raise ValueError('Unrecognized sbatch output')
        except Exception as exc:
            writej(journal, {**report, 'status': 'SUBMISSION_UNCERTAIN', 'reason': str(exc)})
            raise
        report.update({'status': 'SUBMITTED_NOT_RUNNING', 'job_id': job})
        writej(journal, report)
        return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execute', action='store_true')
    print(json.dumps(submit(ROOT, execute=p.parse_args().execute), ensure_ascii=False, indent=2))
