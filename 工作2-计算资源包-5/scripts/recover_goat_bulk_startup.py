"""Restart only nodes whose pool driver exited BEFORE starting any QM.

Uses already allocated nodes and a guarded --overlap step; never cancels the
original step or job. The compute-side guard rejects any existing worker record
or surviving pool process. No quantum task is retried.
"""
import argparse
import json
import os
from pathlib import Path
import shlex
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from clusterflow.common import readj, writej, sha
from clusterflow.runtime import site_config
from scripts.packed_worker import queue_lock
from scripts.goat_bulk_worker import allocation_deadline


def recover(root, execute=False):
    base = root/'operations/goat_bulk_20261015'
    submitted = readj(base/'submission.json'); array = submitted['job_id']
    scope = base/'scope.json'; site = site_config()
    if sha(scope) != submitted['scope_sha256']:
        raise ValueError('Submitted scope changed')
    queue = subprocess.run(['squeue', '-h', '-r', '-j', array, '-o', '%i|%T'],
                           capture_output=True, text=True, timeout=60, check=True)
    live = {p[0] for p in (line.split('|') for line in queue.stdout.splitlines()) if p[1] == 'RUNNING'}
    reports = []
    for task in range(submitted['array_tasks']):
        if array+'_'+str(task) not in live:
            continue
        error_log = base/('slurm-'+array+'_'+str(task)+'.err')
        if not error_log.exists() or 'end_time' not in error_log.read_text():
            continue
        # Resolve each actual numeric allocation once, not once per node.
        result = subprocess.run(['scontrol', 'show', 'job', array+'_'+str(task), '-o'],
                                capture_output=True, text=True, timeout=60, check=True)
        fields = dict(v.split('=', 1) for v in result.stdout.split() if '=' in v)
        if fields.get('ArrayJobId') != array or fields.get('ArrayTaskId') != str(task) or fields.get('JobState') != 'RUNNING':
            raise ValueError('Array allocation identity changed')
        job = fields['JobId']; candidates = []
        for path in (root/'operations/mpi_probes').glob(job+'_*.json'):
            probe = readj(path); hostname = probe['hostname']
            if probe['status'] == 'PASS' and not (base/'nodes'/(job+'_'+hostname+'.json')).exists():
                candidates.append(hostname)
        if not candidates:
            continue
        journal = base/'startup_recovery'/(job+'.json')
        with queue_lock(base/'startup_recovery'):
            if journal.exists():
                reports.append({'job_id': job, 'status': 'EXISTING_RECOVERY_RECORD', 'record': readj(journal)})
                continue
            record = {'job_id': job, 'array_job_id': array, 'array_task_id': task,
                      'nodes': sorted(candidates), 'created_epoch': time.time(), 'status': 'DRY_RUN'}
            if not execute:
                reports.append(record)
                continue
            # This scheduler identity was checked above; subsequent workers read
            # the same cached EndTime under a per-allocation lock.
            from datetime import datetime
            writej(base/'allocations'/job/'deadline.json', {
                'job_id': job, 'deadline_epoch': datetime.fromisoformat(fields['EndTime']).timestamp(),
                'array_job_id': array, 'array_task_id': str(task),
                'scheduler_fields': fields, 'checked_epoch': time.time()})
            env = os.environ.copy()
            env.update({'WORK2_PACKAGE_ROOT': str(root), 'WORK2_RUNTIME_EVIDENCE': str(root/'operations/runtime_124440267.json'),
                        'WORK2_PACKED_PROFILE': 'balanced', 'WORK2_STAGED_RELEASE': '1',
                        'SLURM_ARRAY_JOB_ID': array, 'SLURM_ARRAY_TASK_ID': str(task)})
            memory = fields.get('MinMemoryNode', '')
            if not re.fullmatch(r'\d+M', memory) or int(memory[:-1]) < 114000:
                raise ValueError('Cannot verify allocated memory for supplemental step')
            env['SLURM_MEM_PER_NODE'] = memory[:-1]
            script = 'cd '+shlex.quote(str(root))+'\nsource scripts/scnet_environment.sh\nexec '+shlex.join([
                site['python'], 'scripts/goat_bulk_worker.py', '--scope', str(scope),
                '--scope-sha256', sha(scope), '--recover-startup'])
            command = ['srun', '--jobid='+job, '--overlap', '--kill-on-bad-exit=0',
                       '--nodes='+str(len(candidates)), '--ntasks='+str(len(candidates)),
                       '--ntasks-per-node=1', '--cpus-per-task=32', '--mem=114000M', '--nodelist='+','.join(sorted(candidates)),
                       'bash', '-lc', script]
            writej(journal, {**record, 'status': 'STARTUP_RECOVERY_INTENT'})
            # The original batch shell knows only its original srun. Protect it
            # before creating supplemental steps that could outlive that srun.
            originals = [readj(p) for p in (base/'nodes').glob(job+'_*.json')]
            originals = [n for n in originals if n.get('status') == 'RUNNING'
                         and n['hostname'] not in candidates and time.time()-n['updated_epoch'] < 60]
            if not originals:
                raise ValueError('No proven original supervisor available for lifetime guard')
            guard_host = originals[0]['hostname']
            guard_command = ['srun', '--jobid='+job, '--overlap', '--nodes=1', '--ntasks=1',
                             '--cpus-per-task=1', '--nodelist='+guard_host,
                             site['python'], str(root/'scripts/goat_bulk_parent_guard.py')]
            with (base/'startup_recovery'/(job+'_guard.log')).open('ab') as output:
                guard_process = subprocess.Popen(guard_command, env=env, stdin=subprocess.DEVNULL,
                                                 stdout=output, stderr=output, start_new_session=True)
            guard_path = base/'parent_guards'/(job+'.json')
            expires = time.monotonic()+90
            while time.monotonic() < expires:
                if guard_path.exists() and readj(guard_path).get('status') == 'PROTECTING_RECOVERY_STEPS':
                    break
                if guard_process.poll() is not None:
                    raise ValueError('Parent guard exited before protection')
                time.sleep(1)
            else:
                raise ValueError('Parent guard did not confirm protection; no new QM launched')
            with (base/'startup_recovery'/(job+'.log')).open('ab') as output:
                process = subprocess.Popen(command, env=env, stdin=subprocess.DEVNULL,
                                           stdout=output, stderr=output, start_new_session=True)
            record.update({'status': 'STARTUP_RECOVERY_SRUN_STARTED', 'pid': process.pid,
                           'guard_srun_pid': guard_process.pid, 'guard_hostname': guard_host,
                           'quantum_retries': 0, 'new_slurm_allocations': 0})
            writej(journal, record)
            reports.append(record)
    return reports


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--execute', action='store_true')
    print(json.dumps(recover(ROOT, execute=p.parse_args().execute), ensure_ascii=False, indent=2))
