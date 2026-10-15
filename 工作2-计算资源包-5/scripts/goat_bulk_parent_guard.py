"""Keep an original batch step alive until its supplemental node steps finish.

Only the Python supervisor (heartbeat writer) is paused; its independent GOAT
worker children and all ORCA processes continue normally. This prevents an
original batch shell from exiting and killing a still-running recovery step.
Unexpected guard failures keep the supervisor held for explicit recovery.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import signal
import socket
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from clusterflow.common import readj, writej
from clusterflow.runtime import process_start, owner_record


def supervisor(node, proc=Path('/proc')):
    parents = set()
    for worker in node['workers']:
        if not worker['alive']:
            continue
        try:
            fields = (proc/str(worker['pid'])/'stat').read_text().rsplit(')', 1)[1].split()
        except FileNotFoundError:
            continue
        if fields[0] != 'Z':
            parents.add(int(fields[1]))
    if len(parents) != 1:
        raise ValueError('Cannot prove one live supervisor for the original workers')
    pid = parents.pop()
    command = (proc/str(pid)/'cmdline').read_bytes().split(b'\0')
    if not any(arg.endswith(b'/goat_bulk_worker.py') for arg in command) or b'--recover-startup' in command:
        raise ValueError('Refusing to pause anything except the original pool supervisor')
    return pid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replace-guard', action='store_true')
    args = parser.parse_args()
    job = os.environ['SLURM_JOB_ID']; base = ROOT/'operations/goat_bulk_20261015'
    recovery = readj(base/'startup_recovery'/(job+'.json'))
    if socket.gethostname() in recovery['nodes']:
        raise ValueError('Guard must run on an original worker node')
    path = base/'parent_guards'/(job+'.json'); path.parent.mkdir(parents=True, exist_ok=True)
    previous = readj(path) if path.exists() else {}
    if args.replace_guard and previous:
        old_pid = previous.get('pid')
        if previous.get('hostname') != socket.gethostname() or previous.get('slurm_job_id') != job:
            raise ValueError('Guard replacement must stay on the same node and allocation')
        if process_start(old_pid) == previous.get('process_start'):
            command = Path('/proc')/str(old_pid)/'cmdline'
            if not any(arg.endswith(b'/goat_bulk_parent_guard.py') for arg in command.read_bytes().split(b'\0')):
                raise ValueError('Previous process is not the expected guard')
            # Replace only this bookkeeping guard, preserving the held supervisor
            # and all quantum children. SIGTERM would resume the old supervisor.
            os.kill(old_pid, signal.SIGKILL)
            time.sleep(1)
    lock = path.with_suffix('.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    node = readj(base/'nodes'/(job+'_'+socket.gethostname()+'.json'))
    pid = previous.get('supervisor_pid') if args.replace_guard else None
    adopt = pid and process_start(pid) == previous.get('supervisor_start')
    if adopt:
        command = (Path('/proc')/str(pid)/'cmdline').read_bytes().split(b'\0')
        if not any(arg.endswith(b'/goat_bulk_worker.py') for arg in command):
            raise ValueError('Held supervisor identity changed')
        started = previous['supervisor_start']
    else:
        if time.time()-node['updated_epoch'] > 90:
            raise ValueError('Original supervisor evidence is stale; do not signal')
        pid = supervisor(node); started = process_start(pid)
    if not started or pid == os.getpid():
        raise ValueError('Unknown supervisor identity')
    record = {**owner_record('parent-guard'), 'supervisor_pid': pid, 'supervisor_start': started,
              'recovery_nodes': recovery['nodes'], 'status': 'GUARD_INTENT',
              'notice': 'Supervisor heartbeat writer only; quantum worker children are not stopped.'}
    writej(path, record)
    resume = False
    def interrupted(signum, frame):
        nonlocal resume
        resume = True
        raise SystemExit('Guard signal '+str(signum))
    for sig in (signal.SIGTERM, signal.SIGHUP, signal.SIGINT):
        signal.signal(sig, interrupted)
    paused = False
    try:
        os.kill(pid, signal.SIGSTOP); paused = True
        while True:
            if process_start(pid) != started:
                raise RuntimeError('Original supervisor disappeared')
            try:
                pending = []
                for hostname in recovery['nodes']:
                    p = base/'nodes'/(job+'_'+hostname+'.json')
                    if not p.exists() or readj(p).get('status') != 'EXITED':
                        pending.append(hostname)
                # Shared-file replacement can briefly return ENOENT on this
                # filesystem. Never interpret a transient read as completion.
                active = ['PENDING_NODE_WORKERS']
                if not pending:
                    from scripts.packed_worker import queue_lock
                    with queue_lock(base):
                        queue = readj(base/'queue.json')
                    active = [v['job_id'] for v in queue['active'].values() if v['hostname'] in recovery['nodes']]
            except (OSError, json.JSONDecodeError):
                time.sleep(2)
                continue
            record.update({'status': 'PROTECTING_RECOVERY_STEPS', 'updated_epoch': time.time(),
                           'pending_nodes': pending, 'recovery_active_claims': len(active)})
            writej(path, record)
            if not pending and not active:
                resume = True
                break
            time.sleep(30)
    finally:
        if resume and paused and process_start(pid) == started:
            os.kill(pid, signal.SIGCONT)
        record.update({'status': 'SUPERVISOR_RESUMED' if resume else 'GUARD_FAILED_SUPERVISOR_HELD',
                       'updated_epoch': time.time()})
        writej(path, record)


if __name__ == '__main__':
    main()
