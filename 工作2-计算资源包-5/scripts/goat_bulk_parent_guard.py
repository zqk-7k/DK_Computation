"""Keep an original batch step alive until its supplemental node steps finish.

Only the Python supervisor (heartbeat writer) is paused; its independent GOAT
worker children and all ORCA processes continue normally. This prevents an
original batch shell from exiting and killing a still-running recovery step.
The supervisor is always resumed on normal guard exit or SIGTERM/SIGHUP.
"""
import fcntl
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
    job = os.environ['SLURM_JOB_ID']; base = ROOT/'operations/goat_bulk_20261015'
    recovery = readj(base/'startup_recovery'/(job+'.json'))
    if socket.gethostname() in recovery['nodes']:
        raise ValueError('Guard must run on an original worker node')
    path = base/'parent_guards'/(job+'.json'); path.parent.mkdir(parents=True, exist_ok=True)
    lock = path.with_suffix('.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    node = readj(base/'nodes'/(job+'_'+socket.gethostname()+'.json'))
    if time.time()-node['updated_epoch'] > 90:
        raise ValueError('Original supervisor evidence is stale; do not signal')
    pid = supervisor(node); started = process_start(pid)
    if not started or pid == os.getpid():
        raise ValueError('Unknown supervisor identity')
    record = {**owner_record('parent-guard'), 'supervisor_pid': pid, 'supervisor_start': started,
              'recovery_nodes': recovery['nodes'], 'status': 'GUARD_INTENT',
              'notice': 'Supervisor heartbeat writer only; quantum worker children are not stopped.'}
    writej(path, record)
    def interrupted(signum, frame):
        raise SystemExit('Guard signal '+str(signum))
    for sig in (signal.SIGTERM, signal.SIGHUP, signal.SIGINT):
        signal.signal(sig, interrupted)
    paused = False
    try:
        os.kill(pid, signal.SIGSTOP); paused = True
        while True:
            if process_start(pid) != started:
                raise RuntimeError('Original supervisor disappeared')
            pending = []
            for hostname in recovery['nodes']:
                p = base/'nodes'/(job+'_'+hostname+'.json')
                if not p.exists() or readj(p).get('status') != 'EXITED':
                    pending.append(hostname)
            queue = readj(base/'queue.json')
            active = [v['job_id'] for v in queue['active'].values() if v['hostname'] in recovery['nodes']]
            record.update({'status': 'PROTECTING_RECOVERY_STEPS', 'updated_epoch': time.time(),
                           'pending_nodes': pending, 'recovery_active_claims': len(active)})
            writej(path, record)
            if not pending and not active:
                break
            time.sleep(30)
    finally:
        if paused and process_start(pid) == started:
            os.kill(pid, signal.SIGCONT)
        record.update({'status': 'SUPERVISOR_RESUMED', 'updated_epoch': time.time()})
        writej(path, record)


if __name__ == '__main__':
    main()
