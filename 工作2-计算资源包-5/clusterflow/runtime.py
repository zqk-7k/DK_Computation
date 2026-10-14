"""SCNET 6.1.0 adapter; no executable is probed implicitly.

Site evidence is generated explicitly on a compute node. Unknown liveness is
never equated with RUNNING or permission to delete a lock/checkpoint.
"""
from contextlib import contextmanager
from pathlib import Path
import json
import os
import re
import shutil
import signal
import socket
import subprocess
import tempfile
import threading
import time
import uuid
from .common import readj, sha, writej

SITE = Path(__file__).resolve().parents[1] / 'config/site_scnet.json'


def site_config():
    site = readj(Path(os.environ.get('WORK2_SITE_CONFIG', str(SITE))))
    profile = os.environ.get('WORK2_PACKED_PROFILE')
    if profile:
        if profile not in site.get('packed_profiles', {}):
            raise ValueError('Unknown packed resource profile: ' + profile)
        site['packed'] = {**site['packed'], **site['packed_profiles'][profile]}
    return site


def wrapper_preflight(executable, xtb=None):
    site = site_config()
    executable = Path(executable).resolve()
    if executable != Path(site['orca_wrapper']).resolve():
        return None  # Native executable path retains the old binary checks.
    reasons = []
    evidence_path = os.environ.get('WORK2_RUNTIME_EVIDENCE')
    evidence = readj(evidence_path) if evidence_path and Path(evidence_path).is_file() else {}
    if not executable.is_file() or not os.access(executable, os.X_OK):
        reasons.append('SCNET wrapper missing/not executable')
    if (not evidence or evidence.get('version') != '6.1.0'
            or evidence.get('wrapper_sha256') != (sha(executable) if executable.is_file() else None)):
        reasons.append('Run explicit site-check on a compute node; current wrapper/version evidence missing')
    for field in ('container', 'xtb'):
        path = Path(site[field])
        if not path.is_file() or evidence.get(field + '_sha256') != (sha(path) if path.is_file() else None):
            reasons.append('Runtime evidence mismatch: ' + field)
    if evidence.get('site_sha256') != sha(Path(os.environ.get('WORK2_SITE_CONFIG', str(SITE)))):
        reasons.append('Site profile changed since runtime verification')
    if xtb is not None and Path(xtb).resolve() != Path(site['xtb']).resolve():
        reasons.append('xTB differs from the explicitly verified SCNET route')
    return {'passed': not reasons, 'reasons': reasons, 'version': '6.1.0',
            'evidence': evidence_path, 'runtime_validation': 'VERSION_ONLY_REAL_QM_PILOT_STILL_REQUIRED',
            'xtb': site['xtb']}


def owner_record(token):
    return {'token': token, 'pid': os.getpid(), 'hostname': socket.gethostname(),
            'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
            'slurm_array_task_id': os.environ.get('SLURM_ARRAY_TASK_ID'),
            'process_start': process_start(os.getpid()), 'started_epoch': time.time()}


def process_start(pid):
    try:
        return Path(f'/proc/{int(pid)}/stat').read_text().rsplit(')', 1)[1].split()[19]
    except (OSError, ValueError, IndexError):
        return None


def lock_state(path):
    """Read-only, no shell, no guessing from mtime and no auto-unlink."""
    try:
        record = readj(path)
    except (OSError, ValueError):
        return 'UNKNOWN_REVIEW'
    if record.get('hostname') == socket.gethostname() and record.get('process_start'):
        return ('RUNNING' if process_start(record.get('pid')) == record['process_start']
                else 'INTERRUPTED_REVIEW')
    job = str(record.get('slurm_job_id') or '')
    array = record.get('slurm_array_task_id')
    if array is not None:
        job += '_' + str(array)
    if not re.fullmatch(r'\d+(?:_\d+)?', job):
        return 'UNKNOWN_REVIEW'
    try:
        result = subprocess.run(['sacct', '-n', '-P', '-j', job, '--format=JobIDRaw,State'],
                                capture_output=True, text=True, timeout=20, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return 'UNKNOWN_REVIEW'
    if result.returncode:
        return 'UNKNOWN_REVIEW'
    rows = [line.split('|') for line in result.stdout.splitlines()]
    states = [r[1].split()[0].split('+')[0] for r in rows if len(r)>1 and r[0] == job]
    if len(states) != 1:
        return 'UNKNOWN_REVIEW'
    if states[0] in {'COMPLETED','CANCELLED','FAILED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED'}:
        return 'INTERRUPTED_REVIEW'
    # A running parent does NOT prove a remote packed child is alive.
    return 'REMOTE_OWNER_CHECK_REQUIRED'


class StorageAdmissionError(ValueError):
    """A capacity wait is distinct from a broken node or corrupted result."""
    def __init__(self, category, path):
        self.category, self.path = category, str(path)
        super().__init__(category + ': ' + str(path))


def check_storage(path, minimum_gib, minimum_inodes=10000, write_probe=True):
    path = Path(path).resolve(strict=True)
    stat = os.statvfs(path)
    free = min(shutil.disk_usage(path).free, stat.f_bavail * stat.f_frsize)
    if free < minimum_gib * 1024**3:
        raise StorageAdmissionError('DISK_ADMISSION_BLOCKED', path)
    if stat.f_files and stat.f_favail < minimum_inodes:
        raise StorageAdmissionError('INODE_ADMISSION_BLOCKED', path)
    if stat.f_flag & getattr(os, 'ST_RDONLY', 1):
        raise StorageAdmissionError('READ_ONLY_NODE_OR_MOUNT', path)
    if write_probe:
        try:
            with tempfile.TemporaryFile(dir=path) as stream:
                stream.write(b'work2-write-probe'); stream.flush(); os.fsync(stream.fileno())
        except OSError as exc:
            import errno
            category = {errno.ENOSPC: 'DISK_ADMISSION_BLOCKED', errno.EDQUOT: 'QUOTA_ADMISSION_BLOCKED',
                        errno.EROFS: 'READ_ONLY_NODE_OR_MOUNT'}.get(exc.errno)
            if category:
                raise StorageAdmissionError(category, path) from exc
            raise
    return {'path': str(path), 'free_bytes': free, 'free_inodes': stat.f_favail}


@contextmanager
def execution_environment(executable, scratch):
    """Private wrapper sandbox per invocation prevents packed-job collisions.

    Keep the scientific scratch on /public; only sandbox/MPI runtime are local.
    This context is used by process-isolated workers, not concurrent threads.
    """
    site = site_config()
    if Path(executable).resolve() != Path(site['orca_wrapper']).resolve():
        yield
        return
    if not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('SCNET quantum execution requires a compute-node allocation')
    scratch = Path(scratch).resolve()
    if not scratch.is_relative_to(Path(site['shared_scratch']).resolve()):
        raise ValueError('SCNET scientific scratch must use the dedicated shared scratch root')
    check_storage(scratch, site['disk']['admit_free_gib'])
    local = Path(site['local_runtime_root'])
    check_storage(local, site['disk']['local_admit_gib'])
    check_storage('/dev/shm', site['disk']['shm_min_gib'], minimum_inodes=100)
    work = Path(tempfile.mkdtemp(prefix='work2-runtime-', dir=local))
    env = dict(site['environment'])
    env.update({'SLURM_TMPDIR': str(work), 'TMPDIR': str(work),
                'PROGLIST': site['xtb'], 'XTB_ROOT': str(Path(site['xtb']).parent.parent)})
    for name in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        env[name] = '1'
    old = {key: os.environ.get(key) for key in env}
    os.environ.update(env)
    stop = threading.Event()
    def heartbeat():
        while not stop.wait(60):
            try:
                writej(scratch / 'runtime_heartbeat.json',
                       {**owner_record('runtime'), 'updated_epoch': time.time(), 'private_runtime': str(work)})
                # Only owned directories, not a recursive touch of matrix files.
                for p in [work, *work.glob('orca610_sandbox_*')]:
                    if p.is_dir() and not p.is_symlink(): os.utime(p, None)
                free = shutil.disk_usage(scratch).free / 1024**3
                if free < site['disk']['stop_new_claims_gib']:
                    writej(scratch / 'LOW_SPACE_STOP_NEW_CLAIMS.json', {'free_gib': free, 'time': time.time()})
            except OSError:
                pass  # Main process/stage-out still enforce integrity; never fake completion.
    thread = threading.Thread(target=heartbeat, daemon=True); thread.start()
    try:
        yield
    finally:
        stop.set(); thread.join(timeout=2)
        for key, value in old.items():
            if value is None: os.environ.pop(key, None)
            else: os.environ[key] = value
        # No recursive cleanup here: wrapper owns sandbox cleanup, retained
        # runtime directories are audited alongside interrupted attempts.


def classify_failure(text, scheduler_state=''):
    text = text.lower()
    if 'read-only file system' in text or 'errno 30' in text: return 'NODE_READ_ONLY'
    if 'no space left' in text or 'disk quota exceeded' in text: return 'DISK_OR_QUOTA'
    if scheduler_state == 'OUT_OF_MEMORY' or any(v in text for v in ('out of memory','not enough memory','oom-kill')): return 'OOM'
    if scheduler_state in ('TIMEOUT','PREEMPTED','CANCELLED'): return 'INTERRUPTED'
    if 'squashfs' in text: return 'CONTAINER_IO'
    if any(v in text for v in ('mpi_win_allocate_shared','mpi_abort','segmentation fault')): return 'MPI_OR_SEGFAULT_REVIEW'
    if 'tmatrixcontainers' in text: return 'MATRIX_IO_REVIEW'
    if 'no valid' in text and ('structure' in text or 'conformer' in text): return 'GOAT_NO_VALID_ENSEMBLE'
    if 'startup' in text or 'propint' in text: return 'MODULE_FAILURE_REVIEW'
    return 'UNKNOWN_REVIEW'


@contextmanager
def termination_handoff():
    """SIGTERM/HUP become a receipt-producing interruption, not false success."""
    if threading.current_thread() is not threading.main_thread():
        yield
        return
    def interrupted(signum, frame):
        os.environ['WORK2_DRAIN_REQUESTED'] = '1'
        raise KeyboardInterrupt('Scheduler signal ' + str(signum))
    names = [signal.SIGTERM] + ([signal.SIGHUP] if hasattr(signal, 'SIGHUP') else [])
    old = {sig: signal.signal(sig, interrupted) for sig in names}
    try:
        yield
    finally:
        for sig, handler in old.items(): signal.signal(sig, handler)
