"""Admission budgets are estimates, never short per-calculation kill timers."""
from contextlib import contextmanager
import fcntl
import math
import os
from pathlib import Path
import time
import uuid

from .common import readj, writej


class DeferredWork(RuntimeError):
    def __init__(self, status, **details):
        self.report = {'status': status, 'quantum_jobs_launched': 0, **details}
        super().__init__(status)


def duration_budget(root, stage, species_id):
    from .cache import package_root
    from .runtime import site_config
    root = package_root(root)
    site = site_config()
    policy = site.get('walltime', {})
    rows = readj(root/'registry/batch_plan_v5.json').get('species', [])
    row = next((r for r in rows if r['species_id'] == species_id), {})
    complexity = row.get('complexity_class', 'X')
    estimates = policy.get('unmeasured_hours', {})
    hours = estimates.get(stage, {}).get(complexity, policy.get('fallback_hours', 240))
    seconds = float(hours) * 3600
    source = 'UNMEASURED_CONFIGURED_BUDGET_NOT_RUNTIME_PREDICTION'
    samples = []
    for path in (root/'operations/timings'/species_id/stage).glob('*.json'):
        sample = readj(path)
        value = float(sample['elapsed_seconds'])
        if math.isfinite(value) and value > 0:
            samples.append(value)
    if samples:
        # Do not extrapolate a fast small molecule to an unmeasured large one.
        seconds = max(seconds, max(samples) * float(policy.get('observed_safety_factor', 1.5)))
        source = 'MAX_CONFIGURED_AND_SAME_SPECIES_OBSERVED_WITH_MARGIN'
    if not math.isfinite(seconds) or seconds <= 0:
        raise ValueError('Invalid positive walltime admission budget')
    return seconds, source


def admit_quantum(root, stage, species_id):
    value = os.environ.get('WORK2_DEADLINE_EPOCH')
    if value is None:
        if os.environ.get('SLURM_JOB_ID'):
            from scripts.packed_worker import end_time
            deadline = end_time()
        else:
            return  # Portable non-SCNET runners retain their explicit timeout option.
    else:
        deadline = float(value)
    from .runtime import site_config
    budget, source = duration_budget(root, stage, species_id)
    margin = float(site_config().get('walltime', {}).get('stageout_margin_seconds', 3600))
    remaining = deadline - time.time()
    if not math.isfinite(remaining) or remaining < budget + margin:
        raise DeferredWork('WAITING_WALLTIME', stage=stage, species_id=species_id,
                           remaining_seconds=remaining, required_seconds=budget + margin,
                           estimate_source=source,
                           notice='No QM started; use a parent with sufficient time. No short QM timeout installed.')


def record_timing(root, stage, species_id, elapsed):
    from .cache import package_root
    writej(package_root(root)/'operations/timings'/species_id/stage/(uuid.uuid4().hex+'.json'),
           {'elapsed_seconds': elapsed, 'recorded_epoch': time.time()})


def _reconcile_cache_owner(root, digest, marker, directory):
    """An unlocked flock after a crash does not prove orphan QM has stopped."""
    from .operations import scheduler_terminal
    from .runtime import lock_state
    old=readj(marker)
    wait={'cache_key':digest,'previous_owner':old,
          'notice':'Reconcile the original execution and stage-out before another copy may run.'}
    if lock_state(marker)=='RUNNING' or not scheduler_terminal(old):
        raise DeferredWork('WAITING_SHARED_RESULT_RECOVERY',**wait)
    source=Path(old['source_directory']).resolve()
    if source!=directory:
        # A different workspace can only reuse a recovered, accepted execution;
        # it must never silently restart an interrupted original.
        if (not source.is_relative_to(root) or not (source/'receipt.json').is_file()
                or any((source/name).exists() for name in ('.orca-run.lock','.clusterflow-run.lock','.stageout-recovery.lock'))):
            raise DeferredWork('WAITING_SHARED_RESULT_RECOVERY',**wait)
        from .cache import key, publish
        from .common import config
        from .gates import inspect_job
        receipt=readj(source/'receipt.json')
        workspace=next(p for p in source.parents if (p/'config/protocol.json').is_file())
        if (receipt.get('status')!='completed' or receipt.get('reuse')
                or key(source,receipt.get('orca_executable_sha256'),receipt.get('runtime_evidence_sha256'))!=digest
                or not inspect_job(source,config(workspace))['accepted']):
            raise DeferredWork('WAITING_SHARED_RESULT_RECOVERY',**wait)
        publish(workspace,source)
        if not (root/'operations/cache'/(digest+'.json')).is_file():
            raise DeferredWork('WAITING_SHARED_RESULT_RECOVERY',**wait)
    # The original runner holds its job locks already. Reaching here for the
    # same directory requires prior owner reconciliation at that directory.
    archive=marker.parent/'owner_history';archive.mkdir(exist_ok=True)
    marker.rename(archive/(digest+'_'+uuid.uuid4().hex+'.json'))


@contextmanager
def single_execution(root, digest, directory=None):
    """Cross-process/cache-key ownership, released by the OS even on SIGKILL.

    The caller rechecks completed cache entries AFTER taking this lock and
    holds it through publishing. A busy dependency defers, never blocks a slot.
    DFT also records a durable owner so an orphan child cannot be duplicated
    from another workspace merely because its driver's flock was released.
    """
    from .cache import package_root
    root=package_root(root)
    path = root/'operations/cache_locks'/(digest+'.lock')
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise DeferredWork('WAITING_SHARED_RESULT', cache_key=digest)
        marker=path.with_suffix('.owner.json')
        token=None
        try:
            if directory is not None:
                from .runtime import owner_record
                directory=Path(directory).resolve()
                if not directory.is_relative_to(root):raise ValueError('Cache owner escapes package')
                if marker.exists():_reconcile_cache_owner(root,digest,marker,directory)
                token=uuid.uuid4().hex
                writej(marker,{**owner_record(token),'source_directory':str(directory)})
            yield
        finally:
            if token is not None and marker.is_file() and readj(marker).get('token')==token:
                marker.unlink()
            fcntl.flock(stream, fcntl.LOCK_UN)
