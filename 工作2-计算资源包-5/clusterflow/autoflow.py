"""Finite, idempotent GOAT -> exploratory NMR waves; no automatic QM by default.

An invocation releases/prepares work and executes at most ONE job only when
``execute=True``.  It never terminates a GOAT run on finding N conformers.
Numerical readiness is advice to a scientist, not a convergence certificate.
"""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import uuid
import time

from .common import config, hash_protocol, inside, readj, sha, writej
from . import pipeline

LOCKS = ('.goat-run.lock', '.orca-run.lock', '.clusterflow-run.lock', '.stageout-recovery.lock')
SCENARIO = 'DMSO_02_W01_trace_water_sensitivity'


def _safe(root, relative):
    """Confinement plus explicit rejection of links/junctions at every level."""
    root = Path(root).resolve()
    path = root / relative
    if Path(relative).is_absolute() or '..' in Path(relative).parts:
        raise ValueError('Unsafe automation-relative path')
    for entry in (path, *path.parents):
        if entry == root:
            break
        if entry.is_symlink() or (hasattr(entry, 'is_junction') and entry.is_junction()):
            raise ValueError('Linked automation/source paths are forbidden')
    return inside(root, relative)


@contextmanager
def _lock(directory):
    path = directory / '.autoflow.lock'
    token = uuid.uuid4().hex
    with path.open('x', encoding='utf-8') as stream:
        from .runtime import owner_record
        json.dump(owner_record(token), stream)
    try:
        yield
    finally:
        if path.is_file() and readj(path).get('token') == token:
            path.unlink()


def _plan(root, species_id):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,180}', species_id) or '..' in species_id:
        raise ValueError('Unsafe species identifier')
    plan = readj(_safe(root, 'registry/batch_plan_v5.json'))
    rows = plan['species']
    found = [r for r in rows if r['species_id'] == species_id]
    if len(found) != 1:
        raise ValueError('Exactly one registered batch plan required')
    row = found[0]
    for key in ('stable_conformer_target', 'initial_candidate_batch_size', 'initial_goat_starts', 'addition_goat_starts'):
        if type(row[key]) is not int or row[key] < 1:
            raise ValueError('Invalid positive batch plan count: ' + key)
    pool = row['goat_start_ids']
    if not pool or len(set(pool)) != len(pool):
        raise ValueError('Missing/duplicate GOAT pool entries')
    jobs = readj(_safe(root, 'registry/jobs.json'))
    indexed = {j['job_id']: j for j in jobs}
    if len(indexed) != len(jobs):
        raise ValueError('Duplicate GOAT registry identity')
    for job_id in pool:
        if (job_id not in indexed or indexed[job_id]['species_id'] != species_id
                or indexed[job_id]['scenario_id'] != SCENARIO):
            raise ValueError('Batch pool differs from fixed-composition GOAT registry')
    return row, indexed


def _job_status(root, row, *, goat):
    directory = _safe(root, row['job_directory'] if goat else row['job_dir'])
    locks = [directory / name for name in LOCKS if (directory / name).exists()]
    if locks:
        from .runtime import lock_state
        states = [lock_state(path) for path in locks]
        status = 'RUNNING' if 'RUNNING' in states else ('INTERRUPTED_REVIEW' if all(s == 'INTERRUPTED_REVIEW' for s in states) else 'OWNER_CHECK_REQUIRED')
        return {'status': status, 'terminal': status == 'INTERRUPTED_REVIEW',
                'directory': str(directory), 'lock_evidence': dict(zip(map(str, locks), states))}
    receipts = [directory / n for n in (('receipt.json', 'run_receipt.json') if goat else ('receipt.json',))
                if (directory / n).is_file()]
    outputs = list(directory.glob('goat.out' if goat else 'job.out'))
    if not receipts:
        return {'status': 'INCOMPLETE_REVIEW' if outputs else 'PENDING',
                'terminal': bool(outputs), 'directory': str(directory),
                'reason': 'Output without a terminal receipt is not a valid completion' if outputs else None}
    if len(receipts) != 1:
        return {'status': 'INVALID_REVIEW', 'terminal': True, 'reason': 'Multiple receipts'}
    try:
        receipt = readj(receipts[0])
        value = str(receipt.get('status', '')).upper()
        if value in ('RUNNING', 'STARTED'):
            return {'status': 'INCOMPLETE_REVIEW', 'terminal': True,
                    'reason': 'Nonterminal receipt without an ownership lock; reconcile scheduler/scratch instead of claiming RUNNING'}
        if value not in ('COMPLETED', 'COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW', 'FAILED', 'TIMEOUT', 'TIMED_OUT'):
            return {'status': 'INVALID_REVIEW', 'terminal': True, 'reason': 'Unknown terminal status'}
        if value in ('FAILED', 'TIMEOUT', 'TIMED_OUT'):
            return {'status': 'FAILED_REVIEW', 'terminal': True,
                    'receipt_sha256': sha(receipts[0]), 'reason': 'Inspect failure; never silently retry'}
        if goat:
            pipeline._goat_source(root, row['job_id'])
            return {'status': 'COMPLETE', 'terminal': True, 'receipt_sha256': sha(receipts[0])}
        from .gates import inspect_job
        checked = inspect_job(directory, config(root))
        if checked['accepted']:
            return {'status': 'COMPLETE', 'terminal': True, 'receipt_sha256': sha(receipts[0])}
        return {'status': 'SCIENTIFIC_REVIEW', 'terminal': True,
                'errors': checked['errors'], 'needs_review': checked.get('needs_review')}
    except (ValueError, OSError, KeyError, TypeError) as exc:
        return {'status': 'INVALID_REVIEW', 'terminal': True, 'reason': str(exc)}


def _copy_bound(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    expected = sha(source)
    if destination.exists():
        if sha(destination) != expected:
            raise ValueError('Refusing to replace a different frozen workspace artifact')
    else:
        with source.open('rb') as src, destination.open('xb') as dst:
            shutil.copyfileobj(src, dst)
    if sha(source) != expected or sha(destination) != expected:
        raise ValueError('Source changed during immutable workspace copy')
    return expected


def _workspace(root, wave, reports, indexed):
    workspace = _safe(root, wave['workspace'])
    workspace.mkdir(parents=True, exist_ok=True)
    bindings = {}
    for name in ('config/protocol.json', 'config/automation.json', 'registry/species.json', 'registry/batch_plan_v5.json'):
        bindings[name] = _copy_bound(_safe(root, name), _safe(workspace, name))
    copied_jobs = []
    for job_id in sorted({j for report in reports.values() for j in report['job_ids']}):
        row = indexed[job_id]
        source = pipeline._goat_source(root, job_id)
        directory = _safe(root, row['job_directory'])
        if any((directory / n).exists() for n in LOCKS):
            raise ValueError('GOAT source became active while snapshotting')
        receipt = source['receipt']
        files = set(receipt.get('input_sha256', {})) | set(receipt.get('output_sha256', {}))
        files |= {'goat.inp', 'start.xyz', 'cluster.mol', 'atom_map.json', 'provenance.json', source['receipt_path'].name}
        if (directory / 'isotope_map.json').is_file():
            files.add('isotope_map.json')
        for name in sorted(files):
            relative = (Path(row['job_directory']) / name).as_posix()
            bindings[relative] = _copy_bound(_safe(root, relative), _safe(workspace, relative))
        copied_jobs.append(row)
    target = workspace / 'registry/jobs.json'
    if target.exists() and readj(target) != copied_jobs:
        raise ValueError('Frozen wave GOAT registry changed')
    if not target.exists():
        writej(target, copied_jobs)
    snapshot = {'original_files_sha256': bindings, 'workspace_jobs_sha256': sha(target),
                'source_goat_ids': [j['job_id'] for j in copied_jobs],
                'notice': 'Cumulative independent workspace; no cached DFT reused across changed selections'}
    snapshot_path = workspace.parent / 'snapshot.json'
    if snapshot_path.exists() and readj(snapshot_path) != snapshot:
        raise ValueError('Frozen wave source snapshot changed')
    if not snapshot_path.exists():
        writej(snapshot_path, snapshot)
    return workspace


def _verify_workspace(root, wave):
    workspace = _safe(root, wave['workspace'])
    snapshot = readj(workspace.parent / 'snapshot.json')
    for relative, digest in snapshot['original_files_sha256'].items():
        source = _safe(root, relative)
        if sha(source) != digest or sha(_safe(workspace, relative)) != digest:
            raise ValueError('Frozen wave source binding changed: ' + relative)
    if sha(workspace / 'registry/jobs.json') != snapshot['workspace_jobs_sha256']:
        raise ValueError('Frozen wave registry changed')
    for row in readj(workspace / 'registry/jobs.json'):
        if any((_safe(root, row['job_directory']) / name).exists() for name in LOCKS):
            raise ValueError('Original frozen GOAT source is being modified/rerun')
    return workspace


def compare_predictions(before, after, *, h_ppm=0.05, c_ppm=0.5):
    """Advisory differences only: unchanged output cannot prove physical accuracy."""
    if not all(math.isfinite(v) and v > 0 for v in (h_ppm, c_ppm)):
        raise ValueError('Advisory comparison limits must be finite positive values')
    for key in ('species_id', 'scenario_id', 'temperature_K', 'j_method', 'protocol_hash'):
        if before['metadata'].get(key) != after['metadata'].get(key):
            raise ValueError('Cannot compare different NMR model identities: ' + key)
    def atoms(value):
        return {r['atom_index_zero_based']: r for r in value['atom_rows']}
    a, b = atoms(before), atoms(after)
    if set(a) != set(b):
        raise ValueError('NMR wave atom identities differ')
    deltas = {'1H': [], '13C': []}
    for idx in sorted(a):
        if a[idx]['nucleus'] != b[idx]['nucleus'] or a[idx].get('exchangeable_H') != b[idx].get('exchangeable_H'):
            raise ValueError('NMR nucleus/exchangeability identity differs')
        if a[idx]['nucleus'] == '1H' and a[idx].get('exchangeable_H'):
            continue
        delta = abs(float(a[idx]['shift_ppm']) - float(b[idx]['shift_ppm']))
        if not math.isfinite(delta):
            raise ValueError('Nonfinite chemical shift difference')
        deltas[a[idx]['nucleus']].append(delta)
    def pairs(value):
        return {(r['atom_i_zero_based'], r['atom_j_zero_based']): r for r in value['j_pairs']}
    aj, bj = pairs(before), pairs(after)
    if set(aj) != set(bj):
        raise ValueError('J-pair identities differ between waves')
    js = [abs(float(aj[k]['J_Hz']) - float(bj[k]['J_Hz'])) for k in aj
          if not aj[k].get('touches_exchangeable_H')]
    if not all(math.isfinite(v) for v in js):
        raise ValueError('Nonfinite J difference')
    maxima = {n: max(v) if v else None for n, v in deltas.items()}
    within = all(maxima[n] is not None and maxima[n] <= limit
                 for n, limit in (('1H', h_ppm), ('13C', c_ppm)))
    return {'status': 'ADVISORY_ONLY', 'max_nonexchangeable_1H_change_ppm': maxima['1H'],
            'max_13C_change_ppm': maxima['13C'], 'max_nonexchangeable_J_change_Hz': max(js) if js else None,
            'reference_limits_ppm': {'1H': h_ppm, '13C': c_ppm}, 'shifts_within_reference_limits': within,
            'J_pattern_requires_review': True, 'scientific_convergence_approved': False,
            'notice': 'Limits guide review; they are neither accuracy guarantees nor mandatory rejection cutoffs.'}


def _prediction(root, wave, species):
    workspace = _verify_workspace(root, wave)
    record = wave['prediction']
    path = _safe(workspace, record['path'])
    if sha(path) != record['sha256']:
        raise ValueError('Recorded wave prediction changed')
    saved = readj(path)
    from .collect import collect_species
    actual = collect_species(workspace, species, SCENARIO, exploratory=True, write=False)
    if any(saved[k] != actual[k] for k in ('atom_rows', 'j_pairs')):
        raise ValueError('Completed wave no longer recomputes from actual bound sources')
    return saved


def _release(state, plans, *, first=False):
    added = []
    for sid, plan in plans.items():
        current = state['released'][sid]
        held = [j for j in plan['goat_start_ids'] if j not in current]
        count = plan['initial_goat_starts'] if first else plan['addition_goat_starts']
        next_ids = held[:count]
        current.extend(next_ids)
        added.extend(next_ids)
    return added


def _execute_goat(root, row, orca, xtb, scratch_root):
    from scripts.goat_portable import run_one
    from .scheduling import single_execution, admit_quantum, record_timing, DeferredWork
    digest='goat_'+hashlib.sha256(row['job_id'].encode()).hexdigest()
    with single_execution(root,digest):
        status=_job_status(root,row,goat=True)
        if status['status']=='COMPLETE':return {'status':'REUSED_SHARED_GOAT','quantum_jobs_launched':0}
        if status['status']!='PENDING':raise DeferredWork('WAITING_SHARED_GOAT',job_status=status)
        admit_quantum(root,'goat',row['species_id'])
        started=time.monotonic()
        result=run_one(root, row['job_id'], Path(orca), Path(xtb), Path(scratch_root), acknowledge_models=True)
        record_timing(root,'goat',row['species_id'],time.monotonic()-started)
        return result


def _execute_dft(workspace, row, orca, scratch_root):
    from .runner import run_one
    from .scheduling import record_timing
    result=run_one(workspace, row['job_dir'], orca, acknowledge=True, scratch_root=scratch_root)
    if result.get('elapsed_seconds') and not result.get('reuse') and not result.get('skipped'):
        meta=readj(Path(workspace)/row['job_dir']/'job.json')
        record_timing(workspace,meta['stage'],meta['species_id'],result['elapsed_seconds'])
    return result


def _parallel_dft(workspace,row,orca,scratch_root,cpus,memory):
    if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,cpus)
    os.environ['SLURM_CPUS_PER_TASK']=str(len(cpus))
    os.environ['WORK2_SLOT_MEMORY_MB']=str(memory)
    from .scheduling import DeferredWork
    try:return _execute_dft(workspace,row,orca,scratch_root)
    except DeferredWork as exc:return exc.report


def _run_pending(workspace,rows,orca,scratch_root):
    """Bounded independent conformer parallelism in separate processes."""
    maximum=int(os.environ.get('WORK2_BRANCH_PARALLEL_JOBS','1'))
    if maximum<1:raise ValueError('Positive branch parallelism required')
    chosen=rows[:maximum]
    if len(chosen)==1:
        executed=_execute_dft(workspace,chosen[0],orca,scratch_root)
        return {'status':'REUSED_ONE_DFT_RESULT' if executed.get('reuse') else 'EXECUTED_ONE_DFT_JOB',
                'execution':executed,'quantum_jobs_launched':executed.get('quantum_jobs_launched',1)}
    from concurrent.futures import ProcessPoolExecutor
    available=sorted(os.sched_getaffinity(0))
    memory=int(os.environ['WORK2_SLOT_MEMORY_MB'])//len(chosen)
    width=len(available)//len(chosen)
    for row in chosen:
        meta=readj(Path(workspace)/row['job_dir']/'job.json')
        if meta['nprocs']>width or meta['nprocs']*meta['maxcore_mb']>.75*memory:
            raise ValueError('Branch parallelism exceeds its isolated CPU/memory budget')
    with ProcessPoolExecutor(max_workers=len(chosen)) as pool:
        futures=[pool.submit(_parallel_dft,workspace,row,orca,scratch_root,
                            available[i*width:(i+1)*width],memory) for i,row in enumerate(chosen)]
        results=[future.result() for future in futures]
    launched=sum(r.get('quantum_jobs_launched',1) for r in results)
    if launched:status='EXECUTED_DFT_BATCH'
    elif any(r.get('status')=='WAITING_WALLTIME' for r in results):status='WAITING_WALLTIME'
    elif any(r.get('status','').startswith('WAITING_SHARED_RESULT') for r in results):
        status=next(r['status'] for r in results if r.get('status','').startswith('WAITING_SHARED_RESULT'))
    else:status='REUSED_DFT_BATCH'
    return {'status':status,'executions':results,'quantum_jobs_launched':launched}


def _stage_step(workspace, plans, *, execute, orca, scratch_root):
    """Advance all accepted members, never silently discard a failed conformer."""
    from .inputs import STAGES
    ph = hash_protocol(config(workspace))
    settings=readj(Path(workspace)/'config/automation.json')
    if settings.get('incremental_downstream') and execute:
        pending=[]
        for sid in plans:
            pipeline.prepare_prescreen(workspace,sid,SCENARIO)
            pipeline.prepare_healthy_downstream(workspace,sid,SCENARIO)
            for stage in ('optfreq','shielding','j_PBE_pcJ3','j_PBE0_pcJ2','prescreen'):
                for row in pipeline.manifest(workspace,stage,sid,SCENARIO):
                    if _job_status(workspace,row,goat=False)['status']=='PENDING':pending.append(row)
        if pending and execute:return _run_pending(workspace,pending,orca,scratch_root)
    for sid in plans:
        pipeline._validated_imports(workspace, sid, SCENARIO)
        pipeline.prepare_prescreen(workspace, sid, SCENARIO)
        for stage in ('prescreen', 'optfreq'):
            rows = pipeline.manifest(workspace, stage, sid, SCENARIO)
            statuses = [_job_status(workspace, row, goat=False) for row in rows]
            pending = [r for r, s in zip(rows, statuses) if s['status'] == 'PENDING']
            if pending and execute:
                return _run_pending(workspace,pending,orca,scratch_root)
            if any(not s['terminal'] for s in statuses):
                return {'status': 'WAITING_DFT_BATCH', 'stage': stage, 'species_id': sid,
                        'jobs': rows, 'job_statuses': statuses}
            if any(s['status'] != 'COMPLETE' for s in statuses):
                return {'status': 'AWAITING_SCIENTIFIC_REVIEW', 'stage': stage, 'species_id': sid,
                        'job_statuses': statuses, 'notice': 'Diagnose and inspect; no automatic retry or scientific approval.'}
            if stage == 'prescreen':
                result = pipeline.advance(workspace, sid, SCENARIO, stage)
                if result['status'] != 'advanced':
                    return {'status': 'AWAITING_SCIENTIFIC_REVIEW', 'stage_result': result}
    # Do NOT call advance(optfreq) before this gate: it creates expensive NMR jobs.
    stable_reports = {sid: pipeline.final_stable_candidates(workspace, sid, SCENARIO) for sid in plans}
    if any(not report['all_prepared_optfreq_pass'] or report['issues'] for report in stable_reports.values()):
        return {'status': 'AWAITING_SCIENTIFIC_REVIEW', 'stable_reports': stable_reports,
                'notice': 'Unverified or human-reviewed approximate minima cannot enter this automatic stable-minimum workflow.'}
    from .review_control import target_for
    effective_targets={sid:target_for(workspace,sid,SCENARIO,plans[sid]['stable_conformer_target'])[0] for sid in plans}
    short = [sid for sid in plans if stable_reports[sid]['eligible_count'] < effective_targets[sid]]
    if short:
        return {'status': 'DFT_STABLE_MINIMA_SHORTAGE', 'short_species': short,
                'stable_reports': stable_reports,
                'minimum_for_first_NMR': {sid: p['stable_conformer_target'] for sid, p in plans.items()},
                'notice': 'Only distinct checked low-energy final DFT minima count; independent properties may exist, but final ensemble is not approved.'}
    for sid in plans:
        result = pipeline.advance(workspace, sid, SCENARIO, 'optfreq')
        if result['status'] != 'advanced':
            return {'status': 'AWAITING_SCIENTIFIC_REVIEW', 'stage_result': result}
        rows = [r for stage in STAGES if stage not in ('prescreen', 'optfreq')
                for r in pipeline.manifest(workspace, stage, sid, SCENARIO)]
        statuses = [_job_status(workspace, row, goat=False) for row in rows]
        pending = [r for r, s in zip(rows, statuses) if s['status'] == 'PENDING']
        if pending and execute:
            return _run_pending(workspace,pending,orca,scratch_root)
        if any(not s['terminal'] for s in statuses):
            return {'status': 'WAITING_PROPERTY_BATCH', 'species_id': sid, 'jobs': rows, 'job_statuses': statuses}
        if any(s['status'] != 'COMPLETE' for s in statuses):
            return {'status': 'AWAITING_SCIENTIFIC_REVIEW', 'stage': 'properties', 'job_statuses': statuses}
        pipeline.validate_ensemble_sources(workspace, sid, SCENARIO)
    return {'status': 'WAVE_DFT_COMPLETE', 'stable_reports': stable_reports, 'quantum_jobs_launched': 0}


def _expand_candidates(state, plans, species_ids):
    """Consume existing reserves before releasing more GOAT sources."""
    wave = state['waves'][-1]
    changes = {}
    for sid in species_ids:
        plan = plans[sid]
        available = wave['selections'][sid].get('energy_eligible_count', wave['selections'][sid]['available_unique_count'])
        current = len(wave['selections'][sid]['selected'])
        if available > current:
            state['targets'][sid] = min(available, current + plan['initial_candidate_batch_size'])
            changes[sid] = {'kind': 'EXISTING_GOAT_RESERVE', 'candidate_target': state['targets'][sid]}
        else:
            held = [j for j in plan['goat_start_ids'] if j not in state['released'][sid]]
            added = held[:plan['addition_goat_starts']]
            if added:
                state['released'][sid].extend(added)
                state['targets'][sid] = current + plan['initial_candidate_batch_size']
                changes[sid] = {'kind': 'ADDITIONAL_WHOLE_GOAT_BATCH', 'newly_released_ids': added,
                                'candidate_target': state['targets'][sid]}
    if changes:
        state['phase'] = 'SAMPLING'
    return changes


def step(root, species_id, *, execute=False, orca=None, xtb=None, scratch_root=None):
    """Prepare/update one species' finite sampling controller, optionally run one job.

    Repeated calls (CLI watch) implement automatic handoff. Ordinary calls never
    execute ORCA. Existing failed jobs are never automatically overwritten.
    """
    root = Path(root).resolve()
    target, indexed = _plan(root, species_id)
    settings = readj(_safe(root, 'config/automation.json'))
    if settings.get('automatic_handoff', {}).get('enabled') is not True:
        raise ValueError('Automatic handoff is not explicitly enabled')
    if execute and not all((orca, xtb, scratch_root)):
        raise ValueError('Explicit execution requires ORCA, xTB and dedicated ASCII scratch paths')
    plans = {species_id: target}
    if species_id != 'TMS_DMSO':
        plans['TMS_DMSO'] = _plan(root, 'TMS_DMSO')[0]
    base = _safe(root, 'automation/' + species_id)
    base.mkdir(parents=True, exist_ok=True)
    state_path = base / 'state.json'
    identities = {name: sha(_safe(root, name)) for name in
                  ('registry/batch_plan_v5.json', 'registry/jobs.json', 'registry/species.json', 'config/protocol.json', 'config/automation.json')}
    with _lock(base):
        if state_path.exists():
            state = readj(state_path)
            if state.get('identity_sha256') != identities or state.get('species_id') != species_id:
                raise ValueError('Controller plan/protocol changed; preserve this run and start a new reviewed run')
        else:
            state = {'schema_version': '5.0', 'species_id': species_id, 'identity_sha256': identities,
                     'released': {sid: [] for sid in plans}, 'waves': [], 'phase': 'SAMPLING',
                     'targets': {sid: p['initial_candidate_batch_size'] for sid, p in plans.items()},
                     'scientific_convergence_approved': False}
            _release(state, plans, first=True)
            writej(state_path, state)

        def finish(status, **details):
            state['last_status'] = status
            state['updated_utc'] = datetime.now(timezone.utc).isoformat()
            writej(state_path, state)
            report = {'status': status, 'species_id': species_id, 'phase': state['phase'],
                      'released_goat_ids': state['released'], 'current_dft_targets': state['targets'],
                      'completed_waves': sum('prediction' in w for w in state['waves']),
                      'scientific_convergence_approved': False, 'quantum_jobs_launched': 0, **details}
            writej(base / 'latest_status.json', report)
            return report

        for previous in state['waves']:
            if (root / previous['workspace']).parent.joinpath('snapshot.json').is_file():
                previous_workspace = _verify_workspace(root, previous)
                if previous.get('stable_reports'):
                    fresh_stable = {sid: pipeline.final_stable_candidates(previous_workspace, sid, SCENARIO)
                                    for sid in previous['stable_reports']}
                    if fresh_stable != previous['stable_reports']:
                        raise ValueError('Frozen completed stable-minimum audit changed')
            if 'prediction' in previous:
                _prediction(root, previous, species_id)
        if state['phase'] == 'REVIEW':
            return finish('AWAITING_SCIENTIFIC_REVIEW', reason=state.get('review_reason'), advice=state.get('advice'),
                          stable_reports=state['waves'][-1].get('stable_reports') if state['waves'] else None)

        if state['phase'] == 'SAMPLING':
            checks = {sid: {j: _job_status(root, indexed[j], goat=True) for j in ids}
                      for sid, ids in state['released'].items()}
            writej(base / 'released_goat_manifest.json', [indexed[j] for ids in state['released'].values() for j in ids])
            if any(not s['terminal'] for by_id in checks.values() for s in by_id.values()):
                pending = [indexed[j] for by_id in checks.values() for j, s in by_id.items() if s['status'] == 'PENDING']
                if execute and pending:
                    result = _execute_goat(root, pending[0], orca, xtb, scratch_root)
                    return finish('EXECUTED_ONE_GOAT_JOB', execution=result, quantum_jobs_launched=result.get('quantum_jobs_launched',1))
                return finish('WAITING_WHOLE_GOAT_BATCH', job_statuses=checks)
            from .autoselect import screen_batch, import_selection
            reports = {}
            for sid in plans:
                good = [j for j, s in checks[sid].items() if s['status'] == 'COMPLETE']
                previous = state['waves'][-1]['selections'][sid]['selected'] if state['waves'] else []
                report = screen_batch(root, good, state['targets'][sid],
                                      previous_selected_ids=[r['conformer_id'] for r in previous], write=False) if good else {
                                          'target_reached': False, 'available_unique_count': 0, 'selected': []}
                reports[sid] = {**report, 'job_ids': good}
            # Candidate batch size is an upper work allocation, not the requested
            # stable-minimum quota. If a complete pool yields fewer candidates,
            # actual DFT may still meet the final quota; do not invent duplicates.
            for sid, report in list(reports.items()):
                available = report.get('energy_eligible_count', report['available_unique_count'])
                no_held = len(state['released'][sid]) == len(plans[sid]['goat_start_ids'])
                if (not report['target_reached'] and available
                        and (no_held or available >= plans[sid]['stable_conformer_target'])):
                    previous = state['waves'][-1]['selections'][sid]['selected'] if state['waves'] else []
                    reports[sid] = screen_batch(root, report['job_ids'], available,
                                               previous_selected_ids=[r['conformer_id'] for r in previous], write=False)
                    state['targets'][sid] = available
            if any(not r['target_reached'] for r in reports.values()):
                added = _release(state, plans)
                if added:
                    return finish('RELEASED_ADDITIONAL_GOAT_BATCH', newly_released_ids=added, selections=reports, job_statuses=checks)
                state['phase'] = 'REVIEW'
                state['review_reason'] = 'Finite GOAT pool exhausted below requested independent-candidate target; never pad duplicates.'
                return finish('AWAITING_SCIENTIFIC_REVIEW', selections=reports, job_statuses=checks, reason=state['review_reason'])
            if state['waves']:
                expanded = False
                for sid, report in reports.items():
                    previous_ids = {r['conformer_id'] for r in state['waves'][-1]['selections'][sid]['selected']}
                    mapped = {report.get('previous_to_current_representatives', {}).get(cid, cid) for cid in previous_ids}
                    if {r['conformer_id'] for r in report['selected']} - mapped:
                        expanded = True
                if not expanded:
                    added = _release(state, plans)
                    if added:
                        return finish('NO_NEW_DISTINCT_CANDIDATE_RELEASED_NEXT_GOAT_BATCH', newly_released_ids=added, selections=reports)
                    state['phase'] = 'REVIEW'
                    state['review_reason'] = 'No additional distinct eligible candidates remain; do not present a repeated identical DFT ensemble as an independent expansion test.'
                    return finish('AWAITING_SCIENTIFIC_REVIEW', reason=state['review_reason'], selections=reports)
            wave_number = len(state['waves']) + 1
            wave = {'number': wave_number,
                    'workspace': f'automation/{species_id}/wave_{wave_number:03d}/workspace',
                    'selections': reports, 'targets': dict(state['targets']), 'goat_job_statuses': checks}
            workspace = _workspace(root, wave, reports, indexed)
            # The report is recomputed inside the snapshot, not edited into approval.
            for sid in plans:
                clean = screen_batch(workspace, reports[sid]['job_ids'], state['targets'][sid],
                                     previous_selected_ids=[r['conformer_id'] for r in (state['waves'][-1]['selections'][sid]['selected'] if state['waves'] else [])], write=False)
                import_selection(workspace, clean)
            state['waves'].append(wave)
            state['phase'] = 'DFT'
            writej(state_path, state)

        wave = state['waves'][-1]
        workspace = _verify_workspace(root, wave)
        stage = _stage_step(workspace, plans, execute=execute, orca=orca, scratch_root=scratch_root)
        if stage['status'] == 'DFT_STABLE_MINIMA_SHORTAGE':
            wave['stable_reports'] = stage['stable_reports']
            changes = _expand_candidates(state, plans, stage['short_species'])
            if changes:
                return finish('DFT_STABLE_MINIMA_SHORTAGE_AUTO_REPLENISH', changes=changes, stable_reports=stage['stable_reports'])
            state['phase'] = 'REVIEW'
            state['review_reason'] = 'Final stable-minimum target unmet after available reserve/GOAT pool exhaustion; independent property work may exist, but final ensemble remains unapproved.'
            return finish('AWAITING_SCIENTIFIC_REVIEW', reason=state['review_reason'], stable_reports=stage['stable_reports'])
        if stage['status'] != 'WAVE_DFT_COMPLETE':
            return finish(stage.pop('status'), workspace=str(workspace), **stage)
        wave['stable_reports'] = stage['stable_reports']
        if species_id == 'TMS_DMSO':
            state['phase'] = 'REVIEW'
            state['review_reason'] = 'TMS reference DFT wave complete; product controllers bind matching reference ensembles independently.'
            return finish('TMS_REFERENCE_WAVE_COMPLETE_REVIEW_REQUIRED', workspace=str(workspace))
        from .collect import collect_species
        prediction = collect_species(workspace, species_id, SCENARIO, exploratory=True, write=True)
        wave['prediction'] = {'path': prediction['result_path'], 'sha256': sha(_safe(workspace, prediction['result_path']))}
        from .spectra import make_spectrum, export_spectrum
        directory = _safe(workspace, prediction['output_directory']) / 'spectra'
        wave['spectra'] = {}
        for nucleus in ('1H', '13C'):
            spectrum = make_spectrum(prediction.get('spectrum_atom_rows', prediction['atom_rows']), prediction['j_pairs'], nucleus, config(workspace)['spectrum'])
            wave['spectra'][nucleus] = export_spectrum(spectrum, directory / nucleus, species_id + ' | exploratory wave ' + str(wave['number']))
        completed_before = [w for w in state['waves'][:-1] if 'prediction' in w]
        if completed_before:
            previous_nmr = completed_before[-1]
            before = _prediction(root, previous_nmr, species_id)
            limits = settings.get('stability_advice', {})
            advice = compare_predictions(before, prediction,
                                         h_ppm=float(limits.get('nonexchangeable_H_ppm', .05)),
                                         c_ppm=float(limits.get('C_ppm', .5)))
            advice['selection_reports'] = {sid: {'before': previous_nmr['selections'][sid], 'after': wave['selections'][sid]} for sid in plans}
            state['advice'] = advice
            writej(base / 'nmr_stability_advice.json', advice)
            # Stop for inspection, not PASS: motifs, J splitting, solvent adequacy
            # and unseen conformational basins remain scientific judgments.
            if advice['shifts_within_reference_limits']:
                state['phase'] = 'REVIEW'
                state['review_reason'] = 'Two real cumulative NMR waves have small shift changes; inspect J patterns, motifs and rediscovery before any scientific approval.'
                return finish('AWAITING_SCIENTIFIC_REVIEW', reason=state['review_reason'], advice=advice)
        changes = _expand_candidates(state, plans, plans)
        if not changes:
            state['phase'] = 'REVIEW'
            state['review_reason'] = 'Finite GOAT pool exhausted; available waves and coverage require scientific review, not automatic convergence.'
            return finish('AWAITING_SCIENTIFIC_REVIEW', reason=state['review_reason'], advice=state.get('advice'))
        return finish('FIRST_OR_EXPANDED_NMR_WAVE_COMPLETE_NEXT_BATCH_RELEASED', changes=changes,
                      prediction=wave['prediction'], notice='Next wave recomputes a cumulative larger selected ensemble; no final conformer cap or accuracy assertion.')
