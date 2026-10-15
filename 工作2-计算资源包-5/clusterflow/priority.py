"""Scheduling estimates only; never alter chemical targets or acceptance gates."""
from datetime import datetime
from pathlib import Path
import statistics
import time

from .common import readj, writej

SUCCESS = {'COMPLETED', 'COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW'}
BENCHMARKS = ['TMS_DMSO', 'F3_N01_source__DMSO', 'F3_N15_source__DMSO',
              'F2_E12_rep_R__DMSO', 'F2_E01_rep_R__DMSO', 'F3_B01_source__DMSO']


def successful_hours(receipt):
    if str(receipt.get('status', '')).upper() not in SUCCESS or receipt.get('reuse'):
        return None
    try:
        value = receipt.get('elapsed_seconds')
        if value is None:
            start = receipt.get('started_utc') or receipt.get('started_at')
            end = receipt.get('finished_utc') or receipt.get('finished_at')
            value = (datetime.fromisoformat(end) - datetime.fromisoformat(start)).total_seconds()
        value = float(value) / 3600
        return value if 0 < value < 1e6 else None
    except (ValueError, TypeError):
        return None


def estimate(root, row, jobs):
    """Remaining first-wave work, revised using successful same-species receipts.

    Descriptor scaling is a provisional ranking heuristic, not a timing claim.
    Shared TMS work is counted once, under its own branch.
    """
    root = Path(root); sid = row['species_id']
    state_path = root/'automation'/sid/'state.json'
    state = readj(state_path) if state_path.is_file() else {}
    descriptors = row.get('descriptors', {})
    scale = max(.65, min(3., (descriptors.get('heavy_atoms', 12) + 11) / 23))
    flex = 1 + .06 * descriptors.get('rotatable_bonds_strict', 0)
    defaults = {'goat': 39.65 * scale * flex, 'prescreen': 5.68 * scale,
                'optfreq': 72.34 * scale, 'shielding': 2.19 * scale,
                'j_PBE_pcJ3': 4.5 * scale, 'j_PBE0_pcJ2': 4.5 * scale}
    counts = {}; samples = {stage: [] for stage in defaults}
    released = state.get('released', {}).get(sid, row.get('goat_start_ids', [])[:row['initial_goat_starts']])
    done = 0
    for job_id in released:
        job = jobs.get(job_id)
        if not job: continue
        receipt = root/job['job_directory']/'receipt.json'
        if receipt.is_file():
            data = readj(receipt)
            if str(data.get('status', '')).upper() in SUCCESS:
                done += 1
                hours = successful_hours(data)
                if hours is not None: samples['goat'].append(hours)
    counts['goat'] = max(0, len(released) - done) if state.get('phase', 'SAMPLING') == 'SAMPLING' else 0
    target = state.get('targets', {}).get(sid, row['initial_candidate_batch_size'])
    counts.update({stage: target for stage in defaults if stage != 'goat'})
    if sid == 'TMS_DMSO': counts['j_PBE_pcJ3'] = counts['j_PBE0_pcJ2'] = 0
    elif row.get('J_method') == 'PBE_pcJ3': counts['j_PBE0_pcJ2'] = 0
    waves = state.get('waves', [])
    if waves:
        workspace = root/waves[-1]['workspace']
        completed = {}
        for path in (workspace/'jobs').glob('*/*/*/*/*/job.json'):
            meta = readj(path)
            if meta.get('species_id') != sid: continue
            stage = meta.get('stage'); receipt = path.with_name('receipt.json')
            if stage not in samples or not receipt.is_file(): continue
            data = readj(receipt)
            if str(data.get('status', '')).upper() in SUCCESS:
                completed[stage] = completed.get(stage, 0) + 1
                hours = successful_hours(data)
                if hours is not None: samples[stage].append(hours)
        for stage, number in completed.items(): counts[stage] = max(0, counts[stage] - number)
    per_job = {stage: statistics.median(samples[stage]) if samples[stage] else value
               for stage, value in defaults.items()}
    hours = sum(counts[stage] * per_job[stage] for stage in counts)
    return {'species_id': sid, 'complexity_class': row['complexity_class'],
            'remaining_serial_hours_estimate': round(hours, 2), 'pending_counts_estimate': counts,
            'per_job_hours_estimate': per_job, 'completed_waves': len(waves),
            'phase': state.get('phase', 'NOT_STARTED'),
            'notice': 'Scheduling heuristic, not a completion date or scientific acceptance.'}


def ranked_species(root, species, *, now=None, write=False):
    root = Path(root); now = time.time() if now is None else now
    rows = {r['species_id']: r for r in readj(root/'registry/batch_plan_v5.json')['species']}
    jobs = {r['job_id']: r for r in readj(root/'registry/jobs.json')}
    policy_path = root/'config/rollout_20261015.json'
    enqueued = readj(policy_path).get('created_epoch', now) if policy_path.is_file() else now
    reports = []
    for sid in species:
        report = estimate(root, rows[sid], jobs)
        latest = root/'operations/packed/latest'/(sid+'.json')
        last = readj(latest).get('returned_epoch', enqueued) if latest.is_file() else enqueued
        report['waiting_hours'] = max(0, (now-last)/3600)
        # A finite starvation guard; no preemption of active ORCA processes.
        report['priority_group'] = (0 if sid in BENCHMARKS else
                                    1 if report['waiting_hours'] >= 48 else 2)
        reports.append(report)
    reports.sort(key=lambda r: (r['priority_group'], -r['remaining_serial_hours_estimate'], r['species_id']))
    if write:
        writej(root/'operations/priority_latest.json', {'updated_epoch': now, 'branches': reports})
    return [r['species_id'] for r in reports]
