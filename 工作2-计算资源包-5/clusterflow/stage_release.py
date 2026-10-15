"""Release runtime stages from actual validated output, not Slurm success."""
from pathlib import Path
import time

from .common import readj, writej, sha, config


def audit(root):
    root = Path(root)
    evidence = {}; errors = []
    for row in readj(root/'registry/jobs.json'):
        receipt = root/row['job_directory']/'receipt.json'
        if not receipt.is_file(): continue
        if str(readj(receipt).get('status', '')).upper() not in {'COMPLETED', 'COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW'}: continue
        try:
            # Polls with no finished chemistry need only the standard library.
            # Loading NumPy/RDKit from shared storage can be slow on login nodes.
            from .pipeline import _goat_source
            _goat_source(root, row['job_id'])
            evidence['goat'] = {'receipt': str(receipt.relative_to(root)), 'sha256': sha(receipt), 'job_id': row['job_id']}
            break
        except (ValueError, OSError, KeyError) as exc:
            errors.append({'receipt': str(receipt.relative_to(root)), 'reason': str(exc)})
    # Only the unrestricted original pilot can validate new DFT stages.
    pilot = readj(root/'config/pilot_scope_20261014.json')['species_ids']
    for sid in pilot:
        for workspace in (root/'automation'/sid).glob('wave_*/workspace'):
            for meta_path in (workspace/'jobs').glob('*/*/*/*/*/job.json'):
                meta = readj(meta_path); stage = meta.get('stage')
                receipt = meta_path.with_name('receipt.json')
                if stage in evidence or not receipt.is_file(): continue
                data = readj(receipt)
                if str(data.get('status', '')).upper() != 'COMPLETED' or data.get('reuse'): continue
                try:
                    from .gates import inspect_job
                    checked = inspect_job(meta_path.parent, config(workspace))
                    if checked['accepted']:
                        evidence[stage] = {'receipt': str(receipt.relative_to(root)), 'sha256': sha(receipt),
                                           'species_id': meta['species_id']}
                except (ValueError, OSError, KeyError) as exc:
                    errors.append({'receipt': str(receipt.relative_to(root)), 'reason': str(exc)})
    report = {'updated_epoch': time.time(), 'allowed_stages': ['goat'] + sorted(set(evidence)-{'goat'}),
              'validated_goat': 'goat' in evidence, 'evidence': evidence, 'errors': errors,
              'scientific_project_accepted': False,
              'notice': 'Stage execution readiness only. Per-job scientific gates remain mandatory.'}
    writej(root/'operations/stage_release.json', report)
    return report


def check(root, stage):
    if stage == 'goat': return
    from .cache import package_root
    from .scheduling import DeferredWork
    root = package_root(root); path = root/'operations/stage_release.json'
    report = readj(path) if path.is_file() else {}
    proof = report.get('evidence', {}).get(stage)
    valid = stage in report.get('allowed_stages', []) and proof
    if valid:
        receipt = (root/proof['receipt']).resolve()
        valid = receipt.is_relative_to(root) and receipt.is_file() and sha(receipt) == proof['sha256']
    if not valid:
        raise DeferredWork('WAITING_STAGE_VALIDATION', stage=stage,
                           notice='Continue independent GOAT work; representative runtime stage not yet validated.')
