"""Run at most one GOAT job in an ASCII scratch directory, retaining provenance.

No work is submitted in bulk. --dry-run does not create files or launch ORCA.
Scratch is deliberately retained, never automatically deleted or retried.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import uuid

try:
    from . import goat_job as runner
except ImportError:
    import goat_job as runner

ROOT = Path(__file__).resolve().parents[1]
IMMUTABLE_FILES = ('goat.inp', 'start.xyz', 'cluster.mol', 'atom_map.json', 'isotope_map.json', 'provenance.json')


def ascii_path(value: Path) -> Path:
    path = Path(value).resolve()
    if not str(path).isascii() or any(c.isspace() for c in str(path)):
        raise ValueError('Scratch root must resolve to an ASCII path without spaces')
    if path == Path(path.anchor):
        raise ValueError('Use a dedicated scratch directory, not a drive/filesystem root')
    if path.exists() and not path.is_dir():
        raise ValueError('Scratch root is not a directory')
    return path


def _hashes(directory: Path) -> dict:
    values = {}
    for name in IMMUTABLE_FILES:
        path = directory / name
        if not path.is_file() or path.is_symlink() or not path.resolve().is_relative_to(directory.resolve()):
            raise ValueError('Missing/unsafe immutable input: ' + name)
        values[name] = runner.digest(path)
    return values


def _new_copy(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open('rb') as src, target.open('xb') as dst:
        shutil.copyfileobj(src, dst, 1024 * 1024)
    if runner.digest(source) != runner.digest(target):
        raise ValueError('Stage-out copy hash mismatch: ' + str(target))


def _receipt(stage_job: Path, job_id: str, source_hashes: dict) -> dict:
    path = stage_job / 'receipt.json'
    if not path.is_file():
        raise ValueError('Staged runner did not produce receipt.json; scratch retained')
    receipt = runner.read_json(path)
    if receipt.get('job_id', receipt.get('id')) != job_id or receipt.get('status') not in (
            'completed', 'COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW', 'FAILED'):
        raise ValueError('Staged receipt is not a terminal receipt for this job')
    inputs = receipt.get('input_sha256')
    outputs = receipt.get('output_sha256')
    if not isinstance(inputs, dict) or not isinstance(outputs, dict):
        raise ValueError('Staged receipt lacks input/output hashes')
    for name in runner.INPUT_HASH_FIELDS:
        if inputs.get(name) != source_hashes[name]:
            raise ValueError('Staged receipt input binding mismatch: ' + name)
    if inputs.get('provenance.json') != source_hashes['provenance.json']:
        raise ValueError('Staged receipt provenance binding mismatch')
    for name, expected in outputs.items():
        value = (stage_job / name).resolve()
        if not value.is_relative_to(stage_job.resolve()) or not value.is_file() or runner.digest(value) != expected:
            raise ValueError('Staged receipt output hash mismatch: ' + name)
    if receipt.get('status') in ('completed', 'COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW') and any(
            name not in outputs for name in ('goat.out', 'goat.finalensemble.xyz')):
        raise ValueError('Completed staged receipt does not bind output and ensemble')
    return receipt


def run_one(root: Path, job_id: str, orca: Path, xtb: Path, scratch_root: Path,
            acknowledge_models: bool = False, dry_run: bool = False,
            timeout: float | None = None) -> dict:
    root = Path(root).resolve()
    scratch_root = ascii_path(scratch_root)
    if not acknowledge_models:
        raise ValueError('Supply --acknowledge-models; this is not scientific validation')
    if timeout is not None and (not math.isfinite(timeout) or timeout <= 0):
        raise ValueError('Timeout must be finite and positive')
    job = runner.validate_job(root, job_id)
    original = job['directory']
    if (original / 'portable_stageout.json').exists():
        raise ValueError('Existing portable stage-out record; no automatic retry')
    original_hashes = _hashes(original)
    runner.check_allocation(job['nprocs'], job['maxcore_mb'])
    preflight = runner.preflight(orca, xtb, job['nprocs'])
    if not preflight['passed']:
        raise ValueError('Preflight blocked: ' + '; '.join(preflight['reasons']))
    identity = json.dumps({'root': str(root), 'job_id': job_id, 'inputs': original_hashes,
                           'manifest': job['manifest_sha256']}, sort_keys=True).encode()
    stage_root = scratch_root / ('work2goat_' + hashlib.sha256(identity).hexdigest()[:20])
    stage_job = stage_root / 'job'
    if stage_root.exists():
        raise ValueError('Scratch leaf already exists; inspect it, no automatic overwrite/retry: ' + str(stage_root))
    if stage_root.is_relative_to(original) or original.is_relative_to(stage_root):
        raise ValueError('Scratch and original job directories must be separate')
    plan = {'job_id': job_id, 'original_job_directory': str(original), 'stage_root': str(stage_root),
            'stage_job_directory': str(stage_job), 'original_manifest_sha256': job['manifest_sha256'],
            'original_input_sha256': original_hashes, 'preflight': preflight,
            'runtime_goat_smoke': 'UNVALIDATED_NOT_EXECUTED', 'scratch_cleanup': 'RETAIN_FOR_AUDIT_NO_AUTOMATIC_DELETE'}
    if dry_run:
        return {**plan, 'status': 'DRY_RUN_NO_EXECUTION_NO_FILES_CREATED', 'quantum_jobs_launched': 0}
    lock = original / '.goat-run.lock'
    token = uuid.uuid4().hex
    with lock.open('x', encoding='utf-8') as stream:
        json.dump({'token': token, 'pid': os.getpid(), 'started_utc': runner.utc(), 'portable_stage_root': str(stage_root)}, stream)
    try:
        fresh = runner.validate_job(root, job_id, allow_own_lock=True)
        if fresh['manifest_sha256'] != job['manifest_sha256'] or _hashes(original) != original_hashes:
            raise ValueError('Original inputs/manifest changed before staging')
        stage_root.mkdir(parents=True, exist_ok=False)
        stage_job.mkdir()
        for name in IMMUTABLE_FILES:
            _new_copy(original / name, stage_job / name)
        mini_row = dict(job['row'])
        mini_row['job_directory'] = 'job'
        (stage_root / 'registry').mkdir()
        runner.atomic_json(stage_root / 'registry' / 'jobs.json', [mini_row])
        mini_manifest_sha = runner.digest(stage_root / 'registry' / 'jobs.json')
        failure = None
        try:
            runner.run_one(stage_root, job_id, orca, xtb, acknowledge_models=True, dry_run=False, timeout=timeout)
        except Exception as exc:
            # Preserve valid terminal failure receipts/logs too, never invent one.
            failure = exc
        receipt = _receipt(stage_job, job_id, original_hashes)
        if receipt.get('manifest_sha256') != mini_manifest_sha:
            raise ValueError('Staged receipt is not bound to its actual mini-registry')
        if _hashes(original) != original_hashes or _hashes(stage_job) != original_hashes:
            raise ValueError('Original/staged immutable inputs changed; no stage-out, scratch retained')
        if runner.digest(root / 'registry' / 'jobs.json') != job['manifest_sha256']:
            raise ValueError('Original manifest changed; no stage-out, scratch retained')
        copies = []
        for path in sorted(stage_job.rglob('*')):
            if path.is_symlink() or not path.resolve().is_relative_to(stage_job.resolve()):
                raise ValueError('Unsafe staged output path')
            if not path.is_file():
                continue
            relative = path.relative_to(stage_job)
            if relative.as_posix() in IMMUTABLE_FILES or relative.as_posix() == 'receipt.json':
                continue
            if path.name == '.goat-run.lock':
                raise ValueError('Staged runner still owns its lock; no stage-out')
            destination = original / relative
            if destination.exists():
                raise ValueError('Refusing to overwrite existing output: ' + str(destination))
            copies.append((path, destination, relative.as_posix()))
        if (original / 'receipt.json').exists() or (original / 'portable_stageout.json').exists():
            raise ValueError('Existing original receipt or stage-out record')
        for source, destination, _ in copies:
            _new_copy(source, destination)
        if _hashes(original) != original_hashes or _hashes(stage_job) != original_hashes:
            raise ValueError('Immutable inputs changed during stage-out; receipt not committed')
        if _receipt(stage_job, job_id, original_hashes) != receipt:
            raise ValueError('Staged receipt changed during stage-out')
        for name, expected in receipt['output_sha256'].items():
            if runner.digest(original / name) != expected:
                raise ValueError('Copied output differs from terminal receipt: ' + name)
        evidence = {**plan, 'status': 'STAGEOUT_EVIDENCE_RECEIPT_COMMITTED_LAST',
                    'staged_manifest_sha256': mini_manifest_sha,
                    'staged_receipt_sha256': runner.digest(stage_job / 'receipt.json'),
                    'copied_output_sha256': {relative: runner.digest(destination) for _, destination, relative in copies},
                    'stageout_utc': runner.utc(), 'receipt_job_directory_is_actual_ascii_scratch': True,
                    'notice': 'Original receipt is copied byte-for-byte; its manifest hash belongs to the staged mini-registry.'}
        # This is new evidence, not a mutation of the hash-bound runner receipt.
        with (original / 'portable_stageout.json').open('x', encoding='utf-8') as stream:
            json.dump(evidence, stream, indent=2, ensure_ascii=False)
        _new_copy(stage_job / 'receipt.json', original / 'receipt.json')
        return {**plan, 'status': receipt['status'], 'quantum_jobs_launched': receipt.get('quantum_jobs_launched', 0),
                'receipt': str(original / 'receipt.json'), 'stageout_record': str(original / 'portable_stageout.json'),
                'runner_exception': str(failure) if failure else None}
    finally:
        if lock.is_file() and runner.read_json(lock).get('token') == token:
            lock.unlink()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('preflight', 'run'):
        command = sub.add_parser(name)
        command.add_argument('--orca', type=Path, required=True)
        command.add_argument('--xtb', type=Path, required=True)
        command.add_argument('--scratch-root', type=Path, required=True)
        if name == 'run':
            command.add_argument('--job-id', required=True)
            command.add_argument('--acknowledge-models', action='store_true')
            command.add_argument('--dry-run', action='store_true')
            command.add_argument('--timeout', type=float)
    args = parser.parse_args(argv)
    try:
        if args.command == 'preflight':
            result = {**runner.preflight(args.orca, args.xtb), 'scratch_root': str(ascii_path(args.scratch_root))}
            code = 0 if result['passed'] else 2
        else:
            result = run_one(ROOT, args.job_id, args.orca, args.xtb, args.scratch_root,
                             args.acknowledge_models, args.dry_run, args.timeout)
            code = 3 if result['status'] == 'FAILED' else 0
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return code
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print('STOP: ' + str(exc), file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
