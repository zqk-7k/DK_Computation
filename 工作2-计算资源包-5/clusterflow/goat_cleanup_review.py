"""Narrow technical review of ORCA success followed by sandbox cleanup failure.

The original FAILED receipt is immutable. This does not approve structures,
sampling convergence, or downstream scientific gates.
"""
import re
from pathlib import Path
from .common import readj, sha


def evidence(directory, receipt_path):
    directory = Path(directory)
    receipt = readj(receipt_path)
    if (receipt.get('status') != 'FAILED' or receipt.get('returncode') != 1
            or receipt.get('errors') != ['ORCA nonzero or absent exit code']):
        raise ValueError('Not an isolated wrapper exit-code failure')
    hashes = {}
    for name in ('goat.out', 'goat.finalensemble.xyz', 'goat.stderr'):
        path = directory / name
        if not path.is_file() or path.is_symlink():
            raise ValueError('Missing or unsafe cleanup-review evidence')
        digest = sha(path)
        if receipt.get('output_sha256', {}).get(name) != digest:
            raise ValueError('Cleanup-review output binding changed: ' + name)
        hashes[name] = digest
    output = (directory / 'goat.out').read_text(errors='replace')
    if ('ORCA TERMINATED NORMALLY' not in output
            or re.search(r'GOAT ERROR|FATAL ERROR|error termination', output, re.I)
            or set(re.findall(r'Program Version\s+(\d+\.\d+\.\d+)', output)) != {'6.1.0'}):
        raise ValueError('ORCA calculation did not pass the normal-output gate')
    cleanup = 0
    for line in (directory / 'goat.stderr').read_text(errors='replace').splitlines():
        if not line.strip():
            continue
        if (line.startswith('rm: cannot remove ') and 'Read-only file system' in line
                and re.search(r'/var/tmp/work2-runtime-[^/]+/orca610_sandbox_', line)):
            cleanup += 1
            continue
        if re.fullmatch(r'INFO:\s+(Starting build\.\.\.|Creating sandbox directory\.\.\.)', line):
            continue
        if re.fullmatch(r'INFO:\s+Verifying bootstrap image /public/home/qikai/apprepo/orca/6\.1\.0-glibc229/app/orca-openmpi418\.sif', line):
            continue
        if re.fullmatch(r'INFO:\s+Build complete: /var/tmp/work2-runtime-[^/]+/orca610_sandbox_[^/]+\.partial', line):
            continue
        raise ValueError('Unreviewed stderr diagnostic: ' + line[:200])
    if not cleanup:
        raise ValueError('No readonly sandbox-cleanup evidence')
    return {'source_receipt_sha256': sha(receipt_path), 'outputs': hashes,
            'cleanup_error_lines': cleanup}


def validate(directory, receipt_path):
    directory = Path(directory)
    record = readj(directory / 'goat_cleanup_review.json')
    if record.get('status') != 'NORMAL_ORCA_OUTPUT_WITH_FAILED_SANDBOX_CLEANUP':
        raise ValueError('No approved technical cleanup review')
    actual = evidence(directory, receipt_path)
    if record.get('evidence') != actual:
        raise ValueError('Technical cleanup-review evidence changed')
    return sha(directory / 'goat_cleanup_review.json')
