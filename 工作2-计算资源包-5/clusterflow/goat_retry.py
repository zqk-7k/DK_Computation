"""Bind an operator-reviewed node-failure retry to retained prior evidence."""
import re
from pathlib import Path
from .common import readj, sha, inside


def retry_token(root, directory, job_id, manifest_sha256, inputs):
    path = Path(directory)/'goat_retry_authorization.json'
    if not path.exists():
        return None
    record = readj(path)
    if (record.get('status') != 'PREPARED_VERIFIED_NODE_FAIL'
            or record.get('job_id') != job_id
            or record.get('manifest_sha256') != manifest_sha256
            or record.get('inputs') != inputs
            or not re.fullmatch(r'[0-9a-f]{32}', record.get('attempt_token', ''))):
        raise ValueError('Retry authorization identity changed')
    proof_path = inside(root, record['proof_path'])
    if not proof_path.is_relative_to(Path(directory).resolve()/'recovery_history'):
        raise ValueError('Retry proof must be retained with this job')
    if sha(proof_path) != record.get('proof_sha256'):
        raise ValueError('Retry proof hash changed')
    proof = readj(proof_path)
    if proof.get('job_id') != job_id or proof.get('scheduler_state') != 'NODE_FAIL':
        raise ValueError('Retry requires this job\'s verified NODE_FAIL allocation')
    previous = Path(proof['stage'])/'receipt.json'
    if sha(previous) != proof['receipt_sha256']:
        raise ValueError('Prior attempt changed after retry review')
    if readj(previous).get('status') != 'RUNNING':
        raise ValueError('Prior terminal result requires review instead of retry')
    return record['attempt_token']
