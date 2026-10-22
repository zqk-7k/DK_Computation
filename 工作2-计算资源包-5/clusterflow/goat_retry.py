"""Bind an operator-reviewed node-failure retry to retained prior evidence."""
import re
from pathlib import Path
from .common import readj, sha, inside


def retry_token(root, directory, job_id, manifest_sha256, inputs):
    path = Path(directory)/'goat_retry_authorization.json'
    if not path.exists():
        return None
    record = readj(path)
    if (record.get('status') not in ('PREPARED_VERIFIED_NODE_FAIL', 'PREPARED_VERIFIED_MPI_ABORT')
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
    if proof.get('job_id') != job_id:
        raise ValueError('Retry proof belongs to another job')
    stage = Path(proof['stage'])
    previous = stage/'receipt.json'
    if sha(previous) != proof['receipt_sha256']:
        raise ValueError('Prior attempt changed after retry review')
    if record['status'] == 'PREPARED_VERIFIED_NODE_FAIL':
        if proof.get('scheduler_state') != 'NODE_FAIL' or readj(previous).get('status') != 'RUNNING':
            raise ValueError('Prior terminal result requires review instead of NODE_FAIL retry')
    else:
        exited = proof.get('producer_exit_proof', {})
        if (proof.get('failure_class') != 'MPI_ABORT_BEFORE_FIRST_GLOBAL_ITERATION'
                or exited.get('matches') != [] or exited.get('failures') != []
                or readj(previous).get('status') != 'FAILED'
                or (stage/'goat.finalensemble.xyz').exists()
                or sha(stage/'goat.stderr') != proof.get('mpi_stderr_sha256')):
            raise ValueError('MPI retry requires retained terminal failure and proven producer exit')
    return record['attempt_token']
