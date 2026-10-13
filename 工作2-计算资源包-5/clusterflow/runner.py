"""Explicit one-job quantum runner with method, resource and path gates."""
from __future__ import annotations
import math
import os
from pathlib import Path
import re
import shutil
import hashlib
import json
import threading
import uuid
from contextlib import ExitStack, contextmanager
from .common import config,hash_protocol,inside,readj,sha
from .orcaio import run_job, _input_hashes

IMMUTABLE_FILES = ('job.inp', 'job.json', 'input.xyz', 'input_structure.mol', 'atom_map.json')
_ENV_LOCK = threading.Lock()


def _safe_file(directory, name, *, existing=True):
    relative = Path(name)
    if (not isinstance(name, str) or not name or relative.is_absolute()
            or '..' in relative.parts or relative.as_posix() != name):
        raise ValueError('Unsafe job-relative file path: ' + str(name))
    path = directory / relative
    if not path.resolve().is_relative_to(directory.resolve()):
        raise ValueError('File path escapes job directory: ' + name)
    for item in (path, *path.parents):
        if item == directory:
            break
        if item.is_symlink() or (hasattr(item, 'is_junction') and item.is_junction()):
            raise ValueError('Linked job file/directory is not permitted: ' + name)
    if existing and not path.is_file():
        raise ValueError('Missing regular job file: ' + name)
    return path


def _inputs(directory):
    for name in IMMUTABLE_FILES:
        _safe_file(directory, name)
    values = _input_hashes(directory)
    for name in values:
        _safe_file(directory, name)
        if name in ('receipt.json', 'portable_stageout.json') or name.startswith(('.', 'attempts/')):
            raise ValueError('Input reference uses a reserved output name: ' + name)
    return values


@contextmanager
def _lock(directory, name):
    path = directory / name
    token = uuid.uuid4().hex
    with path.open('x', encoding='utf-8') as stream:
        json.dump({'token': token, 'pid': os.getpid()}, stream)
    try:
        yield
    finally:
        try:
            if readj(path).get('token') == token:
                path.unlink()
        except (OSError, ValueError, TypeError):
            pass  # Never remove another owner's lock or an unreadable lock.


def _copy_new(source, destination, expected):
    """Exclusive creation: an output appearing after our check is not overwritten."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open('rb') as incoming, destination.open('xb') as outgoing:
        shutil.copyfileobj(incoming, outgoing, 1024 * 1024)
    if sha(source) != expected or sha(destination) != expected:
        raise ValueError('Copy hash mismatch; partial stage-out retained: ' + str(destination))


def _commit_receipt(source, destination, expected, validate):
    """Publish a verified complete receipt atomically, never replacing a path.

    A same-directory hard link has no overwrite semantics on Windows/POSIX.
    Filesystems without hard-link support fail closed, retaining scratch.
    """
    temporary=destination.with_name('.receipt.'+uuid.uuid4().hex+'.tmp')
    try:
        _copy_new(source,temporary,expected)
        validate()
        os.link(temporary,destination)
    finally:
        if temporary.exists():temporary.unlink()


def _terminal_receipt(directory, metadata, inputs, executable_hash):
    path = _safe_file(directory, 'receipt.json')
    receipt = readj(path)
    if (receipt.get('id') != metadata['id'] or receipt.get('stage') != metadata['stage']
            or receipt.get('status') not in ('completed', 'failed')):
        raise ValueError('Receipt is not a terminal receipt for this job')
    if receipt.get('input_sha256') != inputs:
        raise ValueError('Receipt input binding mismatch')
    if receipt.get('orca_executable_sha256') != executable_hash:
        raise ValueError('Receipt executable binding mismatch')
    outputs = receipt.get('output_sha256')
    if not isinstance(outputs, dict) or 'job.out' not in outputs:
        raise ValueError('Receipt lacks output hashes')
    for name, expected in outputs.items():
        if name in inputs or name in ('receipt.json', 'portable_stageout.json') or name.startswith('.'):
            raise ValueError('Receipt output uses a reserved name')
        if sha(_safe_file(directory, name)) != expected:
            raise ValueError('Receipt output hash mismatch: ' + name)
    if receipt['status'] == 'completed':
        for name in metadata.get('required_outputs', []):
            if name not in outputs or not _safe_file(directory, name).stat().st_size:
                raise ValueError('Completed receipt lacks required output: ' + name)
    return receipt


def _existing_outputs(directory, inputs):
    """Recognize previous runner products, not unrelated user notes."""
    names = set()
    for path in directory.iterdir():
        if path.name not in inputs and (path.name.startswith('job.') or path.name in (
                'receipt.json', 'propertyjob.property.txt', 'portable_stageout.json')):
            names.add(path.name)
    for record, field in (('receipt.json', 'output_sha256'), ('portable_stageout.json', 'copied_output_sha256')):
        if (directory / record).is_file():
            values = readj(directory / record).get(field, {})
            if not isinstance(values, dict):
                raise ValueError('Malformed prior output inventory')
            names.update(values)
    for name in names:
        if name in inputs:
            raise ValueError('Prior output inventory overlaps immutable inputs')
        _safe_file(directory, name, existing=False)
    return names


def _archive_outputs(directory, names):
    existing = [name for name in sorted(names) if (directory / name).exists()]
    if not existing:
        return
    archive = directory / 'attempts' / ('before_explicit_retry_' + uuid.uuid4().hex)
    _safe_file(directory, archive.relative_to(directory).as_posix(), existing=False)
    archive.mkdir(parents=True, exist_ok=False)
    for name in existing:
        source = _safe_file(directory, name)
        target = archive / name
        target.parent.mkdir(parents=True, exist_ok=True)
        source.rename(target)

def binary_versions(path):
    versions=set(); tail=b''
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):
            data=tail+b
            versions.update(x.decode() for x in re.findall(rb'Program Version\s+(\d+\.\d+\.\d+)(?![\d.])',data))
            tail=data[-128:]
    return sorted(versions)

def static_preflight(root,orca):
    p=config(root); exe=Path(orca).resolve(); issues=[]
    versions=binary_versions(exe) if exe.is_file() else []
    if not exe.is_file():issues.append('ORCA executable not found')
    if versions!=[p['orca_version']]:issues.append('Exact binary version marker not confirmed')
    mpi=shutil.which('mpiexec') or shutil.which('mpirun')
    if not mpi:issues.append('MPI launcher not in PATH')
    return {'passed':not issues,'issues':issues,'orca':str(exe),'binary_versions':versions,
            'orca_sha256':sha(exe) if exe.is_file() else None,'mpi_launcher':mpi,
            'quantum_jobs_launched':0,'runtime_validation':'NOT_PERFORMED_BY_STATIC_CHECK',
            'external_xtb_required_for_DFT':False,'external_xtb_required_for_GOAT':True}

def allocation(nprocs,maxcore):
    if type(nprocs) is not int or type(maxcore) is not int or min(nprocs,maxcore)<1:
        raise ValueError('Invalid resource specification')
    if int(os.environ.get('SLURM_NTASKS','1'))!=1:
        raise ValueError('Request one driver task, not multiple ORCA drivers')
    cores=int(os.environ.get('SLURM_CPUS_PER_TASK',nprocs))
    if cores<nprocs:raise ValueError('Insufficient allocated CPUs')
    mb=os.environ.get('SLURM_MEM_PER_NODE')
    if not mb and os.environ.get('SLURM_MEM_PER_CPU'):mb=int(os.environ['SLURM_MEM_PER_CPU'])*cores
    if mb and nprocs*maxcore>0.75*int(mb):raise ValueError('Insufficient allocated memory headroom')

def validate_input(root,relative):
    root=Path(root).resolve(); d=inside(root,relative)
    if not d.is_relative_to(root/'jobs'):raise ValueError('Only downstream jobs/ paths are executable here')
    m=readj(d/'job.json'); p=config(root)
    if m.get('protocol_hash')!=hash_protocol(p):raise ValueError('Protocol changed since input preparation')
    if m.get('stage') not in ['prescreen','optfreq','shielding','j_PBE_pcJ3','j_PBE0_pcJ2']:
        raise ValueError('Stage is not in the approved no-revDSD workflow')
    if sha(d/'job.inp')!=m.get('input_sha256'):raise ValueError('Input content hash changed')
    if 'atom_map_sha256' in m and sha(d/'atom_map.json')!=m['atom_map_sha256']:
        raise ValueError('Cluster atom map changed')
    allocation(m['nprocs'],m['maxcore_mb'])
    return d,m

def run_one(root,relative,orca,*,acknowledge=False,dry_run=False,timeout=None,retry=False,scratch_root=None):
    if not acknowledge:raise ValueError('Explicit --acknowledge-models required')
    if timeout is not None and (not math.isfinite(timeout) or timeout<=0):raise ValueError('Invalid timeout')
    root=Path(root).resolve()
    d,m=validate_input(root,relative)
    original_inputs=_inputs(d)
    check=static_preflight(root,orca)
    if not check['passed']:raise ValueError('Preflight failed: '+'; '.join(check['issues']))
    plan={'job':relative,'stage':m['stage'],'command':[check['orca'],'job.inp'],
          'nprocs':m['nprocs'],'maxcore_mb_per_rank':m['maxcore_mb'],'preflight':check}
    stage_dir=None
    if scratch_root:
        scratch=Path(scratch_root).resolve()
        if not str(scratch).isascii() or any(c.isspace() for c in str(scratch)):
            raise ValueError('ORCA scratch path must be ASCII without whitespace')
        if scratch == Path(scratch.anchor) or (scratch.exists() and not scratch.is_dir()):
            raise ValueError('Use a dedicated scratch directory, not a root or regular file')
        identity=json.dumps({'job_directory':str(d),'inputs':original_inputs,
                             'executable_sha256':check['orca_sha256']},sort_keys=True).encode()
        key=hashlib.sha256(identity).hexdigest()[:24]
        stage_dir=scratch/('work2_'+key+('_retry_'+uuid.uuid4().hex[:12] if retry else ''))
        if stage_dir.is_relative_to(d) or d.is_relative_to(stage_dir):
            raise ValueError('Scratch and original job directories must be separate')
        plan['ascii_scratch_directory']=str(stage_dir)
    elif not str(d).isascii() or any(c.isspace() for c in str(d)):
        raise ValueError('Keep the Chinese package name, but use --scratch-root with an ASCII no-space directory for the ORCA process')
    if dry_run:return {**plan,'status':'DRY_RUN_NO_EXECUTION','quantum_jobs_launched':0}
    # The original-job lock covers validation, execution and receipt-last commit.
    # A staged run also owns the low-level lock in the ORIGINAL directory, so a
    # direct orcaio.run_job cannot race a run using a different scratch root.
    with ExitStack() as locks:
        locks.enter_context(_lock(d,'.clusterflow-run.lock'))
        if stage_dir is not None:
            locks.enter_context(_lock(d,'.orca-run.lock'))
        def unchanged():
            fresh_d,fresh=validate_input(root,relative)
            if fresh_d!=d or fresh!=m or _inputs(d)!=original_inputs:
                raise ValueError('Original immutable inputs changed; no valid stage-out')
            if sha(check['orca'])!=check['orca_sha256']:
                raise ValueError('ORCA executable changed since preflight')
        unchanged()
        previous_outputs=_existing_outputs(d,original_inputs)
        if (d/'receipt.json').exists() and not retry:
            previous=_terminal_receipt(d,m,original_inputs,check['orca_sha256'])
            if previous['status']!='completed':
                raise ValueError('Prior incomplete/failed run requires review and explicit --retry')
            return {**previous,'skipped':True,'resume_reason':'matching_inputs_and_outputs'}
        if not (d/'receipt.json').exists() and previous_outputs:
            raise ValueError('Unreceipted existing output requires a separate import audit')
        actual=d
        if stage_dir is not None:
            # Partial scratch work is never silently resumed or overwritten.
            stage_dir.mkdir(parents=True,exist_ok=False)
            for name,expected in original_inputs.items():
                _copy_new(_safe_file(d,name),_safe_file(stage_dir,name,existing=False),expected)
            if _inputs(stage_dir)!=original_inputs:
                raise ValueError('Staged immutable inputs differ')
            unchanged()
            actual=stage_dir
        # In-place retry archival belongs to run_job while it owns its own
        # .orca-run.lock. Never move in-place products before that lock exists.
        if not _ENV_LOCK.acquire(blocking=False):
            raise ValueError('Another in-process ORCA runner owns the execution environment')
        oldpath=os.environ.get('PATH'); oldomp=os.environ.get('OMP_NUM_THREADS')
        try:
            os.environ['PATH']=str(Path(check['orca']).parent)+os.pathsep+(oldpath or '')
            os.environ['OMP_NUM_THREADS']='1'
            returned=run_job(actual,Path(check['orca']),timeout_seconds=timeout,resume=not retry)
        finally:
            if oldpath is None:os.environ.pop('PATH',None)
            else:os.environ['PATH']=oldpath
            if oldomp is None:os.environ.pop('OMP_NUM_THREADS',None)
            else:os.environ['OMP_NUM_THREADS']=oldomp
            _ENV_LOCK.release()
        unchanged()
        if _inputs(actual)!=original_inputs:
            raise ValueError('Staged immutable inputs changed; scratch retained')
        receipt=_terminal_receipt(actual,m,original_inputs,check['orca_sha256'])
        if returned!=receipt:
            raise ValueError('Returned receipt differs from the on-disk terminal receipt')
        if actual==d:return receipt
        receipt_hash=sha(actual/'receipt.json')
        copies={}
        for path in sorted(actual.rglob('*')):
            name=path.relative_to(actual).as_posix()
            # Reject symlinks/junctions including directories before recursion.
            _safe_file(actual,name,existing=False)
            if path.is_dir():continue
            if name in original_inputs or name=='receipt.json':continue
            if name in ('.orca-run.lock','.clusterflow-run.lock','portable_stageout.json'):
                raise ValueError('Unexpected staged lock/stage-out record')
            copies[name]=sha(_safe_file(actual,name))
        if retry:
            # Archive all old products together, not only names emitted again.
            # This prevents stale old property files from contaminating a retry.
            _archive_outputs(d,previous_outputs)
        for name in copies:
            if _safe_file(d,name,existing=False).exists():
                raise ValueError('Refusing to overwrite existing destination output: '+name)
        if (d/'receipt.json').exists() or (d/'portable_stageout.json').exists():
            raise ValueError('Existing destination receipt/stage-out record')
        for name,expected in copies.items():
            _copy_new(_safe_file(actual,name),_safe_file(d,name,existing=False),expected)
        unchanged()
        if _inputs(actual)!=original_inputs or sha(actual/'receipt.json')!=receipt_hash:
            raise ValueError('Staged inputs/receipt changed during stage-out')
        if _terminal_receipt(actual,m,original_inputs,check['orca_sha256'])!=receipt:
            raise ValueError('Staged receipt changed during stage-out')
        for name,expected in copies.items():
            if sha(_safe_file(actual,name))!=expected or sha(_safe_file(d,name))!=expected:
                raise ValueError('Copied output hash mismatch: '+name)
        evidence={**plan,'original_input_sha256':original_inputs,'copied_output_sha256':copies,
                  'staged_receipt_sha256':receipt_hash,'receipt_copied_byte_for_byte':True,
                  'status':'STAGEOUT_EVIDENCE_RECEIPT_COMMITTED_LAST',
                  'scratch_cleanup':'RETAIN_FOR_AUDIT_NO_AUTOMATIC_DELETE'}
        with (d/'portable_stageout.json').open('x',encoding='utf-8') as stream:
            json.dump(evidence,stream,ensure_ascii=False,indent=2)
        # Receipt is the final commit marker. No input/output is silently fixed.
        def final_check():
            unchanged()
            if _inputs(actual)!=original_inputs or sha(actual/'receipt.json')!=receipt_hash:
                raise ValueError('Staged inputs/receipt changed during receipt commit')
            for name,expected in copies.items():
                if sha(_safe_file(actual,name))!=expected or sha(_safe_file(d,name))!=expected:
                    raise ValueError('Copied output hash mismatch during receipt commit: '+name)
        _commit_receipt(actual/'receipt.json',d/'receipt.json',receipt_hash,final_check)
        return receipt
