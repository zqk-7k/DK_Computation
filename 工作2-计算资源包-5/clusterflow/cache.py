"""Exact-input, source-bound DFT reuse across immutable waves and TMS models.

A reuse receipt is explicitly DERIVED, never claimed to be a new ORCA run.
The original execution and all its hashes must remain available for validation.
"""
import hashlib
import json
import os
from pathlib import Path
from .common import readj,writej,sha,config


def package_root(root):
    root=Path(root).resolve()
    for parent in reversed([root,*root.parents]):
        if (parent/'registry/batch_plan_v5.json').is_file() and (parent/'config/protocol.json').is_file():return parent
    return root


def runtime_identity():
    path=os.environ.get('WORK2_RUNTIME_EVIDENCE')
    return sha(path) if path and Path(path).is_file() else None


def key(directory,executable_hash,runtime_evidence):
    from .orcaio import _input_hashes
    inputs=_input_hashes(Path(directory));inputs.pop('job.json',None)
    metadata=readj(Path(directory)/'job.json')
    fields=('stage','method_fingerprint','species_id','scenario_id','charge','multiplicity','protocol_hash')
    payload={'inputs':inputs,'model':{k:metadata.get(k) for k in fields},
             'executable_sha256':executable_hash,'runtime_evidence_sha256':runtime_evidence}
    return hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()


def publish(root,directory):
    if not config(root).get('execution',{}).get('verified_exact_input_reuse',False):return
    from .gates import inspect_job
    directory=Path(directory);receipt=readj(directory/'receipt.json')
    if receipt.get('reuse') or not inspect_job(directory,config(root))['accepted']:return
    exe=receipt['orca_executable_sha256'];evidence=receipt.get('runtime_evidence_sha256')
    digest=key(directory,exe,evidence)
    target=package_root(root)/'operations/cache'/(digest+'.json')
    if not target.exists():
        writej(target,{'source':str(directory.resolve()),'source_receipt_sha256':sha(directory/'receipt.json'),
                      'source_protocol':str((Path(root)/'config/protocol.json').resolve()),'key':digest})


def validate_reuse(directory,receipt,p):
    from .gates import inspect_job
    directory=Path(directory).resolve();reuse=receipt['reuse'];source=Path(reuse['source']).resolve()
    if source==directory or not source.is_relative_to(package_root(directory)):
        raise ValueError('Reuse source must remain inside this production package')
    if sha(source/'receipt.json')!=reuse['source_receipt_sha256']:raise ValueError('Reuse source receipt changed')
    old=readj(source/'receipt.json')
    if old.get('reuse'):raise ValueError('Chained reused receipts are not original quantum executions')
    if old.get('status')!='completed' or not inspect_job(source,p)['accepted']:raise ValueError('Reuse source no longer scientifically accepted')
    source_key=key(source,old['orca_executable_sha256'],old.get('runtime_evidence_sha256'))
    if source_key!=reuse['key'] or key(directory,receipt['orca_executable_sha256'],receipt.get('runtime_evidence_sha256'))!=source_key:
        raise ValueError('Reuse inputs/method/runtime differ')
    if receipt.get('output_sha256')!=old.get('output_sha256'):raise ValueError('Reuse output inventory differs')


def try_reuse(root,directory,executable_hash):
    if not config(root).get('execution',{}).get('verified_exact_input_reuse',False):return None
    from .gates import inspect_job
    from .runner import _safe_file,_copy_new
    from .orcaio import _input_hashes
    directory=Path(directory);digest=key(directory,executable_hash,runtime_identity())
    index=package_root(root)/'operations/cache'/(digest+'.json')
    if not index.is_file():return None
    record=readj(index);source=Path(record['source']).resolve()
    if source==directory:return None
    if not source.is_relative_to(package_root(root)):raise ValueError('Cache source escapes package')
    if sha(source/'receipt.json')!=record['source_receipt_sha256']:raise ValueError('Cached execution changed; audit before recomputation')
    old=readj(source/'receipt.json')
    if old.get('reuse') or not inspect_job(source,config(root))['accepted']:return None
    if key(source,old['orca_executable_sha256'],old.get('runtime_evidence_sha256'))!=digest:return None
    for name,expected in old['output_sha256'].items():
        src=_safe_file(source,name);dst=_safe_file(directory,name,existing=False)
        if name in _input_hashes(directory) or name=='receipt.json':raise ValueError('Unsafe cache output')
        if dst.exists():raise ValueError('Do not overwrite partially copied cache outputs')
        if sha(src)!=expected:raise ValueError('Cached output changed')
    meta=readj(directory/'job.json')
    receipt={**old,'id':meta['id'],'stage':meta['stage'],'input_sha256':_input_hashes(directory),
             'quantum_jobs_launched':0,'receipt_kind':'VERIFIED_REUSE_NOT_NEW_ORCA_EXECUTION',
             'reuse':{'source':str(source),'source_receipt_sha256':record['source_receipt_sha256'],'key':digest}}
    for name,expected in old['output_sha256'].items():
        _copy_new(_safe_file(source,name),_safe_file(directory,name,existing=False),expected)
    validate_reuse(directory,receipt,config(root))
    writej(directory/'receipt.json',receipt)
    return receipt
