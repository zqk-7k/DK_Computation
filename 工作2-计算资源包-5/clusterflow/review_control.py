"""Explicit, source-bound recovery and limited-scope shortfall decisions."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import uuid

from .common import readj, writej, sha, hash_protocol, config, inside


def _digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()


def draft_shortfall(workspace,sid,scenario):
    from .pipeline import final_stable_candidates
    from .cache import package_root
    from .autoflow import _plan, _job_status
    workspace=Path(workspace).resolve();root=package_root(workspace)
    plan,index=_plan(root,sid)
    audit=final_stable_candidates(workspace,sid,scenario)
    sampling={jid:_job_status(root,index[jid],goat=True) for jid in plan['goat_start_ids']}
    return {'status':'DRAFT_SHORTFALL_SCOPE_NOT_APPROVED','species_id':sid,'scenario_id':scenario,
            'workspace':str(workspace.relative_to(root)), 'protocol_hash':hash_protocol(config(workspace)),
            'original_target':plan['stable_conformer_target'],'actual_stable_count':audit['eligible_count'],
            'stable_audit_sha256':_digest(audit),'sampling_sha256':_digest(sampling),
            'all_pool_goat_verified':all(v['status']=='COMPLETE' for v in sampling.values()),
            'all_final_jobs_pass':audit['all_prepared_optfreq_pass'] and not audit['issues'],
            'reviewer':'','reason':'','limitations':'','evidence':[],
            'checks':{'rediscovery_and_motifs_reviewed':False,'available_property_stability_reviewed':False,
                      'actual_count_and_limited_scope_accepted':False},
            'notice':'Original target remains unchanged. No final sampling/experimental accuracy approval.'}


def _validate_shortfall(workspace,sid,scenario,review):
    from .cache import package_root
    current=draft_shortfall(workspace,sid,scenario)
    for field in ('species_id','scenario_id','workspace','protocol_hash','original_target','actual_stable_count',
                  'stable_audit_sha256','sampling_sha256'):
        if review.get(field)!=current[field]:raise ValueError('Shortfall review source changed: '+field)
    if review.get('status')!='APPROVED_LIMITED_SHORTFALL_SCOPE':raise ValueError('Shortfall review is not approved')
    if not current['all_pool_goat_verified'] or not current['all_final_jobs_pass']:
        raise ValueError('Complete verified sampling pool and clean final minima required')
    if not 0<current['actual_stable_count']<current['original_target']:
        raise ValueError('Shortfall decision requires a positive actual count below the unchanged target')
    for key in ('reviewer','reason','limitations'):
        if not isinstance(review.get(key),str) or not review[key].strip():raise ValueError('Review requires '+key)
    if review['reviewer'].lower().strip() in ('auto','automatic','codex','agent','system','ai','unknown'):
        raise ValueError('Do not fabricate a human shortfall decision')
    if any(review.get('checks',{}).get(k) is not True for k in current['checks']):raise ValueError('Shortfall checks incomplete')
    evidence=review.get('evidence',[])
    if not evidence:raise ValueError('Source-bound sampling/property review evidence required')
    for row in evidence:
        if sha(inside(package_root(workspace),row['path']))!=row['sha256']:raise ValueError('Shortfall evidence changed')
    return current


def install_shortfall(workspace,sid,scenario,review_path,execute=False):
    from .cache import package_root
    workspace=Path(workspace).resolve();review=readj(review_path)
    current=_validate_shortfall(workspace,sid,scenario,review)
    target=package_root(workspace)/'operations/target_exceptions'/(_digest(current['workspace']+'|'+sid)+'.json')
    if execute:
        if target.exists():raise ValueError('Shortfall decision already exists; preserve its history')
        writej(target,review)
    return {'status':'SHORTFALL_SCOPE_INSTALLED' if execute else 'SHORTFALL_SCOPE_PLAN',
            'path':str(target),'original_target':current['original_target'],'effective_count':current['actual_stable_count']}


def target_for(workspace,sid,scenario,original):
    from .cache import package_root
    workspace=Path(workspace).resolve();root=package_root(workspace)
    relative=str(workspace.relative_to(root))
    path=root/'operations/target_exceptions'/(_digest(relative+'|'+sid)+'.json')
    if not path.is_file():return original,None
    review=readj(path);current=_validate_shortfall(workspace,sid,scenario,review)
    return current['actual_stable_count'],{'path':str(path),'sha256':sha(path),'original_target':original,
           'effective_count':current['actual_stable_count'],'limitations':review['limitations']}


def resume_branch(root,sid,reason,execute=False,expand=False):
    from . import autoflow
    from .operations import retire_terminal_lock
    from scripts.packed_worker import queue_lock
    root=Path(root).resolve()
    target,_=autoflow._plan(root,sid)
    if not reason.strip():raise ValueError('Recovery reason required')
    base=root/'automation'/sid;path=base/'state.json'
    report={'status':'RESUME_PLAN','species_id':sid,'reason':reason,'scientific_convergence_approved':False}
    if not execute:
        report['state']=readj(path) if path.exists() else None
        return report
    with queue_lock(root/'operations/packed'):
        claim=root/'operations/packed/claims'/(sid+'.json')
        if claim.exists() and not retire_terminal_lock(claim):raise ValueError('Branch is claimed; reconcile its owner before resuming')
        lock=base/'.autoflow.lock'
        if lock.exists() and not retire_terminal_lock(lock):raise ValueError('Branch owner still active/unknown')
        with autoflow._lock(base):
            state=readj(path)
            for name,digest in state['identity_sha256'].items():
                if sha(root/name)!=digest:raise ValueError('Frozen protocol/registry changed; do not resume this run')
            for wave in state['waves']:autoflow._verify_workspace(root,wave)
            old=readj(path)
            if expand:
                if not state['waves']:raise ValueError('No existing wave to expand')
                plans={sid:target}
                if sid!='TMS_DMSO':plans['TMS_DMSO']=autoflow._plan(root,'TMS_DMSO')[0]
                if not autoflow._expand_candidates(state,plans,plans):raise ValueError('No eligible reserve/start remains')
            else:
                if state['waves'] and 'prediction' in state['waves'][-1]:
                    raise ValueError('Completed spectrum wave requires --expand or final review, not replay')
                state['phase']='DFT' if state['waves'] else 'SAMPLING'
            decision=base/'resume_history'/(uuid.uuid4().hex+'.json')
            writej(decision,{**report,'previous_state':old,'new_phase':state['phase'],'time':datetime.now(timezone.utc).isoformat()})
            state['review_reason']=None;state['last_status']='EXPLICIT_RECHECK_REQUESTED'
            writej(path,state)
            blocked=root/'operations/packed/blocked'/(sid+'.json')
            if blocked.exists():
                archive=root/'operations/packed/block_history';archive.mkdir(exist_ok=True)
                blocked.rename(archive/(sid+'_'+uuid.uuid4().hex+'.json'))
    return {**report,'status':'RESUMED_FOR_REVALIDATION_NOT_SCIENTIFIC_APPROVAL','phase':state['phase']}
