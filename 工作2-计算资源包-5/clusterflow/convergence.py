"""Source-bound convergence review; successful jobs alone never imply convergence.

Two successive sampling waves must be compared. Reported H/C shifts are
recalculated from real output-bound ensembles before the numerical stability
test. Motif coverage and basin rediscovery still require documented review.
"""
from pathlib import Path
from datetime import datetime
import math
import numpy as np
from .common import config,hash_protocol,inside,readj,sha

def _bound(root,record):
    p=inside(root,record['path'])
    if sha(p)!=record['sha256']:raise ValueError('Convergence evidence changed: '+record['path'])
    return p

def _recompute_snapshot(root,result,species,scenario,method):
    from .collect import validate_atom_map,_shielding_ensemble,_validate_tms,methyl_average
    from .gates import inspect_job,boltzmann
    from .pipeline import _review_snapshot, _validate_stage_source
    p=config(root); ph=hash_protocol(p); meta=result['metadata']
    for k,v in [('species_id',species),('scenario_id',scenario),('j_method',method),('protocol_hash',ph)]:
        if meta.get(k)!=v:raise ValueError('Convergence snapshot scope mismatch: '+k)
    snap=meta['ensemble_snapshot']; computed={}; fingerprints=[]
    for role,expected in [('target',species),('reference','TMS_DMSO')]:
        e=snap[role]
        frozen=_bound(root,meta.get('frozen_ensembles',{})[role])
        if readj(frozen)!=e:raise ValueError('Embedded ensemble differs from frozen reviewed snapshot')
        if (e['species_id'],e['scenario_id'],e['protocol_hash'])!=(expected,scenario,ph):
            raise ValueError('Snapshot target/reference scope mismatch')
        atoms=validate_atom_map(readj(inside(root,e['atom_map'])))
        evid=meta['source' if role=='target' else 'reference_source']
        if sha(inside(root,e['atom_map']))!=evid['atom_map_sha256']:raise ValueError('Snapshot map changed')
        gs=[]
        for c in e['conformers']:
            check=inspect_job(inside(root,c['optfreq_job']),p)
            if not check['accepted']:raise ValueError('Snapshot Opt/Freq no longer accepted')
            if _review_snapshot(check)!=c.get('scientific_review'):
                raise ValueError('Snapshot scientific review changed; regenerate the downstream evidence')
            if check['job']['conformer_id']!=c['id'] or check['job']['species_id']!=expected or check['job']['scenario_id']!=scenario:
                raise ValueError('Snapshot conformer source mismatch')
            actual_g=check['parsed']['gibbs_free_energy_Eh']
            if not math.isclose(c['G_Eh'],actual_g,rel_tol=0,abs_tol=1e-10):
                raise ValueError('Snapshot stored Gibbs energy differs from accepted output')
            property_dir=inside(root,c['shielding_job'])
            property_meta=readj(property_dir/'job.json')
            if property_meta.get('source',{}).get('scientific_review')!=_review_snapshot(check):
                raise ValueError('Snapshot shielding source scientific review changed')
            # A property's own receipt does not bind its geometry/review to
            # the accepted Opt/Freq source. Check that parent linkage as well.
            _validate_stage_source(root,property_dir,property_meta,e.get('import_index_snapshot',{}))
            gs.append(check['parsed']['gibbs_free_energy_Eh'])
        weights=boltzmann(gs,p['temperature_K'])
        if not np.allclose(weights,[c['weight'] for c in e['conformers']],rtol=1e-10,atol=1e-12):
            raise ValueError('Snapshot weights do not match real final Gibbs energies')
        sigma,_,_,fingerprint=_shielding_ensemble(root,e,atoms,weights,p,evid)
        fingerprints.append(fingerprint)
        if role=='target':sigma,_,_=methyl_average(atoms,sigma,{})
        else:_validate_tms(atoms)
        computed[role]=(sigma,atoms)
    if len(set(fingerprints))!=1:raise ValueError('Snapshot shielding/reference method mismatch')
    ts,ta=computed['target']; rs,ra=computed['reference']
    reference={el:float(np.mean([v for i,v in rs.items() if ra[i]['element']==el])) for el in ['H','C']}
    shifts={i:reference[ta[i]['element']]-v for i,v in ts.items()}
    rows={r['atom_index_zero_based']:r for r in result['atom_rows']}
    if set(rows)!=set(shifts) or any(not math.isclose(rows[i]['shift_ppm'],v,abs_tol=1e-8) for i,v in shifts.items()):
        raise ValueError('Snapshot reported shifts differ from bound outputs')
    return shifts,ta

def _completed_sampling(root,ids,species,scenario):
    from .pipeline import draft_goat_review
    known={j['job_id']:j for j in readj(Path(root)/'registry/jobs.json')}
    counts={species:0,'TMS_DMSO':0}
    if len(ids)!=len(set(ids)):raise ValueError('Repeated GOAT IDs in convergence wave')
    for jid in ids:
        if jid not in known:raise ValueError('Unknown released GOAT job')
        j=known[jid]
        if j['species_id'] not in counts or j['scenario_id']!=scenario:raise ValueError('Unrelated job in convergence wave')
        # Draft creation verifies actual completed output and receipt bindings;
        # generating the draft does not itself approve scientific convergence.
        draft_goat_review(root,jid)
        counts[j['species_id']]+=1
    required=5 if '_W01_' in scenario else 3 if '_02_' in scenario else 2
    advisory=config(root).get('conformer_policy',{}).get('convergence_thresholds_are_advisory') is True
    if min(counts.values())<1:raise ValueError('No successful target/reference sampling evidence')
    if min(counts.values())<required and not advisory:
        raise ValueError('Too few successful target/reference sampling starts')
    return counts

def validate_gate(root,gate_path,species_id,scenario_id,j_method):
    root=Path(root).resolve(); path=Path(gate_path)
    if not path.is_absolute():path=inside(root,str(path))
    g=readj(path)
    if g.get('status')!='REVIEWED_CONVERGENCE':raise ValueError('Convergence review is not approved')
    if not all(isinstance(g.get(k),str) and g[k].strip() for k in ['reviewer','reviewed_utc','reason']):
        raise ValueError('Convergence reviewer/time/reason missing')
    if g['reviewer'].strip().lower() in ['auto','automatic','codex','agent','system','unknown','todo']:
        raise ValueError('Software cannot approve scientific convergence')
    datetime.fromisoformat(g['reviewed_utc'].replace('Z','+00:00'))
    for k in ['basin_rediscovery','no_new_important_basin','motif_coverage','peak_order_stable']:
        if g.get('checks',{}).get(k) is not True:raise ValueError('Required scientific convergence review not passed: '+k)
    evidence=g.get('review_evidence',[])
    if not evidence:raise ValueError('Documented basin/motif review evidence required')
    for e in evidence:_bound(root,e)
    waves=g.get('waves',[])
    if len(waves)<2:raise ValueError('At least two successive sampling-wave snapshots required')
    previous_ids=set(); previous=None; seen_results=set(); maxdelta={'nonexchangeable_1H_ppm':0.0,'exchangeable_H_ppm':0.0,'13C_ppm':0.0}
    for wave in waves:
        identity=(wave['result']['path'],wave['result']['sha256'])
        if identity in seen_results or any(identity[1]==old[1] for old in seen_results):
            raise ValueError('Repeated prediction cannot establish sampling convergence')
        seen_results.add(identity)
        result=readj(_bound(root,wave['result']))
        ids=set(wave['released_goat_jobs'])
        counts=_completed_sampling(root,wave['released_goat_jobs'],species_id,scenario_id)
        if previous_ids and not previous_ids<ids:raise ValueError('Later wave must include additional real sampling starts')
        used_ids=set()
        for role in ['target','reference']:
            e=result['metadata']['ensemble_snapshot'][role]
            imported=e.get('import_index_snapshot')
            if not isinstance(imported,dict):raise ValueError('Missing frozen import/source snapshot')
            source_ids={s['goat_job_id'] for s in imported['sources']}
            selected_ids={s['goat_job_id'] for s in imported['conformers']}
            if source_ids!=selected_ids:raise ValueError('A released start lacks selected-frame source evidence')
            used_ids.update(source_ids)
        if used_ids!=ids:
            raise ValueError('Wave GOAT IDs not exactly bound to this result import selection')
        current,atoms=_recompute_snapshot(root,result,species_id,scenario_id,j_method)
        if previous is not None:
            if set(previous)!=set(current):raise ValueError('Atom identities changed between sampling waves')
            maxdelta={k:0.0 for k in maxdelta}
            for i,v in current.items():
                key='13C_ppm' if atoms[i]['element']=='C' else 'exchangeable_H_ppm' if atoms[i].get('exchangeable_H') else 'nonexchangeable_1H_ppm'
                maxdelta[key]=max(maxdelta[key],abs(v-previous[i]))
        previous,current_meta,previous_ids=current,result['metadata'],ids
    # Last snapshot must represent the current, complete target/reference result.
    for k in ['source','reference_source']:
        ev=current_meta[k]
        if sha(inside(root,ev['ensemble_path']))!=ev['ensemble_sha256']:
            raise ValueError('Convergence gate is stale: current ensemble differs from final snapshot')
        role='target' if k=='source' else 'reference'
        if readj(inside(root,ev['ensemble_path']))!=current_meta['ensemble_snapshot'][role]:
            raise ValueError('Final convergence snapshot is not the complete current ensemble')
    limits={'nonexchangeable_1H_ppm':0.05,'exchangeable_H_ppm':0.2,'13C_ppm':0.5}
    advisory=config(root).get('conformer_policy',{}).get('convergence_thresholds_are_advisory') is True
    exceeded={k:{'observed':maxdelta[k],'reference':limits[k]} for k in limits if maxdelta[k]>limits[k]}
    if exceeded and not advisory:raise ValueError('Last-wave calculated-shift stability gate failed')
    return {'status':'SOURCE_BOUND_REVIEWED_CONVERGENCE','gate_sha256':sha(path),
            'maximum_last_wave_shift_change':maxdelta,'engineering_limits_not_accuracy_guarantees':limits,
            'numerical_references_are_advisory':advisory,'advisory_exceedances':exceeded,
            'notice':'Numerical stability plus documented review, not a proof of exhaustive global sampling.'}
