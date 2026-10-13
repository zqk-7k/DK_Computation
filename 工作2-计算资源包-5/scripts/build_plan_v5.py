"""Mechanically derive the v5 stable-minimum targets; does not run chemistry."""
from pathlib import Path
import json
from collections import Counter

ROOT = Path(__file__).resolve().parents[1]

def read(path):
    return json.loads((ROOT / path).read_text(encoding='utf-8-sig'))

def save(path, value):
    (ROOT / path).write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')

def main():
    if (ROOT / 'automation').exists():
        raise RuntimeError('Do not regenerate a plan after a run has started')
    old = read('registry/sampling_plan_v4.json')
    jobs = read('registry/jobs.json')
    lookup = {8: ('R', 3), 12: ('M', 5), 16: ('H', 10), 24: ('X', 15), 32: ('X', 15)}
    rows = []
    for previous in old['species']:
        row = {k: v for k, v in previous.items() if not k.startswith('suggested_unique') and k not in
               ('advisory_no_quota', 'unique_count_meaning', 'n_starts_per_species')}
        sid = row['species_id']
        category, quota = lookup[row['n_starts']]
        pool = [j for j in jobs if j['species_id'] == sid]
        if len(pool) != row['n_starts']:
            raise ValueError('GOAT pool does not match existing model: '+sid)
        row.update(complexity_class=category, stable_conformer_target=quota,
                   stable_target_meaning='MINIMUM_DISTINCT_FINAL_DFT_MINIMA_BEFORE_FIRST_NMR_NOT_GOAT_STARTS_NOT_ACCURACY_PROOF',
                   initial_candidate_batch_size=2*quota,
                   initial_goat_starts=min(7,len(pool)), addition_goat_starts=4,
                   goat_start_ids=[j['job_id'] for j in pool],
                   temperature_status='USER_CONFIRMED_COMPUTATIONAL_298.15_K_NOT_EXPERIMENTAL_MEASUREMENT',
                   actual_stable_conformer_count=None,
                   count_decision='USER_TARGET_BY_COMPLEXITY_FINAL_DFT_STABLE_MINIMA')
        for index, job in enumerate(pool):
            job['initial_release'] = index < row['initial_goat_starts']
        rows.append(row)
    products = []
    for previous in old['products']:
        product = {k:v for k,v in previous.items() if not k.startswith('suggested_unique')}
        branches = [r for r in rows if r['compound_id'] == product['compound_id']]
        product.update(complexity_class=branches[0]['complexity_class'],
                       stable_target_per_species=branches[0]['stable_conformer_target'],
                       stable_target_all_branches=sum(r['stable_conformer_target'] for r in branches),
                       initial_candidate_batch_size_per_species=branches[0]['initial_candidate_batch_size'],
                       initial_goat_starts_per_species=branches[0]['initial_goat_starts'],
                       temperature_K=298.15,
                       actual_stable_conformer_count=None)
        products.append(product)
    summary = dict(old['summary'])
    summary.update(initially_released_goat_starts=sum(j['initial_release'] for j in jobs),
                   held_goat_starts=sum(not j['initial_release'] for j in jobs),
                   stable_target_product_branches=sum(r['stable_conformer_target'] for r in rows if r['compound_id']!='TMS'),
                   stable_target_TMS=next(r['stable_conformer_target'] for r in rows if r['compound_id']=='TMS'),
                   product_complexity_counts=dict(Counter(p['complexity_class'] for p in products)),
                   targets_are_not_guaranteed_outputs=True)
    plan = dict(schema_version='5.0', package_name=ROOT.name,
                scenario_id=old['scenario_id'], summary=summary, species=rows, products=products,
                interpretation='Targets count distinct eligible FINAL DFT minima after optimization, frequencies and deduplication, not starting geometries or GOAT candidates. Different structural branches remain separate.',
                solvent_note='Uniform 2 DMSO-d6 + 1 ordinary H2O cluster with DMSO continuum is a model assumption, not a proven universal solvation shell or measured water content.',
                status='PLANNED_NOT_QUANTUM_CALCULATED')
    save('registry/jobs.json', jobs)
    save('registry/batch_plan_v5.json', plan)
    protocol = read('config/protocol.json')
    protocol.update(schema_version='5.0', package_name=ROOT.name,
                    temperature_status='USER_CONFIRMED_COMPUTATIONAL_298.15_K_NOT_EXPERIMENTAL_MEASUREMENT')
    protocol['conformer_policy'].update(convergence_thresholds_are_advisory=True,
        stable_minimum_target_is_not_convergence=True,
        final_stable_target_plan='registry/batch_plan_v5.json')
    save('config/protocol.json',protocol)
    save('config/automation.json',{
        'schema_version':'5.0',
        'automatic_handoff':{'enabled':True,'explicit_execute_required':True,
            'complete_entire_GOAT_batch_before_selection':True,
            'stable_target_scope':'FINAL_DFT_UNIQUE_MINIMA_BEFORE_NMR',
            'minimum_targets_by_complexity':{'R':3,'M':5,'H':10,'X':15}},
        'selection':{'rmsd_A':0.10,'representative_energy_window_kcal_mol':6.0,
                     'retain_high_energy_reserve_but_do_not_force_quota':True},
        'stable_minima':{'energy_window_kcal_mol':6.0,'negative_frequencies_eligible':False,
                         'human_reviewed_approximations_count_toward_target':False},
        'stability_advice':{'nonexchangeable_H_ppm':0.05,'C_ppm':0.5,
                            'thresholds_are_advisory_not_accuracy':True},
        'scientific_convergence_may_be_auto_approved':False,
        'finite_pool_exhaustion':'REPORT_SHORTFALL_NEVER_PAD_OR_SILENTLY_LOWER_TARGET',
        'failed_DFT':'PAUSE_FOR_DIAGNOSIS_DO_NOT_SILENTLY_EXCLUDE_LOW_ENERGY_CANDIDATE',
        'notice':'Software planning; no actual GOAT, DFT, or spectra were calculated during package creation.'})
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
