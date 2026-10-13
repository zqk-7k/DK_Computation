"""Generate planned main-scenario starting clusters; never run GOAT or DFT."""
from __future__ import annotations
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shutil
from rdkit import Chem
import micro64 as m

ROOT=Path(__file__).resolve().parents[1]

def build_one(s,plan,old):
    sid=s['species_id'];scenario=plan['scenario_id'];source=ROOT/'solute_models'/sid
    result=[]
    for i in range(1,plan['n_starts']+1):
        jobid=f'{sid}__{scenario}__s{i:02d}'
        if jobid in old:
            row=old[jobid]
            d=ROOT/row['job_directory']
            for name,key in [('goat.inp','input_sha256'),('start.xyz','xyz_sha256'),('atom_map.json','atom_map_sha256'),('cluster.mol','cluster_mol_sha256'),('provenance.json','provenance_sha256')]:
                if m.digest(d/name)!=row[key]:raise ValueError('Existing starting input changed: '+jobid)
            result.append(dict(row));continue
        analyte,prep=m.new_analyte(source,i,m.seed_for(f'{sid}|{i}|v4_51ref'))
        seed=m.seed_for(jobid+'|v4')
        cluster,components,motifs=m.add_environment(analyte,2,1,seed,i)
        amap=m.mapping(cluster,components)
        d=ROOT/'01_goat_explicit'/sid/scenario/f'start_{i:02d}'
        m.put(d/'goat.inp',m.goat_input(s['charge'],s['multiplicity']))
        m.put(d/'start.xyz',m.xyz(cluster,jobid+'; INITIAL MODEL; NOT QM OPTIMIZED'))
        m.put(d/'cluster.mol',Chem.MolToMolBlock(cluster))
        m.putj(d/'atom_map.json',amap);m.putj(d/'isotope_map.json',[r for r in amap if r['isotope']])
        m.putj(d/'provenance.json',{'job_id':jobid,'source_model_sha256':m.digest(source/'model.mol'),
            'source_species':s,'analyte_start':prep,'environment_seed':seed,'motifs':motifs,
            'environment_atoms_excluded_from_analyte_nmr':True,'cluster_optimization':'NOT_RUN','orca_calculation':'NOT_RUN',
            'stereo_geometry_issues':m.geometry_issues(Chem.MolFromMolBlock((source/'model.mol').read_text(encoding='utf-8'),removeHs=False),analyte),
            'composition':{'DMSO_d6':2,'H2O':1},'reference_commit':m.SHA,
            'assumption':'2DMSO-d6+1H2O local model with DMSO continuum; no measured water concentration or converged shell implied.'})
        result.append({'job_id':jobid,'species_id':sid,'compound_id':s['compound_id'],'scenario_id':scenario,
            'start_index':i,'initial_release':True,'status':'PREPARED_NOT_RUN','job_directory':d.relative_to(ROOT).as_posix(),
            'input_sha256':m.digest(d/'goat.inp'),'xyz_sha256':m.digest(d/'start.xyz'),
            'atom_map_sha256':m.digest(d/'atom_map.json'),'cluster_mol_sha256':m.digest(d/'cluster.mol'),
            'provenance_sha256':m.digest(d/'provenance.json'),'atom_count':cluster.GetNumAtoms(),
            'solute_atom_count':analyte.GetNumAtoms(),'charge':s['charge'],'multiplicity':s['multiplicity'],
            'DMSO_d6_count':2,'H2O_count':1,'orca_runtime_validated':False})
    print(f"Prepared {sid}: {len(result)} starts; no QM",flush=True)
    return result

def main(args):
    plan=m.readj(ROOT/'registry/sampling_plan_v4.json');species=m.readj(ROOT/'registry/species.json')
    byid={r['species_id']:r for r in plan['species']}
    if set(byid)!={s['species_id'] for s in species}:raise ValueError('Sampling plan species mismatch')
    hist=ROOT/'history/v3/registry';hist.mkdir(parents=True,exist_ok=True)
    for name in ['jobs.json','scenarios.json','jobs_subset.json','scenarios_subset.json']:
        src=ROOT/'registry'/name;dst=hist/name
        if src.exists() and not dst.exists():shutil.copy2(src,dst)
    old={r['job_id']:r for r in m.readj(hist/'jobs.json')}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        groups=list(pool.map(lambda s:build_one(s,byid[s['species_id']],old),species))
    jobs=sorted([j for group in groups for j in group],key=lambda j:j['job_id'])
    counts=Counter(j['species_id'] for j in jobs)
    if any(counts[k]!=v['n_starts'] for k,v in byid.items()):raise ValueError('Prepared counts mismatch')
    def write(path,value):
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    scenarios=[{**s,'scenario_id':plan['scenario_id'],'explicit_DMSO_d6':2,'explicit_H2O':1,
        'initial_starts':byid[s['species_id']]['n_starts'],'suggested_unique_min':byid[s['species_id']]['suggested_unique_min'],
        'suggested_unique_max':byid[s['species_id']]['suggested_unique_max'],'advisory_no_quota':True,
        'temperature_K':298.15,'temperature_provisional':True,'pooling_key':s['species_id']+'|'+plan['scenario_id'],
        'reference_species_id':'TMS_DMSO'} for s in species]
    write(ROOT/'registry/jobs.json',jobs);write(ROOT/'registry/scenarios.json',scenarios)
    report={'package_name':ROOT.name,'status':'INITIAL_STRUCTURES_PREPARED_NO_QM',**plan['summary'],
        'jobs':len(jobs),'active_composition_scenarios':len(scenarios),'sampling_converged':False,
        'ORCA_GOAT_run':False,'DFT_run':False,'NMR_spectra_generated':False,
        'reused_initial_jobs':sum(j['job_id'] in old for j in jobs),'new_initial_jobs':sum(j['job_id'] not in old for j in jobs),
        'inactive_old_control_directories':'Retained for provenance; absent from active registry/jobs.json',
        'jobs_manifest_sha256':m.digest(ROOT/'registry/jobs.json'),'sampling_plan_sha256':m.digest(ROOT/'registry/sampling_plan_v4.json')}
    write(ROOT/'reports/preparation.json',report)
    with (ROOT/'structures_main_v4.sdf').open('w',encoding='utf-8',newline='\n') as out:
        for j in jobs:
            block=(ROOT/j['job_directory']/'cluster.mol').read_text(encoding='utf-8')
            out.write(block.rstrip()+'\n>  <JOB_ID>\n'+j['job_id']+'\n\n>  <STATUS>\nINITIAL_NOT_QM\n\n$$$$\n')
    print(json.dumps(report,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workers',type=int,default=4)
    main(p.parse_args())
