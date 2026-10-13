"""Prepare, never execute, 18O/DMSO-d6 GOAT starting clusters.

The 51-package composition/sampling contract is retained. Its isotope and
distance-inferred bond maps are NOT reused: authoritative 64-product MOL
topologies and atom ordering are retained throughout.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Lipinski

ROOT = Path(__file__).resolve().parents[1]
SHA = '933cb6159294813a79f443921920127032fb01d4'
M18 = 17.99915961286
MD = 2.01410177812
SCENARIOS = (
    ('DMSO_00_W00_continuum', 0, 0, 3, 2),
    ('DMSO_02_W00_primary', 2, 0, 4, 3),
    ('DMSO_02_W01_trace_water_sensitivity', 2, 1, 7, 5),
)
CEILINGS = {'R': (4, 6, 10), 'M': (6, 10, 16), 'H': (8, 12, 20), 'X': (10, 16, 24)}

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def readj(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def put(path, content):
    """Idempotent preparation; never overwrite a different existing artifact."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = content.encode('utf-8') if isinstance(content, str) else content
    if path.exists() and path.read_bytes() != data:
        raise RuntimeError(f'Refusing to overwrite changed artifact: {path}')
    if not path.exists():
        path.write_bytes(data)

def putj(path, obj):
    put(path, json.dumps(obj, ensure_ascii=False, indent=2) + '\n')

def seed_for(text):
    return int(hashlib.sha256(text.encode()).hexdigest()[:8], 16) & 0x7fffffff

def canon(m):
    return Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(m)), isomericSmiles=True)

def geometry_issues(template, m):
    issues = []
    if [(a.GetSymbol(), a.GetIsotope()) for a in m.GetAtoms()] != [(a.GetSymbol(), a.GetIsotope()) for a in template.GetAtoms()]:
        return ['atom_order_or_isotope_mismatch']
    q = Chem.Mol(m)
    for b in q.GetBonds():
        b.SetBondDir(Chem.BondDir.NONE)
    Chem.AssignStereochemistryFrom3D(q, confId=0, replaceExistingTags=True)
    if canon(q) != canon(template):
        issues.append('coordinate_derived_stereochemistry_mismatch')
    xyz = m.GetConformer().GetPositions()
    pt = Chem.GetPeriodicTable()
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        r = sum(pt.GetRcovalent(m.GetAtomWithIdx(k).GetAtomicNum()) for k in (i, j))
        d = np.linalg.norm(xyz[i] - xyz[j])
        if not 0.57*r < d < 1.65*r:
            issues.append(f'bond_distance_{i}_{j}')
    for i in range(m.GetNumAtoms()):
        for j in range(i):
            if m.GetBondBetweenAtoms(i, j):
                continue
            r = sum(pt.GetRcovalent(m.GetAtomWithIdx(k).GetAtomicNum()) for k in (i, j))
            if np.linalg.norm(xyz[i] - xyz[j]) < 0.55*r:
                issues.append(f'nonbonded_clash_{i}_{j}')
    return issues

def set_xyz(m, xyz):
    q = Chem.Mol(m)
    q.RemoveAllConformers()
    c = Chem.Conformer(q.GetNumAtoms())
    c.Set3D(True)
    for i, p in enumerate(xyz):
        c.SetAtomPosition(i, tuple(float(x) for x in p))
    q.AddConformer(c)
    return q

def new_analyte(source_dir, start_index, seed):
    template = Chem.MolFromMolBlock((source_dir/'model.mol').read_text(encoding='utf-8'), removeHs=False)
    if template is None:
        raise ValueError('Unreadable model MOL')
    m = Chem.Mol(template)
    m.RemoveAllConformers()
    failure = []
    for trial in range(2):
        p = AllChem.ETKDGv3()
        p.randomSeed = (seed + trial*701) & 0x7fffffff
        p.numThreads = 1
        p.enforceChirality = True
        p.useSmallRingTorsions = True
        p.useMacrocycleTorsions = True
        p.useRandomCoords = bool(trial)
        m.RemoveAllConformers()
        if AllChem.EmbedMolecule(m, p) != 0:
            failure.append('ETKDG_embedding_failed')
            continue
        ffname, converged = 'NONE', False
        try:
            if AllChem.MMFFHasAllMoleculeParams(m):
                ffname = 'MMFF94s'
                converged = AllChem.MMFFOptimizeMolecule(m, mmffVariant='MMFF94s', maxIters=400) == 0
            elif AllChem.UFFHasAllMoleculeParams(m):
                ffname = 'UFF'
                converged = AllChem.UFFOptimizeMolecule(m, maxIters=400) == 0
        except (ValueError, RuntimeError) as exc:
            failure.append(type(exc).__name__)
        issues = geometry_issues(template, m)
        if not issues:
            return m, {'method': 'INDEPENDENT_ETKDGv3_WITH_FIXED_ATOM_ORDER', 'seed': seed,
                       'forcefield': ffname, 'forcefield_converged': converged,
                       'rejected_attempts': failure, 'global_sampling_converged': False}
        failure.extend(issues)
    # Difficult cages: retain audited geometry, not a chemically wrong embedding.
    files = sorted(source_dir.glob('c[0-9][0-9][0-9].mol')) or [source_dir/'model.mol']
    chosen = files[(start_index-1) % len(files)]
    base = Chem.MolFromMolBlock(chosen.read_text(encoding='utf-8'), removeHs=False)
    rng = np.random.default_rng(seed)
    for _ in range(20):
        m = set_xyz(base, base.GetConformer().GetPositions()+rng.normal(0, 0.012, (base.GetNumAtoms(), 3)))
        if not geometry_issues(template, m):
            return m, {'method': 'AUDITED_SEED_SMALL_PERTURBATION_NOT_INDEPENDENT_CONFORMER_DISCOVERY',
                       'seed': seed, 'source': chosen.name, 'source_sha256': digest(chosen),
                       'perturbation_sigma_A': 0.012, 'rejected_attempts': failure,
                       'forcefield_converged': False, 'global_sampling_converged': False}
    raise RuntimeError('No stereochemically valid starting geometry')

def solvent(smiles, deuterated=False):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3()
    p.randomSeed = 12478
    if AllChem.EmbedMolecule(m, p) != 0:
        raise RuntimeError('Solvent embedding failed')
    if AllChem.MMFFHasAllMoleculeParams(m):
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    else:
        AllChem.UFFOptimizeMolecule(m, maxIters=500)
    if deuterated:
        for a in m.GetAtoms():
            if a.GetSymbol() == 'H':
                a.SetIsotope(2)
    return m

def rotation(rng):
    u, s, v = np.linalg.svd(rng.normal(size=(3, 3)))
    r = u @ v
    if np.linalg.det(r) < 0:
        u[:, -1] *= -1
        r = u @ v
    return r

def unit(rng):
    v = rng.normal(size=3)
    return v / np.linalg.norm(v)

def acceptable(existing, fragment):
    x, y = existing.GetConformer().GetPositions(), fragment.GetConformer().GetPositions()
    pt = Chem.GetPeriodicTable()
    z1 = [a.GetAtomicNum() for a in existing.GetAtoms()]
    z2 = [a.GetAtomicNum() for a in fragment.GetAtoms()]
    radii1 = np.array([pt.GetRvdw(z) for z in z1])
    radii2 = np.array([pt.GetRvdw(z) for z in z2])
    limits = 0.64*(radii1[:, None] + radii2[None, :])
    base = np.where((np.array(z1)[:, None] == 1) | (np.array(z2)[None, :] == 1), 1.15, 1.65)
    limits = np.maximum(limits, base)
    return bool(np.all(np.linalg.norm(x[:, None, :]-y[None, :, :], axis=2) >= limits))

def donor_targets(analyte):
    x = analyte.GetConformer().GetPositions()
    targets = []
    for a in analyte.GetAtoms():
        if a.GetSymbol() != 'H' or a.GetDegree() != 1:
            continue
        n = a.GetNeighbors()[0]
        if n.GetSymbol() not in ('O', 'N', 'S'):
            continue
        v = x[a.GetIdx()]-x[n.GetIdx()]
        targets.append((a.GetIdx(), x[a.GetIdx()] + 1.90*v/np.linalg.norm(v)))
    return targets

def place_one(existing, analyte, template, kind, rng, index, start_index):
    tx = template.GetConformer().GetPositions()
    oi = next(a.GetIdx() for a in template.GetAtoms() if a.GetSymbol() == 'O')
    donors = donor_targets(analyte)
    rng.shuffle(donors)
    aims = []
    if kind == 'DMSO_d6':
        aims = [(oi, point, f'analyte_H{idx+1}_to_DMSO_O_start') for idx, point in donors]
    else:
        # Water donates to N/O acceptors, or accepts from an analyte X-H donor.
        acc = sorted({i for match in Lipinski._HAcceptors(analyte) for i in match
                      if analyte.GetAtomWithIdx(i).GetSymbol() in ('O', 'N', 'S')})
        rng.shuffle(acc)
        hi = next(a.GetIdx() for a in template.GetAtoms() if a.GetSymbol() == 'H')
        x = analyte.GetConformer().GetPositions()
        for idx in acc:
            aims.append((hi, (idx, x[idx]), f'water_H_to_analyte_acceptor_{idx+1}_start'))
        if start_index % 3 == 0:
            aims = [(oi, point, f'analyte_H{idx+1}_to_water_O_start') for idx, point in donors] + aims
    for anchor, point, motif in aims[:4]:
        for _ in range(160):
            r = rotation(rng)
            if isinstance(point, tuple):
                # Put water H toward the acceptor, with water O pointing away.
                _, acceptor = point
                h_to_o = (tx[oi]-tx[anchor]) @ r
                target = acceptor + 1.95*h_to_o/np.linalg.norm(h_to_o)
            else:
                target = point
            coords = (tx-tx[anchor]) @ r + target
            candidate = set_xyz(template, coords)
            if acceptable(existing, candidate):
                return candidate, motif
    # Random first-shell surface starts are hypotheses, not assigned H bonds.
    ax = analyte.GetConformer().GetPositions()
    heavies = [a.GetIdx() for a in analyte.GetAtoms() if a.GetAtomicNum() > 1]
    centre = ax[heavies].mean(axis=0)
    for trial in range(2400):
        idx = int(rng.choice(heavies))
        v = ax[idx]-centre + 1.5*unit(rng)
        v /= np.linalg.norm(v)
        target = ax[idx] + rng.uniform(2.7, 3.8)*v
        r = rotation(rng)
        coords = (tx-tx.mean(axis=0)) @ r + target
        candidate = set_xyz(template, coords)
        if acceptable(existing, candidate):
            return candidate, 'surface_orientation_candidate_no_specific_Hbond_assignment'
    raise RuntimeError(f'Could not place {kind} without interfragment clashes')

def add_environment(analyte, ndmso, nwater, seed, start_index):
    rng = np.random.default_rng(seed)
    cluster = Chem.Mol(analyte)
    components = ['analyte'] * analyte.GetNumAtoms()
    motifs = []
    for k, (kind, template) in enumerate([('DMSO_d6', solvent('CS(C)=O', True))]*ndmso + [('H2O', solvent('O'))]*nwater):
        fragment, motif = place_one(cluster, analyte, template, kind, rng, k, start_index)
        label = f'{kind}_{k+1:02d}'
        cluster = Chem.CombineMols(cluster, fragment)
        components.extend([label]*fragment.GetNumAtoms())
        motifs.append({'component': label, 'initial_placement': motif})
    return cluster, components, motifs

def mapping(m, components):
    rows = []
    for a in m.GetAtoms():
        i = a.GetIdx()
        environment = components[i] != 'analyte'
        iso = a.GetIsotope()
        mass = M18 if iso == 18 and a.GetSymbol() == 'O' else MD if iso == 2 and a.GetSymbol() == 'H' else None
        if iso and mass is None:
            raise ValueError('Unhandled isotope')
        rows.append({'index_zero_based': i, 'index_one_based': i+1, 'element': a.GetSymbol(),
                     'isotope': iso, 'isotope_mass_u': mass, 'component': components[i],
                     'environment_atom': environment, 'source_atom_index_zero_based': None if environment else i,
                     'formal_charge': a.GetFormalCharge(), 'neighbours_zero_based': [n.GetIdx() for n in a.GetNeighbors()],
                     'include_in_analyte_nmr': not environment and a.GetSymbol() in ('H', 'C'),
                     'exchangeable_H': a.GetSymbol() == 'H' and any(n.GetSymbol() in ('N','O','S') for n in a.GetNeighbors())})
    return rows

def xyz(m, label):
    lines = [str(m.GetNumAtoms()), label]
    for a, p in zip(m.GetAtoms(), m.GetConformer().GetPositions()):
        lines.append(f'{a.GetSymbol():2s} {p[0]: .10f} {p[1]: .10f} {p[2]: .10f}')
    return '\n'.join(lines)+'\n'

def goat_input(charge, multiplicity):
    return f'''! GOAT-ENTROPY XTB ALPB(DMSO)
%pal nprocs 8 end
%maxcore 2000
%goat
 NWorkers 8
 RandomSeed false
 AutoWall true
 Align true
 FreezeAmides false
 FreezeCisTrans true
 MaxEn 12.0
 KeepWorkerData false
 ConfTemp 298.15
end
* xyzfile {charge} {multiplicity} start.xyz
'''

def tier_of(s):
    r = s.get('rotatable_bonds', 0)
    return 'X' if r >= 10 else 'H' if r >= 6 or s['compound_id'].startswith('F4_') else 'M' if r >= 3 else 'R'

def build_species(s, source, full):
    sid = s['species_id']
    sd = source/'structures'/sid
    tier = tier_of(s)
    maxstart = max(CEILINGS[tier]) if full else 7
    analytes = [new_analyte(sd, i, seed_for(f'{sid}|{i}|v2_51ref')) for i in range(1, maxstart+1)]
    for name in ('model.mol', 'model_2D.mol', 'atom_map.json', 'species.json', 'structure.svg'):
        if (sd/name).exists():
            put(ROOT/'solute_models'/sid/name, (sd/name).read_bytes())
    jobs, scenarios = [], []
    for pos, (scenario, nd, nw, release, minimum) in enumerate(SCENARIOS):
        count = CEILINGS[tier][pos] if full else release
        scenarios.append({**s, 'scenario_id': scenario, 'explicit_DMSO_d6': nd, 'explicit_H2O': nw,
                          'initial_starts': release, 'minimum_success_not_convergence': minimum,
                          'maximum_candidates': CEILINGS[tier][pos], 'complexity_tier': tier,
                          'tier_basis': 'NEW_TARGET_HEURISTIC_NOT_OLD_COMPOUND_CLASSIFICATION',
                          'temperature_K': 298.15, 'temperature_provisional': True,
                          'pooling_key': f'{sid}|{scenario}', 'reference_species_id': 'TMS_DMSO'})
        for i in range(1, count+1):
            jobid = f'{sid}__{scenario}__s{i:02d}'
            seed = seed_for(jobid)
            analyte, prep = analytes[i-1]
            cluster, components, motifs = add_environment(analyte, nd, nw, seed, i)
            amap = mapping(cluster, components)
            folder = ROOT/'01_goat_explicit'/sid/scenario/f'start_{i:02d}'
            inp = goat_input(s['charge'], s['multiplicity'])
            put(folder/'goat.inp', inp)
            put(folder/'start.xyz', xyz(cluster, jobid+'; INITIAL MODEL; NOT QM OPTIMIZED'))
            put(folder/'cluster.mol', Chem.MolToMolBlock(cluster))
            putj(folder/'atom_map.json', amap)
            putj(folder/'isotope_map.json', [r for r in amap if r['isotope']])
            putj(folder/'provenance.json', {'job_id': jobid, 'source_model_sha256': digest(sd/'model.mol'),
                  'source_species': s, 'analyte_start': prep, 'environment_seed': seed, 'motifs': motifs,
                  'environment_atoms_excluded_from_analyte_nmr': True,
                  'cluster_optimization': 'NOT_RUN', 'orca_calculation': 'NOT_RUN',
                  'stereo_geometry_issues': geometry_issues(Chem.MolFromMolBlock((sd/'model.mol').read_text(encoding='utf-8'), removeHs=False), analyte),
                  'composition': {'DMSO_d6': nd, 'H2O': nw}, 'reference_commit': SHA,
                  'assumption': 'Two explicit DMSO is a prescribed sensitivity model, not a proven complete solvation shell.'})
            jobs.append({'job_id': jobid, 'species_id': sid, 'compound_id': s['compound_id'],
                         'scenario_id': scenario, 'start_index': i, 'initial_release': i <= release,
                         'status': 'PREPARED_NOT_RUN', 'job_directory': folder.relative_to(ROOT).as_posix(),
                         'input_sha256': digest(folder/'goat.inp'), 'xyz_sha256': digest(folder/'start.xyz'),
                         'atom_map_sha256': digest(folder/'atom_map.json'), 'cluster_mol_sha256': digest(folder/'cluster.mol'),
                         'provenance_sha256': digest(folder/'provenance.json'),
                         'atom_count': cluster.GetNumAtoms(), 'solute_atom_count': analyte.GetNumAtoms(),
                         'charge': s['charge'], 'multiplicity': s['multiplicity'],
                         'DMSO_d6_count': nd, 'H2O_count': nw, 'orca_runtime_validated': False})
    print(f'Prepared {sid}: {len(jobs)} initial candidates', flush=True)
    return jobs, scenarios

def build(args):
    source, reference = Path(args.source), Path(args.reference)
    products = readj(source/'registry/products.json')
    allspecies = readj(source/'registry/species.json')
    active = [s for s in allspecies if s['solvent'] == 'DMSO' and s.get('default_run_enabled')]
    if len(products) != 64 or len(active) != 67:
        raise RuntimeError('Unexpected source registry scope')
    files = ['README.md', 'config/stage_methods.json', 'config/adaptive_goat_policy.json',
             'config/qrrho_and_conformer_policy.json', 'config/reference_validation_policy.json',
             'scripts/prepare_explicit_goat_jobs.py']
    provenance = []
    for name in files:
        put(ROOT/'reference_51'/name, (reference/name).read_bytes())
        provenance.append({'path': name, 'sha256': digest(reference/name)})
    putj(ROOT/'registry/source_provenance.json', {'reference_commit': SHA,
         'repository': 'https://github.com/zqk-7k/DK_Computation', 'reference_files': provenance,
         'source_v1': str(source), 'source_species_sha256': digest(source/'registry/species.json'),
         'source_products_sha256': digest(source/'registry/products.json')})
    putj(ROOT/'registry/products.json', products)
    putj(ROOT/'registry/species.json', active)
    putj(ROOT/'registry/held_species.json', [{'species': s, 'reason': 'Optional sample-state/solvent branch not silently activated'} for s in allspecies if s not in active])
    for name in ('held_models.json', 'stereochemistry_evidence.json', 'geometry_evidence.json'):
        if (source/'registry'/name).exists():
            put(ROOT/'registry'/name, (source/'registry'/name).read_bytes())
    plan = readj(reference/'config/stage_methods.json')
    # Work 2 v3 explicitly excludes this branch; the original reference remains
    # a read-only historical source, never the active execution protocol.
    plan['stages'].pop('revdsd_validation', None)
    plan.update({'name': '64_18O_ADAPTED_FROM_51_DESIGN', 'source_commit': SHA,
                 'implementation_status': 'GOAT_STARTING_INPUTS_AND_SINGLE_JOB_RUNNER_ONLY',
                 'downstream_status': 'NOT_COMMISSIONED_FOR_EXPLICIT_18O_CLUSTERS',
                 'original_51_code_configuration_drift': 'See docs/流程差异与未完成事项.txt',
                 'minimum_required_adaptations': ['18O_and_D_exact_mass', 'analyte_only_NMR', 'fixed_composition_separation',
                    'same_method_reference', 'real_GOAT_ensemble_selection', 'no_blind_v1_cluster_ingestion']})
    putj(ROOT/'config/stage_methods_51_aligned.json', plan)
    putj(ROOT/'config/microsolvation_policy.json', {'scenarios': SCENARIOS,
         'first_wave_starts': [3,4,7], 'minimum_success_not_convergence': [2,3,5],
         'max_candidates_by_tier': CEILINGS, 'DMSO_H_isotope': 2, 'trace_water': 'ordinary_H2_16O',
         'analyte_exchangeable_H': 'retained_as_H; no_unjustified_automatic_HD_exchange',
         'temperature_K': 298.15, 'temperature_status': 'PROVISIONAL',
         'cross_composition_Boltzmann_pooling': False, 'TMS': 'separate_reference_per_scenario_and_method',
         'explicit_2_DMSO_d6_is_not': ['measured_coordination_number','complete_solvation_shell','accuracy_guarantee'],
         'one_water_is_not': 'measured_water_concentration', 'quantum_calculations_started': False})
    if args.only:
        active = [s for s in active if s['species_id'] in args.only or s['compound_id'] in args.only]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(lambda s: build_species(s, source, args.full), active))
    jobs = sorted([j for js, _ in results for j in js], key=lambda j: j['job_id'])
    scenarios = [s for _, ss in results for s in ss]
    suffix = '_subset' if args.only else ''
    putj(ROOT/f'registry/jobs{suffix}.json', jobs)
    putj(ROOT/f'registry/scenarios{suffix}.json', scenarios)
    report = {'status': 'INITIAL_STRUCTURES_PREPARED_NO_QM', 'products': len({s['compound_id'] for s in active if s['compound_id']!='TMS'}),
              'solute_model_branches_excluding_reference': sum(s['compound_id']!='TMS' for s in active),
              'reference_models': sum(s['compound_id']=='TMS' for s in active),
              'fixed_composition_scenarios_including_TMS': len(scenarios), 'jobs': len(jobs),
              'wave1_jobs': sum(j['initial_release'] for j in jobs), 'held_jobs': sum(not j['initial_release'] for j in jobs),
              'sampling_converged': False, 'ORCA_GOAT_run': False, 'DFT_run': False, 'NMR_spectra_generated': False}
    putj(ROOT/f'reports/preparation{suffix}.json', report)
    lines = ['64个产物显式溶剂配置清单（计算假设；未做ORCA优化）',
             '每个模型：DMSO-d6连续介质；0DMSO/0水、2DMSO/0水、2DMSO/1普通水三个独立分支。',
             '初始作业数：3、4、7；不是已收敛构象数。温度暂定298.15 K。', '']
    names = {p['compound_id']: p['name'] for p in products}
    for s in active:
        lines.append(f"{s['compound_id']} | {names.get(s['compound_id'], 'TMS reference')} | {s['species_id']} | q={s['charge']} | 0/0;2/0;2/1")
    put(ROOT/f'64个产物_溶剂配置清单{suffix}.txt', '\n'.join(lines)+'\n')
    print(json.dumps(report, indent=2), flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', required=True)
    p.add_argument('--reference', required=True)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--full', action='store_true', help='Prepare held candidate starts too; never executes them')
    p.add_argument('--only', nargs='+')
    build(p.parse_args())
