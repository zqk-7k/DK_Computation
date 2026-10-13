"""Reviewed GOAT ensembles -> explicit-cluster DFT stages, without execution.

Every scenario is independently addressable. Continuum-only calculations are
not prerequisites for a two-DMSO/one-water branch. No initial seed is promoted.
"""
from __future__ import annotations
from collections import defaultdict
import hashlib
import itertools
import json
from pathlib import Path
import re
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from .common import config, hash_protocol, readj, writej, sha, inside
from .geometry import canonical, parse_xyz_frames, with_xyz, xyz_text, validate_geometry, validate_atom_map
from .inputs import STAGES, J_STAGES, analyte_nuclei, selected_j_pairs, render_input, method_fingerprint
from .gates import inspect_job, boltzmann, receipt_errors
from .orcaio import parse_output

REVIEW_CHECKS = ('identity', 'stereochemistry', 'composition', 'proton_transfer_absent',
                 'cluster_bound', 'unique_basins_across_all_selected_starts')


def _token(value):
    if not isinstance(value, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,180}', value) or '..' in value:
        raise ValueError('Unsafe species/scenario/conformer identifier')
    return value


def _relative(root, path):
    return Path(path).resolve().relative_to(Path(root).resolve()).as_posix()


def _put(path, data):
    path = Path(path)
    payload = data.encode('utf-8') if isinstance(data, str) else data
    if path.exists():
        if path.read_bytes() != payload:
            raise ValueError('Refusing to overwrite different immutable artifact: ' + str(path))
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def _putj(path, value):
    _put(path, json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def _selection(root, species_id, scenario_id):
    _token(species_id)
    _token(scenario_id)
    p = config(root)
    if scenario_id not in p['independent_scenarios']:
        raise ValueError('Unknown independent scenario')
    matches = [s for s in readj(Path(root) / 'registry/species.json') if s['species_id'] == species_id]
    if len(matches) != 1:
        raise ValueError('Use one exact registered species_id')
    return matches[0], p


def _goat_source(root, job_id):
    root = Path(root).resolve()
    rows = readj(root / 'registry/jobs.json')
    matches = [r for r in rows if r.get('job_id') == job_id]
    if len(matches) != 1:
        raise ValueError('GOAT job ID missing or duplicated')
    row = matches[0]
    species, _ = _selection(root, row['species_id'], row['scenario_id'])
    if any(row.get(k) != species.get(k) for k in ('compound_id', 'charge', 'multiplicity')):
        raise ValueError('GOAT registry chemical-state identity mismatch')
    directory = inside(root, row['job_directory'])
    if any((directory / name).exists() for name in ('.goat-run.lock', '.goat-portable.lock', '.run.lock')):
        raise ValueError('GOAT source has an active runner lock; completed batch required')
    receipt_paths = [directory / name for name in ('run_receipt.json', 'receipt.json') if (directory / name).is_file()]
    if len(receipt_paths) != 1:
        raise ValueError('Exactly one completed GOAT run receipt required')
    receipt_path = receipt_paths[0]
    receipt = readj(receipt_path)
    if receipt.get('status') not in ('completed', 'COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW'):
        raise ValueError('GOAT receipt is not completed')
    if receipt.get('job_id', receipt.get('id')) != job_id:
        raise ValueError('GOAT receipt job identity mismatch')
    required_inputs = ['goat.inp', 'start.xyz', 'atom_map.json', 'cluster.mol', 'provenance.json']
    errors = receipt_errors(directory, receipt, required_inputs, ['goat.out', 'goat.finalensemble.xyz'])
    fields = {'goat.inp': 'input_sha256', 'start.xyz': 'xyz_sha256', 'atom_map.json': 'atom_map_sha256',
              'cluster.mol': 'cluster_mol_sha256', 'provenance.json': 'provenance_sha256'}
    for name, field in fields.items():
        if not isinstance(row.get(field), str) or not (directory / name).is_file() or sha(directory / name) != row[field].lower():
            errors.append('GOAT registry hash mismatch:' + name)
    if errors:
        raise ValueError('GOAT source binding failed: ' + str(errors))
    output_path = directory / 'goat.out'
    output = output_path.read_text(encoding='utf-8', errors='replace')
    parsed = parse_output(output)
    actual = re.findall(r'Program Version\s+(\d+\.\d+\.\d+)', output)
    if not parsed['normal_termination'] or not actual or set(actual) != {config(root)['orca_version']}:
        raise ValueError('GOAT output normal termination/version gate failed')
    ensemble = directory / 'goat.finalensemble.xyz'
    frames = parse_xyz_frames(ensemble.read_text(encoding='utf-8'))
    template = Chem.MolFromMolBlock((directory / 'cluster.mol').read_text(encoding='utf-8'), removeHs=False)
    amap = readj(directory / 'atom_map.json')
    issues = validate_atom_map(template, amap)
    if issues:
        raise ValueError('GOAT immutable atom map invalid: ' + str(issues))
    composition = re.fullmatch(r'DMSO_(\d+)_W(\d+)_.+', row['scenario_id'])
    if not composition:
        raise ValueError('Unrecognized fixed-composition scenario identity')
    components = {r['component'] for r in amap if r['environment_atom']}
    dmso_count = sum(name.startswith('DMSO_d6_') for name in components)
    water_count = sum(name.startswith('H2O_') for name in components)
    if (dmso_count, water_count) != tuple(map(int, composition.groups())) or dmso_count + water_count != len(components):
        raise ValueError('GOAT cluster composition differs from selected scenario')
    if row.get('DMSO_d6_count') != dmso_count or row.get('H2O_count') != water_count:
        raise ValueError('GOAT registry composition counts disagree with actual cluster')
    fragments = Chem.GetMolFrags(template, asMols=True, sanitizeFrags=True)
    if canonical(fragments[0]) != canonical(Chem.AddHs(Chem.MolFromSmiles(species['smiles']))):
        raise ValueError('GOAT analyte graph differs from registered chemical species')
    symbols = [a.GetSymbol() for a in template.GetAtoms()]
    if any(frame[0] != symbols for frame in frames):
        raise ValueError('GOAT ensemble element/order/composition mismatch')
    bindings = {'goat_job_id': job_id, 'species_id': row['species_id'], 'scenario_id': row['scenario_id'],
                'source_receipt_sha256': sha(receipt_path), 'goat_output_sha256': sha(output_path),
                'ensemble_sha256': sha(ensemble), 'atom_map_sha256': sha(directory / 'atom_map.json'),
                'cluster_mol_sha256': sha(directory / 'cluster.mol')}
    return {'row': row, 'directory': directory, 'receipt_path': receipt_path, 'receipt': receipt,
            'frames': frames, 'template': template, 'atom_map': amap, 'bindings': bindings}


def draft_goat_review(root, job_id):
    """Return a review draft only; requires real hash-bound completed output."""
    source = _goat_source(root, job_id)
    p = config(root)
    frame_checks = []
    for index, frame in enumerate(source['frames'], start=1):
        candidate = with_xyz(source['template'], frame[2])
        frame_checks.append({'frame_index_one_based': index,
                             'automatic_geometry_errors': validate_geometry(source['template'], candidate, source['atom_map'], p)})
    return {'schema_version': '3.0', 'status': 'DRAFT_REQUIRES_REVIEW', **source['bindings'],
            'reviewer': '', 'reviewed_utc': '', 'reason': '', 'selected_frames': [],
            'checks': {key: False for key in REVIEW_CHECKS},
            'sampling_convergence_established': False, 'frame_checks': frame_checks,
            'notice': 'Select unique physical basins across ALL imported starts of this exact species/scenario; no input seed promotion.'}


def _index_path(root, species_id, scenario_id):
    _selection(root, species_id, scenario_id)
    return Path(root) / 'imports' / species_id / scenario_id / 'index.json'


def import_goat(root, job_id, review_path):
    root = Path(root).resolve()
    source = _goat_source(root, job_id)
    review_path = Path(review_path)
    if not review_path.is_absolute():
        review_path = root / review_path
    review = readj(review_path)
    if review.get('status') != 'APPROVED_FOR_PRESCREEN':
        raise ValueError('A separate approved scientific selection review is required')
    if any(review.get(k) != v for k, v in source['bindings'].items()):
        raise ValueError('Scientific review no longer matches GOAT source hashes/identity')
    if any(review.get('checks', {}).get(k) is not True for k in REVIEW_CHECKS):
        raise ValueError('Scientific identity/composition/stereo/bound-cluster/unique-basin checks must all pass')
    if any(not isinstance(review.get(k), str) or not review[k].strip() for k in ('reviewer', 'reviewed_utc', 'reason')):
        raise ValueError('Reviewer, date and reason are required')
    selected = review.get('selected_frames')
    if not isinstance(selected, list) or not selected:
        raise ValueError('No reviewed GOAT frames selected')
    row = source['row']
    species_id, scenario_id = row['species_id'], row['scenario_id']
    for path in (root / 'jobs').glob('*/' + species_id + '/' + scenario_id + '/*/*/job.json'):
        raise ValueError('Import additional starts before preparing downstream jobs; use a new package for a new sampling wave')
    index_path = _index_path(root, species_id, scenario_id)
    index = readj(index_path) if index_path.is_file() else {'species_id': species_id, 'scenario_id': scenario_id, 'sources': [], 'conformers': []}
    if any(s['goat_job_id'] == job_id for s in index['sources']):
        raise ValueError('GOAT source already imported; imports are immutable')
    amap_digest = source['bindings']['atom_map_sha256']
    if index['sources'] and any(s['atom_map_sha256'] != amap_digest for s in index['sources']):
        raise ValueError('Cannot mix atom maps/compositions/atom orders between starts')
    p = config(root)
    planned, seen_frames, seen_ids = [], set(), {c['id'] for c in index['conformers']}
    source_tag = hashlib.sha256(job_id.encode()).hexdigest()[:10]
    for selection in selected:
        number = selection.get('frame_index_one_based')
        if not isinstance(number, int) or isinstance(number, bool) or not 1 <= number <= len(source['frames']) or number in seen_frames:
            raise ValueError('Invalid/duplicate selected frame index')
        seen_frames.add(number)
        cid = selection.get('conformer_id') or f'g{source_tag}_f{number:04d}'
        _token(cid)
        if len(cid) > 64 or cid in seen_ids:
            raise ValueError('Duplicate/overlong conformer identifier')
        seen_ids.add(cid)
        frame_text = source['frames'][number - 1][2]
        mol = with_xyz(source['template'], frame_text)
        errors = validate_geometry(source['template'], mol, source['atom_map'], p)
        if errors:
            raise ValueError(f'Selected GOAT frame {number} fails geometric gate: {errors}')
        planned.append((cid, number, mol, frame_text))
    base = index_path.parent
    _put(base / 'atom_map.json', (source['directory'] / 'atom_map.json').read_bytes())
    source_base = base / ('source_' + source_tag)
    _put(source_base / 'review.json', review_path.read_bytes())
    imported = []
    for cid, number, mol, frame_text in planned:
        _put(base / (cid + '.xyz'), frame_text)
        _put(base / (cid + '.mol'), Chem.MolToMolBlock(mol))
        imported.append({'id': cid, 'goat_job_id': job_id, 'frame_index_one_based': number,
                         'xyz': _relative(root, base / (cid + '.xyz')), 'mol': _relative(root, base / (cid + '.mol')),
                         'xyz_sha256': sha(base / (cid + '.xyz')), 'mol_sha256': sha(base / (cid + '.mol'))})
    index['sources'].append({**source['bindings'], 'review': _relative(root, source_base / 'review.json'),
                             'review_sha256': sha(source_base / 'review.json'),
                             'sampling_review_claim': review.get('sampling_convergence_established') is True})
    index['conformers'] += imported
    index['atom_map'] = _relative(root, base / 'atom_map.json')
    index['atom_map_sha256'] = amap_digest
    index['sampling_convergence_established'] = False
    index['sampling_status'] = 'NOT_ESTABLISHED_BY_REVIEW_BOOLEAN_OR_GEOMETRY_IMPORT'
    writej(index_path, index)
    return {'status': 'IMPORTED_REVIEWED_GOAT_FRAMES', 'species_id': species_id, 'scenario_id': scenario_id,
            'imported': len(imported), 'conformer_ids': [c['id'] for c in imported], 'quantum_jobs_launched': 0}


def _validated_imports(root, species_id, scenario_id):
    root = Path(root).resolve()
    index = readj(_index_path(root, species_id, scenario_id))
    if index.get('species_id') != species_id or index.get('scenario_id') != scenario_id:
        raise ValueError('Import index identity mismatch')
    amap_path = inside(root, index['atom_map'])
    if sha(amap_path) != index['atom_map_sha256']:
        raise ValueError('Imported atom map changed')
    automatic = None
    if index.get('automatic_selection') is not None:
        from .autoselect import validate_selection
        automatic = validate_selection(root, readj(inside(root, index['automatic_selection'])))
        if automatic['selection_sha256'] != index.get('automatic_selection_sha256'):
            raise ValueError('Automatic import selection binding changed')
        expected = {(r['goat_job_id'], r['frame_index_one_based'], r['conformer_id']) for r in automatic['selected']}
        actual = {(r['goat_job_id'], r['frame_index_one_based'], r['id']) for r in index['conformers']}
        if expected != actual or not automatic['target_reached']:
            raise ValueError('Automatic import does not match complete machine selection')
    sources = {}
    for record in index['sources']:
        source = _goat_source(root, record['goat_job_id'])
        if any(source['bindings'][k] != record.get(k) for k in source['bindings']):
            raise ValueError('Imported GOAT source changed')
        review_path = inside(root, record['review'])
        if sha(review_path) != record['review_sha256']:
            raise ValueError('Imported scientific review changed')
        review = readj(review_path)
        if automatic is not None:
            from .autoselect import STATUS
            selected = [{'frame_index_one_based': r['frame_index_one_based'], 'conformer_id': r['conformer_id']}
                        for r in automatic['selected'] if r['goat_job_id'] == record['goat_job_id']]
            if (review.get('status') != STATUS or review.get('review_type') != 'AUTOMATIC_MACHINE_SCREEN_NOT_HUMAN_APPROVAL'
                    or review.get('human_review_completed') is not False or review.get('review_required') is not True
                    or review.get('selection_sha256') != automatic['selection_sha256']
                    or review.get('selected_frames') != selected
                    or inside(root, review.get('selection_report', '')) != inside(root, index['automatic_selection'])):
                raise ValueError('Automatic screening record changed; it is not human approval')
        elif review.get('status') != 'APPROVED_FOR_PRESCREEN' or any(review.get('checks', {}).get(k) is not True for k in REVIEW_CHECKS):
            raise ValueError('Imported scientific review no longer approved')
        if any(review.get(k) != source['bindings'][k] for k in source['bindings']):
            raise ValueError('Imported review source bindings changed')
        if source['bindings']['atom_map_sha256'] != index['atom_map_sha256']:
            raise ValueError('Imported sources disagree on fixed cluster atom map')
        source['approved_selections'] = review['selected_frames']
        sources[record['goat_job_id']] = source
    if not index['conformers'] or len({c['id'] for c in index['conformers']}) != len(index['conformers']):
        raise ValueError('Missing/duplicate imported conformers')
    for record in index['conformers']:
        for key in ('xyz', 'mol'):
            if sha(inside(root, record[key])) != record[key + '_sha256']:
                raise ValueError('Imported coordinate/topology changed')
        source = sources[record['goat_job_id']]
        approved = [r for r in source['approved_selections'] if r['frame_index_one_based'] == record['frame_index_one_based']]
        if len(approved) != 1 or (approved[0].get('conformer_id') and approved[0]['conformer_id'] != record['id']):
            raise ValueError('Import index frame was not selected by its scientific review')
        if (inside(root, record['xyz'])).read_text(encoding='utf-8') != source['frames'][record['frame_index_one_based'] - 1][2]:
            raise ValueError('Imported geometry differs from actual selected GOAT frame')
    return index


def write_job(root, species, scenario_id, stage, cid, mol, amap, source):
    root = Path(root).resolve()
    p = config(root)
    ph = hash_protocol(p)
    _token(cid)
    species_id = species['species_id']
    issues = validate_geometry(mol, mol, amap, p)
    if issues:
        raise ValueError('Job input geometry gate failed: ' + str(issues))
    if canonical(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)[0]) != canonical(Chem.AddHs(Chem.MolFromSmiles(species['smiles']))):
        raise ValueError('Analyte graph differs from registered species')
    solvent = species.get('solvent', p.get('solvent', 'DMSO'))
    text, pairs = render_input(mol, stage, p, solvent, species['charge'], species['multiplicity'], amap)
    folder = root / 'jobs' / ph / species_id / scenario_id / stage / cid
    _put(folder / 'job.inp', text)
    _put(folder / 'input.xyz', xyz_text(mol))
    _put(folder / 'input_structure.mol', Chem.MolToMolBlock(mol))
    # Copy byte-identical source atom map so all stages share its immutable hash.
    source_map = root / 'imports' / species_id / scenario_id / 'atom_map.json'
    if readj(source_map) != amap:
        raise ValueError('Scenario atom map changed')
    _put(folder / 'atom_map.json', source_map.read_bytes())
    nuclei = analyte_nuclei(amap, ('H', 'C') if stage == 'shielding' else ('H', 'C', 'F'))
    metadata = {'id': _relative(root, folder), 'species_id': species_id, 'compound_id': species['compound_id'],
                'scenario_id': scenario_id, 'stage': stage, 'conformer_id': cid, 'protocol_hash': ph,
                'solvent': solvent, 'charge': species['charge'], 'multiplicity': species['multiplicity'],
                'nprocs': p[stage]['nprocs'], 'maxcore_mb': p[stage]['maxcore_mb'],
                'input_sha256': sha(folder / 'job.inp'), 'atom_map_sha256': sha(folder / 'atom_map.json'),
                'expected_j_pairs_zero_based': pairs,
                'expected_nuclei_zero_based': nuclei if stage == 'shielding' or stage in J_STAGES else [],
                'method_fingerprint': method_fingerprint(p, stage, solvent, has_br=any(a.GetSymbol() == 'Br' for a in mol.GetAtoms())), 'source': source,
                'full_cluster_atom_count': mol.GetNumAtoms(), 'orca_version_required': p['orca_version'],
                'sampling_converged': False, 'scientific_scope': 'ONE_FIXED_COMPOSITION_CLUSTER',
                'bromine_mixed_J_basis_requires_validation': stage in J_STAGES and any(a.GetSymbol() == 'Br' for a in mol.GetAtoms())}
    _putj(folder / 'job.json', metadata)
    return folder


def prepare_prescreen(root, species_id, scenario_id):
    root = Path(root).resolve()
    species, p = _selection(root, species_id, scenario_id)
    index = _validated_imports(root, species_id, scenario_id)
    amap = readj(inside(root, index['atom_map']))
    paths = []
    for record in index['conformers']:
        mol = Chem.MolFromMolBlock(inside(root, record['mol']).read_text(encoding='utf-8'), removeHs=False)
        mol = with_xyz(mol, inside(root, record['xyz']).read_text(encoding='utf-8'))
        source = {'stage': 'automatic_screened_goat' if index.get('automatic_selection') else 'reviewed_goat', 'goat_job_id': record['goat_job_id'],
                  'frame_index_one_based': record['frame_index_one_based'], 'xyz_sha256': record['xyz_sha256'],
                  'import_index_sha256': sha(_index_path(root, species_id, scenario_id))}
        paths.append(_relative(root, write_job(root, species, scenario_id, 'prescreen', record['id'], mol, amap, source)))
    return {'status': 'PREPARED', 'stage': 'prescreen', 'species_id': species_id, 'scenario_id': scenario_id,
            'jobs': paths, 'job_count': len(paths), 'quantum_jobs_launched': 0}


def _review_snapshot(check):
    review=check.get('scientific_review')
    if not review or not review.get('accepted'):return None
    return {'file':'scientific_review.json','sha256':review['sha256'],
            'reviewer':review['review']['reviewer'],'limitations':review['review']['limitations'],
            'scientific_status':check['scientific_status']}


def _validate_stage_source(root, directory, metadata, index):
    """Bind each stage's input coordinates to its actual accepted upstream job."""
    root, directory = Path(root).resolve(), Path(directory)
    source = metadata.get('source', {})
    if source.get('stage') in ('reviewed_goat', 'automatic_screened_goat'):
        if source.get('stage') == 'automatic_screened_goat' and not index.get('automatic_selection'):
            raise ValueError('Automatic stage source lacks source-bound machine screening record')
        if source.get('import_index_sha256') != sha(_index_path(root, metadata['species_id'], metadata['scenario_id'])):
            raise ValueError('Prescreen import-index binding changed')
        records = [c for c in index['conformers'] if c['id'] == metadata['conformer_id']]
        if len(records) != 1:
            raise ValueError('Prescreen conformer absent from approved imports')
        record = records[0]
        for key in ('goat_job_id', 'frame_index_one_based', 'xyz_sha256'):
            if source.get(key) != record[key]:
                raise ValueError('Prescreen source selection differs from reviewed GOAT frame')
        target = inside(root, record['xyz'])
    elif source.get('stage') in ('prescreen', 'optfreq'):
        upstream = inside(root, source.get('job_dir', ''))
        intended = root / 'jobs' / metadata['protocol_hash'] / metadata['species_id'] / metadata['scenario_id'] / source['stage'] / metadata['conformer_id']
        if upstream != intended:
            raise ValueError('Stage handoff source path/identity mismatch')
        for filename, field in (('job.inp', 'input_sha256'), ('job.out', 'output_sha256'),
                                ('receipt.json', 'receipt_sha256'), ('job.xyz', 'geometry_sha256')):
            if not (upstream / filename).is_file() or sha(upstream / filename) != source.get(field):
                raise ValueError('Stage handoff source hash changed: ' + filename)
        if source.get('scientific_review') is not None:
            review_path=upstream/'scientific_review.json'
            if not review_path.is_file() or sha(review_path)!=source['scientific_review']['sha256']:
                raise ValueError('Stage handoff scientific review changed')
        target = upstream / 'job.xyz'
    else:
        raise ValueError('Stage lacks an approved upstream geometry source')
    current_frame = parse_xyz_frames((directory / 'input.xyz').read_text(encoding='utf-8'))
    source_frame = parse_xyz_frames(target.read_text(encoding='utf-8'))
    if len(current_frame) != 1 or len(source_frame) != 1 or current_frame[0][0] != source_frame[0][0] or not np.allclose(current_frame[0][1], source_frame[0][1], rtol=0, atol=1e-7):
        raise ValueError('Stage input coordinates differ from reviewed upstream output')


def cluster_rms(mol_a, mol_b, atom_map):
    """All-atom metric; only graph-identical whole solvent components permute.

    Analyte indices (including diastereotopic H) never permute. Internal solvent
    graph automorphisms are deliberately not inferred; ambiguous duplicate
    rotamers still require selection review rather than manufactured weights.
    """
    groups = defaultdict(list)
    for row in atom_map:
        groups[row['component']].append(row['index_zero_based'])
    equivalent = defaultdict(list)
    for component, indices in groups.items():
        if component == 'analyte':
            continue
        local = {idx: pos for pos, idx in enumerate(indices)}
        signature = tuple((atom_map[i]['element'], atom_map[i]['isotope'], atom_map[i]['formal_charge'],
                           tuple(sorted(local[j] for j in atom_map[i]['neighbours_zero_based']))) for i in indices)
        equivalent[signature].append(component)
    options = [list(itertools.permutations(names)) for names in equivalent.values()]
    count = np.prod([len(choice) for choice in options], dtype=np.int64) if options else 1
    if count > 1000:
        raise ValueError('Too many component permutations; reviewed duplicate audit required')
    best = float('inf')
    for permutations in itertools.product(*options):
        mapping = {i: i for i in range(mol_a.GetNumAtoms())}
        for original, permuted in zip(equivalent.values(), permutations):
            for left, right in zip(original, permuted):
                mapping.update(zip(groups[left], groups[right]))
        rms, _ = rdMolAlign.GetAlignmentTransform(mol_a, mol_b, atomMap=sorted(mapping.items()), reflect=False)
        best = min(best, rms)
    return best


def deduplicate(accepted, p, atom_map):
    unique, duplicates = [], []
    policy = p.get('geometry', {})
    geometry_only = bool(p.get('conformer_policy', {}).get('final_stable_target_plan'))
    for entry in sorted(accepted, key=lambda a: a['parsed']['gibbs_free_energy_Eh']):
        match = None
        for old in unique:
            delta = abs(entry['parsed']['gibbs_free_energy_Eh'] - old['parsed']['gibbs_free_energy_Eh']) * 627.5094740631
            if (geometry_only or delta <= policy.get('dedup_energy_kcal', .10)) and cluster_rms(entry['mol'], old['mol'], atom_map) <= policy.get('dedup_rms_A', .10):
                match = old['job']['conformer_id']
                break
        if match:
            duplicate = {'removed': entry['job']['conformer_id'], 'retained': match,
                         'policy': 'full_cluster_RMSD_whole_solvent_permutations_no_energy_count_inflation' if geometry_only else 'full_cluster_RMSD_plus_G_whole_solvent_component_permutations'}
            if geometry_only:
                duplicate.update({'energy_difference_kcal_mol': delta,
                                  'energy_review_recommended': delta > policy.get('dedup_energy_kcal', .10),
                                  'large_energy_review_required': delta > policy.get('duplicate_energy_review_kcal', 1.0)})
            duplicates.append(duplicate)
        else:
            unique.append(entry)
    return unique, duplicates


def final_stable_candidates(root, species_id, scenario_id):
    """Read-only final DFT minimum quota audit; never launches property jobs.

    A human-reviewed imaginary-mode approximation may remain useful elsewhere,
    but never counts as a verified stable minimum here. The energy window is an
    explicit resource policy, not proof of conformational completeness.
    """
    root = Path(root).resolve()
    _, p = _selection(root, species_id, scenario_id)
    index = _validated_imports(root, species_id, scenario_id)
    ph = hash_protocol(p)
    base = root / 'jobs' / ph / species_id / scenario_id
    settings = readj(root / 'config/automation.json')
    stable_policy = settings.get('stable_minima', {})
    maximum = float(stable_policy.get('energy_window_kcal_mol', stable_policy.get('maximum_relative_G_kcal_mol', 6.0)))
    if not np.isfinite(maximum) or maximum <= 0:
        raise ValueError('Invalid final-stable-minimum relative-G window')
    accepted, issues = [], []
    paths = sorted((base / 'optfreq').glob('*/job.json'))
    expected_gate = base / 'prescreen_gate.json'
    expected_ids = set()
    if expected_gate.is_file():
        previous = readj(expected_gate)
        if previous.get('status') == 'advanced':
            expected_ids = set(previous['retained_ids'])
    if {path.parent.name for path in paths} != expected_ids or not expected_ids:
        issues.append({'conformer_id': None, 'errors': ['optfreq_jobs_do_not_match_completed_prescreen_gate']})
    for path in paths:
        check = inspect_job(path.parent, p)
        errors = list(check['errors'])
        if any(check['job'].get(k) != v for k, v in {'stage': 'optfreq', 'species_id': species_id,
               'scenario_id': scenario_id, 'protocol_hash': ph, 'conformer_id': path.parent.name}.items()):
            errors.append('stage_job_identity_mismatch')
        try:
            _validate_stage_source(root, path.parent, check['job'], index)
        except (ValueError, OSError, KeyError) as exc:
            errors.append('source:' + str(exc))
        frequencies = check.get('parsed', {}).get('frequencies_cm1', [])
        if any(value < 0 for value in frequencies):
            errors.append('negative_frequency_not_verified_stable_even_if_human_reviewed_approximation')
        if check.get('scientific_status') == 'HUMAN_REVIEWED_APPROXIMATION':
            errors.append('human_reviewed_approximation_not_verified_stable_minimum')
        if errors:
            issues.append({'conformer_id': path.parent.name, 'errors': sorted(set(errors))})
        else:
            accepted.append(check)
    amap = readj(inside(root, index['atom_map']))
    unique, duplicates = deduplicate(accepted, p, amap)
    large_discrepancies = [row for row in duplicates if row.get('large_energy_review_required')]
    if large_discrepancies:
        issues.append({'conformer_id': None, 'errors': ['large_G_discrepancy_within_same_geometric_basin_requires_scientific_review'],
                       'duplicates': large_discrepancies})
    minimum = min((check['parsed']['gibbs_free_energy_Eh'] for check in unique), default=None)
    rows = []
    for check in unique:
        cid = check['job']['conformer_id']
        g = check['parsed']['gibbs_free_energy_Eh']
        delta = (g - minimum) * 627.5094740631
        directory = base / 'optfreq' / cid
        rows.append({'conformer_id': cid, 'G_Eh': g, 'relative_G_kcal_mol': delta,
                     'eligible': delta <= maximum, 'status': 'DFT_OPTIMIZED_FREQUENCY_CHECKED_UNIQUE_MINIMUM',
                     'output_sha256': sha(directory / 'job.out'), 'receipt_sha256': sha(directory / 'receipt.json'),
                     'geometry_sha256': sha(directory / 'job.xyz'), 'scientific_review': _review_snapshot(check)})
    eligible = [row for row in rows if row['eligible']]
    return {'schema_version': '5.0', 'status': 'FINAL_DFT_MINIMUM_QUOTA_AUDIT_NOT_CONVERGENCE',
            'species_id': species_id, 'scenario_id': scenario_id, 'eligible_count': len(eligible),
            'eligible_ids': [row['conformer_id'] for row in eligible], 'all_unique_ids': [row['conformer_id'] for row in rows],
            'unique_stable_count': len(rows), 'conformers': rows, 'duplicates': duplicates, 'issues': issues,
            'energy_discrepancy_advisories': [row for row in duplicates if row.get('energy_review_recommended')],
            'human_reviewed_approximations': [row['conformer_id'] for row in issues if 'human_reviewed_approximation_not_verified_stable_minimum' in row['errors']],
            'maximum_relative_G_kcal_mol': maximum, 'all_prepared_optfreq_pass': bool(paths) and not issues,
            'sampling_convergence_established': False, 'NMR_stability_established': False,
            'quantum_jobs_launched': 0}


def advance(root, species_id, scenario_id, from_stage):
    root = Path(root).resolve()
    species, p = _selection(root, species_id, scenario_id)
    ph = hash_protocol(p)
    if from_stage not in ('prescreen', 'optfreq'):
        raise ValueError('Advance only from prescreen or optfreq')
    index = _validated_imports(root, species_id, scenario_id)
    amap = readj(inside(root, index['atom_map']))
    base = root / 'jobs' / ph / species_id / scenario_id
    if from_stage == 'optfreq' and (root / 'config/automation.json').is_file() and (root / 'registry/batch_plan_v5.json').is_file():
        automation = readj(root / 'config/automation.json')
        if automation.get('automatic_handoff', {}).get('enabled') is True:
            plan = readj(root / 'registry/batch_plan_v5.json')
            targets = [row for row in plan['species'] if row['species_id'] == species_id and row['scenario_id'] == scenario_id]
            if len(targets) != 1 or type(targets[0].get('stable_conformer_target')) is not int or targets[0]['stable_conformer_target'] < 1:
                raise ValueError('One exact final stable-minimum target is required by the v5 plan')
            stable_audit = final_stable_candidates(root, species_id, scenario_id)
            target = targets[0]['stable_conformer_target']
            if stable_audit['issues'] or stable_audit['eligible_count'] < target:
                return {'status': 'BLOCKED_FINAL_STABLE_TARGET', 'species_id': species_id, 'scenario_id': scenario_id,
                        'from_stage': from_stage, 'stable_conformer_target': target,
                        'eligible_stable_count': stable_audit['eligible_count'], 'stable_minimum_audit': stable_audit,
                        'quantum_jobs_launched': 0, 'prepared_jobs': []}
    if from_stage == 'prescreen':
        expected = {r['id'] for r in index['conformers']}
    else:
        previous = readj(base / 'prescreen_gate.json')
        if previous.get('status') != 'advanced':
            raise ValueError('Accepted prescreen gate required')
        expected = set(previous['retained_ids'])
        # A completed prior gate is rechecked, not treated as permanent approval.
        for cid in previous['evaluated_ids']:
            upstream_check=inspect_job(base / 'prescreen' / cid, p)
            if not upstream_check['accepted']:
                raise ValueError('Upstream prescreen source no longer passes')
            if previous.get('evaluated_scientific_reviews',{}).get(cid)!=_review_snapshot(upstream_check):
                raise ValueError('Upstream prescreen scientific review changed')
    paths = sorted((base / from_stage).glob('*/job.json'))
    if not expected or {path.parent.name for path in paths} != expected:
        raise ValueError('Prepared stage conformers differ from declared upstream selection')
    checks = [inspect_job(path.parent, p) for path in paths]
    for path, check in zip(paths, checks):
        intended = {'species_id': species_id, 'scenario_id': scenario_id, 'protocol_hash': ph,
                    'stage': from_stage, 'conformer_id': path.parent.name}
        if any(check['job'].get(k) != v for k, v in intended.items()):
            check['errors'].append('stage_job_identity_mismatch')
        try:
            _validate_stage_source(root, path.parent, check['job'], index)
        except (ValueError, OSError, KeyError) as exc:
            check['errors'].append('source:' + str(exc))
        check['accepted'] = not check['errors']
    bad = [{'conformer': c['job']['conformer_id'], 'errors': c['errors']} for c in checks if not c['accepted']]
    if bad:
        report = {'status': 'BLOCKED_ALL_PREPARED_CONFORMERS_MUST_PASS', 'species_id': species_id,
                  'scenario_id': scenario_id, 'from_stage': from_stage, 'details': bad}
        writej(base / (from_stage + '_gate.json'), report)
        return report
    unique, duplicates = deduplicate(checks, p, amap)
    outputs = []
    weights = boltzmann([c['parsed']['gibbs_free_energy_Eh'] for c in unique], p['temperature_K']) if from_stage == 'optfreq' else []
    for offset, check in enumerate(unique):
        cid = check['job']['conformer_id']
        parent = base / from_stage / cid
        source = {'stage': from_stage, 'job_dir': _relative(root, parent),
                  'input_sha256': check['job']['input_sha256'], 'output_sha256': sha(parent / 'job.out'),
                  'receipt_sha256': sha(parent / 'receipt.json'), 'geometry_sha256': sha(parent / 'job.xyz'),
                  'scientific_review':_review_snapshot(check)}
        if from_stage == 'prescreen':
            outputs.append(_relative(root, write_job(root, species, scenario_id, 'optfreq', cid, check['mol'], amap, source)))
        else:
            shielding = write_job(root, species, scenario_id, 'shielding', cid, check['mol'], amap, source)
            pairlist = selected_j_pairs(check['mol'], amap)
            j_jobs = {}
            from .policy import j_methods_for
            for branch in j_methods_for(p, species['compound_id']):
                j_jobs[branch] = (_relative(root, write_job(root, species, scenario_id, 'j_' + branch, cid, check['mol'], amap, source))
                                  if species['compound_id'] != 'TMS' and pairlist else None)
            outputs.append({'id': cid, 'weight': weights[offset], 'G_Eh': check['parsed']['gibbs_free_energy_Eh'],
                            'scientific_review':_review_snapshot(check),
                            'optfreq_job': _relative(root, parent), 'shielding_job': _relative(root, shielding),
                            'j_jobs': j_jobs, 'no_j_policy': 'TMS_REFERENCE' if species['compound_id'] == 'TMS' else ('NO_ELIGIBLE_ANALYTE_PAIRS' if not pairlist else None)})
    report = {'status': 'advanced', 'species_id': species_id, 'scenario_id': scenario_id,
              'protocol_hash': ph, 'from_stage': from_stage, 'evaluated_ids': sorted(expected),
              'retained_ids': [c['job']['conformer_id'] for c in unique], 'duplicates': duplicates,
              'input_import_index_sha256': sha(_index_path(root, species_id, scenario_id)),
              'evaluated_scientific_reviews':{c['job']['conformer_id']:_review_snapshot(c) for c in checks},
              'evaluated_receipt_sha256': {c['job']['conformer_id']: sha(base / from_stage / c['job']['conformer_id'] / 'receipt.json') for c in checks}}
    writej(base / (from_stage + '_gate.json'), report)
    if from_stage == 'optfreq':
        ensemble = {'species_id': species_id, 'compound_id': species['compound_id'], 'scenario_id': scenario_id,
                    'protocol_hash': ph, 'solvent': species.get('solvent', p.get('solvent', 'DMSO')),
                    'charge': species['charge'], 'multiplicity': species['multiplicity'], 'temperature_K': p['temperature_K'],
                    'atom_map': index['atom_map'], 'atom_map_sha256': index['atom_map_sha256'], 'conformers': outputs,
                    'optfreq_gate': _relative(root, base / 'optfreq_gate.json'), 'optfreq_gate_sha256': sha(base / 'optfreq_gate.json'),
                    'import_index_sha256': sha(_index_path(root, species_id, scenario_id)),
                    'import_index_snapshot': index,
                    'reference_ensemble': f'ensembles/{ph}/TMS_DMSO/{scenario_id}/ensemble.json',
                    'sampling_convergence_established': False,
                    'sampling_review_claims': [s['sampling_review_claim'] for s in index['sources']],
                    'sampling_review_claim': all(s['sampling_review_claim'] for s in index['sources']),
                    'selection_review_type': 'AUTOMATIC_MACHINE_SCREEN_NOT_HUMAN_APPROVAL' if index.get('automatic_selection') else 'HUMAN_SCIENTIFIC_REVIEW',
                    'human_selection_review_completed': not bool(index.get('automatic_selection')),
                    'weight_model': 'Native final ORCA qRRHO G; fixed composition; screened unique basins; degeneracy=1',
                    'duplicates': duplicates,
                    'warnings': ['Sampling convergence not established by import or successful DFT.',
                                 'Internal solvent automorphisms and conformer degeneracy require scientific review.']}
        _putj(root / 'ensembles' / ph / species_id / scenario_id / 'ensemble.json', ensemble)
    return {**report, 'prepared_jobs': outputs, 'quantum_jobs_launched': 0}


def manifest(root, stage, species_id=None, scenario_id=None):
    if stage not in STAGES:
        raise ValueError('Unknown stage')
    root = Path(root).resolve()
    ph = hash_protocol(config(root))
    rows = []
    for path in sorted((root / 'jobs' / ph).glob('*/*/' + stage + '/*/job.json')):
        row = readj(path)
        if species_id is not None and row['species_id'] != species_id:
            continue
        if scenario_id is not None and row['scenario_id'] != scenario_id:
            continue
        rows.append({'job_dir': _relative(root, path.parent), **row})
    return rows


def validate_ensemble_sources(root, species_id, scenario_id):
    root = Path(root).resolve()
    _, p = _selection(root, species_id, scenario_id)
    ph = hash_protocol(p)
    index = _validated_imports(root, species_id, scenario_id)
    ensemble_path = root / 'ensembles' / ph / species_id / scenario_id / 'ensemble.json'
    ensemble = readj(ensemble_path)
    if any(ensemble.get(k) != v for k, v in {'species_id': species_id, 'scenario_id': scenario_id, 'protocol_hash': ph,
                                            'temperature_K': p['temperature_K']}.items()):
        raise ValueError('Ensemble identity/protocol/temperature mismatch')
    if sha(_index_path(root, species_id, scenario_id)) != ensemble['import_index_sha256']:
        raise ValueError('Ensemble import selection changed')
    if ensemble.get('import_index_snapshot') != index:
        raise ValueError('Embedded import snapshot differs from frozen import index')
    if ensemble['atom_map'] != index['atom_map'] or ensemble['atom_map_sha256'] != index['atom_map_sha256']:
        raise ValueError('Ensemble atom map differs from fixed imported cluster')
    gate_path = inside(root, ensemble['optfreq_gate'])
    if sha(gate_path) != ensemble['optfreq_gate_sha256']:
        raise ValueError('Optimization gate changed')
    gate = readj(gate_path)
    if gate.get('status') != 'advanced' or {c['id'] for c in ensemble['conformers']} != set(gate['retained_ids']):
        raise ValueError('Ensemble differs from accepted optimization selection')
    base = root / 'jobs' / ph / species_id / scenario_id
    for cid in gate['evaluated_ids']:
        check = inspect_job(base / 'optfreq' / cid, p)
        if not check['accepted'] or sha(base / 'optfreq' / cid / 'receipt.json') != gate['evaluated_receipt_sha256'][cid]:
            raise ValueError('Optimization source no longer passes or receipt changed: ' + cid)
        if gate.get('evaluated_scientific_reviews',{}).get(cid)!=_review_snapshot(check):
            raise ValueError('Optimization scientific review changed: '+cid)
        _validate_stage_source(root, base / 'optfreq' / cid, check['job'], index)
    gibbs = []
    for record in ensemble['conformers']:
        directory = inside(root, record['optfreq_job'])
        if directory != base / 'optfreq' / record['id']:
            raise ValueError('Optimization source path/identity mismatch')
        check = inspect_job(directory, p)
        if record.get('scientific_review')!=_review_snapshot(check):
            raise ValueError('Ensemble scientific-review snapshot changed')
        g = check['parsed']['gibbs_free_energy_Eh']
        if abs(g - record['G_Eh']) > 1e-10:
            raise ValueError('Ensemble G differs from actual accepted final output')
        gibbs.append(g)
        for path in [record['shielding_job'], *record['j_jobs'].values()]:
            if path is None:
                continue
            prop = inside(root, path)
            meta = readj(prop / 'job.json')
            if any(meta.get(k) != v for k, v in {'species_id': species_id, 'scenario_id': scenario_id, 'protocol_hash': ph,
                                                'conformer_id': record['id']}.items()):
                raise ValueError('Property job identity mismatch')
            if meta['atom_map_sha256'] != ensemble['atom_map_sha256'] or sha(prop / 'atom_map.json') != ensemble['atom_map_sha256']:
                raise ValueError('Property atom map differs from ensemble')
            if meta.get('source', {}).get('receipt_sha256') != sha(directory / 'receipt.json'):
                raise ValueError('Property source receipt no longer matches optimization')
            if meta.get('source',{}).get('scientific_review')!=_review_snapshot(check):
                raise ValueError('Property scientific-review source changed')
            _validate_stage_source(root, prop, meta, index)
            template = Chem.MolFromMolBlock((prop / 'input_structure.mol').read_text(encoding='utf-8'), removeHs=False)
            actual = with_xyz(template, (prop / 'input.xyz').read_text(encoding='utf-8'))
            if canonical(actual) != canonical(check['mol']) or not np.allclose(actual.GetConformer().GetPositions(), check['mol'].GetConformer().GetPositions(), rtol=0, atol=1e-7):
                raise ValueError('Property coordinates differ from final accepted optimization')
    weights = boltzmann(gibbs, p['temperature_K'])
    if not np.allclose(weights, [c['weight'] for c in ensemble['conformers']], rtol=1e-10, atol=1e-12):
        raise ValueError('Ensemble weights differ from actual native qRRHO Gibbs model')
    return {'accepted': True, 'species_id': species_id, 'scenario_id': scenario_id,
            'conformers': len(ensemble['conformers']), 'ensemble': ensemble}
