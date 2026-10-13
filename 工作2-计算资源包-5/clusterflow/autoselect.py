"""Source-bound automatic GOAT screening, not scientific approval or convergence.

No quantum execution. Every completed source is revalidated before import and
each subsequent stage. Candidates remain candidates until DFT checks complete.
"""
from __future__ import annotations
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from .common import config, readj, sha, writej
from .geometry import validate_geometry, with_xyz

STATUS = 'AUTOMATIC_SCREENED_CANDIDATE_FOR_DFT'
EH_TO_KCAL = 627.5094740631
NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                     separators=(',', ':'), allow_nan=False).encode('utf-8')).hexdigest()


def _settings(root):
    path = Path(root) / 'config/automation.json'
    value = readj(path)
    if value.get('automatic_handoff', {}).get('enabled') is not True:
        raise ValueError('Automatic handoff is not explicitly enabled')
    supplied = value.get('selection', {})
    result = {'rmsd_A': float(supplied.get('rmsd_A', .10)),
              'representative_energy_window_kcal_mol': float(supplied.get('representative_energy_window_kcal_mol', 6.0)),
              'algorithm_version': '5.0.1', 'energy_source': 'ACTUAL_SOURCE_BOUND_GOAT_ENSEMBLE'}
    if any(not math.isfinite(result[k]) or result[k] <= 0 for k in ('rmsd_A', 'representative_energy_window_kcal_mol')):
        raise ValueError('Invalid automatic selection setting')
    return result


def parse_frame_energies(frames, output):
    """Require explicit Eh in each comment, or corroborate native plain values.

    Native bare XYZ energies are accepted ONLY if the completed GOAT final
    summary declares an Eh minimum and kcal/mol relative energies for EVERY
    corresponding frame. Units are never inferred from magnitude or order.
    """
    comments = [frame[2].splitlines()[1].strip() for frame in frames]
    energies, bare = [], False
    for comment in comments:
        matches = re.findall(rf'(?<![\w.])({NUMBER})\s*(?:Eh|hartree|a\.u\.)(?!\w)', comment, re.I)
        if len(matches) == 1:
            value = float(matches[0].replace('D', 'E').replace('d', 'e'))
        elif not matches and re.fullmatch(rf'(?:Energy\s*[:=]\s*)?({NUMBER})', comment, re.I):
            bare_match = re.fullmatch(rf'(?:Energy\s*[:=]\s*)?({NUMBER})', comment, re.I)
            value = float(bare_match.group(1).replace('D', 'E').replace('d', 'e'))
            bare = True
        else:
            raise ValueError('Ambiguous/missing per-frame GOAT energy or explicit units: ' + comment)
        if not math.isfinite(value):
            raise ValueError('Nonfinite GOAT frame energy')
        energies.append(value)
    if bare:
        minima = list(re.finditer(rf'Lowest energy conformer\s*:\s*({NUMBER})\s+Eh', output, re.I))
        if not minima:
            raise ValueError('Bare XYZ energies require a unit-explicit final GOAT summary')
        match = minima[-1]
        minimum = float(match.group(1).replace('D', 'E').replace('d', 'e'))
        # Only the final ensemble table immediately before the final minimum.
        prefix = output[:match.start()]
        headers = list(re.finditer(r'(?:Erel|E\s*rel|Relative\s+energy).*?kcal\s*/\s*mol[^\n]*', prefix, re.I))
        native_headers = list(re.finditer(r'#\s*Final ensemble info\s*#\s*\n[^\n]*Conformer\s+Energy\s+Degen\.[^\n]*\n\s*\(kcal/mol\)', prefix, re.I))
        headers += native_headers
        if not headers:
            raise ValueError('Final GOAT table lacks explicit relative-energy units')
        table = prefix[max(headers, key=lambda item: item.end()).end():]
        relative = {}
        for line in table.splitlines():
            row = re.fullmatch(rf'\s*(\d+)\s+({NUMBER})\s+\d+\s+{NUMBER}\s+{NUMBER}\s*', line)
            if row:
                index = int(row.group(1))
                if index in relative:
                    raise ValueError('Ambiguous duplicate frame in final GOAT energy table')
                relative[index] = float(row.group(2).replace('D', 'E').replace('d', 'e'))
        if set(relative) == set(range(len(frames))):
            relative = {index + 1: energy for index, energy in relative.items()}
        elif set(relative) != set(range(1, len(frames) + 1)):
            raise ValueError('Incomplete final GOAT energy table; cannot bind bare frame energies')
        for index, energy in enumerate(energies, 1):
            # Final ORCA summary often prints the minimum to six decimals.
            if abs((energy - minimum) * EH_TO_KCAL - relative[index]) > .003:
                raise ValueError('XYZ frame energy disagrees with unit-explicit final GOAT table')
    return energies


def _motif(mol, amap):
    """Geometric proxy only: fixed-order solute torsions + solvent contact sites.

    Solvent component labels are sorted away. This is NOT a complete physical
    hydrogen-bond or ring-pucker classifier, so coverage always remains advisory.
    """
    xyz = mol.GetConformer().GetPositions()
    analyte = [i for i, row in enumerate(amap) if not row['environment_atom'] and row['element'] != 'H']
    groups = defaultdict(list)
    for i, row in enumerate(amap):
        if row['environment_atom'] and row['element'] != 'H':
            groups[row['component']].append(i)
    contacts = []
    for name, ids in groups.items():
        kind = 'DMSO_d6' if name.startswith('DMSO') else 'H2O'
        local = []
        for i in ids:
            if amap[i]['element'] not in ('N', 'O', 'S'):
                continue
            distances = np.linalg.norm(xyz[analyte] - xyz[i], axis=1)
            j = analyte[int(np.argmin(distances))]
            local.append([amap[i]['element'], j, int(float(min(distances)) / .75)])
        contacts.append([kind, sorted(local)])
    torsions = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i not in analyte or j not in analyte:
            continue
        left = sorted(a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors() if a.GetIdx() != j and a.GetIdx() in analyte)
        right = sorted(a.GetIdx() for a in mol.GetAtomWithIdx(j).GetNeighbors() if a.GetIdx() != i and a.GetIdx() in analyte)
        if left and right and left[0] != right[0]:
            angle = rdMolTransforms.GetDihedralDeg(mol.GetConformer(), left[0], i, j, right[0])
            if math.isfinite(angle):
                torsions.append([i, j, int((round(angle, 6) + 180.0) % 360.0 / 60.0)])
    detail = {'solute_torsion_bins': torsions, 'solvent_contact_sites': sorted(contacts)}
    return _digest(detail)[:20], detail


def screen_batch(root, job_ids, target_count, *, previous_selected_ids=None, write=False):
    from .pipeline import _goat_source, cluster_rms
    root = Path(root).resolve()
    policy = _settings(root)
    if not isinstance(target_count, int) or isinstance(target_count, bool) or target_count < 1:
        raise ValueError('Positive integer candidate target required')
    if not isinstance(job_ids, (list, tuple)) or not job_ids or len(set(job_ids)) != len(job_ids):
        raise ValueError('Nonempty unique completed GOAT batch required')
    previous = sorted(set(previous_selected_ids or []))
    sources = [_goat_source(root, job_id) for job_id in sorted(job_ids)]
    first = sources[0]
    species, scenario = first['row']['species_id'], first['row']['scenario_id']
    if any((s['row']['species_id'], s['row']['scenario_id'], s['bindings']['atom_map_sha256']) !=
           (species, scenario, first['bindings']['atom_map_sha256']) for s in sources):
        raise ValueError('Cannot mix species, solvent scenarios, stereoisomers or atom maps')
    rows, molecules, excluded = [], {}, []
    for source in sources:
        job_id = source['row']['job_id']
        energies = parse_frame_energies(source['frames'], (source['directory'] / 'goat.out').read_text(encoding='utf-8', errors='replace'))
        for index, (frame, energy) in enumerate(zip(source['frames'], energies), 1):
            cid = 'g' + hashlib.sha256(job_id.encode()).hexdigest()[:10] + f'_f{index:04d}'
            mol = with_xyz(source['template'], frame[2])
            errors = validate_geometry(source['template'], mol, source['atom_map'], config(root))
            row = {'conformer_id': cid, 'goat_job_id': job_id, 'frame_index_one_based': index,
                   'energy_Eh': energy, 'energy_kind': 'GOAT_POTENTIAL_ENERGY_NOT_DFT_GIBBS',
                   'xyz_sha256': hashlib.sha256(frame[2].encode('utf-8')).hexdigest(),
                   'geometry_sha256': _digest(frame[1].tolist())}
            if errors:
                excluded.append({**row, 'geometry_errors': errors, 'disposition': 'REQUIRES_REVIEW_NOT_AUTOMATICALLY_PROMOTED'})
                continue
            motif_id, motif = _motif(mol, source['atom_map'])
            row.update({'motif_id': motif_id, 'motif_proxy': motif, 'status': STATUS,
                        'human_review_completed': False, 'review_required': True, 'DFT_minimum_verified': False})
            rows.append(row)
            molecules[cid] = mol
    ordered = sorted(rows, key=lambda row: (row['energy_Eh'], row['conformer_id']))
    unique, duplicates, representatives = [], [], {}
    for row in ordered:
        duplicate = None
        for old in unique:
            rms = cluster_rms(molecules[row['conformer_id']], molecules[old['conformer_id']], first['atom_map'])
            if rms <= policy['rmsd_A']:
                duplicate = old
                duplicates.append({'removed': row['conformer_id'], 'retained': old['conformer_id'],
                                   'removed_goat_job_id': row['goat_job_id'], 'retained_goat_job_id': old['goat_job_id'],
                                   'rmsd_A': round(float(rms), 8), 'energy_difference_kcal_mol': round((row['energy_Eh'] - old['energy_Eh']) * EH_TO_KCAL, 8),
                                   'policy': 'FULL_CLUSTER_FIXED_ANALYTE_INDICES_WHOLE_SOLVENT_PERMUTATIONS'})
                break
        representative = duplicate or row
        representatives[row['conformer_id']] = representative['conformer_id']
        if duplicate is None:
            unique.append(row)
    if any(cid not in representatives for cid in previous):
        raise ValueError('Previous selected candidate is absent or no longer passes screening')
    by_id = {row['conformer_id']: row for row in unique}
    selected_ids = []
    for cid in previous:
        mapped = representatives[cid]
        if mapped not in selected_ids:
            selected_ids.append(mapped)
    if len(selected_ids) > target_count:
        raise ValueError('Cumulative target cannot discard previously selected basins')
    if unique and unique[0]['conformer_id'] not in selected_ids and len(selected_ids) >= target_count:
        raise ValueError('Increase cumulative candidate target to retain previous basins AND the new lowest-energy basin')
    minimum = unique[0]['energy_Eh'] if unique else None
    for row in unique:
        row['relative_energy_kcal_mol'] = round((row['energy_Eh'] - minimum) * EH_TO_KCAL, 8)
        row['automatic_DFT_energy_eligible'] = row['relative_energy_kcal_mol'] <= policy['representative_energy_window_kcal_mol']
    if unique and unique[0]['conformer_id'] not in selected_ids and len(selected_ids) < target_count:
        selected_ids.append(unique[0]['conformer_id'])
    covered = {by_id[cid]['motif_id'] for cid in selected_ids}
    for row in unique:
        if len(selected_ids) >= target_count:
            break
        if row['motif_id'] not in covered and row['relative_energy_kcal_mol'] <= policy['representative_energy_window_kcal_mol']:
            selected_ids.append(row['conformer_id'])
            covered.add(row['motif_id'])
    for row in unique:
        if len(selected_ids) >= target_count:
            break
        if row['conformer_id'] not in selected_ids and row['automatic_DFT_energy_eligible']:
            selected_ids.append(row['conformer_id'])
    selected = [by_id[cid] for cid in selected_ids]
    reserve = [row for row in unique if row['conformer_id'] not in selected_ids]
    discovered_by = defaultdict(set)
    for row in rows:
        discovered_by[representatives[row['conformer_id']]].add(row['goat_job_id'])
    covered = {row['motif_id'] for row in selected}
    uncovered = sorted({row['motif_id'] for row in unique} - covered)
    report = {'schema_version': '5.0', 'status': STATUS, 'species_id': species, 'scenario_id': scenario,
              'job_ids': sorted(job_ids), 'target_count': target_count, 'previous_selected_ids': previous,
              'target_kind': 'GOAT_CANDIDATE_BATCH_SIZE_NOT_DFT_STABLE_MINIMUM_QUOTA',
              'policy': policy, 'source_bindings': [s['bindings'] for s in sources],
              'available_unique_count': len(unique), 'target_reached': len(selected) >= target_count,
              'energy_eligible_count': sum(row['automatic_DFT_energy_eligible'] for row in unique),
              'candidates': unique, 'selected': selected, 'reserve': reserve, 'duplicates': duplicates,
              'excluded': excluded, 'raw_frame_count': sum(len(s['frames']) for s in sources),
              'previous_to_current_representatives': {cid: representatives[cid] for cid in previous},
              'rediscovery': [{'conformer_id': cid, 'independent_goat_jobs': sorted(jobs), 'job_count': len(jobs)} for cid, jobs in sorted(discovered_by.items())],
              'motif_coverage': {'observed_proxy_count': len({row['motif_id'] for row in unique}),
                                 'selected_proxy_count': len(covered), 'uncovered_proxy_ids': uncovered,
                                 'complete_physical_motif_coverage_established': False},
              'advisories': ['Candidate target is an initial DFT budget, not a stopping or reliability criterion.',
                            'No GOAT task is terminated when the candidate target is reached.',
                            'GOAT potential energies are not final DFT Gibbs energies; beyond-window candidates remain in reserve for manual review, never forced to fill a quota.',
                            'Contact/torsion motif proxies are incomplete and require scientific inspection.',
                            'Internal solvent automorphisms are not exhaustively deduplicated; inspect possible methyl/water orientation duplicates.'],
              'human_review_completed': False, 'sampling_convergence_established': False,
              'NMR_stability_established': False, 'quantum_jobs_launched': 0}
    if uncovered:
        report['advisories'].append('Some observed geometric motif proxies are not in this wave; they remain in reserve.')
    if any(row['relative_energy_kcal_mol'] > policy['representative_energy_window_kcal_mol'] for row in selected):
        report['advisories'].append('Previously selected basin moved outside the GOAT energy focus after new sampling; retained as cumulative comparison, not newly promoted to fill the quota.')
    report['selection_sha256'] = _digest(report)
    if write:
        from .pipeline import _putj
        _putj(root / 'reports' / ('automatic_selection_' + report['selection_sha256'] + '.json'), report)
    return report


def validate_selection(root, report):
    if isinstance(report, (str, Path)):
        report = readj(report)
    if not isinstance(report, dict) or report.get('status') != STATUS:
        raise ValueError('Expected automatic screening record, not human approval')
    fresh = screen_batch(root, report['job_ids'], report['target_count'], previous_selected_ids=report.get('previous_selected_ids'), write=False)
    if fresh != report:
        raise ValueError('Automatic selection/source/policy changed or report was tampered with')
    return fresh


def import_selection(root, report):
    """Prepare immutable imports only after full re-screen; never forge approval."""
    from .pipeline import _goat_source, _index_path, _put, _relative
    root = Path(root).resolve()
    report = validate_selection(root, report)
    if not report['target_reached'] or not report['selected']:
        raise ValueError('Candidate target not reached; continue the next completed GOAT batch or record a manual scientific decision')
    species, scenario = report['species_id'], report['scenario_id']
    index_path = _index_path(root, species, scenario)
    if index_path.is_file():
        old = readj(index_path)
        if old.get('automatic_selection_sha256') != report['selection_sha256']:
            raise ValueError('Immutable imports already exist; use a new cumulative-wave workspace')
        from .pipeline import _validated_imports
        _validated_imports(root, species, scenario)
        return {'status': 'AUTOMATIC_SELECTION_ALREADY_IMPORTED', 'imported': len(old['conformers']), 'conformer_ids': [r['id'] for r in old['conformers']], 'quantum_jobs_launched': 0}
    if list((root / 'jobs').glob('*/' + species + '/' + scenario + '/*/*/job.json')):
        raise ValueError('Downstream jobs already exist; use a new cumulative-wave workspace')
    base = index_path.parent
    report_path = base / 'automatic_selection.json'
    payloads = {report_path: (json.dumps(report, ensure_ascii=False, indent=2) + '\n').encode('utf-8')}
    selected_by_job = defaultdict(list)
    for row in report['selected']:
        selected_by_job[row['goat_job_id']].append(row)
    records, source_records = [], []
    for job_id, selections in sorted(selected_by_job.items()):
        source = _goat_source(root, job_id)
        payloads[base / 'atom_map.json'] = (source['directory'] / 'atom_map.json').read_bytes()
        source_tag = hashlib.sha256(job_id.encode()).hexdigest()[:10]
        review_path = base / ('source_' + source_tag) / 'review.json'
        review = {'schema_version': '5.0', 'status': STATUS, 'review_type': 'AUTOMATIC_MACHINE_SCREEN_NOT_HUMAN_APPROVAL',
                  **source['bindings'], 'human_review_completed': False, 'review_required': True,
                  'sampling_convergence_established': False, 'selection_sha256': report['selection_sha256'],
                  'selection_report': _relative(root, report_path),
                  'selected_frames': [{'frame_index_one_based': r['frame_index_one_based'], 'conformer_id': r['conformer_id']} for r in selections]}
        review_payload = (json.dumps(review, ensure_ascii=False, indent=2) + '\n').encode('utf-8')
        payloads[review_path] = review_payload
        source_records.append({**source['bindings'], 'review': _relative(root, review_path),
                               'review_sha256': hashlib.sha256(review_payload).hexdigest(), 'sampling_review_claim': False,
                               'review_type': 'AUTOMATIC_MACHINE_SCREEN_NOT_HUMAN_APPROVAL'})
        for row in selections:
            cid, number = row['conformer_id'], row['frame_index_one_based']
            frame_text = source['frames'][number - 1][2]
            mol = with_xyz(source['template'], frame_text)
            xyz_payload = frame_text.encode('utf-8')
            mol_payload = Chem.MolToMolBlock(mol).encode('utf-8')
            payloads[base / (cid + '.xyz')] = xyz_payload
            payloads[base / (cid + '.mol')] = mol_payload
            records.append({'id': cid, 'goat_job_id': job_id, 'frame_index_one_based': number,
                            'xyz': _relative(root, base / (cid + '.xyz')), 'mol': _relative(root, base / (cid + '.mol')),
                            'xyz_sha256': hashlib.sha256(xyz_payload).hexdigest(), 'mol_sha256': hashlib.sha256(mol_payload).hexdigest(),
                            'selection_status': STATUS, 'human_review_completed': False, 'review_required': True})
    index = {'species_id': species, 'scenario_id': scenario, 'sources': source_records, 'conformers': records,
             'atom_map': _relative(root, base / 'atom_map.json'), 'atom_map_sha256': report['source_bindings'][0]['atom_map_sha256'],
             'automatic_selection': _relative(root, report_path), 'automatic_selection_sha256': report['selection_sha256'],
             'sampling_convergence_established': False, 'sampling_status': 'NOT_ESTABLISHED_BY_AUTOMATIC_SCREENING',
             'human_review_completed': False}
    # Collision/integrity preflight before ANY artifact is written.
    payloads[index_path] = (json.dumps(index, ensure_ascii=False, indent=2) + '\n').encode('utf-8')
    for path, payload in payloads.items():
        if path.exists() and path.read_bytes() != payload:
            raise ValueError('Refusing partial/overwriting automatic import: ' + str(path))
    for path, payload in payloads.items():
        _put(path, payload)
    return {'status': 'IMPORTED_AUTOMATIC_SCREENED_CANDIDATES', 'species_id': species, 'scenario_id': scenario,
            'imported': len(records), 'conformer_ids': [r['id'] for r in records], 'quantum_jobs_launched': 0}


machine_select = screen_batch
import_machine_selection = import_selection
