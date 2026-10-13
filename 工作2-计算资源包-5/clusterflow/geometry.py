"""Fixed-order full-cluster coordinates and fail-closed geometric checks.

Distance gates detect likely graph changes; they do not infer a new chemical
species or prove solution speciation. Solvent molecules remain in the model.
"""
from __future__ import annotations
from collections import defaultdict
import math
import numpy as np
from rdkit import Chem

MASS_18O = 17.99915961286
MASS_D = 2.01410177812


def canonical(mol):
    copy = Chem.RemoveHs(Chem.Mol(mol))
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(copy, isomericSmiles=True)


def parse_xyz_frames(text):
    lines = text.splitlines()
    frames = []
    offset = 0
    while offset < len(lines):
        if not lines[offset].strip():
            offset += 1
            continue
        try:
            count = int(lines[offset])
        except ValueError as exc:
            raise ValueError('Invalid XYZ atom count') from exc
        if count < 1 or offset + count + 2 > len(lines):
            raise ValueError('Incomplete/empty XYZ frame')
        elements, positions = [], []
        for row in lines[offset + 2:offset + count + 2]:
            fields = row.split()
            if len(fields) != 4:
                raise ValueError('Plain four-column XYZ required; isotope masses belong in inline ORCA coordinates')
            elements.append(fields[0])
            positions.append([float(value) for value in fields[1:]])
        positions = np.asarray(positions, dtype=float)
        if not np.isfinite(positions).all():
            raise ValueError('Nonfinite XYZ coordinates')
        frames.append((elements, positions, '\n'.join(lines[offset:offset + count + 2]) + '\n'))
        offset += count + 2
    if not frames:
        raise ValueError('No XYZ frames')
    return frames


def parse_xyz(text):
    frames = parse_xyz_frames(text)
    if len(frames) != 1:
        raise ValueError('Exactly one XYZ frame required')
    return frames[0][:2]


def with_xyz(template, text):
    if template is None:
        raise ValueError('Missing MOL topology')
    elements, positions = parse_xyz(text)
    if elements != [a.GetSymbol() for a in template.GetAtoms()]:
        raise ValueError('Atom count/order/element mismatch; no automatic remapping')
    copy = Chem.Mol(template)
    copy.RemoveAllConformers()
    conformer = Chem.Conformer(copy.GetNumAtoms())
    for index, xyz in enumerate(positions):
        conformer.SetAtomPosition(index, xyz)
    conformer.Set3D(True)
    copy.AddConformer(conformer)
    return copy


def xyz_text(mol, comment='REVIEWED ORCA GEOMETRY; NOT EXPERIMENTAL DATA'):
    positions = mol.GetConformer().GetPositions()
    lines = [str(mol.GetNumAtoms()), comment]
    for atom, xyz in zip(mol.GetAtoms(), positions):
        lines.append(f'{atom.GetSymbol():2s} {xyz[0]: .10f} {xyz[1]: .10f} {xyz[2]: .10f}')
    return '\n'.join(lines) + '\n'


def isotope_mass(atom):
    number = atom.GetIsotope()
    if not number:
        return None
    if atom.GetSymbol() == 'O' and number == 18:
        return MASS_18O
    if atom.GetSymbol() == 'H' and number == 2:
        return MASS_D
    mass = Chem.GetPeriodicTable().GetMassForIsotope(atom.GetAtomicNum(), number)
    if not math.isfinite(mass) or mass <= 0:
        raise ValueError('Unsupported isotope')
    return mass


def validate_atom_map(mol, rows):
    if mol is None or not isinstance(rows, list) or len(rows) != mol.GetNumAtoms():
        return ['atom_map_length_or_template_invalid']
    issues = []
    for index, (atom, row) in enumerate(zip(mol.GetAtoms(), rows)):
        if not isinstance(row, dict):
            issues.append(f'atom_map_row_invalid:{index}')
            continue
        expected = {'index_zero_based': index, 'index_one_based': index + 1,
                    'element': atom.GetSymbol(), 'isotope': atom.GetIsotope(),
                    'formal_charge': atom.GetFormalCharge()}
        if any(row.get(k) != v for k, v in expected.items()):
            issues.append(f'atom_map_identity_mismatch:{index}')
        if sorted(row.get('neighbours_zero_based', [])) != sorted(n.GetIdx() for n in atom.GetNeighbors()):
            issues.append(f'atom_map_topology_mismatch:{index}')
        if not isinstance(row.get('environment_atom'), bool) or not isinstance(row.get('include_in_analyte_nmr'), bool):
            issues.append(f'atom_map_scope_missing:{index}')
        if not isinstance(row.get('exchangeable_H'), bool) or not row.get('component'):
            issues.append(f'atom_map_component_or_exchange_missing:{index}')
        if row.get('environment_atom') and row.get('include_in_analyte_nmr'):
            issues.append(f'environment_nmr_leak:{index}')
        if row.get('component') == 'analyte' and row.get('environment_atom'):
            issues.append(f'analyte_scope_mismatch:{index}')
        environment = row.get('component') != 'analyte'
        if row.get('environment_atom') is not environment:
            issues.append(f'component_environment_scope_mismatch:{index}')
        if row.get('source_atom_index_zero_based') != (None if environment else index):
            issues.append(f'source_atom_index_mismatch:{index}')
        if row.get('include_in_analyte_nmr') is not (not environment and atom.GetSymbol() in ('H', 'C')):
            issues.append(f'analyte_nmr_scope_mismatch:{index}')
        expected_mass = isotope_mass(atom)
        actual_mass = row.get('isotope_mass_u')
        if expected_mass is None:
            if actual_mass is not None:
                issues.append(f'unexpected_mass:{index}')
        elif not isinstance(actual_mass, (int, float)) or not math.isclose(actual_mass, expected_mass, abs_tol=5e-11, rel_tol=0):
            issues.append(f'isotope_mass_mismatch:{index}')
    # Component labels must describe exactly the immutable disconnected graph.
    for fragment in Chem.GetMolFrags(mol):
        components = {rows[i].get('component') for i in fragment if isinstance(rows[i], dict)}
        if len(components) != 1:
            issues.append('covalent_graph_crosses_component_boundary')
    component_count = len({row.get('component') for row in rows if isinstance(row, dict)})
    if component_count != len(Chem.GetMolFrags(mol)):
        issues.append('component_labels_do_not_match_graph_fragments')
    return issues


def validate_geometry(template, candidate, atom_map, config=None):
    issues = validate_atom_map(template, atom_map)
    if issues:
        return issues
    if candidate is None or candidate.GetNumConformers() != 1:
        return ['candidate_geometry_missing']
    if [(a.GetSymbol(), a.GetIsotope(), a.GetFormalCharge()) for a in template.GetAtoms()] != [
            (a.GetSymbol(), a.GetIsotope(), a.GetFormalCharge()) for a in candidate.GetAtoms()]:
        return ['atom_order_isotope_or_charge_mismatch']
    if canonical(candidate) != canonical(template):
        issues.append('immutable_graph_mismatch')
    copy = Chem.Mol(candidate)
    for bond in copy.GetBonds():
        bond.SetBondDir(Chem.BondDir.NONE)
    Chem.AssignStereochemistryFrom3D(copy, confId=0, replaceExistingTags=True)
    if canonical(copy) != canonical(template):
        issues.append('3D_stereochemistry_mismatch')
    positions = candidate.GetConformer().GetPositions()
    if not np.isfinite(positions).all():
        return issues + ['nonfinite_geometry']
    table = Chem.GetPeriodicTable()
    distances = np.linalg.norm(positions[:, None, :] - positions[None, :, :], axis=2)
    for bond in template.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        radius = sum(table.GetRcovalent(template.GetAtomWithIdx(k).GetAtomicNum()) for k in (i, j))
        if not 0.57 * radius < distances[i, j] < 1.65 * radius:
            issues.append(f'covalent_bond_distance:{i}:{j}:{distances[i,j]:.3f}')
    for i in range(candidate.GetNumAtoms()):
        atom = candidate.GetAtomWithIdx(i)
        for j in range(i):
            if template.GetBondBetweenAtoms(i, j):
                continue
            radius = table.GetRcovalent(atom.GetAtomicNum()) + table.GetRcovalent(candidate.GetAtomWithIdx(j).GetAtomicNum())
            if distances[i, j] < 0.55 * radius:
                issues.append(f'nonbonded_clash:{i}:{j}')
            if atom.GetSymbol() != 'H' and candidate.GetAtomWithIdx(j).GetSymbol() != 'H' and distances[i, j] < 1.10 * radius:
                scope = 'intercomponent' if atom_map[i]['component'] != atom_map[j]['component'] else 'intracomponent'
                issues.append(f'possible_{scope}_covalent_bond:{i}:{j}')
        if atom.GetSymbol() == 'H' and atom.GetDegree() == 1:
            parent = atom.GetNeighbors()[0].GetIdx()
            old_distance = distances[i, parent]
            for other in candidate.GetAtoms():
                j = other.GetIdx()
                if j == parent or other.GetSymbol() not in {'N', 'O', 'S', 'P'}:
                    continue
                if distances[i, j] < 1.30 and distances[i, j] + 0.10 < old_distance:
                    issues.append(f'possible_proton_transfer:{i}:{parent}->{j}')
    groups = defaultdict(list)
    for i, row in enumerate(atom_map):
        if row['element'] != 'H':
            groups[row['component']].append(i)
    if 'analyte' not in groups:
        issues.append('analyte_component_missing')
    else:
        policy = (config or {}).get('geometry', {})
        cutoff = float(policy.get('max_component_contact_A', 4.5))
        reached = {'analyte'}
        while True:
            newly = {name for name, indices in groups.items() if name not in reached and any(
                np.min(distances[np.ix_(indices, groups[seen])]) <= cutoff for seen in reached)}
            if not newly:
                break
            reached |= newly
        issues.extend('dissociated_component:' + name for name in groups if name not in reached)
    return sorted(set(issues))


def isotope_inline_coordinates(mol, atom_map=None):
    if atom_map is not None:
        issues = validate_atom_map(mol, atom_map)
        if issues:
            raise ValueError('Atom/isotope map gate: ' + str(issues))
    lines = []
    for atom, xyz in zip(mol.GetAtoms(), mol.GetConformer().GetPositions()):
        line = f'{atom.GetSymbol():2s} {xyz[0]: .10f} {xyz[1]: .10f} {xyz[2]: .10f}'
        mass = isotope_mass(atom)
        if mass is not None:
            line += f' M = {mass:.11f}'
        lines.append(line)
    return '\n'.join(lines)
