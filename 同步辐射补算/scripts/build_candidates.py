#!/usr/bin/env python3
"""Build explicit, unoptimized oxygen-coordinated W hypothesis structures.

The original adatom CIFs are audit references. No Co/Fe neighbours are silently
changed into oxygen. Two distinct substitutional hypotheses are constructed on
the pristine LDH and must pass unconstrained local relaxation before FEFF.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
import warnings

import numpy as np
from pymatgen.core import Structure
from pymatgen.io.cif import CifParser, CifWriter

ROOT = Path(__file__).resolve().parent.parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def safe_child(root, value):
    path = (Path(root) / value).resolve()
    if not path.is_relative_to(Path(root).resolve()):
        raise ValueError('Configured path leaves workflow root')
    return path


def load_cif(path):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        s = CifParser(str(path), frac_tolerance=0, occupancy_tolerance=1.0).parse_structures(primitive=False)[0]
    if not s.is_ordered:
        raise ValueError('Partial occupancy must be resolved before model construction')
    return s, [str(w.message) for w in caught]


def unit(vector):
    vector = np.asarray(vector, dtype=float)
    norm = np.linalg.norm(vector)
    if norm < 1e-10:
        raise ValueError('Degenerate geometry direction')
    return vector / norm


def displacement(lattice, xyz_a, xyz_b):
    fa, fb = lattice.get_fractional_coords(xyz_a), lattice.get_fractional_coords(xyz_b)
    _, image = lattice.get_distance_and_image(fa, fb)
    return lattice.get_cartesian_coords(fb + image - fa)


def distance(lattice, a, b):
    return float(np.linalg.norm(displacement(lattice, a['xyz'], b['xyz'])))


def get_atom(atoms, ident):
    return next(a for a in atoms if a['id'] == ident)


def outward_oh(lattice, atoms, metal_id):
    metal = get_atom(atoms, metal_id)
    candidates = []
    for oxygen in atoms:
        if oxygen['element'] != 'O':
            continue
        delta = displacement(lattice, metal['xyz'], oxygen['xyz'])
        if not (1.35 < np.linalg.norm(delta) < 2.65 and delta[2] > 0.25):
            continue
        hydrogens = [h for h in atoms if h['element'] == 'H' and distance(lattice, oxygen, h) < 1.25]
        if len(hydrogens) == 1:
            candidates.append((oxygen['id'], hydrogens[0]['id'], float(np.linalg.norm(delta))))
    return sorted(candidates, key=lambda x: (x[2], x[0]))


def choose_sites(lattice, atoms, w_id, co_id):
    w_sites, co_sites = outward_oh(lattice, atoms, w_id), outward_oh(lattice, atoms, co_id)
    options = []
    for w, co, dep in itertools.product(w_sites, co_sites, co_sites):
        if len({w[0], co[0], dep[0]}) != 3 or len({w[1], co[1], dep[1]}) != 3:
            continue
        options.append((sum(x[2] for x in [w, co, dep]), w, co, dep))
    if not options:
        raise ValueError('No three distinct outward OH groups for paired W/Co ligand exchange and Co deprotonation')
    _, w, co, dep = min(options)
    return {'eg_at_w': w, 'eg_at_co': co, 'co_deprotonation': dep}


def perpendicular(vector, preferred=None):
    vector = unit(vector)
    preferred = np.array([0., 0., 1.]) if preferred is None else np.asarray(preferred)
    candidate = preferred - np.dot(preferred, vector) * vector
    if np.linalg.norm(candidate) < 1e-7:
        preferred = np.array([1., 0., 0.])
        candidate = preferred - np.dot(preferred, vector) * vector
    return unit(candidate)


def carbon_hydrogens(center, bond1, bond2):
    b1, b2 = unit(bond1), unit(bond2)
    axis = unit(np.cross(b1, b2))
    common = -0.5 * (b1 + b2)
    scale = math.sqrt(max(0., 1. - float(np.dot(common, common))))
    return [np.asarray(center) + 1.09 * unit(common + sign * scale * axis) for sign in [-1, 1]]


def eg_geometry(lattice, oxygen, metal, phi):
    """Retain the donor O and build O-CH2-CH2-OH; donor proton is removed."""
    o1 = np.asarray(metal['xyz']) + displacement(lattice, metal['xyz'], oxygen['xyz'])
    outward = unit(o1 - np.asarray(metal['xyz']))
    v1 = unit(outward + np.array([0., 0., 1.]))
    p, q = perpendicular(v1), unit(np.cross(v1, perpendicular(v1)))
    v2 = unit(v1 / 3. + math.sqrt(8. / 9.) * (math.cos(phi) * p + math.sin(phi) * q))
    c1, c2 = o1 + 1.43 * v1, o1 + 1.43 * v1 + 1.52 * v2
    v3 = unit(v2 / 3. + math.sqrt(8. / 9.) * perpendicular(v2))
    o2 = c2 + 1.43 * v3
    h_oxygen = o2 + 0.97 * unit(v3 / 3. + math.sqrt(8. / 9.) * perpendicular(v3))
    points = [('C', c1), ('C', c2), ('O', o2)]
    points += [('H', h) for h in carbon_hydrogens(c1, o1 - c1, c2 - c1)]
    points += [('H', h) for h in carbon_hydrogens(c2, c1 - c2, o2 - c2)]
    points += [('H', h_oxygen)]
    return [dict(id=f'EG_{i+1}', element=el, xyz=list(map(float, xyz)), origin='constructed_EG_alkoxide')
            for i, (el, xyz) in enumerate(points)]


def pair_threshold(a, b):
    if a['element'] == b['element'] == 'H':
        return 0.55
    if 'H' in {a['element'], b['element']}:
        return 0.65
    return 0.85


def clash_report(lattice, atoms):
    violations, minimum = [], float('inf')
    fractional = lattice.get_fractional_coords(np.array([a['xyz'] for a in atoms]))
    matrix = lattice.get_all_distances(fractional, fractional)
    for i, a in enumerate(atoms):
        for j in range(i+1, len(atoms)):
            b = atoms[j]
            d = float(matrix[i, j])
            minimum = min(minimum, d)
            if d < pair_threshold(a, b):
                violations.append(dict(atom1=a['id'], atom2=b['id'], distance_A=d, minimum_A=pair_threshold(a, b)))
    return dict(passed=not violations, minimum_pair_distance_A=minimum, violations=violations)


def attach_eg(lattice, atoms, donor, metal_id):
    oxygen_id, hydrogen_id, _ = donor
    oxygen, metal = get_atom(atoms, oxygen_id), get_atom(atoms, metal_id)
    retained = [a.copy() for a in atoms if a['id'] != hydrogen_id]
    trials = []
    for angle in np.linspace(0, 2*math.pi, 24, endpoint=False):
        added = eg_geometry(lattice, oxygen, metal, angle)
        candidate = retained + added
        clash = clash_report(lattice, candidate)
        # Pick an unclashed trial with the greatest adsorbate-to-host separation.
        clearance = float(np.min(lattice.get_all_distances(
            lattice.get_fractional_coords([a['xyz'] for a in added]),
            lattice.get_fractional_coords([a['xyz'] for a in retained if a['id'] != oxygen_id]))))
        trials.append((len(clash['violations']), -clearance, float(angle), candidate, clash))
    _, _, angle, candidate, clash = min(trials, key=lambda x: x[:3])
    return candidate, dict(donor_oxygen_id=oxygen_id, removed_H_id=hydrogen_id,
                          retained_original_O=True, added_formula='C2H5O',
                          replacing_ligand='OH -> OCH2CH2OH (mono-deprotonated EG)',
                          torsion_radians=angle, clash_check=clash)


def coordination(lattice, atoms, w_id, cutoff=2.6):
    w = get_atom(atoms, w_id)
    neighbours = [dict(id=a['id'], element=a['element'], distance_A=distance(lattice, w, a))
                  for a in atoms if a['id'] != w_id and distance(lattice, w, a) <= cutoff]
    neighbours.sort(key=lambda x: x['distance_A'])
    heavy = [a for a in neighbours if a['element'] != 'H']
    return dict(cutoff_A=cutoff, neighbours=neighbours,
                all_heavy_neighbours_are_O=bool(heavy) and all(a['element'] == 'O' for a in heavy),
                oxygen_count=sum(a['element'] == 'O' for a in heavy))


def build(config_path):
    config_path = Path(config_path).resolve(); root = config_path.parent
    config = read_json(config_path); options = config.get('candidates', {})
    source = safe_child(root, options.get('source_cif', 'structures_original/CoFe-LDH.cif'))
    structure, parser_warnings = load_cif(source)
    base = [dict(id=f'orig_{i}', element=site.specie.symbol, xyz=list(map(float, site.coords)),
                 origin='source_CIF', original_index=i, original_label=site.label)
            for i, site in enumerate(structure)]
    first_fe = next(a for a in base if a['element'] == 'Fe')
    first_co = min((a for a in base if a['element'] == 'Co'), key=lambda a: (distance(structure.lattice, first_fe, a), a['id']))
    anchors = [('W_sub_Fe', first_fe), ('W_sub_Co', first_co)]
    if options.get('anchors'):
        anchors = [(a['id'], base[int(a['source_index'])]) for a in options['anchors']]
    manifest = []
    for anchor_name, original in anchors:
        if original['element'] not in {'Fe', 'Co'}:
            raise ValueError('Substitution anchor must be an actual Co or Fe atom')
        anchor = [a.copy() for a in base]
        get_atom(anchor, original['id'])['element'] = 'W'
        response_co = min((a for a in anchor if a['element'] == 'Co'),
                          key=lambda a: (distance(structure.lattice, original, a), a['id']))
        sites = choose_sites(structure.lattice, anchor, original['id'], response_co['id'])
        for background, adsorbate in itertools.product(['hydroxylated', 'Co_deprotonated'], ['none', 'EG_at_W', 'EG_at_Co']):
            ident = f'{anchor_name}__{background}__{adsorbate}'
            if not all(c.isalnum() or c in '_-' for c in ident):
                raise ValueError('Unsafe model identifier')
            atoms = [a.copy() for a in anchor]; removed = []
            if background == 'Co_deprotonated':
                removed.append(sites['co_deprotonation'][1])
                atoms = [a for a in atoms if a['id'] not in removed]
            exchange = None
            if adsorbate != 'none':
                key = 'eg_at_w' if adsorbate == 'EG_at_W' else 'eg_at_co'
                metal = original['id'] if adsorbate == 'EG_at_W' else response_co['id']
                atoms, exchange = attach_eg(structure.lattice, atoms, sites[key], metal)
                removed.append(exchange['removed_H_id'])
            clashes = clash_report(structure.lattice, atoms)
            coord = coordination(structure.lattice, atoms, original['id'])
            if not clashes['passed'] or not coord['all_heavy_neighbours_are_O']:
                raise ValueError(f'Candidate {ident} has clashes or non-O initial W first shell; inspect construction')
            output = safe_child(root, f'candidates/{ident}'); output.mkdir(parents=True, exist_ok=True)
            candidate = Structure(structure.lattice, [a['element'] for a in atoms], [a['xyz'] for a in atoms],
                                  coords_are_cartesian=True, labels=[a['id'] for a in atoms])
            cif = output/'structure.cif'
            candidate_text = str(CifWriter(candidate, symprec=None))
            if cif.exists() and cif.read_text(encoding='utf-8') != candidate_text:
                raise ValueError(f'Changed candidate would overwrite {cif}; use a new project or model identifier')
            cif.write_text(candidate_text, encoding='utf-8', newline='\n')
            # CifParser sorts species. Retain stable labels and reorder the map
            # to the actual parsed input order used by QE, not construction order.
            parsed_candidate, _ = load_cif(cif)
            by_id = {a['id']: a for a in atoms}
            if set(site.label for site in parsed_candidate) != set(by_id):
                raise ValueError('CIF round trip lost the stable atom labels')
            atoms = [dict(by_id[site.label], xyz=list(map(float, site.coords))) for site in parsed_candidate]
            metadata = dict(id=ident, status='candidate_not_DFT_optimized', source_cif=source.relative_to(root).as_posix(),
                            source_sha256=sha(source), structure_sha256=sha(cif), anchor=anchor_name,
                            substitution=dict(atom_id=original['id'], from_element=original['element'], to_element='W'),
                            background=background, adsorbate=adsorbate, response_co_id=response_co['id'],
                            ligand_sites=sites, removed_original_atoms=removed, ligand_exchange=exchange,
                            original_composition=dict(Counter(a['element'] for a in base)),
                            composition=dict(Counter(a['element'] for a in atoms)), atom_map=atoms,
                            atom_map_order='CifParser input order; stable atom labels retained',
                            initial_coordination=coord, initial_clash_check=clashes, parser_warnings=parser_warnings,
                            potential_assignment=None, fixed_lattice=True, geometry_constraints='none on atomic positions',
                            scientific_scope='Alternative substitutional hypothesis; not a correction of the original optimized adatom structure. Background and adsorption factors are crossed; no spectral response is imposed.')
            write_json(output/'metadata.json', metadata)
            manifest.append(dict(id=ident, path=output.relative_to(root).as_posix(), structure_sha256=sha(cif),
                                 status=metadata['status'], anchor=anchor_name, background=background, adsorbate=adsorbate,
                                 w_atom_id=original['id'], response_co_id=response_co['id']))
    write_json(root/'candidates/manifest.json', manifest)
    print(f'Prepared {len(manifest)} explicit candidate structures; DFT optimization is still required.')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default=str(ROOT/'config.json'))
    args = parser.parse_args(); build(args.config)


if __name__ == '__main__':
    main()
