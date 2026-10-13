"""Inspect a cluster XYZ and render an isotope-aware coordinate block ONLY.

This is NOT a GOAT-selection approval, DFT input generator, or NMR executor.
It never reads optimized coordinates from a different atom ordering and never
adds environment atoms to the analyte H/C masks.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from rdkit import Chem

def load_cluster(job_directory, xyz_path):
    d = Path(job_directory)
    m = Chem.MolFromMolBlock((d/'cluster.mol').read_text(encoding='utf-8'), removeHs=False)
    amap = json.loads((d/'atom_map.json').read_text(encoding='utf-8'))
    lines = Path(xyz_path).read_text(encoding='utf-8').splitlines()
    n = int(lines[0])
    if n != m.GetNumAtoms() or len(amap) != n or len(lines) != n+2:
        raise ValueError('Require exactly one complete cluster XYZ with identical atom count')
    coords = []
    for i, (a, row) in enumerate(zip(m.GetAtoms(), lines[2:])):
        p = row.split()
        if len(p) != 4 or p[0] != a.GetSymbol():
            raise ValueError('Element/order/XYZ column mismatch')
        record = amap[i]
        if record['index_zero_based'] != i or record['element'] != a.GetSymbol() or record['isotope'] != a.GetIsotope():
            raise ValueError('Map/known topology mismatch')
        coords.append([float(v) for v in p[1:]])
    if not np.isfinite(coords).all():
        raise ValueError('Nonfinite coordinates')
    return m, amap, np.asarray(coords)

def coordinate_block(m, amap, coords, charge=0, multiplicity=1):
    if charge != sum(a.GetFormalCharge() for a in m.GetAtoms()):
        raise ValueError('Charge differs from frozen cluster')
    rows = [f'* xyz {charge} {multiplicity}']
    for a, record, p in zip(m.GetAtoms(), amap, coords):
        line = f'{a.GetSymbol():2s} {p[0]: .10f} {p[1]: .10f} {p[2]: .10f}'
        isotope = a.GetIsotope()
        expected = 17.99915961286 if a.GetSymbol() == 'O' and isotope == 18 else 2.01410177812 if a.GetSymbol() == 'H' and isotope == 2 else None
        if isotope and expected is None:
            raise ValueError('Unsupported isotope')
        if record['isotope_mass_u'] != expected:
            raise ValueError('Isotope mass corrupted')
        if expected is not None:
            line += f' M={expected:.11f}'
        rows.append(line)
    return '\n'.join(rows+['*'])+'\n'

def analyte_masks(amap):
    masks = {e: [] for e in ('H','C','F')}
    for r in amap:
        if r['environment_atom']:
            if r['include_in_analyte_nmr']:
                raise ValueError('Solvent atom must not appear in analyte NMR')
            continue
        if r['element'] in masks:
            masks[r['element']].append(r['index_zero_based'])
    pairs = []
    for e1, e2 in (('H','H'), ('H','F'), ('C','F')):
        pairs.extend([i, j] for i in masks[e1] for j in masks[e2] if e1 != e2 or i < j)
    return {'zero_based_analyte_nuclei': masks, 'HH_HF_CF_pairs': pairs,
            'no_environment_pairs': True, 'note': 'Masks only; do not imply a commissioned property/collect pipeline.'}

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('job_directory')
    p.add_argument('xyz')
    p.add_argument('--charge', type=int, default=0)
    p.add_argument('--multiplicity', type=int, default=1)
    args = p.parse_args()
    mol, amap, coords = load_cluster(args.job_directory, args.xyz)
    print('# COORDINATE/MASK INSPECTION ONLY; NOT APPROVED PRODUCTION INPUT')
    print(coordinate_block(mol, amap, coords, args.charge, args.multiplicity))
    print(json.dumps(analyte_masks(amap), indent=2))
