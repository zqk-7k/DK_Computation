"""Reconstruct one complete alpha-PW12O40 anion from published coordinates.

Source: Deblonde and Colliard, Acta Cryst. E 80, 721-724 (2024).
DOI 10.1107/S2056989024005565, CCDC 2361833.
Coordinates and cell are transcribed from the article HTML table, NOT the CIF.
https://journals.iucr.org/e/issues/2024/07/00/ox2004/index.html
The crystal contains Cs/Tl/Cl counterions; these are excluded deliberately.
No geometry optimization is claimed. The 30 A CIF box is storage only.
"""
from pathlib import Path
import json
import numpy as np
from pymatgen.core import Lattice, Structure, Molecule
from pymatgen.io.cif import CifWriter

ROOT = Path(__file__).resolve().parent
SITES = [
    ('P', 'P1', [.25, .75, .75]),
    ('W', 'W1', [.35861, .64139, .74467]),
    ('O', 'O2', [.3797, .7015, .7985]),
    ('O', 'O3', [.3161, .6010, .6839]),
    ('O', 'O1', [.4101, .5899, .7606]),
    ('O', 'O4', [.2910, .7090, .7090]),
]

def main():
    crystal = Structure.from_spacegroup(
        225, Lattice.cubic(22.9963), [r[0] for r in SITES], [r[2] for r in SITES]
    )
    p = next(site for site in crystal if site.specie.symbol == 'P')
    rows = [('P', np.zeros(3))] + [
        (n.specie.symbol, np.array(n.coords) - p.coords)
        for n in crystal.get_neighbors(p, 6.2) if n.specie.symbol in ('W', 'O')
    ]
    rows.sort(key=lambda r: ({'P': 0, 'W': 1, 'O': 2}[r[0]], np.linalg.norm(r[1]), *r[1]))
    mol = Molecule([r[0] for r in rows], [r[1] for r in rows], charge=-3)
    assert len(mol) == 53
    assert mol.composition.as_dict() == {'P': 1., 'W': 12., 'O': 40.}
    (ROOT / 'alpha_PW12O40_crystallographic.xyz').write_text(mol.to(fmt='xyz'), encoding='utf-8')
    boxed = Structure(Lattice.cubic(30), mol.species, mol.cart_coords + 15, coords_are_cartesian=True)
    CifWriter(boxed).write_file(ROOT / 'alpha_PW12O40_crystallographic.cif')
    audit = []
    for i, site in enumerate(mol):
        if site.specie.symbol != 'W':
            continue
        ns = sorted((float(np.linalg.norm(t.coords-site.coords)), t.specie.symbol, j)
                    for j, t in enumerate(mol) if j != i)
        first = [x for x in ns if x[0] < 2.65]
        assert len(first) == 6 and all(x[1] == 'O' for x in first)
        audit.append({
            'W_index_zero_based': i,
            'first_shell': [{'element': x[1], 'distance_A': round(x[0],6), 'index_zero_based': x[2]} for x in first],
            'nearest_other_W_A': round(min(x[0] for x in ns if x[1]=='W'), 6),
        })
    info = {
        'doi': '10.1107/S2056989024005565',
        'url': 'https://journals.iucr.org/e/issues/2024/07/00/ox2004/index.html',
        'ccdc': '2361833',
        'method': 'Reconstructed from published HTML fractional coordinate table and Fm-3m symmetry; not a downloaded CIF. Cs/Tl/Cl excluded. Salt-derived crystallographic geometry; not optimized catalyst geometry.',
        'cell_a_angstrom': 22.9963,
        'space_group_number': 225,
        'asymmetric_sites': [{'element': r[0], 'label': r[1], 'frac': r[2]} for r in SITES],
        'formula': 'PW12O40', 'heavy_atoms': 53, 'formal_anion_charge': -3,
        'periodic_box_is_storage_only_A': 30, 'no_hydrogen_or_counterions': True,
        'first_shell_audit': audit,
    }
    (ROOT/'alpha_PW12O40_provenance.json').write_text(json.dumps(info, indent=2), encoding='utf-8')
    # A geometric truncation, NOT a ligand-complete molecule or DFT input.
    centre_i = audit[0]['W_index_zero_based']
    fragment_ids = [centre_i] + [r['index_zero_based'] for r in audit[0]['first_shell']]
    fragment_species = [mol[i].specie.symbol for i in fragment_ids]
    fragment_coords = [mol[i].coords - mol[centre_i].coords for i in fragment_ids]
    xyz = ['7', 'Geometric WO6 fragment from crystallographic POM; diagnostic only; no physical charge assigned']
    xyz.extend(f'{el} {v[0]:.10f} {v[1]:.10f} {v[2]:.10f}' for el,v in zip(fragment_species, fragment_coords))
    (ROOT/'W1O6_first_shell_DIAGNOSTIC.xyz').write_text('\n'.join(xyz)+'\n', encoding='utf-8')
    fragbox = Structure(Lattice.cubic(30), fragment_species, np.array(fragment_coords)+15, coords_are_cartesian=True)
    CifWriter(fragbox).write_file(ROOT/'W1O6_first_shell_DIAGNOSTIC.cif')

    # Verify element-resolved distances to every atom, not only the first shell.
    radial = []
    for row in audit:
        c = mol[row['W_index_zero_based']].coords
        radial.append({el: sorted(float(np.linalg.norm(x.coords-c)) for x in mol if x.specie.symbol==el)
                       for el in ('P','W','O')})
    max_deviation = max(max(abs(a-b) for a,b in zip(radial[0][el],r[el]))
                        for r in radial for el in ('P','W','O'))
    assert max_deviation < 1e-8
    edges = [(row['W_index_zero_based'], o['index_zero_based']) for row in audit for o in row['first_shell']]
    oxygen_degrees = {i:sum(i==b for _,b in edges) for i,s in enumerate(mol) if s.specie.symbol=='O'}
    degree_histogram = {str(n):sum(d==n for d in oxygen_degrees.values()) for n in sorted(set(oxygen_degrees.values()))}
    assert degree_histogram == {'1':12,'2':24,'3':4}
    manifest = {
        'model_scope':'Isolated full POM reference geometry plus first-shell geometric diagnostic; neither is the supported W single-atom catalyst.',
        'coordinate_order':'XYZ indices zero-based. CIF writer may reorder sites; locate absorber by coordinates/labels.',
        'nonperiodic':True,
        'pbc_forbidden_for_feff':True,
        'full_model':{'xyz':'alpha_PW12O40_crystallographic.xyz','cif':'alpha_PW12O40_crystallographic.cif',
                      'n_atoms':53,'composition':{'P':1,'W':12,'O':40},'formal_ion_charge':-3,
                      'W_indices_zero_based':[r['W_index_zero_based'] for r in audit],
                      'representative_W_index_zero_based':centre_i,
                      'element_resolved_full_radial_equivalence_max_deviation_A':max_deviation,
                      'all_W_in_single_crystallographic_orbit':True,
                      'W_O_edges':edges,'O_number_of_W_neighbors_histogram':degree_histogram},
        'diagnostic_model':{'xyz':'W1O6_first_shell_DIAGNOSTIC.xyz','cif':'W1O6_first_shell_DIAGNOSTIC.cif',
                            'n_atoms':7,'W_index_zero_based':0,'parent_atom_indices_zero_based':fragment_ids,
                            'physical_charge':None,'relaxation_allowed':False,
                            'limitation':'Removed bridging partners leave dangling O ligands; not a stable molecular model. Only geometric/scattering-path diagnostic. Do not use for quantitative XANES or assign as a catalyst structure.'}
    }
    (ROOT/'reference_model_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print('Validated PW12O40: 53 atoms, all 12 W have six O first-shell neighbors.')

if __name__ == '__main__':
    main()
