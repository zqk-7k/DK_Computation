#!/usr/bin/env python3
"""Optional four-model local-relaxation input preparation; no DFT is run.

Model CIFs are unchanged unoptimized hypotheses. Common host labels use the
same movable mask across the four controls. The outer host stays fixed.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import warnings

import numpy as np
from pymatgen.core import Element
from pymatgen.io.cif import CifParser

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False)+'\n', encoding='utf-8', newline='\n')


def load(path):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return CifParser(path, frac_tolerance=0, occupancy_tolerance=1.0).parse_structures(primitive=False)[0]


def inputs(root):
    result = []
    for record in read(root/'manifest.json'):
        directory = root/record['path']; path = directory/'structure.cif'
        if sha(path) != record['source_sha256']:
            raise ValueError('Source candidate CIF changed: '+record['id'])
        structure = load(path)
        labels = [s.label for s in structure]
        if len(set(labels)) != len(labels):
            raise ValueError('Stable labels must be unique')
        result.append((record, structure, read(directory/'source_metadata.json')))
    return result


def common_mask(models, radius):
    if radius < 3.5:
        raise ValueError('A radius below 3.5 A omits the adjacent metal shell')
    movable_host = set()
    for record, structure, metadata in models:
        labels = [s.label for s in structure]
        centers = [metadata['substitution']['atom_id'], metadata['response_co_id']]
        indices = [labels.index(label) for label in centers]
        distances = structure.lattice.get_all_distances(structure.frac_coords, structure.frac_coords[indices])
        movable_host.update(labels[i] for i in range(len(structure))
                            if labels[i].startswith('orig_') and float(np.min(distances[i])) <= radius)
    # Move the proton attached to each mobile host oxygen so that a mobile OH
    # group is not artificially pinned by its H. Iterate to the common union.
    for _, structure, _ in models:
        labels = [s.label for s in structure]
        distance = structure.distance_matrix
        oxygen = [i for i, site in enumerate(structure) if site.specie.symbol == 'O' and labels[i] in movable_host]
        for i, site in enumerate(structure):
            if site.specie.symbol == 'H' and labels[i].startswith('orig_') and oxygen and np.min(distance[i, oxygen]) < 1.25:
                movable_host.add(labels[i])
    return movable_host


def dipole_window(structure):
    matrix = np.asarray(structure.lattice.matrix)
    if np.max(np.abs(matrix[:2, 2])) > 1e-7 or np.max(np.abs(matrix[2, :2])) > 1e-7:
        raise ValueError('Automatic dipole placement requires an xy slab with c along z')
    c = float(matrix[2, 2]); positions = np.sort(np.mod(structure.frac_coords[:, 2], 1.))
    endpoints = np.r_[0., positions, 1.]
    low, high = max(zip(endpoints[:-1], endpoints[1:]), key=lambda pair: pair[1]-pair[0])
    available = float((high-low)*c-4.)
    if available < .75:
        raise ValueError('Insufficient empty vacuum for a dipole window with 2 A nuclear margins')
    width = min(1.5, available)/c; center = float(low+high)/2
    return dict(emaxpos=center-width/2, eopreg=width, eamp=0., edir=3, nuclear_margin_A=2.)


def prepare(root, output, pseudo_dir, radius=4., masks_only=False, allow_missing=False):
    models = inputs(root); mobile = common_mask(models, radius)
    definitions = read(root/'docs/SSSP_1.3.0_PBE_efficiency.json')
    output.mkdir(parents=True, exist_ok=True)
    records = []
    for source, structure, metadata in models:
        target = output/source['id']; target.mkdir(parents=True, exist_ok=True)
        if (target/'local_mask.json').exists():
            previous = read(target/'local_mask.json')
            if previous.get('radius_A') != radius or previous.get('source_sha256') != source['source_sha256']:
                raise ValueError('Existing local mask uses different geometry/radius; choose a new --output directory')
        labels = [s.label for s in structure]
        flags = [label in mobile or label.startswith('EG_') for label in labels]
        if not any(flags) or all(flags):
            raise ValueError('Local-relax mask must contain both mobile and fixed atoms')
        atoms = [dict(index_zero=i, index_one=i+1, label=site.label, element=site.specie.symbol,
                      cartesian_A=list(map(float, site.coords)), movable=bool(flags[i]),
                      qe_if_pos=[1, 1, 1] if flags[i] else [0, 0, 0])
                 for i, site in enumerate(structure)]
        mask = dict(id=source['id'], source_cif=source['path']+'/structure.cif', source_sha256=source['source_sha256'],
                    status='candidate_not_optimized', radius_A=radius,
                    centers=[metadata['substitution']['atom_id'], metadata['response_co_id']],
                    rule='Union of host atoms within radius of W or response Co across all four models; attached OH hydrogens and all EG atoms mobile; remaining host fixed.',
                    movable_host_labels=sorted(mobile), movable_count=sum(flags), fixed_count=len(flags)-sum(flags),
                    composition=dict(Counter(site.specie.symbol for site in structure)), atoms=atoms,
                    limitation='Constrained local relaxation only. Fixed atom residual forces are not a convergence test on mobile atoms. Validate a larger mobile region and full relaxation before quantitative publication claims.')
        write(target/'local_mask.json', mask)
        header='index_zero\tindex_one\tlabel\telement\tx_A\ty_A\tz_A\tmove_x\tmove_y\tmove_z\n'
        rows=['\t'.join([str(a['index_zero']),str(a['index_one']),a['label'],a['element']]+[f'{v:.10f}' for v in a['cartesian_A']]+[str(v) for v in a['qe_if_pos']]) for a in atoms]
        (target/'atom_constraints.txt').write_text(header+'\n'.join(rows)+'\n',encoding='utf-8',newline='\n')
        xyz=[str(len(structure)),f"Unoptimized candidate {source['id']}; original source SHA256 {source['source_sha256']}"]
        xyz += [a['element']+' '+' '.join(f'{v:.10f}' for v in a['cartesian_A']) for a in atoms]
        (target/'structure.xyz').write_text('\n'.join(xyz)+'\n',encoding='utf-8',newline='\n')
        record=dict(id=source['id'],path=target.relative_to(root).as_posix(),movable_count=sum(flags),fixed_count=len(flags)-sum(flags),source_sha256=source['source_sha256'])
        if not masks_only:
            elements=list(dict.fromkeys(site.specie.symbol for site in structure)); pseudo_records={}
            for element in elements:
                definition=definitions[element];file=pseudo_dir/definition['filename']
                if not file.is_file() and not allow_missing:
                    raise FileNotFoundError(f'Missing pseudopotential: {file}; install the pinned SSSP files or use --allow-missing-pseudos for inspection only')
                if file.is_file() and hashlib.md5(file.read_bytes()).hexdigest()!=definition['md5']:
                    raise ValueError('Pseudopotential digest mismatch: '+file.name)
                pseudo_records[element]=dict(filename=definition['filename'],md5=definition['md5'],installed=file.is_file(),sha256=sha(file) if file.is_file() else None)
            dipole=dipole_window(structure)
            relpseudo=os.path.relpath(pseudo_dir,target).replace(os.sep,'/')
            lines=['&CONTROL'," calculation='relax'"," restart_mode='from_scratch'"," prefix='local_pilot'",f" pseudo_dir='{relpseudo}'"," outdir='./scratch'",' nstep=150',' tprnfor=.true.',' tstress=.true.',' etot_conv_thr=1.0d-5',' forc_conv_thr=1.0d-3',' tefield=.true.',' dipfield=.true.','/','&SYSTEM',' ibrav=0',f' nat={len(structure)}',f' ntyp={len(elements)}',' ecutwfc=90',' ecutrho=1080',' nspin=2',' tot_charge=0'," occupations='smearing'"," smearing='mv'",' degauss=0.02',' nosym=.true.',' noinv=.true.',' edir=3',' eamp=0.0',f" emaxpos={dipole['emaxpos']:.12g}",f" eopreg={dipole['eopreg']:.12g}"]
            for i,element in enumerate(elements,1):lines.append(f" starting_magnetization({i})={dict(Co=.5,Fe=.6).get(element,0.):.8g}")
            lines+=['/','&ELECTRONS',' conv_thr=1.0d-8',' electron_maxstep=250',' mixing_beta=0.2'," diagonalization='david'",'/','&IONS'," ion_dynamics='bfgs'",'/','ATOMIC_SPECIES']
            lines += [f"{el} {float(Element(el).atomic_mass):.10g} {definitions[el]['filename']}" for el in elements]
            lines+=['CELL_PARAMETERS angstrom']+[' '.join(f'{v:.12f}' for v in row) for row in structure.lattice.matrix]+['ATOMIC_POSITIONS angstrom']
            lines += [a['element']+' '+' '.join(f'{v:.12f}' for v in a['cartesian_A'])+' '+' '.join(str(v) for v in a['qe_if_pos']) for a in atoms]
            lines+=['K_POINTS automatic','2 2 1 0 0 0']
            inp=target/'qe.in';text='\n'.join(lines)+'\n'
            if inp.exists() and inp.read_text(encoding='utf-8')!=text:
                raise ValueError('Changed input would overwrite an existing job; choose a new --output directory')
            inp.write_text(text,encoding='utf-8',newline='\n')
            record.update(input_sha256=sha(inp),pseudopotentials=pseudo_records,dipole=dipole,
                          status='input_inspection_only_missing_pseudos' if not all(x['installed'] for x in pseudo_records.values()) else 'prepared_not_run')
            write(target/'job.json',record)
        records.append(record)
    write(output/'local_relax_manifest.json',records)
    return records


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',default='local_relax_R4')
    parser.add_argument('--pseudo-dir',default='../software/pseudo')
    parser.add_argument('--radius',type=float,default=4.)
    parser.add_argument('--masks-only',action='store_true')
    parser.add_argument('--allow-missing-pseudos',action='store_true')
    args=parser.parse_args()
    output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT.resolve()):raise ValueError('--output must remain inside catalyst_pilot')
    records=prepare(ROOT,output,(ROOT/args.pseudo_dir).resolve(),args.radius,args.masks_only,args.allow_missing_pseudos)
    for r in records:print(f"{r['id']}: {r['movable_count']} mobile / {r['fixed_count']} fixed")


if __name__=='__main__':main()
