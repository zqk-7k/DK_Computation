#!/usr/bin/env python3
"""Prepare, execute and validate fixed-cell Quantum ESPRESSO relaxations.

Candidate construction is not optimization. A relaxed CIF is exported only
after a normal QE ending, electronic and ionic convergence, and force checks.
No electrode potential, oxidation state, or active-site conclusion is assigned.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

import numpy as np
from pymatgen.core import Element, Structure
from pymatgen.io.cif import CifWriter
from build_candidates import (ROOT, clash_report, coordination, load_cif, read_json,
                              safe_child, sha, write_json)

BOHR_A = 0.529177210903
RY_EV = 13.605693122994
NUMBER = r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?'


def value(text):
    return float(str(text).replace('D', 'E').replace('d', 'e'))


def now():
    return datetime.now(timezone.utc).isoformat()


def fingerprint(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def settings(config):
    data = dict(command=['pw.x', '-in', 'qe.in'], pseudo_dir='software/pseudo',
                ecutwfc_Ry=90., ecutrho_Ry=1080., kpoints=[2, 2, 1], nspin=2,
                starting_magnetization={'Co': .5, 'Fe': .6, 'W': 0., 'O': 0., 'C': 0., 'H': 0.},
                hubbard_U_eV={}, conv_thr_Ry=1e-8, etot_conv_thr_Ry=1e-5,
                forc_conv_thr_Ry_bohr=.001, nstep=200, electron_maxstep=250,
                mixing_beta=.2, degauss_Ry=.02, timeout_seconds=172800, total_charge=0.)
    data.update(dipole_correction=True, dipole_vacuum_margin_A=2., dipole_region_width_A=1.5)
    data.update(config.get('dft', {}))
    if int(data['nspin']) != 2:
        raise ValueError('This workflow requires an explicit spin-polarized baseline (nspin=2)')
    for key in ['ecutwfc_Ry', 'ecutrho_Ry', 'conv_thr_Ry', 'etot_conv_thr_Ry', 'forc_conv_thr_Ry_bohr', 'timeout_seconds']:
        if not np.isfinite(float(data[key])) or float(data[key]) <= 0:
            raise ValueError('Invalid positive DFT setting: '+key)
    if float(data['ecutrho_Ry']) < float(data['ecutwfc_Ry']):
        raise ValueError('ecutrho must be >= ecutwfc')
    if len(data['kpoints']) != 3 or any(int(k) != k or k < 1 for k in data['kpoints']):
        raise ValueError('kpoints must be three positive integers')
    if not isinstance(data['command'], list) or not data['command'] or not all(isinstance(x, str) for x in data['command']):
        raise ValueError('dft.command must be a nonempty argv list, not a shell command')
    return data


def dipole_region(structure, dft):
    """Put the sawtooth return region in a nucleus-free vacuum interval.

    This geometric check does not replace inspection of planar-averaged charge
    and potential, especially for charged or electron-spilling slabs.
    """
    if not dft.get('dipole_correction', True):
        return dict(enabled=False)
    matrix = np.asarray(structure.lattice.matrix)
    if (np.max(np.abs(matrix[:2, 2])) > 1e-7 or np.max(np.abs(matrix[2, :2])) > 1e-7
            or matrix[2, 2] <= 0):
        raise ValueError('Automatic dipole placement requires an xy slab and positive z-aligned c vector')
    length = float(matrix[2, 2]); margin = float(dft.get('dipole_vacuum_margin_A', 2.))
    preferred = float(dft.get('dipole_region_width_A', 1.5))
    if margin < 1.5 or preferred <= 0:
        raise ValueError('Dipole vacuum margin must be >=1.5 A and region width positive')
    positions = np.sort(np.mod(structure.frac_coords[:, 2], 1.))
    endpoints = np.r_[0., positions, 1.]
    intervals = sorted(zip(endpoints[:-1], endpoints[1:]), key=lambda x: x[1]-x[0], reverse=True)
    start, stop = intervals[0]
    available = (stop-start)*length-2*margin
    if available < .75:
        raise ValueError('No adequate vacuum region for dipole correction; inspect/expand the slab vacuum before calculation')
    width = min(preferred, available); center = .5*(start+stop)
    result = dict(enabled=True, edir=3, eamp=0., emaxpos=float(center-.5*width/length),
                  eopreg=float(width/length), vacuum_interval_fractional=[float(start), float(stop)],
                  required_nuclear_margin_A=margin, cell_c_A=length,
                  interpretation='Zero applied sawtooth amplitude plus slab dipole correction; not an electrode-potential model')
    if not (0 < result['emaxpos'] < result['emaxpos']+result['eopreg'] < 1):
        raise ValueError('Dipole region leaves the primary cell')
    return result


def dipole_clearance(structure, xyz, region):
    if not region.get('enabled'):
        return dict(passed=True, enabled=False)
    positions = np.mod(structure.lattice.get_fractional_coords(xyz)[:, 2], 1.)
    start, stop = region['emaxpos'], region['emaxpos']+region['eopreg']
    nearest = float('inf')
    for f in positions:
        for shift in [-1., 0., 1.]:
            point = f+shift
            separation = start-point if point < start else point-stop if point > stop else 0.
            nearest = min(nearest, separation*region['cell_c_A'])
    return dict(enabled=True, passed=bool(nearest+1e-6 >= region['required_nuclear_margin_A']),
                minimum_nuclear_clearance_A=float(nearest), required_nuclear_margin_A=region['required_nuclear_margin_A'])


def pseudo_records(root, config, elements, allow_missing=False):
    dft = settings(config)
    definitions = dft.get('pseudopotentials')
    if not definitions:
        meta = root/'docs/SSSP_1.3.0_PBE_efficiency.json'
        if not meta.exists():
            raise ValueError('Supply dft.pseudopotentials or the official SSSP metadata')
        definitions = read_json(meta)
    directory = safe_child(root, dft['pseudo_dir'])
    records = {}
    for element in elements:
        record = definitions.get(element)
        if not isinstance(record, dict):
            raise ValueError('Missing pseudopotential for '+element)
        filename = record.get('file') or record.get('filename')
        if not filename or Path(filename).name != filename:
            raise ValueError('Pseudopotential filename must be a single safe name')
        path = directory/filename
        checksum = record.get('sha256') or record.get('md5')
        algorithm = 'sha256' if record.get('sha256') else 'md5'
        if not checksum:
            raise ValueError('A published pseudopotential digest is required for '+element)
        exists = path.is_file()
        if not exists and not allow_missing:
            raise FileNotFoundError(f'Pseudopotential missing: {path}; run the installer first')
        if exists:
            actual = hashlib.new(algorithm, path.read_bytes()).hexdigest()
            if actual.lower() != checksum.lower():
                raise ValueError('Pseudopotential digest mismatch: '+filename)
        records[element] = dict(file=filename, digest_algorithm=algorithm, published_digest=checksum,
                                sha256=sha(path) if exists else None, installed=exists,
                                suggested_cutoff_Ry=record.get('cutoff_wfc'),
                                suggested_charge_cutoff_Ry=record.get('cutoff_rho'))
    return records


def qe_input(structure, dft, pseudos, pseudo_path):
    elements = list(dict.fromkeys(site.specie.symbol for site in structure))
    dipole = dipole_region(structure, dft)
    lines = ['&CONTROL', " calculation='relax'", " restart_mode='from_scratch'", " prefix='xas_relax'",
             f" pseudo_dir='{pseudo_path.replace(chr(39), chr(39)*2)}'", " outdir='./scratch'",
             " verbosity='high'", ' tprnfor=.true.', ' tstress=.true.',
             f" nstep={int(dft['nstep'])}", f" etot_conv_thr={float(dft['etot_conv_thr_Ry']):.12g}",
             f" forc_conv_thr={float(dft['forc_conv_thr_Ry_bohr']):.12g}", '/', '&SYSTEM',
             ' ibrav=0', f' nat={len(structure)}', f' ntyp={len(elements)}',
             f" ecutwfc={float(dft['ecutwfc_Ry']):.12g}", f" ecutrho={float(dft['ecutrho_Ry']):.12g}",
             ' nspin=2', f" tot_charge={float(dft['total_charge']):.12g}",
             " occupations='smearing'", " smearing='mv'", f" degauss={float(dft['degauss_Ry']):.12g}",
             ' nosym=.true.', ' noinv=.true.']
    if dipole['enabled']:
        control_end = lines.index('/')
        lines[control_end:control_end] = [' tefield=.true.', ' dipfield=.true.']
        lines += [' edir=3', ' eamp=0.0', f" emaxpos={dipole['emaxpos']:.12g}", f" eopreg={dipole['eopreg']:.12g}"]
    for i, element in enumerate(elements, 1):
        lines.append(f" starting_magnetization({i})={float(dft['starting_magnetization'].get(element, 0.)):.8g}")
    lines += ['/', '&ELECTRONS', f" conv_thr={float(dft['conv_thr_Ry']):.12g}",
              f" electron_maxstep={int(dft['electron_maxstep'])}", f" mixing_beta={float(dft['mixing_beta']):.12g}",
              " diagonalization='david'", '/', '&IONS', " ion_dynamics='bfgs'", '/', 'ATOMIC_SPECIES']
    for element in elements:
        lines.append(f'{element} {float(Element(element).atomic_mass):.10g} {pseudos[element]["file"]}')
    lines.append('CELL_PARAMETERS angstrom')
    lines += [' '.join(f'{x:.12f}' for x in vector) for vector in structure.lattice.matrix]
    lines.append('ATOMIC_POSITIONS angstrom')
    lines += [f'{site.specie.symbol} ' + ' '.join(f'{x:.12f}' for x in site.coords) for site in structure]
    lines += ['K_POINTS automatic', ' '.join(str(int(k)) for k in dft['kpoints'])+' 0 0 0']
    positive_u = {key: float(v) for key, v in dft.get('hubbard_U_eV', {}).items() if float(v) > 0 and key in elements}
    if positive_u:
        orbitals = {'Co': '3d', 'Fe': '3d', 'W': '5d'}
        lines.append('HUBBARD (ortho-atomic)')
        for element, u in positive_u.items():
            if element not in orbitals:
                raise ValueError('Explicit Hubbard orbital mapping needed for '+element)
            lines.append(f'U {element}-{orbitals[element]} {u:.8g}')
    return '\n'.join(lines)+'\n'


def prepare(config_path, allow_missing=False):
    path = Path(config_path).resolve(); root = path.parent; config = read_json(path); dft = settings(config)
    candidates = read_json(root/'candidates/manifest.json')
    jobs = []
    for item in candidates:
        source = safe_child(root, item['path'])/'structure.cif'
        if sha(source) != item['structure_sha256']:
            raise ValueError('Candidate structure changed: '+item['id'])
        structure, _ = load_cif(source)
        elements = list(dict.fromkeys(site.specie.symbol for site in structure))
        pseudos = pseudo_records(root, config, elements, allow_missing)
        directory = safe_child(root, 'dft/'+item['id']); directory.mkdir(parents=True, exist_ok=True)
        pseudo_path = os.path.relpath(safe_child(root, dft['pseudo_dir']), directory).replace(os.sep, '/')
        text = qe_input(structure, dft, pseudos, pseudo_path)
        inp = directory/'qe.in'; metadata_path = directory/'job.json'
        contract = dict(dft=dft, pseudos=pseudos, structure_sha256=sha(source))
        contract_sha = fingerprint(contract)
        if inp.exists() and inp.read_text(encoding='utf-8') != text:
            raise ValueError('Refusing to overwrite changed QE input: '+str(inp)+'; use a separate project/profile')
        if metadata_path.exists() and read_json(metadata_path).get('settings_sha256') != contract_sha:
            raise ValueError('QE provenance/settings changed: '+str(metadata_path)+'; preserve this job and create a new run directory')
        inp.write_text(text, encoding='utf-8', newline='\n')
        job = dict(id=item['id'], path=directory.relative_to(root).as_posix(),
                   candidate_path=item['path'], structure_sha256=sha(source), input_sha256=sha(inp),
                   settings_sha256=contract_sha, pseudopotentials=pseudos,
                   dipole_region=dipole_region(structure, dft),
                   settings=dft, status='prepared_not_calculated',
                   all_pseudopotentials_installed=all(p['installed'] for p in pseudos.values()))
        write_json(metadata_path, job); jobs.append(job)
    write_json(root/'dft/jobs.json', jobs)
    print(f'Prepared {len(jobs)} QE relaxation inputs; none is marked optimized.')
    return jobs


def final_positions(text, expected_symbols, lattice):
    lines = text.splitlines(); blocks = []
    alat_matches = re.findall(r'lattice parameter\s*\(alat\)\s*=\s*('+NUMBER+r')', text, re.I)
    alat = value(alat_matches[-1])*BOHR_A if alat_matches else None
    for i, line in enumerate(lines):
        if not re.match(r'\s*ATOMIC_POSITIONS\b', line):
            continue
        unit_name = re.sub(r'[(){}]', ' ', line).split()[1:]
        unit_name = unit_name[0].lower() if unit_name else 'alat'
        positions, symbols = [], []
        for atom_line in lines[i+1:i+1+len(expected_symbols)]:
            parts = atom_line.split()
            if len(parts) < 4 or not re.fullmatch(r'[A-Z][a-z]?', parts[0]):
                break
            try:
                positions.append([value(x) for x in parts[1:4]]); symbols.append(parts[0])
            except ValueError:
                break
        if len(positions) == len(expected_symbols):
            blocks.append((unit_name, symbols, np.asarray(positions)))
    if not blocks:
        raise ValueError('No complete final ATOMIC_POSITIONS block')
    unit_name, symbols, xyz = blocks[-1]
    if symbols != expected_symbols:
        raise ValueError('Final atom identities/order differ from the input')
    if unit_name == 'crystal':
        xyz = lattice.get_cartesian_coords(xyz)
    elif unit_name == 'bohr':
        xyz *= BOHR_A
    elif unit_name == 'alat' and alat is not None:
        xyz *= alat
    elif unit_name != 'angstrom':
        raise ValueError('Unknown ATOMIC_POSITIONS unit: '+unit_name)
    if not np.isfinite(xyz).all():
        raise ValueError('Non-finite final coordinates')
    return xyz


def parse_output(text, structure, dft):
    failures = []
    if 'JOB DONE.' not in text:
        failures.append('QE normal completion marker missing')
    electronic = re.findall(r'convergence (?:has been achieved|NOT achieved)[^\n]*', text, re.I)
    if not electronic or 'has been achieved' not in electronic[-1].lower():
        failures.append('Final SCF convergence not confirmed')
    ionic = bool(re.search(r'bfgs converged in', text, re.I)) and 'End of BFGS Geometry Optimization' in text
    if not ionic:
        failures.append('Ionic BFGS convergence not confirmed')
    energies = re.findall(r'!\s*total energy\s*=\s*('+NUMBER+r')\s*Ry', text, re.I)
    if not energies:
        failures.append('Final total energy missing')
    starts = [match.end() for match in re.finditer(r'Forces acting on atoms', text)]
    forces = []
    if starts:
        segment = text[starts[-1]:]
        for match in re.finditer(r'atom\s+\d+\s+type\s+\d+\s+force\s*=\s*('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')', segment):
            forces.append([value(match.group(j)) for j in [1, 2, 3]])
    force_max = float(np.max(np.linalg.norm(forces, axis=1))) if len(forces) == len(structure) else None
    if force_max is None:
        failures.append('Complete final atomic-force block missing')
    elif force_max > float(dft['forc_conv_thr_Ry_bohr'])*1.01:
        failures.append('Final maximum atomic force exceeds configured threshold')
    try:
        xyz = final_positions(text, [s.specie.symbol for s in structure], structure.lattice)
    except ValueError as exc:
        failures.append(str(exc)); xyz = None
    if re.search(r'Error in routine|MPI_ABORT|SIGSEGV', text, re.I):
        failures.append('QE/MPI error marker present')
    metrics = dict(normal_completion='JOB DONE.' in text, electronic_convergence=bool(electronic) and 'has been achieved' in electronic[-1].lower(),
                   ionic_convergence=ionic, final_energy_Ry=value(energies[-1]) if energies else None,
                   final_energy_eV=value(energies[-1])*RY_EV if energies else None,
                   max_atomic_force_Ry_bohr=force_max,
                   max_atomic_force_eV_A=force_max*RY_EV/BOHR_A if force_max is not None else None,
                   force_threshold_Ry_bohr=float(dft['forc_conv_thr_Ry_bohr']), failures=failures)
    return xyz, metrics


def collect_one(root, job, extra=None):
    directory = safe_child(root, job['path']); source = safe_child(root, job['candidate_path'])/'structure.cif'
    if sha(directory/'qe.in') != job['input_sha256'] or sha(source) != job['structure_sha256']:
        raise ValueError('DFT input/candidate changed after preparation')
    output = directory/'qe.out'
    base = dict(id=job['id'], status='failed', collected_utc=now(), input_sha256=job['input_sha256'],
                structure_sha256=job['structure_sha256'], settings_sha256=job['settings_sha256'],
                settings=job['settings'], pseudopotentials=job['pseudopotentials'])
    if extra:
        base.update(extra)
    if not output.is_file():
        base['failures'] = ['qe.out missing']; write_json(directory/'result.json', base); return base
    structure, _ = load_cif(source)
    xyz, metrics = parse_output(output.read_text(encoding='utf-8', errors='replace'), structure, job['settings'])
    base.update(metrics=metrics, output_sha256=sha(output))
    if extra and extra.get('returncode') not in {None, 0}:
        metrics['failures'].append('QE process returned a nonzero exit status')
    if metrics['failures']:
        write_json(directory/'result.json', base); return base
    metadata = read_json(safe_child(root, job['candidate_path'])/'metadata.json')
    atoms = [dict(a, xyz=list(map(float, pos))) for a, pos in zip(metadata['atom_map'], xyz)]
    dipole = job.get('dipole_region', dipole_region(structure, job['settings']))
    dipole_check = dipole_clearance(structure, xyz, dipole)
    base.update(dipole_region=dipole, dipole_clearance_check=dipole_check)
    if not dipole_check['passed']:
        metrics['failures'].append('Relaxed atoms entered the dipole correction vacuum safety region')
        write_json(directory/'result.json', base); return base
    clashes = clash_report(structure.lattice, atoms)
    if not clashes['passed']:
        metrics['failures'].append('Relaxed geometry has severe interatomic overlaps')
        base['clash_check'] = clashes; write_json(directory/'result.json', base); return base
    if [a['element'] for a in metadata['atom_map']] != [s.specie.symbol for s in structure]:
        raise ValueError('Candidate atom map order does not match QE input atom order')
    relaxed = Structure(structure.lattice, [s.specie.symbol for s in structure], xyz,
                        coords_are_cartesian=True, labels=[a['id'] for a in metadata['atom_map']])
    relaxed_path = directory/'relaxed.cif'; relaxed_path.write_text(str(CifWriter(relaxed, symprec=None)), encoding='utf-8', newline='\n')
    coord = coordination(structure.lattice, atoms, metadata['substitution']['atom_id'])
    base.update(status='optimized', relaxed_cif=relaxed_path.relative_to(root).as_posix(), relaxed_sha256=sha(relaxed_path),
                composition=dict(Counter(s.specie.symbol for s in relaxed)), relaxed_coordination=coord, clash_check=clashes,
                oxygen_coordination_hypothesis_retained=coord['all_heavy_neighbours_are_O'],
                warnings=[] if coord['all_heavy_neighbours_are_O'] else ['Relaxation changed W first-shell chemistry; retain and report this unfavorable outcome. Do not alter the geometry to pass.'],
                physical_interpretation='PBE spin-polarized fixed-cell candidate relaxation; no applied-potential assignment and no inherited VASP electronic state')
    write_json(directory/'result.json', base)
    return base


def run(config_path, task):
    config_path = Path(config_path).resolve(); root = config_path.parent
    jobs = read_json(root/'dft/jobs.json')
    if task < 0 or task >= len(jobs):
        raise ValueError('DFT task index out of range')
    job = jobs[task]; directory = safe_child(root, job['path'])
    config = read_json(config_path)
    if settings(config) != job['settings']:
        raise ValueError('DFT settings changed since preparation; create a separate run')
    pseudo_records(root, config, job['pseudopotentials'].keys())
    existing = directory/'result.json'
    if existing.exists():
        result = read_json(existing)
        if (result.get('status') == 'optimized' and result.get('input_sha256') == sha(directory/'qe.in')
                and result.get('settings_sha256') == job['settings_sha256']
                and (directory/'qe.out').is_file() and result.get('output_sha256') == sha(directory/'qe.out')
                and (directory/'relaxed.cif').is_file() and result.get('relaxed_sha256') == sha(directory/'relaxed.cif')):
            print('Already optimized:', job['id']); return result
    lock = directory/'.run.lock'
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError as exc:
        raise RuntimeError(f'Job lock exists: {lock}; inspect scheduler/process before clearing a stale lock') from exc
    with os.fdopen(descriptor, 'w', encoding='utf-8') as stream:
        stream.write(json.dumps(dict(pid=os.getpid(), started_utc=now())))
    started = time.monotonic(); command = []
    try:
        if sha(directory/'qe.in') != job['input_sha256']:
            raise ValueError('QE input changed after preparation')
        for name in ['qe.out', 'qe.err']:
            old = directory/name
            if old.exists():
                tag = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
                old.rename(directory/(name+'.previous_'+tag))
        command = [arg.format(root=str(root), job_dir=str(directory), input='qe.in') for arg in job['settings']['command']]
        with (directory/'qe.out').open('w', encoding='utf-8') as out, (directory/'qe.err').open('w', encoding='utf-8') as err:
            process = subprocess.run(command, cwd=directory, stdout=out, stderr=err,
                                     timeout=float(job['settings']['timeout_seconds']), check=False)
        result = collect_one(root, job, dict(command=command, returncode=process.returncode,
                                             duration_seconds=time.monotonic()-started, execution_started_by_runner=True))
        if result['status'] != 'optimized':
            raise RuntimeError('QE job failed validation: '+job['id'])
        return result
    except (OSError, subprocess.TimeoutExpired) as exc:
        write_json(directory/'result.json', dict(id=job['id'], status='failed', input_sha256=job['input_sha256'],
                   settings_sha256=job['settings_sha256'], duration_seconds=time.monotonic()-started,
                   command=command, error=str(exc)))
        raise
    finally:
        lock.unlink(missing_ok=True)


def collect(config_path):
    root = Path(config_path).resolve().parent; jobs = read_json(root/'dft/jobs.json'); results = []
    for job in jobs:
        existing = safe_child(root, job['path'])/'result.json'
        previous = read_json(existing) if existing.exists() else {}
        extra = {key: previous[key] for key in ['command', 'returncode', 'duration_seconds', 'execution_started_by_runner'] if key in previous}
        results.append(collect_one(root, job, extra))
    report = dict(optimized=sum(r['status'] == 'optimized' for r in results), expected=len(jobs), results=results,
                  all_complete=all(r['status'] == 'optimized' for r in results))
    write_json(root/'dft/collection_report.json', report)
    if not report['all_complete']:
        raise RuntimeError(f"Only {report['optimized']}/{len(jobs)} DFT jobs are optimized; inspect result.json and raw outputs")
    print(f"Validated {len(jobs)} optimized structures. Inspect oxygen coordination and model assumptions before interpretation.")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'run', 'collect'])
    parser.add_argument('--config', default=str(ROOT/'config.json'))
    parser.add_argument('--task', type=int)
    parser.add_argument('--allow-missing-pseudos', action='store_true', help='Input inspection only, not execution')
    args = parser.parse_args()
    if args.action == 'prepare':
        prepare(args.config, args.allow_missing_pseudos)
    elif args.action == 'run':
        if args.task is None:
            parser.error('run requires --task INDEX')
        run(args.config, args.task)
    else:
        collect(args.config)


if __name__ == '__main__':
    main()
