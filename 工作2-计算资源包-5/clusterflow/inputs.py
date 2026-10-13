"""Cluster-preserving ORCA inputs. Only analyte H/C/F nuclei are requested."""
import hashlib
import json
import math
import numpy as np
from .geometry import isotope_inline_coordinates, validate_atom_map

STAGES = ('prescreen', 'optfreq', 'shielding', 'j_PBE_pcJ3', 'j_PBE0_pcJ2')
J_STAGES = ('j_PBE_pcJ3', 'j_PBE0_pcJ2')
ALLOWED_SOLVENTS = {'DMSO', 'Methanol', 'Water', 'Chloroform', 'Acetonitrile', 'Acetone'}


def analyte_nuclei(atom_map, elements=('H', 'C')):
    return [r['index_zero_based'] for r in atom_map
            if (r.get('include_in_analyte_nmr') is True or r['element'] == 'F')
            and r.get('component') == 'analyte' and r.get('environment_atom') is False
            and r['element'] in elements and not (r['element'] == 'H' and r.get('isotope') == 2)]


def selected_j_pairs(mol, atom_map, cutoff=1000.0):
    # The whitelist is topology/scope based, not conformer-distance based. All
    # eligible pairs are retained so missing long-range values never mean zero.
    if not math.isfinite(float(cutoff)) or float(cutoff) < 1000.0:
        raise ValueError('This fixed-pair workflow requires distance_cutoff_A >= 1000')
    indices = analyte_nuclei(atom_map, ('H', 'C', 'F'))
    pairs = []
    for offset, i in enumerate(indices):
        for j in indices[:offset]:
            if sorted((atom_map[i]['element'], atom_map[j]['element'])) in [['H', 'H'], ['F', 'H'], ['C', 'F']]:
                pairs.append([min(i, j), max(i, j)])
    return sorted(pairs)


def method_fingerprint(config, stage, solvent, has_br=False):
    settings = {k: v for k, v in config[stage].items() if k not in ('nprocs', 'maxcore_mb')}
    payload = {'orca_version': config['orca_version'], 'stage': stage, 'settings': settings,
               'solvent': solvent, 'gauge': 'GIAO' if stage == 'shielding' else None,
               'nuclei_scope': 'ANALYTE_ONLY_EXPLICIT_CLUSTER_RETAINED',
               'actual_Br_mixed_J_basis': config[stage].get('bromine_basis') if has_br and stage in J_STAGES else None}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def render_input(mol, stage, config, solvent, charge, multiplicity=1, atom_map=None):
    if stage not in STAGES or solvent not in ALLOWED_SOLVENTS:
        raise ValueError('Unsupported stage or solvent')
    if atom_map is None:
        raise ValueError('Explicit immutable cluster atom_map required')
    issues = validate_atom_map(mol, atom_map)
    if issues:
        raise ValueError('Invalid cluster atom map: ' + str(issues))
    if sum(a.GetFormalCharge() for a in mol.GetAtoms()) != int(charge):
        raise ValueError('Cluster charge does not match fixed topology')
    if int(multiplicity) < 1:
        raise ValueError('Invalid multiplicity')
    settings = config[stage]
    keywords = settings['keywords']
    if any(word in keywords.lower() for word in ('revdsd', 'new_job', 'xyzfile')) or '\n' in keywords:
        raise ValueError('Unsafe or prohibited method keywords')
    lines = [f'! {keywords} CPCM({solvent}) MASS2016']
    if stage in ('prescreen', 'optfreq'):
        lines[0] += ' Freq'
    if stage == 'shielding':
        lines[0] += ' NMR'
    ranks, memory = int(settings['nprocs']), int(settings['maxcore_mb'])
    if ranks < 1 or memory < 1:
        raise ValueError('Invalid resource settings')
    lines += [f'%maxcore {memory}', f'%pal nprocs {ranks} end']
    if stage in ('prescreen', 'optfreq'):
        lines += ['%freq', '  QuasiRRHO true',
                  f"  QRRHORefFreq {float(config['qrrho_reference_frequency_cm1'])}",
                  f"  CutOffFreq {float(config['qrrho_cutoff_frequency_cm1'])}",
                  f"  Temp {float(config['temperature_K'])}", '  Pressure 1.0', 'end']
    if stage in J_STAGES and any(a.GetSymbol() == 'Br' for a in mol.GetAtoms()):
        if settings.get('bromine_basis') != 'def2-TZVPP':
            raise ValueError('Br requires explicitly declared def2-TZVPP mixed J basis')
        lines += ['# Br mixed-basis J model: domain validation still required', '%basis',
                  '  NewGTO Br "def2-TZVPP" end', 'end']
    lines += [f'* xyz {int(charge)} {int(multiplicity)}', isotope_inline_coordinates(mol, atom_map), '*']
    pairs = []
    if stage == 'shielding' or stage in J_STAGES:
        indices = analyte_nuclei(atom_map, ('H', 'C') if stage == 'shielding' else ('H', 'C', 'F'))
        if not indices:
            raise ValueError('No eligible analyte nuclei')
        lines += ['%eprnmr']
        if stage == 'shielding':
            lines += ['  Ori = GIAO']
        for element, isotope in (('H', 1), ('C', 13), ('F', 19)):
            selected = [i + 1 for i in indices if atom_map[i]['element'] == element]
            if selected:
                request = 'shift' if stage == 'shielding' else 'ssall'
                lines += [f"  Nuclei = {','.join(map(str, selected))} {{ {request}, ist = {isotope} }}"]
        if stage in J_STAGES:
            cutoff = float(settings.get('distance_cutoff_A', 1000.0))
            pairs = selected_j_pairs(mol, atom_map, cutoff)
            if not pairs:
                raise ValueError('No eligible J pairs; record explicit no-J policy instead')
            lines += ['  SpinSpinAtomPairs ' + ' '.join('{' + str(i) + ' ' + str(j) + '}' for i, j in pairs),
                      f'  SpinSpinRThresh {cutoff}']
        lines += ['end']
    return '\n'.join(lines) + '\n', pairs
