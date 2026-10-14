"""Execution receipts and full-cluster scientific gates are independent."""
from pathlib import Path
import math
import re
import numpy as np
from rdkit import Chem
from .common import readj, sha, hash_protocol
from .geometry import with_xyz, validate_geometry
from .inputs import render_input, method_fingerprint, analyte_nuclei, J_STAGES
from .orcaio import parse_output
from .issue_review import user_review_policy, negative_modes, validate_scientific_review


def receipt_errors(directory, receipt, required_inputs, required_outputs):
    directory = Path(directory).resolve()
    errors = []
    for field, required in (('input_sha256', required_inputs), ('output_sha256', required_outputs)):
        hashes = receipt.get(field)
        if not isinstance(hashes, dict) or not hashes:
            errors.append('receipt_hash_section_missing:' + field)
            continue
        for name in required:
            if name not in hashes:
                errors.append('receipt_binding_missing:' + name)
        for name, expected in hashes.items():
            path = (directory / name).resolve()
            if not path.is_relative_to(directory) or not path.is_file() or not isinstance(expected, str) or sha(path).lower() != expected.lower():
                errors.append('receipt_hash_mismatch:' + name)
    return errors


def inspect_job(job_dir, config, require_geometry=True, review_path=None):
    directory = Path(job_dir)
    meta = readj(directory / 'job.json')
    errors = []
    warnings = []
    scientific_review = None
    needs_review = False
    if meta.get('protocol_hash') != hash_protocol(config):
        errors.append('protocol_hash_mismatch')
    output = directory / 'job.out'
    if not output.exists():
        return {'accepted': False, 'errors': ['not_run'], 'job': meta, 'parsed': {}, 'mol': None,
                'warnings': [], 'needs_review': False, 'scientific_review': None,
                'scientific_status': 'INCOMPLETE_NOT_A_CALCULATED_RESULT'}
    receipt_path = directory / 'receipt.json'
    if not receipt_path.is_file():
        errors.append('receipt_missing')
    else:
        receipt = readj(receipt_path)
        if receipt.get('status') != 'completed':
            errors.append('receipt_not_completed')
        if receipt.get('reuse'):
            try:
                from .cache import validate_reuse
                validate_reuse(directory,receipt,config)
            except (ValueError,OSError,KeyError) as exc:
                errors.append('reuse_evidence_invalid:' + str(exc))
        required = ['job.out']
        for name in ('job.xyz', 'job.property.txt'):
            if (directory / name).exists():
                required.append(name)
        errors += receipt_errors(directory, receipt, ['job.inp', 'job.json', 'input.xyz', 'input_structure.mol', 'atom_map.json'], required)
    text = output.read_text(encoding='utf-8', errors='replace')
    property_file = directory / 'job.property.txt'
    parsed = parse_output(text, property_file.read_text(encoding='utf-8', errors='replace') if property_file.exists() else None)
    if not parsed.get('normal_termination'):
        errors.append('not_normally_terminated')
    versions = re.findall(r'Program Version\s+(\d+\.\d+\.\d+)', text)
    if not versions or set(versions) != {config['orca_version']}:
        errors.append('orca_version_mismatch_or_missing')
    for filename, key in (('job.inp', 'input_sha256'), ('atom_map.json', 'atom_map_sha256')):
        if not (directory / filename).is_file() or sha(directory / filename) != meta.get(key):
            errors.append('metadata_hash_mismatch:' + filename)
    template = Chem.MolFromMolBlock((directory / 'input_structure.mol').read_text(encoding='utf-8'), removeHs=False)
    stage = meta['stage']
    try:
        atom_map = readj(directory / 'atom_map.json')
        input_mol = with_xyz(template, (directory / 'input.xyz').read_text(encoding='utf-8'))
        expected_input, expected_pairs = render_input(input_mol, stage, config, meta['solvent'],
                                                    meta['charge'], meta['multiplicity'], atom_map)
        if expected_input != (directory / 'job.inp').read_text(encoding='utf-8'):
            errors.append('input_does_not_match_frozen_protocol_and_coordinates')
        if expected_pairs != meta.get('expected_j_pairs_zero_based', []):
            errors.append('expected_J_pairs_do_not_match_analyte_scope')
        expected_nuclei = analyte_nuclei(atom_map, ('H', 'C') if stage == 'shielding' else ('H', 'C', 'F')) if stage == 'shielding' or stage in J_STAGES else []
        if expected_nuclei != meta.get('expected_nuclei_zero_based', []):
            errors.append('expected_nuclei_do_not_match_analyte_scope')
        actual_br = any(a.GetSymbol() == 'Br' for a in input_mol.GetAtoms())
        if meta.get('method_fingerprint') != method_fingerprint(config, stage, meta['solvent'], has_br=actual_br):
            errors.append('method_fingerprint_mismatch')
    except (ValueError, RuntimeError, KeyError, OSError) as exc:
        errors.append('input_protocol_gate:' + str(exc))
    if stage in ('prescreen', 'optfreq'):
        if not parsed.get('optimization_converged'):
            errors.append('optimization_not_converged')
        frequencies = parsed.get('frequencies_cm1', [])
        if not frequencies:
            errors.append('frequencies_missing')
        if template is None or len(frequencies) != 3 * template.GetNumAtoms():
            errors.append('frequency_count_not_3N')
        if user_review_policy(config):
            # Frequency magnitude is advisory, never an arbitrary approval cutoff.
            if negative_modes(parsed):
                scientific_review = validate_scientific_review(directory, config, parsed, review_path)
                warnings.append('negative_frequencies_observed_requires_documented_mode_review')
                if not scientific_review['accepted']:
                    needs_review = True
                    errors.append('imaginary_frequencies_require_mode_review')
                    errors.extend(scientific_review['errors'])
                else:
                    warnings.append('human_reviewed_approximation_not_a_certified_physical_minimum')
        elif any(f < -float(config['imaginary_frequency_tolerance_cm1']) for f in frequencies):
            errors.append('imaginary_frequencies_require_mode_review')
        g = parsed.get('gibbs_free_energy_Eh')
        if g is None or not math.isfinite(g):
            errors.append('G_missing_or_nonfinite')
        temperature = parsed.get('temperature_K')
        if temperature is None or abs(temperature - config['temperature_K']) > .05:
            errors.append('temperature_mismatch')
        if not parsed.get('qrrho_enabled'):
            errors.append('qrrho_not_confirmed')
        q = parsed.get('qrrho_reference_frequency_cm1')
        if q is None or abs(q - config['qrrho_reference_frequency_cm1']) > .01:
            errors.append('qrrho_reference_mismatch')
    elif stage == 'shielding':
        rows = parsed.get('shielding', [])
        actual = [r['atom_index_zero_based'] for r in rows]
        expected = meta.get('expected_nuclei_zero_based', [])
        if len(set(actual)) != len(actual) or set(actual) != set(expected):
            errors.append('shielding_nucleus_whitelist_mismatch')
    elif stage.startswith('j_'):
        actual = [(r['atom_i_zero_based'], r['atom_j_zero_based']) for r in parsed.get('j_pairs', [])]
        expected = [tuple(pair) for pair in meta.get('expected_j_pairs_zero_based', [])]
        if not expected or len(set(actual)) != len(actual) or set(actual) != set(expected):
            errors.append('J_pair_whitelist_mismatch')
    else:
        errors.append('unknown_stage')
    mol = None
    if require_geometry:
        coordinate = directory / ('job.xyz' if stage in ('prescreen', 'optfreq') else 'input.xyz')
        try:
            mol = with_xyz(template, coordinate.read_text(encoding='utf-8'))
            errors += validate_geometry(template, mol, readj(directory / 'atom_map.json'), config)
        except (OSError, ValueError, RuntimeError) as exc:
            errors.append('geometry:' + str(exc))
    reviewed_approximation = bool(scientific_review and scientific_review['accepted'])
    return {'accepted': not errors, 'errors': sorted(set(errors)), 'job': meta, 'parsed': parsed, 'mol': mol,
            'warnings': sorted(set(warnings)), 'needs_review': needs_review,
            'scientific_review': scientific_review,
            'scientific_status': ('NEEDS_REVIEW_OR_RECALCULATION' if errors else
                                  'HUMAN_REVIEWED_APPROXIMATION' if reviewed_approximation else
                                  'PASSED_CONFIGURED_CHECKS_NOT_EXPERIMENTAL_VALIDATION')}


def boltzmann(gibbs, temperature):
    values = np.asarray(gibbs, dtype=float)
    if not len(values) or not np.isfinite(values).all() or not math.isfinite(temperature) or temperature <= 0:
        raise ValueError('Invalid Gibbs energies/temperature')
    weights = np.exp(-(values - values.min()) * 627.5094740631 / (0.00198720425864083 * temperature))
    return (weights / weights.sum()).tolist()
