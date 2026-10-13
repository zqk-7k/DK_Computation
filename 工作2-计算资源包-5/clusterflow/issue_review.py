"""Advisory failure diagnosis and source-bound human review, never auto-repair.

An exception may acknowledge negative modes in a COMPLETE opt/frequency job.
It does not turn a failed or incomplete calculation into a valid result and does
not certify a stationary structure as a physical minimum.
"""
from datetime import datetime
from pathlib import Path
import math
import re

from .common import readj, sha, hash_protocol
from .orcaio import parse_output


REVIEW_CHECKS = (
    'all_negative_modes_inspected',
    'solute_connectivity_and_stereochemistry_checked',
    'solvent_contacts_and_proton_transfer_checked',
    'reported_thermochemistry_and_weighting_limitations_acknowledged',
    'alternative_reoptimization_or_displacement_considered',
)
SOURCE_NAMES = ('job.json', 'job.inp', 'input.xyz', 'input_structure.mol',
                'atom_map.json', 'receipt.json', 'job.out', 'job.xyz',
                'job.property.txt', 'job.hess')


def user_review_policy(config):
    value = config.get('scientific_issue_policy', 'strict')
    return (value.get('mode') if isinstance(value, dict) else value) == 'user_review'


def negative_modes(parsed):
    return [{'mode_index_zero_based': i, 'frequency_cm1': float(value)}
            for i, value in enumerate(parsed.get('frequencies_cm1', []))
            if math.isfinite(value) and value < 0]


def source_hashes(job_dir):
    directory = Path(job_dir)
    return {name: sha(directory / name) for name in SOURCE_NAMES
            if (directory / name).is_file()}


def _read_outputs(job_dir):
    directory = Path(job_dir)
    out = directory / 'job.out'
    prop = directory / 'job.property.txt'
    text = out.read_text(encoding='utf-8', errors='replace') if out.is_file() else ''
    parsed = parse_output(text, prop.read_text(encoding='utf-8', errors='replace')
                          if prop.is_file() else None)
    return text, parsed


def diagnose_job(job_dir, config):
    """Return recommendations even for missing/truncated outputs; execute nothing."""
    directory = Path(job_dir)
    meta = readj(directory / 'job.json') if (directory / 'job.json').is_file() else {}
    receipt = readj(directory / 'receipt.json') if (directory / 'receipt.json').is_file() else {}
    text, parsed = _read_outputs(directory)
    findings = []
    if not (directory / 'job.out').is_file():
        findings.append({'code': 'not_run_or_output_missing', 'reviewable_exception': False,
                         'suggestion': 'Locate the actual output or explicitly run the job; do not infer results from the initial geometry.'})
    else:
        if receipt.get('timed_out') or receipt.get('failure_reason') == 'timeout':
            findings.append({'code': 'timeout', 'reviewable_exception': False,
                             'suggestion': 'Inspect the last completed iteration and checkpoint. Decide whether to resume/restart in a new attempt; preserve the original output. A timeout is not a complete result.'})
        if not parsed.get('normal_termination'):
            findings.append({'code': 'not_normally_terminated', 'reviewable_exception': False,
                             'suggestion': 'Read the final error and resource/checkpoint information. Correct the identified cause, then obtain a normally terminated actual output.'})
        if meta.get('stage') in ('prescreen', 'optfreq') and not parsed.get('optimization_converged'):
            findings.append({'code': 'optimization_not_converged', 'reviewable_exception': False,
                             'suggestion': 'Inspect energy/gradient history and geometry; consider a better initial structure, smaller steps, or an updated Hessian. Choose settings after review and rerun; do not relabel the old geometry as converged.'})
        if re.search(r'(SCF\s+(?:NOT\s+CONVERGED|CONVERGENCE\s+FAIL)|SCF.*DID\s+NOT\s+CONVERGE)', text, re.I):
            findings.append({'code': 'scf_difficulty', 'reviewable_exception': False,
                             'suggestion': 'Review charge, multiplicity, contacts, and the SCF trace. Consider an appropriate initial guess or convergence aids; verify a new actual result.'})
        modes = negative_modes(parsed)
        if modes:
            findings.append({'code': 'negative_frequencies', 'reviewable_exception': True,
                             'modes': modes,
                             'suggestion': 'Animate every negative mode and distinguish numerical/solvent motions from real solute or reaction coordinates. Consider displacements in both directions and reoptimization. No frequency-magnitude cutoff automatically accepts or rejects a complete job; record a source-bound human decision.'})
    return {'status': 'ADVISORY_ONLY_NOT_SCIENTIFIC_APPROVAL', 'job_dir': str(directory),
            'protocol_hash': hash_protocol(config), 'job': meta, 'source_sha256': source_hashes(directory),
            'findings': findings, 'negative_modes': negative_modes(parsed),
            'normal_termination_observed': bool(parsed.get('normal_termination')),
            'optimization_convergence_observed': bool(parsed.get('optimization_converged')),
            'automatic_retry': False, 'automatic_deletion': False,
            'physical_minimum_certified': False}


def draft_scientific_review(job_dir, config):
    """Create an unapproved review object. The caller may save it for the user."""
    diagnosis = diagnose_job(job_dir, config)
    return {'schema_version': '1.0', 'status': 'DRAFT_REQUIRES_HUMAN_REVIEW',
            'scope': 'negative_frequency_exception_only',
            'protocol_hash': hash_protocol(config),
            'stage': diagnosis['job'].get('stage'),
            'source_sha256': diagnosis['source_sha256'],
            'negative_modes': diagnosis['negative_modes'],
            'reviewer': '', 'reviewed_utc': '', 'rationale': '',
            'limitations': '', 'checks': {key: False for key in REVIEW_CHECKS},
            'mode_review_evidence': [], 'physical_minimum_certified': False,
            'note': 'No fields in this draft constitute approval. This cannot override failed/incomplete jobs, missing data, geometry failures, or source tampering.'}


def validate_scientific_review(job_dir, config, parsed, review_path=None):
    """Validate a manual negative-mode exception; all other gates stay separate."""
    directory = Path(job_dir).resolve()
    path = Path(review_path) if review_path is not None else directory / 'scientific_review.json'
    errors = []
    if not path.is_file():
        return {'accepted': False, 'errors': ['scientific_review_missing'], 'review': None}
    try:
        review = readj(path)
        if review.get('status') != 'ACCEPT_AS_REVIEWED_APPROXIMATION':
            errors.append('scientific_review_not_approved')
        if review.get('scope') != 'negative_frequency_exception_only':
            errors.append('scientific_review_scope_invalid')
        if review.get('protocol_hash') != hash_protocol(config):
            errors.append('scientific_review_protocol_mismatch')
        meta = readj(directory / 'job.json')
        if review.get('stage') != meta.get('stage') or meta.get('stage') not in ('prescreen', 'optfreq'):
            errors.append('scientific_review_stage_mismatch')
        if review.get('source_sha256') != source_hashes(directory):
            errors.append('scientific_review_source_hash_mismatch')
        if review.get('negative_modes') != negative_modes(parsed) or not negative_modes(parsed):
            errors.append('scientific_review_mode_list_mismatch_or_empty')
        reviewer = review.get('reviewer', '').strip()
        if not reviewer or reviewer.lower() in {'ai', 'assistant', 'codex', 'claude', 'chatgpt', 'system', 'automatic', 'auto', 'software'}:
            errors.append('scientific_review_requires_named_human')
        try:
            reviewed = datetime.fromisoformat(review.get('reviewed_utc', '').replace('Z', '+00:00'))
            if reviewed.tzinfo is None:
                raise ValueError('timezone missing')
        except (TypeError, ValueError):
            errors.append('scientific_review_timestamp_invalid')
        for field in ('rationale', 'limitations'):
            if not isinstance(review.get(field), str) or not review[field].strip():
                errors.append('scientific_review_' + field + '_missing')
        if any(review.get('checks', {}).get(key) is not True for key in REVIEW_CHECKS):
            errors.append('scientific_review_checks_incomplete')
        if review.get('physical_minimum_certified') is not False:
            errors.append('scientific_review_cannot_certify_physical_minimum')
        evidence = review.get('mode_review_evidence', [])
        if not isinstance(evidence, list) or not evidence:
            errors.append('scientific_review_mode_evidence_missing')
        else:
            for row in evidence:
                candidate = (directory / row.get('path', '')).resolve()
                if (not candidate.is_relative_to(directory) or not candidate.is_file()
                        or not row.get('description', '').strip()
                        or sha(candidate) != row.get('sha256')):
                    errors.append('scientific_review_mode_evidence_invalid')
    except (OSError, ValueError, TypeError, AttributeError) as exc:
        return {'accepted': False, 'errors': ['scientific_review_invalid:' + str(exc)], 'review': None}
    return {'accepted': not errors, 'errors': sorted(set(errors)), 'review': review,
            'path': str(path), 'sha256': sha(path), 'physical_minimum_certified': False}
