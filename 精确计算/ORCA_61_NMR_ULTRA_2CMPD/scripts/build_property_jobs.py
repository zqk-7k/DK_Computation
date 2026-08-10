"""Build analyte-only NMR shielding/J inputs from APPROVED_FOR_NMR gates.

The complete explicit cluster remains in the electronic model.  Only nuclei
marked ``ANALYTE`` in the hash-locked per-conformer map are requested.  The
builder is fail-closed for unapproved reference scales, J benchmarks and the
F3_006 iodine/X2C protocol.  It never launches ORCA.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

from draft_minimum_boltzmann_gate import normalized_weights, parse_orca_qrrho
from workflow_schema import (
    NOTICE,
    WorkflowValidationError,
    atomic_write_bytes,
    atomic_write_csv,
    ensure_within,
    normalize_sha256,
    orca_input_bundle_sha256,
    parse_xyz_frames,
    read_csv,
    read_json,
    read_text_auto,
    sha256_bytes,
    sha256_file,
    validate_completed_opt_freq_state,
)


PROJECT = Path(__file__).resolve().parents[1]
MINIMUM_APPROVED = "APPROVED_FOR_NMR"
PROTOCOL_APPROVED = "APPROVED_FOR_NMR_PROPERTIES"
PASS = "PASS"
GENERIC_J_FINGERPRINT = (
    "PBE|pcJ-3|AutoAux|CPCM(actual_solvent)|VeryTightSCF|DefGrid3|"
    "EPRNMR_Tol=1e-10|Z_Tol=1e-12|SpinSpinRThresh=6.05A|ORCA_6.1.1"
)
ISOTOPE = {"H": 1, "C": 13, "F": 19, "P": 31}
NUCLEUS = {"H": "1H", "C": "13C", "F": "19F", "P": "31P"}
ALLOWED_J_CLASSES = {
    tuple(sorted(pair))
    for pair in (("H", "H"), ("H", "F"), ("F", "F"), ("C", "F"), ("H", "P"), ("P", "P"), ("C", "P"))
}


@dataclass
class Plan:
    directory: Path
    files: dict[str, bytes]
    row: dict[str, object]


def safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    if not token or token in {".", ".."}:
        raise WorkflowValidationError(f"Unsafe identifier: {value!r}")
    return token


def continuum_token(value: str) -> str:
    token = str(value).strip().upper()
    aliases = {"DMSO": "DMSO", "DMSO-D6": "DMSO", "WATER": "WATER", "D2O": "WATER"}
    if token not in aliases:
        raise WorkflowValidationError(f"Unsupported CPCM domain: {value!r}")
    return aliases[token]


def solvent_domain(primary: str, continuum: str) -> str:
    primary_upper = str(primary).upper()
    return "D2O" if "D2O" in primary_upper or continuum_token(continuum) == "WATER" else "DMSO-d6"


def shielding_branches(continuum: str) -> list[dict[str, str]]:
    solvent = continuum_token(continuum)
    return [
        {
            "branch": "production_revdsd_pcsseg3",
            "fingerprint": (
                f"revDSD-PBEP86-D4/2021|pcSseg-3|AutoAux|RIJCOSX|GIAO|NoFrozenCore|"
                f"CPCM({solvent})|ExtremeSCF|DefGrid3|ORCA_6.1.1"
            ),
            "keywords": (
                "revDSD-PBEP86-D4/2021 pcSseg-3 AutoAux NMR NoFrozenCore "
                "CPCM({CONTINUUM}) RIJCOSX DefGrid3 ExtremeSCF"
            ),
            "basis": "",
            "extra_blocks": "%mp2\n  Density Relaxed\nend",
        },
        {
            "branch": "crosscheck_pbe0_pcsseg3",
            "fingerprint": (
                f"PBE0|pcSseg-3|AutoAux|RIJCOSX|GIAO|CPCM({solvent})|"
                "ExtremeSCF|DefGrid3|ORCA_6.1.1"
            ),
            "keywords": "PBE0 pcSseg-3 AutoAux NMR CPCM({CONTINUUM}) RIJCOSX DefGrid3 ExtremeSCF",
            "basis": "",
            "extra_blocks": "",
        },
    ]


def shielding_fingerprint(continuum: str) -> str:
    """Backward-compatible name for the locked production branch fingerprint."""
    return shielding_branches(continuum)[0]["fingerprint"]


def _csv_bytes(rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(fields), extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return ("\ufeff" + stream.getvalue()).encode("utf-8")


def analyte_index_map(path: Path) -> dict[str, list[int]]:
    result: dict[str, list[int]] = defaultdict(list)
    for row in read_csv(path):
        if row.get("nmr_scope", "") != "ANALYTE":
            continue
        element = row.get("element", "").strip().upper()
        if element in ISOTOPE:
            result[element].append(int(row["atom_index_1based"]))
    return {element: sorted(indices) for element, indices in result.items()}


def j_pair_allowlist(
    index_map: Mapping[str, Sequence[int]], optimized_xyz: Path, maximum_distance_A: float = 6.0
) -> list[dict[str, object]]:
    frame = parse_xyz_frames(optimized_xyz)
    if len(frame) != 1:
        raise WorkflowValidationError(f"Expected one optimized XYZ frame: {optimized_xyz}")
    xyz = frame[0]
    selected = sorted((index, element) for element, values in index_map.items() for index in values)
    rows: list[dict[str, object]] = []
    for left_position, (left, left_element) in enumerate(selected):
        for right, right_element in selected[left_position + 1 :]:
            pair_class = tuple(sorted((left_element, right_element)))
            if pair_class not in ALLOWED_J_CLASSES:
                continue
            a = xyz.coordinates[left - 1]
            b = xyz.coordinates[right - 1]
            distance = math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b, strict=True)))
            if distance > maximum_distance_A:
                continue
            rows.append(
                {
                    "left_atom_index_1based_nuclei": left,
                    "right_atom_index_1based_nuclei": right,
                    "left_atom_index_0based_spinspin": left - 1,
                    "right_atom_index_0based_spinspin": right - 1,
                    "left_element": left_element,
                    "right_element": right_element,
                    "pair_class": f"{left_element}-{right_element}",
                    "distance_A": f"{distance:.6f}",
                    "scope": "ANALYTE_ONLY",
                }
            )
    return rows


def render_shielding_input(
    charge: int,
    multiplicity: int,
    continuum: str,
    index_map: Mapping[str, Sequence[int]],
    nprocs: int,
    maxcore_mb: int,
    simple_keywords: str | None = None,
    basis_block: str = "",
    extra_blocks: str = "",
) -> str:
    solvent = continuum_token(continuum)
    keywords = (
        simple_keywords.replace("{CONTINUUM}", solvent)
        if simple_keywords
        else f"PBE0 pcSseg-2 AutoAux RIJCOSX TightSCF DefGrid3 CPCM({solvent})"
    )
    if f"CPCM({solvent})" not in keywords.upper():
        raise WorkflowValidationError("Shielding keywords do not use the actual scenario CPCM solvent")
    lines = [
        f"# {NOTICE}",
        "# Complete explicit cluster; NMR request is analyte-only by reviewed atom indices.",
        f"! {keywords}",
        f"%pal nprocs {nprocs} end",
        f"%maxcore {maxcore_mb}",
    ]
    if basis_block.strip():
        lines.extend((basis_block.strip(), ""))
    if extra_blocks.strip():
        lines.extend((extra_blocks.strip(), ""))
    # ORCA 6.1 requires the EPR/NMR block after the coordinate specification.
    lines.extend((f"* xyzfile {charge} {multiplicity} optimized.xyz", "%eprnmr", "  Ori = GIAO"))
    for element in ("H", "C", "F", "P"):
        indices = list(index_map.get(element, ()))
        if indices:
            lines.append(
                f"  Nuclei = {','.join(str(value) for value in indices)} "
                f"{{ shift, ist = {ISOTOPE[element]} }};"
            )
    lines.extend(("end", ""))
    return "\n".join(lines)


def render_j_input(
    charge: int,
    multiplicity: int,
    continuum: str,
    index_map: Mapping[str, Sequence[int]],
    pairs: Sequence[Mapping[str, object]],
    nprocs: int,
    maxcore_mb: int,
    simple_keywords: str | None = None,
    basis_block: str = "",
) -> str:
    if not pairs:
        raise WorkflowValidationError("Cannot render a J job with an empty pair allowlist")
    solvent = continuum_token(continuum)
    keywords = (
        simple_keywords.replace("{CONTINUUM}", solvent)
        if simple_keywords
        else f"PBE pcJ-3 AutoAux VeryTightSCF DefGrid3 CPCM({solvent})"
    )
    if f"CPCM({solvent})" not in keywords.upper():
        raise WorkflowValidationError("J keywords do not use the actual scenario CPCM solvent")
    involved = {
        int(row[field]) + 1
        for row in pairs
        for field in ("left_atom_index_0based_spinspin", "right_atom_index_0based_spinspin")
    }
    lines = [
        f"# {NOTICE}",
        "# Complete explicit cluster; Nuclei are 1-based, SpinSpinAtomPairs are 0-based.",
        f"! {keywords}",
        f"%pal nprocs {nprocs} end",
        f"%maxcore {maxcore_mb}",
        "%method Z_Tol 1e-12 end",
    ]
    if basis_block.strip():
        lines.extend((basis_block.strip(), ""))
    lines.extend((f"* xyzfile {charge} {multiplicity} optimized.xyz", "%eprnmr", "  Tol = 1e-10"))
    for element in ("H", "C", "F", "P"):
        indices = [value for value in index_map.get(element, ()) if value in involved]
        if indices:
            lines.append(
                f"  Nuclei = {','.join(str(value) for value in indices)} "
                f"{{ ssall, ist = {ISOTOPE[element]} }};"
            )
    pair_text = " ".join(
        f"{{{row['left_atom_index_0based_spinspin']} {row['right_atom_index_0based_spinspin']}}}"
        for row in pairs
    )
    lines.extend((f"  SpinSpinAtomPairs {pair_text}", "  SpinSpinRThresh 6.05", "end", ""))
    return "\n".join(lines)


def _approved_registry(path: Path, project: Path) -> tuple[dict, list[dict[str, object]]]:
    registry = read_json(path)
    if registry.get("status") != PROTOCOL_APPROVED:
        raise WorkflowValidationError(
            f"Reference registry is not {PROTOCOL_APPROVED}; production property build is blocked: {path}"
        )
    entries = registry.get("reference_jobs")
    if not isinstance(entries, list) or not entries:
        raise WorkflowValidationError(f"Approved reference registry has no entries: {path}")
    validated: list[dict[str, object]] = []
    for entry in entries:
        if not isinstance(entry, dict) or entry.get("approval_status") != PROTOCOL_APPROVED:
            continue
        input_path = ensure_within(Path(str(entry.get("input_file", ""))), project, "reference input")
        output_path = ensure_within(Path(str(entry.get("output_file", ""))), project, "reference output")
        if sha256_file(input_path) != normalize_sha256(str(entry.get("input_sha256", ""))):
            raise WorkflowValidationError(f"Reference input hash mismatch: {entry.get('reference_id')}")
        if sha256_file(output_path) != normalize_sha256(str(entry.get("output_sha256", ""))):
            raise WorkflowValidationError(f"Reference output hash mismatch: {entry.get('reference_id')}")
        if "ORCA TERMINATED NORMALLY" not in read_text_auto(output_path):
            raise WorkflowValidationError(
                f"Reference output lacks normal ORCA termination: {entry.get('reference_id')}"
            )
        validated.append(entry)
    return registry, validated


def _find_reference(
    entries: Sequence[Mapping[str, object]], nucleus: str, domain: str, fingerprint: str
) -> Mapping[str, object]:
    matches = [
        entry
        for entry in entries
        if entry.get("nucleus") == nucleus
        and entry.get("solvent_domain") == domain
        and entry.get("method_fingerprint") == fingerprint
    ]
    if len(matches) != 1:
        raise WorkflowValidationError(
            f"Expected exactly one approved reference for {nucleus}/{domain}/{fingerprint}; found {len(matches)}"
        )
    return matches[0]


def _approved_j_gate(path: Path, project: Path) -> dict:
    gate = read_json(path)
    if gate.get("approval_status") != PROTOCOL_APPROVED:
        raise WorkflowValidationError(f"J benchmark gate is not {PROTOCOL_APPROVED}: {path}")
    for field in (
        "representative_set_coverage_check",
        "numerical_agreement_review",
        "production_method_lock_check",
        "pair_screen_sensitivity_review",
    ):
        if gate.get(field) != PASS:
            raise WorkflowValidationError(f"J benchmark {field} is not PASS: {path}")
    if gate.get("production_method_fingerprint") != GENERIC_J_FINGERPRINT:
        raise WorkflowValidationError(f"J production method fingerprint mismatch: {path}")
    if not str(gate.get("reviewer", "")).strip() or not str(gate.get("reviewed_at", "")).strip():
        raise WorkflowValidationError(f"J benchmark gate lacks reviewer/time: {path}")
    records = gate.get("representative_case_results")
    if not isinstance(records, list) or not records:
        raise WorkflowValidationError(f"J benchmark gate has no representative evidence: {path}")
    for record in records:
        if not isinstance(record, dict) or record.get("approval_status") != PROTOCOL_APPROVED:
            raise WorkflowValidationError(f"Unapproved J benchmark record: {path}")
        for prefix in ("production", "benchmark"):
            input_path = ensure_within(
                Path(str(record.get(f"{prefix}_input_file", ""))), project, f"J {prefix} input"
            )
            output_path = ensure_within(
                Path(str(record.get(f"{prefix}_output_file", ""))), project, f"J {prefix} output"
            )
            if sha256_file(input_path) != normalize_sha256(str(record.get(f"{prefix}_input_sha256", ""))):
                raise WorkflowValidationError(f"J benchmark {prefix} input hash mismatch")
            if sha256_file(output_path) != normalize_sha256(str(record.get(f"{prefix}_output_sha256", ""))):
                raise WorkflowValidationError(f"J benchmark {prefix} output hash mismatch")
            if "ORCA TERMINATED NORMALLY" not in read_text_auto(output_path):
                raise WorkflowValidationError(f"J benchmark {prefix} output lacks normal termination")
    return gate


def _iodine_branches(path: Path, project: Path) -> tuple[list[dict[str, str]], dict[str, str]]:
    gate = read_json(path)
    if gate.get("approval_status") != PROTOCOL_APPROVED:
        raise WorkflowValidationError(f"F3_006 iodine/X2C gate is not approved: {path}")
    for field in (
        "basis_coverage_check",
        "picture_change_and_relativistic_review",
        "x2c_sensitivity_review",
        "j_basis_compatibility_review",
    ):
        if gate.get(field) != PASS:
            raise WorkflowValidationError(f"F3_006 iodine gate {field} is not PASS")
    if not str(gate.get("reviewer", "")).strip() or not str(gate.get("reviewed_at", "")).strip():
        raise WorkflowValidationError("F3_006 iodine gate lacks reviewer/time")
    records = gate.get("syntax_test_records")
    if not isinstance(records, list) or not records:
        raise WorkflowValidationError("F3_006 iodine gate lacks syntax-test evidence")
    for record in records:
        if not isinstance(record, dict) or record.get("approval_status") != PROTOCOL_APPROVED:
            raise WorkflowValidationError("F3_006 syntax-test record is not approved")
        for role in ("input", "output"):
            evidence_path = ensure_within(
                Path(str(record.get(f"{role}_file", ""))), project, f"iodine syntax {role}"
            )
            if sha256_file(evidence_path) != normalize_sha256(str(record.get(f"{role}_sha256", ""))):
                raise WorkflowValidationError(f"F3_006 syntax {role} hash mismatch")
            if role == "output" and "ORCA TERMINATED NORMALLY" not in read_text_auto(evidence_path):
                raise WorkflowValidationError("F3_006 syntax output lacks normal termination")
    branches = []
    for tag, prefix in (("iodine_primary", "primary_shielding"), ("iodine_x2c", "x2c_sensitivity")):
        fingerprint = str(gate.get(f"{prefix}_method_fingerprint", "")).strip()
        keywords = str(gate.get(f"{prefix}_simple_keywords", "")).strip()
        basis = str(gate.get(f"{prefix}_basis_block", "")).strip()
        if not fingerprint or fingerprint.startswith("TO_") or not keywords or not basis:
            raise WorkflowValidationError(f"F3_006 {prefix} method is not fully locked")
        branches.append(
            {"branch": tag, "fingerprint": fingerprint, "keywords": keywords, "basis": basis}
        )
    j = {
        "fingerprint": str(gate.get("production_j_method_fingerprint", "")).strip(),
        "keywords": str(gate.get("production_j_simple_keywords", "")).strip(),
        "basis": str(gate.get("production_j_basis_block", "")).strip(),
    }
    if not all(j.values()) or j["fingerprint"].startswith("TO_"):
        raise WorkflowValidationError("F3_006 iodine J method is not fully locked")
    return branches, j


def _gate_states(
    project: Path,
    gate: Mapping[str, object],
    gate_path: Path,
    figure_rows: list[dict[str, str]],
) -> dict[str, dict[str, object]]:
    records = gate.get("upstream_state_hashes")
    if not isinstance(records, list):
        raise WorkflowValidationError(f"Minimum gate lacks upstream states: {gate_path}")
    rows_by_job = {row["job_id"]: row for row in figure_rows}
    states: dict[str, dict[str, object]] = {}
    for record in records:
        if not isinstance(record, dict):
            raise WorkflowValidationError(f"Malformed upstream state record: {gate_path}")
        job_id = str(record.get("job_id", ""))
        if job_id not in rows_by_job or job_id in states:
            raise WorkflowValidationError(f"Unknown/duplicate opt/freq state in gate: {job_id}")
        state_path = ensure_within(Path(str(record.get("state_file", ""))), project, "minimum-gate state")
        if sha256_file(state_path) != normalize_sha256(str(record.get("sha256", ""))):
            raise WorkflowValidationError(f"Minimum-gate state hash mismatch: {job_id}")
        evidence = validate_completed_opt_freq_state(rows_by_job[job_id], state_path, project)
        states[job_id] = evidence
    if set(states) != set(rows_by_job):
        raise WorkflowValidationError("Approved minimum gate must bind every Figure opt/freq state exactly once")
    return states


def _validate_approved_minimum_gate(
    project: Path,
    gate_path: Path,
    dft_rows: list[dict[str, str]],
) -> tuple[dict, dict[str, dict[str, object]], dict[str, dict[str, object]], dict[tuple[str, str], dict]]:
    gate = read_json(gate_path)
    if gate.get("approval_status") != MINIMUM_APPROVED:
        raise WorkflowValidationError(f"Minimum gate is not {MINIMUM_APPROVED}: {gate_path}")
    if gate.get("same_composition_only_check") != PASS:
        raise WorkflowValidationError(f"Minimum gate same-composition check is not PASS: {gate_path}")
    if not str(gate.get("reviewer", "")).strip() or not str(gate.get("reviewed_at", "")).strip():
        raise WorkflowValidationError(f"Minimum gate lacks reviewer/time: {gate_path}")
    figure_id = str(gate.get("figure_id", ""))
    figure_rows = [row for row in dft_rows if row["figure_id"] == figure_id]
    if not figure_rows:
        raise WorkflowValidationError(f"No DFT rows for minimum gate Figure {figure_id}")
    states = _gate_states(project, gate, gate_path, figure_rows)
    rows_by_job = {row["job_id"]: row for row in figure_rows}
    decisions_payload = gate.get("conformer_decisions")
    if not isinstance(decisions_payload, list):
        raise WorkflowValidationError(f"Minimum gate lacks conformer decisions: {gate_path}")
    decisions: dict[str, dict[str, object]] = {}
    approved: dict[str, dict[str, object]] = {}
    for decision in decisions_payload:
        if not isinstance(decision, dict):
            raise WorkflowValidationError(f"Malformed conformer decision: {gate_path}")
        job_id = str(decision.get("upstream_dft_job_id", ""))
        if job_id not in rows_by_job or job_id in decisions:
            raise WorkflowValidationError(f"Unknown/duplicate conformer decision: {job_id}")
        if decision.get("condition_scenario_id") != rows_by_job[job_id]["scenario_id"]:
            raise WorkflowValidationError(f"Conformer decision crosses scenario: {job_id}")
        evidence = states[job_id]
        if Path(str(decision.get("opt_freq_state_file", ""))).resolve() != Path(evidence["state_path"]).resolve():
            raise WorkflowValidationError(f"Conformer state path mismatch: {job_id}")
        if normalize_sha256(str(decision.get("opt_freq_state_sha256", ""))) != evidence["state_sha256"]:
            raise WorkflowValidationError(f"Conformer state hash mismatch: {job_id}")
        optimized = Path(evidence["artifacts"]["optimized_xyz"]["path"])
        if Path(str(decision.get("optimized_coordinate_file", ""))).resolve() != optimized.resolve():
            raise WorkflowValidationError(f"Conformer optimized XYZ path mismatch: {job_id}")
        if normalize_sha256(str(decision.get("optimized_coordinate_sha256", ""))) != evidence["artifacts"]["optimized_xyz"]["sha256"]:
            raise WorkflowValidationError(f"Conformer optimized XYZ hash mismatch: {job_id}")
        if Path(str(decision.get("hessian_file", ""))).resolve() != Path(
            evidence["artifacts"]["hess"]["path"]
        ).resolve():
            raise WorkflowValidationError(f"Conformer Hessian path mismatch: {job_id}")
        if normalize_sha256(str(decision.get("hessian_sha256", ""))) != evidence["artifacts"]["hess"]["sha256"]:
            raise WorkflowValidationError(f"Conformer Hessian hash mismatch: {job_id}")
        if str(decision.get("conformer_id", "")) != job_id:
            raise WorkflowValidationError(f"Conformer ID is not the bound opt/freq job ID: {job_id}")
        map_path = Path(evidence["job_directory"]) / "nmr_atom_map_v2.csv"
        if Path(str(decision.get("nmr_atom_map_file", ""))).resolve() != map_path.resolve():
            raise WorkflowValidationError(f"Conformer NMR map path mismatch: {job_id}")
        if normalize_sha256(str(decision.get("nmr_atom_map_sha256", ""))) != sha256_file(map_path):
            raise WorkflowValidationError(f"Conformer NMR map hash mismatch: {job_id}")
        by_nucleus = decision.get("thermochemistry_by_nucleus")
        if not isinstance(by_nucleus, dict) or not by_nucleus:
            raise WorkflowValidationError(
                f"Conformer lacks authoritative nucleus-specific thermochemistry: {job_id}"
            )
        parsed_by_nucleus: dict[str, dict[str, object]] = {}
        for nucleus, recorded in by_nucleus.items():
            if not isinstance(recorded, dict):
                raise WorkflowValidationError(f"Malformed {nucleus} thermochemistry: {job_id}")
            temperature = float(recorded.get("temperature_K"))
            parsed = parse_orca_qrrho(Path(evidence["output_path"]), temperature)
            for field in (
                "electronic_energy_hartree",
                "qrrho_gibbs_free_energy_hartree",
                "temperature_K",
            ):
                if abs(float(recorded.get(field)) - float(parsed[field])) > 1e-8:
                    raise WorkflowValidationError(
                        f"Conformer {nucleus}/{field} mismatch: {job_id}"
                    )
            parsed_by_nucleus[str(nucleus)] = parsed
        primary_nucleus = str(decision.get("primary_thermochemistry_nucleus_for_legacy_fields", ""))
        if primary_nucleus not in parsed_by_nucleus:
            raise WorkflowValidationError(f"Invalid primary thermochemistry nucleus: {job_id}")
        thermo = parsed_by_nucleus[primary_nucleus]
        for field in ("electronic_energy_hartree", "qrrho_gibbs_free_energy_hartree", "temperature_K"):
            if abs(float(decision.get(field)) - float(thermo[field])) > 1e-8:
                raise WorkflowValidationError(f"Conformer legacy {field} mismatch: {job_id}")
        if int(decision.get("imaginary_frequency_count", -1)) != int(thermo["imaginary_frequency_count"]):
            raise WorkflowValidationError(f"Conformer imaginary-frequency count mismatch: {job_id}")
        status = decision.get("approval_status")
        if status == MINIMUM_APPROVED:
            if (
                int(thermo["imaginary_frequency_count"]) != 0
                or decision.get("stationary_point_check") != PASS
                or decision.get("cluster_integrity_check") != PASS
            ):
                raise WorkflowValidationError(f"Approved NMR conformer failed minimum checks: {job_id}")
            approved[job_id] = decision
        elif status == "REJECTED_AFTER_OPT_FREQ":
            if not str(decision.get("decision_reason", "")).strip():
                raise WorkflowValidationError(f"Rejected conformer lacks reason: {job_id}")
        else:
            raise WorkflowValidationError(f"Unresolved conformer status in approved gate: {job_id}")
        decisions[job_id] = decision
    if set(decisions) != set(rows_by_job):
        raise WorkflowValidationError("Approved gate must decide every completed opt/freq conformer")
    for scenario in {row["scenario_id"] for row in figure_rows}:
        if not any(rows_by_job[job]["scenario_id"] == scenario for job in approved):
            raise WorkflowValidationError(f"No approved NMR conformer remains for scenario {scenario}")

    expected_weight_keys: set[tuple[str, str]] = set()
    nuclei_by_scenario: dict[str, set[str]] = defaultdict(set)
    index_signature_by_scenario: dict[str, dict[str, list[int]]] = {}
    for job_id in approved:
        row = rows_by_job[job_id]
        index_map = analyte_index_map(Path(states[job_id]["job_directory"]) / "nmr_atom_map_v2.csv")
        scenario_id = row["scenario_id"]
        if scenario_id in index_signature_by_scenario and index_signature_by_scenario[scenario_id] != index_map:
            raise WorkflowValidationError(f"Analyte nucleus/index map differs within scenario {scenario_id}")
        index_signature_by_scenario[scenario_id] = index_map
        nuclei_by_scenario[scenario_id].update(NUCLEUS[element] for element in index_map)
    for scenario, nuclei in nuclei_by_scenario.items():
        expected_weight_keys.update((scenario, nucleus) for nucleus in nuclei)
    weight_sets_payload = gate.get("boltzmann_weight_sets")
    if not isinstance(weight_sets_payload, list):
        raise WorkflowValidationError("Approved gate lacks Boltzmann weight sets")
    weight_sets: dict[tuple[str, str], dict] = {}
    for weight_set in weight_sets_payload:
        if not isinstance(weight_set, dict) or weight_set.get("approval_status") != MINIMUM_APPROVED:
            raise WorkflowValidationError("Boltzmann weight set is not approved")
        key = (str(weight_set.get("condition_scenario_id", "")), str(weight_set.get("nucleus", "")))
        if key in weight_sets:
            raise WorkflowValidationError(f"Duplicate Boltzmann weight set: {key}")
        temperature = float(weight_set.get("temperature_K"))
        approved_jobs = sorted(
            job for job in approved if rows_by_job[job]["scenario_id"] == key[0]
        )
        supplied = weight_set.get("weights")
        if not isinstance(supplied, list):
            raise WorkflowValidationError(f"Malformed weight records: {key}")
        supplied_map = {str(item["conformer_id"]): float(item["weight"]) for item in supplied}
        if set(supplied_map) != set(approved_jobs):
            raise WorkflowValidationError(f"Weight conformer set differs from approved set: {key}")
        energies: dict[str, float] = {}
        for job in approved_jobs:
            thermochemistry = approved[job].get("thermochemistry_by_nucleus")
            if not isinstance(thermochemistry, dict) or key[1] not in thermochemistry:
                raise WorkflowValidationError(
                    f"Weight set {key} lacks matching nucleus thermochemistry for {job}"
                )
            record = thermochemistry[key[1]]
            if not isinstance(record, dict):
                raise WorkflowValidationError(f"Malformed thermochemistry record for {job}/{key[1]}")
            if abs(float(record.get("temperature_K")) - temperature) > 0.02:
                raise WorkflowValidationError(f"Weight temperature lacks matching thermochemistry: {key}")
            energies[job] = float(record["qrrho_gibbs_free_energy_hartree"])
        expected = normalized_weights(energies, temperature)
        if not math.isclose(math.fsum(supplied_map.values()), 1.0, rel_tol=0.0, abs_tol=1e-9):
            raise WorkflowValidationError(f"Weights do not sum to one: {key}")
        if any(abs(supplied_map[job] - expected[job]) > 1e-8 for job in approved_jobs):
            raise WorkflowValidationError(f"Weights do not reproduce QRRHO Boltzmann values: {key}")
        weight_sets[key] = weight_set
    if set(weight_sets) != expected_weight_keys:
        raise WorkflowValidationError(
            f"Weight-set coverage mismatch; expected={sorted(expected_weight_keys)}, got={sorted(weight_sets)}"
        )
    return gate, states, approved, weight_sets


def _resource_values(path: Path) -> tuple[dict[str, int], str]:
    policy = read_json(path)
    cap = int(policy["execution"]["absolute_mpi_rank_cap"])
    execution = policy["execution"]
    legacy_serial = execution.get("one_compound_at_a_time") is True
    hpc_array = (
        execution.get("execution_model") == "one_orca_calculation_per_scheduler_array_task"
        and execution.get("one_orca_calculation_per_array_task") is True
    )
    if cap > 32 or not (legacy_serial or hpc_array):
        raise WorkflowValidationError(
            "Property resource policy must use legacy single-instance execution or one ORCA "
            "calculation per scheduler-array task, with <=32 MPI ranks"
        )
    result = {}
    for stage in ("shielding", "j_coupling"):
        values = policy["stage_resources"][stage]
        ranks = int(values["mpi_ranks"])
        if ranks > cap:
            raise WorkflowValidationError(f"{stage} ranks exceed resource cap")
        result[f"{stage}_ranks"] = ranks
        result[f"{stage}_maxcore"] = int(values["maxcore_mb_per_rank"])
    return result, sha256_file(path)


def build_property_jobs(
    project: Path,
    gate_paths: list[Path],
    dft_manifest_path: Path,
    reference_registry_path: Path,
    j_gate_path: Path,
    iodine_gate_path: Path,
    resource_policy_path: Path,
    output_root: Path,
    shielding_manifest_path: Path,
    j_manifest_path: Path,
    weights_path: Path,
    replace_unstarted: bool = False,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    project = Path(project).resolve()
    if not gate_paths:
        raise WorkflowValidationError("No approved minimum/Boltzmann gates were supplied")
    dft_manifest_path = ensure_within(dft_manifest_path, project, "DFT manifest")
    reference_registry_path = ensure_within(reference_registry_path, project, "reference registry")
    resource_policy_path = ensure_within(resource_policy_path, project, "resource policy")
    output_root = ensure_within(output_root, project, "property output root")
    shielding_manifest_path = ensure_within(shielding_manifest_path, project, "shielding manifest")
    j_manifest_path = ensure_within(j_manifest_path, project, "J manifest")
    weights_path = ensure_within(weights_path, project, "weight table")
    gate_paths = [ensure_within(path, project, "minimum gate") for path in gate_paths]
    dft_rows = read_csv(dft_manifest_path)
    _registry, reference_entries = _approved_registry(reference_registry_path, project)
    resources, resource_hash = _resource_values(resource_policy_path)
    dft_manifest_hash = sha256_file(dft_manifest_path)
    reference_registry_hash = sha256_file(reference_registry_path)

    gate_records = []
    seen_figures: set[str] = set()
    any_j = False
    need_iodine = False
    for gate_path in gate_paths:
        gate, states, approved, weight_sets = _validate_approved_minimum_gate(
            project, gate_path, dft_rows
        )
        figure_id = str(gate["figure_id"])
        if figure_id in seen_figures:
            raise WorkflowValidationError(f"Duplicate approved minimum gate for {figure_id}")
        seen_figures.add(figure_id)
        rows_by_job = {row["job_id"]: row for row in dft_rows if row["figure_id"] == figure_id}
        for job_id in approved:
            index_map = analyte_index_map(Path(states[job_id]["job_directory"]) / "nmr_atom_map_v2.csv")
            pairs = j_pair_allowlist(index_map, Path(states[job_id]["artifacts"]["optimized_xyz"]["path"]))
            any_j = any_j or bool(pairs)
            frame = parse_xyz_frames(Path(states[job_id]["artifacts"]["optimized_xyz"]["path"]))[0]
            if "I" in {symbol.upper() for symbol in frame.symbols}:
                need_iodine = True
                if figure_id != "F3_006":
                    raise WorkflowValidationError(f"Unexpected iodine outside F3_006: {job_id}")
        gate_records.append((gate_path, gate, states, approved, weight_sets, rows_by_job))
    resolved_j_gate = ensure_within(j_gate_path, project, "J benchmark gate") if any_j else None
    resolved_iodine_gate = (
        ensure_within(iodine_gate_path, project, "F3_006 iodine gate") if need_iodine else None
    )
    # A draft J benchmark gate blocks J generation only.  It must not block
    # independently referenced shielding jobs that are already approved.
    j_gate = None
    if resolved_j_gate is not None:
        raw_j_gate = read_json(resolved_j_gate)
        if raw_j_gate.get("approval_status") == PROTOCOL_APPROVED:
            j_gate = _approved_j_gate(resolved_j_gate, project)
    iodine_branches = (
        _iodine_branches(resolved_iodine_gate, project)
        if resolved_iodine_gate is not None
        else None
    )
    j_gate_hash = sha256_file(resolved_j_gate) if j_gate is not None else ""
    iodine_gate_hash = sha256_file(resolved_iodine_gate) if iodine_branches is not None else ""

    shielding_plans: list[Plan] = []
    j_plans: list[Plan] = []
    weight_rows: list[dict[str, object]] = []
    for gate_path, gate, states, approved, weight_sets, rows_by_job in gate_records:
        gate_hash = sha256_file(gate_path)
        for key, weight_set in sorted(weight_sets.items()):
            for item in weight_set["weights"]:
                weight_rows.append(
                    {
                        "figure_id": gate["figure_id"],
                        "scenario_id": key[0],
                        "nucleus": key[1],
                        "temperature_K": weight_set["temperature_K"],
                        "conformer_id": item["conformer_id"],
                        "qrrho_gibbs_free_energy_hartree": item.get(
                            "qrrho_gibbs_free_energy_hartree",
                            approved[item["conformer_id"]]["qrrho_gibbs_free_energy_hartree"],
                        ),
                        "weight": item["weight"],
                        "approved_gate_sha256": gate_hash,
                        "notice": NOTICE,
                    }
                )
        for job_id, decision in sorted(approved.items()):
            row = rows_by_job[job_id]
            evidence = states[job_id]
            scenario_id = row["scenario_id"]
            optimized_path = Path(evidence["artifacts"]["optimized_xyz"]["path"])
            optimized_bytes = optimized_path.read_bytes()
            optimized_hash = evidence["artifacts"]["optimized_xyz"]["sha256"]
            map_path = Path(evidence["job_directory"]) / "nmr_atom_map_v2.csv"
            map_bytes = map_path.read_bytes()
            map_hash = sha256_file(map_path)
            isotope_map_path = Path(evidence["job_directory"]) / "source_isotope_map.csv"
            if not isotope_map_path.is_file():
                raise WorkflowValidationError(
                    f"Missing source isotope map for property handoff: {job_id}"
                )
            isotope_map_bytes = isotope_map_path.read_bytes()
            isotope_map_hash = sha256_file(isotope_map_path)
            index_map = analyte_index_map(map_path)
            if not index_map:
                continue
            pairs = j_pair_allowlist(index_map, optimized_path)
            continuum = continuum_token(row["continuum_keyword"])
            domain = solvent_domain(row.get("primary_solvent", ""), continuum)
            try:
                charge = int(row["charge"])
                multiplicity = int(row["multiplicity"])
            except (TypeError, ValueError) as exc:
                raise WorkflowValidationError(f"Invalid property charge/multiplicity: {job_id}") from exc
            contains_iodine = "I" in {
                symbol.upper() for symbol in parse_xyz_frames(optimized_path)[0].symbols
            }
            if contains_iodine:
                assert iodine_branches is not None
                shield_branches, iodine_j = iodine_branches
            else:
                shield_branches = shielding_branches(continuum)
                iodine_j = None
            weights_json = {
                NUCLEUS[element]: float(
                    next(
                        item["weight"]
                        for item in weight_sets[(scenario_id, NUCLEUS[element])]["weights"]
                        if item["conformer_id"] == job_id
                    )
                )
                for element in index_map
            }
            request_rows = [
                {
                    "element": element,
                    "nucleus": NUCLEUS[element],
                    "isotope": ISOTOPE[element],
                    "atom_indices_1based_nuclei": ",".join(str(value) for value in indices),
                    "scope": "ANALYTE_ONLY",
                }
                for element, indices in sorted(index_map.items())
            ]
            for branch in shield_branches:
                references = {
                    NUCLEUS[element]: _find_reference(
                        reference_entries, NUCLEUS[element], domain, branch["fingerprint"]
                    )
                    for element in index_map
                }
                text = render_shielding_input(
                    charge,
                    multiplicity,
                    continuum,
                    index_map,
                    resources["shielding_ranks"],
                    resources["shielding_maxcore"],
                    branch["keywords"] or None,
                    branch["basis"],
                    branch.get("extra_blocks", ""),
                )
                input_bytes = text.encode("utf-8")
                bundle_hash = orca_input_bundle_sha256(
                    "shielding.inp", input_bytes, (("optimized.xyz", optimized_bytes),)
                )
                property_job_id = f"{job_id}__shielding__{safe_token(branch['branch'])}"
                directory = (
                    output_root
                    / safe_token(str(gate["figure_id"]))
                    / safe_token(scenario_id)
                    / safe_token(job_id)
                    / f"shielding_{safe_token(branch['branch'])}"
                )
                provenance = {
                    "schema_version": "2.0",
                    "notice": NOTICE,
                    "job_id": property_job_id,
                    "figure_id": gate["figure_id"],
                    "scenario_id": scenario_id,
                    "source_opt_freq_job_id": job_id,
                    "approved_minimum_gate_file": str(gate_path),
                    "approved_minimum_gate_sha256": gate_hash,
                    "source_opt_freq_state_sha256": evidence["state_sha256"],
                    "optimized_coordinate_sha256": optimized_hash,
                    "nmr_atom_map_sha256": map_hash,
                    "source_isotope_map_sha256": isotope_map_hash,
                    "requested_indices_1based": index_map,
                    "environment_nuclei_requested": False,
                    "method_fingerprint": branch["fingerprint"],
                    "solvent_domain": domain,
                    "continuum_keyword": continuum,
                    "weights_by_nucleus": weights_json,
                    "reference_jobs": references,
                    "reference_registry_sha256": reference_registry_hash,
                    "source_dft_manifest_sha256": dft_manifest_hash,
                    "resource_policy_sha256": resource_hash,
                    "iodine_gate_sha256": iodine_gate_hash if contains_iodine else "",
                }
                files = {
                    "optimized.xyz": optimized_bytes,
                    "nmr_atom_map_v2.csv": map_bytes,
                    "source_isotope_map.csv": isotope_map_bytes,
                    "shielding_request_map.csv": _csv_bytes(
                        request_rows,
                        (
                            "element",
                            "nucleus",
                            "isotope",
                            "atom_indices_1based_nuclei",
                            "scope",
                        ),
                    ),
                    "shielding.inp": input_bytes,
                    "provenance.json": (
                        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n"
                    ).encode("utf-8"),
                    "NOTICE.txt": (NOTICE + "\n").encode("utf-8"),
                }
                shielding_plans.append(
                    Plan(
                        directory,
                        files,
                        {
                            "job_id": property_job_id,
                            "figure_id": gate["figure_id"],
                            "scenario_id": scenario_id,
                            "conformer_id": job_id,
                            "property_branch": branch["branch"],
                            "primary_solvent": row.get("primary_solvent", ""),
                            "continuum_keyword": continuum,
                            "job_directory": str(directory.resolve()),
                            "input_file": "shielding.inp",
                            "output_file": "shielding.out",
                            "required_artifacts": json.dumps(
                                ["shielding.gbw", "shielding.property.txt"], separators=(",", ":")
                            ),
                            "input_sha256": sha256_bytes(input_bytes),
                            "expected_input_bundle_sha256": bundle_hash,
                            "conformer_xyz_sha256": optimized_hash,
                            "nmr_atom_map_sha256": map_hash,
                            "source_isotope_map_sha256": isotope_map_hash,
                            "approved_minimum_gate_sha256": gate_hash,
                            "source_opt_freq_state_sha256": evidence["state_sha256"],
                            "method_fingerprint": branch["fingerprint"],
                            "reference_ids": json.dumps(
                                {key: value["reference_id"] for key, value in references.items()},
                                separators=(",", ":"),
                            ),
                            "weights_by_nucleus": json.dumps(weights_json, separators=(",", ":")),
                            "enabled": "true",
                            "readiness": "ready",
                            "approval_status": PROTOCOL_APPROVED,
                            "notice": NOTICE,
                        },
                    )
                )
            if pairs and j_gate is not None:
                if contains_iodine:
                    assert iodine_j is not None
                    j_fingerprint = iodine_j["fingerprint"]
                    j_keywords = iodine_j["keywords"]
                    j_basis = iodine_j["basis"]
                else:
                    j_fingerprint = GENERIC_J_FINGERPRINT
                    j_keywords = ""
                    j_basis = ""
                text = render_j_input(
                    charge,
                    multiplicity,
                    continuum,
                    index_map,
                    pairs,
                    resources["j_coupling_ranks"],
                    resources["j_coupling_maxcore"],
                    j_keywords or None,
                    j_basis,
                )
                input_bytes = text.encode("utf-8")
                bundle_hash = orca_input_bundle_sha256(
                    "j_coupling.inp", input_bytes, (("optimized.xyz", optimized_bytes),)
                )
                property_job_id = f"{job_id}__j_coupling"
                directory = (
                    output_root
                    / safe_token(str(gate["figure_id"]))
                    / safe_token(scenario_id)
                    / safe_token(job_id)
                    / "j_coupling"
                )
                pair_fields = (
                    "left_atom_index_1based_nuclei",
                    "right_atom_index_1based_nuclei",
                    "left_atom_index_0based_spinspin",
                    "right_atom_index_0based_spinspin",
                    "left_element",
                    "right_element",
                    "pair_class",
                    "distance_A",
                    "scope",
                )
                provenance = {
                    "schema_version": "2.0",
                    "notice": NOTICE,
                    "job_id": property_job_id,
                    "figure_id": gate["figure_id"],
                    "scenario_id": scenario_id,
                    "source_opt_freq_job_id": job_id,
                    "approved_minimum_gate_file": str(gate_path),
                    "approved_minimum_gate_sha256": gate_hash,
                    "source_opt_freq_state_sha256": evidence["state_sha256"],
                    "optimized_coordinate_sha256": optimized_hash,
                    "nmr_atom_map_sha256": map_hash,
                    "source_isotope_map_sha256": isotope_map_hash,
                    "method_fingerprint": j_fingerprint,
                    "solvent_domain": domain,
                    "continuum_keyword": continuum,
                    "pair_count": len(pairs),
                    "pair_index_policy": "Nuclei=1-based; SpinSpinAtomPairs=0-based (ORCA 6.1 manual)",
                    "environment_nuclei_requested": False,
                    "j_benchmark_gate_sha256": j_gate_hash,
                    "iodine_gate_sha256": iodine_gate_hash if contains_iodine else "",
                    "weights_by_nucleus": weights_json,
                    "resource_policy_sha256": resource_hash,
                    "source_dft_manifest_sha256": dft_manifest_hash,
                }
                files = {
                    "optimized.xyz": optimized_bytes,
                    "nmr_atom_map_v2.csv": map_bytes,
                    "source_isotope_map.csv": isotope_map_bytes,
                    "j_pair_allowlist.csv": _csv_bytes(pairs, pair_fields),
                    "j_coupling.inp": input_bytes,
                    "provenance.json": (
                        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n"
                    ).encode("utf-8"),
                    "NOTICE.txt": (NOTICE + "\n").encode("utf-8"),
                }
                j_plans.append(
                    Plan(
                        directory,
                        files,
                        {
                            "job_id": property_job_id,
                            "figure_id": gate["figure_id"],
                            "scenario_id": scenario_id,
                            "conformer_id": job_id,
                            "primary_solvent": row.get("primary_solvent", ""),
                            "continuum_keyword": continuum,
                            "job_directory": str(directory.resolve()),
                            "input_file": "j_coupling.inp",
                            "output_file": "j_coupling.out",
                            "required_artifacts": json.dumps(
                                ["j_coupling.gbw", "j_coupling.property.txt"], separators=(",", ":")
                            ),
                            "input_sha256": sha256_bytes(input_bytes),
                            "expected_input_bundle_sha256": bundle_hash,
                            "conformer_xyz_sha256": optimized_hash,
                            "nmr_atom_map_sha256": map_hash,
                            "source_isotope_map_sha256": isotope_map_hash,
                            "approved_minimum_gate_sha256": gate_hash,
                            "source_opt_freq_state_sha256": evidence["state_sha256"],
                            "method_fingerprint": j_fingerprint,
                            "j_pair_count": len(pairs),
                            "j_benchmark_gate_sha256": j_gate_hash,
                            "weights_by_nucleus": json.dumps(weights_json, separators=(",", ":")),
                            "enabled": "true",
                            "readiness": "ready",
                            "approval_status": PROTOCOL_APPROVED,
                            "notice": NOTICE,
                        },
                    )
                )

    for path in (shielding_manifest_path, j_manifest_path, weights_path):
        if path.exists() and not replace_unstarted:
            raise WorkflowValidationError(f"Output table exists; refusing overwrite: {path}")
    for plan in (*shielding_plans, *j_plans):
        if plan.directory.exists():
            started = list(plan.directory.rglob("*.out")) + list(plan.directory.rglob("*.gbw"))
            if started:
                raise WorkflowValidationError(f"Property job already started: {started[0]}")
            if not replace_unstarted:
                raise WorkflowValidationError(f"Property input directory exists: {plan.directory}")
    for plan in (*shielding_plans, *j_plans):
        for relative, data in plan.files.items():
            atomic_write_bytes(plan.directory / relative, data)
    shielding_rows = [dict(plan.row, queue_order=index) for index, plan in enumerate(shielding_plans, 1)]
    j_rows = [dict(plan.row, queue_order=index) for index, plan in enumerate(j_plans, 1)]
    shielding_fields = list(shielding_rows[0]) if shielding_rows else ["queue_order", "job_id", "notice"]
    j_fields = list(j_rows[0]) if j_rows else ["queue_order", "job_id", "notice"]
    weight_fields = list(weight_rows[0]) if weight_rows else [
        "figure_id",
        "scenario_id",
        "nucleus",
        "temperature_K",
        "conformer_id",
        "weight",
        "notice",
    ]
    atomic_write_csv(shielding_manifest_path, shielding_rows, shielding_fields)
    atomic_write_csv(j_manifest_path, j_rows, j_fields)
    atomic_write_csv(weights_path, weight_rows, weight_fields)
    return shielding_rows, j_rows, weight_rows


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--project", type=Path, default=PROJECT)
    value.add_argument("--gate", type=Path, action="append", dest="gates")
    value.add_argument("--dft-manifest", type=Path)
    value.add_argument("--reference-registry", type=Path)
    value.add_argument("--j-gate", type=Path)
    value.add_argument("--iodine-gate", type=Path)
    value.add_argument("--resource-policy", type=Path)
    value.add_argument("--output-root", type=Path)
    value.add_argument("--shielding-manifest", type=Path)
    value.add_argument("--j-manifest", type=Path)
    value.add_argument("--weights", type=Path)
    value.add_argument("--replace-unstarted", action="store_true")
    return value


def main() -> None:
    args = parser().parse_args()
    project = args.project.resolve()
    gates = args.gates or sorted(
        (project / "jobs").glob("F3_*/minimum_boltzmann_audit/minimum_boltzmann_approved.json")
    )
    shielding, coupling, weights = build_property_jobs(
        project=project,
        gate_paths=gates,
        dft_manifest_path=args.dft_manifest or project / "data" / "DFT_explicit_job_manifest.csv",
        reference_registry_path=args.reference_registry or project / "config" / "reference_job_registry.json",
        j_gate_path=args.j_gate or project / "jobs" / "j_method_benchmark_gate.json",
        iodine_gate_path=args.iodine_gate or project / "jobs" / "F3_006_iodine_x2c_gate.json",
        resource_policy_path=args.resource_policy or project / "config" / "resource_policy.json",
        output_root=args.output_root or project / "03_nmr_properties_explicit",
        shielding_manifest_path=args.shielding_manifest or project / "data" / "NMR_shielding_job_manifest.csv",
        j_manifest_path=args.j_manifest or project / "data" / "NMR_j_coupling_job_manifest.csv",
        weights_path=args.weights or project / "data" / "NMR_conformer_weights.csv",
        replace_unstarted=args.replace_unstarted,
    )
    print(
        f"Built {len(shielding)} shielding and {len(coupling)} J inputs with {len(weights)} weights; "
        "ORCA was not started."
    )


if __name__ == "__main__":
    main()
