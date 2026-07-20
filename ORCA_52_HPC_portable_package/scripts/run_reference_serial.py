"""Validate and, only after review approval, run NMR reference jobs serially.

The default mode is a read-only dry run.  ``--execute`` is required before an
ORCA process can be created, and even then the protocol gate, every model
decision, the manifest rows and the model declarations must all carry the
explicit status ``APPROVED_FOR_REFERENCE_CALCULATION``.  Optimization/frequency
and shielding stages are run one at a time under a project-wide lock.

The calculated raw shielding reference is deliberately kept separate from the
experimental processing/reference convention.  In particular, the declared
microsolvated H3PO4 models are not treated as 85 wt% H3PO4 external standards.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import time
import uuid
from datetime import datetime
from pathlib import Path
from typing import Iterable, Mapping

from draft_minimum_boltzmann_gate import parse_orca_qrrho
from prepare_reference_jobs import intercomponent_contact_audit
from run_opt_freq_serial import (
    ExclusiveLock,
    available_ram_gib,
    cpu_percent,
    running_process_names,
    utc_now,
)
from workflow_schema import (
    DEUTERIUM_MASS_U,
    NOTICE,
    WorkflowValidationError,
    atomic_write_json,
    ensure_within,
    inline_orca_input_bundle_sha256,
    normalize_sha256,
    orca_input_bundle_sha256,
    parse_xyz_frames,
    read_csv,
    read_json,
    sha256_bytes,
    sha256_file,
)


PROJECT = Path(__file__).resolve().parents[1]
APPROVED = "APPROVED_FOR_REFERENCE_CALCULATION"
REJECTED = "REJECTED_REFERENCE_MODEL"
NORMAL_MARKER = "ORCA TERMINATED NORMALLY"
ISOTOPE_BY_ELEMENT = {"H": 1, "C": 13, "F": 19, "P": 31}
NUCLEUS_BY_ELEMENT = {"H": "1H", "C": "13C", "F": "19F", "P": "31P"}


def _read_manifest(path: Path) -> list[dict[str, str]]:
    path = Path(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle)]
    if not rows:
        raise WorkflowValidationError(f"Empty reference manifest: {path}")
    model_ids = [row.get("model_id", "").strip() for row in rows]
    if any(not value for value in model_ids) or len(model_ids) != len(set(model_ids)):
        raise WorkflowValidationError("Reference manifest has blank or duplicate model IDs")
    try:
        orders = [int(row.get("queue_order", "")) for row in rows]
    except ValueError as exc:
        raise WorkflowValidationError("Reference manifest has a non-integer queue order") from exc
    if orders != list(range(1, len(rows) + 1)):
        raise WorkflowValidationError("Reference manifest queue order is not consecutive")
    return rows


def _parse_nuclei_requests(template: str) -> dict[str, list[int]]:
    result: dict[str, list[int]] = {}
    isotope_to_element = {value: key for key, value in ISOTOPE_BY_ELEMENT.items()}
    pattern = re.compile(
        r"Nuclei\s*=\s*([0-9,\s]+)\s*\{\s*shift\s*,\s*ist\s*=\s*(\d+)\s*\}\s*;",
        flags=re.IGNORECASE,
    )
    for indices_text, isotope_text in pattern.findall(template):
        isotope = int(isotope_text)
        if isotope not in isotope_to_element:
            raise WorkflowValidationError(f"Unsupported isotope in reference template: {isotope}")
        element = isotope_to_element[isotope]
        if element in result:
            raise WorkflowValidationError(f"Duplicate Nuclei request for {element}")
        indices = [int(value.strip()) for value in indices_text.split(",") if value.strip()]
        if not indices or len(indices) != len(set(indices)) or any(value < 1 for value in indices):
            raise WorkflowValidationError(f"Invalid one-based Nuclei request for {element}")
        result[element] = indices
    if not result:
        raise WorkflowValidationError("Reference shielding template has no Nuclei request")
    return result


def _validate_atom_map_and_inputs(
    row: Mapping[str, str], job_dir: Path, declaration: Mapping[str, object]
) -> dict[str, object]:
    start_path = job_dir / "start.xyz"
    map_path = job_dir / "reference_atom_map.csv"
    opt_path = job_dir / str(row.get("opt_input_file", ""))
    template_path = job_dir / str(row.get("shielding_template_file", ""))
    frames = parse_xyz_frames(start_path)
    if len(frames) != 1:
        raise WorkflowValidationError(f"Reference start XYZ must contain one frame: {start_path}")
    frame = frames[0]
    atom_rows = read_csv(map_path)
    if len(atom_rows) != len(frame.symbols):
        raise WorkflowValidationError(f"Reference atom-map count differs from XYZ: {row['model_id']}")
    requested_from_map: dict[str, list[int]] = {}
    mass_count = 0
    for expected_index, (atom_row, symbol, coordinates) in enumerate(
        zip(atom_rows, frame.symbols, frame.coordinates, strict=True), start=1
    ):
        if int(atom_row.get("atom_index_1based", "0")) != expected_index:
            raise WorkflowValidationError(f"Non-consecutive reference atom map: {row['model_id']}")
        if atom_row.get("element", "").strip() != symbol:
            raise WorkflowValidationError(f"Reference atom symbol differs from XYZ: {row['model_id']}")
        mapped_xyz = tuple(float(atom_row[key]) for key in ("x_A", "y_A", "z_A"))
        if any(abs(left - right) > 2e-8 for left, right in zip(mapped_xyz, coordinates, strict=True)):
            raise WorkflowValidationError(f"Reference atom coordinates differ from XYZ: {row['model_id']}")
        if atom_row.get("isotope_mass_u", "").strip():
            if symbol != "H" or not math.isclose(
                float(atom_row["isotope_mass_u"]), DEUTERIUM_MASS_U, rel_tol=0.0, abs_tol=1e-11
            ):
                raise WorkflowValidationError(f"Unsupported reference isotope mass: {row['model_id']}")
            mass_count += 1
        nucleus = atom_row.get("requested_nucleus", "").strip()
        site = atom_row.get("reference_site", "").strip().lower() == "true"
        environment = atom_row.get("environment_atom", "").strip().lower() == "true"
        if site != bool(nucleus):
            raise WorkflowValidationError(f"Reference-site/nucleus map disagreement: {row['model_id']}")
        if site:
            if environment:
                raise WorkflowValidationError(f"Environment atom requested as reference site: {row['model_id']}")
            expected_nucleus = NUCLEUS_BY_ELEMENT.get(symbol)
            if nucleus != expected_nucleus:
                raise WorkflowValidationError(f"Reference-site element/isotope mismatch: {row['model_id']}")
            requested_from_map.setdefault(symbol, []).append(expected_index)

    opt_text = opt_path.read_text(encoding="utf-8")
    if "xyzfile" in opt_text.casefold() or "* xyz " not in opt_text.casefold():
        raise WorkflowValidationError(f"Reference opt/freq input must contain inline coordinates: {opt_path}")
    if opt_text.count(f"M={DEUTERIUM_MASS_U:.11f}") != mass_count:
        raise WorkflowValidationError(f"Reference deuterium annotations differ from map: {row['model_id']}")
    template_text = template_path.read_text(encoding="utf-8")
    if "optimized.xyz" not in template_text or "* xyzfile" not in template_text.casefold():
        raise WorkflowValidationError(f"Reference shielding template lacks optimized.xyz: {template_path}")
    requested_from_template = _parse_nuclei_requests(template_text)
    if requested_from_template != requested_from_map:
        raise WorkflowValidationError(f"Reference Nuclei indices differ from atom map: {row['model_id']}")
    requested_nuclei = {NUCLEUS_BY_ELEMENT[element] for element in requested_from_map}
    declared_nuclei = {str(value) for value in declaration.get("requested_nuclei", [])}
    row_nuclei = {value for value in row.get("requested_nuclei", "").split(";") if value}
    if requested_nuclei != declared_nuclei or requested_nuclei != row_nuclei:
        raise WorkflowValidationError(f"Reference requested-nucleus sets differ: {row['model_id']}")
    expected_fingerprint = (
        f"PBE0|pcSseg-2|AutoAux|RIJCOSX|GIAO|"
        f"CPCM({declaration['continuum_keyword']})|TightSCF|DefGrid3|ORCA_6.1.1"
    )
    if declaration.get("method_fingerprint") != expected_fingerprint:
        raise WorkflowValidationError(f"Reference method fingerprint is not canonical: {row['model_id']}")
    required_template_tokens = (
        "PBE0",
        "pcSseg-2",
        "AutoAux",
        "RIJCOSX",
        "TightSCF",
        "DefGrid3",
        f"CPCM({declaration['continuum_keyword']})",
        "Ori = GIAO",
    )
    if any(token.casefold() not in template_text.casefold() for token in required_template_tokens):
        raise WorkflowValidationError(f"Reference template and method fingerprint differ: {row['model_id']}")
    contact_audit = intercomponent_contact_audit(atom_rows)
    declared_contact = declaration.get("initial_cluster_contact_audit")
    if declared_contact != contact_audit:
        raise WorkflowValidationError(f"Reference cluster contact audit differs: {row['model_id']}")
    if contact_audit["automatic_intercomponent_collision_check"] == "FAIL":
        raise WorkflowValidationError(f"Reference cluster has an intercomponent collision: {row['model_id']}")
    if row.get("initial_collision_check") != contact_audit[
        "automatic_intercomponent_collision_check"
    ]:
        raise WorkflowValidationError(f"Reference manifest collision status differs: {row['model_id']}")
    for field in ("minimum_cross_component_distance_A", "minimum_nonbonded_clearance_ratio"):
        declared_value = contact_audit[field]
        row_value = row.get(field, "").strip()
        if declared_value is None:
            if row_value:
                raise WorkflowValidationError(f"Unexpected reference contact metric: {row['model_id']}")
        elif not math.isclose(float(row_value), float(declared_value), rel_tol=0.0, abs_tol=1e-6):
            raise WorkflowValidationError(f"Reference contact metric differs ({field}): {row['model_id']}")
    return {
        "atom_count": len(atom_rows),
        "reference_indices_1based": requested_from_map,
        "isotope_mass_annotation_count": mass_count,
    }


def validate_reference_manifest_job(
    row: Mapping[str, str], project: Path, gate_decision: Mapping[str, object]
) -> dict[str, object]:
    project = Path(project).resolve()
    model_id = str(row.get("model_id", "")).strip()
    job_dir = ensure_within(Path(str(row.get("job_directory", ""))), project, "reference job")
    files_and_hashes = (
        ("start.xyz", "start_xyz_sha256"),
        ("reference_atom_map.csv", "reference_atom_map_sha256"),
        ("model_declaration.json", "model_declaration_sha256"),
        (str(row.get("opt_input_file", "")), "opt_input_sha256"),
        (str(row.get("shielding_template_file", "")), "shielding_template_sha256"),
    )
    paths: dict[str, Path] = {}
    for name, hash_field in files_and_hashes:
        if not name or Path(name).name != name:
            raise WorkflowValidationError(f"Unsafe reference file name ({hash_field}): {model_id}")
        path = ensure_within(job_dir / name, project, f"reference {hash_field}")
        if sha256_file(path) != normalize_sha256(str(row.get(hash_field, "")), hash_field):
            raise WorkflowValidationError(f"Reference file hash mismatch ({hash_field}): {model_id}")
        paths[hash_field] = path
    declaration_path = paths["model_declaration_sha256"]
    declaration = read_json(declaration_path)
    if declaration.get("model_id") != model_id or declaration.get("notice") != NOTICE:
        raise WorkflowValidationError(f"Reference model declaration identity differs: {model_id}")
    if str(declaration.get("solvent_domain", "")) != str(row.get("solvent_domain", "")):
        raise WorkflowValidationError(f"Reference solvent domain differs: {model_id}")
    if str(declaration.get("continuum_keyword", "")) != str(row.get("continuum_keyword", "")):
        raise WorkflowValidationError(f"Reference continuum differs: {model_id}")
    fingerprint = str(row.get("shielding_method_fingerprint", ""))
    if str(declaration.get("method_fingerprint", "")) != fingerprint:
        raise WorkflowValidationError(f"Reference shielding method differs: {model_id}")
    fingerprint_hash = sha256_bytes(fingerprint.encode())
    if fingerprint_hash != normalize_sha256(
        str(row.get("shielding_method_fingerprint_sha256", "")), "method fingerprint hash"
    ) or fingerprint_hash != normalize_sha256(
        str(declaration.get("method_fingerprint_sha256", "")), "declaration method hash"
    ):
        raise WorkflowValidationError(f"Reference shielding method hash differs: {model_id}")
    if normalize_sha256(str(declaration.get("structure_sha256", ""))) != normalize_sha256(
        str(row.get("start_xyz_sha256", ""))
    ):
        raise WorkflowValidationError(f"Reference structure hash differs: {model_id}")
    if normalize_sha256(str(declaration.get("opt_input_sha256", ""))) != normalize_sha256(
        str(row.get("opt_input_sha256", ""))
    ):
        raise WorkflowValidationError(f"Reference opt input declaration hash differs: {model_id}")
    if normalize_sha256(
        str(declaration.get("shielding_template_sha256", ""))
    ) != normalize_sha256(str(row.get("shielding_template_sha256", ""))):
        raise WorkflowValidationError(f"Reference shielding template declaration hash differs: {model_id}")
    if declaration.get("raw_shielding_and_processing_convention_separated") is not True:
        raise WorkflowValidationError(f"Reference raw/processing separation is not declared: {model_id}")
    raw_role = str(declaration.get("raw_calculation_reference", "")).strip()
    processing = str(declaration.get("experimental_processing_convention", "")).strip()
    if not raw_role or not processing or raw_role == processing:
        raise WorkflowValidationError(f"Reference raw and processing records are not distinct: {model_id}")
    aggregation = str(declaration.get("raw_shielding_aggregation_policy", ""))
    if not all(token in aggregation for token in ("PRESERVE_EACH_REQUESTED_ATOM_VALUE", "population standard deviation", "never")):
        raise WorkflowValidationError(f"Reference raw-shielding aggregation policy is incomplete: {model_id}")
    requested = {value for value in str(row.get("requested_nuclei", "")).split(";") if value}
    if "31P" in requested:
        if declaration.get("equivalent_to_85wt_percent_H3PO4_external_standard") is not False:
            raise WorkflowValidationError(f"31P model falsely claims 85 wt% H3PO4 equivalence: {model_id}")
        if "MODEL_ONLY_NOT_EQUIVALENT_TO_85WT_PERCENT_H3PO4" not in str(
            declaration.get("model_limitation", "")
        ):
            raise WorkflowValidationError(f"31P model lacks external-standard limitation: {model_id}")
    if gate_decision.get("model_id") != model_id:
        raise WorkflowValidationError(f"Reference gate decision identity differs: {model_id}")
    if ensure_within(
        Path(str(gate_decision.get("model_declaration_file", ""))), project, "gate declaration"
    ) != declaration_path.resolve():
        raise WorkflowValidationError(f"Reference gate declaration path differs: {model_id}")
    if sha256_file(declaration_path) != normalize_sha256(
        str(gate_decision.get("model_declaration_sha256", "")), "gate declaration hash"
    ):
        raise WorkflowValidationError(f"Reference gate declaration hash differs: {model_id}")
    map_evidence = _validate_atom_map_and_inputs(row, job_dir, declaration)
    return {
        "model_id": model_id,
        "job_directory": job_dir,
        "start_path": paths["start_xyz_sha256"],
        "start_xyz_sha256": sha256_file(paths["start_xyz_sha256"]),
        "map_path": paths["reference_atom_map_sha256"],
        "map_sha256": sha256_file(paths["reference_atom_map_sha256"]),
        "declaration_path": declaration_path,
        "declaration": declaration,
        "declaration_sha256": sha256_file(declaration_path),
        "opt_input_path": paths["opt_input_sha256"],
        "opt_input_sha256": sha256_file(paths["opt_input_sha256"]),
        "opt_input_bundle_sha256": inline_orca_input_bundle_sha256(
            paths["opt_input_sha256"].name, paths["opt_input_sha256"].read_bytes()
        ),
        "shielding_template_path": paths["shielding_template_sha256"],
        "shielding_template_sha256": sha256_file(paths["shielding_template_sha256"]),
        "method_fingerprint": fingerprint,
        **map_evidence,
    }


def _approval_check(
    gate: Mapping[str, object], rows: list[dict[str, str]], decisions: Mapping[str, Mapping[str, object]]
) -> set[str]:
    if gate.get("approval_status") != APPROVED:
        raise WorkflowValidationError(
            f"Reference execution requires protocol status {APPROVED}; current status is "
            f"{gate.get('approval_status')!r}"
        )
    for field in (
        "raw_shielding_reference_check",
        "experimental_processing_convention_check",
        "no_false_external_standard_equivalence_check",
    ):
        if gate.get(field) != "PASS":
            raise WorkflowValidationError(f"Reference protocol approval field is not PASS: {field}")
    if not str(gate.get("reviewer", "")).strip() or not str(gate.get("reviewed_at", "")).strip():
        raise WorkflowValidationError("Reference protocol approval lacks reviewer or reviewed_at")
    approved: set[str] = set()
    for row in rows:
        model_id = row["model_id"]
        decision = decisions[model_id]
        status = decision.get("approval_status")
        if status == APPROVED:
            if row.get("approval_status") != APPROVED:
                raise WorkflowValidationError(f"Reference manifest row is not approved: {model_id}")
            for field in (
                "structure_and_speciation_review",
                "reference_site_index_review",
                "method_match_review",
                "experimental_convention_review",
            ):
                if decision.get(field) != "PASS":
                    raise WorkflowValidationError(
                        f"Reference model review is not PASS ({field}): {model_id}"
                    )
            approved.add(model_id)
        elif status == REJECTED:
            if row.get("approval_status") != REJECTED:
                raise WorkflowValidationError(f"Rejected reference manifest row is unresolved: {model_id}")
            if not str(decision.get("decision_reason", "")).strip():
                raise WorkflowValidationError(f"Rejected reference model lacks a reason: {model_id}")
        else:
            raise WorkflowValidationError(f"Reference model decision is unresolved: {model_id}")
    route = gate.get("d2o_1h_13c_reference_route_decision")
    if not isinstance(route, dict) or route.get("status") != "PASS":
        raise WorkflowValidationError("D2O 1H/13C reference route is unresolved")
    selected_route = str(route.get("selected_route", ""))
    allowed = {str(value) for value in route.get("allowed_routes", [])}
    if selected_route not in allowed or not str(
        route.get("actual_standard_or_indirect_metadata", "")
    ).strip() or not str(route.get("decision_reason", "")).strip():
        raise WorkflowValidationError("D2O 1H/13C route lacks a valid selection or metadata")
    molecular_candidates = {"REF_D2O_DSS_NA_6D2O", "REF_D2O_TSP_NA_6D2O"}
    expected_approved = (
        {selected_route} if selected_route in molecular_candidates else set()
    )
    if approved & molecular_candidates != expected_approved:
        raise WorkflowValidationError(
            "Exactly the selected DSS/TSP model, or neither for an indirect route, may be approved"
        )
    if not approved:
        raise WorkflowValidationError("Reference protocol approves no molecular calculation")
    return approved


def validate_reference_protocol(
    project: Path,
    manifest_path: Path,
    gate_path: Path,
    *,
    require_approved: bool = False,
    only: Iterable[str] | None = None,
) -> tuple[dict, list[dict[str, str]], list[dict[str, object]]]:
    project = Path(project).resolve()
    manifest_path = ensure_within(Path(manifest_path), project, "reference manifest")
    gate_path = ensure_within(Path(gate_path), project, "reference protocol gate")
    rows = _read_manifest(manifest_path)
    gate = read_json(gate_path)
    if gate.get("notice") != NOTICE or gate.get("stage") != "reference_model_protocol":
        raise WorkflowValidationError("Reference protocol gate identity or notice differs")
    gate_manifest = ensure_within(
        Path(str(gate.get("reference_manifest_file", ""))), project, "gate reference manifest"
    )
    if gate_manifest != manifest_path.resolve():
        raise WorkflowValidationError("Reference protocol gate points to a different manifest")
    if sha256_file(manifest_path) != normalize_sha256(
        str(gate.get("reference_manifest_sha256", "")), "gate manifest hash"
    ):
        raise WorkflowValidationError("Reference protocol manifest hash differs")
    raw_decisions = gate.get("model_decisions")
    if not isinstance(raw_decisions, list) or not all(isinstance(value, dict) for value in raw_decisions):
        raise WorkflowValidationError("Reference protocol model decisions are missing")
    decisions = {str(value.get("model_id", "")): value for value in raw_decisions}
    if len(decisions) != len(raw_decisions) or set(decisions) != {row["model_id"] for row in rows}:
        raise WorkflowValidationError("Reference protocol decisions do not match manifest models")
    selected = set(only or [])
    unknown = selected - {row["model_id"] for row in rows}
    if unknown:
        raise WorkflowValidationError(f"Unknown reference model selection: {sorted(unknown)}")
    all_evidence = [
        validate_reference_manifest_job(row, project, decisions[row["model_id"]])
        for row in rows
    ]
    if require_approved:
        approved = _approval_check(gate, rows, decisions)
        if selected and not selected <= approved:
            raise WorkflowValidationError(
                f"Requested reference model is not approved: {sorted(selected - approved)}"
            )
        selected_ids = selected or approved
        for item in all_evidence:
            decision_status = decisions[str(item["model_id"])].get("approval_status")
            if item["declaration"].get("approval_status") != decision_status:
                raise WorkflowValidationError(
                    f"Reference declaration status differs from gate: {item['model_id']}"
                )
    else:
        selected_ids = selected or {row["model_id"] for row in rows}
    selected_rows = [row for row in rows if row["model_id"] in selected_ids]
    evidence_by_id = {str(item["model_id"]): item for item in all_evidence}
    evidence = [evidence_by_id[row["model_id"]] for row in selected_rows]
    if not selected_rows:
        raise WorkflowValidationError("No reference jobs selected")
    return gate, selected_rows, evidence


def dry_run_validate(
    project: Path,
    manifest_path: Path,
    gate_path: Path,
    only: Iterable[str] | None = None,
) -> list[dict[str, object]]:
    return validate_reference_protocol(
        project, manifest_path, gate_path, require_approved=False, only=only
    )[2]


def _resource_gate(project: Path, policy: Mapping[str, object], stage_name: str) -> dict[str, object]:
    execution = policy["execution"]
    stage = policy["stage_resources"][stage_name]
    start = policy["start_gates"]
    ranks = int(stage["mpi_ranks"])
    if ranks > int(execution["absolute_mpi_rank_cap"]):
        raise WorkflowValidationError(f"Reference {stage_name} MPI ranks exceed absolute cap")
    blocked_names = {
        str(value).casefold() for value in start["block_if_external_process_names_match"]
    }
    blocked = sorted(running_process_names() & blocked_names)
    if blocked:
        raise WorkflowValidationError(f"External quantum-chemistry process is running: {blocked}")
    sample_seconds = float(start["cpu_sample_seconds"])
    cpu = cpu_percent(sample_seconds)
    ram = available_ram_gib()
    free_disk = shutil.disk_usage(project.anchor).free / (1024**3)
    if cpu > float(start["maximum_average_cpu_percent"]):
        raise WorkflowValidationError(f"Pre-start CPU {cpu:.1f}% exceeds gate")
    if ram < float(stage["minimum_free_ram_gib"]):
        raise WorkflowValidationError(
            f"Free RAM {ram:.1f} GiB is below reference {stage_name} gate"
        )
    required_disk = float(start["minimum_free_space_gib"].get(project.drive.rstrip(":"), 0))
    if free_disk < required_disk:
        raise WorkflowValidationError(f"Free disk {free_disk:.1f} GiB is below gate")
    return {
        "stage": stage_name,
        "cpu_average_percent": cpu,
        "cpu_sample_seconds": sample_seconds,
        "free_ram_gib": ram,
        "free_project_drive_gib": free_disk,
        "mpi_ranks": ranks,
        "maxcore_mb_per_rank": int(stage["maxcore_mb_per_rank"]),
    }


def _orca_environment(work: Path, policy: Mapping[str, object], orca_exe: Path) -> dict[str, str]:
    environment = os.environ.copy()
    environment.update({"OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
    temporary = work / "_tmp"
    temporary.mkdir()
    environment.update({"TEMP": str(temporary), "TMP": str(temporary)})
    mpi_bin = str(policy.get("machine", {}).get("msmpi_bin", "")).strip()
    prefixes = [str(orca_exe.parent)] + ([mpi_bin] if mpi_bin else [])
    environment["PATH"] = os.pathsep.join(prefixes) + os.pathsep + environment.get("PATH", "")
    return environment


def _artifact(path: Path) -> dict[str, object]:
    exists = path.is_file() and path.stat().st_size > 0
    return {
        "path": str(path.resolve()),
        "exists_nonempty": exists,
        "size_bytes": path.stat().st_size if exists else 0,
        "sha256": sha256_file(path) if exists else "",
    }


def _run_stage(
    *,
    model_id: str,
    stage_name: str,
    source_input: Path,
    work: Path,
    policy: Mapping[str, object],
    project: Path,
    orca_exe: Path,
    coordinate_source: Path | None = None,
) -> tuple[dict[str, object], bool]:
    gate_snapshot = _resource_gate(project, policy, stage_name)
    work.mkdir(parents=True, exist_ok=False)
    input_name = "opt_freq.inp" if stage_name == "opt_freq" else "shielding.inp"
    input_path = work / input_name
    shutil.copy2(source_input, input_path)
    if sha256_file(input_path) != sha256_file(source_input):
        raise WorkflowValidationError(f"Copied reference input differs: {model_id}/{stage_name}")
    coordinate_path: Path | None = None
    if coordinate_source is not None:
        coordinate_path = work / "optimized.xyz"
        shutil.copy2(coordinate_source, coordinate_path)
        if sha256_file(coordinate_path) != sha256_file(coordinate_source):
            raise WorkflowValidationError(f"Copied optimized reference geometry differs: {model_id}")
    input_bytes = input_path.read_bytes()
    bundle_hash = (
        inline_orca_input_bundle_sha256(input_path.name, input_bytes)
        if coordinate_path is None
        else orca_input_bundle_sha256(
            input_path.name, input_bytes, (("optimized.xyz", coordinate_path.read_bytes()),)
        )
    )
    output_path = work / f"{input_path.stem}.out"
    stderr_path = work / f"{input_path.stem}.stderr.log"
    record: dict[str, object] = {
        "stage": stage_name,
        "status": "running",
        "started_utc": utc_now(),
        "source_input_file": str(source_input.resolve()),
        "source_input_sha256": sha256_file(source_input),
        "run_input_file": str(input_path.resolve()),
        "run_input_sha256": sha256_file(input_path),
        "input_bundle_sha256": bundle_hash,
        "coordinate_source_file": str(coordinate_source.resolve()) if coordinate_source else "",
        "coordinate_source_sha256": sha256_file(coordinate_source) if coordinate_source else "",
        "run_coordinate_file": str(coordinate_path.resolve()) if coordinate_path else "",
        "run_coordinate_sha256": sha256_file(coordinate_path) if coordinate_path else "",
        "prestart_resource_gate": gate_snapshot,
        "command": [str(orca_exe.resolve()), input_path.name],
        "orca_executable_sha256": sha256_file(orca_exe),
        "output_file": str(output_path.resolve()),
        "stderr_file": str(stderr_path.resolve()),
    }
    environment = _orca_environment(work, policy, orca_exe)
    creationflags = getattr(subprocess, "CREATE_NO_WINDOW", 0) | getattr(
        subprocess, "BELOW_NORMAL_PRIORITY_CLASS", 0
    )
    started = time.monotonic()
    stop_after_current = False
    sustained = 0.0
    heartbeat = float(policy["runtime_policy"]["heartbeat_seconds"])
    high_threshold = float(policy["runtime_policy"]["high_cpu_percent"])
    sustained_limit = float(policy["runtime_policy"]["high_cpu_sustained_seconds"])
    with output_path.open("wb") as stdout_handle, stderr_path.open("wb") as stderr_handle:
        process = subprocess.Popen(
            [str(orca_exe.resolve()), input_path.name],
            cwd=work,
            env=environment,
            stdout=stdout_handle,
            stderr=stderr_handle,
            creationflags=creationflags,
        )
        while process.poll() is None:
            sampled = cpu_percent(heartbeat)
            sustained = sustained + heartbeat if sampled >= high_threshold else 0.0
            if sustained >= sustained_limit:
                stop_after_current = True
        return_code = int(process.returncode)
    normal = output_path.is_file() and NORMAL_MARKER in output_path.read_text(
        encoding="utf-8", errors="replace"
    )
    stem = input_path.stem
    expected = (
        {
            "gbw": work / f"{stem}.gbw",
            "hess": work / f"{stem}.hess",
            "property": work / f"{stem}.property.txt",
            "optimized_xyz": work / f"{stem}.xyz",
        }
        if stage_name == "opt_freq"
        else {
            "gbw": work / f"{stem}.gbw",
            "property": work / f"{stem}.property.txt",
        }
    )
    artifacts = {role: _artifact(path) for role, path in expected.items()}
    source_unchanged = sha256_file(source_input) == record["source_input_sha256"]
    coordinate_unchanged = (
        coordinate_source is None
        or (
            sha256_file(coordinate_source) == record["coordinate_source_sha256"]
            and coordinate_path is not None
            and sha256_file(coordinate_path) == record["run_coordinate_sha256"]
        )
    )
    complete = (
        return_code == 0
        and normal
        and source_unchanged
        and coordinate_unchanged
        and all(bool(value["exists_nonempty"]) for value in artifacts.values())
    )
    record.update(
        {
            "status": "completed" if complete else "failed",
            "ended_utc": utc_now(),
            "elapsed_seconds": round(time.monotonic() - started, 2),
            "return_code": return_code,
            "normal_termination_marker": normal,
            "source_input_unchanged": source_unchanged,
            "coordinate_source_unchanged": coordinate_unchanged,
            "output_sha256": sha256_file(output_path) if output_path.is_file() else "",
            "artifacts": artifacts,
            "stop_after_current_due_to_sustained_high_cpu": stop_after_current,
        }
    )
    return record, stop_after_current


def _validate_optimized_reference_geometry(path: Path, evidence: Mapping[str, object]) -> None:
    optimized = parse_xyz_frames(path)
    source = parse_xyz_frames(Path(evidence["start_path"]))
    if len(optimized) != 1 or len(source) != 1:
        raise WorkflowValidationError(f"Reference optimized XYZ frame count differs: {evidence['model_id']}")
    if optimized[0].symbols != source[0].symbols:
        raise WorkflowValidationError(f"Reference optimized XYZ atom order differs: {evidence['model_id']}")


def _validate_optimum(output_path: Path) -> dict[str, object]:
    thermo = parse_orca_qrrho(output_path)
    if int(thermo["imaginary_frequency_count"]) != 0:
        raise WorkflowValidationError(
            f"Reference optimization has imaginary frequencies: {output_path}"
        )
    if not math.isclose(float(thermo["qrrho_cutoff_cm1"]), 100.0, rel_tol=0.0, abs_tol=0.01):
        raise WorkflowValidationError(f"Reference QRRHO cutoff is not 100.0 cm-1: {output_path}")
    if not math.isclose(float(thermo["temperature_K"]), 298.15, rel_tol=0.0, abs_tol=0.02):
        raise WorkflowValidationError(f"Reference thermochemistry temperature is not 298.15 K: {output_path}")
    return thermo


def _current_sources(evidence: Mapping[str, object]) -> dict[str, str]:
    return {
        "start_xyz_sha256": sha256_file(Path(evidence["start_path"])),
        "reference_atom_map_sha256": sha256_file(Path(evidence["map_path"])),
        "model_declaration_sha256": sha256_file(Path(evidence["declaration_path"])),
        "opt_input_sha256": sha256_file(Path(evidence["opt_input_path"])),
        "shielding_template_sha256": sha256_file(Path(evidence["shielding_template_path"])),
    }


def _protocol_sources_unchanged(
    project: Path,
    protocol_gate: Path,
    expected_gate_hash: str,
    expected_source_hashes: Mapping[str, str],
    evidence: Mapping[str, object],
) -> bool:
    if sha256_file(protocol_gate) != expected_gate_hash:
        return False
    gate = read_json(protocol_gate)
    manifest = ensure_within(
        Path(str(gate.get("reference_manifest_file", ""))), project, "reference manifest"
    )
    if sha256_file(manifest) != normalize_sha256(
        str(gate.get("reference_manifest_sha256", "")), "reference manifest hash"
    ):
        return False
    return _current_sources(evidence) == dict(expected_source_hashes)


def _completed_state_valid(
    state_path: Path, evidence: Mapping[str, object], gate_hash: str, project: Path
) -> bool:
    if not state_path.is_file():
        return False
    try:
        state = read_json(state_path)
        if state.get("status") != "completed_calculated_reference_unapproved_for_processing":
            return False
        if state.get("model_id") != evidence["model_id"] or state.get("protocol_gate_sha256") != gate_hash:
            return False
        if state.get("source_hashes") != _current_sources(evidence):
            return False
        for stage in ("opt_freq", "shielding"):
            record = state["stages"][stage]
            if (
                record.get("status") != "completed"
                or record.get("normal_termination_marker") is not True
                or record.get("source_input_unchanged") is not True
                or record.get("coordinate_source_unchanged") is not True
            ):
                return False
            output = ensure_within(Path(record["output_file"]), project, "reference output")
            if sha256_file(output) != normalize_sha256(str(record["output_sha256"])):
                return False
            if NORMAL_MARKER not in output.read_text(encoding="utf-8", errors="replace"):
                return False
            run_input = ensure_within(
                Path(record["run_input_file"]), project, "reference run input"
            )
            source_input = ensure_within(
                Path(record["source_input_file"]), project, "reference source input"
            )
            if sha256_file(run_input) != normalize_sha256(str(record["run_input_sha256"])):
                return False
            if sha256_file(source_input) != normalize_sha256(
                str(record["source_input_sha256"])
            ):
                return False
            if stage == "opt_freq":
                current_bundle = inline_orca_input_bundle_sha256(
                    run_input.name, run_input.read_bytes()
                )
                if record.get("minimum_and_atom_order_gate") != "PASS":
                    return False
                _validate_optimum(output)
            else:
                coordinate = ensure_within(
                    Path(record["coordinate_source_file"]), project, "reference optimized geometry"
                )
                run_coordinate = ensure_within(
                    Path(record["run_coordinate_file"]), project, "reference run geometry"
                )
                if sha256_file(coordinate) != normalize_sha256(
                    str(record["coordinate_source_sha256"])
                ) or sha256_file(run_coordinate) != normalize_sha256(
                    str(record["run_coordinate_sha256"])
                ) or sha256_file(run_coordinate) != sha256_file(coordinate):
                    return False
                current_bundle = orca_input_bundle_sha256(
                    run_input.name,
                    run_input.read_bytes(),
                    (("optimized.xyz", run_coordinate.read_bytes()),),
                )
            if current_bundle != normalize_sha256(str(record["input_bundle_sha256"])):
                return False
            for artifact in record.get("artifacts", {}).values():
                path = ensure_within(Path(artifact["path"]), project, "reference artifact")
                if not path.is_file() or sha256_file(path) != normalize_sha256(str(artifact["sha256"])):
                    return False
        return True
    except (KeyError, OSError, ValueError, WorkflowValidationError):
        return False


def run_one_reference(
    row: Mapping[str, str],
    evidence: Mapping[str, object],
    *,
    project: Path,
    protocol_gate: Path,
    state_root: Path,
    run_root: Path,
    policy: Mapping[str, object],
    orca_exe: Path,
) -> tuple[dict[str, object], bool]:
    model_id = str(evidence["model_id"])
    state_path = state_root / f"{model_id}.json"
    gate_hash = sha256_file(protocol_gate)
    if _completed_state_valid(state_path, evidence, gate_hash, project):
        return read_json(state_path), False
    attempt_id = f"attempt_{datetime.now().strftime('%Y%m%dT%H%M%S')}_{uuid.uuid4().hex[:8]}"
    attempt_root = run_root / model_id / attempt_id
    state: dict[str, object] = {
        "schema_version": "2.0",
        "notice": NOTICE,
        "model_id": model_id,
        "status": "running",
        "raw_shielding_role": "CALCULATED_REFERENCE_MODEL_NOT_EXPERIMENTAL_STANDARD",
        "experimental_processing_convention_applied": False,
        "protocol_gate_file": str(protocol_gate.resolve()),
        "protocol_gate_sha256": gate_hash,
        "resource_policy_sha256": sha256_file(project / "config" / "resource_policy.json"),
        "source_hashes": _current_sources(evidence),
        "attempt_id": attempt_id,
        "attempt_root": str(attempt_root.resolve()),
        "stages": {},
    }
    state_root.mkdir(parents=True, exist_ok=True)
    atomic_write_json(state_path, state)
    if not _protocol_sources_unchanged(
        project, protocol_gate, gate_hash, state["source_hashes"], evidence
    ):
        state["status"] = "blocked_source_changed_before_opt_freq"
        atomic_write_json(state_path, state)
        return state, True
    opt_record, stop_after_opt = _run_stage(
        model_id=model_id,
        stage_name="opt_freq",
        source_input=Path(evidence["opt_input_path"]),
        work=attempt_root / "opt_freq" / "work",
        policy=policy,
        project=project,
        orca_exe=orca_exe,
    )
    state["stages"]["opt_freq"] = opt_record
    if opt_record["status"] != "completed":
        state["status"] = "failed_opt_freq"
        atomic_write_json(state_path, state)
        return state, True
    optimized_path = Path(opt_record["artifacts"]["optimized_xyz"]["path"])
    try:
        _validate_optimized_reference_geometry(optimized_path, evidence)
        opt_record["qrrho_minimum_evidence"] = _validate_optimum(
            Path(opt_record["output_file"])
        )
        opt_record["minimum_and_atom_order_gate"] = "PASS"
    except (OSError, ValueError, WorkflowValidationError) as exc:
        opt_record["minimum_and_atom_order_gate"] = "FAIL"
        opt_record["minimum_and_atom_order_failure"] = str(exc)
        state["status"] = "failed_minimum_or_atom_order_gate"
        atomic_write_json(state_path, state)
        return state, True
    atomic_write_json(state_path, state)
    if stop_after_opt:
        state["status"] = "stopped_after_opt_freq_due_to_sustained_high_cpu"
        atomic_write_json(state_path, state)
        return state, True
    if not _protocol_sources_unchanged(
        project, protocol_gate, gate_hash, state["source_hashes"], evidence
    ):
        state["status"] = "blocked_source_changed_before_shielding"
        atomic_write_json(state_path, state)
        return state, True
    shielding_record, stop_after_shielding = _run_stage(
        model_id=model_id,
        stage_name="shielding",
        source_input=Path(evidence["shielding_template_path"]),
        coordinate_source=optimized_path,
        work=attempt_root / "shielding" / "work",
        policy=policy,
        project=project,
        orca_exe=orca_exe,
    )
    state["stages"]["shielding"] = shielding_record
    if shielding_record["status"] != "completed":
        state["status"] = "failed_shielding"
        atomic_write_json(state_path, state)
        return state, True
    state["status"] = "completed_calculated_reference_unapproved_for_processing"
    state["completed_utc"] = utc_now()
    state["next_required_action"] = (
        "Extract and independently review raw shieldings, then bind exact input/output hashes in "
        "an approved reference registry; do not infer an experimental reference convention."
    )
    atomic_write_json(state_path, state)
    return state, stop_after_shielding


def run_reference_serial(
    *,
    project: Path,
    manifest_path: Path,
    gate_path: Path,
    execute: bool = False,
    only: Iterable[str] | None = None,
    orca_exe: Path | None = None,
    state_root: Path | None = None,
    run_root: Path | None = None,
) -> list[dict[str, object]]:
    project = Path(project).resolve()
    gate, rows, evidence = validate_reference_protocol(
        project,
        manifest_path,
        gate_path,
        require_approved=execute,
        only=only,
    )
    if not execute:
        return evidence
    policy_path = project / "config" / "resource_policy.json"
    policy = read_json(policy_path)
    executable = Path(
        orca_exe or Path(str(policy.get("machine", {}).get("orca_executable", "")))
    ).resolve()
    if not executable.is_file():
        raise WorkflowValidationError(f"ORCA executable is missing: {executable}")
    if sha256_file(executable) != normalize_sha256(
        str(policy.get("machine", {}).get("orca_sha256", "")), "ORCA executable hash"
    ):
        raise WorkflowValidationError("ORCA executable hash differs from resource policy")
    states = (state_root or project / "logs" / "reference_serial" / "jobs").resolve()
    runs = (run_root or project / "reference_runs").resolve()
    results: list[dict[str, object]] = []
    with ExclusiveLock(project / "logs" / "reference_serial" / "RUNNER.lock"):
        for row, item in zip(rows, evidence, strict=True):
            state, stop = run_one_reference(
                row,
                item,
                project=project,
                protocol_gate=Path(gate_path),
                state_root=states,
                run_root=runs,
                policy=policy,
                orca_exe=executable,
            )
            results.append(state)
            if stop or not str(state["status"]).startswith("completed_"):
                break
    return results


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--project", type=Path, default=PROJECT)
    value.add_argument("--manifest", type=Path)
    value.add_argument("--gate", type=Path)
    value.add_argument("--only", nargs="*")
    value.add_argument("--execute", action="store_true")
    value.add_argument("--orca-exe", type=Path)
    value.add_argument("--state-root", type=Path)
    value.add_argument("--run-root", type=Path)
    return value


def main() -> None:
    args = parser().parse_args()
    project = args.project.resolve()
    manifest = (args.manifest or project / "data" / "reference_job_manifest_DRAFT.csv").resolve()
    gate = (args.gate or project / "jobs" / "reference_model_protocol_DRAFT.json").resolve()
    results = run_reference_serial(
        project=project,
        manifest_path=manifest,
        gate_path=gate,
        execute=args.execute,
        only=args.only,
        orca_exe=args.orca_exe,
        state_root=args.state_root,
        run_root=args.run_root,
    )
    if not args.execute:
        gate_status = read_json(gate).get("approval_status")
        print(
            f"DRY RUN: validated {len(results)} reference models; gate={gate_status}; "
            "ORCA was not started."
        )
    else:
        for state in results:
            print(f"{state['model_id']}: {state['status']}")


if __name__ == "__main__":
    main()
