"""Build final PBE0-D4/def2-TZVPP opt/freq jobs from approved GOAT gates.

This program is intentionally a *builder*, not a runner.  It never launches
ORCA.  A production job can be written only after byte-level provenance,
scenario isolation, reviewed frame extraction, atom ordering and isotope maps
all pass a fail-closed validation.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from workflow_schema import (
    HARTREE_TO_KCAL_MOL,
    NOTICE,
    WorkflowValidationError,
    atomic_write_bytes,
    atomic_write_csv,
    atomic_write_json,
    atomic_write_text,
    downstream_atom_map_rows,
    ensure_within,
    frames_numerically_equal,
    inline_orca_input_bundle_sha256,
    normalize_sha256,
    parse_xyz_frames,
    read_csv,
    read_json,
    render_opt_freq_input,
    scientific_dependency_sha256,
    sha256_bytes,
    sha256_file,
    validate_cluster_maps,
    validate_completed_goat_state,
)


PROJECT = Path(__file__).resolve().parents[1]
APPROVED = "APPROVED_FOR_OPT_FREQ"
PASS = "PASS"
STARTED_PATTERNS = (
    "*.out",
    "*.gbw",
    "*.hess",
    "*.engrad",
    "*.opt",
    "*.property.txt",
    "*property*.txt",
    "*.tmp",
)


@dataclass
class PlannedJob:
    job_id: str
    directory: Path
    files: dict[str, bytes]
    manifest_row: dict[str, object]


def safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    if not token or token in {".", ".."}:
        raise WorkflowValidationError(f"Unsafe identifier: {value!r}")
    return token


def _manifest_index(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    index: dict[str, dict[str, str]] = {}
    for row in rows:
        job_id = row.get("job_id", "").strip()
        if not job_id or job_id in index:
            raise WorkflowValidationError(f"Blank or duplicate GOAT job_id: {job_id!r}")
        index[job_id] = row
    return index


def _registry_index(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    index: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        key = (row.get("figure_id", "").strip(), row.get("scenario_id", "").strip())
        if not all(key) or key in index:
            raise WorkflowValidationError(f"Blank or duplicate scenario key: {key}")
        index[key] = row
    return index


def _resource_values(policy_path: Path) -> tuple[int, int, str]:
    policy = read_json(policy_path)
    execution = policy.get("execution", {})
    stage = policy.get("stage_resources", {}).get("opt_freq", {})
    if not isinstance(execution, dict) or not isinstance(stage, dict):
        raise WorkflowValidationError(f"Malformed opt/freq resource policy: {policy_path}")
    try:
        nprocs = int(stage["mpi_ranks"])
        maxcore = int(stage["maxcore_mb_per_rank"])
        cap = int(execution["absolute_mpi_rank_cap"])
    except (KeyError, TypeError, ValueError) as exc:
        raise WorkflowValidationError(f"Incomplete opt/freq resource policy: {policy_path}") from exc
    legacy_serial = execution.get("one_compound_at_a_time") is True
    hpc_array = (
        execution.get("execution_model") == "one_orca_calculation_per_scheduler_array_task"
        and execution.get("one_orca_calculation_per_array_task") is True
    )
    if nprocs > cap or cap > 32 or not (legacy_serial or hpc_array):
        raise WorkflowValidationError(
            "Resource policy must use either legacy single-instance execution or one ORCA "
            "calculation per scheduler-array task, with <=32 MPI ranks"
        )
    return nprocs, maxcore, sha256_file(policy_path)


def _calibration_references(project: Path) -> list[dict[str, str]]:
    root = project / "data" / "real_nmr_calibration_audit"
    if not root.is_dir():
        return []
    references: list[dict[str, str]] = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        references.append(
            {
                "file": str(path.resolve()),
                "project_relative_file": path.relative_to(project).as_posix(),
                "sha256": sha256_file(path),
                "role": "DOWNSTREAM_LINE_SHAPE_NOISE_AND_REFERENCE_CALIBRATION_ONLY",
            }
        )
    return references


def _locked_acquisition_temperatures(
    rows: list[dict[str, str]], figure_id: str
) -> dict[str, float]:
    """Return the reviewed/explicitly-assumed temperature for each NMR nucleus."""
    result: dict[str, float] = {}
    for row in rows:
        if row.get("figure_id", "").strip() != figure_id:
            continue
        nucleus = row.get("nucleus", "").strip()
        if nucleus not in {"1H", "13C", "19F", "31P"}:
            continue
        if nucleus in result:
            raise WorkflowValidationError(f"Duplicate acquisition condition for {figure_id}/{nucleus}")
        status = row.get("condition_status", "").strip()
        try:
            temperature = float(row["temperature_K"])
        except (KeyError, TypeError, ValueError) as exc:
            raise WorkflowValidationError(
                f"Invalid acquisition temperature for {figure_id}/{nucleus}"
            ) from exc
        if temperature <= 0.0 or not status or status.upper() == "UNKNOWN":
            raise WorkflowValidationError(f"Unlocked acquisition condition for {figure_id}/{nucleus}")
        result[nucleus] = temperature
    if not result:
        raise WorkflowValidationError(f"No locked nucleus acquisition conditions for {figure_id}")
    return result


def _assert_approved_gate_shape(gate: Mapping[str, object], gate_path: Path) -> None:
    if gate.get("approval_status") != APPROVED:
        raise WorkflowValidationError(f"Gate is not {APPROVED}: {gate_path}")
    if gate.get("selection_completeness_check") != PASS:
        raise WorkflowValidationError(f"Selection completeness is not PASS: {gate_path}")
    if gate.get("protonation_and_cluster_integrity_review") != PASS:
        raise WorkflowValidationError(f"Cluster-integrity review is not PASS: {gate_path}")
    if not str(gate.get("reviewer", "")).strip() or not str(gate.get("reviewed_at", "")).strip():
        raise WorkflowValidationError(f"Approved gate lacks reviewer/reviewed_at: {gate_path}")
    selections = gate.get("selected_conformers")
    if not isinstance(selections, list) or not selections:
        raise WorkflowValidationError(f"Approved gate has no selected conformers: {gate_path}")
    for selection in selections:
        if not isinstance(selection, dict):
            raise WorkflowValidationError(f"Malformed selected conformer in {gate_path}")
        if selection.get("approval_status") != APPROVED:
            raise WorkflowValidationError(f"Selected conformer is not approved in {gate_path}")
        for field in (
            "connectivity_stereochemistry_check",
            "cluster_composition_check",
            "manual_extraction_check",
        ):
            if selection.get(field) != PASS:
                raise WorkflowValidationError(f"{field} is not PASS in {gate_path}")


def _resolve_state_records(
    project: Path,
    gate: Mapping[str, object],
    gate_path: Path,
    figure_rows: list[dict[str, str]],
) -> dict[str, dict[str, object]]:
    records = gate.get("upstream_state_hashes")
    if not isinstance(records, list):
        raise WorkflowValidationError(f"Gate lacks upstream_state_hashes: {gate_path}")
    expected_job_ids = {row["job_id"] for row in figure_rows}
    supplied: dict[str, dict[str, object]] = {}
    manifest_by_job = {row["job_id"]: row for row in figure_rows}
    for record in records:
        if not isinstance(record, dict):
            raise WorkflowValidationError(f"Malformed upstream state record: {gate_path}")
        job_id = str(record.get("job_id", "")).strip()
        if not job_id or job_id in supplied or job_id not in expected_job_ids:
            raise WorkflowValidationError(f"Unknown/duplicate upstream state job_id {job_id!r}")
        state_path = ensure_within(
            Path(str(record.get("state_file", ""))), project, "approved-gate GOAT state"
        )
        expected_state_hash = normalize_sha256(str(record.get("sha256", "")), "gate state sha256")
        if sha256_file(state_path) != expected_state_hash:
            raise WorkflowValidationError(f"Gate state hash mismatch: {state_path}")
        evidence = validate_completed_goat_state(manifest_by_job[job_id], state_path, project)
        if evidence["state_sha256"] != expected_state_hash:
            raise WorkflowValidationError(f"Validated state hash changed: {state_path}")
        supplied[job_id] = evidence
    if set(supplied) != expected_job_ids:
        missing = sorted(expected_job_ids - set(supplied))
        extra = sorted(set(supplied) - expected_job_ids)
        raise WorkflowValidationError(
            f"Approved gate must bind every manifest start exactly once; missing={missing}, extra={extra}"
        )
    return supplied


def _safe_existing_directory(directory: Path, replace_unstarted: bool) -> None:
    if not directory.exists():
        return
    started = [path for pattern in STARTED_PATTERNS for path in directory.rglob(pattern)]
    if started:
        raise WorkflowValidationError(
            f"Refusing to alter a DFT directory containing calculation artifacts: {started[0]}"
        )
    if not replace_unstarted:
        raise WorkflowValidationError(
            f"DFT directory already exists; use --replace-unstarted only after inspection: {directory}"
        )


def _plan_gate_jobs(
    project: Path,
    gate_path: Path,
    manifest_rows: list[dict[str, str]],
    manifest_by_job: dict[str, dict[str, str]],
    registry: dict[tuple[str, str], dict[str, str]],
    output_root: Path,
    nprocs: int,
    maxcore_mb: int,
    resource_policy_hash: str,
    manifest_hash: str,
    registry_hash: str,
    calibration_references: list[dict[str, str]],
    nucleus_condition_rows: list[dict[str, str]],
    nucleus_conditions_hash: str,
) -> list[PlannedJob]:
    gate = read_json(gate_path)
    _assert_approved_gate_shape(gate, gate_path)
    gate_hash = sha256_file(gate_path)
    figure_id = str(gate.get("figure_id", "")).strip()
    if not figure_id:
        raise WorkflowValidationError(f"Approved gate has blank figure_id: {gate_path}")
    temperatures_by_nucleus = _locked_acquisition_temperatures(
        nucleus_condition_rows, figure_id
    )
    thermochemistry_temperatures = sorted(set(temperatures_by_nucleus.values()))
    figure_rows = [row for row in manifest_rows if row["figure_id"] == figure_id]
    if not figure_rows:
        raise WorkflowValidationError(f"No manifest rows for gate figure {figure_id}")
    states = _resolve_state_records(project, gate, gate_path, figure_rows)
    expected_scenarios = {row["scenario_id"] for row in figure_rows}
    selections = list(gate["selected_conformers"])
    supplied_scenarios = {str(item.get("condition_scenario_id", "")) for item in selections}
    if supplied_scenarios != expected_scenarios:
        raise WorkflowValidationError(
            f"Approved gate must retain >=1 conformer for each scenario; expected={sorted(expected_scenarios)}, "
            f"supplied={sorted(supplied_scenarios)}"
        )
    try:
        energy_window = float(gate["selection_energy_window_kcal_mol"])
    except (KeyError, TypeError, ValueError) as exc:
        raise WorkflowValidationError(f"Approved gate has invalid energy window: {gate_path}") from exc
    if energy_window <= 0:
        raise WorkflowValidationError(f"Approved gate energy window must be positive: {gate_path}")
    scenario_minimum_hartree: dict[str, float] = {}
    for row in figure_rows:
        frames = parse_xyz_frames(Path(states[row["job_id"]]["ensemble_path"]))
        energies = [
            float(frame.energy_hartree)
            for frame in frames
            if frame.converged is True and frame.energy_hartree is not None
        ]
        if not energies:
            continue
        scenario_id = row["scenario_id"]
        minimum = min(energies)
        scenario_minimum_hartree[scenario_id] = min(
            minimum, scenario_minimum_hartree.get(scenario_id, minimum)
        )
    if set(scenario_minimum_hartree) != expected_scenarios:
        raise WorkflowValidationError(
            f"At least one scenario lacks a converged energy-labelled GOAT frame: {gate_path}"
        )

    ordered = sorted(
        selections,
        key=lambda item: (
            str(item["condition_scenario_id"]),
            float(item.get("relative_energy_kcal_mol", 0.0)),
            str(item.get("conformer_id", "")),
            str(item.get("source_goat_job_id", "")),
        ),
    )
    scenario_ordinals: dict[str, int] = defaultdict(int)
    planned: list[PlannedJob] = []
    seen_selection_keys: set[tuple[str, str, int]] = set()
    for selection in ordered:
        scenario_id = str(selection["condition_scenario_id"])
        source_job_id = str(selection.get("source_goat_job_id", "")).strip()
        if source_job_id not in manifest_by_job or source_job_id not in states:
            raise WorkflowValidationError(f"Selected source GOAT job is not bound by gate: {source_job_id}")
        source_row = manifest_by_job[source_job_id]
        if source_row["figure_id"] != figure_id or source_row["scenario_id"] != scenario_id:
            raise WorkflowValidationError(
                f"Selection crosses figure/scenario boundary: {figure_id}/{scenario_id} <- {source_job_id}"
            )
        scenario_key = (figure_id, scenario_id)
        if scenario_key not in registry:
            raise WorkflowValidationError(f"Scenario missing from registry: {scenario_key}")
        scenario = registry[scenario_key]
        evidence = states[source_job_id]

        state_path = Path(str(selection.get("source_state_file", ""))).resolve()
        if state_path != Path(evidence["state_path"]).resolve():
            raise WorkflowValidationError(f"Selection state path mismatch for {source_job_id}")
        if normalize_sha256(str(selection.get("source_state_sha256", ""))) != evidence["state_sha256"]:
            raise WorkflowValidationError(f"Selection state hash mismatch for {source_job_id}")
        artifact_path = Path(str(selection.get("source_artifact_file", ""))).resolve()
        if artifact_path != Path(evidence["ensemble_path"]).resolve():
            raise WorkflowValidationError(f"Selection ensemble path mismatch for {source_job_id}")
        artifact_hash = normalize_sha256(str(selection.get("source_artifact_sha256", "")))
        if artifact_hash != evidence["ensemble_sha256"] or sha256_file(artifact_path) != artifact_hash:
            raise WorkflowValidationError(f"Selection ensemble hash mismatch for {source_job_id}")
        try:
            frame_index = int(selection["source_frame_index_1based"])
        except (KeyError, TypeError, ValueError) as exc:
            raise WorkflowValidationError(f"Invalid source frame index for {source_job_id}") from exc
        selection_key = (scenario_id, source_job_id, frame_index)
        if selection_key in seen_selection_keys:
            raise WorkflowValidationError(f"Duplicate selected source frame: {selection_key}")
        seen_selection_keys.add(selection_key)
        source_frames = parse_xyz_frames(artifact_path)
        if frame_index < 1 or frame_index > len(source_frames):
            raise WorkflowValidationError(f"Selected frame index is outside ensemble: {selection_key}")
        source_frame = source_frames[frame_index - 1]
        if source_frame.converged is not True or source_frame.energy_hartree is None:
            raise WorkflowValidationError(f"Approved source frame is not converged/energy-labelled: {selection_key}")
        try:
            recorded_energy = float(selection["source_frame_energy_hartree"])
            recorded_relative = float(selection["relative_energy_kcal_mol"])
        except (KeyError, TypeError, ValueError) as exc:
            raise WorkflowValidationError(f"Approved selection lacks numeric energy provenance: {selection_key}") from exc
        recalculated_relative = (
            float(source_frame.energy_hartree) - scenario_minimum_hartree[scenario_id]
        ) * HARTREE_TO_KCAL_MOL
        if abs(recorded_energy - float(source_frame.energy_hartree)) > 1e-9:
            raise WorkflowValidationError(f"Recorded source-frame energy is inconsistent: {selection_key}")
        if abs(recorded_relative - recalculated_relative) > 1e-5:
            raise WorkflowValidationError(f"Recorded relative energy is inconsistent: {selection_key}")
        if recalculated_relative > energy_window + 1e-8:
            raise WorkflowValidationError(
                f"Approved frame lies outside the reviewed {energy_window:g} kcal/mol window: {selection_key}"
            )

        selected_path = ensure_within(
            Path(str(selection.get("selected_coordinate_file", ""))),
            project,
            "approved selected coordinate",
        )
        selected_hash = normalize_sha256(str(selection.get("selected_coordinate_sha256", "")))
        if sha256_file(selected_path) != selected_hash:
            raise WorkflowValidationError(f"Selected-coordinate hash mismatch: {selected_path}")
        selected_frames = parse_xyz_frames(selected_path)
        if len(selected_frames) != 1 or not frames_numerically_equal(source_frame, selected_frames[0]):
            raise WorkflowValidationError(
                f"Selected coordinate is not the declared source ensemble frame: {selected_path}"
            )
        frame = selected_frames[0]

        source_dir = Path(source_row["job_directory"]).resolve()
        atom_map_path = source_dir / "atom_map.csv"
        isotope_map_path = source_dir / "isotope_map.csv"
        if Path(str(selection.get("source_atom_map_file", ""))).resolve() != atom_map_path:
            raise WorkflowValidationError(f"Selection atom-map path is not from its origin start: {source_job_id}")
        if Path(str(selection.get("source_isotope_map_file", ""))).resolve() != isotope_map_path:
            raise WorkflowValidationError(f"Selection isotope-map path is not from its origin start: {source_job_id}")
        atom_map_hash = sha256_file(atom_map_path)
        isotope_map_hash = sha256_file(isotope_map_path)
        if normalize_sha256(str(selection.get("source_atom_map_sha256", ""))) != atom_map_hash:
            raise WorkflowValidationError(f"Selection atom-map hash mismatch: {source_job_id}")
        if normalize_sha256(str(selection.get("source_isotope_map_sha256", ""))) != isotope_map_hash:
            raise WorkflowValidationError(f"Selection isotope-map hash mismatch: {source_job_id}")
        atom_rows, isotope_rows, isotope_masses = validate_cluster_maps(
            atom_map_path, isotope_map_path, frame
        )
        downstream_rows = downstream_atom_map_rows(
            atom_rows,
            isotope_masses,
            figure_id,
            scenario_id,
            scenario.get("species_smiles", ""),
        )
        downstream_fields = list(downstream_rows[0].keys())

        scenario_ordinals[scenario_id] += 1
        ordinal = scenario_ordinals[scenario_id]
        job_id = f"{safe_token(figure_id)}__{safe_token(scenario_id)}__c{ordinal:03d}"
        directory = output_root / safe_token(figure_id) / safe_token(scenario_id) / f"conformer_{ordinal:03d}"
        try:
            charge = int(source_row["charge"])
            multiplicity = int(source_row["multiplicity"])
        except (TypeError, ValueError) as exc:
            raise WorkflowValidationError(f"Invalid charge/multiplicity for {source_job_id}") from exc
        continuum = scenario.get("continuum_keyword", "") or source_row["continuum_keyword"]
        input_text = render_opt_freq_input(
            frame,
            charge,
            multiplicity,
            continuum,
            isotope_masses,
            nprocs=nprocs,
            maxcore_mb=maxcore_mb,
            temperature_K=thermochemistry_temperatures,
        )
        input_bytes = input_text.encode("utf-8")
        input_hash = sha256_bytes(input_bytes)
        input_bundle_hash = inline_orca_input_bundle_sha256("dft_opt_freq.inp", input_bytes)

        atom_map_text_path = directory / "nmr_atom_map_v2.csv"
        # Produce CSV bytes through a deterministic in-memory serialization.
        import csv as _csv
        import io as _io

        stream = _io.StringIO(newline="")
        writer = _csv.DictWriter(stream, fieldnames=downstream_fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(downstream_rows)
        downstream_map_bytes = ("\ufeff" + stream.getvalue()).encode("utf-8")

        selected_bytes = selected_path.read_bytes()
        source_atom_bytes = atom_map_path.read_bytes()
        source_isotope_bytes = isotope_map_path.read_bytes()
        dependency_hash = scientific_dependency_sha256(
            (
                ("approved_gate", gate_hash),
                ("manifest", manifest_hash),
                ("scenario_registry", registry_hash),
                ("resource_policy", resource_policy_hash),
                ("source_state", evidence["state_sha256"]),
                ("source_ensemble", artifact_hash),
                ("selected_coordinate", selected_hash),
                ("source_atom_map", atom_map_hash),
                ("source_isotope_map", isotope_map_hash),
                ("nucleus_conditions", nucleus_conditions_hash),
            )
        )
        provenance = {
            "schema_version": "2.0",
            "notice": NOTICE,
            "job_id": job_id,
            "figure_id": figure_id,
            "scenario_id": scenario_id,
            "condition_isolation_key": [figure_id, scenario_id],
            "approved_gate_file": str(gate_path.resolve()),
            "approved_gate_sha256": gate_hash,
            "approved_conformer_id": selection.get("conformer_id", ""),
            "source_goat_job_id": source_job_id,
            "source_state_file": str(evidence["state_path"]),
            "source_state_sha256": evidence["state_sha256"],
            "source_finalensemble_file": str(artifact_path),
            "source_finalensemble_sha256": artifact_hash,
            "source_frame_index_1based": frame_index,
            "selected_coordinate_source_file": str(selected_path),
            "selected_coordinate_sha256": selected_hash,
            "source_atom_map_file": str(atom_map_path),
            "source_atom_map_sha256": atom_map_hash,
            "source_isotope_map_file": str(isotope_map_path),
            "source_isotope_map_sha256": isotope_map_hash,
            "scenario_species_smiles": scenario.get("species_smiles", ""),
            "primary_solvent": scenario.get("primary_solvent", source_row.get("primary_solvent", "")),
            "continuum_keyword": continuum,
            "explicit_solvent_count": scenario.get(
                "explicit_solvent_count", source_row.get("explicit_solvent_count", "")
            ),
            "charge": charge,
            "multiplicity": multiplicity,
            "cluster_atom_count": len(frame.symbols),
            "deuterium_atom_indices_1based": sorted(isotope_masses),
            "method": "PBE0-D4/def2-TZVPP RIJCOSX CPCM(actual solvent) VeryTightOpt Freq VeryTightSCF DefGrid3",
            "acquisition_temperatures_K_by_nucleus": temperatures_by_nucleus,
            "thermochemistry_temperatures_K": thermochemistry_temperatures,
            "thermochemistry_temperature_policy": (
                "ORCA QRRHO free energies evaluated at every locked nucleus acquisition "
                "temperature from the same isotope-aware Hessian"
            ),
            "nucleus_conditions_sha256": nucleus_conditions_hash,
            "mpi_ranks": nprocs,
            "maxcore_mb_per_rank": maxcore_mb,
            "dft_input_file": "dft_opt_freq.inp",
            "dft_input_sha256": input_hash,
            "dft_input_bundle_sha256": input_bundle_hash,
            "resource_policy_sha256": resource_policy_hash,
            "scientific_dependency_sha256": dependency_hash,
            "downstream_real_nmr_calibration_reference": calibration_references,
            "checks": {
                "approved_gate_and_per_conformer_checks": PASS,
                "all_figure_goat_states_hash_bound": PASS,
                "scenario_composition_isolated": PASS,
                "selected_xyz_matches_declared_source_frame": PASS,
                "explicit_cluster_atom_order_preserved": PASS,
                "origin_atom_and_isotope_maps_hash_verified": PASS,
                "deuterium_inline_mass_exact": PASS,
                "environment_atoms_excluded_only_from_analyte_nmr_map_not_electronic_model": PASS,
            },
        }
        provenance_bytes = (json.dumps(provenance, indent=2, ensure_ascii=False) + "\n").encode("utf-8")
        files = {
            "approved_selected_cluster.xyz": selected_bytes,
            "source_atom_map.csv": source_atom_bytes,
            "source_isotope_map.csv": source_isotope_bytes,
            atom_map_text_path.name: downstream_map_bytes,
            "dft_opt_freq.inp": input_bytes,
            "provenance.json": provenance_bytes,
            "SCIENTIFIC_DEPENDENCY_SHA256.txt": (dependency_hash.upper() + "\n").encode("ascii"),
            "NOTICE.txt": (NOTICE + "\n").encode("utf-8"),
        }
        planned.append(
            PlannedJob(
                job_id=job_id,
                directory=directory,
                files=files,
                manifest_row={
                    "queue_order": 0,
                    "job_id": job_id,
                    "figure_id": figure_id,
                    "compound_number": source_row.get("compound_number", ""),
                    "display_name": source_row.get("display_name", ""),
                    "scenario_id": scenario_id,
                    "primary_solvent": provenance["primary_solvent"],
                    "continuum_keyword": continuum,
                    "explicit_solvent_count": provenance["explicit_solvent_count"],
                    "conformer_ordinal": ordinal,
                    "approved_conformer_id": selection.get("conformer_id", ""),
                    "charge": charge,
                    "multiplicity": multiplicity,
                    "source_goat_job_id": source_job_id,
                    "source_frame_index_1based": frame_index,
                    "job_directory": str(directory.resolve()),
                    "input_file": "dft_opt_freq.inp",
                    "output_file": "dft_opt_freq.out",
                    "required_artifacts": json.dumps(
                        [
                            "dft_opt_freq.gbw",
                            "dft_opt_freq.hess",
                            "dft_opt_freq.property.txt",
                            "dft_opt_freq.xyz",
                        ],
                        separators=(",", ":"),
                    ),
                    "input_sha256": input_hash,
                    "expected_input_bundle_sha256": input_bundle_hash,
                    "approved_gate_sha256": gate_hash,
                    "source_state_sha256": evidence["state_sha256"],
                    "source_finalensemble_sha256": artifact_hash,
                    "selected_coordinate_sha256": selected_hash,
                    "source_atom_map_sha256": atom_map_hash,
                    "source_isotope_map_sha256": isotope_map_hash,
                    "scientific_dependency_sha256": dependency_hash,
                    "acquisition_temperatures_K_by_nucleus": json.dumps(
                        temperatures_by_nucleus, separators=(",", ":"), sort_keys=True
                    ),
                    "thermochemistry_temperatures_K": json.dumps(
                        thermochemistry_temperatures, separators=(",", ":")
                    ),
                    "nucleus_conditions_sha256": nucleus_conditions_hash,
                    "enabled": "true",
                    "readiness": "ready",
                    "approval_status": APPROVED,
                    "notice": NOTICE,
                },
            )
        )
    return planned


def build_dft_jobs(
    project: Path,
    gate_paths: list[Path],
    manifest_path: Path,
    scenario_registry_path: Path,
    resource_policy_path: Path,
    output_root: Path,
    output_manifest: Path,
    replace_unstarted: bool = False,
    nucleus_conditions_path: Path | None = None,
) -> list[dict[str, object]]:
    project = Path(project).resolve()
    if not gate_paths:
        raise WorkflowValidationError("No approved gate files were supplied or discovered")
    manifest_path = ensure_within(Path(manifest_path), project, "GOAT explicit manifest")
    scenario_registry_path = ensure_within(
        Path(scenario_registry_path), project, "explicit scenario registry"
    )
    resource_policy_path = ensure_within(Path(resource_policy_path), project, "resource policy")
    nucleus_conditions_path = ensure_within(
        Path(nucleus_conditions_path or project / "config" / "nucleus_conditions.csv"),
        project,
        "nucleus conditions",
    )
    output_root = ensure_within(Path(output_root), project, "DFT output root")
    output_manifest = ensure_within(Path(output_manifest), project, "DFT output manifest")
    gate_paths = [ensure_within(Path(path), project, "approved GOAT selection gate") for path in gate_paths]
    manifest_rows = read_csv(manifest_path)
    manifest_by_job = _manifest_index(manifest_rows)
    registry = _registry_index(read_csv(scenario_registry_path))
    nprocs, maxcore_mb, resource_policy_hash = _resource_values(resource_policy_path)
    manifest_hash = sha256_file(manifest_path)
    registry_hash = sha256_file(scenario_registry_path)
    nucleus_condition_rows = read_csv(nucleus_conditions_path)
    nucleus_conditions_hash = sha256_file(nucleus_conditions_path)
    calibration_references = _calibration_references(project)
    planned: list[PlannedJob] = []
    seen_figures: set[str] = set()
    for gate_path in sorted(Path(item).resolve() for item in gate_paths):
        gate_figure = str(read_json(gate_path).get("figure_id", ""))
        if gate_figure in seen_figures:
            raise WorkflowValidationError(f"More than one approved gate supplied for {gate_figure}")
        seen_figures.add(gate_figure)
        planned.extend(
            _plan_gate_jobs(
                project,
                gate_path,
                manifest_rows,
                manifest_by_job,
                registry,
                Path(output_root),
                nprocs,
                maxcore_mb,
                resource_policy_hash,
                manifest_hash,
                registry_hash,
                calibration_references,
                nucleus_condition_rows,
                nucleus_conditions_hash,
            )
        )
    job_ids = [job.job_id for job in planned]
    if len(job_ids) != len(set(job_ids)):
        raise WorkflowValidationError("Planned DFT job IDs collide")
    for job in planned:
        _safe_existing_directory(job.directory, replace_unstarted)
    if Path(output_manifest).exists() and not replace_unstarted:
        raise WorkflowValidationError(
            f"Output manifest exists; use --replace-unstarted only after inspection: {output_manifest}"
        )

    rows: list[dict[str, object]] = []
    for queue_order, job in enumerate(
        sorted(
            planned,
            key=lambda item: (int(item.manifest_row["compound_number"]), item.job_id),
        ),
        start=1,
    ):
        for relative, data in job.files.items():
            atomic_write_bytes(job.directory / relative, data)
        row = dict(job.manifest_row)
        row["queue_order"] = queue_order
        rows.append(row)
    fields = list(rows[0].keys()) if rows else []
    if not rows:
        raise WorkflowValidationError("No DFT jobs were planned")
    atomic_write_csv(output_manifest, rows, fields)
    return rows


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--project", type=Path, default=PROJECT)
    value.add_argument("--gate", type=Path, action="append", dest="gates")
    value.add_argument("--manifest", type=Path)
    value.add_argument("--scenario-registry", type=Path)
    value.add_argument("--resource-policy", type=Path)
    value.add_argument("--output-root", type=Path)
    value.add_argument("--output-manifest", type=Path)
    value.add_argument("--replace-unstarted", action="store_true")
    value.add_argument("--nucleus-conditions", type=Path)
    return value


def main() -> None:
    args = parser().parse_args()
    project = args.project.resolve()
    gates = args.gates or sorted(
        (project / "jobs").glob(
            "F3_*/goat_ensemble_selection/goat_ensemble_selection_approved.json"
        )
    )
    rows = build_dft_jobs(
        project=project,
        gate_paths=gates,
        manifest_path=args.manifest or project / "data" / "GOAT_explicit_job_manifest.csv",
        scenario_registry_path=args.scenario_registry or project / "data" / "explicit_scenario_registry.csv",
        resource_policy_path=args.resource_policy or project / "config" / "resource_policy.json",
        output_root=args.output_root or project / "02_dft_opt_freq_explicit",
        output_manifest=args.output_manifest or project / "data" / "DFT_explicit_job_manifest.csv",
        replace_unstarted=args.replace_unstarted,
        nucleus_conditions_path=args.nucleus_conditions,
    )
    print(f"Built {len(rows)} DFT opt/freq inputs; ORCA was not started.")


if __name__ == "__main__":
    main()
