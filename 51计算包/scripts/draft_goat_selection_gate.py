"""Create review-only GOAT ensemble selection gates.

The program never writes an approved gate and never launches ORCA.  It first
validates every declared GOAT start for a figure, then proposes all unique,
converged conformers inside a conservative energy window, independently for
each explicit-solvent scenario.
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from workflow_schema import (
    HARTREE_TO_KCAL_MOL,
    NOTICE,
    WorkflowValidationError,
    atomic_write_json,
    atomic_write_text,
    fixed_order_kabsch_rmsd,
    normalize_sha256,
    parse_xyz_frames,
    read_csv,
    sha256_bytes,
    sha256_file,
    validate_cluster_maps,
    validate_completed_goat_state,
    xyz_text,
)


PROJECT = Path(__file__).resolve().parents[1]


def safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    if not token or token in {".", ".."}:
        raise WorkflowValidationError(f"Unsafe empty identifier: {value!r}")
    return token


def scenario_registry_index(rows: Iterable[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    index: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        key = (row.get("figure_id", "").strip(), row.get("scenario_id", "").strip())
        if not all(key) or key in index:
            raise WorkflowValidationError(f"Blank or duplicate scenario-registry key: {key}")
        index[key] = row
    return index


def _propose_for_figure(
    project: Path,
    figure_id: str,
    manifest_rows: list[dict[str, str]],
    registry: dict[tuple[str, str], dict[str, str]],
    state_root: Path,
    output_root: Path,
    energy_window_kcal_mol: float,
    duplicate_rmsd_A: float,
) -> tuple[Path, dict, list[tuple[Path, str]]]:
    figure_rows = [row for row in manifest_rows if row.get("figure_id") == figure_id]
    if not figure_rows:
        raise WorkflowValidationError(f"No explicit GOAT manifest rows for {figure_id}")
    by_scenario: dict[str, list[dict[str, object]]] = defaultdict(list)
    upstream_states: list[dict[str, object]] = []

    for row in sorted(figure_rows, key=lambda item: item["job_id"]):
        key = (figure_id, row["scenario_id"])
        if key not in registry:
            raise WorkflowValidationError(f"Scenario missing from registry: {key}")
        state_path = state_root / f"{row['job_id']}.json"
        evidence = validate_completed_goat_state(row, state_path, project)
        upstream_states.append(
            {
                "job_id": row["job_id"],
                "scenario_id": row["scenario_id"],
                "state_file": str(evidence["state_path"]),
                "sha256": evidence["state_sha256"],
                "input_bundle_sha256": evidence["input_bundle_sha256"],
                "output_file": str(evidence["output_path"]),
                "output_sha256": evidence["output_sha256"],
                "finalensemble_file": str(evidence["ensemble_path"]),
                "finalensemble_sha256": evidence["ensemble_sha256"],
            }
        )
        frames = parse_xyz_frames(Path(evidence["ensemble_path"]))
        source_dir = Path(evidence["job_directory"])
        atom_map_path = source_dir / "atom_map.csv"
        isotope_map_path = source_dir / "isotope_map.csv"
        for frame in frames:
            validate_cluster_maps(atom_map_path, isotope_map_path, frame)
            if frame.energy_hartree is None or frame.converged is not True:
                continue
            by_scenario[row["scenario_id"]].append(
                {
                    "row": row,
                    "evidence": evidence,
                    "frame": frame,
                    "atom_map_path": atom_map_path,
                    "atom_map_sha256": sha256_file(atom_map_path),
                    "isotope_map_path": isotope_map_path,
                    "isotope_map_sha256": sha256_file(isotope_map_path),
                }
            )

    expected_scenarios = sorted({row["scenario_id"] for row in figure_rows})
    if sorted(by_scenario) != expected_scenarios:
        absent = sorted(set(expected_scenarios) - set(by_scenario))
        raise WorkflowValidationError(f"No converged, energy-labelled GOAT frames for {figure_id}: {absent}")

    gate_dir = output_root / figure_id / "goat_ensemble_selection"
    writes: list[tuple[Path, str]] = []
    proposed: list[dict[str, object]] = []
    scenario_summaries: list[dict[str, object]] = []
    for scenario_id in expected_scenarios:
        raw = sorted(
            by_scenario[scenario_id],
            key=lambda item: (
                float(item["frame"].energy_hartree),
                str(item["row"]["job_id"]),
                int(item["frame"].source_frame_index_1based or 0),
            ),
        )
        minimum = float(raw[0]["frame"].energy_hartree)
        within = [
            item
            for item in raw
            if (float(item["frame"].energy_hartree) - minimum) * HARTREE_TO_KCAL_MOL
            <= energy_window_kcal_mol + 1e-9
        ]
        unique: list[dict[str, object]] = []
        for item in within:
            frame = item["frame"]
            if all(fixed_order_kabsch_rmsd(frame, kept["frame"]) > duplicate_rmsd_A for kept in unique):
                unique.append(item)
        if not unique:
            raise WorkflowValidationError(f"Selection unexpectedly empty for {figure_id}/{scenario_id}")
        scenario_summaries.append(
            {
                "scenario_id": scenario_id,
                "completed_goat_starts": len(
                    {str(item["row"]["job_id"]) for item in by_scenario[scenario_id]}
                ),
                "converged_frames": len(raw),
                "frames_inside_energy_window": len(within),
                "proposed_unique_conformers": len(unique),
                "minimum_energy_hartree": minimum,
            }
        )
        for ordinal, item in enumerate(unique, start=1):
            row = item["row"]
            evidence = item["evidence"]
            frame = item["frame"]
            conformer_id = f"conf{ordinal:03d}"
            coordinate_path = (
                gate_dir
                / "draft_candidates"
                / safe_token(scenario_id)
                / f"{conformer_id}__{safe_token(row['job_id'])}__frame{frame.source_frame_index_1based:04d}.xyz"
            )
            coordinate_text = xyz_text(
                frame,
                (
                    f"{NOTICE}; DRAFT candidate; {figure_id}; {scenario_id}; "
                    f"source={row['job_id']}; source_frame={frame.source_frame_index_1based}; "
                    f"energy_hartree={frame.energy_hartree:.12f}"
                ),
            )
            writes.append((coordinate_path, coordinate_text))
            proposed.append(
                {
                    "condition_scenario_id": scenario_id,
                    "conformer_id": conformer_id,
                    "source_goat_job_id": row["job_id"],
                    "source_goat_stage": row["job_id"],
                    "source_state_file": str(evidence["state_path"]),
                    "source_state_sha256": evidence["state_sha256"],
                    "source_artifact_file": str(evidence["ensemble_path"]),
                    "source_artifact_sha256": evidence["ensemble_sha256"],
                    "source_frame_index_1based": frame.source_frame_index_1based,
                    "source_frame_energy_hartree": frame.energy_hartree,
                    "relative_energy_kcal_mol": (
                        float(frame.energy_hartree) - minimum
                    )
                    * HARTREE_TO_KCAL_MOL,
                    "source_frame_canonical_sha256": sha256_bytes(xyz_text(frame).encode("utf-8")),
                    "selected_coordinate_file": str(coordinate_path.resolve()),
                    "selected_coordinate_sha256": sha256_bytes(coordinate_text.encode("utf-8")),
                    "source_atom_map_file": str(item["atom_map_path"]),
                    "source_atom_map_sha256": item["atom_map_sha256"],
                    "source_isotope_map_file": str(item["isotope_map_path"]),
                    "source_isotope_map_sha256": item["isotope_map_sha256"],
                    "connectivity_stereochemistry_check": "REVIEW_REQUIRED",
                    "cluster_composition_check": "REVIEW_REQUIRED",
                    "manual_extraction_check": "REVIEW_REQUIRED",
                    "approval_status": "DRAFT_NOT_APPROVED",
                }
            )

    gate_path = gate_dir / "goat_ensemble_selection_DRAFT.json"
    gate = {
        "schema_version": "2.0",
        "notice": NOTICE,
        "approval_status": "DRAFT_NOT_APPROVED",
        "figure_id": figure_id,
        "stage": "goat_ensemble_selection",
        "model_id": "all_enabled_explicit_scenarios",
        "selection_energy_window_kcal_mol": energy_window_kcal_mol,
        "duplicate_rmsd_threshold_A": duplicate_rmsd_A,
        "rmsd_policy": "FIXED_ATOM_ORDER_ALL_CLUSTER_ATOMS_KABSCH_NO_REFLECTION",
        "selection_completeness_check": "REVIEW_REQUIRED",
        "protonation_and_cluster_integrity_review": "REVIEW_REQUIRED",
        "scenario_summaries": scenario_summaries,
        "upstream_state_hashes": upstream_states,
        "selected_conformers": proposed,
        "waiver": {"approval_status": "NOT_APPLICABLE", "reason": ""},
        "reviewer": "",
        "reviewed_at": "",
        "decision_notes": [
            "This file is a DRAFT. The generator cannot approve a conformer selection.",
            "Review each explicit-solvent scenario independently; never merge compositions.",
            "Before DFT, set the top-level status to APPROVED_FOR_OPT_FREQ, both top-level checks to PASS, and every retained conformer status/check to its required approved/PASS value.",
            "Remove rejected conformers rather than labelling them approved.",
            "Save the reviewed copy as goat_ensemble_selection_approved.json; do not overwrite upstream GOAT artifacts or maps.",
        ],
    }
    return gate_path, gate, writes


def build_draft_gates(
    project: Path,
    manifest_path: Path,
    scenario_registry_path: Path,
    state_root: Path,
    output_root: Path,
    figure_ids: list[str] | None = None,
    energy_window_kcal_mol: float = 12.0,
    duplicate_rmsd_A: float = 0.25,
    replace_draft: bool = False,
) -> list[Path]:
    project = Path(project).resolve()
    manifest_rows = read_csv(manifest_path)
    registry = scenario_registry_index(read_csv(scenario_registry_path))
    all_figures = sorted({row.get("figure_id", "") for row in manifest_rows})
    selected_figures = sorted(set(figure_ids or all_figures))
    unknown = sorted(set(selected_figures) - set(all_figures))
    if unknown:
        raise WorkflowValidationError(f"Unknown figure IDs: {unknown}")
    if energy_window_kcal_mol <= 0 or duplicate_rmsd_A <= 0:
        raise WorkflowValidationError("Energy window and RMSD threshold must be positive")

    proposals = [
        _propose_for_figure(
            project,
            figure_id,
            manifest_rows,
            registry,
            Path(state_root),
            Path(output_root),
            energy_window_kcal_mol,
            duplicate_rmsd_A,
        )
        for figure_id in selected_figures
    ]
    for gate_path, _gate, _writes in proposals:
        approved = gate_path.with_name("goat_ensemble_selection_approved.json")
        if approved.exists():
            raise WorkflowValidationError(
                f"Approved gate already exists; refusing to create a competing draft: {approved}"
            )
        if gate_path.exists() and not replace_draft:
            raise WorkflowValidationError(f"Draft gate exists; use --replace-draft to refresh it: {gate_path}")

    written: list[Path] = []
    for gate_path, gate, writes in proposals:
        for coordinate_path, coordinate_text in writes:
            atomic_write_text(coordinate_path, coordinate_text)
        atomic_write_json(gate_path, gate)
        written.append(gate_path)
    return written


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--project", type=Path, default=PROJECT)
    value.add_argument("--manifest", type=Path)
    value.add_argument("--scenario-registry", type=Path)
    value.add_argument("--state-root", type=Path)
    value.add_argument("--output-root", type=Path)
    value.add_argument("--figure-id", action="append", dest="figure_ids")
    value.add_argument("--energy-window-kcal-mol", type=float, default=12.0)
    value.add_argument("--duplicate-rmsd-A", type=float, default=0.25)
    value.add_argument("--replace-draft", action="store_true")
    return value


def main() -> None:
    args = parser().parse_args()
    project = args.project.resolve()
    paths = build_draft_gates(
        project=project,
        manifest_path=args.manifest or project / "data" / "GOAT_explicit_job_manifest.csv",
        scenario_registry_path=args.scenario_registry or project / "data" / "explicit_scenario_registry.csv",
        state_root=args.state_root or project / "logs" / "goat_explicit_batch" / "jobs",
        output_root=args.output_root or project / "jobs",
        figure_ids=args.figure_ids,
        energy_window_kcal_mol=args.energy_window_kcal_mol,
        duplicate_rmsd_A=args.duplicate_rmsd_A,
        replace_draft=args.replace_draft,
    )
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
