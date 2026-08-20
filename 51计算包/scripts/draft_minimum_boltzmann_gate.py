"""Create DRAFT-only minimum/QRRHO/Boltzmann review gates.

The script validates every opt/freq state for a Figure, parses ORCA-reported
QRRHO thermochemistry, and proposes normalized weights independently for each
explicit composition, nucleus and locked temperature.  It never approves a
production minimum and never launches ORCA.
"""

from __future__ import annotations

import argparse
import math
import re
from collections import defaultdict
from pathlib import Path

from workflow_schema import (
    NOTICE,
    WorkflowValidationError,
    atomic_write_json,
    parse_xyz_frames,
    read_csv,
    read_text_auto,
    sha256_file,
    validate_cluster_maps,
    validate_completed_opt_freq_state,
)


PROJECT = Path(__file__).resolve().parents[1]
R_KCAL_MOL_K = 0.00198720425864083
SUPPORTED_NUCLEI = {"H": "1H", "C": "13C", "F": "19F", "P": "31P"}


def parse_orca_qrrho_temperatures(path: Path) -> list[dict[str, object]]:
    """Parse every ORCA QRRHO thermochemistry block from one Hessian run.

    ORCA 6.1 accepts up to 16 values in ``%freq Temp`` and prints a complete
    thermochemistry block for each value.  Treating the final block as if it
    applied to every nucleus silently gives wrong conformer populations, so the
    multi-temperature representation is the authoritative one.
    """
    text = read_text_auto(path)

    def last_float(pattern: str, label: str) -> float:
        matches = re.findall(pattern, text, flags=re.IGNORECASE | re.MULTILINE)
        if not matches:
            raise WorkflowValidationError(f"Cannot parse {label} from {path}")
        return float(matches[-1])

    electronic = last_float(
        r"FINAL\s+SINGLE\s+POINT\s+ENERGY\s+([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)",
        "electronic energy",
    )
    cutoff = last_float(
        r"QRRHO[\s\S]{0,180}?reference\s+frequency\s+of\s+([-+]?\d+(?:\.\d+)?)\s+cm-1",
        "QRRHO reference frequency",
    )
    if "Vibrational entropy computed according to the QRRHO of S. Grimme" not in text:
        raise WorkflowValidationError(f"ORCA output does not declare Grimme QRRHO entropy: {path}")
    frequencies = [
        float(value)
        for value in re.findall(
            r"^\s*\d+:\s+([-+]?\d+(?:\.\d+)?)\s+cm\*\*-1",
            text,
            flags=re.MULTILINE,
        )
    ]
    if not frequencies:
        raise WorkflowValidationError(f"Cannot parse vibrational frequencies from {path}")
    imaginary = [value for value in frequencies if value < -0.1]
    temperature_pattern = re.compile(
        r"Temperature\s+\.\.\.\s+([-+]?\d+(?:\.\d+)?)\s+K",
        flags=re.IGNORECASE,
    )
    gibbs_pattern = re.compile(
        r"Final\s+Gibbs\s+free\s+energy\s+\.\.\.\s+"
        r"([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)\s+Eh",
        flags=re.IGNORECASE,
    )
    temperature_matches = list(temperature_pattern.finditer(text))
    if not temperature_matches:
        raise WorkflowValidationError(f"Cannot parse thermochemistry temperature from {path}")
    records: list[dict[str, object]] = []
    seen: list[float] = []
    for position, match in enumerate(temperature_matches):
        start = match.end()
        stop = (
            temperature_matches[position + 1].start()
            if position + 1 < len(temperature_matches)
            else len(text)
        )
        gibbs_matches = gibbs_pattern.findall(text[start:stop])
        if not gibbs_matches:
            raise WorkflowValidationError(
                f"Temperature {match.group(1)} K lacks a matching Gibbs free energy in {path}"
            )
        temperature = float(match.group(1))
        if any(math.isclose(temperature, prior, rel_tol=0.0, abs_tol=1e-8) for prior in seen):
            raise WorkflowValidationError(f"Duplicate thermochemistry temperature {temperature} K in {path}")
        seen.append(temperature)
        records.append(
            {
                "electronic_energy_hartree": electronic,
                "qrrho_gibbs_free_energy_hartree": float(gibbs_matches[-1]),
                "temperature_K": temperature,
                "qrrho_cutoff_cm1": cutoff,
                "frequency_scale_factor": 1.0,
                "frequency_count_including_translation_rotation_zeros": len(frequencies),
                "imaginary_frequencies_cm1": imaginary,
                "imaginary_frequency_count": len(imaginary),
                "thermochemistry_source": "ORCA_6_1_REPORTED_GRIMME_QRRHO_SAME_HESSIAN_MULTI_T",
            }
        )
    return records


def parse_orca_qrrho(path: Path, temperature_K: float | None = None) -> dict[str, object]:
    """Backward-compatible single-block view, optionally selected by temperature."""
    records = parse_orca_qrrho_temperatures(path)
    if temperature_K is None:
        return records[-1]
    matches = [
        record
        for record in records
        if math.isclose(float(record["temperature_K"]), temperature_K, rel_tol=0.0, abs_tol=0.02)
    ]
    if len(matches) != 1:
        raise WorkflowValidationError(
            f"Expected one QRRHO block at {temperature_K} K in {path}; found {len(matches)}"
        )
    return matches[0]


def locked_nucleus_conditions(path: Path, figure_id: str) -> dict[str, dict[str, object]]:
    """Load the per-nucleus acquisition-temperature lock for one compound."""
    rows = [row for row in read_csv(path) if row.get("figure_id", "").strip() == figure_id]
    conditions: dict[str, dict[str, object]] = {}
    for row in rows:
        nucleus = row.get("nucleus", "").strip()
        if nucleus not in set(SUPPORTED_NUCLEI.values()):
            continue
        if nucleus in conditions:
            raise WorkflowValidationError(f"Duplicate {figure_id}/{nucleus} condition in {path}")
        try:
            temperature = float(row["temperature_K"])
        except (KeyError, TypeError, ValueError) as exc:
            raise WorkflowValidationError(f"Invalid {figure_id}/{nucleus} temperature in {path}") from exc
        status = row.get("condition_status", "").strip()
        if temperature <= 0.0 or not status or status.upper() == "UNKNOWN":
            raise WorkflowValidationError(f"Unlocked {figure_id}/{nucleus} acquisition condition in {path}")
        conditions[nucleus] = {
            "temperature_K": temperature,
            "condition_status": status,
            "solvent": row.get("solvent", "").strip(),
            "source_note": row.get("source_note", "").strip(),
        }
    if not conditions:
        raise WorkflowValidationError(f"No locked nucleus conditions for {figure_id} in {path}")
    return conditions


def analyte_nuclei(map_path: Path) -> dict[str, list[int]]:
    result: dict[str, list[int]] = defaultdict(list)
    for row in read_csv(map_path):
        if row.get("nmr_scope", "") != "ANALYTE":
            continue
        nucleus = SUPPORTED_NUCLEI.get(row.get("element", "").strip().upper())
        if nucleus:
            result[nucleus].append(int(row["atom_index_1based"]))
    return {key: sorted(values) for key, values in result.items()}


def normalized_weights(energies_hartree: dict[str, float], temperature_K: float) -> dict[str, float]:
    if not energies_hartree or temperature_K <= 0:
        raise WorkflowValidationError("Cannot calculate Boltzmann weights without energies and T > 0")
    minimum = min(energies_hartree.values())
    factors = {
        key: math.exp(-((value - minimum) * 627.5094740631) / (R_KCAL_MOL_K * temperature_K))
        for key, value in energies_hartree.items()
    }
    partition = math.fsum(factors.values())
    weights = {key: value / partition for key, value in factors.items()}
    if not math.isclose(math.fsum(weights.values()), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise WorkflowValidationError("Boltzmann normalization failed")
    return weights


def _proposal(
    project: Path,
    figure_id: str,
    rows: list[dict[str, str]],
    state_root: Path,
    output_root: Path,
    nucleus_conditions_path: Path,
) -> tuple[Path, dict]:
    figure_rows = [row for row in rows if row.get("figure_id") == figure_id]
    if not figure_rows:
        raise WorkflowValidationError(f"No DFT manifest rows for {figure_id}")
    locked_conditions = locked_nucleus_conditions(nucleus_conditions_path, figure_id)
    decisions: list[dict[str, object]] = []
    upstream: list[dict[str, object]] = []
    by_scenario: dict[str, list[dict[str, object]]] = defaultdict(list)
    nuclei_by_scenario: dict[str, dict[str, list[int]]] = {}
    composition_by_scenario: dict[str, tuple[str, ...]] = {}
    for row in sorted(figure_rows, key=lambda item: item["job_id"]):
        state_path = state_root / f"{row['job_id']}.json"
        evidence = validate_completed_opt_freq_state(row, state_path, project)
        thermochemistry = {
            nucleus: parse_orca_qrrho(
                Path(evidence["output_path"]), float(condition["temperature_K"])
            )
            for nucleus, condition in sorted(locked_conditions.items())
        }
        primary_nucleus = "1H" if "1H" in thermochemistry else sorted(thermochemistry)[0]
        thermo = thermochemistry[primary_nucleus]
        if any(
            not math.isclose(float(item["qrrho_cutoff_cm1"]), 100.0, rel_tol=0.0, abs_tol=0.01)
            for item in thermochemistry.values()
        ):
            raise WorkflowValidationError(
                f"{row['job_id']} used QRRHO cutoff {thermo['qrrho_cutoff_cm1']} cm-1, expected 100.0"
            )
        optimized = Path(evidence["artifacts"]["optimized_xyz"]["path"])
        frames = parse_xyz_frames(optimized)
        if len(frames) != 1:
            raise WorkflowValidationError(f"Optimized XYZ must contain one frame: {optimized}")
        job_dir = Path(evidence["job_directory"])
        validate_cluster_maps(job_dir / "source_atom_map.csv", job_dir / "source_isotope_map.csv", frames[0])
        nuclei = analyte_nuclei(job_dir / "nmr_atom_map_v2.csv")
        scenario_id = row["scenario_id"]
        if scenario_id in nuclei_by_scenario and nuclei_by_scenario[scenario_id] != nuclei:
            raise WorkflowValidationError(f"Analyte nucleus/index map differs inside scenario {scenario_id}")
        nuclei_by_scenario[scenario_id] = nuclei
        composition = (
            str(row.get("primary_solvent", "")),
            str(row.get("continuum_keyword", "")),
            str(row.get("explicit_solvent_count", "")),
            str(row.get("charge", "")),
            str(row.get("multiplicity", "")),
            str(evidence["provenance"].get("scenario_species_smiles", "")),
        )
        if scenario_id in composition_by_scenario and composition_by_scenario[scenario_id] != composition:
            raise WorkflowValidationError(f"Composition differs inside scenario {scenario_id}")
        composition_by_scenario[scenario_id] = composition
        conformer_id = row["job_id"]
        upstream.append(
            {
                "job_id": row["job_id"],
                "scenario_id": scenario_id,
                "state_file": str(evidence["state_path"]),
                "sha256": evidence["state_sha256"],
                "input_bundle_sha256": evidence["input_bundle_sha256"],
                "output_sha256": evidence["output_sha256"],
            }
        )
        decision = {
            "condition_scenario_id": scenario_id,
            "conformer_id": conformer_id,
            "upstream_dft_job_id": row["job_id"],
            "opt_freq_stage": "opt_freq",
            "opt_freq_state_file": str(evidence["state_path"]),
            "opt_freq_state_sha256": evidence["state_sha256"],
            "optimized_coordinate_file": str(optimized),
            "optimized_coordinate_sha256": evidence["artifacts"]["optimized_xyz"]["sha256"],
            "hessian_file": str(evidence["artifacts"]["hess"]["path"]),
            "hessian_sha256": evidence["artifacts"]["hess"]["sha256"],
            "nmr_atom_map_file": str((job_dir / "nmr_atom_map_v2.csv").resolve()),
            "nmr_atom_map_sha256": sha256_file(job_dir / "nmr_atom_map_v2.csv"),
            **thermo,
            "primary_thermochemistry_nucleus_for_legacy_fields": primary_nucleus,
            "thermochemistry_by_nucleus": thermochemistry,
            "automated_zero_imaginary_screen": (
                "PASS" if int(thermo["imaginary_frequency_count"]) == 0 else "FAIL"
            ),
            "stationary_point_check": "REVIEW_REQUIRED",
            "cluster_integrity_check": "REVIEW_REQUIRED",
            "approval_status": "DRAFT_NOT_APPROVED",
            "decision_reason": "",
        }
        decisions.append(decision)
        by_scenario[scenario_id].append(decision)

    weight_sets: list[dict[str, object]] = []
    scenario_summaries: list[dict[str, object]] = []
    for scenario_id in sorted(by_scenario):
        candidates = [
            item for item in by_scenario[scenario_id] if item["automated_zero_imaginary_screen"] == "PASS"
        ]
        if not candidates:
            raise WorkflowValidationError(f"No zero-imaginary candidate for {figure_id}/{scenario_id}")
        energies = {
            str(item["conformer_id"]): float(item["qrrho_gibbs_free_energy_hartree"])
            for item in candidates
        }
        for nucleus in sorted(nuclei_by_scenario[scenario_id]):
            if nucleus not in locked_conditions:
                raise WorkflowValidationError(
                    f"Analyte contains {nucleus}, but {figure_id} lacks its acquisition-temperature lock"
                )
            temperature = float(locked_conditions[nucleus]["temperature_K"])
            energies = {
                str(item["conformer_id"]): float(
                    item["thermochemistry_by_nucleus"][nucleus][
                        "qrrho_gibbs_free_energy_hartree"
                    ]
                )
                for item in candidates
            }
            weights = normalized_weights(energies, temperature)
            minimum = min(energies.values())
            weight_sets.append(
                {
                    "condition_scenario_id": scenario_id,
                    "nucleus": nucleus,
                    "temperature_K": temperature,
                    "condition_status": locked_conditions[nucleus]["condition_status"],
                    "condition_source_note": locked_conditions[nucleus]["source_note"],
                    "energy_model": "ORCA_6_1_REPORTED_GRIMME_QRRHO_GIBBS_FREE_ENERGY_SAME_HESSIAN_AT_NUCLEUS_T",
                    "qrrho_cutoff_cm1": 100.0,
                    "frequency_scale_factor": 1.0,
                    "weights": [
                        {
                            "conformer_id": key,
                            "qrrho_gibbs_free_energy_hartree": energies[key],
                            "relative_qrrho_free_energy_kcal_mol": (
                                energies[key] - minimum
                            )
                            * 627.5094740631,
                            "weight": weights[key],
                        }
                        for key in sorted(weights)
                    ],
                    "normalization_sum": math.fsum(weights.values()),
                    "approval_status": "DRAFT_NOT_APPROVED",
                }
            )
        scenario_summaries.append(
            {
                "scenario_id": scenario_id,
                "composition_lock": composition_by_scenario[scenario_id],
                "completed_opt_freq_conformers": len(by_scenario[scenario_id]),
                "zero_imaginary_candidates": len(candidates),
                "analyte_nuclei": nuclei_by_scenario[scenario_id],
            }
        )

    path = output_root / figure_id / "minimum_boltzmann_audit" / "minimum_boltzmann_DRAFT.json"
    payload = {
        "schema_version": "3.0",
        "notice": NOTICE,
        "approval_status": "DRAFT_NOT_APPROVED",
        "figure_id": figure_id,
        "stage": "minimum_boltzmann_audit",
        "model_id": "all_enabled_explicit_scenarios",
        "thermochemistry_method": "ORCA_6_1_REPORTED_GRIMME_QRRHO",
        "frequency_scale_factor": 1.0,
        "qrrho_cutoff_cm1": 100.0,
        "same_composition_only_check": "REVIEW_REQUIRED",
        "temperature_policy": (
            "Weights are separate by scenario/nucleus/acquisition temperature. Each temperature-specific "
            "ORCA QRRHO free energy was re-evaluated from the same isotope-aware Hessian."
        ),
        "nucleus_conditions_file": str(nucleus_conditions_path.resolve()),
        "nucleus_conditions_sha256": sha256_file(nucleus_conditions_path),
        "locked_nucleus_conditions": locked_conditions,
        "scenario_summaries": scenario_summaries,
        "upstream_state_hashes": upstream,
        "conformer_decisions": decisions,
        "boltzmann_weight_sets": weight_sets,
        "waiver": {"approval_status": "NOT_APPLICABLE", "reason": ""},
        "reviewer": "",
        "reviewed_at": "",
        "decision_notes": [
            "This generator never approves minima or weights.",
            "There must be one decision for every completed opt/freq state.",
            "APPROVED_FOR_NMR requires zero imaginary frequencies plus manual stationary-point and cluster-integrity PASS.",
            "Rejected conformers must use REJECTED_AFTER_OPT_FREQ with a non-empty reason.",
            "After any rejection, recalculate and normalize each scenario/nucleus/temperature weight set using only approved conformers.",
            "Do not combine different explicit solvent counts, protonation states, counterion arrangements or solvents.",
        ],
    }
    return path, payload


def build_draft_minimum_gates(
    project: Path,
    manifest_path: Path,
    state_root: Path,
    output_root: Path,
    figure_ids: list[str] | None = None,
    nucleus_conditions_path: Path | None = None,
    replace_draft: bool = False,
) -> list[Path]:
    project = Path(project).resolve()
    rows = read_csv(manifest_path)
    all_figures = sorted({row["figure_id"] for row in rows})
    selected = sorted(set(figure_ids or all_figures))
    unknown = sorted(set(selected) - set(all_figures))
    if unknown:
        raise WorkflowValidationError(f"Unknown Figure IDs: {unknown}")
    resolved_conditions = Path(nucleus_conditions_path or project / "config" / "nucleus_conditions.csv")
    proposals = [
        _proposal(project, figure_id, rows, Path(state_root), Path(output_root), resolved_conditions)
        for figure_id in selected
    ]
    for path, _payload in proposals:
        approved = path.with_name("minimum_boltzmann_approved.json")
        if approved.exists():
            raise WorkflowValidationError(f"Approved minimum gate already exists: {approved}")
        if path.exists() and not replace_draft:
            raise WorkflowValidationError(f"DRAFT minimum gate already exists: {path}")
    for path, payload in proposals:
        atomic_write_json(path, payload)
    return [path for path, _payload in proposals]


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--project", type=Path, default=PROJECT)
    value.add_argument("--manifest", type=Path)
    value.add_argument("--state-root", type=Path)
    value.add_argument("--output-root", type=Path)
    value.add_argument("--figure-id", action="append", dest="figure_ids")
    value.add_argument("--nucleus-conditions", type=Path)
    value.add_argument("--replace-draft", action="store_true")
    return value


def main() -> None:
    args = parser().parse_args()
    project = args.project.resolve()
    paths = build_draft_minimum_gates(
        project=project,
        manifest_path=args.manifest or project / "data" / "DFT_explicit_job_manifest.csv",
        state_root=args.state_root or project / "logs" / "opt_freq_serial" / "jobs",
        output_root=args.output_root or project / "jobs",
        figure_ids=args.figure_ids,
        nucleus_conditions_path=args.nucleus_conditions,
        replace_draft=args.replace_draft,
    )
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
