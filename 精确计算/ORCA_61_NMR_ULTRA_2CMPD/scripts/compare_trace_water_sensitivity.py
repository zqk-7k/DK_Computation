"""Compare matched dry and trace-water NMR prediction bundles.

The comparison is conditional on a fixed explicit-solvent composition.  It
never mixes dry and wet populations and never selects the branch that happens
to be closer to experiment.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import uuid
from pathlib import Path
from typing import Any, Sequence


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT / "config" / "workflow_profiles.json"
DEFAULT_SOLUTION_EQUIVALENCE = PROJECT / "config" / "solution_equivalence_groups.json"
NOTICE = "CALCULATED WATER-SENSITIVITY COMPARISON — NOT EXPERIMENTAL DATA"
SHIFT_FIELDS = {
    "raw_tms": "raw_tms_reference_shift_ppm",
    "provisional_regression": "provisional_external_regression_shift_ppm",
}


class WaterSensitivityError(RuntimeError):
    """Raised when dry/wet predictions are not scientifically comparable."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def safe_token(value: str) -> str:
    token = "".join(character if character.isalnum() or character in "_.-" else "_" for character in value)
    token = token.strip("_.")
    if not token or token in {".", ".."}:
        raise WaterSensitivityError(f"Unsafe identifier: {value!r}")
    return token


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise WaterSensitivityError(f"Missing input: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise WaterSensitivityError(f"Expected a JSON object: {path}")
    return payload


def atomic_write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{uuid.uuid4().hex}")
    try:
        temporary.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        raise WaterSensitivityError("Cannot write an empty water-sensitivity table")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    temporary = path.with_name(f".{path.name}.tmp.{uuid.uuid4().hex}")
    try:
        with temporary.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def selected_profile(config: dict[str, Any], profile_name: str | None) -> tuple[str, dict[str, Any]]:
    selected = str(profile_name or config.get("default_profile", ""))
    profiles = config.get("profiles", {})
    if not selected or not isinstance(profiles, dict) or selected not in profiles:
        raise WaterSensitivityError(f"Unknown workflow profile: {selected!r}")
    profile = profiles[selected]
    if not isinstance(profile, dict):
        raise WaterSensitivityError(f"Malformed workflow profile: {selected}")
    return selected, profile


def bundle_identity(bundle: dict[str, Any], path: Path) -> tuple[str, str, str]:
    compound = bundle.get("compound")
    nuclei = bundle.get("averaged_nuclei")
    if not isinstance(compound, dict) or not isinstance(nuclei, list):
        raise WaterSensitivityError(f"Malformed prediction bundle: {path}")
    figure_id = str(compound.get("figure_id", ""))
    scenario_id = str(compound.get("scenario_id", ""))
    molecular_identity = str(
        compound.get("inchikey") or compound.get("inchi") or compound.get("isomeric_smiles") or ""
    )
    if not figure_id or not scenario_id:
        raise WaterSensitivityError(f"Bundle lacks figure/scenario identity: {path}")
    return figure_id, scenario_id, molecular_identity


def scientific_signature(bundle: dict[str, Any], path: Path) -> dict[str, Any]:
    methods = bundle.get("methods")
    weighting = bundle.get("weighting")
    references = bundle.get("reference_shieldings")
    reference_approval = bundle.get("reference_registry_approval")
    condition_state = bundle.get("condition_state")
    source_manifest_records = bundle.get("source_manifest_hashes")
    if not isinstance(methods, dict) or not isinstance(weighting, dict):
        raise WaterSensitivityError(f"Bundle lacks method/weighting provenance: {path}")
    if not isinstance(references, dict) or not isinstance(reference_approval, dict):
        raise WaterSensitivityError(f"Bundle lacks reference provenance: {path}")
    if not isinstance(condition_state, dict) or not isinstance(source_manifest_records, dict):
        raise WaterSensitivityError(f"Bundle lacks gate/source-manifest provenance: {path}")
    gate_binding = condition_state.get("weight_gate_binding")
    if not isinstance(gate_binding, dict):
        raise WaterSensitivityError(f"Bundle lacks weight-gate binding: {path}")
    gate_hash = str(gate_binding.get("approved_minimum_gate_sha256", "")).upper()
    if not gate_hash or len(gate_hash) != 64 or any(
        character not in "0123456789ABCDEF" for character in gate_hash
    ):
        raise WaterSensitivityError(f"Invalid approved minimum-gate hash: {path}")
    required_manifests = {
        "dft_manifest",
        "weight_table",
        "nucleus_conditions",
        "shielding_manifest",
        "coupling_manifest",
        "reference_manifest",
        "reference_registry",
    }
    manifest_hashes: dict[str, str] = {}
    for key in sorted(required_manifests):
        record = source_manifest_records.get(key)
        digest = str(record.get("sha256", "")).upper() if isinstance(record, dict) else ""
        if len(digest) != 64 or any(
            character not in "0123456789ABCDEF" for character in digest
        ):
            raise WaterSensitivityError(f"Invalid source-manifest hash for {key}: {path}")
        manifest_hashes[key] = digest
    fingerprint = str(methods.get("shielding_method_fingerprint", ""))
    method = str(methods.get("shielding_method", ""))
    tms_output = str(methods.get("tms_reference_output", ""))
    dft_keywords = str(methods.get("dft_opt_freq_input_keywords", ""))
    qrrho_protocol = methods.get("opt_freq_qrrho_protocol")
    qrrho_fingerprint = str(methods.get("opt_freq_qrrho_protocol_fingerprint", ""))
    shielding_keywords = str(methods.get("shielding_input_keywords", ""))
    weighting_source = str(weighting.get("source", ""))
    if (
        not fingerprint
        or not method
        or not tms_output
        or not dft_keywords
        or not isinstance(qrrho_protocol, dict)
        or not qrrho_fingerprint
        or not shielding_keywords
        or not weighting_source
    ):
        raise WaterSensitivityError(f"Bundle has incomplete shielding/TMS identity: {path}")
    expected_qrrho_fingerprint = hashlib.sha256(
        json.dumps(qrrho_protocol, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest().upper()
    if qrrho_fingerprint.upper() != expected_qrrho_fingerprint:
        raise WaterSensitivityError(f"Invalid opt/freq-QRRHO protocol fingerprint: {path}")
    nucleus_sets = weighting.get("nucleus_weight_sets")
    if not isinstance(nucleus_sets, dict) or not nucleus_sets:
        raise WaterSensitivityError(f"Bundle lacks nucleus-specific weight conditions: {path}")
    conditions: dict[str, dict[str, Any]] = {}
    for nucleus, record in nucleus_sets.items():
        if not isinstance(record, dict) or record.get("temperature_K") is None:
            raise WaterSensitivityError(f"Malformed weight condition for {nucleus}: {path}")
        conditions[str(nucleus)] = {
            "temperature_K": float(record["temperature_K"]),
            "condition_status": str(record.get("condition_status", "")),
            "condition_source_note": str(record.get("condition_source_note", "")),
        }
    return {
        "shielding_method": method,
        "shielding_method_fingerprint": fingerprint,
        "dft_opt_freq_input_keywords": dft_keywords,
        "opt_freq_qrrho_protocol": qrrho_protocol,
        "opt_freq_qrrho_protocol_fingerprint": qrrho_fingerprint.upper(),
        "shielding_input_keywords": shielding_keywords,
        "tms_reference_output": tms_output,
        "reference_shieldings": references,
        "reference_registry_approval": reference_approval,
        "approved_minimum_gate_sha256": gate_hash,
        "source_manifest_hashes": manifest_hashes,
        "empirical_regression_policy": bundle.get("empirical_regression_policy"),
        "weighting_policy": weighting.get("policy"),
        "weighting_source": weighting_source,
        "nucleus_weight_conditions": conditions,
    }


def analyte_nuclei(bundle: dict[str, Any], shift_field: str) -> dict[tuple[str, int], dict[str, Any]]:
    result: dict[tuple[str, int], dict[str, Any]] = {}
    for row in bundle["averaged_nuclei"]:
        element = str(row.get("element", ""))
        if element not in {"H", "C"}:
            continue
        index = int(row["nucleus_index_0based"])
        value = row.get(shift_field)
        if value is None or not math.isfinite(float(value)):
            raise WaterSensitivityError(
                f"Missing finite {shift_field} for {element}{index + 1}"
            )
        key = (element, index)
        if key in result:
            raise WaterSensitivityError(f"Duplicate analyte nucleus: {element}{index + 1}")
        result[key] = row
    if not result:
        raise WaterSensitivityError("Prediction bundle contains no comparable 1H/13C nuclei")
    return result


def solution_equivalent_shifts(
    nuclei: dict[tuple[str, int], dict[str, Any]],
    shift_field: str,
    rules: dict[str, Any],
) -> dict[tuple[str, int], float]:
    shifts = {key: float(row[shift_field]) for key, row in nuclei.items()}
    for group in rules.get("groups", []):
        element = str(group.get("element", ""))
        members = [(element, int(index)) for index in group.get("atom_indices_0based", [])]
        if not members or any(member not in shifts for member in members):
            raise WaterSensitivityError(
                f"Solution-equivalence group is incompatible with analyte map: {group.get('group_id')}"
            )
        mean = sum(shifts[member] for member in members) / len(members)
        for member in members:
            shifts[member] = mean
    return shifts


def exchangeable_order_changed(
    dry_rows: dict[tuple[str, int], dict[str, Any]],
    dry_shifts: dict[tuple[str, int], float],
    wet_shifts: dict[tuple[str, int], float],
) -> bool:
    indices = sorted(
        index
        for (element, index), row in dry_rows.items()
        if element == "H" and int(row.get("exchangeable_H_flag", 0)) == 1
    )
    for position, first in enumerate(indices):
        for second in indices[position + 1 :]:
            dry_difference = dry_shifts[("H", first)] - dry_shifts[("H", second)]
            wet_difference = wet_shifts[("H", first)] - wet_shifts[("H", second)]
            if dry_difference * wet_difference < 0.0:
                return True
    return False


def compare(
    dry_path: Path,
    wet_path: Path,
    config_path: Path,
    output_directory: Path,
    profile_name: str | None = None,
    shift_variant: str = "raw_tms",
    solution_equivalence_path: Path = DEFAULT_SOLUTION_EQUIVALENCE,
) -> dict[str, Any]:
    config = load_json(config_path)
    selected_name, profile = selected_profile(config, profile_name)
    if shift_variant not in SHIFT_FIELDS:
        raise WaterSensitivityError(f"Unsupported shift variant: {shift_variant}")
    shift_field = SHIFT_FIELDS[shift_variant]

    dry_bundle = load_json(dry_path)
    wet_bundle = load_json(wet_path)
    dry_figure, dry_scenario, dry_identity = bundle_identity(dry_bundle, dry_path)
    wet_figure, wet_scenario, wet_identity = bundle_identity(wet_bundle, wet_path)
    if (
        dry_figure != wet_figure
        or not dry_identity
        or not wet_identity
        or dry_identity != wet_identity
    ):
        raise WaterSensitivityError("Dry and wet bundles do not describe the same compound")
    dry_signature = scientific_signature(dry_bundle, dry_path)
    wet_signature = scientific_signature(wet_bundle, wet_path)
    if dry_signature != wet_signature:
        raise WaterSensitivityError(
            "Dry/wet method, TMS/reference, regression, temperature or weighting policy differs"
        )

    equivalence = load_json(solution_equivalence_path)
    compound_rules = equivalence.get("compounds", {}).get(dry_figure)
    if not isinstance(compound_rules, dict):
        raise WaterSensitivityError(f"No solution-equivalence policy for {dry_figure}")
    expected_identity = str(compound_rules.get("inchikey", ""))
    if not expected_identity or dry_identity != expected_identity:
        raise WaterSensitivityError(
            f"Solution-equivalence identity mismatch for {dry_figure}: {dry_identity}"
        )

    policy = profile.get("fixed_composition_policy", {})
    pairs = policy.get("matched_water_delta_pairs", {}) if isinstance(policy, dict) else {}
    expected = pairs.get(dry_figure) if isinstance(pairs, dict) else None
    if not isinstance(expected, dict):
        raise WaterSensitivityError(f"No declared matched dry/wet pair for {dry_figure}")
    if dry_scenario != expected.get("dry") or wet_scenario != expected.get("wet"):
        raise WaterSensitivityError(
            f"Scenario mismatch for {dry_figure}: expected {expected.get('dry')} versus "
            f"{expected.get('wet')}; received {dry_scenario} versus {wet_scenario}"
        )
    if policy.get("cross_stoichiometry_boltzmann_mixing") != "FORBIDDEN":
        raise WaterSensitivityError("Workflow profile does not forbid cross-stoichiometry mixing")

    dry_nuclei = analyte_nuclei(dry_bundle, shift_field)
    wet_nuclei = analyte_nuclei(wet_bundle, shift_field)
    if set(dry_nuclei) != set(wet_nuclei):
        raise WaterSensitivityError("Dry/wet analyte nucleus maps differ")
    dry_shifts = solution_equivalent_shifts(dry_nuclei, shift_field, compound_rules)
    wet_shifts = solution_equivalent_shifts(wet_nuclei, shift_field, compound_rules)

    thresholds = profile.get("adaptive_escalation", {}).get("water_sensitivity_flags", {})
    nonexchange_h = float(thresholds["nonexchangeable_1H_abs_delta_ppm"])
    exchange_h = float(thresholds["exchangeable_1H_abs_delta_ppm"])
    carbon = float(thresholds["13C_abs_delta_ppm"])
    order_change = exchangeable_order_changed(dry_nuclei, dry_shifts, wet_shifts)

    rows: list[dict[str, Any]] = []
    for element, index in sorted(dry_nuclei, key=lambda item: (item[0] != "H", item[1])):
        dry_row = dry_nuclei[(element, index)]
        wet_row = wet_nuclei[(element, index)]
        if str(dry_row.get("element")) != str(wet_row.get("element")):
            raise WaterSensitivityError(f"Element mismatch at analyte index {index}")
        exchangeable = element == "H" and int(dry_row.get("exchangeable_H_flag", 0)) == 1
        if exchangeable != (
            element == "H" and int(wet_row.get("exchangeable_H_flag", 0)) == 1
        ):
            raise WaterSensitivityError(f"Exchangeable-H identity mismatch at H{index + 1}")
        threshold = carbon if element == "C" else exchange_h if exchangeable else nonexchange_h
        dry_shift = dry_shifts[(element, index)]
        wet_shift = wet_shifts[(element, index)]
        delta = wet_shift - dry_shift
        flagged = abs(delta) > threshold or (exchangeable and order_change)
        rows.append(
            {
                "notice": NOTICE,
                "workflow_profile": selected_name,
                "figure_id": dry_figure,
                "dry_scenario_id": dry_scenario,
                "wet_scenario_id": wet_scenario,
                "shift_variant": shift_variant,
                "shift_field": shift_field,
                "nucleus": "1H" if element == "H" else "13C",
                "nucleus_index_0based": index,
                "atom_label_1based": str(dry_row.get("atom_label_1based", f"{element}{index + 1}")),
                "exchangeable_H_flag": int(exchangeable),
                "dry_shift_ppm": dry_shift,
                "wet_shift_ppm": wet_shift,
                "wet_minus_dry_delta_ppm": delta,
                "absolute_delta_ppm": abs(delta),
                "declared_threshold_ppm": threshold,
                "exchangeable_H_order_change": int(exchangeable and order_change),
                "environment_sensitive": int(flagged),
                "interpretation": (
                    "ENVIRONMENT_SENSITIVE" if flagged else "WITHIN_DECLARED_SENSITIVITY_THRESHOLD"
                ),
                "cross_stoichiometry_weight": "NONE",
            }
        )

    status = (
        "ENVIRONMENT_SENSITIVE"
        if any(int(row["environment_sensitive"]) for row in rows)
        else "WITHIN_DECLARED_SENSITIVITY_THRESHOLDS"
    )
    stem = f"{safe_token(dry_figure)}__trace_water_sensitivity__{safe_token(shift_variant)}"
    csv_path = output_directory / f"{stem}.csv"
    json_path = output_directory / f"{stem}.json"
    atomic_write_csv(csv_path, rows)
    report = {
        "schema_version": "1.0",
        "notice": NOTICE,
        "status": status,
        "workflow_profile": selected_name,
        "workflow_profile_config": str(config_path),
        "workflow_profile_config_sha256": sha256(config_path),
        "scientific_signature": dry_signature,
        "solution_equivalence_config": str(solution_equivalence_path),
        "solution_equivalence_config_sha256": sha256(solution_equivalence_path),
        "solution_equivalence_policy": "APPLIED_IDENTICALLY_TO_MATCHED_DRY_AND_WET_BUNDLES",
        "figure_id": dry_figure,
        "dry_scenario_id": dry_scenario,
        "wet_scenario_id": wet_scenario,
        "shift_variant": shift_variant,
        "dry_bundle": str(dry_path),
        "dry_bundle_sha256": sha256(dry_path),
        "wet_bundle": str(wet_path),
        "wet_bundle_sha256": sha256(wet_path),
        "n_compared_nuclei": len(rows),
        "n_environment_sensitive": sum(int(row["environment_sensitive"]) for row in rows),
        "exchangeable_H_order_change": order_change,
        "cross_stoichiometry_boltzmann_mixing": "FORBIDDEN",
        "main_prediction_replacement": "FORBIDDEN",
        "appearance_water_peak_inference": "FORBIDDEN",
        "rows": rows,
    }
    atomic_write_json(json_path, report)
    return {"csv": csv_path, "json": json_path, "report": report}


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--dry-bundle", type=Path, required=True)
    value.add_argument("--wet-bundle", type=Path, required=True)
    value.add_argument("--workflow-profiles", type=Path, default=DEFAULT_CONFIG)
    value.add_argument("--workflow-profile")
    value.add_argument("--shift-variant", choices=sorted(SHIFT_FIELDS), default="raw_tms")
    value.add_argument(
        "--solution-equivalence",
        type=Path,
        default=DEFAULT_SOLUTION_EQUIVALENCE,
    )
    value.add_argument("--output-directory", type=Path, required=True)
    return value


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        result = compare(
            args.dry_bundle.resolve(),
            args.wet_bundle.resolve(),
            args.workflow_profiles.resolve(),
            args.output_directory.resolve(),
            args.workflow_profile,
            args.shift_variant,
            args.solution_equivalence.resolve(),
        )
    except (WaterSensitivityError, OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"Trace-water comparison failed: {exc}")
        return 2
    print(json.dumps({key: str(value) for key, value in result.items() if key != "report"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
