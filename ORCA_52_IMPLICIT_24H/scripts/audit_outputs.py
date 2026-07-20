#!/usr/bin/env python3
"""Strict completeness and provenance audit for the 24 h prediction bundle."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from common import (
    ROOT,
    artifact_manifest_issues,
    frozen_geometry_issues,
    read_csv,
    read_banner_csv,
    read_json,
    sha256,
    write_csv,
    write_json,
)


def line_count(path: Path) -> int:
    with path.open("rb") as handle:
        return sum(1 for _ in handle)


def sum_observed_atoms(peak_path: Path) -> int:
    rows = read_csv(peak_path)
    if not rows or "integral_or_equivalent_atom_count" not in rows[0]:
        return 0
    return sum(int(float(row["integral_or_equivalent_atom_count"])) for row in rows)


def portable_artifact_path(run_dir: Path, declared: Any, inferred: Path) -> Path | None:
    """Resolve a result artifact after the result directory has been moved.

    Canonical paths inferred below the current run directory deliberately take
    precedence over absolute paths recorded on the machine that produced the
    bundle.
    """

    candidates = [inferred]
    text = str(declared or "").strip()
    if text and text != "NOT_APPLICABLE":
        recorded = Path(text)
        if not recorded.is_absolute():
            candidates.append(run_dir / recorded)
        candidates.append(recorded)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def raw_output_path(run_dir: Path, rows: list[dict[str, str]]) -> Path | None:
    """Resolve the raw ORCA output, preferring its run-relative identity."""

    if not rows:
        return None
    first = rows[0]
    relative = first.get("source_output_relative_to_run_dir", "").strip()
    expected_hash = first.get("source_output_sha256", "").strip().upper()
    if relative:
        candidate = run_dir / Path(relative)
        if candidate.is_file():
            return candidate

    # A relative path is normally sufficient.  The hash search is a migration
    # fallback for an older bundle whose run directory was reorganized.
    basename = Path(first.get("source_output", "")).name
    search_root = run_dir / "raw_orca"
    if basename and search_root.is_dir() and expected_hash:
        for candidate in search_root.rglob(basename):
            if candidate.is_file() and sha256(candidate) == expected_hash:
                return candidate
    recorded = Path(first.get("source_output", ""))
    return recorded if recorded.is_file() else None


def parsed_shielding_issues(
    *,
    run_dir: Path,
    parsed_path: Path,
    run_identity: dict[str, Any],
    label: str,
    expected_figure_id: str,
) -> list[str]:
    """Audit parsed shielding provenance back to retained ORCA inputs."""

    issues: list[str] = []
    if not parsed_path.is_file():
        return [f"{label}_shielding_missing"]
    rows = read_csv(parsed_path)
    if not rows:
        return [f"{label}_shielding_empty"]

    expected_uuid = str(run_identity.get("run_uuid", ""))
    expected_fingerprint = str(run_identity.get("workflow_fingerprint", ""))
    required = (
        "run_uuid",
        "workflow_fingerprint",
        "source_output_relative_to_run_dir",
        "source_output_sha256",
        "source_input_sha256",
        "source_geometry_sha256",
    )
    for field in required:
        if any(not row.get(field, "").strip() for row in rows):
            issues.append(f"{label}_parsed_missing:{field}")
    if not expected_uuid or not expected_fingerprint:
        issues.append("run_config_identity_missing")
    if expected_uuid and any(row.get("run_uuid") != expected_uuid for row in rows):
        issues.append(f"{label}_run_uuid_mismatch")
    if expected_fingerprint and any(
        row.get("workflow_fingerprint") != expected_fingerprint for row in rows
    ):
        issues.append(f"{label}_workflow_fingerprint_mismatch")
    if any(row.get("figure_id") != expected_figure_id for row in rows):
        issues.append(f"{label}_figure_id_mismatch")

    # Every row must describe one and the same retained ORCA calculation.
    for field in required:
        values = {row.get(field, "") for row in rows}
        if len(values) != 1:
            issues.append(f"{label}_parsed_inconsistent:{field}")

    raw_output = raw_output_path(run_dir, rows)
    first = rows[0]
    if raw_output is None:
        issues.append(f"{label}_raw_output_missing")
        return issues
    if sha256(raw_output) != first.get("source_output_sha256", "").upper():
        issues.append(f"{label}_raw_output_hash_mismatch")

    record_path = raw_output.parent / "run_record.json"
    if not record_path.is_file():
        issues.append(f"{label}_run_record_missing")
        return issues
    record = read_json(record_path)
    input_hash = first.get("source_input_sha256", "").upper()
    geometry_hash = first.get("source_geometry_sha256", "").upper()
    if str(record.get("input_sha256", "")).upper() != input_hash:
        issues.append(f"{label}_run_record_input_hash_mismatch")
    if str(record.get("geometry_sha256", "")).upper() != geometry_hash:
        issues.append(f"{label}_run_record_geometry_hash_mismatch")
    if str(record.get("status", "")) != "COMPLETE":
        issues.append(f"{label}_run_record_not_complete")

    input_matches = [path for path in raw_output.parent.glob("*.inp") if sha256(path) == input_hash]
    if not input_matches:
        issues.append(f"{label}_retained_input_hash_mismatch")
    geometry_matches = [path for path in raw_output.parent.glob("*.xyz") if sha256(path) == geometry_hash]
    if not geometry_matches:
        issues.append(f"{label}_retained_geometry_hash_mismatch")
    return issues


def peak_line_conservation_issues(folder: Path, nucleus: str) -> tuple[list[str], float]:
    """Check first-order line area and centroid against each predicted signal."""

    issues: list[str] = []
    peaks = read_csv(folder / f"{nucleus}_predicted_peaks.csv")
    lines = read_csv(folder / f"{nucleus}_first_order_lines.csv")
    if not peaks or "signal_index" not in peaks[0]:
        return issues, 0.0
    if not lines or "signal_index" not in lines[0]:
        return [f"line_table_missing_for_analyte:{nucleus}"], sum(
            float(row["integral_or_equivalent_atom_count"]) for row in peaks
        )

    peak_by_signal = {int(row["signal_index"]): row for row in peaks}
    lines_by_signal: dict[int, list[dict[str, str]]] = {}
    for line in lines:
        lines_by_signal.setdefault(int(line["signal_index"]), []).append(line)
    unexpected = sorted(set(lines_by_signal) - set(peak_by_signal))
    if unexpected:
        issues.append(f"unexpected_line_signal:{nucleus}:{','.join(map(str, unexpected))}")

    peak_total = sum(float(row["integral_or_equivalent_atom_count"]) for row in peaks)
    line_total = sum(float(row["relative_area"]) for row in lines)
    area_tolerance = max(1.0e-7, 1.0e-6 * max(1.0, peak_total))
    if not math.isclose(line_total, peak_total, rel_tol=1.0e-6, abs_tol=area_tolerance):
        issues.append(f"line_peak_total_area_mismatch:{nucleus}:{line_total:.9g}!={peak_total:.9g}")

    for signal_id, peak in peak_by_signal.items():
        signal_lines = lines_by_signal.get(signal_id, [])
        if not signal_lines:
            issues.append(f"signal_lines_missing:{nucleus}:{signal_id}")
            continue
        expected_area = float(peak["integral_or_equivalent_atom_count"])
        observed_area = sum(float(line["relative_area"]) for line in signal_lines)
        if not math.isclose(observed_area, expected_area, rel_tol=1.0e-6, abs_tol=1.0e-7):
            issues.append(
                f"signal_area_mismatch:{nucleus}:{signal_id}:{observed_area:.9g}!={expected_area:.9g}"
            )
        if observed_area <= 0:
            issues.append(f"signal_nonpositive_area:{nucleus}:{signal_id}")
            continue
        centroid = sum(
            float(line["ppm"]) * float(line["relative_area"]) for line in signal_lines
        ) / observed_area
        delta = float(peak["delta_ppm"])
        if not math.isclose(centroid, delta, rel_tol=0.0, abs_tol=2.0e-5):
            issues.append(f"signal_centroid_mismatch:{nucleus}:{signal_id}:{centroid:.8f}!={delta:.8f}")
    return issues, peak_total


def audit_compound(run_dir: Path, row: dict[str, str]) -> dict[str, Any]:
    figure_id = row["figure_id"]
    folder = run_dir / "NMR_CALCULATED_NOT_EXPERIMENTAL" / figure_id
    metadata_path = folder / "metadata.json"
    issues: list[str] = []
    if not metadata_path.is_file():
        return {"figure_id": figure_id, "status": "FAIL", "issues": ["metadata_missing"]}
    metadata = read_json(metadata_path)
    run_config_path = run_dir / "run_config.json"
    run_identity = read_json(run_config_path) if run_config_path.is_file() else {}
    if not run_identity:
        issues.append("run_config_missing_or_empty")
    elif not run_identity.get("runtime_environment_identity"):
        issues.append("runtime_environment_identity_missing")
    if not str(metadata.get("status", "")).startswith("COMPLETE"):
        issues.append(f"metadata_status:{metadata.get('status')}")
    if metadata.get("prediction_label") != "CALCULATED/SYNTHETIC - NOT EXPERIMENTAL":
        issues.append("prediction_label_missing_or_changed")

    expected_counts = {
        "1H": int(row["observable_1H_count_after_exchange_policy"]),
        "13C": int(row["observable_13C_count"]),
    }
    expected_points = {"1H": 65536, "13C": 32768}
    for nucleus in ("1H", "13C"):
        required = [
            folder / f"{nucleus}_spectrum_full_points.csv",
            folder / f"{nucleus}_spectrum_full_points.txt",
            folder / f"{nucleus}_predicted_peaks.csv",
            folder / f"{nucleus}_first_order_lines.csv",
            folder / f"{nucleus}_clean.png",
            folder / f"{nucleus}_clean.svg",
            folder / f"{nucleus}_synthetic_appearance.png",
            folder / f"{nucleus}_synthetic_appearance.svg",
        ]
        for path in required:
            if not path.is_file() or path.stat().st_size == 0:
                issues.append(f"missing_or_empty:{path.name}")
        csv_path = folder / f"{nucleus}_spectrum_full_points.csv"
        if csv_path.is_file() and line_count(csv_path) != expected_points[nucleus] + 1:
            issues.append(
                f"point_count:{nucleus}:{line_count(csv_path) - 1}!={expected_points[nucleus]}"
            )
        peak_path = folder / f"{nucleus}_predicted_peaks.csv"
        if peak_path.is_file() and sum_observed_atoms(peak_path) != expected_counts[nucleus]:
            issues.append(
                f"observable_atom_count:{nucleus}:{sum_observed_atoms(peak_path)}!={expected_counts[nucleus]}"
            )
        line_path = folder / f"{nucleus}_first_order_lines.csv"
        if peak_path.is_file() and line_path.is_file():
            conservation, analyte_area = peak_line_conservation_issues(folder, nucleus)
            issues.extend(conservation)
            if analyte_area > 0:
                profile = metadata.get("profiles", {}).get(nucleus, {})
                try:
                    capture = float(profile["clean_area_capture_fraction"])
                except (KeyError, TypeError, ValueError):
                    issues.append(f"clean_area_capture_fraction_missing:{nucleus}")
                else:
                    # A unit-area Lorentzian is sampled on a finite digital grid.
                    # The 13C grid is only about one point per 1 Hz clean FWHM, so
                    # quadrature alone can deviate by roughly +/-2.4%.  A symmetric
                    # 5% guard still catches clipping or a broken normalization.
                    if not 0.95 <= capture <= 1.05:
                        issues.append(f"clean_area_capture_fraction_out_of_range:{nucleus}:{capture:.8g}")
    workbook = folder / "NMR_prediction_data.xlsx"
    if not workbook.is_file() or workbook.stat().st_size == 0:
        issues.append("xlsx_missing_or_empty")

    geometry_value = metadata.get("geometry")
    if geometry_value and geometry_value != "NOT_APPLICABLE":
        geometry_tier = metadata.get("selected_geometry_tier")
        if geometry_tier == "REFINED_GEOMETRY":
            inferred_geometry = run_dir / "geometries" / "refined" / f"{figure_id}.xyz"
        elif geometry_tier == "BASELINE_GEOMETRY":
            inferred_geometry = run_dir / "geometries" / "baseline" / f"{figure_id}.xyz"
        else:
            inferred_geometry = ROOT / "structures" / figure_id / "start.xyz"
        geometry = portable_artifact_path(run_dir, geometry_value, inferred_geometry)
        if geometry is None:
            issues.append("selected_geometry_missing")
        elif sha256(geometry) != metadata.get("geometry_sha256"):
            issues.append("selected_geometry_hash_mismatch")
        else:
            issues.extend(
                frozen_geometry_issues(
                    ROOT / "structures" / figure_id / "start.xyz",
                    geometry,
                    ROOT / "structures" / figure_id / "bonds.csv",
                    row["model_smiles"],
                )
            )

    if expected_counts["1H"] > 0 or expected_counts["13C"] > 0:
        method_id = str(metadata.get("selected_method_id", ""))
        selected_inferred = (
            run_dir / "parsed_shieldings" / "compounds" / figure_id / f"{method_id}.csv"
        )
        shield_path = portable_artifact_path(
            run_dir, metadata.get("selected_shielding_file"), selected_inferred
        )
        if shield_path is None:
            issues.append("selected_shielding_missing")
        elif sha256(shield_path) != metadata.get("selected_shielding_sha256"):
            issues.append("selected_shielding_hash_mismatch")
        else:
            issues.extend(
                parsed_shielding_issues(
                    run_dir=run_dir,
                    parsed_path=shield_path,
                    run_identity=run_identity,
                    label="selected",
                    expected_figure_id=figure_id,
                )
            )
        if not metadata.get("reference_shieldings_ppm"):
            issues.append("reference_shielding_missing")
        reference_id = (
            "TMS_DMSO" if row["orca_continuum_keyword"] == "DMSO" else "TSP_MINUS_WATER"
        )
        reference_inferred = (
            run_dir / "parsed_shieldings" / "references" / reference_id / f"{method_id}.csv"
        )
        reference_path = portable_artifact_path(
            run_dir, metadata.get("reference_shielding_file"), reference_inferred
        )
        if reference_path is None:
            issues.append("reference_shielding_file_missing")
        elif sha256(reference_path) != metadata.get("reference_shielding_sha256"):
            issues.append("reference_shielding_hash_mismatch")
        else:
            issues.extend(
                parsed_shielding_issues(
                    run_dir=run_dir,
                    parsed_path=reference_path,
                    run_identity=run_identity,
                    label="reference",
                    expected_figure_id=reference_id,
                )
            )
        if metadata.get("proton_splitting_model") != "EMPIRICAL GEOMETRY/TOPOLOGY J - NOT ORCA J":
            issues.append("empirical_J_provenance_missing")

    atom_map = read_csv(ROOT / "structures" / figure_id / "atom_map.csv")
    if any(atom["environment_atom"].lower() == "true" or atom["component"] != "analyte" for atom in atom_map):
        issues.append("explicit_environment_atom_found")
    return {
        "figure_id": figure_id,
        "display_name": row["display_name"],
        "method_id": metadata.get("selected_method_id", ""),
        "geometry_tier": metadata.get("selected_geometry_tier", ""),
        "status": "PASS" if not issues else "FAIL",
        "issues": issues,
    }


def audit_hrms(run_dir: Path, expected_ids: set[str]) -> tuple[set[str], list[str]]:
    output = run_dir / "HRMS_CALCULATED_NOT_EXPERIMENTAL"
    summary = output / "Figure3_formula_HRMS_prediction_summary.csv"
    validation_path = output / "VALIDATION_AND_PROVENANCE.json"
    manifest_path = output / "FILE_MANIFEST_SHA256.csv"
    issues: list[str] = []
    if not summary.is_file():
        observed_ids: set[str] = set()
        issues.append("summary_missing")
    else:
        observed_ids = {row["figure_id"] for row in read_banner_csv(summary)}
        if observed_ids != expected_ids:
            issues.append("summary_id_set_mismatch")
    if not validation_path.is_file():
        return observed_ids, [*issues, "validation_missing"]
    audit = read_json(validation_path)
    expected_count = len(expected_ids)
    if audit.get("prediction_status") != "CALCULATED HRMS PREDICTION - NOT EXPERIMENTAL":
        issues.append("prediction_label_invalid")
    if audit.get("expected_compound_count") != expected_count:
        issues.append("expected_compound_count_mismatch")
    if audit.get("observed_compound_count") != expected_count:
        issues.append("observed_compound_count_mismatch")
    if audit.get("missing_ids") or audit.get("extra_ids"):
        issues.append("missing_or_extra_ids")
    if audit.get("all_self_tests_pass") is not True:
        issues.append("element_self_tests_failed")
    crosscheck = audit.get("independent_reference_mass_crosscheck", {})
    if crosscheck.get("pass") is not True:
        issues.append("independent_mass_crosscheck_failed")
    if crosscheck.get("comparison_count") != expected_count:
        issues.append("independent_mass_crosscheck_count_mismatch")
    if audit.get("fragment_ions_generated") is not False:
        issues.append("fragment_scope_flag_invalid")
    if audit.get("experimental_peaks_generated") is not False:
        issues.append("experimental_scope_flag_invalid")
    if audit.get("input_sdf_sha256") != sha256(
        ROOT / "registry" / "Figure3_52_ORCA_ordered_reference_states.sdf"
    ):
        issues.append("input_sdf_hash_mismatch")
    if audit.get("input_registry_sha256") != sha256(
        ROOT / "registry" / "Figure3_52_ORCA_reference_state_registry.csv"
    ):
        issues.append("input_registry_hash_mismatch")
    if audit.get("pinned_engine_sha256") != sha256(
        ROOT / "scripts" / "generate_formula_hrms.py"
    ):
        issues.append("pinned_engine_hash_mismatch")
    issues.extend(
        f"internal_{item}"
        for item in artifact_manifest_issues(
            manifest_path, workflow_root=ROOT, generated_root=output
        )
    )
    return observed_ids, issues


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--only", default="")
    args = parser.parse_args()
    species = read_csv(ROOT / "config" / "primary_species.csv")
    selected = {item.strip() for item in args.only.split(",") if item.strip()}
    species = [row for row in species if not selected or row["figure_id"] in selected]
    reports = [audit_compound(args.run_dir, row) for row in species]
    expected_hrms_ids = {row["figure_id"] for row in species}
    observed_hrms, hrms_global_issues = audit_hrms(args.run_dir, expected_hrms_ids)
    for report in reports:
        if report["figure_id"] not in observed_hrms:
            report["issues"].append("hrms_candidate_ion_missing")
            report["status"] = "FAIL"
        if hrms_global_issues:
            report["issues"].extend(f"hrms_{issue}" for issue in hrms_global_issues)
            report["status"] = "FAIL"

    flattened = []
    for report in reports:
        flattened.append(
            {
                "figure_id": report["figure_id"],
                "display_name": report["display_name"],
                "method_id": report["method_id"],
                "geometry_tier": report["geometry_tier"],
                "status": report["status"],
                "issues": ";".join(report["issues"]),
            }
        )
    write_csv(
        args.run_dir / "FINAL_COMPOUND_AUDIT.csv", flattened,
        ["figure_id", "display_name", "method_id", "geometry_tier", "status", "issues"],
    )
    failures = [report["figure_id"] for report in reports if report["status"] != "PASS"]
    tier_counts: dict[str, int] = {}
    for report in reports:
        tier_counts[report["method_id"]] = tier_counts.get(report["method_id"], 0) + 1
    summary = {
        "schema_version": "1.0",
        "prediction_only": True,
        "experimental_data_claimed": False,
        "expected_compound_count": len(species),
        "pass_count": len(species) - len(failures),
        "failures": failures,
        "method_tier_counts": tier_counts,
        "hrms_global_issues": hrms_global_issues,
        "status": "PASS" if not failures else "FAIL",
    }
    write_json(args.run_dir / "FINAL_AUDIT.json", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
