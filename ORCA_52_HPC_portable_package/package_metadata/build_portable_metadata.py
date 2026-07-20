#!/usr/bin/env python3
"""Normalize copied Figure 3 ORCA inputs for location-independent HPC use.

This script never reads from or writes to the active Windows calculation tree.
It operates only inside the portable package that contains this file.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
GOAT_ROOT = ROOT / "01_goat_explicit"
SOURCE_PROJECT_WIN = (
    r"E:\For_codex\C_N_coupling\Nature_sub\Figure3_finalized_reference_states"
    r"\ORCA_52_NMR_serial_high_quality"
)
SOURCE_STATE_ROOT_WIN = (
    r"E:\For_codex\C_N_coupling\Nature_sub\Figure3_finalized_reference_states"
)
SOURCE_WORKSPACE_WIN = r"E:\For_codex\C_N_coupling"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


def write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def portable_string(value: str) -> str:
    """Convert known Windows paths to package-relative or explicit placeholders."""
    if value.startswith(SOURCE_PROJECT_WIN):
        suffix = value[len(SOURCE_PROJECT_WIN) :].lstrip("\\/")
        return suffix.replace("\\", "/") or "."
    if value.startswith(SOURCE_STATE_ROOT_WIN):
        suffix = value[len(SOURCE_STATE_ROOT_WIN) :].lstrip("\\/")
        known_registry = {
            "Figure3_52_ORCA_reference_state_registry.csv",
            "Figure3_52_ORCA_reference_state_registry.json",
            "Figure3_52_ORCA_ordered_reference_states.sdf",
            "F3_016_017_019_021_023_reference_state_manifest.csv",
            "F3_016_017_019_021_023_reference_state_manifest.json",
        }
        name = Path(suffix.replace("\\", "/")).name
        if name in known_registry:
            return f"registry/{name}"
        return "${LEGACY_REFERENCE_STATE_ROOT_NOT_INCLUDED}/" + suffix.replace("\\", "/")
    if value.startswith(SOURCE_WORKSPACE_WIN):
        suffix = value[len(SOURCE_WORKSPACE_WIN) :].lstrip("\\/")
        return "${LOCAL_SOURCE_NOT_INCLUDED}/" + suffix.replace("\\", "/")
    # Template-only pseudo paths such as E:\...\state\F3_xxx are made relative.
    pseudo = re.match(r"^[A-Za-z]:\\\.\.\.\\(.+)$", value)
    if pseudo:
        return pseudo.group(1).replace("\\", "/")
    if re.match(r"^[A-Za-z]:[\\/]", value):
        return "${EXTERNAL_PATH_NOT_INCLUDED}/" + Path(value.replace("\\", "/")).name
    return value


def normalize_object(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: normalize_object(item) for key, item in value.items()}
    if isinstance(value, list):
        return [normalize_object(item) for item in value]
    if isinstance(value, str):
        return portable_string(value)
    return value


def normalize_json_files() -> int:
    changed = 0
    roots = [
        ROOT / "config",
        ROOT / "data",
        ROOT / "jobs",
        ROOT / "reference_jobs",
        ROOT / "real_spectrum_calibration",
        GOAT_ROOT,
    ]
    for base in roots:
        for path in sorted(base.rglob("*.json")):
            original = json.loads(path.read_text(encoding="utf-8-sig"))
            normalized = normalize_object(original)
            if normalized != original:
                write_json(path, normalized)
                changed += 1
    return changed


def bundle_sha(directory: Path) -> str:
    h = hashlib.sha256()
    for name in (
        "goat.inp",
        "start.xyz",
        "atom_map.csv",
        "isotope_map.csv",
        "structure_validation.json",
    ):
        h.update(name.encode("utf-8"))
        h.update((directory / name).read_bytes())
    return h.hexdigest().upper()


def rewrite_goat_manifest() -> tuple[int, int, int, str]:
    manifest = DATA / "GOAT_explicit_job_manifest.csv"
    with manifest.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1120:
        raise ValueError(f"Expected 1120 GOAT rows, found {len(rows)}")

    output_rows: list[dict[str, str]] = []
    for row in rows:
        figure_id = row["figure_id"].strip()
        scenario_id = row["scenario_id"].strip()
        start_index = int(row["start_index"])
        rel = Path("01_goat_explicit") / figure_id / scenario_id / f"start_{start_index:02d}"
        directory = ROOT / rel
        if not directory.is_dir():
            raise FileNotFoundError(directory)
        required = [
            "goat.inp",
            "start.xyz",
            "atom_map.csv",
            "isotope_map.csv",
            "structure_validation.json",
            "INPUT_BUNDLE_SHA256.txt",
        ]
        missing = [name for name in required if not (directory / name).is_file()]
        if missing:
            raise FileNotFoundError(f"{rel.as_posix()}: {missing}")
        input_hash = sha256(directory / "goat.inp")
        coordinate_hash = sha256(directory / "start.xyz")
        current_bundle = bundle_sha(directory)
        (directory / "INPUT_BUNDLE_SHA256.txt").write_text(
            current_bundle + "\n", encoding="ascii"
        )
        row.update(
            {
                "stage": "goat",
                "job_directory": rel.as_posix(),
                "portable_job_directory": rel.as_posix(),
                "relative_job_directory": rel.as_posix(),
                "input_file": "goat.inp",
                "coordinate_file": "start.xyz",
                "output_file": "goat.out",
                "goat_input_sha256": input_hash,
                "start_xyz_sha256": coordinate_hash,
                "input_bundle_sha256": current_bundle,
                "expected_input_sha256": input_hash,
                "expected_coordinate_sha256": coordinate_hash,
                "expected_input_bundle_sha256": current_bundle,
            }
        )
        output_rows.append(row)

    preferred = [
        "queue_order",
        "job_id",
        "figure_id",
        "compound_number",
        "display_name",
        "stage",
        "scenario_id",
        "primary_solvent",
        "continuum_keyword",
        "explicit_solvent_count",
        "start_index",
        "charge",
        "multiplicity",
        "preparation_seed",
        "job_directory",
        "portable_job_directory",
        "relative_job_directory",
        "input_file",
        "coordinate_file",
        "output_file",
        "goat_input_sha256",
        "start_xyz_sha256",
        "input_bundle_sha256",
        "expected_input_sha256",
        "expected_coordinate_sha256",
        "expected_input_bundle_sha256",
    ]
    extras = [key for key in output_rows[0] if key not in preferred]
    fields = preferred + extras
    for target in (manifest, DATA / "GOAT_explicit_job_manifest_portable.csv"):
        with target.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(output_rows)

    compounds = len({row["figure_id"] for row in output_rows})
    scenarios = len({(row["figure_id"], row["scenario_id"]) for row in output_rows})
    return len(output_rows), compounds, scenarios, sha256(manifest)


def rewrite_reference_manifest() -> None:
    path = DATA / "reference_job_manifest_DRAFT.csv"
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 6:
        raise ValueError(f"Expected six reference models, found {len(rows)}")
    fields = list(rows[0])
    for new_field in ("portable_job_directory", "stage"):
        if new_field not in fields:
            fields.append(new_field)
    for row in rows:
        rel = f"reference_jobs/{row['model_id']}"
        directory = ROOT / rel
        if not directory.is_dir():
            raise FileNotFoundError(directory)
        row["job_directory"] = rel
        row["portable_job_directory"] = rel
        row["stage"] = "reference_opt_freq_then_shielding"
        for field, filename in (
            ("start_xyz_sha256", "start.xyz"),
            ("reference_atom_map_sha256", "reference_atom_map.csv"),
            ("model_declaration_sha256", "model_declaration.json"),
            ("opt_input_sha256", row["opt_input_file"]),
            ("shielding_template_sha256", row["shielding_template_file"]),
        ):
            row[field] = sha256(directory / filename).lower()
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def update_summary(jobs: int, compounds: int, scenarios: int, manifest_hash: str) -> None:
    for name in ("explicit_preparation_summary.json", "GOAT_explicit_preparation_audit.json"):
        path = DATA / name
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
        payload["manifest"] = "data/GOAT_explicit_job_manifest.csv"
        payload["manifest_sha256"] = manifest_hash
        payload["jobs" if "jobs" in payload else "goat_jobs"] = jobs
        payload["compounds"] = compounds
        payload["scenarios"] = scenarios
        if "isotope_mass_policy" in payload:
            payload["isotope_mass_policy"] = "data/ORCA_6_1_1_isotope_mass_policy.json"
        if "isotope_mass_policy_file" in payload:
            payload["isotope_mass_policy_file"] = "data/ORCA_6_1_1_isotope_mass_policy.json"
        write_json(path, payload)


def refresh_reference_registry_hash() -> None:
    protocol = ROOT / "jobs" / "reference_model_protocol_DRAFT.json"
    registry = ROOT / "config" / "reference_job_registry.json"
    payload = json.loads(registry.read_text(encoding="utf-8-sig"))
    payload["reference_model_protocol_gate"] = "jobs/reference_model_protocol_DRAFT.json"
    payload["reference_model_protocol_gate_sha256"] = sha256(protocol).lower()
    write_json(registry, payload)


def rewrite_text_documents() -> int:
    replacements = [
        (SOURCE_PROJECT_WIN, "${PACKAGE_ROOT}"),
        (SOURCE_STATE_ROOT_WIN, "${LEGACY_REFERENCE_STATE_ROOT_NOT_INCLUDED}"),
        (r"D:\ORCA_6.1.1\xtb-6.7.1pre", "${XTB_DIR}"),
        (r"D:\ORCA_6.1.1\orca.exe", "${ORCA_EXE}"),
        (r"D:\ORCA_6.1.1", "${ORCA_ROOT}"),
        (r"C:\Program Files\Microsoft MPI\Bin", "${MPI_BIN}"),
        (SOURCE_WORKSPACE_WIN, "${LOCAL_SOURCE_NOT_INCLUDED}"),
    ]
    changed = 0
    for base in (ROOT / "docs", ROOT / "real_spectrum_calibration"):
        for path in base.rglob("*"):
            if not path.is_file() or path.suffix.lower() not in {".md", ".txt", ".json"}:
                continue
            text = path.read_text(encoding="utf-8-sig")
            new_text = text
            for old, new in replacements:
                new_text = new_text.replace(old, new)
            if new_text != text:
                path.write_text(new_text, encoding="utf-8")
                changed += 1
    return changed


def main() -> int:
    normalized_json = normalize_json_files()
    jobs, compounds, scenarios, manifest_hash = rewrite_goat_manifest()
    rewrite_reference_manifest()
    update_summary(jobs, compounds, scenarios, manifest_hash)
    refresh_reference_registry_hash()
    docs_changed = rewrite_text_documents()
    report = {
        "status": "PASS",
        "notice": "CALCULATED MODELS ONLY - NOT EXPERIMENTAL DATA",
        "package_root_policy": "all production paths are package-root relative",
        "goat_jobs": jobs,
        "compounds": compounds,
        "scenarios": scenarios,
        "goat_manifest": "data/GOAT_explicit_job_manifest.csv",
        "goat_manifest_sha256": sha256(DATA / "GOAT_explicit_job_manifest.csv"),
        "json_files_normalized": normalized_json,
        "documentation_files_rewritten": docs_changed,
        "source_project_modified": False,
    }
    write_json(ROOT / "package_metadata" / "PORTABLE_REWRITE_PROVENANCE.json", report)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
