#!/usr/bin/env python3
"""Read-only scientific and portability audit for the HPC upload package."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "package_metadata"
WINDOWS_ABSOLUTE = re.compile(r"(?i)(?:^|[^A-Za-z0-9_])[A-Za-z]:[\\/]")
REQUIRED_INPUT_FILES = {
    "goat.inp",
    "start.xyz",
    "atom_map.csv",
    "isotope_map.csv",
    "structure_validation.json",
    "INPUT_BUNDLE_SHA256.txt",
}
RUNTIME_TEXT_SUFFIXES = {
    ".csv",
    ".json",
    ".py",
    ".sh",
    ".slurm",
    ".sbatch",
    ".yaml",
    ".yml",
    ".toml",
}
GENERATED_METADATA = {
    "PACKAGE_SHA256SUMS.txt",
    "PACKAGE_INVENTORY.csv",
    "PORTABILITY_AUDIT.json",
    "PORTABILITY_AUDIT.md",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


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


def is_relative_portable(value: str) -> bool:
    if not value or WINDOWS_ABSOLUTE.search(value) or value.startswith(("/", "\\")):
        return False
    parts = Path(value).parts
    return ".." not in parts


def audit_manifest(issues: list[dict[str, Any]]) -> dict[str, Any]:
    path = ROOT / "data" / "GOAT_explicit_job_manifest.csv"
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    required_columns = {
        "job_id",
        "figure_id",
        "scenario_id",
        "stage",
        "job_directory",
        "portable_job_directory",
        "input_file",
        "output_file",
        "expected_input_sha256",
        "expected_coordinate_sha256",
        "expected_input_bundle_sha256",
    }
    missing_columns = sorted(required_columns - set(rows[0] if rows else []))
    if missing_columns:
        issues.append({"scope": "manifest", "problem": "missing columns", "details": missing_columns})
    seen_jobs: set[str] = set()
    compounds: set[str] = set()
    scenarios: set[tuple[str, str]] = set()
    solvents: Counter[str] = Counter()
    for line, row in enumerate(rows, start=2):
        job_id = row.get("job_id", "")
        if job_id in seen_jobs:
            issues.append({"scope": job_id, "problem": "duplicate job_id", "line": line})
        seen_jobs.add(job_id)
        compound = row.get("figure_id", "")
        scenario = row.get("scenario_id", "")
        compounds.add(compound)
        scenarios.add((compound, scenario))
        solvents[row.get("primary_solvent", "UNKNOWN")] += 1
        rel = row.get("portable_job_directory") or row.get("job_directory", "")
        if not is_relative_portable(rel):
            issues.append({"scope": job_id, "problem": "nonportable job directory", "value": rel})
            continue
        directory = ROOT / rel
        if not directory.is_dir():
            issues.append({"scope": job_id, "problem": "job directory missing", "value": rel})
            continue
        names = {path.name for path in directory.iterdir() if path.is_file()}
        missing = sorted(REQUIRED_INPUT_FILES - names)
        if missing:
            issues.append({"scope": job_id, "problem": "input bundle incomplete", "details": missing})
            continue
        checks = {
            "expected_input_sha256": sha256(directory / row.get("input_file", "goat.inp")),
            "expected_coordinate_sha256": sha256(directory / row.get("coordinate_file", "start.xyz")),
            "expected_input_bundle_sha256": bundle_sha(directory),
        }
        for field, actual in checks.items():
            expected = row.get(field, "").upper()
            if expected != actual:
                issues.append(
                    {"scope": job_id, "problem": f"{field} mismatch", "expected": expected, "actual": actual}
                )
        lock_value = (directory / "INPUT_BUNDLE_SHA256.txt").read_text(encoding="ascii").strip().upper()
        if lock_value != checks["expected_input_bundle_sha256"]:
            issues.append({"scope": job_id, "problem": "bundle lock file mismatch"})
        text = (directory / "goat.inp").read_text(encoding="utf-8", errors="replace")
        if "start.xyz" not in text:
            issues.append({"scope": job_id, "problem": "goat.inp does not reference local start.xyz"})
    if len(rows) != 1120:
        issues.append({"scope": "manifest", "problem": "wrong job count", "actual": len(rows), "expected": 1120})
    if len(compounds) != 52:
        issues.append({"scope": "manifest", "problem": "wrong compound count", "actual": len(compounds), "expected": 52})
    if len(scenarios) != 112:
        issues.append({"scope": "manifest", "problem": "wrong scenario count", "actual": len(scenarios), "expected": 112})
    return {
        "jobs": len(rows),
        "unique_job_ids": len(seen_jobs),
        "compounds": len(compounds),
        "scenarios": len(scenarios),
        "solvent_jobs": dict(solvents),
        "manifest_sha256": sha256(path),
    }


def audit_references(issues: list[dict[str, Any]]) -> dict[str, Any]:
    path = ROOT / "data" / "reference_job_manifest_DRAFT.csv"
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        rel = row.get("portable_job_directory") or row.get("job_directory", "")
        if not is_relative_portable(rel):
            issues.append({"scope": row.get("model_id"), "problem": "nonportable reference directory", "value": rel})
            continue
        directory = ROOT / rel
        for name in ("start.xyz", "reference_atom_map.csv", "model_declaration.json", "opt_freq.inp", "shielding.inp.template"):
            if not (directory / name).is_file():
                issues.append({"scope": row.get("model_id"), "problem": "reference file missing", "value": name})
    if len(rows) != 6:
        issues.append({"scope": "references", "problem": "wrong reference model count", "actual": len(rows), "expected": 6})
    return {"reference_models": len(rows), "status": "DRAFT_NOT_APPROVED"}


def audit_runtime_paths(issues: list[dict[str, Any]]) -> dict[str, Any]:
    roots = ["config", "data", "jobs", "reference_jobs", "scripts", "postprocess", "tests", "hpc"]
    scanned = 0
    findings: list[dict[str, Any]] = []
    for name in roots:
        base = ROOT / name
        if not base.exists():
            if name != "hpc":
                issues.append({"scope": name, "problem": "required directory missing"})
            continue
        for path in base.rglob("*"):
            if not path.is_file() or path.suffix.lower() not in RUNTIME_TEXT_SUFFIXES:
                continue
            scanned += 1
            text = path.read_text(encoding="utf-8-sig", errors="replace")
            if WINDOWS_ABSOLUTE.search(text):
                findings.append({"file": path.relative_to(ROOT).as_posix(), "problem": "Windows absolute path"})
    issues.extend({"scope": item["file"], "problem": item["problem"]} for item in findings)
    return {"runtime_text_files_scanned": scanned, "absolute_path_findings": len(findings)}


def audit_forbidden_artifacts(issues: list[dict[str, Any]]) -> dict[str, Any]:
    forbidden_dirs = {"__pycache__", "logs", "scratch", "state", "02_goat_explicit_runs"}
    forbidden_suffixes = {".gbw", ".tmp", ".lock", ".pyc", ".hess", ".engrad"}
    forbidden_names = {"orca.exe", "xtb.exe", "mpiexec.exe"}
    found: list[str] = []
    for path in ROOT.rglob("*"):
        rel = path.relative_to(ROOT)
        if any(part in forbidden_dirs for part in rel.parts):
            found.append(rel.as_posix())
        elif path.is_file() and (path.suffix.lower() in forbidden_suffixes or path.name.lower() in forbidden_names):
            found.append(rel.as_posix())
    for rel in sorted(set(found)):
        issues.append({"scope": rel, "problem": "forbidden cache/output/binary artifact"})
    return {"forbidden_artifacts": len(set(found))}


def write_inventory_and_hashes() -> dict[str, Any]:
    inventory_path = META / "PACKAGE_INVENTORY.csv"
    hashes_path = META / "PACKAGE_SHA256SUMS.txt"
    excluded = {inventory_path.resolve(), hashes_path.resolve()}
    files = [path for path in ROOT.rglob("*") if path.is_file() and path.resolve() not in excluded]
    files.sort(key=lambda path: path.relative_to(ROOT).as_posix())
    with inventory_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["relative_path", "size_bytes", "sha256"], lineterminator="\n")
        writer.writeheader()
        for path in files:
            writer.writerow(
                {
                    "relative_path": path.relative_to(ROOT).as_posix(),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    # Include the inventory itself; exclude only the self-referential checksum file.
    checksum_files = sorted(
        (path for path in ROOT.rglob("*") if path.is_file() and path.resolve() != hashes_path.resolve()),
        key=lambda path: path.relative_to(ROOT).as_posix(),
    )
    hashes_path.write_text(
        "".join(f"{sha256(path)}  {path.relative_to(ROOT).as_posix()}\n" for path in checksum_files),
        encoding="utf-8",
    )
    return {
        "inventory_entries_excluding_generated_inventory_and_checksum": len(files),
        "checksum_entries_excluding_checksum_self": len(checksum_files),
        "total_bytes_excluding_checksum": sum(path.stat().st_size for path in checksum_files),
    }


def main() -> int:
    issues: list[dict[str, Any]] = []
    report = {
        "audit_time_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "notice": "CALCULATED MODELS ONLY - NOT EXPERIMENTAL DATA",
        "audit_mode": "read-only inspection of scientific/runtime files; only audit metadata is written",
        "goat": audit_manifest(issues),
        "references": audit_references(issues),
        "runtime_paths": audit_runtime_paths(issues),
        "artifact_exclusions": audit_forbidden_artifacts(issues),
    }
    report["issue_count"] = len(issues)
    report["issues"] = issues
    report["status"] = "PASS" if not issues else "FAIL"
    META.mkdir(parents=True, exist_ok=True)
    json_path = META / "PORTABILITY_AUDIT.json"
    json_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    md = [
        "# HPC portable package audit",
        "",
        f"- Status: **{report['status']}**",
        f"- GOAT jobs: {report['goat']['jobs']}",
        f"- Compounds: {report['goat']['compounds']}",
        f"- Fixed scenarios: {report['goat']['scenarios']}",
        f"- Reference models: {report['references']['reference_models']} (still DRAFT, not auto-approved)",
        f"- Runtime absolute Windows paths: {report['runtime_paths']['absolute_path_findings']}",
        f"- Forbidden outputs/caches/binaries: {report['artifact_exclusions']['forbidden_artifacts']}",
        f"- Issues: {report['issue_count']}",
        "",
        "This audit does not run ORCA and does not convert calculated data into experimental data.",
    ]
    if issues:
        md.extend(["", "## Issues", ""] + [f"- `{item.get('scope')}`: {item.get('problem')}" for item in issues])
    (META / "PORTABILITY_AUDIT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    inventory = write_inventory_and_hashes()
    print(json.dumps({"status": report["status"], "issues": len(issues), **inventory}, indent=2))
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
