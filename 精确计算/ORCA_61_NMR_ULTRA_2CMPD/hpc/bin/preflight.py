#!/usr/bin/env python3
from __future__ import annotations

import argparse
import fnmatch
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path
import hashlib
from typing import Any

from hpc_common import (
    approval_problem,
    atomic_json,
    dependency_contract,
    enabled_row,
    input_coordinate_references,
    load_json,
    parse_bool,
    parse_orca_resources,
    read_manifest,
    required_artifact_patterns,
    resolve_manifest_path,
    resolve_package_path,
    resolve_required_helpers,
    row_value,
    safe_component,
    sha256_file,
    utc_now,
    verify_materialized_dependency_receipt,
)


def add_check(checks: list[dict[str, Any]], level: str, name: str, message: str, **details: Any) -> None:
    checks.append({"level": level, "name": name, "message": message, **details})


def resolve_input(package_root: Path, row: dict[str, str], input_default: str) -> tuple[Path, Path]:
    raw_job = row_value(
        row,
        "portable_job_directory",
        "relative_job_directory",
        "working_directory",
        "job_directory",
    )
    job_dir = resolve_package_path(package_root, raw_job, field="job_directory")
    raw_input = row_value(row, "portable_input_file", "input_file") or input_default
    token = raw_input.replace("\\", "/")
    if len(token) >= 3 and token[1:3] == ":/":
        token = Path(token).name
    candidate = job_dir / token
    if not candidate.is_file():
        package_candidate = package_root / token
        candidate = package_candidate if package_candidate.is_file() else candidate
    if not candidate.is_file():
        raise FileNotFoundError(candidate)
    return job_dir, candidate.resolve()


def main() -> int:
    parser = argparse.ArgumentParser(description="Static and executable preflight for an ORCA HPC stage.")
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--scheduler", choices=("slurm", "pbs"), required=True)
    parser.add_argument("--site-config", type=Path)
    parser.add_argument("--stage-config", type=Path)
    parser.add_argument("--allow-non-linux", action="store_true", help="Testing only")
    parser.add_argument("--no-execute", action="store_true", help="Do not execute ORCA --version")
    args = parser.parse_args()

    root = args.package_root.resolve()
    site_path = (args.site_config or root / "hpc/config/site_config.json").resolve()
    resources_path = (args.stage_config or root / "hpc/config/stage_resources.json").resolve()
    site = load_json(site_path)
    strict_hash = parse_bool(site.get("strict_hash"), default=True)
    stages = load_json(resources_path).get("stages", {})
    if args.stage not in stages:
        raise SystemExit(f"Unknown stage: {args.stage}")
    stage = stages[args.stage]
    manifest = resolve_manifest_path(root, args.manifest)
    _, rows = read_manifest(manifest)
    checks: list[dict[str, Any]] = []

    if platform.system() == "Linux" or args.allow_non_linux:
        add_check(checks, "PASS", "operating_system", platform.platform())
    else:
        add_check(checks, "ERROR", "operating_system", "Production package requires Linux", actual=platform.platform())

    scheduler_command = "sbatch" if args.scheduler == "slurm" else "qsub"
    if shutil.which(scheduler_command):
        add_check(checks, "PASS", "scheduler_command", scheduler_command, path=shutil.which(scheduler_command))
    else:
        add_check(checks, "ERROR", "scheduler_command", f"{scheduler_command} is not on PATH")

    python_command = str(site.get("python_executable") or "python3")
    if Path(python_command).is_absolute() and Path(python_command).is_file():
        add_check(checks, "PASS", "python", python_command)
    elif shutil.which(python_command):
        add_check(checks, "PASS", "python", python_command, path=shutil.which(python_command))
    else:
        add_check(checks, "ERROR", "python", f"Python executable is unavailable: {python_command}")

    orca = Path(str(site.get("orca_executable") or "")).expanduser()
    if not orca.is_absolute() or not orca.is_file():
        add_check(checks, "ERROR", "orca_executable", "Set an existing absolute Linux path", configured=str(orca))
    elif not os.access(orca, os.X_OK):
        add_check(checks, "ERROR", "orca_executable", "ORCA file is not executable", configured=str(orca))
    else:
        actual_hash = sha256_file(orca)
        expected_hash = str(site.get("orca_sha256") or "").strip().upper()
        if strict_hash and not re.fullmatch(r"[0-9A-F]{64}", expected_hash):
            add_check(
                checks,
                "ERROR",
                "orca_sha256",
                "strict_hash=true requires the exact 64-hex SHA-256 of the Linux ORCA binary",
                actual=actual_hash,
            )
        elif expected_hash and expected_hash != actual_hash:
            add_check(checks, "ERROR", "orca_sha256", "ORCA binary hash mismatch", expected=expected_hash, actual=actual_hash)
        else:
            add_check(
                checks,
                "PASS" if expected_hash else "WARNING",
                "orca_sha256",
                "ORCA executable is hash-bound" if expected_hash else "ORCA hash is not pinned because strict_hash=false",
                path=str(orca),
                sha256=actual_hash,
            )
        if not args.no_execute:
            try:
                result = subprocess.run(
                    [str(orca), "--version"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    timeout=30,
                    check=False,
                )
                text = (result.stdout or "").strip()[-4000:]
                expected_version = str(site.get("orca_expected_version") or "").strip()
                if expected_version and expected_version not in text:
                    add_check(
                        checks,
                        "ERROR",
                        "orca_version",
                        "Configured ORCA version was not confirmed from `orca --version`",
                        expected=expected_version,
                        return_code=result.returncode,
                        output=text,
                    )
                else:
                    add_check(checks, "PASS", "orca_version", expected_version or "not pinned", output=text)
            except Exception as exc:
                add_check(checks, "ERROR", "orca_version", f"Version probe failed: {type(exc).__name__}: {exc}")

    if args.stage == "goat":
        xtb_raw = str(site.get("xtb_executable") or "").strip()
        xtb_candidates = [Path(xtb_raw).expanduser()] if xtb_raw else []
        if orca.is_absolute():
            xtb_candidates.extend(
                (
                    orca.parent / "otool_xtb",
                    orca.parent / "xtb",
                    orca.parent / "xtb-6.7.1pre" / "xtb",
                )
            )
        xtb = next((path for path in xtb_candidates if path.is_file() and os.access(path, os.X_OK)), None)
        if xtb:
            add_check(checks, "PASS", "xtb_executable", str(xtb))
        else:
            add_check(checks, "ERROR", "xtb_executable", "GOAT requires the ORCA-compatible xTB executable")

    enabled_count = 0
    blocked_count = 0
    path_errors = 0
    resource_errors = 0
    hash_absent = 0
    hash_errors = 0
    declared_stage_errors = 0
    seen_job_ids: set[str] = set()
    duplicate_job_ids = 0
    dependency_errors = 0
    helper_errors = 0
    required_contract_errors = 0
    for index, row in enumerate(rows, start=1):
        if not enabled_row(row):
            continue
        enabled_count += 1
        job_id = row_value(row, "job_id", "model_id") or f"row_{index:05d}"
        if job_id in seen_job_ids:
            duplicate_job_ids += 1
        seen_job_ids.add(job_id)
        declared_stage = row_value(row, "stage")
        if declared_stage and declared_stage != args.stage:
            declared_stage_errors += 1
        if approval_problem(row, args.stage):
            blocked_count += 1
        try:
            job_dir, input_path = resolve_input(root, row, str(stage["input_default"]))
            try:
                patterns = required_artifact_patterns(stage, row)
                if not patterns:
                    raise ValueError("No completion artifacts are declared")
                if any(fnmatch.fnmatch(input_path.name, pattern) for pattern in patterns):
                    raise ValueError("An input file can satisfy the completion artifact contract")
            except Exception:
                required_contract_errors += 1
            try:
                resolve_required_helpers(site, orca, row)
            except Exception:
                helper_errors += 1
            try:
                dependency, parent_stage, _ = dependency_contract(row, args.stage)
                if dependency and parent_stage == args.stage:
                    raise ValueError(
                        "HPC arrays do not guarantee same-stage dependency order; "
                        "split and materialize the parent stage first"
                    )
                if dependency:
                    verify_materialized_dependency_receipt(
                        root, row, job_dir, args.stage
                    )
            except Exception:
                dependency_errors += 1
            for coordinate in input_coordinate_references(input_path):
                if not coordinate.is_file():
                    raise FileNotFoundError(coordinate)
            request = parse_orca_resources(input_path)
            if request["required_cpu_slots"] > int(stage["cpus"]):
                resource_errors += 1
            maxcore_total = request["nprocs"] * request["maxcore_mb_per_rank"]
            if maxcore_total and maxcore_total > int(stage["memory_mb"]) * 0.90:
                resource_errors += 1
            expected_hash = row_value(
                row,
                "expected_input_sha256",
                "input_file_sha256",
                "input_sha256",
                "goat_input_sha256",
                "opt_input_sha256",
                "shielding_input_sha256",
                "j_input_sha256",
            )
            if strict_hash and not expected_hash:
                hash_absent += 1
            elif expected_hash and sha256_file(input_path) != expected_hash.upper():
                hash_errors += 1

            expected_coordinate = row_value(
                row,
                "expected_coordinate_sha256",
                "coordinate_sha256",
                "start_xyz_sha256",
                "selected_xyz_sha256",
                "optimized_xyz_sha256",
            )
            if expected_coordinate:
                coordinate_raw = row_value(row, "coordinate_file")
                if coordinate_raw:
                    coordinate = job_dir / coordinate_raw.replace("\\", "/")
                    candidates = [coordinate.resolve()]
                else:
                    candidates = input_coordinate_references(input_path)
                if len(candidates) != 1 or not candidates[0].is_file():
                    hash_errors += 1
                elif sha256_file(candidates[0]) != expected_coordinate.upper():
                    hash_errors += 1

            expected_bundle = row_value(
                row, "expected_input_bundle_sha256", "input_bundle_sha256"
            )
            if expected_bundle:
                bundle_names = (
                    "goat.inp",
                    "start.xyz",
                    "atom_map.csv",
                    "isotope_map.csv",
                    "structure_validation.json",
                )
                if all((job_dir / name).is_file() for name in bundle_names):
                    digest = hashlib.sha256()
                    for name in bundle_names:
                        digest.update(name.encode())
                        digest.update((job_dir / name).read_bytes())
                    if digest.hexdigest().upper() != expected_bundle.upper():
                        hash_errors += 1
                else:
                    recorded = job_dir / "INPUT_BUNDLE_SHA256.txt"
                    if not recorded.is_file() or recorded.read_text(encoding="utf-8").strip().upper() != expected_bundle.upper():
                        hash_errors += 1
        except Exception:
            path_errors += 1

    add_check(checks, "PASS" if enabled_count else "ERROR", "manifest_enabled_rows", str(enabled_count), total_rows=len(rows))
    add_check(
        checks,
        "ERROR" if blocked_count else "PASS",
        "scientific_gates",
        f"{blocked_count} enabled rows are not approved" if blocked_count else "All enabled rows pass manifest gate tokens",
    )
    add_check(
        checks,
        "ERROR" if path_errors else "PASS",
        "portable_input_paths",
        f"{path_errors} enabled rows have missing/non-portable inputs" if path_errors else "All enabled input paths are package-relative and present",
    )
    add_check(
        checks,
        "ERROR" if resource_errors else "PASS",
        "resource_fit",
        f"{resource_errors} rows exceed CPU or 90% memory policy" if resource_errors else "All ORCA %pal/%maxcore requests fit the stage allocation",
    )
    add_check(
        checks,
        "ERROR" if (hash_absent or hash_errors) else "PASS",
        "manifest_input_hashes",
        (
            f"{hash_absent} enabled rows lack an expected input hash; "
            f"{hash_errors} input/coordinate/bundle hashes mismatch"
            if (hash_absent or hash_errors)
            else "Expected input, coordinate and bundle hashes are present and match"
        ),
    )
    add_check(
        checks,
        "ERROR" if (declared_stage_errors or duplicate_job_ids) else "PASS",
        "manifest_identity",
        (
            f"{declared_stage_errors} rows declare a different stage; "
            f"{duplicate_job_ids} duplicate job IDs"
            if (declared_stage_errors or duplicate_job_ids)
            else "Stage declarations and job IDs are consistent"
        ),
    )
    add_check(
        checks,
        "ERROR" if dependency_errors else "PASS",
        "stage_dependencies",
        (
            f"{dependency_errors} enabled rows lack verified materialized dependencies"
            if dependency_errors
            else "All stage dependencies are hash-bound and materialized"
        ),
    )
    add_check(
        checks,
        "ERROR" if helper_errors else "PASS",
        "required_helpers",
        (
            f"{helper_errors} enabled rows require unavailable helper executables"
            if helper_errors
            else "All row-level helper executable requirements are satisfied"
        ),
    )
    add_check(
        checks,
        "ERROR" if required_contract_errors else "PASS",
        "required_output_contracts",
        (
            f"{required_contract_errors} enabled rows have unsafe/empty output contracts"
            if required_contract_errors
            else "Stage and row-level required outputs are safe and cannot match their inputs"
        ),
    )

    scratch_raw = str(site.get("scratch_base") or "").strip()
    if scratch_raw:
        scratch = Path(os.path.expandvars(os.path.expanduser(scratch_raw)))
        if scratch.is_dir() and os.access(scratch, os.W_OK):
            add_check(checks, "PASS", "scratch", str(scratch))
        else:
            add_check(checks, "ERROR", "scratch", "Configured scratch is missing or not writable", path=str(scratch))
    else:
        add_check(checks, "PASS", "scratch", "Scheduler node-local scratch will be auto-selected from SLURM_TMPDIR/TMPDIR")

    errors = sum(check["level"] == "ERROR" for check in checks)
    warnings = sum(check["level"] == "WARNING" for check in checks)
    payload = {
        "schema_version": "1.0",
        "notice": "CALCULATED MODEL WORKFLOW - NOT EXPERIMENTAL DATA",
        "timestamp": utc_now(),
        "status": "PASS" if errors == 0 else "FAIL",
        "errors": errors,
        "warnings": warnings,
        "package_root": str(root),
        "manifest": str(manifest.relative_to(root)),
        "manifest_sha256": sha256_file(manifest),
        "stage": args.stage,
        "scheduler": args.scheduler,
        "checks": checks,
    }
    timestamp = utc_now().replace(":", "").replace("+", "_")
    output = root / "state/preflight" / f"{args.stage}_{args.scheduler}_{timestamp}.json"
    atomic_json(output, payload)
    print(f"Preflight {payload['status']}: errors={errors}, warnings={warnings}; report={output}")
    for check in checks:
        print(f"[{check['level']}] {check['name']}: {check['message']}")
    return 0 if errors == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
