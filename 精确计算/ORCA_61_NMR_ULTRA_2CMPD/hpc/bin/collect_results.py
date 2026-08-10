#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any

from hpc_common import (
    atomic_json,
    enabled_row,
    load_json,
    read_manifest,
    resolve_manifest_path,
    resolve_package_path,
    row_value,
    safe_component,
    sha256_file,
    utc_now,
    write_csv,
)


COMPLETE = {"COMPLETED", "RESUMED_ALREADY_COMPLETE"}


def verified_ledger(attempt: Path, state: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Verify every attempt file against the runner-created SHA-256 ledger.

    Materialization is intentionally stricter than a normal collection pass:
    the top-level task state must bind the ledger hash, and every listed file
    must be a flat, non-empty, byte-identical attempt artifact.
    """

    ledger_path = attempt / "FILE_MANIFEST_SHA256.csv"
    if not ledger_path.is_file():
        raise ValueError(f"Missing attempt ledger: {ledger_path}")
    expected_ledger = str(state.get("file_manifest_sha256", "")).strip().upper()
    actual_ledger = sha256_file(ledger_path)
    if not expected_ledger or expected_ledger != actual_ledger:
        raise ValueError(
            f"Attempt ledger hash is absent or mismatched: expected={expected_ledger!r} actual={actual_ledger}"
        )
    with ledger_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"file", "size_bytes", "sha256"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Malformed attempt ledger header: {ledger_path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Attempt ledger is empty: {ledger_path}")
    verified: dict[str, dict[str, Any]] = {}
    for row in rows:
        name = str(row.get("file", "")).strip()
        candidate = Path(name)
        if not name or candidate.is_absolute() or candidate.name != name or name in {".", ".."}:
            raise ValueError(f"Unsafe attempt-ledger path: {name!r}")
        if name in verified:
            raise ValueError(f"Duplicate attempt-ledger entry: {name}")
        source = attempt / name
        if not source.is_file():
            raise ValueError(f"Attempt artifact is missing: {source}")
        try:
            expected_size = int(str(row.get("size_bytes", "")))
        except ValueError as exc:
            raise ValueError(f"Invalid size in attempt ledger for {name}") from exc
        expected_hash = str(row.get("sha256", "")).strip().upper()
        actual_hash = sha256_file(source)
        if source.stat().st_size != expected_size or not expected_hash or actual_hash != expected_hash:
            raise ValueError(f"Attempt artifact hash/size mismatch: {source}")
        verified[name] = {
            "path": source,
            "size_bytes": expected_size,
            "sha256": actual_hash,
        }
    return verified


def materialize_job_directory(
    root: Path,
    row: dict[str, str],
    verified: dict[str, dict[str, Any]],
    *,
    stage: str,
    job_id: str,
    state: dict[str, Any],
) -> list[dict[str, Any]]:
    raw_job = row_value(
        row,
        "portable_job_directory",
        "relative_job_directory",
        "working_directory",
        "job_directory",
    )
    target = resolve_package_path(root, raw_job, field="job_directory")
    selected = {name: record for name, record in verified.items() if name != "run_metadata.json"}
    attempt = Path(str(state.get("attempt_path", "")))
    if not attempt.is_absolute():
        attempt = root / attempt
    attempt = attempt.resolve()
    try:
        attempt_relative = attempt.relative_to(root.resolve()).as_posix()
    except ValueError as exc:
        raise ValueError(f"Dependency attempt escapes package root: {attempt}") from exc
    receipt_path = target / f".hpc_materialization__{safe_component(stage)}__{safe_component(job_id)}.json"
    receipt = {
        "schema_version": "1.0",
        "status": "HASH_VERIFIED_MATERIALIZATION",
        "stage": stage,
        "job_id": job_id,
        "attempt_path": attempt_relative,
        "file_manifest_sha256": str(state.get("file_manifest_sha256", "")).upper(),
        "artifacts": [
            {
                "file": name,
                "size_bytes": record["size_bytes"],
                "sha256": record["sha256"],
            }
            for name, record in sorted(selected.items())
        ],
    }

    # Check every destination before writing anything.  Existing identical
    # files are accepted; any differing file is a hard conflict.
    for name, record in selected.items():
        destination = target / name
        resolved = destination.resolve()
        try:
            resolved.relative_to(target.resolve())
        except ValueError as exc:
            raise ValueError(f"Materialization target escapes job directory: {destination}") from exc
        if destination.exists():
            if not destination.is_file() or sha256_file(destination) != record["sha256"]:
                raise ValueError(f"Refusing to overwrite conflicting job-local artifact: {destination}")
    if receipt_path.exists() and load_json(receipt_path) != receipt:
        raise ValueError(f"Refusing conflicting materialization receipt: {receipt_path}")

    materialized: list[dict[str, Any]] = []
    for name, record in selected.items():
        source = Path(record["path"])
        destination = target / name
        action = "already_present"
        if not destination.exists():
            fd, temporary_name = tempfile.mkstemp(
                dir=str(target), prefix=f".{name}.", suffix=".materialize.tmp"
            )
            os.close(fd)
            temporary = Path(temporary_name)
            try:
                shutil.copy2(source, temporary)
                if sha256_file(temporary) != record["sha256"]:
                    raise ValueError(f"Materialized temporary hash mismatch: {temporary}")
                os.replace(temporary, destination)
            finally:
                if temporary.exists():
                    temporary.unlink()
            action = "copied"
        materialized.append(
            {
                "file": name,
                "sha256": record["sha256"],
                "size_bytes": record["size_bytes"],
                "action": action,
            }
        )
    if not receipt_path.exists():
        atomic_json(receipt_path, receipt)
    return materialized


def main() -> int:
    parser = argparse.ArgumentParser(description="Create a hash-indexed manifest of completed ORCA attempts.")
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--site-config", type=Path)
    parser.add_argument("--require-all", action="store_true")
    parser.add_argument(
        "--materialize-job-directories",
        action="store_true",
        help="Hash-verify completed attempts and copy artifacts into their package job directories",
    )
    args = parser.parse_args()

    root = args.package_root.resolve()
    manifest = resolve_manifest_path(root, args.manifest)
    site_path = (args.site_config or root / "hpc/config/site_config.json").resolve()
    site = load_json(site_path)
    _, rows = read_manifest(manifest)
    state_dir = root / "state/tasks" / safe_component(args.stage)
    collected: list[dict[str, Any]] = []
    enabled_count = 0
    incomplete: list[str] = []
    materialization_errors: list[dict[str, str]] = []
    for index, row in enumerate(rows, start=1):
        declared_stage = row_value(row, "stage")
        if declared_stage and declared_stage != args.stage:
            continue
        if not enabled_row(row):
            continue
        enabled_count += 1
        job_id = safe_component(row_value(row, "job_id", "model_id") or f"row_{index:05d}")
        state_path = state_dir / f"{job_id}.json"
        if not state_path.is_file():
            incomplete.append(job_id)
            continue
        state = load_json(state_path)
        if state.get("status") not in COMPLETE:
            incomplete.append(job_id)
            continue
        attempt_path = Path(str(state.get("attempt_path", "")))
        if not attempt_path.is_absolute():
            attempt_path = root / attempt_path
        attempt_path = attempt_path.resolve()
        try:
            attempt_relative = attempt_path.relative_to(root)
        except ValueError:
            incomplete.append(job_id)
            continue
        output_file = attempt_path / str(state.get("output_file", ""))
        if not output_file.is_file():
            incomplete.append(job_id)
            continue
        file_manifest = attempt_path / "FILE_MANIFEST_SHA256.csv"
        materialized: list[dict[str, Any]] = []
        if args.materialize_job_directories:
            try:
                ledger = verified_ledger(attempt_path, state)
                materialized = materialize_job_directory(
                    root,
                    row,
                    ledger,
                    stage=args.stage,
                    job_id=job_id,
                    state=state,
                )
            except Exception as exc:
                incomplete.append(job_id)
                materialization_errors.append(
                    {"job_id": job_id, "error": f"{type(exc).__name__}: {exc}"}
                )
                continue
        collected.append(
            {
                "manifest_row": index,
                "job_id": job_id,
                "stage": args.stage,
                "input_sha256": state.get("input_sha256", ""),
                "attempt": state.get("attempt", ""),
                "attempt_directory": str(attempt_relative).replace("\\", "/"),
                "output_file": output_file.name,
                "output_sha256": sha256_file(output_file),
                "normal_termination_marker": state.get("normal_termination_marker", True),
                "required_artifacts_ok": state.get("required_artifacts_ok", True),
                "file_manifest": str(file_manifest.relative_to(root)).replace("\\", "/") if file_manifest.is_file() else "",
                "file_manifest_sha256": sha256_file(file_manifest) if file_manifest.is_file() else "",
                "materialized_file_count": len(materialized),
                "materialization_status": (
                    "HASH_VERIFIED_MATERIALIZED" if args.materialize_job_directories else "NOT_REQUESTED"
                ),
                "finished_at": state.get("finished_at", ""),
                "duration_seconds": state.get("duration_seconds", ""),
                "notice": "CALCULATED MODEL - NOT EXPERIMENTAL DATA",
            }
        )

    configured_results_root = resolve_package_path(
        root,
        str(site.get("results_root") or "results/hpc"),
        field="results_root",
        allow_missing=True,
    )
    results_root = configured_results_root / safe_component(args.stage)
    results_root.mkdir(parents=True, exist_ok=True)
    csv_path = results_root / f"{safe_component(args.stage)}_completed_results_manifest.csv"
    json_path = results_root / f"{safe_component(args.stage)}_completed_results_manifest.json"
    fields = collected[0].keys() if collected else (
        "manifest_row",
        "job_id",
        "stage",
        "input_sha256",
        "attempt",
        "attempt_directory",
        "output_file",
        "output_sha256",
        "notice",
    )
    write_csv(csv_path, fields, collected)
    payload = {
        "schema_version": "1.0",
        "notice": "CALCULATED MODEL WORKFLOW - NOT EXPERIMENTAL DATA",
        "timestamp": utc_now(),
        "stage": args.stage,
        "manifest": str(manifest.relative_to(root)),
        "enabled_tasks": enabled_count,
        "collected_tasks": len(collected),
        "all_complete": len(collected) == enabled_count and enabled_count > 0,
        "incomplete_job_ids": incomplete,
        "materialize_job_directories": args.materialize_job_directories,
        "materialization_errors": materialization_errors,
        "records": collected,
    }
    atomic_json(json_path, payload)
    print(f"Collected {len(collected)}/{enabled_count} completed tasks")
    print(f"CSV: {csv_path}")
    print(f"JSON: {json_path}")
    if materialization_errors:
        return 3
    if args.require_all and incomplete:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
