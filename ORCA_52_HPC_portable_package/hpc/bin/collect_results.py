#!/usr/bin/env python3
from __future__ import annotations

import argparse
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


def main() -> int:
    parser = argparse.ArgumentParser(description="Create a hash-indexed manifest of completed ORCA attempts.")
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--site-config", type=Path)
    parser.add_argument("--require-all", action="store_true")
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
    for index, row in enumerate(rows, start=1):
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
        "records": collected,
    }
    atomic_json(json_path, payload)
    print(f"Collected {len(collected)}/{enabled_count} completed tasks")
    print(f"CSV: {csv_path}")
    print(f"JSON: {json_path}")
    if args.require_all and incomplete:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
