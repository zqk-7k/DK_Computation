#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path
from typing import Any

from hpc_common import (
    atomic_json,
    enabled_row,
    load_json,
    read_manifest,
    resolve_manifest_path,
    row_value,
    safe_component,
    utc_now,
    write_csv,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize task state without querying or modifying the scheduler.")
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()

    root = args.package_root.resolve()
    manifest = resolve_manifest_path(root, args.manifest)
    _, rows = read_manifest(manifest)
    task_state_dir = root / "state/tasks" / safe_component(args.stage)
    lock_dir = root / "state/locks" / safe_component(args.stage)
    records: list[dict[str, Any]] = []
    for index, row in enumerate(rows, start=1):
        job_id = safe_component(row_value(row, "job_id", "model_id") or f"row_{index:05d}")
        state_path = task_state_dir / f"{job_id}.json"
        state = load_json(state_path) if state_path.is_file() else {}
        lock_path = lock_dir / f"{job_id}.lock"
        lock_age = round(time.time() - lock_path.stat().st_mtime, 1) if lock_path.exists() else ""
        records.append(
            {
                "manifest_row": index,
                "job_id": job_id,
                "enabled": enabled_row(row),
                "status": state.get("status", "NOT_STARTED"),
                "attempt": state.get("attempt", ""),
                "started_at": state.get("started_at", ""),
                "finished_at": state.get("finished_at", ""),
                "duration_seconds": state.get("duration_seconds", ""),
                "return_code": state.get("return_code", ""),
                "normal_termination_marker": state.get("normal_termination_marker", ""),
                "required_artifacts_ok": state.get("required_artifacts_ok", ""),
                "attempt_path": state.get("attempt_path", ""),
                "error": state.get("error", state.get("reason", "")),
                "lock_present": lock_path.exists(),
                "lock_age_seconds": lock_age,
                "scheduler": (state.get("scheduler") or {}).get("scheduler", "") if isinstance(state.get("scheduler"), dict) else "",
                "scheduler_job_id": (state.get("scheduler") or {}).get("job_id", "") if isinstance(state.get("scheduler"), dict) else "",
                "scheduler_array_index": (state.get("scheduler") or {}).get("array_index", "") if isinstance(state.get("scheduler"), dict) else "",
            }
        )

    status_counts = Counter(str(record["status"]) for record in records if record["enabled"])
    completed_tokens = {"COMPLETED", "RESUMED_ALREADY_COMPLETE"}
    enabled_records = [record for record in records if record["enabled"]]
    completed = sum(record["status"] in completed_tokens for record in enabled_records)
    all_complete = bool(enabled_records) and completed == len(enabled_records)
    payload = {
        "schema_version": "1.0",
        "notice": "CALCULATED MODEL WORKFLOW - NOT EXPERIMENTAL DATA",
        "timestamp": utc_now(),
        "stage": args.stage,
        "manifest": str(manifest.relative_to(root)),
        "enabled_tasks": len(enabled_records),
        "completed_tasks": completed,
        "all_complete": all_complete,
        "status_counts": dict(sorted(status_counts.items())),
        "records": records,
    }
    output_dir = root / "state/status"
    csv_path = output_dir / f"{safe_component(args.stage)}_status.csv"
    json_path = output_dir / f"{safe_component(args.stage)}_status.json"
    write_csv(csv_path, records[0].keys() if records else ("manifest_row", "job_id", "status"), records)
    atomic_json(json_path, payload)
    print(f"Stage {args.stage}: {completed}/{len(enabled_records)} complete; {dict(status_counts)}")
    print(f"CSV: {csv_path}")
    print(f"JSON: {json_path}")
    if args.require_complete and not all_complete:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
