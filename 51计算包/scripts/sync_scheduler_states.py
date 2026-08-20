#!/usr/bin/env python3
"""Reconcile active manifest jobs with scheduler accounting (dry-run by default)."""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from package_io import ACTIVE_STATES, atomic_write_json, load_json


STATE_MAP = {
    "PENDING": "PENDING", "CONFIGURING": "PENDING", "RUNNING": "RUNNING", "COMPLETING": "COMPLETING",
    "COMPLETED": "COMPLETED", "FAILED": "FAILED", "CANCELLED": "CANCELLED", "TIMEOUT": "TIMEOUT",
    "OUT_OF_MEMORY": "FAILED", "NODE_FAIL": "FAILED", "BOOT_FAIL": "FAILED", "PREEMPTED": "FAILED",
}


def scheduler_state(text: str) -> str:
    tokens = re.findall(r"[A-Z_]+", text.upper())
    for token in tokens:
        if token in STATE_MAP: return STATE_MAP[token]
    return "UNKNOWN"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--site-config", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args(); root = args.package_root.resolve(); site = load_json(args.site_config.resolve())
    with args.manifest.open("r", encoding="utf-8-sig", newline="") as handle: rows = list(csv.DictReader(handle))
    query_template = [str(x) for x in site.get("query_command", [])]
    if not query_template: raise SystemExit("site config query_command is empty")
    results = []
    for row in rows:
        job_id = row.get("job_id", ""); scheduler_id = row.get("scheduler_job_id", "").strip()
        raw_dir = row.get("job_dir") or row.get("job_directory") or row.get("relative_job_directory")
        if not (job_id and scheduler_id and raw_dir): continue
        job_dir = Path(raw_dir); job_dir = job_dir if job_dir.is_absolute() else (root / job_dir).resolve()
        state_file = job_dir / "job_state.json"; current = str(row.get("state", "HELD")).upper()
        if state_file.is_file(): current = str(load_json(state_file).get("state", current)).upper()
        if current not in ACTIVE_STATES: continue
        command = [part.replace("{scheduler_job_id}", scheduler_id) for part in query_template]
        proc = subprocess.run(command, text=True, capture_output=True, check=False)
        observed = scheduler_state(proc.stdout) if proc.returncode == 0 else "QUERY_FAILED"
        record = {"job_id": job_id, "scheduler_job_id": scheduler_id, "prior_state": current, "observed_state": observed,
                  "query_returncode": proc.returncode, "stderr": proc.stderr.strip()}
        if args.apply and observed in STATE_MAP.values() and observed != current:
            if observed == "COMPLETED" and not (job_dir / "attempts").is_dir():
                # A scheduler COMPLETED state without a runner receipt is not a scientific success.
                observed = "FAILED"
                record["note"] = "scheduler completed but runner attempt/receipt is absent"
            atomic_write_json(state_file, {"schema_version": "1.0", "job_id": job_id, "state": observed,
                                          "scheduler_job_id": scheduler_id, "updated_utc": datetime.now(timezone.utc).isoformat(),
                                          "accepted": False, "source": "scheduler_reconciliation"})
            record["applied"] = True
        else: record["applied"] = False
        results.append(record)
    payload = {"schema_version": "1.0", "mode": "APPLY" if args.apply else "DRY_RUN", "results": results}
    print(json.dumps(payload, ensure_ascii=False, indent=2)); return 0


if __name__ == "__main__":
    raise SystemExit(main())
