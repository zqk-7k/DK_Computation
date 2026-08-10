#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import time
from pathlib import Path

from hpc_common import atomic_json, load_json, safe_component, utc_now


def main() -> int:
    parser = argparse.ArgumentParser(description="Remove one stale task lock only after an explicit scheduler check.")
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--confirm-scheduler-job-is-not-running", action="store_true", required=True)
    parser.add_argument("--minimum-age-seconds", type=int, default=3600)
    args = parser.parse_args()
    root = args.package_root.resolve()
    lock = root / "state/locks" / safe_component(args.stage) / f"{safe_component(args.job_id)}.lock"
    if not lock.is_dir():
        print(f"No lock exists: {lock}")
        return 0
    age = time.time() - lock.stat().st_mtime
    if age < args.minimum_age_seconds:
        raise SystemExit(f"Lock age {age:.0f}s is below minimum {args.minimum_age_seconds}s")
    owner_path = lock / "owner.json"
    owner = load_json(owner_path) if owner_path.is_file() else {}
    audit = {
        "schema_version": "1.0",
        "timestamp": utc_now(),
        "action": "STALE_LOCK_REMOVED_AFTER_EXPLICIT_SCHEDULER_CONFIRMATION",
        "stage": args.stage,
        "job_id": args.job_id,
        "lock_age_seconds": round(age, 1),
        "prior_owner": owner,
    }
    audit_path = root / "state/unlock_audit" / f"{safe_component(args.stage)}_{safe_component(args.job_id)}_{int(time.time())}.json"
    atomic_json(audit_path, audit)
    shutil.rmtree(lock)
    print(f"Removed stale lock; audit={audit_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
