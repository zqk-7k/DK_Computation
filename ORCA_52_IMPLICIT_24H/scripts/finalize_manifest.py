#!/usr/bin/env python3
"""Generate the immutable final result-file manifest as the last write step."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import uuid
from pathlib import Path

from common import read_json, sha256, utc_stamp


MANIFEST_NAME = "FINAL_FILE_MANIFEST_SHA256.csv"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    if not run_dir.is_dir():
        raise FileNotFoundError(run_dir)

    # Prepare a future pipeline_final and a manifest that hashes those exact
    # future bytes.  Commit the manifest first and pipeline_final second.  If a
    # walltime kill occurs between the two atomic replaces, pipeline_final stays
    # PENDING and the inconsistency is detectable; it can never falsely claim
    # COMPLETE while the manifest is missing or stale.
    pipeline_final_path = run_dir / "pipeline_final.json"
    if not pipeline_final_path.is_file():
        raise FileNotFoundError(pipeline_final_path)
    pipeline_final = read_json(pipeline_final_path)
    pipeline_final["manifest_generation_status"] = "COMPLETE"
    pipeline_final["manifest_generated_utc"] = utc_stamp()
    pipeline_final["manifest_exclusion"] = MANIFEST_NAME
    future_pipeline_bytes = (
        json.dumps(pipeline_final, ensure_ascii=False, indent=2) + "\n"
    ).encode("utf-8")

    for stale in run_dir.glob(".finalize_*.tmp"):
        stale.unlink(missing_ok=True)

    rows = []
    for path in sorted(run_dir.rglob("*")):
        if (
            not path.is_file()
            or path.name == MANIFEST_NAME
            or path == pipeline_final_path
            or path.name.startswith(".finalize_")
        ):
            continue
        rows.append(
            {
                "path": path.relative_to(run_dir).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    rows.append(
        {
            "path": "pipeline_final.json",
            "size_bytes": len(future_pipeline_bytes),
            "sha256": hashlib.sha256(future_pipeline_bytes).hexdigest().upper(),
        }
    )
    rows.sort(key=lambda row: str(row["path"]))

    token = uuid.uuid4().hex
    manifest_tmp = run_dir / f".finalize_manifest_{token}.tmp"
    pipeline_tmp = run_dir / f".finalize_pipeline_{token}.tmp"
    with manifest_tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "size_bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)
    pipeline_tmp.write_bytes(future_pipeline_bytes)
    os.replace(manifest_tmp, run_dir / MANIFEST_NAME)
    os.replace(pipeline_tmp, pipeline_final_path)
    print(f"Manifested {len(rows)} files under {run_dir}; {MANIFEST_NAME} excludes itself by design.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
