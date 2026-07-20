#!/usr/bin/env python3
"""Read-only comparison of portable GOAT physics inputs with the local source."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT.parent / "ORCA_52_NMR_serial_high_quality"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def normalized_validation(value: Any) -> Any:
    if isinstance(value, dict):
        result = {key: normalized_validation(item) for key, item in value.items()}
        policy = result.get("isotope_mass_policy")
        if isinstance(policy, dict) and "policy_file" in policy:
            policy["policy_file"] = "data/ORCA_6_1_1_isotope_mass_policy.json"
        return result
    if isinstance(value, list):
        return [normalized_validation(item) for item in value]
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    args = parser.parse_args()
    source = args.source.resolve()
    if not (source / "data" / "GOAT_explicit_job_manifest.csv").is_file():
        raise FileNotFoundError(source)
    with (ROOT / "data" / "GOAT_explicit_job_manifest.csv").open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle))
    files_compared = 0
    validation_compared = 0
    mismatches: list[dict[str, str]] = []
    for row in rows:
        rel = Path(row["portable_job_directory"])
        source_dir = source / rel
        package_dir = ROOT / rel
        for name in ("goat.inp", "start.xyz", "atom_map.csv", "isotope_map.csv"):
            files_compared += 1
            source_hash = sha256(source_dir / name)
            package_hash = sha256(package_dir / name)
            if source_hash != package_hash:
                mismatches.append(
                    {"job_id": row["job_id"], "file": name, "source_sha256": source_hash, "package_sha256": package_hash}
                )
        source_validation = normalized_validation(
            json.loads((source_dir / "structure_validation.json").read_text(encoding="utf-8-sig"))
        )
        package_validation = normalized_validation(
            json.loads((package_dir / "structure_validation.json").read_text(encoding="utf-8-sig"))
        )
        validation_compared += 1
        if source_validation != package_validation:
            mismatches.append({"job_id": row["job_id"], "file": "structure_validation.json", "problem": "semantic mismatch"})
    report = {
        "status": "PASS" if not mismatches else "FAIL",
        "mode": "READ_ONLY_SOURCE_COMPARISON",
        "source_location_recorded": "sibling active project; absolute workstation path deliberately omitted",
        "source_files_written": 0,
        "jobs": len(rows),
        "physics_files_compared": files_compared,
        "structure_validation_records_compared": validation_compared,
        "expected_portable_difference": "policy_file is package-relative and bundle hashes were regenerated accordingly",
        "mismatch_count": len(mismatches),
        "mismatches": mismatches,
    }
    output = ROOT / "package_metadata" / "SOURCE_COPY_READ_ONLY_AUDIT.json"
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if not mismatches else 1


if __name__ == "__main__":
    raise SystemExit(main())
