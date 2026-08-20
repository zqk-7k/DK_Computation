#!/usr/bin/env python3
"""Create the effective held GOAT manifest and rendered Slurm launchers."""

from __future__ import annotations

import argparse
import csv
import json
from io import StringIO
from pathlib import Path

from package_io import atomic_write_text
from render_slurm import render


PROFILE_BY_TIER = {"R": "goat_small", "M": "goat_medium", "H": "goat_large", "X": "goat_large"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--source-manifest", type=Path)
    parser.add_argument("--site-config", type=Path, required=True)
    parser.add_argument("--output-manifest", type=Path)
    args = parser.parse_args()
    root = args.package_root.resolve()
    source = (args.source_manifest or root / "data" / "GOAT_explicit_job_manifest.csv").resolve()
    output = (args.output_manifest or root / "hpc" / "manifests" / "goat_adaptive_manifest.csv").resolve()
    site = json.loads(args.site_config.read_text(encoding="utf-8"))
    resources_path = root / "hpc" / "resource_profiles.json"
    resources = json.loads(resources_path.read_text(encoding="utf-8"))["profiles"]
    template_path = root / "hpc" / "slurm" / "run_manifest_task.sbatch.in"
    template = template_path.read_text(encoding="utf-8")
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle); rows = list(reader)
    output_rows = []
    scripts_dir = root / "hpc" / "job_scripts" / "goat"
    scripts_dir.mkdir(parents=True, exist_ok=True)
    for row in rows:
        job_id = row["job_id"]
        tier = row.get("complexity_tier", "R").upper()
        profile_name = PROFILE_BY_TIER.get(tier)
        if not profile_name: raise ValueError(f"unknown complexity tier: {tier}")
        profile = resources[profile_name]
        rel_job_dir = row.get("relative_job_directory") or row.get("portable_job_directory") or row["job_directory"]
        job_dir = (root / rel_job_dir).resolve()
        input_file = job_dir / row.get("input_file", "goat.inp")
        script = scripts_dir / f"{job_id}.sbatch"
        values = {
            "JOB_NAME": job_id, "SCHEDULER_CPUS": profile["scheduler_cpus"],
            "REQUEST_MEMORY_MB": profile["requested_memory_mb"], "WALLTIME": profile["walltime"],
            "JOB_DIR": job_dir.as_posix(), "PYTHON_EXECUTABLE": site["python_executable"],
            "PACKAGE_ROOT": root.as_posix(), "SITE_CONFIG": args.site_config.resolve().as_posix(),
            "RESOURCE_PROFILES": resources_path.resolve().as_posix(), "PROFILE": profile_name,
            "INPUT_FILE": input_file.as_posix(),
        }
        atomic_write_text(script, render(template, values))
        output_rows.append({
            "job_id": job_id, "compound_id": row["figure_id"], "scenario_id": row["scenario_id"],
            "wave": row["wave"], "target_success": row["target_success"], "max_candidates": row["max_candidates"],
            "state": "HELD", "accepted": "false", "job_dir": rel_job_dir.replace("\\", "/"),
            "submit_script": script.relative_to(root).as_posix(), "resource_profile": profile_name,
            "scheduler_job_id": "", "source_approval_status": row.get("approval_status", ""),
            "source_input_bundle_sha256": row.get("input_bundle_sha256", ""),
        })
    fields = list(output_rows[0]) if output_rows else []
    stream = StringIO(newline=""); writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader(); writer.writerows(output_rows)
    atomic_write_text(output, stream.getvalue())
    print(json.dumps({"jobs": len(output_rows), "manifest": str(output), "scripts": str(scripts_dir)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
