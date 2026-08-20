#!/usr/bin/env python3
"""Render a fail-closed Slurm launcher from a resource profile."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from package_io import atomic_write_text


TOKEN_RE = re.compile(r"\{\{([A-Z0-9_]+)\}\}")


def render(template: str, values: dict[str, object]) -> str:
    unknown = sorted(set(TOKEN_RE.findall(template)) - set(values))
    if unknown:
        raise ValueError(f"missing template values: {', '.join(unknown)}")
    result = template
    for key, value in values.items():
        result = result.replace("{{" + key + "}}", str(value))
    remaining = sorted(set(TOKEN_RE.findall(result)))
    if remaining:
        raise ValueError(f"unresolved template values: {', '.join(remaining)}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--site-config", type=Path, required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--job-dir", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--job-name")
    args = parser.parse_args()

    root = args.package_root.resolve()
    site = json.loads(args.site_config.read_text(encoding="utf-8"))
    rp_path = root / "hpc" / "resource_profiles.json"
    resource = json.loads(rp_path.read_text(encoding="utf-8"))
    if args.profile not in resource["profiles"]:
        raise SystemExit(f"unknown resource profile: {args.profile}")
    p = resource["profiles"][args.profile]
    job_dir = args.job_dir.resolve()
    input_file = args.input.resolve()
    if input_file.parent != job_dir:
        raise SystemExit("input must be directly inside job-dir")
    values = {
        "JOB_NAME": args.job_name or job_dir.name,
        "SCHEDULER_CPUS": p["scheduler_cpus"],
        "REQUEST_MEMORY_MB": p["requested_memory_mb"],
        "WALLTIME": p["walltime"],
        "JOB_DIR": job_dir.as_posix(),
        "PYTHON_EXECUTABLE": site["python_executable"],
        "PACKAGE_ROOT": root.as_posix(),
        "SITE_CONFIG": args.site_config.resolve().as_posix(),
        "RESOURCE_PROFILES": rp_path.resolve().as_posix(),
        "PROFILE": args.profile,
        "INPUT_FILE": input_file.as_posix(),
    }
    template = (root / "hpc" / "slurm" / "run_manifest_task.sbatch.in").read_text(encoding="utf-8")
    atomic_write_text(args.output.resolve(), render(template, values))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
