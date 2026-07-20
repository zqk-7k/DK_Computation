#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from hpc_common import atomic_json, load_json, read_manifest, resolve_manifest_path, safe_component, sha256_file, utc_now


TOKENS = (
    "JOB_NAME",
    "CPUS",
    "MEMORY_MB",
    "WALLTIME",
    "ARRAY_SPEC",
    "MAX_CONCURRENT",
    "SCHEDULER_LOG_DIR",
    "OPTIONAL_DIRECTIVES",
    "MODULE_COMMANDS",
    "ORCA_EXE",
    "SITE_CONFIG",
    "PACKAGE_ROOT",
    "MANIFEST",
    "PYTHON_EXE",
    "STAGE",
)


def shell_quote(value: str | Path) -> str:
    return shlex.quote(str(value))


def render(template: str, values: dict[str, str]) -> str:
    result = template
    for token in TOKENS:
        result = result.replace("{{" + token + "}}", values[token])
    leftovers = sorted(set(re.findall(r"\{\{[A-Z0-9_]+\}\}", result)))
    if leftovers:
        raise ValueError(f"Unresolved template tokens: {leftovers}")
    return result


def optional_directives(scheduler: str, site: dict[str, Any], dependency: str) -> list[str]:
    values: list[str] = []
    if scheduler == "slurm":
        config = site.get("slurm", {})
        if config.get("account"):
            values.append(f"#SBATCH --account={config['account']}")
        if config.get("partition"):
            values.append(f"#SBATCH --partition={config['partition']}")
        if config.get("qos"):
            values.append(f"#SBATCH --qos={config['qos']}")
        if dependency:
            values.append(f"#SBATCH --dependency=afterok:{dependency}")
        values.extend(str(item) for item in config.get("extra_directives", []))
    else:
        config = site.get("pbs", {})
        if config.get("project"):
            values.append(f"#PBS -P {config['project']}")
        if config.get("queue"):
            values.append(f"#PBS -q {config['queue']}")
        if dependency:
            values.append(f"#PBS -W depend=afterok:{dependency}")
        values.extend(str(item) for item in config.get("extra_directives", []))
    return values


def run_preflight(
    root: Path,
    manifest: Path,
    stage: str,
    scheduler: str,
    site_path: Path,
    resources_path: Path,
) -> None:
    command = [
        sys.executable,
        str(root / "hpc/bin/preflight.py"),
        "--package-root",
        str(root),
        "--manifest",
        str(manifest),
        "--stage",
        stage,
        "--scheduler",
        scheduler,
        "--site-config",
        str(site_path),
        "--stage-config",
        str(resources_path),
    ]
    result = subprocess.run(command, check=False)
    if result.returncode != 0:
        raise SystemExit("Submission stopped because preflight failed.")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Render or submit a bounded-concurrency SLURM/PBS ORCA array. Default: render only."
    )
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--scheduler", choices=("slurm", "pbs"), required=True)
    parser.add_argument("--site-config", type=Path)
    parser.add_argument("--stage-config", type=Path)
    parser.add_argument("--max-concurrent", type=int)
    parser.add_argument(
        "--array-range",
        help="Optional one-based inclusive pilot range, for example 1-2 or 1-10",
    )
    parser.add_argument("--job-name")
    parser.add_argument("--dependency", default="", help="Scheduler job ID; uses afterok")
    parser.add_argument("--submit", action="store_true", help="Actually call sbatch/qsub")
    parser.add_argument("--skip-preflight", action="store_true", help="Not recommended; recorded")
    args = parser.parse_args()

    root = args.package_root.resolve()
    manifest = resolve_manifest_path(root, args.manifest)
    site_path = (args.site_config or root / "hpc/config/site_config.json").resolve()
    resources_path = (args.stage_config or root / "hpc/config/stage_resources.json").resolve()
    site = load_json(site_path)
    stage_map = load_json(resources_path).get("stages", {})
    if args.stage not in stage_map:
        raise SystemExit(f"Unknown stage {args.stage!r}; choose from {sorted(stage_map)}")
    stage = stage_map[args.stage]
    _, rows = read_manifest(manifest)
    if not rows:
        raise SystemExit("Manifest has no data rows")
    array_spec = f"1-{len(rows)}"
    selected_task_count = len(rows)
    if args.array_range:
        match = re.fullmatch(r"(\d+)(?:-(\d+))?", args.array_range.strip())
        if not match:
            raise SystemExit("--array-range must be one index or one inclusive range such as 1-10")
        start = int(match.group(1))
        end = int(match.group(2) or match.group(1))
        if start < 1 or end < start or end > len(rows):
            raise SystemExit(
                f"--array-range {args.array_range!r} is outside the one-based manifest range 1-{len(rows)}"
            )
        array_spec = str(start) if start == end else f"{start}-{end}"
        selected_task_count = end - start + 1
    max_concurrent = args.max_concurrent or int(stage["max_concurrent"])
    if max_concurrent < 1:
        raise SystemExit("--max-concurrent must be >= 1")

    orca_exe = str(site.get("orca_executable") or "")
    if not Path(orca_exe).is_absolute():
        raise SystemExit("Edit hpc/config/site_config.json: orca_executable must be an absolute Linux path")
    python_exe = str(site.get("python_executable") or "python3")
    module_commands = [str(command) for command in site.get("module_commands", [])]
    optional = optional_directives(args.scheduler, site, args.dependency)

    logs = root / "state/scheduler_logs" / safe_component(args.stage)
    generated = root / "state/generated_submissions"
    logs.mkdir(parents=True, exist_ok=True)
    generated.mkdir(parents=True, exist_ok=True)
    job_name = safe_component(args.job_name or f"orca_{args.stage}")[:64]
    extension = "sbatch" if args.scheduler == "slurm" else "pbs"
    template = root / "hpc/templates" / f"{args.scheduler}_array.{extension}.in"
    if not template.is_file():
        raise SystemExit(f"Missing scheduler template: {template}")
    timestamp = utc_now().replace(":", "").replace("+", "_")
    script_path = generated / f"{job_name}_{args.scheduler}_{timestamp}.{extension}"
    values = {
        "JOB_NAME": job_name,
        "CPUS": str(int(stage["cpus"])),
        "MEMORY_MB": str(int(stage["memory_mb"])),
        "WALLTIME": str(stage["walltime"]),
        "ARRAY_SPEC": array_spec,
        "MAX_CONCURRENT": str(max_concurrent),
        "SCHEDULER_LOG_DIR": str(logs),
        "OPTIONAL_DIRECTIVES": "\n".join(optional) if optional else "# no optional scheduler directives",
        "MODULE_COMMANDS": "\n".join(module_commands) if module_commands else ": # no module commands configured",
        "ORCA_EXE": shell_quote(orca_exe),
        "SITE_CONFIG": shell_quote(site_path),
        "PACKAGE_ROOT": shell_quote(root),
        "MANIFEST": shell_quote(manifest),
        "PYTHON_EXE": shell_quote(python_exe),
        "STAGE": shell_quote(args.stage),
    }
    content = render(template.read_text(encoding="utf-8"), values)
    script_path.write_text(content, encoding="utf-8", newline="\n")
    try:
        script_path.chmod(0o750)
    except OSError:
        pass

    record: dict[str, Any] = {
        "schema_version": "1.0",
        "notice": "CALCULATED MODEL WORKFLOW - NOT EXPERIMENTAL DATA",
        "timestamp": utc_now(),
        "mode": "SUBMIT" if args.submit else "RENDER_ONLY",
        "scheduler": args.scheduler,
        "stage": args.stage,
        "manifest": str(manifest.relative_to(root)),
        "manifest_sha256": sha256_file(manifest),
        "task_count": len(rows),
        "array_spec": array_spec,
        "selected_task_count": selected_task_count,
        "max_concurrent": max_concurrent,
        "cpus_per_task": int(stage["cpus"]),
        "memory_mb_per_task": int(stage["memory_mb"]),
        "walltime": stage["walltime"],
        "dependency": args.dependency,
        "script": str(script_path.relative_to(root)),
        "script_sha256": sha256_file(script_path),
        "preflight_skipped": args.skip_preflight,
    }
    record_path = script_path.with_suffix(script_path.suffix + ".json")

    if not args.submit:
        atomic_json(record_path, record)
        print(f"Rendered only: {script_path}")
        print(f"Array: {array_spec}%{max_concurrent}; selected={selected_task_count}; cpus/task={stage['cpus']}; memory/task={stage['memory_mb']} MB")
        print("Inspect the script, run preflight, then repeat with --submit.")
        return 0

    if not args.skip_preflight:
        run_preflight(root, manifest, args.stage, args.scheduler, site_path, resources_path)
    command = ["sbatch", str(script_path)] if args.scheduler == "slurm" else ["qsub", str(script_path)]
    if shutil.which(command[0]) is None:
        raise SystemExit(f"Scheduler submission command is unavailable: {command[0]}")
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    record.update(
        {
            "submit_command": command,
            "submit_return_code": result.returncode,
            "scheduler_stdout": result.stdout.strip(),
            "scheduler_stderr": result.stderr.strip(),
        }
    )
    atomic_json(record_path, record)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        return result.returncode
    print(result.stdout.strip())
    print(f"Submission provenance: {record_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
