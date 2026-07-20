#!/usr/bin/env python3
"""Fail-fast server and input audit before the 24 h clock is spent."""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

import psutil

from common import (
    ROOT,
    physical_cpu_ids,
    read_csv,
    read_json,
    runtime_environment_identity,
    sha256,
    utc_stamp,
    write_json,
)


def check(condition: bool, name: str, detail: str, critical: bool = True) -> dict[str, object]:
    return {
        "name": name,
        "status": "PASS" if condition else "FAIL" if critical else "WARN",
        "critical": critical,
        "detail": detail,
    }


def orca_version(orca: Path) -> str:
    for arguments in ([str(orca), "-v"], [str(orca), "--version"]):
        try:
            result = subprocess.run(arguments, capture_output=True, text=True, timeout=20)
        except Exception:
            continue
        text = (result.stdout + "\n" + result.stderr).strip()
        if text:
            return text[:4000]
    return "VERSION_QUERY_RETURNED_NO_TEXT"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orca", type=Path, required=True, help="Absolute path to the ORCA driver")
    parser.add_argument("--scratch", type=Path, required=True, help="Node-local NVMe scratch directory")
    parser.add_argument("--production", action="store_true", help="Enforce the 104-core server thresholds")
    parser.add_argument("--report", type=Path, default=ROOT / "results" / "preflight.json")
    args = parser.parse_args()

    hardware = read_json(ROOT / "config" / "hardware_104physical_24h.json")
    cap = int(hardware["resource_limits"]["physical_core_cap"])
    checks: list[dict[str, object]] = []
    checks.append(check(args.orca.is_absolute(), "orca_absolute_path", str(args.orca)))
    checks.append(check(args.orca.is_file(), "orca_executable_exists", str(args.orca)))
    checks.append(check(os.access(args.orca, os.X_OK), "orca_executable_permission", str(args.orca)))
    version = orca_version(args.orca) if args.orca.is_file() else "ORCA_MISSING"
    checks.append(check("6.1" in version, "orca_6_1_version", version, critical=args.production))
    checks.append(check(os.name != "nt", "linux_runtime", platform.platform(), critical=args.production))

    cpu_ids = physical_cpu_ids(cap)
    minimum_cores = 90 if args.production else 1
    checks.append(
        check(
            len(cpu_ids) >= minimum_cores,
            "physical_core_affinity",
            f"selected={len(cpu_ids)} cap={cap} ids={cpu_ids}",
            critical=args.production,
        )
    )
    memory_gib = psutil.virtual_memory().total / (1024**3)
    checks.append(
        check(
            memory_gib >= (700 if args.production else 2),
            "memory_capacity",
            f"detected={memory_gib:.1f} GiB",
            critical=args.production,
        )
    )
    args.scratch.mkdir(parents=True, exist_ok=True)
    free_gib = shutil.disk_usage(args.scratch).free / (1024**3)
    checks.append(
        check(
            free_gib >= (300 if args.production else 1),
            "scratch_free_space",
            f"path={args.scratch.resolve()} free={free_gib:.1f} GiB",
            critical=args.production,
        )
    )
    scratch_device = "UNKNOWN"
    resolved_scratch = args.scratch.resolve()
    matching_mounts = []
    for partition in psutil.disk_partitions(all=True):
        try:
            resolved_scratch.relative_to(Path(partition.mountpoint).resolve())
        except ValueError:
            continue
        matching_mounts.append(partition)
    if matching_mounts:
        selected_mount = max(matching_mounts, key=lambda item: len(str(Path(item.mountpoint))))
        scratch_device = f"device={selected_mount.device} mount={selected_mount.mountpoint} fstype={selected_mount.fstype}"
    checks.append(
        check(
            "nvme" in scratch_device.lower(),
            "scratch_device_nvme_hint",
            scratch_device + "; verify with the site administrator if RAID/LVM hides the NVMe name",
            critical=False,
        )
    )
    checks.append(check(shutil.which("taskset") is not None, "taskset_available", str(shutil.which("taskset")), critical=args.production))

    manifest = read_csv(ROOT / "config" / "input_manifest_sha256.csv")
    mismatches: list[str] = []
    missing: list[str] = []
    for row in manifest:
        path = ROOT / row["path"]
        if not path.is_file():
            missing.append(row["path"])
        elif sha256(path) != row["sha256"]:
            mismatches.append(row["path"])
    checks.append(
        check(
            not missing and not mismatches,
            "frozen_input_hashes",
            f"files={len(manifest)} missing={missing} mismatched={mismatches}",
        )
    )

    species = read_csv(ROOT / "config" / "primary_species.csv")
    checks.append(check(len(species) == 52, "compound_count", f"observed={len(species)} expected=52"))
    environment_issues: list[str] = []
    for row in species:
        atom_map = read_csv(ROOT / "structures" / row["figure_id"] / "atom_map.csv")
        if any(atom["environment_atom"].lower() == "true" or atom["component"] != "analyte" for atom in atom_map):
            environment_issues.append(row["figure_id"])
    checks.append(
        check(
            not environment_issues,
            "implicit_only_structures",
            f"structures_with_environment_atoms={environment_issues}",
        )
    )

    report = {
        "schema_version": "1.0",
        "generated_utc": utc_stamp(),
        "production_mode": args.production,
        "python": sys.version,
        "platform": platform.platform(),
        "runtime_environment_identity": runtime_environment_identity(),
        "orca_path": str(args.orca.resolve()),
        "orca_version_query": version,
        "scratch": str(args.scratch.resolve()),
        "selected_physical_cpu_ids": cpu_ids,
        "checks": checks,
        "pass": all(item["status"] != "FAIL" for item in checks),
    }
    write_json(args.report, report)
    for item in checks:
        print(f"[{item['status']}] {item['name']}: {item['detail']}")
    print(f"Preflight report: {args.report}")
    return 0 if report["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
