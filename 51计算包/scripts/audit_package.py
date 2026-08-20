#!/usr/bin/env python3
"""Static, fail-closed audit for the portable 51-compound ORCA package."""

from __future__ import annotations

import argparse
import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path

from package_io import atomic_write_json, sha256_file


EXPECTED_IDS = {f"F3_{i:03d}" for i in range(1, 53)} - {"F3_005"}


class Audit:
    def __init__(self): self.errors: list[str] = []; self.warnings: list[str] = []; self.checks: list[dict[str, object]] = []
    def require(self, condition: bool, label: str, detail: str = ""):
        self.checks.append({"label": label, "pass": bool(condition), "detail": detail})
        if not condition: self.errors.append(label + (f": {detail}" if detail else ""))
    def warn(self, condition: bool, message: str):
        if not condition: self.warnings.append(message)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle: return list(csv.DictReader(handle))


def audit_resources(root: Path, audit: Audit):
    path = root / "hpc" / "resource_profiles.json"
    audit.require(path.is_file(), "resource profiles exist")
    if not path.is_file(): return
    data = json.loads(path.read_text(encoding="utf-8")); fraction = float(data["policy"]["memory_fraction_limit"])
    audit.require(0 < fraction <= 0.75, "memory fraction is conservative", str(fraction))
    for name, p in data["profiles"].items():
        required = {"scheduler_cpus", "orca_nprocs", "maxcore_mb", "requested_memory_mb", "scratch_mb", "walltime"}
        audit.require(required <= set(p), f"profile {name} has required fields")
        if required <= set(p):
            ratio = int(p["orca_nprocs"]) * int(p["maxcore_mb"]) / int(p["requested_memory_mb"])
            audit.require(ratio <= fraction + 1e-12, f"profile {name} passes memory inequality", f"ratio={ratio:.3f}, limit={fraction:.3f}")
            audit.require(int(p["scheduler_cpus"]) >= int(p["orca_nprocs"]), f"profile {name} scheduler CPUs cover ORCA ranks")
    for name in ("revdsd_pilot", "revdsd_medium", "revdsd_large_4rank"):
        audit.require(name in data["profiles"], f"required profile {name} exists")


def audit_manifest(root: Path, audit: Audit, full_hash: bool):
    path = root / "data" / "GOAT_explicit_job_manifest.csv"
    audit.require(path.is_file(), "GOAT source manifest exists")
    if not path.is_file(): return
    rows = read_csv(path); audit.require(bool(rows), "GOAT source manifest is nonempty")
    ids = {r.get("figure_id", "") for r in rows}
    audit.require(ids == EXPECTED_IDS, "manifest contains exactly 51 expected compounds", f"missing={sorted(EXPECTED_IDS-ids)}, extra={sorted(ids-EXPECTED_IDS)}")
    audit.require("F3_005" not in ids, "benzamide F3_005 is excluded")
    job_ids = [r.get("job_id", "") for r in rows]
    audit.require(len(job_ids) == len(set(job_ids)) and all(job_ids), "manifest job IDs are unique and nonblank")
    required = {"figure_id", "scenario_id", "target_success", "max_candidates", "wave", "job_directory", "input_file", "coordinate_file"}
    audit.require(required <= set(rows[0]), "manifest has adaptive/HPC fields", f"missing={sorted(required-set(rows[0]))}")
    to_check = rows if full_hash else rows[: min(100, len(rows))]
    for r in to_check:
        rel = r.get("relative_job_directory") or r.get("portable_job_directory") or r["job_directory"]
        job_dir = (root / rel).resolve()
        try: job_dir.relative_to(root)
        except ValueError: audit.errors.append(f"job path escapes package: {rel}"); continue
        inp, xyz = job_dir / r["input_file"], job_dir / r["coordinate_file"]
        if not inp.is_file() or not xyz.is_file(): audit.errors.append(f"missing input bundle: {r['job_id']}"); continue
        expected_inp = (r.get("expected_input_sha256") or r.get("goat_input_sha256") or "").lower()
        expected_xyz = (r.get("expected_coordinate_sha256") or r.get("start_xyz_sha256") or "").lower()
        if expected_inp and sha256_file(inp).lower() != expected_inp: audit.errors.append(f"input hash mismatch: {r['job_id']}")
        if expected_xyz and sha256_file(xyz).lower() != expected_xyz: audit.errors.append(f"coordinate hash mismatch: {r['job_id']}")
    audit.checks.append({"label": "input bundle existence/hash scan", "pass": not any("hash mismatch" in x or "missing input bundle" in x for x in audit.errors), "detail": f"checked {len(to_check)}/{len(rows)} jobs"})
    initial_counts = {}
    for r in rows:
        key = (r["figure_id"], r["scenario_id"]); initial_counts.setdefault(key, 0)
        if str(r.get("initial_release", "")).lower() in {"1", "true", "yes"}: initial_counts[key] += 1
    for key, count in initial_counts.items():
        target = int(next(r["target_success"] for r in rows if (r["figure_id"], r["scenario_id"]) == key))
        maximum = int(next(r["max_candidates"] for r in rows if (r["figure_id"], r["scenario_id"]) == key))
        if not target <= count <= maximum: audit.errors.append(f"unsafe initial release count {key}: target={target}, release={count}, max={maximum}")


def audit_code(root: Path, audit: Audit):
    required = [
        "scripts/adaptive_goat_controller.py", "scripts/run_orca_job.py", "scripts/classify_failure.py",
        "scripts/write_package_checksums.py", "scripts/verify_receipt.py", "scripts/prepare_hpc_jobs.py",
        "scripts/sync_scheduler_states.py",
        "hpc/slurm/run_manifest_task.sbatch.in", "hpc/site_config.example.json",
    ]
    for rel in required: audit.require((root / rel).is_file(), f"required control file exists: {rel}")
    forbidden = re.compile(r"\b(?:mpirun|mpiexec|srun)\b[^\n]{0,160}\borca\b", re.IGNORECASE)
    for path in list((root / "scripts").glob("*.py")) + list((root / "hpc" / "slurm").glob("*")):
        if path.is_file():
            source = path.read_text(encoding="utf-8", errors="replace")
            audit.require(not bool(forbidden.search(source)), f"no outer MPI ORCA invocation in {path.relative_to(root)}")
            if path.name in {"adaptive_goat_controller.py", "sync_scheduler_states.py"}:
                audit.require("scancel" not in source.lower(), f"controller code does not invoke scheduler cancellation: {path.relative_to(root)}")
    runner = (root / "scripts" / "run_orca_job.py").read_text(encoding="utf-8", errors="replace") if (root / "scripts" / "run_orca_job.py").is_file() else ""
    for marker in ("Heartbeat", "publish_results", "sha256_file", "node-local scratch", "os.replace"):
        audit.require(marker in runner, f"runner contains {marker} safeguard")


def audit_site(path: Path, audit: Audit, live: bool):
    audit.require(path.is_file(), "site config exists", str(path))
    if not path.is_file(): return
    data = json.loads(path.read_text(encoding="utf-8")); orca_text = str(data.get("orca_executable", "")); orca = Path(orca_text)
    # A portable package may be audited on Windows while carrying a POSIX cluster path.
    audit.require(orca.is_absolute() or orca_text.startswith("/"), "ORCA executable is an absolute path", orca_text)
    audit.require("6.1.1" in str(data.get("orca_version_required", "")), "ORCA version lock is 6.1.1")
    if live: audit.require(orca.is_file(), "live ORCA executable exists", str(orca))
    env = data.get("environment", {})
    for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        audit.require(str(env.get(name)) == "1", f"site config fixes {name}=1")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--site-config", type=Path)
    parser.add_argument("--live-site", action="store_true")
    parser.add_argument("--full-hash", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(); root = args.package_root.resolve(); audit = Audit()
    for directory in ("config", "data", "01_goat_explicit", "hpc", "scripts", "gates", "package_metadata"):
        audit.require((root / directory).is_dir(), f"required directory exists: {directory}")
    audit_resources(root, audit); audit_manifest(root, audit, args.full_hash); audit_code(root, audit)
    site = args.site_config or (root / "hpc" / "site_config.example.json")
    audit_site(site.resolve(), audit, args.live_site)
    for schema in (root / "package_metadata" / "schemas").glob("*.json"):
        try: json.loads(schema.read_text(encoding="utf-8")); audit.checks.append({"label": f"valid JSON: {schema.name}", "pass": True})
        except Exception as exc: audit.errors.append(f"invalid JSON schema {schema.name}: {exc}")
    result = {"schema_version": "1.0", "created_utc": datetime.now(timezone.utc).isoformat(), "package_root": str(root),
              "status": "PASS" if not audit.errors else "FAIL", "errors": audit.errors, "warnings": audit.warnings, "checks": audit.checks,
              "scope": "static package audit; no ORCA calculation was executed"}
    output = args.output or (root / "package_metadata" / "audit_report.json")
    atomic_write_json(output, result); print(json.dumps(result, ensure_ascii=False, indent=2))
    return 1 if audit.errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
