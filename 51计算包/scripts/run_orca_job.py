#!/usr/bin/env python3
"""Run one ORCA job with local scratch, heartbeat and hash-verified stage-out."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import socket
import subprocess
import threading
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

from package_io import atomic_write_json, load_json, sha256_file


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def validate_profile(inp: Path, profile: dict[str, object], fraction: float) -> None:
    text = inp.read_text(encoding="utf-8", errors="replace")
    nprocs = int(profile["orca_nprocs"])
    maxcore = int(profile["maxcore_mb"])
    memory = int(profile["requested_memory_mb"])
    if nprocs * maxcore > fraction * memory:
        raise ValueError(f"unsafe memory profile: {nprocs}*{maxcore} > {fraction}*{memory}")
    pal = re.search(r"%pal\s+.*?nprocs\s+(\d+)", text, re.IGNORECASE | re.DOTALL)
    one_line = re.search(r"!.*?PAL(\d+)\b", text, re.IGNORECASE)
    actual_nprocs = int(pal.group(1)) if pal else (int(one_line.group(1)) if one_line else None)
    if actual_nprocs != nprocs:
        raise ValueError(f"ORCA input nprocs={actual_nprocs}; profile requires {nprocs}")
    match = re.search(r"%maxcore\s+(\d+)", text, re.IGNORECASE)
    if not match or int(match.group(1)) != maxcore:
        raise ValueError(f"ORCA input maxcore does not equal profile value {maxcore} MB")


def choose_scratch(site: dict[str, object]) -> Path:
    for name in site.get("node_local_scratch_environment_variables", []):
        value = os.environ.get(str(name), "").strip()
        if value:
            root = Path(os.path.expandvars(value)).resolve()
            if root.is_dir() and os.access(root, os.W_OK):
                return root
    if site.get("allow_fallback_scratch") is True:
        root = Path(os.path.expandvars(str(site["fallback_scratch_root"]))).resolve()
        root.mkdir(parents=True, exist_ok=True)
        if root.is_dir() and os.access(root, os.W_OK):
            return root
    raise RuntimeError("no writable node-local scratch was found; shared-filesystem fallback is disabled")


def input_files(job_dir: Path, inp: Path) -> list[Path]:
    bundle = job_dir / "input_bundle.json"
    if bundle.is_file():
        names = load_json(bundle).get("files", [])
        files = [(job_dir / str(name)).resolve() for name in names]
    else:
        allowed = {".inp", ".xyz", ".gbw", ".hess", ".engrad", ".pc", ".pot"}
        files = [p.resolve() for p in job_dir.iterdir() if p.is_file() and p.suffix.lower() in allowed]
    if inp.resolve() not in files:
        files.append(inp.resolve())
    for path in files:
        if path.parent != job_dir.resolve() or not path.is_file():
            raise ValueError(f"input bundle has missing or non-local file: {path}")
    return sorted(set(files))


def acquire_lock(job_dir: Path, site: dict[str, object], clearance: Path | None) -> tuple[Path, str]:
    lock = job_dir / ".runner.lock.json"
    token = uuid.uuid4().hex
    payload = {"schema_version": "1.0", "token": token, "created_utc": utc_now(), "hostname": socket.gethostname(),
               "pid": os.getpid(), "scheduler_job_id": os.environ.get("SLURM_JOB_ID", "")}
    try:
        fd = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL)
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2); handle.write("\n"); handle.flush(); os.fsync(handle.fileno())
        return lock, token
    except FileExistsError:
        old = load_json(lock)
        heartbeat = job_dir / "heartbeat.json"
        age = time.time() - heartbeat.stat().st_mtime if heartbeat.exists() else time.time() - lock.stat().st_mtime
        if age <= int(site.get("stale_lock_seconds", 900)):
            raise RuntimeError(f"live lock exists (heartbeat age {age:.0f}s); refusing concurrent run")
        if clearance is None or not clearance.is_file():
            raise RuntimeError("stale lock found; a manual stale-lock clearance receipt is required")
        approval = load_json(clearance)
        terminal = {str(x).upper() for x in site.get("terminal_scheduler_states", [])}
        if approval.get("approved") is not True or str(approval.get("prior_scheduler_state", "")).upper() not in terminal:
            raise RuntimeError("stale-lock clearance is not approved or does not prove a terminal scheduler state")
        if approval.get("prior_lock_token") != old.get("token"):
            raise RuntimeError("stale-lock clearance token does not match the lock")
        archive = job_dir / f".runner.lock.stale.{int(time.time())}.json"
        os.replace(lock, archive)
        return acquire_lock(job_dir, site, None)


class Heartbeat:
    def __init__(self, path: Path, token: str, seconds: int):
        self.path, self.token, self.seconds = path, token, max(10, seconds)
        self.stop_event = threading.Event()
        self.thread = threading.Thread(target=self._run, daemon=True)

    def _run(self):
        while not self.stop_event.is_set():
            atomic_write_json(self.path, {"schema_version": "1.0", "token": self.token, "updated_utc": utc_now(),
                                          "hostname": socket.gethostname(), "pid": os.getpid(),
                                          "scheduler_job_id": os.environ.get("SLURM_JOB_ID", "")})
            self.stop_event.wait(self.seconds)

    def __enter__(self):
        self.thread.start(); return self

    def __exit__(self, *_):
        self.stop_event.set(); self.thread.join(timeout=self.seconds + 2)


def publish_results(scratch_job: Path, job_dir: Path, receipt: dict[str, object]) -> Path:
    attempts = job_dir / "attempts"
    attempts.mkdir(exist_ok=True)
    attempt_id = str(receipt["attempt_id"])
    temporary = attempts / f".{attempt_id}.stageout"
    final = attempts / attempt_id
    if temporary.exists() or final.exists():
        raise FileExistsError(f"attempt destination already exists: {attempt_id}")
    temporary.mkdir()
    hashes = {}
    for source in sorted(p for p in scratch_job.rglob("*") if p.is_file()):
        rel = source.relative_to(scratch_job)
        destination = temporary / rel
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        source_hash, dest_hash = sha256_file(source), sha256_file(destination)
        if source_hash != dest_hash:
            raise IOError(f"stage-out hash mismatch: {rel}")
        hashes[rel.as_posix()] = source_hash
    receipt["files_sha256"] = hashes
    atomic_write_json(temporary / "runner_receipt.json", receipt)
    os.replace(temporary, final)
    return final


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--site-config", type=Path, required=True)
    parser.add_argument("--resource-profiles", type=Path, required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--job-dir", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--stale-lock-clearance", type=Path)
    args = parser.parse_args()
    root, job_dir, inp = args.package_root.resolve(), args.job_dir.resolve(), args.input.resolve()
    job_dir.relative_to(root)
    if inp.parent != job_dir or not inp.is_file():
        raise SystemExit("input must be an existing file directly inside job-dir")
    site = load_json(args.site_config.resolve())
    resources = load_json(args.resource_profiles.resolve())
    if args.profile not in resources["profiles"]:
        raise SystemExit(f"unknown profile {args.profile}")
    profile = resources["profiles"][args.profile]
    validate_profile(inp, profile, float(resources["policy"]["memory_fraction_limit"]))
    orca = Path(str(site["orca_executable"]))
    if not orca.is_absolute() or not orca.is_file():
        raise SystemExit("site_config.orca_executable must be the existing absolute ORCA driver path")
    lock, token = acquire_lock(job_dir, site, args.stale_lock_clearance.resolve() if args.stale_lock_clearance else None)
    attempt_id = f"attempt_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}_{uuid.uuid4().hex[:8]}"
    scratch_job = choose_scratch(site) / f"orca_nmr_{attempt_id}"
    scratch_job.mkdir(parents=False, exist_ok=False)
    receipt: dict[str, object] = {"schema_version": "1.0", "attempt_id": attempt_id, "started_utc": utc_now(),
                                  "profile": args.profile, "orca_executable": str(orca), "orca_version_required": site.get("orca_version_required"),
                                  "command": [str(orca), inp.name], "outer_mpi_launcher": False, "scratch_path": str(scratch_job),
                                  "scheduler_job_id": os.environ.get("SLURM_JOB_ID", "")}
    state_path = job_dir / "job_state.json"
    try:
        for source in input_files(job_dir, inp):
            shutil.copy2(source, scratch_job / source.name)
        env = os.environ.copy()
        for key, value in site.get("environment", {}).items(): env[str(key)] = str(value)
        for key, value in site.get("fabric_environment", {}).items(): env[str(key)] = str(value)
        for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
            env[key] = "1"
        atomic_write_json(state_path, {"schema_version": "1.0", "state": "RUNNING", "attempt_id": attempt_id, "updated_utc": utc_now(), "accepted": False})
        stdout_file = scratch_job / f"{inp.stem}.out"
        with Heartbeat(job_dir / "heartbeat.json", token, int(site.get("heartbeat_seconds", 60))):
            with stdout_file.open("w", encoding="utf-8", errors="replace") as output:
                proc = subprocess.run([str(orca), inp.name], cwd=scratch_job, env=env, stdout=output, stderr=subprocess.STDOUT, check=False)
        text = stdout_file.read_text(encoding="utf-8", errors="replace")
        normal = "ORCA TERMINATED NORMALLY" in text
        receipt.update(finished_utc=utc_now(), returncode=proc.returncode, normal_termination_marker=normal,
                       status="COMPLETED" if proc.returncode == 0 and normal else "FAILED")
        final = publish_results(scratch_job, job_dir, receipt)
        atomic_write_json(state_path, {"schema_version": "1.0", "state": receipt["status"], "attempt_id": attempt_id,
                                      "result_dir": str(final), "updated_utc": utc_now(), "accepted": False,
                                      "note": "Scientific acceptance is set only by downstream gates."})
        return 0 if receipt["status"] == "COMPLETED" else 1
    except Exception as exc:
        receipt.update(finished_utc=utc_now(), status="FAILED_RUNNER", exception=repr(exc))
        try:
            publish_results(scratch_job, job_dir, receipt)
        except Exception as stage_exc:
            atomic_write_json(job_dir / f"runner_emergency_receipt_{attempt_id}.json", {**receipt, "stageout_exception": repr(stage_exc)})
        atomic_write_json(state_path, {"schema_version": "1.0", "state": "FAILED", "attempt_id": attempt_id,
                                      "updated_utc": utc_now(), "accepted": False, "runner_exception": repr(exc)})
        raise
    finally:
        if lock.exists():
            try:
                current = load_json(lock)
                if current.get("token") == token: lock.unlink()
            except Exception:
                pass
        if scratch_job.exists():
            # The exact per-attempt directory was created above under a verified scratch root.
            shutil.rmtree(scratch_job, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
