"""Hash-locked, single-instance runner for approved explicit-cluster opt/freq jobs.

Dry-run validation is the default.  ORCA is started only with ``--execute``.
One job is run at a time, with pre-start CPU/RAM/disk/process gates and a
finish-current-then-stop policy if sustained CPU load exceeds the reviewed
threshold.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import json
import os
import shutil
import subprocess
import time
import uuid
from ctypes import wintypes
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from workflow_schema import (
    NOTICE,
    WorkflowValidationError,
    atomic_write_json,
    inline_orca_input_bundle_sha256,
    normalize_sha256,
    read_json,
    sha256_file,
    validate_completed_opt_freq_state,
)


PROJECT = Path(__file__).resolve().parents[1]
NORMAL_MARKER = "ORCA TERMINATED NORMALLY"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


class MEMORYSTATUSEX(ctypes.Structure):
    _fields_ = [
        ("dwLength", wintypes.DWORD),
        ("dwMemoryLoad", wintypes.DWORD),
        ("ullTotalPhys", ctypes.c_ulonglong),
        ("ullAvailPhys", ctypes.c_ulonglong),
        ("ullTotalPageFile", ctypes.c_ulonglong),
        ("ullAvailPageFile", ctypes.c_ulonglong),
        ("ullTotalVirtual", ctypes.c_ulonglong),
        ("ullAvailVirtual", ctypes.c_ulonglong),
        ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
    ]


def available_ram_gib() -> float:
    status = MEMORYSTATUSEX()
    status.dwLength = ctypes.sizeof(status)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
        raise OSError("GlobalMemoryStatusEx failed")
    return status.ullAvailPhys / (1024**3)


def _filetime_value(value: wintypes.FILETIME) -> int:
    return (int(value.dwHighDateTime) << 32) + int(value.dwLowDateTime)


def _system_times() -> tuple[int, int, int]:
    idle = wintypes.FILETIME()
    kernel = wintypes.FILETIME()
    user = wintypes.FILETIME()
    if not ctypes.windll.kernel32.GetSystemTimes(
        ctypes.byref(idle), ctypes.byref(kernel), ctypes.byref(user)
    ):
        raise OSError("GetSystemTimes failed")
    return _filetime_value(idle), _filetime_value(kernel), _filetime_value(user)


def cpu_percent(sample_seconds: float) -> float:
    before = _system_times()
    time.sleep(max(0.2, sample_seconds))
    after = _system_times()
    idle = after[0] - before[0]
    total = (after[1] - before[1]) + (after[2] - before[2])
    return 0.0 if total <= 0 else max(0.0, min(100.0, 100.0 * (1.0 - idle / total)))


def running_process_names() -> set[str]:
    completed = subprocess.run(
        ["tasklist", "/FO", "CSV", "/NH"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=20,
        check=False,
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
    )
    if completed.returncode != 0:
        raise WorkflowValidationError(f"tasklist failed: {completed.stderr.strip()}")
    names: set[str] = set()
    for row in csv.reader(completed.stdout.splitlines()):
        if row:
            names.add(Path(row[0]).stem.casefold())
    return names


class ExclusiveLock:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.descriptor: int | None = None

    def __enter__(self) -> "ExclusiveLock":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        try:
            self.descriptor = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError as exc:
            raise WorkflowValidationError(f"Another opt/freq runner lock exists: {self.path}") from exc
        os.write(self.descriptor, f"pid={os.getpid()} started={utc_now()}\n".encode("ascii"))
        return self

    def __exit__(self, *_args: object) -> None:
        if self.descriptor is not None:
            os.close(self.descriptor)
        self.path.unlink(missing_ok=True)


def _read_manifest(path: Path) -> list[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle)]
    if not rows:
        raise WorkflowValidationError(f"Empty opt/freq manifest: {path}")
    job_ids = [row.get("job_id", "") for row in rows]
    if any(not value for value in job_ids) or len(job_ids) != len(set(job_ids)):
        raise WorkflowValidationError("Opt/freq manifest has blank or duplicate job IDs")
    return rows


def validate_manifest_job(row: dict[str, str], project: Path) -> dict[str, object]:
    job_id = row.get("job_id", "")
    job_dir = Path(row.get("job_directory", "")).resolve()
    try:
        job_dir.relative_to(project.resolve())
    except ValueError as exc:
        raise WorkflowValidationError(f"Opt/freq job escapes project: {job_id}") from exc
    input_path = job_dir / row.get("input_file", "")
    if not input_path.is_file():
        raise WorkflowValidationError(f"Missing opt/freq input: {input_path}")
    input_bytes = input_path.read_bytes()
    if b"xyzfile" in input_bytes.lower():
        raise WorkflowValidationError(f"Opt/freq hand-off requires inline, map-locked coordinates: {input_path}")
    input_hash = sha256_file(input_path)
    expected_input = normalize_sha256(row.get("input_sha256", ""), "opt/freq input_sha256")
    if input_hash != expected_input:
        raise WorkflowValidationError(f"Opt/freq input hash mismatch: {job_id}")
    bundle_hash = inline_orca_input_bundle_sha256(input_path.name, input_bytes)
    expected_bundle = normalize_sha256(
        row.get("expected_input_bundle_sha256", ""), "opt/freq expected_input_bundle_sha256"
    )
    if bundle_hash != expected_bundle:
        raise WorkflowValidationError(f"Opt/freq input bundle hash mismatch: {job_id}")
    provenance_path = job_dir / "provenance.json"
    provenance = read_json(provenance_path)
    if provenance.get("job_id") != job_id or provenance.get("scenario_id") != row.get("scenario_id"):
        raise WorkflowValidationError(f"Opt/freq provenance identity mismatch: {job_id}")
    if normalize_sha256(str(provenance.get("dft_input_sha256", ""))) != input_hash:
        raise WorkflowValidationError(f"Opt/freq provenance input hash mismatch: {job_id}")
    if normalize_sha256(str(provenance.get("dft_input_bundle_sha256", ""))) != bundle_hash:
        raise WorkflowValidationError(f"Opt/freq provenance bundle hash mismatch: {job_id}")
    dependency_pairs = (
        ("approved_gate_file", "approved_gate_sha256"),
        ("source_state_file", "source_state_sha256"),
        ("source_finalensemble_file", "source_finalensemble_sha256"),
        ("selected_coordinate_source_file", "selected_coordinate_sha256"),
        ("source_atom_map_file", "source_atom_map_sha256"),
        ("source_isotope_map_file", "source_isotope_map_sha256"),
    )
    for path_field, hash_field in dependency_pairs:
        dependency = Path(str(provenance.get(path_field, ""))).resolve()
        try:
            dependency.relative_to(project.resolve())
        except ValueError as exc:
            raise WorkflowValidationError(
                f"Opt/freq provenance dependency escapes project ({path_field}): {job_id}"
            ) from exc
        if sha256_file(dependency) != normalize_sha256(str(provenance.get(hash_field, "")), hash_field):
            raise WorkflowValidationError(f"Opt/freq provenance dependency hash mismatch ({path_field}): {job_id}")
    for copied_name, hash_field in (
        ("approved_selected_cluster.xyz", "selected_coordinate_sha256"),
        ("source_atom_map.csv", "source_atom_map_sha256"),
        ("source_isotope_map.csv", "source_isotope_map_sha256"),
    ):
        if sha256_file(job_dir / copied_name) != normalize_sha256(
            str(provenance.get(hash_field, "")), hash_field
        ):
            raise WorkflowValidationError(f"Copied opt/freq dependency differs ({copied_name}): {job_id}")
    dependency_text = (job_dir / "SCIENTIFIC_DEPENDENCY_SHA256.txt").read_text(
        encoding="ascii"
    ).strip().lower()
    if dependency_text != normalize_sha256(
        str(provenance.get("scientific_dependency_sha256", "")),
        "scientific_dependency_sha256",
    ):
        raise WorkflowValidationError(f"Scientific-dependency hash record differs: {job_id}")
    return {
        "job_id": job_id,
        "job_directory": job_dir,
        "input_path": input_path,
        "input_sha256": input_hash,
        "input_bundle_sha256": bundle_hash,
        "provenance_path": provenance_path,
        "provenance_sha256": sha256_file(provenance_path),
    }


def _artifact_records(work: Path, stem: str) -> dict[str, dict[str, object]]:
    names = {
        "gbw": f"{stem}.gbw",
        "hess": f"{stem}.hess",
        "property": f"{stem}.property.txt",
        "optimized_xyz": f"{stem}.xyz",
    }
    records: dict[str, dict[str, object]] = {}
    for role, name in names.items():
        path = work / name
        exists = path.is_file() and path.stat().st_size > 0
        records[role] = {
            "path": str(path.resolve()),
            "exists_nonempty": exists,
            "size_bytes": path.stat().st_size if exists else 0,
            "sha256": sha256_file(path) if exists else "",
        }
    return records


def _completed_state_valid(
    state_path: Path,
    evidence: dict[str, object],
    row: dict[str, str],
    project: Path,
) -> bool:
    if not state_path.is_file():
        return False
    try:
        completed = validate_completed_opt_freq_state(row, state_path, project)
    except (OSError, ValueError, WorkflowValidationError):
        return False
    return completed["input_bundle_sha256"] == evidence["input_bundle_sha256"]


def prestart_gate(project: Path, resource_policy: dict) -> dict[str, object]:
    execution = resource_policy["execution"]
    stage = resource_policy["stage_resources"]["opt_freq"]
    start = resource_policy["start_gates"]
    blocked_names = {str(value).casefold() for value in start["block_if_external_process_names_match"]}
    running = running_process_names()
    blocked = sorted(running & blocked_names)
    if blocked:
        raise WorkflowValidationError(f"External quantum-chemistry process is already running: {blocked}")
    sample_seconds = float(start["cpu_sample_seconds"])
    cpu = cpu_percent(sample_seconds)
    ram = available_ram_gib()
    free_disk = shutil.disk_usage(project.anchor).free / (1024**3)
    if cpu > float(start["maximum_average_cpu_percent"]):
        raise WorkflowValidationError(f"Pre-start CPU {cpu:.1f}% exceeds gate")
    if ram < float(stage["minimum_free_ram_gib"]):
        raise WorkflowValidationError(f"Free RAM {ram:.1f} GiB is below opt/freq gate")
    required_disk = float(start["minimum_free_space_gib"].get(project.drive.rstrip(":"), 0))
    if free_disk < required_disk:
        raise WorkflowValidationError(f"Free disk {free_disk:.1f} GiB is below gate")
    return {
        "cpu_average_percent": cpu,
        "cpu_sample_seconds": sample_seconds,
        "free_ram_gib": ram,
        "free_project_drive_gib": free_disk,
        "mpi_ranks": int(stage["mpi_ranks"]),
        "maxcore_mb_per_rank": int(stage["maxcore_mb_per_rank"]),
    }


def run_one(
    row: dict[str, str],
    evidence: dict[str, object],
    project: Path,
    state_root: Path,
    run_root: Path,
    resource_policy: dict,
    orca_exe: Path,
) -> tuple[dict, bool]:
    job_id = str(evidence["job_id"])
    state_path = state_root / f"{job_id}.json"
    if _completed_state_valid(state_path, evidence, row, project):
        return read_json(state_path), False
    gate_snapshot = prestart_gate(project, resource_policy)
    attempt_id = f"attempt_{datetime.now().strftime('%Y%m%dT%H%M%S')}_{uuid.uuid4().hex[:8]}"
    work = run_root / job_id / attempt_id / "work"
    work.mkdir(parents=True, exist_ok=False)
    input_path = work / "dft_opt_freq.inp"
    shutil.copy2(Path(evidence["input_path"]), input_path)
    if sha256_file(input_path) != evidence["input_sha256"]:
        raise WorkflowValidationError(f"Copied opt/freq input hash changed before launch: {job_id}")
    output_path = work / "dft_opt_freq.out"
    stderr_path = work / "dft_opt_freq.stderr.log"
    previous = read_json(state_path) if state_path.is_file() else {}
    attempts = list(previous.get("attempts", [])) if isinstance(previous.get("attempts", []), list) else []
    attempt = {
        "attempt_no": len(attempts) + 1,
        "status": "running",
        "started_utc": utc_now(),
        "input_sha256": evidence["input_sha256"],
        "input_bundle_sha256": evidence["input_bundle_sha256"],
        "run_directory": str(work.resolve()),
        "output_file": str(output_path.resolve()),
        "stderr_file": str(stderr_path.resolve()),
        "prestart_resource_gate": gate_snapshot,
        "command": [str(orca_exe.resolve()), input_path.name],
        "orca_executable_sha256": sha256_file(orca_exe),
    }
    attempts.append(attempt)
    state = {
        "schema_version": 2,
        "notice": NOTICE,
        "job_id": job_id,
        "figure_id": row.get("figure_id", ""),
        "scenario_id": row.get("scenario_id", ""),
        "source_job_directory": str(evidence["job_directory"]),
        "status": "running",
        "input_bundle_sha256": evidence["input_bundle_sha256"],
        "provenance_sha256": evidence["provenance_sha256"],
        "resource_policy_sha256": sha256_file(project / "config" / "resource_policy.json"),
        "attempt_count": len(attempts),
        "attempts": attempts,
    }
    atomic_write_json(state_path, state)
    environment = os.environ.copy()
    environment.update({"OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
    temporary = work / "_tmp"
    temporary.mkdir()
    environment.update({"TEMP": str(temporary), "TMP": str(temporary)})
    mpi_bin = str(resource_policy.get("machine", {}).get("msmpi_bin", "")).strip()
    path_prefix = [str(orca_exe.parent)]
    if mpi_bin:
        path_prefix.append(mpi_bin)
    environment["PATH"] = os.pathsep.join(path_prefix) + os.pathsep + environment.get("PATH", "")
    creationflags = getattr(subprocess, "CREATE_NO_WINDOW", 0) | getattr(
        subprocess, "BELOW_NORMAL_PRIORITY_CLASS", 0
    )
    start_time = time.monotonic()
    stop_after_current = False
    sustained = 0.0
    heartbeat = float(resource_policy["runtime_policy"]["heartbeat_seconds"])
    high_threshold = float(resource_policy["runtime_policy"]["high_cpu_percent"])
    sustained_limit = float(resource_policy["runtime_policy"]["high_cpu_sustained_seconds"])
    with output_path.open("wb") as stdout_handle, stderr_path.open("wb") as stderr_handle:
        process = subprocess.Popen(
            [str(orca_exe.resolve()), input_path.name],
            cwd=work,
            env=environment,
            stdout=stdout_handle,
            stderr=stderr_handle,
            creationflags=creationflags,
        )
        while process.poll() is None:
            sampled = cpu_percent(heartbeat)
            sustained = sustained + heartbeat if sampled >= high_threshold else 0.0
            if sustained >= sustained_limit:
                stop_after_current = True
        return_code = int(process.returncode)
    elapsed = time.monotonic() - start_time
    artifacts = _artifact_records(work, input_path.stem)
    normal = output_path.is_file() and NORMAL_MARKER in output_path.read_text(
        encoding="utf-8", errors="replace"
    )
    source_unchanged = (
        sha256_file(Path(evidence["input_path"])) == evidence["input_sha256"]
    )
    complete = return_code == 0 and normal and source_unchanged and all(
        bool(record["exists_nonempty"]) for record in artifacts.values()
    )
    attempt.update(
        {
            "status": "completed" if complete else "failed",
            "ended_utc": utc_now(),
            "elapsed_seconds": round(elapsed, 2),
            "return_code": return_code,
            "normal_termination_marker": normal,
            "source_input_unchanged": source_unchanged,
            "output_sha256": sha256_file(output_path) if output_path.is_file() else "",
            "artifacts": artifacts,
            "stop_after_current_due_to_sustained_high_cpu": stop_after_current,
        }
    )
    state["status"] = "completed" if complete else "failed"
    state["reason"] = "" if complete else "return code, normal termination or required artifacts failed"
    atomic_write_json(state_path, state)
    return state, stop_after_current


def dry_run_validate(
    project: Path, manifest_path: Path, only: Iterable[str] | None = None
) -> list[dict[str, object]]:
    rows = _read_manifest(manifest_path)
    selected = set(only or [])
    if selected:
        rows = [row for row in rows if row["job_id"] in selected or row.get("figure_id") in selected]
    if not rows:
        raise WorkflowValidationError("No opt/freq jobs selected")
    return [validate_manifest_job(row, project) for row in rows]


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--project", type=Path, default=PROJECT)
    value.add_argument("--manifest", type=Path)
    value.add_argument("--only", nargs="*")
    value.add_argument("--execute", action="store_true")
    value.add_argument("--orca-exe", type=Path)
    value.add_argument("--state-root", type=Path)
    value.add_argument("--run-root", type=Path)
    return value


def main() -> None:
    args = parser().parse_args()
    project = args.project.resolve()
    manifest = (args.manifest or project / "data" / "DFT_explicit_job_manifest.csv").resolve()
    rows = _read_manifest(manifest)
    selected = set(args.only or [])
    if selected:
        rows = [row for row in rows if row["job_id"] in selected or row.get("figure_id") in selected]
    if not rows:
        raise WorkflowValidationError("No opt/freq jobs selected")
    evidence = [validate_manifest_job(row, project) for row in rows]
    if not args.execute:
        print(f"DRY RUN: validated {len(evidence)} opt/freq jobs; ORCA was not started.")
        return
    policy_path = project / "config" / "resource_policy.json"
    policy = read_json(policy_path)
    orca_exe = (args.orca_exe or Path(policy["machine"]["orca_executable"])).resolve()
    if not orca_exe.is_file():
        raise WorkflowValidationError(f"ORCA executable is missing: {orca_exe}")
    if sha256_file(orca_exe) != normalize_sha256(str(policy["machine"]["orca_sha256"])):
        raise WorkflowValidationError("ORCA executable hash differs from resource policy")
    state_root = (args.state_root or project / "logs" / "opt_freq_serial" / "jobs").resolve()
    run_root = (args.run_root or project / "02_dft_opt_freq_runs").resolve()
    with ExclusiveLock(project / "logs" / "opt_freq_serial" / "RUNNER.lock"):
        for row, item in zip(rows, evidence, strict=True):
            state, stop = run_one(row, item, project, state_root, run_root, policy, orca_exe)
            print(f"{row['job_id']}: {state['status']}")
            if state["status"] != "completed" or stop:
                break


if __name__ == "__main__":
    main()
