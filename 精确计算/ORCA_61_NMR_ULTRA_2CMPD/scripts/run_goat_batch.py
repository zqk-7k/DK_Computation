#!/usr/bin/env python3
"""Run the prepared ORCA GOAT jobs safely and recoverably.

The runner never launches ORCA in a source job directory.  It fingerprints and
copies ``goat.inp`` plus its local XYZ files into a unique attempt directory,
then executes the copy.  Per-job state and aggregate CSV/JSON files are written
atomically, so an interrupted batch can be restarted without repeating valid,
completed jobs.

The defaults deliberately favour calculation quality and recoverability over
throughput: one job at a time, the input's own parallel settings unchanged, and
a generous seven-day per-attempt timeout.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import hashlib
import io
import json
import os
import re
import shutil
import signal
import socket
import subprocess
import sys
import threading
import time
import uuid
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCRIPT = Path(__file__).resolve()
PROJECT = SCRIPT.parents[1]
DEFAULT_MANIFEST = PROJECT / "data" / "GOAT_explicit_job_manifest.csv"
DEFAULT_RUN_ROOT = PROJECT / "02_goat_explicit_runs"
DEFAULT_STATE_ROOT = PROJECT / "logs" / "goat_explicit_batch"
# HPC/package portability: command-line arguments and environment variables
# take precedence; these relative fallbacks are then resolved from PATH.
DEFAULT_ORCA = Path(os.environ.get("ORCA_EXE", "orca"))
DEFAULT_XTB_DIR = Path(os.environ.get("XTB_DIR", "xtb"))
NORMAL_MARKER = "ORCA TERMINATED NORMALLY"
SCHEMA_VERSION = 1
HOSTNAME = socket.gethostname()
STOP_AFTER_CURRENT = PROJECT / "logs" / "STOP_AFTER_CURRENT_STAGE"
GOAT_MINIMUM_FREE_RAM_GIB = 24.0
LEGACY_TWO_FILE_PILOT_JOB_ID = "F3_001_DMSO_01_primary_s01"
PROVENANCE_INPUT_NAMES = (
    "atom_map.csv",
    "isotope_map.csv",
    "structure_validation.json",
)
XYZFILE_LINE = re.compile(
    r"^\s*\*\s*xyzfile\s+[+-]?\d+\s+\d+\s+"
    r"(?:\"(?P<double>[^\"]+)\"|'(?P<single>[^']+)'|(?P<bare>\S+))",
    flags=re.IGNORECASE,
)

AGGREGATE_FIELDS = [
    "job_id",
    "figure_id",
    "queue_order",
    "compound_number",
    "display_name",
    "scenario_id",
    "primary_solvent",
    "explicit_solvent_count",
    "start_index",
    "charge",
    "multiplicity",
    "status",
    "attempt_count",
    "latest_attempt",
    "started_utc",
    "ended_utc",
    "elapsed_seconds",
    "return_code",
    "reason",
    "canonical_state_input_bundle_sha256",
    "latest_attempt_input_bundle_sha256",
    "source_job_directory",
    "latest_run_directory",
    "output_file",
    "globalminimum_file",
    "finalensemble_file",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def parse_utc(value: str | None) -> datetime | None:
    if not value:
        return None
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None


def safe_job_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    if not cleaned:
        raise ValueError(f"Unsafe or empty job id: {value!r}")
    return cleaned


class FILETIME(ctypes.Structure):
    _fields_ = [("dwLowDateTime", ctypes.c_uint32), ("dwHighDateTime", ctypes.c_uint32)]


class MEMORYSTATUSEX(ctypes.Structure):
    _fields_ = [
        ("dwLength", ctypes.c_uint32),
        ("dwMemoryLoad", ctypes.c_uint32),
        ("ullTotalPhys", ctypes.c_uint64),
        ("ullAvailPhys", ctypes.c_uint64),
        ("ullTotalPageFile", ctypes.c_uint64),
        ("ullAvailPageFile", ctypes.c_uint64),
        ("ullTotalVirtual", ctypes.c_uint64),
        ("ullAvailVirtual", ctypes.c_uint64),
        ("ullAvailExtendedVirtual", ctypes.c_uint64),
    ]


def _filetime_value(value: FILETIME) -> int:
    return (int(value.dwHighDateTime) << 32) | int(value.dwLowDateTime)


def system_cpu_times() -> tuple[int, int]:
    if os.name != "nt":
        values = os.times()
        return int(values.elapsed * 10_000_000), int((values.elapsed - values.system - values.user) * 10_000_000)
    idle, kernel, user = FILETIME(), FILETIME(), FILETIME()
    if not ctypes.windll.kernel32.GetSystemTimes(ctypes.byref(idle), ctypes.byref(kernel), ctypes.byref(user)):
        raise OSError("GetSystemTimes failed")
    return _filetime_value(kernel) + _filetime_value(user), _filetime_value(idle)


def cpu_percent_between(previous: tuple[int, int], current: tuple[int, int]) -> float:
    total_delta = current[0] - previous[0]
    idle_delta = current[1] - previous[1]
    if total_delta <= 0:
        return 0.0
    return max(0.0, min(100.0, 100.0 * (total_delta - idle_delta) / total_delta))


def sample_cpu(seconds: float, stop_event: threading.Event) -> float | None:
    interval = 2.0
    samples: list[float] = []
    previous = system_cpu_times()
    deadline = time.monotonic() + max(seconds, interval)
    while time.monotonic() < deadline:
        if stop_event.wait(min(interval, max(0.0, deadline - time.monotonic()))):
            return None
        current = system_cpu_times()
        samples.append(cpu_percent_between(previous, current))
        previous = current
    return round(sum(samples) / len(samples), 2) if samples else 0.0


def free_ram_gib() -> float:
    if os.name != "nt":
        return float("inf")
    status = MEMORYSTATUSEX()
    status.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
        raise OSError("GlobalMemoryStatusEx failed")
    return status.ullAvailPhys / (1024.0**3)


def external_compute_processes() -> list[str]:
    if os.name != "nt":
        return []
    completed = subprocess.run(
        ["tasklist.exe", "/FO", "CSV", "/NH"],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=30,
        check=False,
    )
    names: list[str] = []
    for raw in completed.stdout.splitlines():
        try:
            fields = next(csv.reader([raw]))
        except Exception:
            continue
        if fields and fields[0].casefold() in {"orca.exe", "mpiexec.exe", "xtb.exe"}:
            names.append(raw)
    return names


def atomic_write_text(path: Path, text: str, encoding: str = "utf-8") -> None:
    """Durably replace *path* with text without exposing a partial file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}.{uuid.uuid4().hex}")
    try:
        with temporary.open("w", encoding=encoding, newline="") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        try:
            temporary.unlink(missing_ok=True)
        except OSError:
            pass


def atomic_write_json(path: Path, value: Any) -> None:
    atomic_write_text(path, json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def load_json(path: Path) -> dict[str, Any] | None:
    try:
        with path.open(encoding="utf-8") as handle:
            value = json.load(handle)
        return value if isinstance(value, dict) else None
    except (OSError, json.JSONDecodeError):
        return None


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_input_files(job_directory: Path, input_name: str) -> list[Path]:
    """Return the immutable GOAT inputs that are staged for an attempt.

    Only coordinate files explicitly referenced by an ORCA
    ``* xyzfile charge multiplicity filename`` line are included.  Every
    reference must resolve to an existing regular file inside the source job
    directory.  This prevents prior/unrelated XYZ outputs from contaminating a
    new attempt and rejects path traversal (including traversal through a
    symlink or junction).
    """
    try:
        root = job_directory.resolve(strict=True)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Missing source job directory: {job_directory}") from exc
    if not root.is_dir():
        raise NotADirectoryError(f"Source job directory is not a directory: {root}")

    try:
        input_path = (root / input_name).resolve(strict=True)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Missing ORCA input: {root / input_name}") from exc
    try:
        input_path.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"ORCA input escapes source job directory: {input_path}") from exc
    if not input_path.is_file():
        raise FileNotFoundError(f"ORCA input is not a regular file: {input_path}")

    referenced: list[Path] = []
    for line_number, line in enumerate(
        input_path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1
    ):
        match = XYZFILE_LINE.match(line)
        if match is None:
            continue
        filename = next(value for value in match.group("double", "single", "bare") if value is not None)
        raw_path = Path(filename).expanduser()
        candidate = raw_path if raw_path.is_absolute() else root / raw_path
        try:
            resolved = candidate.resolve(strict=True)
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                f"ORCA xyzfile reference on line {line_number} does not exist: {filename!r} "
                f"(resolved from {input_path})"
            ) from exc
        try:
            resolved.relative_to(root)
        except ValueError as exc:
            raise ValueError(
                f"ORCA xyzfile reference on line {line_number} escapes source job directory: "
                f"{filename!r} -> {resolved}"
            ) from exc
        if not resolved.is_file():
            raise FileNotFoundError(
                f"ORCA xyzfile reference on line {line_number} is not a regular file: {resolved}"
            )
        referenced.append(resolved)

    provenance_files: list[Path] = []
    for name in PROVENANCE_INPUT_NAMES:
        candidate = (root / name).resolve(strict=True)
        try:
            candidate.relative_to(root)
        except ValueError as exc:
            raise ValueError(f"Provenance input escapes source job directory: {candidate}") from exc
        if not candidate.is_file():
            raise FileNotFoundError(f"Missing required provenance input: {candidate}")
        provenance_files.append(candidate)

    # This order intentionally matches prepare_explicit_goat_jobs.py:
    # goat.inp, start.xyz, atom_map.csv, isotope_map.csv,
    # structure_validation.json.  It is part of the locked bundle definition.
    files = [input_path, *referenced, *provenance_files]
    unique: dict[str, Path] = {}
    for path in files:
        relative = path.relative_to(root)
        key = os.path.normcase(str(relative))
        unique.setdefault(key, path)
    return list(unique.values())


def input_bundle_sha256(files: Sequence[Path]) -> str:
    if not files:
        raise ValueError("Cannot fingerprint an empty input bundle")
    root = files[0].resolve().parent
    digest = hashlib.sha256()
    for path in files:
        relative = path.resolve().relative_to(root).as_posix()
        digest.update(relative.encode("utf-8"))
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest()


def stage_inputs(files: Sequence[Path], work_directory: Path) -> None:
    if not files:
        raise ValueError("Cannot stage an empty input bundle")
    source_root = files[0].resolve().parent
    work_directory.mkdir(parents=True, exist_ok=False)
    for source in files:
        relative = source.resolve().relative_to(source_root)
        destination = work_directory / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def path_prepend(environment: dict[str, str], entries: Iterable[Path]) -> None:
    current = environment.get("PATH", "")
    parts: list[str] = []
    seen: set[str] = set()
    for entry in entries:
        text = str(entry)
        key = os.path.normcase(os.path.abspath(text))
        if entry.is_dir() and key not in seen:
            parts.append(text)
            seen.add(key)
    for text in current.split(os.pathsep):
        if not text:
            continue
        key = os.path.normcase(os.path.abspath(text))
        if key not in seen:
            parts.append(text)
            seen.add(key)
    environment["PATH"] = os.pathsep.join(parts)


def resolve_orca(value: str | None) -> Path:
    candidates: list[Path] = []
    if value:
        candidates.append(Path(value).expanduser())
    if os.environ.get("ORCA_EXE"):
        candidates.append(Path(os.environ["ORCA_EXE"]).expanduser())
    candidates.append(DEFAULT_ORCA)
    located = shutil.which("orca.exe") or shutil.which("orca")
    if located:
        candidates.append(Path(located))
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    attempted = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"Could not find ORCA executable; tried: {attempted}")


def resolve_xtb_dir(value: str | None) -> Path | None:
    candidates: list[Path] = []
    if value:
        candidates.append(Path(value).expanduser())
    if os.environ.get("XTB_DIR"):
        candidates.append(Path(os.environ["XTB_DIR"]).expanduser())
    candidates.append(DEFAULT_XTB_DIR)
    located = shutil.which("xtb.exe") or shutil.which("xtb")
    if located:
        candidates.append(Path(located).parent)
    for candidate in candidates:
        executable = candidate if candidate.is_file() else candidate / ("xtb.exe" if os.name == "nt" else "xtb")
        if executable.is_file():
            return executable.parent.resolve()
    return None


def resolve_msmpi_bin(value: str | None) -> Path | None:
    candidates: list[Path] = []
    if value:
        candidates.append(Path(value).expanduser())
    if os.environ.get("MSMPI_BIN"):
        candidates.append(Path(os.environ["MSMPI_BIN"]).expanduser())
    program_files = os.environ.get("ProgramFiles")
    if program_files:
        candidates.append(Path(program_files) / "Microsoft MPI" / "Bin")
    located = shutil.which("mpiexec.exe") or shutil.which("mpiexec")
    if located:
        candidates.append(Path(located).parent)
    for candidate in candidates:
        executable = candidate if candidate.is_file() else candidate / ("mpiexec.exe" if os.name == "nt" else "mpiexec")
        if executable.is_file():
            return executable.parent.resolve()
    return None


def parse_env_assignments(values: Sequence[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"--set-env requires KEY=VALUE, got {value!r}")
        key, assigned = value.split("=", 1)
        key = key.strip()
        if not key or "\x00" in key or "=" in key:
            raise ValueError(f"Invalid environment variable name: {key!r}")
        result[key] = assigned
    return result


def build_environment(
    orca_exe: Path,
    xtb_dir: Path | None,
    msmpi_bin: Path | None,
    omp_threads: int,
    additions: Mapping[str, str],
) -> dict[str, str]:
    environment = os.environ.copy()
    path_prepend(environment, [orca_exe.parent, *(p for p in (xtb_dir, msmpi_bin) if p)])
    # Avoid nested BLAS/OpenMP oversubscription around the parallelism already
    # requested in the unmodified ORCA input.
    environment["OMP_NUM_THREADS"] = str(omp_threads)
    environment["MKL_NUM_THREADS"] = str(omp_threads)
    environment["OPENBLAS_NUM_THREADS"] = str(omp_threads)
    environment.update(additions)
    if xtb_dir is not None:
        xtb_executable = (xtb_dir / ("xtb.exe" if os.name == "nt" else "xtb")).resolve()
        if not xtb_executable.is_file():
            raise FileNotFoundError(f"xTB executable not found: {xtb_executable}")
        # ORCA's external xTB interface consumes XTBEXE.  Set it after generic
        # additions so the validated executable selected by --xtb-dir/default
        # cannot be accidentally replaced by a stale inherited value.
        environment["XTBEXE"] = str(xtb_executable)
    return environment


def input_requests_parallel_orca(path: Path) -> bool:
    text = path.read_text(encoding="utf-8", errors="replace")
    values = [int(match) for match in re.findall(r"(?i)\bnprocs\s+(\d+)", text)]
    values.extend(int(match) for match in re.findall(r"(?i)\bPAL(\d+)\b", text))
    values.extend(int(match) for match in re.findall(r"(?i)\bNWorkers\s+(\d+)", text))
    return any(value > 1 for value in values)


def pid_is_alive(pid: int, hostname: str | None) -> bool:
    if pid <= 0 or (hostname and hostname.casefold() != HOSTNAME.casefold()):
        return False
    if pid == os.getpid():
        return True
    if os.name == "nt":
        process_query_limited_information = 0x1000
        still_active = 259
        handle = ctypes.windll.kernel32.OpenProcess(process_query_limited_information, False, pid)
        if not handle:
            return False
        try:
            exit_code = ctypes.c_ulong()
            if not ctypes.windll.kernel32.GetExitCodeProcess(handle, ctypes.byref(exit_code)):
                return False
            return exit_code.value == still_active
        finally:
            ctypes.windll.kernel32.CloseHandle(handle)
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


class ExclusiveFileLock:
    """Small cross-process lock based on atomic O_EXCL file creation."""

    def __init__(self, path: Path, stale_hours: float, wait_seconds: float = 0.0):
        self.path = path
        self.stale_hours = stale_hours
        self.wait_seconds = wait_seconds
        self.token = uuid.uuid4().hex
        self.acquired = False

    def _existing_is_stale(self) -> bool:
        record = load_json(self.path) or {}
        pid = int(record.get("pid") or -1)
        hostname = str(record.get("hostname") or "")
        if hostname.casefold() == HOSTNAME.casefold():
            # A live same-host owner is authoritative regardless of lock age.
            # Long production batches can legitimately outlive stale_hours.
            return not pid_is_alive(pid, hostname)
        created = parse_utc(str(record.get("created_utc") or ""))
        if created is None:
            try:
                age_seconds = time.time() - self.path.stat().st_mtime
            except OSError:
                return True
        else:
            age_seconds = (datetime.now(timezone.utc) - created).total_seconds()
        return age_seconds > self.stale_hours * 3600.0

    def acquire(self) -> bool:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        deadline = time.monotonic() + self.wait_seconds
        payload = json.dumps(
            {
                "token": self.token,
                "pid": os.getpid(),
                "hostname": HOSTNAME,
                "created_utc": utc_now(),
            },
            indent=2,
        ) + "\n"
        while True:
            try:
                descriptor = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
                    handle.write(payload)
                    handle.flush()
                    os.fsync(handle.fileno())
                self.acquired = True
                return True
            except FileExistsError:
                if self._existing_is_stale():
                    stale = self.path.with_name(f"{self.path.name}.stale.{int(time.time())}.{uuid.uuid4().hex[:8]}")
                    try:
                        os.replace(self.path, stale)
                    except OSError:
                        pass
                    continue
                if time.monotonic() >= deadline:
                    return False
                time.sleep(0.20)

    def release(self) -> None:
        if not self.acquired:
            return
        record = load_json(self.path) or {}
        if record.get("token") == self.token:
            try:
                self.path.unlink(missing_ok=True)
            except OSError:
                pass
        self.acquired = False

    def __enter__(self) -> "ExclusiveFileLock":
        if not self.acquire():
            raise TimeoutError(f"Could not acquire lock: {self.path}")
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        self.release()


@dataclass(frozen=True)
class Job:
    row: dict[str, str]
    job_id: str
    figure_id: str
    source_directory: Path


def load_manifest(path: Path) -> list[Job]:
    if not path.is_file():
        raise FileNotFoundError(
            f"Manifest not found: {path}. Run scripts/prepare_high_quality.py first."
        )
    jobs: list[Job] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row_number, row in enumerate(csv.DictReader(handle), start=2):
            job_id = (row.get("job_id") or "").strip()
            figure_id = (row.get("figure_id") or "").strip()
            directory = (row.get("job_directory") or "").strip()
            if not job_id or not figure_id or not directory:
                raise ValueError(f"Manifest row {row_number} lacks job_id, figure_id, or job_directory")
            if job_id in seen:
                raise ValueError(f"Duplicate job_id in manifest: {job_id}")
            seen.add(job_id)
            jobs.append(
                Job(
                    row=dict(row),
                    job_id=job_id,
                    figure_id=figure_id,
                    source_directory=Path(directory).expanduser().resolve(),
                )
            )
    if not jobs:
        raise ValueError(f"Manifest is empty: {path}")
    return jobs


def normalized_selectors(values: Sequence[str] | None) -> set[str]:
    selectors: set[str] = set()
    for raw in values or []:
        for token in raw.split(","):
            token = token.strip().upper()
            if not token:
                continue
            if token.isdigit():
                token = f"F3_{int(token):03d}"
            selectors.add(token)
    return selectors


def select_jobs(jobs: Sequence[Job], selectors: set[str]) -> list[Job]:
    if not selectors:
        return list(jobs)
    selected = [
        job for job in jobs
        if job.job_id.upper() in selectors or job.figure_id.upper() in selectors
    ]
    known = {job.job_id.upper() for job in jobs} | {job.figure_id.upper() for job in jobs}
    missing = sorted(selectors - known)
    if missing:
        raise ValueError(f"Unknown --only selector(s): {', '.join(missing)}")
    return selected


class StateStore:
    def __init__(self, state_root: Path, project: Path, jobs: Sequence[Job]):
        self.state_root = state_root
        self.jobs_root = state_root / "jobs"
        self.locks_root = state_root / "locks"
        self.aggregate_csv = project / "data" / "GOAT_batch_status.csv"
        self.aggregate_json = project / "data" / "GOAT_batch_status.json"
        self.jobs = list(jobs)
        self.thread_lock = threading.RLock()

    def state_path(self, job: Job) -> Path:
        return self.jobs_root / f"{safe_job_id(job.job_id)}.json"

    def job_lock_path(self, job: Job) -> Path:
        return self.locks_root / f"{safe_job_id(job.job_id)}.lock"

    def load(self, job: Job) -> dict[str, Any]:
        return load_json(self.state_path(job)) or {
            "schema_version": SCHEMA_VERSION,
            "job_id": job.job_id,
            "figure_id": job.figure_id,
            "source_job_directory": str(job.source_directory),
            "status": "pending",
            "attempt_count": 0,
            "attempts": [],
        }

    def save(self, job: Job, state: dict[str, Any]) -> None:
        state["schema_version"] = SCHEMA_VERSION
        state["job_id"] = job.job_id
        state["figure_id"] = job.figure_id
        state["source_job_directory"] = str(job.source_directory)
        with self.thread_lock:
            atomic_write_json(self.state_path(job), state)
            self.write_aggregate()

    def _row_for(self, job: Job) -> dict[str, Any]:
        state = self.load(job)
        attempts = state.get("attempts") if isinstance(state.get("attempts"), list) else []
        latest = attempts[-1] if attempts and isinstance(attempts[-1], dict) else {}
        artifacts = latest.get("artifacts") if isinstance(latest.get("artifacts"), dict) else {}
        return {
            "job_id": job.job_id,
            "figure_id": job.figure_id,
            "queue_order": job.row.get("queue_order", ""),
            "compound_number": job.row.get("compound_number", ""),
            "display_name": job.row.get("display_name", ""),
            "scenario_id": job.row.get("scenario_id", ""),
            "primary_solvent": job.row.get("primary_solvent", ""),
            "explicit_solvent_count": job.row.get("explicit_solvent_count", ""),
            "start_index": job.row.get("start_index", ""),
            "charge": job.row.get("charge", ""),
            "multiplicity": job.row.get("multiplicity", ""),
            "status": state.get("status", "pending"),
            "attempt_count": state.get("attempt_count", len(attempts)),
            "latest_attempt": latest.get("attempt_no", ""),
            "started_utc": latest.get("started_utc", ""),
            "ended_utc": latest.get("ended_utc", ""),
            "elapsed_seconds": latest.get("elapsed_seconds", ""),
            "return_code": latest.get("return_code", ""),
            "reason": latest.get("reason", state.get("reason", "")),
            "canonical_state_input_bundle_sha256": state.get("input_bundle_sha256", ""),
            "latest_attempt_input_bundle_sha256": latest.get("input_bundle_sha256", ""),
            "source_job_directory": str(job.source_directory),
            "latest_run_directory": latest.get("run_directory", ""),
            "output_file": latest.get("output_file", ""),
            "globalminimum_file": artifacts.get("globalminimum", ""),
            "finalensemble_file": artifacts.get("finalensemble", ""),
        }

    def write_aggregate(self) -> None:
        """Serialize aggregate writers across threads and runner processes."""
        with self.thread_lock:
            lock = ExclusiveFileLock(
                self.state_root / "aggregate.lock", stale_hours=1.0, wait_seconds=30.0
            )
            if not lock.acquire():
                raise TimeoutError("Timed out waiting for the aggregate status lock")
            try:
                rows = [self._row_for(job) for job in self.jobs]
                payload = {
                    "schema_version": SCHEMA_VERSION,
                    "updated_utc": utc_now(),
                    "manifest_jobs": len(self.jobs),
                    "jobs": rows,
                }
                atomic_write_json(self.aggregate_json, payload)
                stream = io.StringIO(newline="")
                writer = csv.DictWriter(stream, fieldnames=AGGREGATE_FIELDS, extrasaction="ignore")
                writer.writeheader()
                writer.writerows(rows)
                atomic_write_text(self.aggregate_csv, stream.getvalue(), encoding="utf-8-sig")
            finally:
                lock.release()


def read_tail(path: Path, maximum_bytes: int = 2 * 1024 * 1024) -> str:
    if not path.is_file():
        return ""
    with path.open("rb") as handle:
        size = handle.seek(0, os.SEEK_END)
        handle.seek(max(0, size - maximum_bytes))
        return handle.read().decode("utf-8", errors="replace")


def completed_artifacts_are_valid(state: Mapping[str, Any]) -> bool:
    if state.get("status") != "completed":
        return False
    attempts = state.get("attempts")
    if not isinstance(attempts, list) or not attempts or not isinstance(attempts[-1], dict):
        return False
    latest = attempts[-1]
    output = Path(str(latest.get("output_file") or ""))
    artifacts = latest.get("artifacts") if isinstance(latest.get("artifacts"), dict) else {}
    globalminimum = Path(str(artifacts.get("globalminimum") or ""))
    finalensemble = Path(str(artifacts.get("finalensemble") or ""))
    return (
        latest.get("status") == "completed"
        and latest.get("return_code") == 0
        and output.is_file()
        and NORMAL_MARKER in read_tail(output)
        and globalminimum.is_file()
        and globalminimum.stat().st_size > 0
        and finalensemble.is_file()
        and finalensemble.stat().st_size > 0
    )


def completed_state_is_valid(state: Mapping[str, Any], current_hash: str) -> bool:
    return (
        str(state.get("input_bundle_sha256") or "").casefold() == current_hash.casefold()
        and completed_artifacts_are_valid(state)
    )


def verified_legacy_pilot_migration(
    job: Job,
    state: Mapping[str, Any],
    current_files: Sequence[Path],
    current_hash: str,
    run_root: Path,
) -> tuple[dict[str, Any] | None, str]:
    """Verify the one known two-file pilot before adopting a five-file hash.

    This is deliberately not a general hash migration mechanism.  A completed
    result may be reused only for the named compatibility pilot, and only when
    the staged ``goat.inp`` and ``start.xyz`` are byte-identical (SHA-256) to
    the current, manifest-locked physics inputs.  Provenance maps do not affect
    the GOAT potential-energy surface, but their absence from the legacy work
    directory is recorded rather than hidden.
    """
    if job.job_id != LEGACY_TWO_FILE_PILOT_JOB_ID:
        return None, "legacy migration is restricted to the named compatibility pilot"
    if state.get("status") != "completed":
        return None, "state is not completed"
    attempts = state.get("attempts")
    if not isinstance(attempts, list) or not attempts or not isinstance(attempts[-1], dict):
        return None, "completed state lacks a latest attempt record"
    latest = attempts[-1]
    if latest.get("status") != "completed":
        return None, "latest attempt is not completed"

    source_root = job.source_directory.resolve()
    current_by_relative: dict[str, Path] = {}
    try:
        for path in current_files:
            relative = path.resolve(strict=True).relative_to(source_root).as_posix()
            current_by_relative[relative] = path.resolve(strict=True)
    except (OSError, ValueError) as exc:
        return None, f"current canonical input path validation failed: {exc}"
    physics_names = ("goat.inp", "start.xyz")
    if any(name not in current_by_relative for name in physics_names):
        return None, "current canonical bundle does not contain goat.inp and start.xyz"

    raw_run_directory = str(latest.get("run_directory") or "").strip()
    if not raw_run_directory:
        return None, "latest attempt lacks run_directory"
    try:
        work_directory = Path(raw_run_directory).resolve(strict=True)
        expected_job_root = (run_root.resolve() / safe_job_id(job.job_id)).resolve()
        work_directory.relative_to(expected_job_root)
    except (OSError, ValueError) as exc:
        return None, f"legacy work directory is missing or outside the job run root: {exc}"
    if not work_directory.is_dir():
        return None, "legacy run_directory is not a directory"

    verified_files: list[dict[str, str]] = []
    for name in physics_names:
        staged = work_directory / name
        if not staged.is_file():
            return None, f"legacy work directory lacks required staged physics input: {name}"
        source_sha = sha256_file(current_by_relative[name])
        staged_sha = sha256_file(staged)
        if source_sha.casefold() != staged_sha.casefold():
            return None, f"legacy staged {name} differs from the current manifest-locked source"
        verified_files.append(
            {
                "relative_path": name,
                "current_source_sha256": source_sha,
                "legacy_staged_sha256": staged_sha,
                "verification": "SHA256_MATCH",
            }
        )

    provenance_presence = {
        name: (work_directory / name).is_file() for name in PROVENANCE_INPUT_NAMES
    }
    return (
        {
            "migrated_utc": utc_now(),
            "migration_kind": "VERIFIED_LEGACY_TWO_FILE_COMPATIBILITY_PILOT",
            "reason": "expanded runtime bundle to include immutable atom/isotope/validation provenance",
            "legacy_runtime_bundle_sha256": state.get("input_bundle_sha256", ""),
            "legacy_attempt_bundle_sha256": latest.get("input_bundle_sha256", ""),
            "canonical_manifest_bundle_sha256": current_hash,
            "legacy_run_directory": str(work_directory),
            "verified_physics_inputs": verified_files,
            "legacy_staged_provenance_presence": provenance_presence,
            "orca_physics_inputs_unchanged": True,
            "verification_basis": "PER_FILE_SHA256_MATCH",
        },
        "verified legacy pilot physics inputs",
    )


def priority_class_name(value: int | None) -> str:
    if value is None:
        return "UNKNOWN"
    names = {
        0x00000040: "IDLE_PRIORITY_CLASS",
        0x00004000: "BELOW_NORMAL_PRIORITY_CLASS",
        0x00000020: "NORMAL_PRIORITY_CLASS",
        0x00008000: "ABOVE_NORMAL_PRIORITY_CLASS",
        0x00000080: "HIGH_PRIORITY_CLASS",
        0x00000100: "REALTIME_PRIORITY_CLASS",
    }
    return names.get(int(value), f"UNKNOWN_0x{int(value):08X}")


def ensure_windows_below_normal_priority(
    pid: int, kernel32: Any | None = None
) -> dict[str, Any]:
    """Best-effort verification/enforcement of the ORCA parent priority.

    Failure is non-destructive: the healthy calculation is left running and a
    warning is returned for the immutable attempt record.
    """
    target = int(getattr(subprocess, "BELOW_NORMAL_PRIORITY_CLASS", 0x00004000))
    result: dict[str, Any] = {
        "requested": priority_class_name(target),
        "requested_value": target,
        "before": "UNKNOWN",
        "before_value": None,
        "after": "UNKNOWN",
        "after_value": None,
        "set_attempted": False,
        "set_succeeded": False,
        "verified": False,
        "warning": "",
    }
    if os.name != "nt" and kernel32 is None:
        result["warning"] = "Windows priority verification is not applicable on this host"
        return result
    api = kernel32 if kernel32 is not None else ctypes.windll.kernel32
    process_set_information = 0x0200
    process_query_limited_information = 0x1000
    handle = api.OpenProcess(
        process_set_information | process_query_limited_information, False, int(pid)
    )
    if not handle:
        result["warning"] = "OpenProcess failed; ORCA priority could not be verified"
        return result
    try:
        before = int(api.GetPriorityClass(handle) or 0)
        result["before_value"] = before or None
        result["before"] = priority_class_name(before or None)
        if before != target:
            result["set_attempted"] = True
            result["set_succeeded"] = bool(api.SetPriorityClass(handle, target))
        after = int(api.GetPriorityClass(handle) or 0)
        result["after_value"] = after or None
        result["after"] = priority_class_name(after or None)
        result["verified"] = after == target
        if not result["verified"]:
            result["warning"] = (
                "ORCA parent process is not verified at BELOW_NORMAL priority; "
                "the calculation was left running"
            )
    except Exception as exc:
        result["warning"] = f"Priority verification failed: {type(exc).__name__}: {exc}"
    finally:
        try:
            api.CloseHandle(handle)
        except Exception:
            pass
    return result


def terminate_process_tree(process: subprocess.Popen[bytes]) -> None:
    if process.poll() is not None:
        return
    if os.name == "nt":
        try:
            subprocess.run(
                ["taskkill", "/PID", str(process.pid), "/T", "/F"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                timeout=30,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired):
            pass
    else:
        try:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=15)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
        except (OSError, ProcessLookupError):
            pass
    if process.poll() is None:
        try:
            process.kill()
        except OSError:
            pass


def artifact_paths(work_directory: Path, input_stem: str) -> dict[str, str]:
    def choose(pattern: str, preferred: str) -> Path | None:
        preferred_path = work_directory / preferred
        if preferred_path.is_file():
            return preferred_path
        candidates = sorted(work_directory.glob(pattern))
        return candidates[0] if candidates else None

    globalminimum = choose("*.globalminimum.xyz", f"{input_stem}.globalminimum.xyz")
    finalensemble = choose("*.finalensemble.xyz", f"{input_stem}.finalensemble.xyz")
    return {
        "globalminimum": str(globalminimum.resolve()) if globalminimum else "",
        "finalensemble": str(finalensemble.resolve()) if finalensemble else "",
    }


def classify_result(
    return_code: int | None,
    output_file: Path,
    artifacts: Mapping[str, str],
    timed_out: bool,
    interrupted: bool,
    source_unchanged: bool,
) -> tuple[str, str]:
    if interrupted:
        return "interrupted", "Runner interruption requested; process tree terminated"
    if timed_out:
        return "timeout", "Per-attempt wall-clock timeout exceeded"
    if not source_unchanged:
        return "failed", "Source input bundle changed while the attempt was running"
    if return_code != 0:
        return "failed", f"ORCA returned non-zero exit code {return_code}"
    if NORMAL_MARKER not in read_tail(output_file):
        return "failed", f"Output lacks normal termination marker: {NORMAL_MARKER}"
    missing = [name for name, value in artifacts.items() if not value or not Path(value).is_file()]
    if missing:
        return "failed", "Normal termination reported but required GOAT artifact(s) missing: " + ", ".join(missing)
    return "completed", "Return code, normal termination marker, and GOAT artifacts verified"


@dataclass
class RunnerConfig:
    orca_exe: Path
    environment: dict[str, str]
    input_name: str
    run_root: Path
    timeout_seconds: float | None
    max_attempts: int
    force: bool
    stale_lock_hours: float
    stop_event: threading.Event
    cpu_sample_seconds: float
    maximum_launch_cpu_percent: float
    minimum_free_ram_gib: float
    minimum_free_project_gib: float
    resource_wait_seconds: float
    runtime_high_cpu_percent: float
    runtime_high_cpu_seconds: float


def wait_for_resource_gate(config: RunnerConfig) -> tuple[bool, str]:
    """Wait until a new single-instance job can start without crowding the PC."""
    while not config.stop_event.is_set():
        if STOP_AFTER_CURRENT.exists():
            config.stop_event.set()
            return False, f"Stop-after-current flag exists: {STOP_AFTER_CURRENT}"
        external = external_compute_processes()
        if external:
            reason = "External ORCA/MPI/xTB process detected; waiting: " + "; ".join(external[:4])
            print(reason, flush=True)
            if config.stop_event.wait(config.resource_wait_seconds):
                return False, "Runner stop requested during external-process wait"
            continue
        average_cpu = sample_cpu(config.cpu_sample_seconds, config.stop_event)
        if average_cpu is None:
            return False, "Runner stop requested during CPU launch sample"
        ram = free_ram_gib()
        project_free = shutil.disk_usage(PROJECT).free / (1024.0**3)
        conditions = [
            average_cpu <= config.maximum_launch_cpu_percent,
            ram >= config.minimum_free_ram_gib,
            project_free >= config.minimum_free_project_gib,
        ]
        if all(conditions):
            return True, (
                f"Launch gate passed: CPU average {average_cpu:.1f}%, free RAM {ram:.1f} GiB, "
                f"project volume free {project_free:.1f} GiB"
            )
        reason = (
            f"Launch gate waiting: CPU {average_cpu:.1f}% (max {config.maximum_launch_cpu_percent:.1f}%), "
            f"RAM {ram:.1f} GiB (min {config.minimum_free_ram_gib:.1f}), "
            f"project free {project_free:.1f} GiB (min {config.minimum_free_project_gib:.1f})"
        )
        print(reason, flush=True)
        if config.stop_event.wait(config.resource_wait_seconds):
            return False, "Runner stop requested during resource wait"
    return False, "Runner stop requested"


def mark_interrupted_attempt(state: dict[str, Any]) -> None:
    attempts = state.get("attempts")
    if state.get("status") != "running" or not isinstance(attempts, list) or not attempts:
        return
    latest = attempts[-1]
    if not isinstance(latest, dict) or latest.get("status") != "running":
        return
    latest["status"] = "interrupted"
    latest["ended_utc"] = utc_now()
    latest["reason"] = "Recovered stale running state from an interrupted runner"
    state["status"] = "interrupted"
    state["reason"] = latest["reason"]


def run_job(job: Job, store: StateStore, config: RunnerConfig) -> str:
    lock = ExclusiveFileLock(
        store.job_lock_path(job), stale_hours=config.stale_lock_hours, wait_seconds=0.0
    )
    if not lock.acquire():
        return "locked"
    try:
        files = source_input_files(job.source_directory, config.input_name)
        source_hash_before = input_bundle_sha256(files)
        manifest_hash = str(job.row.get("input_bundle_sha256") or "").strip()
        if not manifest_hash:
            raise ValueError(f"Manifest lacks input_bundle_sha256 for {job.job_id}")
        if manifest_hash.casefold() != source_hash_before.casefold():
            raise ValueError(
                f"Canonical input bundle hash mismatch for {job.job_id}: "
                f"manifest={manifest_hash}, current={source_hash_before}"
            )
        state = store.load(job)
        # The first compatibility pilot was launched before the runner began
        # staging the three provenance maps.  Reuse is allowed only for that
        # named pilot and only after its staged goat.inp/start.xyz have been
        # independently verified against the current manifest-locked sources.
        # The original attempt hash remains untouched as an audit trail.
        if (
            completed_artifacts_are_valid(state)
            and str(state.get("input_bundle_sha256") or "").casefold()
            != source_hash_before.casefold()
        ):
            migration, migration_reason = verified_legacy_pilot_migration(
                job=job,
                state=state,
                current_files=files,
                current_hash=source_hash_before,
                run_root=config.run_root,
            )
            if migration is not None:
                migrations = state.setdefault("input_bundle_hash_migrations", [])
                if not isinstance(migrations, list):
                    migrations = []
                    state["input_bundle_hash_migrations"] = migrations
                migrations.append(migration)
                state["input_bundle_sha256"] = source_hash_before
                store.save(job, state)
            else:
                print(
                    f"[{job.job_id}] completed-state hash migration rejected: "
                    f"{migration_reason}; a new attempt is required",
                    file=sys.stderr,
                    flush=True,
                )
        if completed_state_is_valid(state, source_hash_before) and not config.force:
            return "skipped_completed"
        mark_interrupted_attempt(state)
        attempt_count = int(state.get("attempt_count") or 0)
        if attempt_count >= config.max_attempts and not config.force:
            state["status"] = "exhausted"
            state["reason"] = f"Maximum {config.max_attempts} attempts reached; use --force to retry"
            state["input_bundle_sha256"] = source_hash_before
            store.save(job, state)
            return "exhausted"
        if config.stop_event.is_set():
            return "interrupted"

        launch_allowed, launch_reason = wait_for_resource_gate(config)
        if not launch_allowed:
            return "interrupted"
        print(f"{job.job_id}: {launch_reason}", flush=True)

        attempt_no = attempt_count + 1
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        attempt_directory = (
            config.run_root
            / safe_job_id(job.job_id)
            / f"attempt_{attempt_no:03d}_{stamp}_{uuid.uuid4().hex[:8]}"
        )
        work_directory = attempt_directory / "work"
        stage_inputs(files, work_directory)
        copied_input = work_directory / config.input_name
        output_file = work_directory / f"{copied_input.stem}.out"
        error_file = work_directory / f"{copied_input.stem}.stderr.log"
        command = [str(config.orca_exe), copied_input.name]
        attempt: dict[str, Any] = {
            "attempt_no": attempt_no,
            "status": "running",
            "started_utc": utc_now(),
            "ended_utc": "",
            "elapsed_seconds": "",
            "return_code": "",
            "reason": "",
            "input_bundle_sha256": source_hash_before,
            "run_directory": str(work_directory.resolve()),
            "output_file": str(output_file.resolve()),
            "stderr_file": str(error_file.resolve()),
            "command": command,
            "orca_executable": str(config.orca_exe),
            "runner_pid": os.getpid(),
            "hostname": HOSTNAME,
            "artifacts": {},
        }
        attempts = state.setdefault("attempts", [])
        if not isinstance(attempts, list):
            attempts = []
            state["attempts"] = attempts
        attempts.append(attempt)
        state["attempt_count"] = attempt_no
        state["status"] = "running"
        state["reason"] = ""
        state["input_bundle_sha256"] = source_hash_before
        store.save(job, state)

        start = time.monotonic()
        timed_out = False
        interrupted = False
        return_code: int | None = None
        high_cpu_seconds = 0.0
        previous_cpu_times = system_cpu_times()
        next_cpu_check = time.monotonic() + 10.0
        creationflags = (
            subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.BELOW_NORMAL_PRIORITY_CLASS
            if os.name == "nt"
            else 0
        )
        popen_kwargs: dict[str, Any] = {}
        if os.name != "nt":
            popen_kwargs["start_new_session"] = True
        try:
            with output_file.open("wb") as stdout_handle, error_file.open("wb") as stderr_handle:
                process = subprocess.Popen(
                    command,
                    cwd=work_directory,
                    env=config.environment,
                    stdout=stdout_handle,
                    stderr=stderr_handle,
                    stdin=subprocess.DEVNULL,
                    shell=False,
                    creationflags=creationflags,
                    **popen_kwargs,
                )
                attempt["orca_pid"] = process.pid
                priority_record = ensure_windows_below_normal_priority(process.pid)
                attempt["orca_parent_process_priority"] = priority_record
                if priority_record.get("warning"):
                    attempt["priority_warning"] = priority_record["warning"]
                    print(
                        f"[{job.job_id}] priority warning: {priority_record['warning']}",
                        file=sys.stderr,
                        flush=True,
                    )
                store.save(job, state)
                while process.poll() is None:
                    if config.stop_event.wait(1.0):
                        interrupted = True
                        terminate_process_tree(process)
                        break
                    if (
                        config.timeout_seconds is not None
                        and time.monotonic() - start > config.timeout_seconds
                    ):
                        timed_out = True
                        terminate_process_tree(process)
                        break
                    if time.monotonic() >= next_cpu_check:
                        current_cpu_times = system_cpu_times()
                        current_cpu = cpu_percent_between(previous_cpu_times, current_cpu_times)
                        previous_cpu_times = current_cpu_times
                        next_cpu_check = time.monotonic() + 10.0
                        if current_cpu >= config.runtime_high_cpu_percent:
                            high_cpu_seconds += 10.0
                        else:
                            high_cpu_seconds = 0.0
                        attempt["latest_system_cpu_percent"] = round(current_cpu, 2)
                        attempt["high_cpu_sustained_seconds"] = high_cpu_seconds
                        if (
                            high_cpu_seconds >= config.runtime_high_cpu_seconds
                            and not STOP_AFTER_CURRENT.exists()
                        ):
                            atomic_write_text(
                                STOP_AFTER_CURRENT,
                                "High CPU sustained; finish current ORCA job, then stop.\n",
                            )
                            attempt["stop_after_current_reason"] = "sustained high system CPU"
                        store.save(job, state)
                try:
                    return_code = process.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    terminate_process_tree(process)
                    return_code = process.poll()
        except OSError as exc:
            attempt["launch_error"] = f"{type(exc).__name__}: {exc}"

        elapsed = round(time.monotonic() - start, 3)
        artifacts = artifact_paths(work_directory, copied_input.stem)
        try:
            source_hash_after = input_bundle_sha256(source_input_files(job.source_directory, config.input_name))
        except (OSError, FileNotFoundError):
            source_hash_after = ""
        status, reason = classify_result(
            return_code=return_code,
            output_file=output_file,
            artifacts=artifacts,
            timed_out=timed_out,
            interrupted=interrupted,
            source_unchanged=(source_hash_after == source_hash_before),
        )
        attempt.update(
            {
                "status": status,
                "ended_utc": utc_now(),
                "elapsed_seconds": elapsed,
                "return_code": return_code,
                "reason": reason,
                "source_input_unchanged": source_hash_after == source_hash_before,
                "source_input_bundle_sha256_after": source_hash_after,
                "output_sha256": sha256_file(output_file) if output_file.is_file() else "",
                "artifacts": artifacts,
            }
        )
        state["status"] = status
        state["reason"] = reason
        state["input_bundle_sha256"] = source_hash_before
        store.save(job, state)
        if STOP_AFTER_CURRENT.exists():
            # The current healthy job has finished.  Prevent every queued future
            # from starting another ORCA process; do not retroactively kill this job.
            config.stop_event.set()
        return status
    except Exception as exc:
        state = store.load(job)
        state["status"] = "failed"
        state["reason"] = f"Runner exception: {type(exc).__name__}: {exc}"
        store.save(job, state)
        return "failed"
    finally:
        lock.release()


def mark_environment_block(
    jobs: Sequence[Job], store: StateStore, reason: str, input_name: str
) -> None:
    for job in jobs:
        state = store.load(job)
        if state.get("status") == "completed":
            continue
        state["status"] = "blocked_environment"
        state["reason"] = reason
        try:
            files = source_input_files(job.source_directory, input_name)
            state["input_bundle_sha256"] = input_bundle_sha256(files)
        except OSError:
            pass
        store.save(job, state)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run prepared ORCA GOAT jobs in isolated attempt directories with "
            "atomic status files and restart-safe completed-job skipping."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--only",
        nargs="+",
        metavar="ID",
        help="Run only Figure IDs or job IDs; accepts spaces/commas and bare numbers (for example F3_001 or 16)",
    )
    parser.add_argument("-j", "--concurrency", type=int, default=1, help="Number of ORCA jobs launched concurrently")
    parser.add_argument("--timeout-hours", type=float, default=168.0, help="Wall-clock limit per attempt; 0 disables the timeout")
    parser.add_argument("--max-attempts", type=int, default=3, help="Maximum recorded attempts per job before requiring --force")
    parser.add_argument("--force", action="store_true", help="Rerun even a valid completed or attempt-exhausted job")
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Launch ORCA after validation; without this flag the runner is dry-run only",
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST, help="Prepared GOAT CSV manifest")
    parser.add_argument("--input-name", default="goat.inp", help="ORCA input filename inside each manifest job directory")
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT, help="Root for immutable per-attempt working directories")
    parser.add_argument("--state-root", type=Path, default=DEFAULT_STATE_ROOT, help="Root for per-job JSON and lock files")
    parser.add_argument("--orca-exe", help="ORCA executable; then ORCA_EXE and PATH are tried")
    parser.add_argument("--xtb-dir", help="Directory containing xtb.exe; prepended to PATH (then XTB_DIR/default/PATH are tried)")
    parser.add_argument("--msmpi-bin", help="Directory containing mpiexec.exe; prepended to PATH (then MSMPI_BIN/default/PATH are tried)")
    parser.add_argument("--allow-missing-msmpi", action="store_true", help="Explicitly allow parallel ORCA inputs when mpiexec is not discoverable")
    parser.add_argument("--omp-threads", type=int, default=1, help="OMP/MKL/OpenBLAS threads per ORCA process")
    parser.add_argument("--set-env", action="append", default=[], metavar="KEY=VALUE", help="Additional child-process environment assignment; repeat as needed")
    parser.add_argument("--stale-lock-hours", type=float, default=192.0, help="Age after which a foreign-host lock may be recovered")
    parser.add_argument("--cpu-sample-seconds", type=float, default=30.0, help="System-CPU averaging window before every ORCA launch")
    parser.add_argument("--maximum-launch-cpu-percent", type=float, default=45.0, help="Do not launch a new ORCA job above this average system load")
    parser.add_argument("--minimum-free-ram-gib", type=float, default=GOAT_MINIMUM_FREE_RAM_GIB)
    parser.add_argument("--minimum-free-project-gib", type=float, default=100.0)
    parser.add_argument("--resource-wait-seconds", type=float, default=60.0)
    parser.add_argument("--runtime-high-cpu-percent", type=float, default=88.0, help="Soft stop-after-current threshold")
    parser.add_argument("--runtime-high-cpu-seconds", type=float, default=90.0, help="Sustained duration for the soft stop")
    parser.add_argument("--clear-stop-after-current", action="store_true", help="Clear a prior soft-stop flag before resuming")
    return parser


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if args.concurrency != 1:
        parser.error("This high-quality workstation runner requires --concurrency 1")
    if args.timeout_hours < 0:
        parser.error("--timeout-hours cannot be negative")
    if args.max_attempts < 1:
        parser.error("--max-attempts must be at least 1")
    if args.omp_threads < 1:
        parser.error("--omp-threads must be at least 1")
    if args.stale_lock_hours <= 0:
        parser.error("--stale-lock-hours must be positive")
    if not args.input_name or Path(args.input_name).name != args.input_name:
        parser.error("--input-name must be a filename, not a path")
    if args.cpu_sample_seconds < 2 or args.resource_wait_seconds < 1:
        parser.error("CPU sampling must be at least 2 s and resource wait at least 1 s")
    for name in (
        "maximum_launch_cpu_percent",
        "runtime_high_cpu_percent",
    ):
        value = float(getattr(args, name))
        if value <= 0 or value >= 100:
            parser.error(f"--{name.replace('_', '-')} must be between 0 and 100")
    if args.runtime_high_cpu_seconds < 10:
        parser.error("--runtime-high-cpu-seconds must be at least 10")


def effective_minimum_free_ram_gib(requested: float) -> float:
    """Apply the non-weakenable workstation RAM floor."""
    return max(float(requested), GOAT_MINIMUM_FREE_RAM_GIB)


def execute_runner(args: argparse.Namespace) -> int:
    if args.clear_stop_after_current:
        STOP_AFTER_CURRENT.unlink(missing_ok=True)
    try:
        additions = parse_env_assignments(args.set_env)
        jobs = load_manifest(args.manifest.expanduser().resolve())
        selected = select_jobs(jobs, normalized_selectors(args.only))
        if not selected:
            raise ValueError("No jobs selected")
        orca_exe = resolve_orca(args.orca_exe)
        xtb_dir = resolve_xtb_dir(args.xtb_dir)
        msmpi_bin = resolve_msmpi_bin(args.msmpi_bin)
        environment = build_environment(
            orca_exe=orca_exe,
            xtb_dir=xtb_dir,
            msmpi_bin=msmpi_bin,
            omp_threads=args.omp_threads,
            additions=additions,
        )
        parallel_inputs = [
            job for job in selected
            if input_requests_parallel_orca(job.source_directory / args.input_name)
        ]
    except (OSError, ValueError) as exc:
        print(f"Configuration error: {exc}", file=sys.stderr)
        return 2

    print(f"Project:       {PROJECT}")
    print(f"Manifest:      {args.manifest.expanduser().resolve()}")
    print(f"Selected jobs: {len(selected)} / {len(jobs)}")
    print(f"Concurrency:   {args.concurrency}")
    print(f"ORCA:          {orca_exe}")
    print(f"xTB directory: {xtb_dir or 'not found (ORCA may use its native XTB implementation)'}")
    print(f"MS-MPI bin:    {msmpi_bin or 'not found'}")
    print(f"Run root:      {args.run_root.expanduser().resolve()}")
    print(f"State root:    {args.state_root.expanduser().resolve()}")

    if not args.execute:
        for job in selected:
            print(f"  {job.job_id:<14} {job.figure_id:<8} {job.source_directory}")
        if parallel_inputs and msmpi_bin is None and not args.allow_missing_msmpi:
            print(
                "Dry-run warning: selected input(s) request parallel workers but mpiexec was not found.",
                file=sys.stderr,
            )
        return 0

    store = StateStore(args.state_root.expanduser().resolve(), PROJECT, jobs)
    if parallel_inputs and msmpi_bin is None and not args.allow_missing_msmpi:
        reason = (
            "Selected ORCA input requests multiple processes/workers, but mpiexec was not found. "
            "Install/configure Microsoft MPI or pass --msmpi-bin; use --allow-missing-msmpi only "
            "after independently verifying the runtime."
        )
        mark_environment_block(selected, store, reason, args.input_name)
        print(f"Blocked: {reason}", file=sys.stderr)
        return 2

    store.write_aggregate()
    stop_event = threading.Event()
    config = RunnerConfig(
        orca_exe=orca_exe,
        environment=environment,
        input_name=args.input_name,
        run_root=args.run_root.expanduser().resolve(),
        timeout_seconds=None if args.timeout_hours == 0 else args.timeout_hours * 3600.0,
        max_attempts=args.max_attempts,
        force=args.force,
        stale_lock_hours=args.stale_lock_hours,
        stop_event=stop_event,
        cpu_sample_seconds=args.cpu_sample_seconds,
        maximum_launch_cpu_percent=args.maximum_launch_cpu_percent,
        minimum_free_ram_gib=args.minimum_free_ram_gib,
        minimum_free_project_gib=args.minimum_free_project_gib,
        resource_wait_seconds=args.resource_wait_seconds,
        runtime_high_cpu_percent=args.runtime_high_cpu_percent,
        runtime_high_cpu_seconds=args.runtime_high_cpu_seconds,
    )

    counts: dict[str, int] = {}
    futures: dict[Future[str], Job] = {}
    interrupted_by_user = False
    with ThreadPoolExecutor(max_workers=args.concurrency, thread_name_prefix="goat") as executor:
        try:
            for job in selected:
                futures[executor.submit(run_job, job, store, config)] = job
            for future in as_completed(futures):
                job = futures[future]
                try:
                    status = future.result()
                except Exception as exc:  # defensive: run_job normally records its own failures
                    status = "failed"
                    print(f"[{job.job_id}] runner exception: {type(exc).__name__}: {exc}", file=sys.stderr)
                counts[status] = counts.get(status, 0) + 1
                print(f"[{job.job_id}] {status}", flush=True)
        except KeyboardInterrupt:
            # Catch inside the executor context: setting the event here lets
            # worker loops kill their process trees before __exit__ waits.
            print("Interrupt received; terminating active ORCA process trees safely...", file=sys.stderr)
            interrupted_by_user = True
            stop_event.set()
            for future in futures:
                future.cancel()
    stop_event.set()
    store.write_aggregate()

    if interrupted_by_user:
        return 130

    summary = ", ".join(f"{key}={counts[key]}" for key in sorted(counts))
    print(f"Batch summary: {summary or 'no completed futures'}")
    print(f"Aggregate CSV:  {store.aggregate_csv}")
    print(f"Aggregate JSON: {store.aggregate_json}")
    acceptable = {"completed", "skipped_completed"}
    return 0 if counts and set(counts).issubset(acceptable) else 1


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)
    # This floor is part of the machine-specific locked GOAT policy.  A caller
    # may request a stricter threshold, but cannot silently weaken it.
    args.minimum_free_ram_gib = effective_minimum_free_ram_gib(args.minimum_free_ram_gib)
    if not args.execute:
        return execute_runner(args)

    state_root = args.state_root.expanduser().resolve()
    runner_lock = ExclusiveFileLock(
        state_root / "runner_global.lock",
        stale_hours=args.stale_lock_hours,
        wait_seconds=0.0,
    )
    if not runner_lock.acquire():
        print(
            f"Blocked: another Python GOAT runner owns {runner_lock.path}",
            file=sys.stderr,
        )
        return 3
    try:
        return execute_runner(args)
    finally:
        runner_lock.release()


if __name__ == "__main__":
    raise SystemExit(main())
