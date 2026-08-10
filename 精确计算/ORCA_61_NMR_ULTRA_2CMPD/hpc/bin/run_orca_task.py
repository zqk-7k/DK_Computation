#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import uuid
from pathlib import Path
from typing import Any, Mapping

from hpc_common import (
    NORMAL_MARKER,
    approval_problem,
    atomic_json,
    canonical_row_hash,
    dependency_contract,
    enabled_row,
    glob_nonempty,
    input_coordinate_references,
    load_json,
    parse_bool,
    parse_orca_resources,
    read_manifest,
    required_artifact_patterns,
    resolve_manifest_path,
    resolve_package_path,
    resolve_required_helpers,
    row_value,
    safe_component,
    sha256_file,
    utc_now,
    verify_materialized_dependency_receipt,
    write_csv,
)


INPUT_HASH_FIELDS = (
    "expected_input_sha256",
    "input_file_sha256",
    "input_sha256",
    "goat_input_sha256",
    "opt_input_sha256",
    "shielding_input_sha256",
    "j_input_sha256",
    "coupling_input_sha256",
)
COORDINATE_HASH_FIELDS = (
    "expected_coordinate_sha256",
    "coordinate_sha256",
    "start_xyz_sha256",
    "selected_xyz_sha256",
    "optimized_xyz_sha256",
    "conformer_xyz_sha256",
)
LEGACY_GOAT_BUNDLE_FILES = (
    "goat.inp",
    "start.xyz",
    "atom_map.csv",
    "isotope_map.csv",
    "structure_validation.json",
)


class TaskFailure(RuntimeError):
    pass


def scheduler_context() -> dict[str, str]:
    return {
        "scheduler": "slurm" if os.environ.get("SLURM_JOB_ID") else (
            "pbs" if os.environ.get("PBS_JOBID") else "direct"
        ),
        "job_id": os.environ.get("SLURM_JOB_ID") or os.environ.get("PBS_JOBID", ""),
        "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID") or os.environ.get("PBS_ARRAY_ID", ""),
        "array_index": os.environ.get("SLURM_ARRAY_TASK_ID") or os.environ.get("PBS_ARRAY_INDEX", ""),
        "hostname": os.environ.get("HOSTNAME", ""),
    }


def allocated_cpus(explicit: int | None) -> int:
    if explicit is not None:
        return explicit
    for name in ("SLURM_CPUS_PER_TASK", "PBS_NCPUS", "NCPUS"):
        value = os.environ.get(name, "").strip()
        if value:
            return int(value)
    return 1


def job_directory(package_root: Path, row: Mapping[str, str]) -> Path:
    raw = row_value(
        row,
        "portable_job_directory",
        "relative_job_directory",
        "working_directory",
        "job_directory",
    )
    return resolve_package_path(package_root, raw, field="job_directory")


def input_path_for(
    package_root: Path,
    job_dir: Path,
    row: Mapping[str, str],
    stage_config: Mapping[str, Any],
) -> Path:
    raw = row_value(row, "portable_input_file", "input_file") or str(stage_config["input_default"])
    token = raw.replace("\\", "/")
    if len(token) >= 3 and token[1:3] == ":/":
        token = Path(token).name
    candidate = Path(token)
    if not candidate.is_absolute():
        direct = job_dir / candidate
        if direct.is_file():
            return direct.resolve()
        package_candidate = package_root / candidate
        if package_candidate.is_file():
            return package_candidate.resolve()
    raise FileNotFoundError(f"Input file is missing for {row_value(row, 'job_id')}: {raw}")


def output_name_for(row: Mapping[str, str], input_path: Path) -> str:
    raw = row_value(row, "output_file", "orca_output_file")
    if raw:
        name = Path(raw.replace("\\", "/")).name
    else:
        name = input_path.with_suffix(".out").name
    if not name:
        raise ValueError("Empty ORCA output name")
    return name


def verify_materialized_dependency(
    package_root: Path, row: Mapping[str, str], job_dir: Path, stage: str
) -> None:
    try:
        parent_id, parent_stage, _ = dependency_contract(row, stage)
        if parent_id and parent_stage != stage:
            verify_materialized_dependency_receipt(package_root, row, job_dir, stage)
    except (OSError, ValueError) as exc:
        raise TaskFailure(str(exc)) from exc


def legacy_goat_bundle_hash(job_dir: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    for name in LEGACY_GOAT_BUNDLE_FILES:
        path = job_dir / name
        if not path.is_file():
            raise FileNotFoundError(f"Legacy GOAT bundle file missing: {path}")
        digest.update(name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest().upper()


def verify_hashes(row: Mapping[str, str], job_dir: Path, input_path: Path) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    expected_input = row_value(row, *INPUT_HASH_FIELDS).upper()
    actual_input = sha256_file(input_path)
    if expected_input:
        checks.append(
            {
                "kind": "input",
                "path": str(input_path),
                "expected": expected_input,
                "actual": actual_input,
                "matches": expected_input == actual_input,
            }
        )

    coordinate_paths = input_coordinate_references(input_path)
    expected_coordinate = row_value(row, *COORDINATE_HASH_FIELDS).upper()
    if expected_coordinate:
        if len(coordinate_paths) != 1:
            checks.append(
                {
                    "kind": "coordinate",
                    "expected": expected_coordinate,
                    "actual": "",
                    "matches": False,
                    "reason": f"Expected exactly one xyzfile/gzmtfile reference, found {len(coordinate_paths)}",
                }
            )
        else:
            coordinate = coordinate_paths[0]
            if not coordinate.is_file():
                raise FileNotFoundError(f"Referenced coordinate file is missing: {coordinate}")
            actual_coordinate = sha256_file(coordinate)
            checks.append(
                {
                    "kind": "coordinate",
                    "path": str(coordinate),
                    "expected": expected_coordinate,
                    "actual": actual_coordinate,
                    "matches": expected_coordinate == actual_coordinate,
                }
            )

    expected_bundle = row_value(row, "expected_input_bundle_sha256", "input_bundle_sha256").upper()
    recorded_bundle_path = job_dir / "INPUT_BUNDLE_SHA256.txt"
    if expected_bundle:
        recorded = recorded_bundle_path.read_text(encoding="utf-8").strip().upper() if recorded_bundle_path.is_file() else ""
        checks.append(
            {
                "kind": "recorded_bundle",
                "path": str(recorded_bundle_path),
                "expected": expected_bundle,
                "actual": recorded,
                "matches": bool(recorded) and expected_bundle == recorded,
            }
        )
        if all((job_dir / name).is_file() for name in LEGACY_GOAT_BUNDLE_FILES):
            actual_bundle = legacy_goat_bundle_hash(job_dir)
            checks.append(
                {
                    "kind": "recomputed_legacy_goat_bundle",
                    "expected": expected_bundle,
                    "actual": actual_bundle,
                    "matches": expected_bundle == actual_bundle,
                }
            )

    return {
        "strict_checks_present": bool(checks),
        "matches": all(bool(check["matches"]) for check in checks),
        "checks": checks,
        "actual_input_sha256": actual_input,
    }


def copy_source_to_scratch(source: Path, target: Path) -> None:
    def ignore(directory: str, names: list[str]) -> set[str]:
        ignored: set[str] = set()
        for name in names:
            lower = name.casefold()
            if lower.endswith((".out", ".stderr.log")) or lower.startswith("attempt_"):
                ignored.add(name)
            if lower in {"worker_tmp", "scratch", "tmp"}:
                ignored.add(name)
        return ignored

    shutil.copytree(source, target, ignore=ignore)


def next_attempt_directory(job_result_root: Path) -> tuple[int, Path, Path]:
    existing: list[int] = []
    if job_result_root.is_dir():
        for child in job_result_root.iterdir():
            if child.is_dir() and child.name.startswith("attempt_") and ".partial." not in child.name:
                try:
                    existing.append(int(child.name.split("_", 1)[1]))
                except ValueError:
                    pass
    number = max(existing, default=0) + 1
    final = job_result_root / f"attempt_{number:03d}"
    partial = job_result_root / f"attempt_{number:03d}.partial.{uuid.uuid4().hex[:10]}"
    return number, partial, final


def copy_back(
    scratch_work: Path,
    partial: Path,
    final: Path,
    stage_config: Mapping[str, Any],
    metadata: dict[str, Any],
    required_patterns: tuple[str, ...] = (),
) -> list[dict[str, Any]]:
    partial.mkdir(parents=True, exist_ok=False)
    copied: set[Path] = set()
    copy_patterns = list(stage_config.get("copy_back_globs", ["*.out", "*.stderr.log"]))
    for pattern in required_patterns:
        if pattern not in copy_patterns:
            copy_patterns.append(pattern)
    for pattern in copy_patterns:
        for source in scratch_work.glob(str(pattern)):
            if not source.is_file() or source in copied:
                continue
            shutil.copy2(source, partial / source.name)
            copied.add(source)
    atomic_json(partial / "run_metadata.json", metadata)
    records: list[dict[str, Any]] = []
    for path in sorted(item for item in partial.iterdir() if item.is_file()):
        records.append(
            {
                "file": path.name,
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    write_csv(partial / "FILE_MANIFEST_SHA256.csv", ("file", "size_bytes", "sha256"), records)
    os.replace(partial, final)
    return records


def completed_state_is_valid(
    state: Mapping[str, Any],
    expected_input_hash: str,
    package_root: Path,
    required_patterns: tuple[str, ...],
    orca_executable_sha256: str,
    required_helpers: list[dict[str, str]],
) -> bool:
    if state.get("status") not in {"COMPLETED", "RESUMED_ALREADY_COMPLETE"}:
        return False
    if str(state.get("input_sha256", "")).upper() != expected_input_hash.upper():
        return False
    if str(state.get("orca_executable_sha256", "")).upper() != orca_executable_sha256.upper():
        return False
    if state.get("required_helpers") != required_helpers:
        return False
    attempt_path = Path(str(state.get("attempt_path", "")))
    if not attempt_path.is_absolute():
        attempt_path = package_root / attempt_path
    output_name = str(state.get("output_file", ""))
    output = attempt_path / output_name
    if not output.is_file() or NORMAL_MARKER not in output.read_text(encoding="utf-8", errors="replace"):
        return False
    ledger_path = attempt_path / "FILE_MANIFEST_SHA256.csv"
    expected_ledger_hash = str(state.get("file_manifest_sha256", "")).strip().upper()
    if not ledger_path.is_file() or not expected_ledger_hash:
        return False
    if sha256_file(ledger_path) != expected_ledger_hash:
        return False
    try:
        with ledger_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or not {"file", "size_bytes", "sha256"}.issubset(reader.fieldnames):
                return False
            ledger_rows = list(reader)
        ledger: dict[str, str] = {}
        for row in ledger_rows:
            name = str(row.get("file", "")).strip()
            token = Path(name)
            if not name or token.is_absolute() or token.name != name or name in ledger:
                return False
            artifact = attempt_path / name
            if not artifact.is_file():
                return False
            if artifact.stat().st_size != int(str(row.get("size_bytes", ""))):
                return False
            digest = str(row.get("sha256", "")).strip().upper()
            if not digest or sha256_file(artifact) != digest:
                return False
            ledger[name] = digest
        if output_name not in ledger:
            return False
        for pattern in required_patterns:
            matches = glob_nonempty(attempt_path, str(pattern))
            if not matches or any(path.name not in ledger for path in matches):
                return False
    except (OSError, TypeError, ValueError):
        return False
    return True


def terminate_process(process: subprocess.Popen[bytes], grace_seconds: float = 30.0) -> None:
    if process.poll() is not None:
        return
    try:
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGTERM)
        else:
            process.terminate()
    except ProcessLookupError:
        return
    deadline = time.monotonic() + grace_seconds
    while process.poll() is None and time.monotonic() < deadline:
        time.sleep(0.25)
    if process.poll() is None:
        try:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGKILL)
            else:
                process.kill()
        except ProcessLookupError:
            pass


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run one hash-verified ORCA manifest row in node-local scratch."
    )
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--index", type=int, required=True, help="One-based CSV data-row index")
    parser.add_argument("--site-config", type=Path)
    parser.add_argument("--stage-config", type=Path)
    parser.add_argument("--allocated-cpus", type=int, help="Testing/direct-run override")
    args = parser.parse_args()

    package_root = args.package_root.resolve()
    if not package_root.is_dir():
        raise SystemExit(f"Package root does not exist: {package_root}")
    site_config_from_env = os.environ.get("HPC_SITE_CONFIG", "").strip()
    site_config_path = (
        args.site_config
        or (Path(site_config_from_env) if site_config_from_env else None)
        or package_root / "hpc/config/site_config.json"
    )
    if not site_config_path.is_absolute():
        site_config_path = package_root / site_config_path
    stage_config_path = args.stage_config or package_root / "hpc/config/stage_resources.json"
    if not stage_config_path.is_absolute():
        stage_config_path = package_root / stage_config_path
    site_config = load_json(site_config_path.resolve())
    all_stages = load_json(stage_config_path.resolve()).get("stages", {})
    if args.stage not in all_stages:
        raise SystemExit(f"Unknown stage {args.stage!r}; choose from {sorted(all_stages)}")
    stage_config = all_stages[args.stage]

    manifest = resolve_manifest_path(package_root, args.manifest)
    _, rows = read_manifest(manifest)
    if args.index < 1 or args.index > len(rows):
        raise SystemExit(f"Index {args.index} is outside 1..{len(rows)}")
    row = rows[args.index - 1]
    job_id = safe_component(row_value(row, "job_id", "model_id") or f"row_{args.index:05d}")
    declared_stage = row_value(row, "stage")

    results_raw = str(site_config.get("results_root") or "results/hpc")
    results_root = resolve_package_path(package_root, results_raw, field="results_root", allow_missing=True)
    state_root = package_root / "state"
    state_file = state_root / "tasks" / safe_component(args.stage) / f"{job_id}.json"
    lock_dir = state_root / "locks" / safe_component(args.stage) / f"{job_id}.lock"

    if not enabled_row(row):
        state = {
            "schema_version": "1.0",
            "notice": "CALCULATED MODEL WORKFLOW - NOT EXPERIMENTAL DATA",
            "status": "SKIPPED_DISABLED",
            "timestamp": utc_now(),
            "stage": args.stage,
            "job_id": job_id,
            "manifest": str(manifest.relative_to(package_root)),
            "manifest_row": args.index,
        }
        atomic_json(state_file, state)
        return 0

    if not declared_stage or declared_stage != args.stage:
        raise SystemExit(
            f"Enabled manifest row {args.index} declares stage={declared_stage!r}, not {args.stage!r}"
        )

    approval_issue = approval_problem(row, args.stage)
    if approval_issue:
        atomic_json(
            state_file,
            {
                "schema_version": "1.0",
                "notice": "CALCULATED MODEL WORKFLOW - NOT EXPERIMENTAL DATA",
                "status": "BLOCKED_SCIENTIFIC_GATE",
                "timestamp": utc_now(),
                "stage": args.stage,
                "job_id": job_id,
                "reason": approval_issue,
            },
        )
        print(f"BLOCKED: {job_id}: {approval_issue}", file=sys.stderr)
        return 3

    source_job_dir = job_directory(package_root, row)
    input_path = input_path_for(package_root, source_job_dir, row, stage_config)
    verify_materialized_dependency(package_root, row, source_job_dir, args.stage)
    required_patterns = required_artifact_patterns(stage_config, row)
    output_name = output_name_for(row, input_path)
    hashes = verify_hashes(row, source_job_dir, input_path)
    strict_hash = parse_bool(site_config.get("strict_hash"), default=True)
    if strict_hash and (not hashes["strict_checks_present"] or not hashes["matches"]):
        atomic_json(
            state_file,
            {
                "schema_version": "1.0",
                "notice": "CALCULATED MODEL WORKFLOW - NOT EXPERIMENTAL DATA",
                "status": "BLOCKED_HASH_VALIDATION",
                "timestamp": utc_now(),
                "stage": args.stage,
                "job_id": job_id,
                "hash_validation": hashes,
            },
        )
        print(f"BLOCKED: {job_id}: required input hashes are absent or mismatched", file=sys.stderr)
        return 4

    current_input_hash = str(hashes["actual_input_sha256"])
    try:
        orca_exe = Path(
            os.environ.get("ORCA_EXE") or str(site_config.get("orca_executable", ""))
        ).expanduser()
        if (
            not orca_exe.is_absolute()
            or not orca_exe.is_file()
            or not os.access(orca_exe, os.X_OK)
        ):
            raise TaskFailure(
                f"ORCA_EXE must be an existing executable absolute path: {orca_exe}"
            )
        orca_exe = orca_exe.resolve()
        actual_orca_hash = sha256_file(orca_exe)
        expected_orca_hash = str(site_config.get("orca_sha256") or "").strip().upper()
        valid_expected_orca_hash = len(expected_orca_hash) == 64 and all(
            character in "0123456789ABCDEF" for character in expected_orca_hash
        )
        if strict_hash and not valid_expected_orca_hash:
            raise TaskFailure(
                "strict_hash=true requires orca_sha256 to contain the exact 64-hex Linux ORCA binary hash"
            )
        if expected_orca_hash and expected_orca_hash != actual_orca_hash:
            raise TaskFailure(
                f"ORCA executable hash mismatch: expected={expected_orca_hash} actual={actual_orca_hash}"
            )
        required_helpers = resolve_required_helpers(site_config, orca_exe, row)
    except (OSError, ValueError, TaskFailure) as exc:
        atomic_json(
            state_file,
            {
                "schema_version": "1.0",
                "notice": "CALCULATED MODEL WORKFLOW - NOT EXPERIMENTAL DATA",
                "status": "BLOCKED_EXECUTABLE_IDENTITY",
                "timestamp": utc_now(),
                "stage": args.stage,
                "job_id": job_id,
                "reason": str(exc),
            },
        )
        print(f"BLOCKED: {job_id}: {exc}", file=sys.stderr)
        return 5
    if state_file.is_file():
        old_state = load_json(state_file)
        if completed_state_is_valid(
            old_state,
            current_input_hash,
            package_root,
            required_patterns,
            actual_orca_hash,
            required_helpers,
        ):
            old_state["status"] = "RESUMED_ALREADY_COMPLETE"
            old_state["resume_checked_at"] = utc_now()
            atomic_json(state_file, old_state)
            print(f"RESUME: {job_id} was already complete and hash-valid")
            return 0

    lock_dir.parent.mkdir(parents=True, exist_ok=True)
    try:
        lock_dir.mkdir()
    except FileExistsError:
        age = time.time() - lock_dir.stat().st_mtime
        print(
            f"LOCKED: {job_id}; lock age={age:.0f}s. Inspect with hpc/bin/status.py. "
            "Never delete a lock until the scheduler confirms that the task is no longer running.",
            file=sys.stderr,
        )
        return 75

    process: subprocess.Popen[bytes] | None = None
    interrupted_signal: int | None = None

    def on_signal(signum: int, _frame: Any) -> None:
        nonlocal interrupted_signal
        interrupted_signal = signum
        if process is not None:
            terminate_process(process)

    old_handlers = {
        sig: signal.signal(sig, on_signal)
        for sig in (signal.SIGTERM, signal.SIGINT)
    }

    scratch_dir: Path | None = None
    partial: Path | None = None
    try:
        xtb_exe: Path | None = None
        if args.stage == "goat":
            xtb_raw = str(site_config.get("xtb_executable") or "").strip()
            xtb_candidates = [Path(xtb_raw).expanduser()] if xtb_raw else []
            xtb_candidates.extend(
                (
                    orca_exe.parent / "otool_xtb",
                    orca_exe.parent / "xtb",
                    orca_exe.parent / "xtb-6.7.1pre" / "xtb",
                )
            )
            xtb_exe = next(
                (path.resolve() for path in xtb_candidates if path.is_file() and os.access(path, os.X_OK)),
                None,
            )
            if xtb_exe is None or not xtb_exe.is_absolute():
                raise TaskFailure(
                    "GOAT requires an absolute xTB executable path. Set xtb_executable "
                    "in hpc/config/site_config.json."
                )

        allocation = allocated_cpus(args.allocated_cpus)
        resource_request = parse_orca_resources(input_path)
        if resource_request["required_cpu_slots"] > allocation:
            raise TaskFailure(
                "ORCA input requests more parallel workers than the scheduler allocated: "
                f"required={resource_request['required_cpu_slots']} allocated={allocation}"
            )

        scratch_base_raw = str(site_config.get("scratch_base") or "").strip()
        if scratch_base_raw:
            scratch_base = Path(os.path.expandvars(os.path.expanduser(scratch_base_raw)))
        else:
            scratch_base = Path(
                os.environ.get("SLURM_TMPDIR")
                or os.environ.get("TMPDIR")
                or tempfile.gettempdir()
            )
        scratch_base.mkdir(parents=True, exist_ok=True)
        scheduler = scheduler_context()
        scratch_dir = Path(
            tempfile.mkdtemp(
                prefix=f"orca_{safe_component(scheduler['job_id'] or 'direct')}_{args.index}_{job_id}_",
                dir=str(scratch_base),
            )
        )
        scratch_work = scratch_dir / "work"
        copy_source_to_scratch(source_job_dir, scratch_work)
        scratch_input = scratch_work / input_path.relative_to(source_job_dir)
        if not scratch_input.is_file():
            raise TaskFailure(f"Input was not copied to node-local scratch: {scratch_input}")

        job_result_root = results_root / safe_component(args.stage) / job_id
        attempt_number, partial, attempt_final = next_attempt_directory(job_result_root)
        running_state: dict[str, Any] = {
            "schema_version": "1.0",
            "notice": "CALCULATED MODEL WORKFLOW - NOT EXPERIMENTAL DATA",
            "status": "RUNNING",
            "started_at": utc_now(),
            "stage": args.stage,
            "job_id": job_id,
            "manifest": str(manifest.relative_to(package_root)),
            "manifest_sha256": sha256_file(manifest),
            "manifest_row": args.index,
            "manifest_row_sha256": canonical_row_hash(row),
            "source_job_directory": str(source_job_dir.relative_to(package_root)),
            "input_file": str(input_path.relative_to(package_root)),
            "input_sha256": current_input_hash,
            "output_file": output_name,
            "attempt": attempt_number,
            "attempt_path": str(attempt_final.relative_to(package_root)).replace("\\", "/"),
            "scratch_directory": str(scratch_dir),
            "scheduler": scheduler,
            "allocated_cpus": allocation,
            "orca_resources": resource_request,
            "orca_executable": str(orca_exe),
            "orca_executable_sha256": actual_orca_hash,
            "required_helpers": required_helpers,
            "required_artifact_patterns": list(required_patterns),
            "xtb_executable": str(xtb_exe or ""),
            "hash_validation": hashes,
        }
        atomic_json(state_file, running_state)
        atomic_json(lock_dir / "owner.json", running_state)

        environment = os.environ.copy()
        environment.update(
            {
                "OMP_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1",
                "BLIS_NUM_THREADS": "1",
                "VECLIB_MAXIMUM_THREADS": "1",
                "NUMEXPR_NUM_THREADS": "1",
                "OMP_MAX_ACTIVE_LEVELS": "1",
                "TMPDIR": str(scratch_dir),
            }
        )
        executable_directories = [str(orca_exe.parent)]
        if xtb_exe is not None:
            executable_directories.append(str(xtb_exe.parent))
            environment["XTB_EXE"] = str(xtb_exe)
        environment["PATH"] = os.pathsep.join(executable_directories + [environment.get("PATH", "")])
        if interrupted_signal is not None:
            raise TaskFailure("Task received a termination signal before ORCA launch")
        stdout_path = scratch_work / output_name
        stderr_path = scratch_work / f"{Path(output_name).stem}.stderr.log"
        start_monotonic = time.monotonic()
        with stdout_path.open("wb") as stdout_handle, stderr_path.open("wb") as stderr_handle:
            process = subprocess.Popen(
                [str(orca_exe), scratch_input.name],
                cwd=str(scratch_work),
                env=environment,
                stdout=stdout_handle,
                stderr=stderr_handle,
                start_new_session=True,
            )
            return_code = process.wait()
        duration_seconds = time.monotonic() - start_monotonic

        output_text = stdout_path.read_text(encoding="utf-8", errors="replace") if stdout_path.is_file() else ""
        normal_termination = NORMAL_MARKER in output_text
        required: dict[str, list[str]] = {}
        required_ok = True
        for pattern in required_patterns:
            matches = glob_nonempty(scratch_work, str(pattern))
            required[str(pattern)] = [path.name for path in matches]
            required_ok = required_ok and bool(matches)

        status = "COMPLETED" if return_code == 0 and normal_termination and required_ok else "FAILED"
        if interrupted_signal is not None:
            status = "INTERRUPTED"
        completed_metadata = {
            **running_state,
            "status": status,
            "finished_at": utc_now(),
            "duration_seconds": round(duration_seconds, 3),
            "return_code": return_code,
            "normal_termination_marker": normal_termination,
            "required_artifacts": required,
            "required_artifacts_ok": required_ok,
            "interrupted_signal": interrupted_signal,
        }
        records = copy_back(
            scratch_work,
            partial,
            attempt_final,
            stage_config,
            completed_metadata,
            required_patterns,
        )
        partial = None
        completed_metadata["copied_file_count"] = len(records)
        completed_metadata["attempt_path"] = str(attempt_final.relative_to(package_root)).replace("\\", "/")
        completed_metadata["file_manifest_sha256"] = sha256_file(
            attempt_final / "FILE_MANIFEST_SHA256.csv"
        )
        atomic_json(state_file, completed_metadata)

        keep_failed = parse_bool(site_config.get("keep_failed_scratch"), default=True)
        if status == "COMPLETED" or not keep_failed:
            shutil.rmtree(scratch_dir, ignore_errors=True)
            scratch_dir = None
        if status != "COMPLETED":
            print(
                f"FAILED: {job_id}: rc={return_code}, normal={normal_termination}, "
                f"required_artifacts={required_ok}",
                file=sys.stderr,
            )
            return 1
        print(f"COMPLETED: {job_id} -> {attempt_final}")
        return 0
    except BaseException as exc:
        status = "INTERRUPTED" if interrupted_signal is not None else "FAILED_RUNNER"
        failure = {
            "schema_version": "1.0",
            "notice": "CALCULATED MODEL WORKFLOW - NOT EXPERIMENTAL DATA",
            "status": status,
            "finished_at": utc_now(),
            "stage": args.stage,
            "job_id": job_id,
            "manifest": str(manifest.relative_to(package_root)),
            "manifest_row": args.index,
            "input_sha256": locals().get("current_input_hash", ""),
            "scratch_directory": str(scratch_dir or ""),
            "interrupted_signal": interrupted_signal,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "scheduler": scheduler_context(),
        }
        atomic_json(state_file, failure)
        print(f"{status}: {job_id}: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
    finally:
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)
        if partial is not None and partial.exists():
            shutil.rmtree(partial, ignore_errors=True)
        if lock_dir.exists():
            shutil.rmtree(lock_dir, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
