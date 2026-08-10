from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence


PROJECT = Path(__file__).resolve().parents[1]
HPC_BIN = PROJECT / "hpc" / "bin"
if str(HPC_BIN) not in sys.path:
    sys.path.insert(0, str(HPC_BIN))

from hpc_common import (
    approval_problem,
    dependency_contract,
    enabled_row,
    load_json,
    manifest_name_list,
    resolve_package_path,
    resolve_required_helpers,
    row_value,
    safe_component,
    verify_materialized_dependency_receipt,
)

ORCA_EXE = Path(os.environ.get("ORCA_EXE", "orca"))
MPI_DIR = Path(os.environ.get("MPI_DIR", ""))
NORMAL_MARKER = "****ORCA TERMINATED NORMALLY****"
EXPECTED_BUNDLE_HASH_FIELDS = (
    "expected_input_bundle_sha256",
    "input_bundle_sha256",
)
EXPECTED_INPUT_HASH_FIELDS = (
    "expected_input_sha256",
    "input_file_sha256",
    "input_sha256",
)
EXPECTED_COORDINATE_HASH_FIELDS = (
    "expected_coordinate_sha256",
    "coordinate_sha256",
    "coordinates_sha256",
    "xyzfile_sha256",
    "start_xyz_sha256",
    "conformer_xyz_sha256",
)
XYZFILE_PATTERN = re.compile(
    r"(?im)^\s*\*\s*(?:xyzfile|gzmtfile)\s+[-+]?\d+\s+\d+\s+"
    r'(?:"([^"]+)"|\'([^\']+)\'|(\S+))'
)


def validated_job_directory(row: Mapping[str, str]) -> Path:
    raw = row_value(
        row,
        "portable_job_directory",
        "relative_job_directory",
        "working_directory",
        "job_directory",
    )
    return resolve_package_path(PROJECT, raw, field="job_directory")


def verify_cross_stage_dependency(row: Mapping[str, str]) -> dict[str, Any]:
    child_stage = row_value(row, "stage")
    parent_id, parent_stage, required = dependency_contract(row, child_stage)
    if not parent_id or parent_stage == child_stage:
        raise ValueError("Expected a declared cross-stage dependency")
    job_dir = validated_job_directory(row)
    receipt_error = ""
    try:
        verified = verify_materialized_dependency_receipt(
            PROJECT, row, job_dir, child_stage
        )
        return {"source": "hpc_materialization_receipt", **verified}
    except (OSError, ValueError) as exc:
        receipt_error = str(exc)

    status_path = PROJECT / "logs" / safe_component(parent_stage) / (
        f"{safe_component(parent_id)}.json"
    )
    try:
        status = load_json(status_path)
        if (
            status.get("status") not in {"complete", "already_complete"}
            or status.get("stage") != parent_stage
            or status.get("job_id") != parent_id
            or status.get("required_artifacts_complete") is not True
            or not status.get("input_bundle_sha256")
            or not status.get("orca_executable_sha256")
            or not status.get("orca_version")
        ):
            raise ValueError("local dependency state identity/provenance is incomplete")
        output = Path(str(status.get("output", ""))).resolve()
        try:
            output.relative_to(job_dir.resolve())
        except ValueError as exc:
            raise ValueError("local dependency output is outside the child job directory") from exc
        if (
            not output.is_file()
            or sha256_file(output) != str(status.get("output_sha256", "")).lower()
            or not terminated_normally(output)
        ):
            raise ValueError("local dependency output is absent, altered, or incomplete")
        records: dict[str, Mapping[str, Any]] = {}
        for record in status.get("required_artifacts", []):
            if not isinstance(record, Mapping):
                continue
            name = str(record.get("name", "")).strip()
            if name and name not in records:
                records[name] = record
        for name in required:
            expected_path = (job_dir / name).resolve()
            record = records.get(name)
            if record is None:
                raise ValueError(f"local dependency state does not bind {name}")
            recorded_path = Path(str(record.get("path", ""))).resolve()
            if recorded_path != expected_path:
                raise ValueError(f"local dependency path mismatch for {name}")
            digest = str(record.get("sha256", "")).strip().lower()
            if (
                not digest
                or not expected_path.is_file()
                or sha256_file(expected_path) != digest
            ):
                raise ValueError(f"local dependency artifact is absent or altered: {name}")
        return {
            "source": "local_hash_bound_completion_state",
            "parent_job_id": parent_id,
            "parent_stage": parent_stage,
            "status_file": str(status_path),
            "artifacts": list(required),
        }
    except (OSError, ValueError) as exc:
        raise ValueError(
            f"Unverified cross-stage dependency {parent_stage}/{parent_id}; "
            f"materialization receipt: {receipt_error}; local state: {exc}"
        ) from exc


def topological_rows(rows: Sequence[dict[str, str]]) -> list[dict[str, str]]:
    by_id = {row_value(row, "job_id", "model_id"): row for row in rows}
    if "" in by_id or len(by_id) != len(rows):
        raise ValueError("Enabled manifest rows require unique non-empty job_id values")
    dependencies: dict[str, set[str]] = {}
    for job_id, row in by_id.items():
        dependency = row_value(row, "depends_on_job_id")
        dependencies[job_id] = {dependency} if dependency and dependency in by_id else set()
        if dependency and dependency not in by_id:
            child_stage = row_value(row, "stage")
            _, parent_stage, _ = dependency_contract(row, child_stage)
            if not child_stage or parent_stage == child_stage:
                raise ValueError(
                    f"{job_id} depends on unselected or disabled same-stage job: {dependency}"
                )
            verify_cross_stage_dependency(row)
    ordered: list[dict[str, str]] = []
    remaining = dict(dependencies)
    while remaining:
        ready = sorted(job_id for job_id, deps in remaining.items() if not deps)
        if not ready:
            raise ValueError(f"Dependency cycle in enabled manifest rows: {sorted(remaining)}")
        for job_id in ready:
            ordered.append(by_id[job_id])
            remaining.pop(job_id)
        completed = set(ready)
        for deps in remaining.values():
            deps.difference_update(completed)
    return ordered


def utc_now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8", newline="\n")
    temporary.replace(path)


def atomic_json(path: Path, payload: object) -> None:
    atomic_text(path, json.dumps(payload, ensure_ascii=False, indent=2) + "\n")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def referenced_coordinate_files(input_path: Path) -> list[tuple[str, Path]]:
    text = input_path.read_text(encoding="utf-8", errors="replace")
    references: list[tuple[str, Path]] = []
    seen: set[str] = set()
    for match in XYZFILE_PATTERN.finditer(text):
        token = next(value for value in match.groups() if value is not None).strip()
        candidate = Path(token)
        if not candidate.is_absolute():
            candidate = input_path.parent / candidate
        resolved = candidate.resolve()
        key = os.path.normcase(str(resolved))
        if key not in seen:
            references.append((token, resolved))
            seen.add(key)
    return references


def build_input_bundle(input_path: Path) -> dict[str, Any]:
    """Fingerprint the ORCA input and coordinate files named by xyzfile/gzmtfile."""

    resolved_input = input_path.resolve()
    files: list[dict[str, Any]] = [
        {
            "role": "input",
            "logical_name": input_path.name,
            "path": str(resolved_input),
        }
    ]
    for token, path in referenced_coordinate_files(resolved_input):
        files.append(
            {
                "role": "coordinate",
                "logical_name": token.replace("\\", "/"),
                "path": str(path),
            }
        )
    missing = [item["path"] for item in files if not Path(item["path"]).is_file()]
    if missing:
        raise FileNotFoundError("Missing input-bundle file(s): " + ", ".join(missing))
    for item in files:
        path = Path(item["path"])
        item["size_bytes"] = path.stat().st_size
        item["sha256"] = sha256_file(path)
    identity_payload = [
        {
            "role": item["role"],
            "logical_name": item["logical_name"],
            "size_bytes": item["size_bytes"],
            "sha256": item["sha256"],
        }
        for item in files
    ]
    bundle_hash = hashlib.sha256(
        json.dumps(
            identity_payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
    ).hexdigest()
    return {
        "algorithm": "sha256(canonical-json(role,logical_name,size,sha256))",
        "sha256": bundle_hash,
        "files": files,
    }


def _expected_value(row: Mapping[str, str], fields: Sequence[str]) -> tuple[str, str]:
    lowered = {str(key).casefold(): str(value).strip() for key, value in row.items()}
    for field in fields:
        value = lowered.get(field.casefold(), "")
        if value:
            return field, value.casefold()
    return "", ""


def validate_manifest_hashes(
    row: Mapping[str, str], bundle: Mapping[str, Any]
) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    bundle_field, expected_bundle = _expected_value(row, EXPECTED_BUNDLE_HASH_FIELDS)
    if expected_bundle:
        checks.append(
            {
                "field": bundle_field,
                "expected": expected_bundle,
                "actual": str(bundle["sha256"]).casefold(),
                "matches": expected_bundle == str(bundle["sha256"]).casefold(),
            }
        )
    input_record = next(
        item for item in bundle["files"] if item.get("role") == "input"
    )
    input_field, expected_input = _expected_value(row, EXPECTED_INPUT_HASH_FIELDS)
    if expected_input:
        checks.append(
            {
                "field": input_field,
                "expected": expected_input,
                "actual": str(input_record["sha256"]).casefold(),
                "matches": expected_input == str(input_record["sha256"]).casefold(),
            }
        )
    coordinate_records = [
        item for item in bundle["files"] if item.get("role") == "coordinate"
    ]
    coordinate_field, expected_coordinate = _expected_value(
        row, EXPECTED_COORDINATE_HASH_FIELDS
    )
    if expected_coordinate:
        if len(coordinate_records) == 1:
            actual_coordinate = str(coordinate_records[0]["sha256"]).casefold()
            checks.append(
                {
                    "field": coordinate_field,
                    "expected": expected_coordinate,
                    "actual": actual_coordinate,
                    "matches": expected_coordinate == actual_coordinate,
                }
            )
        else:
            checks.append(
                {
                    "field": coordinate_field,
                    "expected": expected_coordinate,
                    "actual": "",
                    "matches": False,
                    "reason": (
                        "A single coordinate hash was supplied but the input references "
                        f"{len(coordinate_records)} coordinate files"
                    ),
                }
            )
    return {
        "checks": checks,
        "matches": all(bool(check["matches"]) for check in checks),
    }


def manifest_required_artifacts(row: Mapping[str, str]) -> list[str]:
    for field in ("required_artifacts", "required_artifact_files"):
        raw = str(row.get(field, "")).strip()
        if not raw:
            continue
        if raw.startswith("["):
            parsed = json.loads(raw)
            if not isinstance(parsed, list):
                raise ValueError(f"{field} must contain a JSON list or delimited names")
            return [str(value).strip() for value in parsed if str(value).strip()]
        return [value.strip() for value in re.split(r"[;,|]", raw) if value.strip()]
    return []


def infer_required_artifacts(
    row: Mapping[str, str], input_name: str, output_name: str, stage: str
) -> list[str]:
    explicit = manifest_required_artifacts(row)
    output_globs = list(
        manifest_name_list(row, "required_output_globs", allow_globs=True)
    )
    if explicit or output_globs:
        return list(dict.fromkeys([*explicit, *output_globs]))
    stem = Path(output_name).stem or Path(input_name).stem
    context = " ".join(
        (
            stage,
            str(row.get("job_id", "")),
            input_name,
            output_name,
        )
    ).casefold()
    if "dft_opt_freq" in context or ("opt" in context and "freq" in context):
        return [f"{stem}.gbw", f"{stem}.hess", f"{stem}.property.txt", f"{stem}.xyz"]
    if "shield" in context or "coupling" in context or "jhh" in context:
        return [f"{stem}.gbw", f"{stem}.property.txt"]
    return []


def inspect_artifacts(job_dir: Path, names: Sequence[str]) -> tuple[list[dict[str, Any]], bool]:
    records: list[dict[str, Any]] = []
    complete = True
    for name in names:
        token = Path(name)
        if token.is_absolute() or token.name != name or name in {".", ".."}:
            raise ValueError(f"Unsafe required artifact pattern: {name!r}")
        candidates = (
            sorted(job_dir.glob(name))
            if any(character in name for character in "*?[]")
            else [job_dir / name]
        )
        matches = [
            candidate.resolve()
            for candidate in candidates
            if candidate.is_file() and candidate.stat().st_size > 0
        ]
        if not matches:
            complete = False
            records.append(
                {
                    "name": name,
                    "required_pattern": name,
                    "path": "",
                    "exists_nonempty": False,
                    "size_bytes": 0,
                    "mtime_ns": 0,
                    "sha256": "",
                }
            )
            continue
        for candidate in matches:
            records.append(
                {
                    "name": candidate.name,
                    "required_pattern": name,
                    "path": str(candidate),
                    "exists_nonempty": True,
                    "size_bytes": candidate.stat().st_size,
                    "mtime_ns": candidate.stat().st_mtime_ns,
                    "sha256": sha256_file(candidate),
                }
            )
    return records, complete


def detect_orca_version(executable: Path, environment: Mapping[str, str]) -> str:
    try:
        completed = subprocess.run(
            [str(executable), "--version"],
            env=dict(environment),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=15,
            check=False,
        )
        text = completed.stdout + "\n" + completed.stderr
        match = re.search(r"Program\s+Version\s+([0-9][0-9A-Za-z.+_-]*)", text, re.I)
        if match:
            return match.group(1)
        match = re.search(r"Python\s+([0-9][0-9A-Za-z.+_-]*)", text, re.I)
        if match:
            return "Python " + match.group(1)
    except Exception:
        pass
    match = re.search(r"ORCA[_-]([0-9]+(?:\.[0-9]+)+)", str(executable), re.I)
    return match.group(1) if match else "unknown"


def identify_orca(executable: Path, environment: Mapping[str, str]) -> dict[str, Any]:
    resolved = executable.resolve()
    return {
        "path": str(resolved),
        "version": detect_orca_version(resolved, environment),
        "sha256": sha256_file(resolved),
        "size_bytes": resolved.stat().st_size,
        "mtime_ns": resolved.stat().st_mtime_ns,
    }


def load_status(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def completed_state_is_valid(
    previous: Mapping[str, Any],
    input_bundle: Mapping[str, Any],
    orca_identity: Mapping[str, Any],
    output_path: Path,
    artifact_records: Sequence[Mapping[str, Any]],
    artifacts_complete: bool,
    required_helpers: Sequence[Mapping[str, str]],
) -> tuple[bool, list[str]]:
    reasons: list[str] = []
    if previous.get("status") not in {"complete", "already_complete"}:
        reasons.append("no trusted completed status")
    if previous.get("input_bundle_sha256") != input_bundle.get("sha256"):
        reasons.append("input bundle hash changed or was not previously recorded")
    if previous.get("orca_executable_sha256") != orca_identity.get("sha256"):
        reasons.append("ORCA executable hash changed or was not previously recorded")
    if str(previous.get("orca_version", "")) != str(orca_identity.get("version", "")):
        reasons.append("ORCA version changed or was not previously recorded")
    if str(previous.get("output", "")) != str(output_path):
        reasons.append("output path changed or was not previously recorded")
    if not terminated_normally(output_path):
        reasons.append("output lacks the normal-termination marker")
    current_output_hash = sha256_file(output_path) if output_path.is_file() else ""
    if not previous.get("output_sha256") or previous.get("output_sha256") != current_output_hash:
        reasons.append("output hash changed or was not previously recorded")
    if not artifacts_complete:
        reasons.append("one or more required artifacts are missing or empty")
    if list(previous.get("required_helpers", [])) != list(required_helpers):
        reasons.append("required helper identity/hash changed or was not previously recorded")

    def signatures(
        records: Sequence[Mapping[str, Any]],
    ) -> list[tuple[str, str, int, str]]:
        return sorted(
            (
                str(record.get("required_pattern", record.get("name", ""))),
                str(record.get("path", "")),
                int(record.get("size_bytes", 0)),
                str(record.get("sha256", "")),
            )
            for record in records
        )

    previous_records = [
        record
        for record in previous.get("required_artifacts", [])
        if isinstance(record, Mapping)
    ]
    if signatures(previous_records) != signatures(artifact_records):
        reasons.append("required artifact set/hash changed or was not previously recorded")
    return not reasons, reasons


def terminated_normally(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    with path.open("rb") as handle:
        handle.seek(max(0, path.stat().st_size - 65536))
        tail = handle.read().decode("utf-8", errors="replace")
    return NORMAL_MARKER in tail


def process_environment(orca_executable: Path, mpi_directory: Path = MPI_DIR) -> dict[str, str]:
    environment = os.environ.copy()
    prefix = os.pathsep.join((str(orca_executable.resolve().parent), str(mpi_directory)))
    environment["PATH"] = prefix + os.pathsep + environment.get("PATH", "")
    return environment


def run_one(
    row: dict[str, str],
    timeout_seconds: float | None,
    force: bool,
    stage: str,
    status_dir: Path,
    environment: dict[str, str],
    orca_executable: Path,
    orca_identity: Mapping[str, Any],
) -> dict[str, object]:
    job_id = row["job_id"]
    if not enabled_row(row):
        raise ValueError(f"Refusing to execute disabled manifest row: {job_id}")
    declared_stage = row_value(row, "stage")
    if declared_stage != stage:
        raise ValueError(f"Refusing stage mismatch for {job_id}: {declared_stage!r} != {stage!r}")
    gate_issue = approval_problem(row, stage)
    if gate_issue:
        raise ValueError(f"Refusing unapproved manifest row {job_id}: {gate_issue}")
    job_dir = validated_job_directory(row)
    required_helpers = resolve_required_helpers({}, orca_executable, row)
    input_name = row.get("input_file", "").strip()
    output_name = row.get("output_file", "").strip()
    if not input_name or not output_name:
        raise ValueError(f"Manifest row {job_id} lacks input_file/output_file")
    input_path = job_dir / input_name
    output_path = job_dir / output_name
    error_path = job_dir / (Path(output_name).stem + ".err")
    status_path = status_dir / f"{job_id}.json"
    orca_fields = {
        "orca_executable": str(orca_identity["path"]),
        "orca_version": str(orca_identity["version"]),
        "orca_executable_sha256": str(orca_identity["sha256"]),
    }
    if not input_path.exists():
        result = {
            "job_id": job_id,
            "figure_id": row.get("figure_id", ""),
            "stage": stage,
            "status": "missing_input",
            "input": str(input_path),
            "output": str(output_path),
            "checked_at": utc_now(),
            **orca_fields,
        }
        atomic_json(status_path, result)
        return result

    try:
        input_bundle = build_input_bundle(input_path)
    except (OSError, ValueError) as exc:
        result = {
            "job_id": job_id,
            "figure_id": row.get("figure_id", ""),
            "stage": stage,
            "status": "missing_input_bundle",
            "error": str(exc),
            "input": str(input_path),
            "output": str(output_path),
            "checked_at": utc_now(),
            **orca_fields,
        }
        atomic_json(status_path, result)
        return result

    manifest_hash_validation = validate_manifest_hashes(row, input_bundle)
    required_names = infer_required_artifacts(
        row, input_name=input_name, output_name=output_name, stage=stage
    )
    artifact_records, artifacts_complete = inspect_artifacts(job_dir, required_names)
    provenance = {
        "input_bundle_sha256": input_bundle["sha256"],
        "input_bundle": input_bundle,
        "manifest_hash_validation": manifest_hash_validation,
        "required_artifact_names": required_names,
        "required_helpers": required_helpers,
        **orca_fields,
    }
    if not manifest_hash_validation["matches"]:
        result = {
            "job_id": job_id,
            "figure_id": row.get("figure_id", ""),
            "stage": stage,
            "status": "manifest_hash_mismatch",
            "input": str(input_path),
            "output": str(output_path),
            "normal_termination_marker": terminated_normally(output_path),
            "required_artifacts": artifact_records,
            "required_artifacts_complete": artifacts_complete,
            "checked_at": utc_now(),
            **provenance,
        }
        atomic_json(status_path, result)
        return result

    previous = load_status(status_path)
    can_skip, invalidation_reasons = completed_state_is_valid(
        previous=previous,
        input_bundle=input_bundle,
        orca_identity=orca_identity,
        output_path=output_path,
        artifact_records=artifact_records,
        artifacts_complete=artifacts_complete,
        required_helpers=required_helpers,
    )
    if not force and can_skip:
        result = {
            "job_id": job_id,
            "figure_id": row.get("figure_id", ""),
            "stage": stage,
            "status": "already_complete",
            "return_code": 0,
            "normal_termination_marker": True,
            "input": str(input_path),
            "output": str(output_path),
            "output_sha256": sha256_file(output_path),
            "required_artifacts": artifact_records,
            "required_artifacts_complete": True,
            "completion_validation": "input, executable, output and artifacts matched",
            "checked_at": utc_now(),
            **provenance,
        }
        atomic_json(status_path, result)
        return result

    job_dir.mkdir(parents=True, exist_ok=True)
    started = utc_now()
    start_clock = time.monotonic()
    start_wall_ns = time.time_ns()
    result: dict[str, object]
    with output_path.open("w", encoding="utf-8", newline="\n") as stdout, error_path.open(
        "w", encoding="utf-8", newline="\n"
    ) as stderr:
        try:
            completed = subprocess.run(
                [str(orca_executable), input_name],
                cwd=job_dir,
                env=environment,
                stdout=stdout,
                stderr=stderr,
                text=True,
                timeout=timeout_seconds,
                check=False,
            )
            return_code = int(completed.returncode)
            normal = terminated_normally(output_path)
            artifact_records, artifacts_complete = inspect_artifacts(job_dir, required_names)
            artifacts_from_this_run = all(
                int(record.get("mtime_ns", 0)) >= start_wall_ns
                for record in artifact_records
            )
            status = (
                "complete"
                if return_code == 0
                and normal
                and artifacts_complete
                and artifacts_from_this_run
                else "failed"
            )
            failure_reasons: list[str] = []
            if return_code != 0:
                failure_reasons.append(f"ORCA return code {return_code}")
            if not normal:
                failure_reasons.append("normal-termination marker absent")
            if not artifacts_complete:
                failure_reasons.append("required artifact missing or empty")
            if artifacts_complete and not artifacts_from_this_run:
                failure_reasons.append("required artifact was not produced by this run")
            result = {
                "job_id": job_id,
                "figure_id": row.get("figure_id", ""),
                "stage": stage,
                "status": status,
                "reason": "; ".join(failure_reasons),
                "return_code": return_code,
                "normal_termination_marker": normal,
                "started_at": started,
                "finished_at": utc_now(),
                "elapsed_seconds": round(time.monotonic() - start_clock, 3),
                "input": str(input_path),
                "output": str(output_path),
                "output_sha256": sha256_file(output_path) if output_path.is_file() else "",
                "stderr": str(error_path),
                "required_artifacts": artifact_records,
                "required_artifacts_complete": artifacts_complete,
                "required_artifacts_produced_this_run": artifacts_from_this_run,
                "rerun_reasons": ["--force"] if force else invalidation_reasons,
                **provenance,
            }
        except subprocess.TimeoutExpired:
            artifact_records, artifacts_complete = inspect_artifacts(job_dir, required_names)
            result = {
                "job_id": job_id,
                "figure_id": row.get("figure_id", ""),
                "stage": stage,
                "status": "timeout",
                "return_code": None,
                "normal_termination_marker": False,
                "started_at": started,
                "finished_at": utc_now(),
                "elapsed_seconds": round(time.monotonic() - start_clock, 3),
                "timeout_seconds": timeout_seconds,
                "input": str(input_path),
                "output": str(output_path),
                "output_sha256": sha256_file(output_path) if output_path.is_file() else "",
                "stderr": str(error_path),
                "required_artifacts": artifact_records,
                "required_artifacts_complete": artifacts_complete,
                "rerun_reasons": ["--force"] if force else invalidation_reasons,
                **provenance,
            }
        except Exception as exc:
            artifact_records, artifacts_complete = inspect_artifacts(job_dir, required_names)
            result = {
                "job_id": job_id,
                "figure_id": row.get("figure_id", ""),
                "stage": stage,
                "status": "runner_exception",
                "error": repr(exc),
                "started_at": started,
                "finished_at": utc_now(),
                "elapsed_seconds": round(time.monotonic() - start_clock, 3),
                "input": str(input_path),
                "output": str(output_path),
                "output_sha256": sha256_file(output_path) if output_path.is_file() else "",
                "stderr": str(error_path),
                "required_artifacts": artifact_records,
                "required_artifacts_complete": artifacts_complete,
                "rerun_reasons": ["--force"] if force else invalidation_reasons,
                **provenance,
            }
    atomic_json(status_path, result)
    return result


def write_aggregate(rows: list[dict[str, object]], json_path: Path, csv_path: Path) -> None:
    ordered = sorted(rows, key=lambda row: str(row.get("job_id", "")))
    atomic_json(json_path, ordered)
    fields: list[str] = []
    for row in ordered:
        for key in row:
            if key not in fields:
                fields.append(key)
    temporary = csv_path.with_suffix(csv_path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(ordered)
    temporary.replace(csv_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate an ORCA CSV manifest. Default is dry-run; --execute is required to launch ORCA."
    )
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--only", nargs="*", help="Optional job IDs or Figure3 IDs")
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--timeout-hours", type=float, default=0.0, help="0 disables the timeout")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--stage", required=True)
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Launch ORCA after all path, hash, stage, dependency and scientific-gate checks pass",
    )
    parser.add_argument(
        "--orca-exe",
        type=Path,
        default=ORCA_EXE,
        help="ORCA executable (default: ORCA_EXE environment variable, then PATH)",
    )
    parser.add_argument(
        "--status-root",
        type=Path,
        help="Optional status directory; default is PROJECT/logs/<stage>.",
    )
    args = parser.parse_args()
    if args.max_workers < 1:
        raise ValueError("--max-workers must be at least 1")
    manifest = args.manifest.resolve()
    if not manifest.exists():
        raise FileNotFoundError(manifest)
    try:
        manifest.relative_to(PROJECT.resolve())
    except ValueError as exc:
        raise ValueError("Manifest must remain inside the package root") from exc
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    only = set(args.only or [])
    if only:
        rows = [row for row in rows if row.get("job_id") in only or row.get("figure_id") in only]
    selected: list[dict[str, str]] = []
    for row in rows:
        declared_stage = row_value(row, "stage")
        if declared_stage != args.stage:
            continue
        if not enabled_row(row):
            continue
        gate_issue = approval_problem(row, args.stage)
        if gate_issue:
            raise ValueError(
                f"Scientific gate rejected {row_value(row, 'job_id', 'model_id')!r}: {gate_issue}"
            )
        selected.append(row)
    if not selected:
        raise RuntimeError("No enabled, stage-matching, explicitly approved jobs selected")
    rows = topological_rows(selected)
    has_dependencies = any(row_value(row, "depends_on_job_id") for row in rows)
    if has_dependencies and args.max_workers != 1:
        raise ValueError("Manifests with dependencies require --max-workers 1 in this local runner")

    planned: list[dict[str, Any]] = []
    for row in rows:
        job_id = row_value(row, "job_id", "model_id")
        job_dir = validated_job_directory(row)
        input_name = row_value(row, "input_file")
        if not input_name or Path(input_name).is_absolute() or Path(input_name).name != input_name:
            raise ValueError(f"Unsafe or missing input_file for {job_id}: {input_name!r}")
        input_path = (job_dir / input_name).resolve()
        try:
            input_path.relative_to(job_dir.resolve())
        except ValueError as exc:
            raise ValueError(f"Input escapes job directory for {job_id}") from exc
        bundle = build_input_bundle(input_path)
        hash_validation = validate_manifest_hashes(row, bundle)
        if not hash_validation["checks"] or not hash_validation["matches"]:
            raise ValueError(f"Missing or mismatched manifest hashes for {job_id}")
        helper_records = (
            resolve_required_helpers({}, args.orca_exe.resolve(), row)
            if manifest_name_list(row, "required_helpers")
            else []
        )
        planned.append(
            {
                "job_id": job_id,
                "stage": args.stage,
                "job_directory": str(job_dir.relative_to(PROJECT)).replace("\\", "/"),
                "input_file": input_name,
                "depends_on_job_id": row_value(row, "depends_on_job_id"),
                "input_bundle_sha256": bundle["sha256"],
                "required_helpers": helper_records,
            }
        )
    print(
        json.dumps(
            {
                "mode": "EXECUTE" if args.execute else "DRY_RUN",
                "stage": args.stage,
                "job_count": len(planned),
                "jobs": planned,
            },
            indent=2,
        )
    )
    if not args.execute:
        return

    orca_executable = args.orca_exe.resolve()
    if not orca_executable.is_file():
        raise FileNotFoundError(orca_executable)

    status_root = (
        args.status_root.resolve()
        if args.status_root is not None
        else PROJECT / "logs" / args.stage
    )
    status_root.mkdir(parents=True, exist_ok=True)
    aggregate_json = status_root / "aggregate_status.json"
    aggregate_csv = status_root / "aggregate_status.csv"
    timeout_seconds = args.timeout_hours * 3600.0 if args.timeout_hours > 0 else None
    environment = process_environment(orca_executable)
    orca_identity = identify_orca(orca_executable, environment)
    results: list[dict[str, object]] = []
    lock = threading.Lock()
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = {
            executor.submit(
                run_one,
                row,
                timeout_seconds,
                args.force,
                args.stage,
                status_root,
                environment,
                orca_executable,
                orca_identity,
            ): row["job_id"]
            for row in rows
        }
        for future in as_completed(futures):
            try:
                result = future.result()
            except Exception as exc:
                result = {
                    "job_id": futures[future],
                    "stage": args.stage,
                    "status": "unhandled_runner_exception",
                    "error": repr(exc),
                    "finished_at": utc_now(),
                }
            with lock:
                results.append(result)
                write_aggregate(results, aggregate_json, aggregate_csv)
            print(f"{result.get('job_id')}: {result.get('status')}", flush=True)
    failures = [row for row in results if row.get("status") not in {"complete", "already_complete"}]
    print(f"Finished {len(results)} jobs; failures={len(failures)}")
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()
