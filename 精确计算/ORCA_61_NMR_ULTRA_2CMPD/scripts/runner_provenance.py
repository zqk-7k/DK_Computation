"""Validate ORCA output hashes against local or portable-HPC runner state."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


class RunnerProvenanceError(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def canonical_row_hash(row: Mapping[str, Any]) -> str:
    payload = json.dumps(dict(row), sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest().upper()


def _inside(path: Path, root: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    try:
        resolved.relative_to(root.resolve())
    except ValueError as exc:
        raise RunnerProvenanceError(f"{label} escapes package root: {resolved}") from exc
    return resolved


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise RunnerProvenanceError(f"Missing runner state: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise RunnerProvenanceError(f"Runner state is not a JSON object: {path}")
    return payload


def _validate_hpc_state(
    state_path: Path,
    state: Mapping[str, Any],
    row: Mapping[str, str],
    stage: str,
    package_root: Path,
) -> dict[str, Any]:
    if state.get("status") not in {"COMPLETED", "RESUMED_ALREADY_COMPLETE"}:
        raise RunnerProvenanceError(f"HPC runner state is not complete: {state_path}")
    if state.get("job_id") != row.get("job_id") or state.get("stage") != stage:
        raise RunnerProvenanceError(f"HPC runner identity mismatch: {state_path}")
    if str(state.get("manifest_row_sha256", "")).upper() != canonical_row_hash(row):
        raise RunnerProvenanceError(f"HPC manifest-row hash mismatch: {state_path}")
    job_directory = Path(str(row["job_directory"]))
    if not job_directory.is_absolute():
        job_directory = package_root / job_directory
    input_path = _inside(job_directory / row["input_file"], package_root, "input")
    if sha256_file(input_path) != str(state.get("input_sha256", "")).upper():
        raise RunnerProvenanceError(f"HPC input hash mismatch: {state_path}")
    hash_validation = state.get("hash_validation")
    if not isinstance(hash_validation, Mapping) or not hash_validation.get("strict_checks_present"):
        raise RunnerProvenanceError(f"HPC runner lacks strict manifest hash checks: {state_path}")
    if hash_validation.get("matches") is not True:
        raise RunnerProvenanceError(f"HPC runner manifest hash validation failed: {state_path}")
    if state.get("normal_termination_marker") is not True or state.get("required_artifacts_ok") is not True:
        raise RunnerProvenanceError(f"HPC runner completion/artifact checks failed: {state_path}")
    attempt = Path(str(state.get("attempt_path", "")))
    if not attempt.is_absolute():
        attempt = package_root / attempt
    attempt = _inside(attempt, package_root, "attempt_path")
    output = _inside(attempt / str(state.get("output_file", "")), attempt, "HPC output")
    file_manifest = attempt / "FILE_MANIFEST_SHA256.csv"
    if not output.is_file() or not file_manifest.is_file():
        raise RunnerProvenanceError(f"HPC output/file manifest is missing: {attempt}")
    expected_manifest_hash = str(state.get("file_manifest_sha256", "")).upper()
    if not expected_manifest_hash or sha256_file(file_manifest) != expected_manifest_hash:
        raise RunnerProvenanceError(f"HPC file-manifest hash mismatch: {file_manifest}")
    with file_manifest.open("r", encoding="utf-8-sig", newline="") as handle:
        records = list(csv.DictReader(handle))
    by_name = {record.get("file", ""): record for record in records}
    for name, record in by_name.items():
        if not name:
            raise RunnerProvenanceError(f"HPC file manifest has a blank artifact name: {file_manifest}")
        artifact = _inside(attempt / name, attempt, "HPC recorded artifact")
        if not artifact.is_file():
            raise RunnerProvenanceError(f"HPC recorded artifact is missing: {artifact}")
        if sha256_file(artifact) != str(record.get("sha256", "")).upper():
            raise RunnerProvenanceError(f"HPC artifact hash mismatch: {artifact}")
    if output.name not in by_name:
        raise RunnerProvenanceError(f"HPC output absent from file hash manifest: {output}")
    text = output.read_text(encoding="utf-8", errors="replace")
    if "****ORCA TERMINATED NORMALLY****" not in text:
        raise RunnerProvenanceError(f"HPC output lacks ORCA normal termination: {output}")
    metadata_path = attempt / "run_metadata.json"
    metadata = _read_json(metadata_path)
    if metadata.get("job_id") != row.get("job_id") or metadata.get("stage") != stage:
        raise RunnerProvenanceError(f"HPC run_metadata identity mismatch: {metadata_path}")
    return {
        "runner_kind": "portable_hpc",
        "runner_state": str(state_path),
        "runner_state_sha256": sha256_file(state_path),
        "output_path": output,
        "output_sha256": sha256_file(output),
        "input_sha256": sha256_file(input_path),
        "attempt_path": str(attempt),
        "orca_version": state.get("orca_version", metadata.get("orca_version", "")),
    }


def _validate_local_state(
    state_path: Path,
    state: Mapping[str, Any],
    row: Mapping[str, str],
    stage: str,
    package_root: Path,
) -> dict[str, Any]:
    if state.get("status") not in {"complete", "already_complete"}:
        raise RunnerProvenanceError(f"Local runner state is not complete: {state_path}")
    if state.get("job_id") != row.get("job_id") or state.get("stage") != stage:
        raise RunnerProvenanceError(f"Local runner identity mismatch: {state_path}")
    validation = state.get("manifest_hash_validation")
    if not isinstance(validation, Mapping) or validation.get("matches") is not True:
        raise RunnerProvenanceError(f"Local runner manifest hash validation failed: {state_path}")
    checks = validation.get("checks")
    if not isinstance(checks, list) or not checks or any(
        not isinstance(check, Mapping) or check.get("matches") is not True
        for check in checks
    ):
        raise RunnerProvenanceError(
            f"Local runner lacks non-empty, successful manifest hash checks: {state_path}"
        )
    job_directory = Path(str(row["job_directory"]))
    if not job_directory.is_absolute():
        job_directory = package_root / job_directory
    input_path = _inside(job_directory / row["input_file"], package_root, "input")
    expected_input = str(
        row.get("input_sha256", "") or row.get("expected_input_sha256", "")
    ).upper()
    if not expected_input or sha256_file(input_path) != expected_input:
        raise RunnerProvenanceError(f"Local manifest/input hash mismatch: {state_path}")
    bundle = state.get("input_bundle")
    bundle_files = bundle.get("files") if isinstance(bundle, Mapping) else None
    if not isinstance(bundle_files, list):
        raise RunnerProvenanceError(f"Local runner lacks input-bundle provenance: {state_path}")
    input_records = [
        item for item in bundle_files
        if isinstance(item, Mapping) and item.get("role") == "input"
    ]
    if len(input_records) != 1 or str(input_records[0].get("sha256", "")).upper() != expected_input:
        raise RunnerProvenanceError(f"Local runner input-bundle hash mismatch: {state_path}")
    output = _inside(Path(str(state.get("output", ""))), package_root, "local output")
    if not output.is_file() or sha256_file(output) != str(state.get("output_sha256", "")).upper():
        raise RunnerProvenanceError(f"Local runner output hash mismatch: {state_path}")
    if "****ORCA TERMINATED NORMALLY****" not in output.read_text(
        encoding="utf-8", errors="replace"
    ):
        raise RunnerProvenanceError(f"Local runner output lacks ORCA normal termination: {output}")
    if state.get("required_artifacts_complete") is not True:
        raise RunnerProvenanceError(f"Local runner required artifacts are incomplete: {state_path}")
    artifact_records = state.get("required_artifacts")
    if not isinstance(artifact_records, list):
        raise RunnerProvenanceError(f"Local runner lacks required-artifact provenance: {state_path}")
    for record in artifact_records:
        if not isinstance(record, Mapping):
            raise RunnerProvenanceError(f"Malformed local artifact record: {state_path}")
        artifact = _inside(Path(str(record.get("path", ""))), package_root, "local artifact")
        expected = str(record.get("sha256", "")).upper()
        if not artifact.is_file() or not expected or sha256_file(artifact) != expected:
            raise RunnerProvenanceError(f"Local runner artifact hash mismatch: {artifact}")
    return {
        "runner_kind": "local_manifest_runner",
        "runner_state": str(state_path),
        "runner_state_sha256": sha256_file(state_path),
        "output_path": output,
        "output_sha256": sha256_file(output),
        "input_sha256": expected_input,
        "attempt_path": str(output.parent),
        "orca_version": state.get("orca_version", ""),
    }


def validate_runner_output(
    row: Mapping[str, str],
    stage: str,
    package_root: Path,
    state_root: Path | None = None,
) -> dict[str, Any]:
    """Resolve a trusted output only through a hash-verified runner state."""
    job_id = str(row.get("job_id", "")).strip()
    if not job_id:
        raise RunnerProvenanceError("Manifest row lacks job_id")
    candidates: list[Path] = []
    if state_root is not None:
        candidates.append(state_root.expanduser().resolve() / f"{job_id}.json")
    else:
        candidates.extend(
            (
                package_root / "state" / "tasks" / stage / f"{job_id}.json",
                package_root / "logs" / stage / f"{job_id}.json",
            )
        )
    state_path = next((path for path in candidates if path.is_file()), None)
    if state_path is None:
        raise RunnerProvenanceError(
            f"No trusted runner state for {job_id}; tried: {', '.join(str(path) for path in candidates)}"
        )
    state = _read_json(state_path)
    if str(state.get("status", "")).isupper():
        return _validate_hpc_state(state_path, state, row, stage, package_root)
    return _validate_local_state(state_path, state, row, stage, package_root)
