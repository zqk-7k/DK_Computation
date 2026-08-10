#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Mapping

from hpc_common import (
    NORMAL_MARKER,
    atomic_json,
    enabled_row,
    input_coordinate_references,
    load_json,
    read_manifest,
    resolve_manifest_path,
    resolve_package_path,
    row_value,
    safe_component,
    sha256_file,
    utc_now,
    write_csv,
)


COMPLETE = {"COMPLETED", "RESUMED_ALREADY_COMPLETE"}
GOAT_PREPARATION_BUNDLE = (
    "goat.inp",
    "start.xyz",
    "atom_map.csv",
    "isotope_map.csv",
    "structure_validation.json",
)


def source_job_directory(root: Path, row: Mapping[str, str]) -> Path:
    return resolve_package_path(
        root,
        row_value(
            row,
            "portable_job_directory",
            "relative_job_directory",
            "working_directory",
            "job_directory",
        ),
        field="job_directory",
    )


def input_path(job_dir: Path, row: Mapping[str, str], default: str) -> Path:
    raw = row_value(row, "portable_input_file", "input_file") or default
    path = job_dir / Path(raw.replace("\\", "/")).name
    if not path.is_file():
        path = job_dir / raw.replace("\\", "/")
    if not path.is_file():
        raise FileNotFoundError(path)
    return path.resolve()


def goat_runner_bundle_hash(input_file: Path) -> str:
    files = [input_file, *input_coordinate_references(input_file)]
    root = input_file.parent.resolve()
    unique = {os.path.normcase(str(path.resolve())): path.resolve() for path in files}
    digest = hashlib.sha256()
    for path in sorted(unique.values(), key=lambda value: os.path.normcase(str(value))):
        relative = path.relative_to(root).as_posix().encode("utf-8")
        digest.update(len(relative).to_bytes(4, "big"))
        digest.update(relative)
        digest.update(path.read_bytes())
    return digest.hexdigest().upper()


def goat_preparation_bundle_hash(job_dir: Path) -> str:
    digest = hashlib.sha256()
    for name in GOAT_PREPARATION_BUNDLE:
        path = job_dir / name
        if not path.is_file():
            raise FileNotFoundError(path)
        digest.update(name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest().upper()


def inline_input_bundle_hash(path: Path) -> str:
    data = path.read_bytes()
    identity = [
        {
            "role": "input",
            "logical_name": path.name,
            "size_bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
        }
    ]
    return hashlib.sha256(
        json.dumps(identity, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    ).hexdigest().upper()


def attempt_directory(root: Path, state: Mapping[str, Any]) -> Path:
    raw = str(state.get("attempt_path", ""))
    path = Path(raw)
    if not path.is_absolute():
        path = root / path
    resolved = path.resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"Attempt directory escapes package root: {resolved}") from exc
    if not resolved.is_dir():
        raise FileNotFoundError(resolved)
    return resolved


def verified_ledger(attempt: Path) -> dict[str, str]:
    ledger_path = attempt / "FILE_MANIFEST_SHA256.csv"
    if not ledger_path.is_file():
        raise FileNotFoundError(ledger_path)
    with ledger_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    ledger: dict[str, str] = {}
    for row in rows:
        name = str(row.get("file", "")).strip()
        expected = str(row.get("sha256", "")).strip().upper()
        if not name or name in ledger or len(expected) != 64:
            raise ValueError(f"Malformed/duplicate FILE_MANIFEST row: {row}")
        path = attempt / name
        if not path.is_file() or sha256_file(path) != expected:
            raise ValueError(f"Attempt artifact hash mismatch: {path}")
        ledger[name] = expected
    return ledger


def ledger_file(attempt: Path, ledger: Mapping[str, str], name: str) -> Path:
    path = attempt / name
    if name not in ledger or not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Required attempt file is absent from verified ledger: {path}")
    return path.resolve()


def choose_artifact(attempt: Path, ledger: Mapping[str, str], preferred: str, pattern: str) -> Path:
    if preferred and preferred in ledger:
        return ledger_file(attempt, ledger, preferred)
    matches = [path for path in attempt.glob(pattern) if path.name in ledger and path.stat().st_size > 0]
    if len(matches) != 1:
        raise ValueError(f"Expected one {pattern} artifact in {attempt}, found {[p.name for p in matches]}")
    return matches[0].resolve()


def export_goat(
    root: Path,
    row: dict[str, str],
    hpc_state: dict[str, Any],
    attempt: Path,
    ledger: Mapping[str, str],
) -> tuple[dict[str, Any], dict[str, str]]:
    job_dir = source_job_directory(root, row)
    inp = input_path(job_dir, row, "goat.inp")
    current_input_hash = sha256_file(inp)
    expected_input_hash = row_value(row, "expected_input_sha256", "goat_input_sha256", "input_sha256").upper()
    if not expected_input_hash or current_input_hash != expected_input_hash:
        raise ValueError(f"GOAT source input changed: {row['job_id']}")
    expected_preparation = row_value(
        row, "expected_input_bundle_sha256", "input_bundle_sha256"
    ).upper()
    if expected_preparation and all((job_dir / name).is_file() for name in GOAT_PREPARATION_BUNDLE):
        if goat_preparation_bundle_hash(job_dir) != expected_preparation:
            raise ValueError(f"GOAT preparation bundle changed: {row['job_id']}")
    runner_bundle = goat_runner_bundle_hash(inp)
    output_name = str(hpc_state.get("output_file") or row_value(row, "output_file") or "goat.out")
    output = ledger_file(attempt, ledger, Path(output_name).name)
    if NORMAL_MARKER not in output.read_text(encoding="utf-8", errors="replace"):
        raise ValueError(f"GOAT output lacks normal termination: {output}")
    ensemble = choose_artifact(attempt, ledger, "goat.finalensemble.xyz", "*.finalensemble.xyz")
    globalminimum = choose_artifact(attempt, ledger, "goat.globalminimum.xyz", "*.globalminimum.xyz")
    attempt_record = {
        "status": "completed",
        "return_code": 0,
        "source_input_unchanged": True,
        "input_bundle_sha256": runner_bundle.lower(),
        "source_input_bundle_sha256_after": runner_bundle.lower(),
        "output_file": str(output),
        "output_sha256": sha256_file(output).lower(),
        "normal_termination_marker": True,
        "artifacts": {
            "finalensemble": str(ensemble),
            "globalminimum": str(globalminimum),
        },
        "hpc_attempt_directory": str(attempt),
        "hpc_file_manifest_sha256": sha256_file(attempt / "FILE_MANIFEST_SHA256.csv").lower(),
    }
    legacy = {
        "schema_version": "HPC_BRIDGE_TO_LEGACY_GATE_STATE_1.0",
        "notice": "CALCULATED MODEL WORKFLOW - NOT EXPERIMENTAL DATA",
        "status": "completed",
        "job_id": row["job_id"],
        "figure_id": row["figure_id"],
        "scenario_id": row.get("scenario_id", ""),
        "input_bundle_sha256": runner_bundle.lower(),
        "attempts": [attempt_record],
        "hpc_source_state": str(root / "state/tasks/goat" / f"{safe_component(row['job_id'])}.json"),
        "hpc_source_state_sha256": sha256_file(root / "state/tasks/goat" / f"{safe_component(row['job_id'])}.json").lower(),
        "exported_at": utc_now(),
    }
    compat = dict(row)
    compat["job_directory"] = str(job_dir.relative_to(root)).replace("\\", "/")
    compat["input_bundle_sha256"] = runner_bundle.lower()
    return legacy, compat


def export_optfreq(
    root: Path,
    row: dict[str, str],
    hpc_state: dict[str, Any],
    attempt: Path,
    ledger: Mapping[str, str],
) -> tuple[dict[str, Any], dict[str, str]]:
    job_dir = source_job_directory(root, row)
    inp = input_path(job_dir, row, "dft_opt_freq.inp")
    input_hash = sha256_file(inp)
    expected_input = row_value(row, "expected_input_sha256", "input_sha256").upper()
    if not expected_input or input_hash != expected_input:
        raise ValueError(f"Opt/freq source input changed: {row['job_id']}")
    bundle_hash = inline_input_bundle_hash(inp)
    expected_bundle = row_value(row, "expected_input_bundle_sha256", "input_bundle_sha256").upper()
    if not expected_bundle or bundle_hash != expected_bundle:
        raise ValueError(f"Opt/freq input bundle changed: {row['job_id']}")
    provenance = job_dir / "provenance.json"
    if not provenance.is_file():
        raise FileNotFoundError(provenance)
    output_name = str(hpc_state.get("output_file") or row_value(row, "output_file") or "dft_opt_freq.out")
    output = ledger_file(attempt, ledger, Path(output_name).name)
    if NORMAL_MARKER not in output.read_text(encoding="utf-8", errors="replace"):
        raise ValueError(f"Opt/freq output lacks normal termination: {output}")
    stem = inp.stem
    roles = {
        "gbw": choose_artifact(attempt, ledger, f"{stem}.gbw", "*.gbw"),
        "hess": choose_artifact(attempt, ledger, f"{stem}.hess", "*.hess"),
        "property": choose_artifact(attempt, ledger, f"{stem}.property.txt", "*.property.txt"),
        "optimized_xyz": choose_artifact(attempt, ledger, f"{stem}.xyz", "*.xyz"),
    }
    artifact_records = {
        role: {
            "path": str(path),
            "sha256": sha256_file(path).lower(),
            "exists_nonempty": True,
        }
        for role, path in roles.items()
    }
    resource_policy = root / "config/resource_policy.json"
    attempt_record = {
        "status": "completed",
        "return_code": 0,
        "normal_termination_marker": True,
        "source_input_unchanged": True,
        "input_sha256": input_hash.lower(),
        "input_bundle_sha256": bundle_hash.lower(),
        "output_file": str(output),
        "output_sha256": sha256_file(output).lower(),
        "artifacts": artifact_records,
        "hpc_attempt_directory": str(attempt),
        "hpc_file_manifest_sha256": sha256_file(attempt / "FILE_MANIFEST_SHA256.csv").lower(),
        "orca_executable_sha256": "",
    }
    legacy = {
        "schema_version": "HPC_BRIDGE_TO_LEGACY_GATE_STATE_1.0",
        "notice": "CALCULATED MODEL WORKFLOW - NOT EXPERIMENTAL DATA",
        "status": "completed",
        "job_id": row["job_id"],
        "figure_id": row["figure_id"],
        "scenario_id": row["scenario_id"],
        "input_bundle_sha256": bundle_hash.lower(),
        "provenance_sha256": sha256_file(provenance).lower(),
        "resource_policy_sha256": sha256_file(resource_policy).lower() if resource_policy.is_file() else "",
        "attempts": [attempt_record],
        "hpc_source_state": str(root / "state/tasks/optfreq" / f"{safe_component(row['job_id'])}.json"),
        "hpc_source_state_sha256": sha256_file(root / "state/tasks/optfreq" / f"{safe_component(row['job_id'])}.json").lower(),
        "exported_at": utc_now(),
    }
    compat = dict(row)
    compat["job_directory"] = str(job_dir.relative_to(root)).replace("\\", "/")
    compat["input_sha256"] = input_hash.lower()
    compat["expected_input_bundle_sha256"] = bundle_hash.lower()
    return legacy, compat


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Export hash-verified HPC states to the legacy scientific-gate schema."
    )
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--stage", choices=("goat", "optfreq"), required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-state-root", type=Path)
    parser.add_argument("--compat-manifest", type=Path)
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()

    root = args.package_root.resolve()
    manifest = resolve_manifest_path(root, args.manifest)
    fields, rows = read_manifest(manifest)
    default_state = "logs/goat_explicit_batch/jobs" if args.stage == "goat" else "logs/opt_freq_serial/jobs"
    state_output = args.output_state_root or Path(default_state)
    if not state_output.is_absolute():
        state_output = root / state_output
    state_output = state_output.resolve()
    state_output.relative_to(root)
    compat_output = args.compat_manifest or Path(f"data/{args.stage}_manifest_gate_compatible.csv")
    if not compat_output.is_absolute():
        compat_output = root / compat_output
    compat_output = compat_output.resolve()
    compat_output.relative_to(root)
    if compat_output.exists() and not args.replace:
        raise SystemExit(f"Compatibility manifest exists; inspect it or use --replace: {compat_output}")

    converted: list[dict[str, str]] = []
    state_hashes: list[dict[str, str]] = []
    for row in rows:
        if not enabled_row(row):
            converted.append(dict(row))
            continue
        job_id = row_value(row, "job_id", "model_id")
        source_state_path = root / "state/tasks" / args.stage / f"{safe_component(job_id)}.json"
        if not source_state_path.is_file():
            raise SystemExit(f"HPC state is missing: {source_state_path}")
        hpc_state = load_json(source_state_path)
        if hpc_state.get("status") not in COMPLETE or hpc_state.get("job_id") != job_id:
            raise SystemExit(f"HPC task is not completed and identity-bound: {source_state_path}")
        attempt = attempt_directory(root, hpc_state)
        ledger = verified_ledger(attempt)
        if args.stage == "goat":
            legacy, compat = export_goat(root, row, hpc_state, attempt, ledger)
        else:
            legacy, compat = export_optfreq(root, row, hpc_state, attempt, ledger)
        output_state = state_output / f"{safe_component(job_id)}.json"
        if output_state.exists() and not args.replace:
            raise SystemExit(f"Legacy-compatible state exists; use --replace after inspection: {output_state}")
        atomic_json(output_state, legacy)
        state_hashes.append(
            {
                "job_id": job_id,
                "source_hpc_state_sha256": sha256_file(source_state_path),
                "exported_state_sha256": sha256_file(output_state),
            }
        )
        converted.append(compat)

    compat_fields = list(fields)
    for field in ("job_directory", "input_bundle_sha256", "input_sha256", "expected_input_bundle_sha256"):
        if field not in compat_fields:
            compat_fields.append(field)
    write_csv(compat_output, compat_fields, converted)
    audit = {
        "schema_version": "1.0",
        "notice": "CALCULATED MODEL WORKFLOW - NOT EXPERIMENTAL DATA",
        "timestamp": utc_now(),
        "stage": args.stage,
        "source_manifest": str(manifest.relative_to(root)).replace("\\", "/"),
        "source_manifest_sha256": sha256_file(manifest),
        "compatibility_manifest": str(compat_output.relative_to(root)).replace("\\", "/"),
        "compatibility_manifest_sha256": sha256_file(compat_output),
        "state_output_root": str(state_output.relative_to(root)).replace("\\", "/"),
        "exported_enabled_states": len(state_hashes),
        "state_hashes": state_hashes,
    }
    atomic_json(compat_output.with_suffix(compat_output.suffix + ".bridge_audit.json"), audit)
    print(f"Exported {len(state_hashes)} {args.stage} states -> {state_output}")
    print(f"Gate-compatible manifest -> {compat_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
