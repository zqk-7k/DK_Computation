#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
from typing import Any

from hpc_common import (
    atomic_json,
    input_coordinate_references,
    is_windows_absolute,
    load_json,
    read_manifest,
    row_value,
    sha256_file,
    utc_now,
    write_csv,
)


GOAT_BUNDLE_NAMES = (
    "goat.inp",
    "start.xyz",
    "atom_map.csv",
    "isotope_map.csv",
    "structure_validation.json",
)


def inside_root(root: Path, raw: str, field: str) -> Path:
    token = raw.strip()
    if not token:
        raise ValueError(f"Missing {field}")
    if is_windows_absolute(token) and os.name != "nt":
        raise ValueError(f"{field} is a Windows path on Linux and cannot be mapped safely: {token}")
    candidate = Path(token.replace("\\", os.sep))
    if not candidate.is_absolute():
        candidate = root / candidate
    resolved = candidate.resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"{field} is outside the package root: {resolved}") from exc
    return resolved


def legacy_goat_bundle_hash(directory: Path) -> str:
    digest = hashlib.sha256()
    for name in GOAT_BUNDLE_NAMES:
        path = directory / name
        if not path.is_file():
            raise FileNotFoundError(path)
        digest.update(name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest().upper()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert a generated ORCA manifest to package-relative, hash-locked HPC paths."
    )
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True, help="Generated CSV manifest")
    parser.add_argument("--output", type=Path, required=True, help="Portable CSV manifest inside package")
    parser.add_argument("--stage", required=True)
    parser.add_argument("--stage-config", type=Path)
    args = parser.parse_args()

    root = args.package_root.resolve()
    input_manifest = args.input.resolve()
    output_manifest = args.output
    if not output_manifest.is_absolute():
        output_manifest = root / output_manifest
    output_manifest = output_manifest.resolve()
    try:
        input_manifest.relative_to(root)
        output_manifest.relative_to(root)
    except ValueError as exc:
        raise SystemExit("Both input and output manifests must be inside the package root") from exc
    if input_manifest == output_manifest:
        raise SystemExit("Refusing in-place normalization; choose a new output path")

    resources_path = (args.stage_config or root / "hpc/config/stage_resources.json").resolve()
    stages = load_json(resources_path).get("stages", {})
    if args.stage not in stages:
        raise SystemExit(f"Unknown stage: {args.stage}")
    input_default = str(stages[args.stage]["input_default"])

    original_fields, rows = read_manifest(input_manifest)
    additions = [
        "stage",
        "portable_job_directory",
        "relative_job_directory",
        "job_directory",
        "input_file",
        "coordinate_file",
        "output_file",
        "expected_input_sha256",
        "expected_coordinate_sha256",
        "expected_input_bundle_sha256",
        "enabled",
    ]
    fields = list(original_fields)
    for field in additions:
        if field not in fields:
            fields.append(field)

    normalized: list[dict[str, str]] = []
    seen_job_ids: set[str] = set()
    for index, source_row in enumerate(rows, start=1):
        row = dict(source_row)
        job_id = row_value(row, "job_id", "model_id") or f"row_{index:05d}"
        if job_id in seen_job_ids:
            raise SystemExit(f"Duplicate job_id: {job_id}")
        seen_job_ids.add(job_id)
        declared_stage = row_value(row, "stage")
        if declared_stage and declared_stage != args.stage:
            raise SystemExit(
                f"Row {index} declares stage={declared_stage!r}, expected {args.stage!r}"
            )
        raw_job = row_value(
            row,
            "portable_job_directory",
            "relative_job_directory",
            "working_directory",
            "job_directory",
        )
        directory = inside_root(root, raw_job, f"row {index} job_directory")
        if not directory.is_dir():
            raise SystemExit(f"Row {index} job directory is missing: {directory}")
        relative_directory = str(directory.relative_to(root)).replace("\\", "/")

        raw_input = row_value(row, "portable_input_file", "input_file") or input_default
        input_candidate = Path(raw_input.replace("\\", os.sep))
        if input_candidate.is_absolute() or is_windows_absolute(raw_input):
            input_path = inside_root(root, raw_input, f"row {index} input_file")
        else:
            input_path = (directory / input_candidate).resolve()
        try:
            input_relative_to_job = input_path.relative_to(directory)
        except ValueError as exc:
            raise SystemExit(f"Row {index} input file is not inside its job directory: {input_path}") from exc
        if not input_path.is_file():
            raise SystemExit(f"Row {index} input file is missing: {input_path}")

        coordinate_raw = row_value(row, "coordinate_file")
        coordinate_path: Path | None = None
        if coordinate_raw:
            coordinate_candidate = Path(coordinate_raw.replace("\\", os.sep))
            if coordinate_candidate.is_absolute() or is_windows_absolute(coordinate_raw):
                coordinate_path = inside_root(root, coordinate_raw, f"row {index} coordinate_file")
            else:
                coordinate_path = (directory / coordinate_candidate).resolve()
        else:
            references = input_coordinate_references(input_path)
            if len(references) == 1:
                coordinate_path = references[0]
        coordinate_relative = ""
        coordinate_hash = ""
        if coordinate_path is not None:
            try:
                coordinate_relative = str(coordinate_path.relative_to(directory)).replace("\\", "/")
            except ValueError as exc:
                raise SystemExit(
                    f"Row {index} coordinate file is not inside its job directory: {coordinate_path}"
                ) from exc
            if not coordinate_path.is_file():
                raise SystemExit(f"Row {index} coordinate file is missing: {coordinate_path}")
            coordinate_hash = sha256_file(coordinate_path)

        output_name = Path(row_value(row, "output_file") or input_path.with_suffix(".out").name).name
        row.update(
            {
                "stage": args.stage,
                "portable_job_directory": relative_directory,
                "relative_job_directory": relative_directory,
                "job_directory": relative_directory,
                "input_file": str(input_relative_to_job).replace("\\", "/"),
                "coordinate_file": coordinate_relative,
                "output_file": output_name,
                "expected_input_sha256": sha256_file(input_path),
                "expected_coordinate_sha256": coordinate_hash,
                "enabled": row_value(row, "enabled") or "false",
            }
        )
        if all((directory / name).is_file() for name in GOAT_BUNDLE_NAMES):
            row["expected_input_bundle_sha256"] = legacy_goat_bundle_hash(directory)
        normalized.append(row)

    write_csv(output_manifest, fields, normalized)
    audit = {
        "schema_version": "1.0",
        "notice": "CALCULATED MODEL WORKFLOW - NOT EXPERIMENTAL DATA",
        "timestamp": utc_now(),
        "stage": args.stage,
        "row_count": len(normalized),
        "input_manifest": str(input_manifest.relative_to(root)).replace("\\", "/"),
        "input_manifest_sha256": sha256_file(input_manifest),
        "output_manifest": str(output_manifest.relative_to(root)).replace("\\", "/"),
        "output_manifest_sha256": sha256_file(output_manifest),
        "policy": "All job/input/coordinate paths are inside package_root and output paths are package-relative; hashes were recomputed.",
    }
    atomic_json(output_manifest.with_suffix(output_manifest.suffix + ".normalization.json"), audit)
    print(f"Normalized {len(normalized)} rows -> {output_manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
