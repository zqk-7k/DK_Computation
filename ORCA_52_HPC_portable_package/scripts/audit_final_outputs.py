#!/usr/bin/env python3
"""Read-only integrity and completeness audit for the 51-compound workflow.

The script never launches ORCA and never changes a production calculation,
manifest, status file, prediction bundle, spectrum, or HRMS result.  Its only
writes are three audit reports (JSON, CSV and Markdown) under ``--output-dir``;
use ``--stdout-only`` to suppress even those report writes.

Exit codes
----------
0  COMPLETE: every strict gate passed.
2  INCOMPLETE: one or more required jobs/results are not available yet.
3  INVALID: an integrity, provenance, strict-method, or semantic gate failed.
4  AUDITOR_ERROR: the auditor itself raised an unexpected exception.

An explicitly labelled external H/C regression is allowed only as a secondary
view.  Raw shifts from same-level calculated reference compounds must be
present and are treated as the primary formal baseline.  Partial property
coverage and every provisional imaginary-frequency/incomplete-DFT override
are rejected.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCRIPT = Path(__file__).resolve()
DEFAULT_PROJECT = SCRIPT.parents[1]
EXPECTED_IDS = tuple(f"F3_{index:03d}" for index in range(1, 53) if index != 5)
EXPECTED_ID_SET = set(EXPECTED_IDS)
EXPECTED_GOAT_JOBS = 103
NO_CARBON_IDS = {"F3_035", "F3_043", "F3_044"}
VALID_HASH = re.compile(r"^[0-9a-f]{64}$", re.IGNORECASE)
FIGURE_ID = re.compile(r"^F3_\d{3}$")
XYZFILE_PATTERN = re.compile(
    r"(?im)^\s*\*\s*(?:xyzfile|gzmtfile)\s+[-+]?\d+\s+\d+\s+"
    r'(?:"([^"]+)"|\'([^\']+)\'|(\S+))'
)
GOAT_XYZFILE_PATTERN = re.compile(
    r"^\s*\*\s*xyzfile\s+[+-]?\d+\s+\d+\s+"
    r'(?:"([^"]+)"|\'([^\']+)\'|(\S+))',
    flags=re.IGNORECASE,
)
ORCA_NORMAL_MARKERS = ("****ORCA TERMINATED NORMALLY****", "ORCA TERMINATED NORMALLY")
COMPLETE_ORCA_STATUSES = {"complete", "already_complete"}
CALCULATED_NOTICE_WORDS = ("calculated", "not experimental")
FLOAT_TOL = 2.0e-8


@dataclass
class Check:
    check_id: str
    section: str
    subject: str
    status: str
    message: str
    expected: Any = None
    observed: Any = None
    paths: list[str] | None = None


class Auditor:
    def __init__(self, project: Path):
        self.project = project.resolve()
        self.checks: list[Check] = []
        self.hash_cache: dict[tuple[str, int, int], str] = {}

    def add(
        self,
        check_id: str,
        section: str,
        subject: str,
        status: str,
        message: str,
        *,
        expected: Any = None,
        observed: Any = None,
        paths: Iterable[Path | str] = (),
    ) -> None:
        if status not in {"PASS", "INFO", "INCOMPLETE", "INVALID"}:
            raise ValueError(f"Unsupported audit status: {status}")
        self.checks.append(
            Check(
                check_id=check_id,
                section=section,
                subject=subject,
                status=status,
                message=message,
                expected=expected,
                observed=observed,
                paths=[str(value) for value in paths] or None,
            )
        )

    def pass_(self, check_id: str, section: str, subject: str, message: str, **kwargs: Any) -> None:
        self.add(check_id, section, subject, "PASS", message, **kwargs)

    def info(self, check_id: str, section: str, subject: str, message: str, **kwargs: Any) -> None:
        self.add(check_id, section, subject, "INFO", message, **kwargs)

    def incomplete(
        self, check_id: str, section: str, subject: str, message: str, **kwargs: Any
    ) -> None:
        self.add(check_id, section, subject, "INCOMPLETE", message, **kwargs)

    def invalid(self, check_id: str, section: str, subject: str, message: str, **kwargs: Any) -> None:
        self.add(check_id, section, subject, "INVALID", message, **kwargs)

    def sha256(self, path: Path) -> str:
        resolved = path.resolve()
        stat = resolved.stat()
        key = (os.path.normcase(str(resolved)), stat.st_size, stat.st_mtime_ns)
        cached = self.hash_cache.get(key)
        if cached is not None:
            return cached
        digest = hashlib.sha256()
        with resolved.open("rb") as handle:
            for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
                digest.update(block)
        value = digest.hexdigest()
        self.hash_cache[key] = value
        return value

    def summary(self) -> dict[str, Any]:
        counts = Counter(check.status for check in self.checks)
        if counts["INVALID"]:
            overall, exit_code = "INVALID", 3
        elif counts["INCOMPLETE"]:
            overall, exit_code = "INCOMPLETE", 2
        else:
            overall, exit_code = "COMPLETE", 0
        section_counts: dict[str, dict[str, int]] = {}
        for section in sorted({check.section for check in self.checks}):
            values = Counter(check.status for check in self.checks if check.section == section)
            section_counts[section] = {
                name: values[name] for name in ("PASS", "INFO", "INCOMPLETE", "INVALID")
            }
        return {
            "overall_status": overall,
            "exit_code": exit_code,
            "counts": {name: counts[name] for name in ("PASS", "INFO", "INCOMPLETE", "INVALID")},
            "sections": section_counts,
        }


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path: Path, comments: bool = False) -> list[dict[str, str]]:
    text = path.read_text(encoding="utf-8-sig")
    if comments:
        text = "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("#"))
    reader = csv.DictReader(io.StringIO(text))
    if reader.fieldnames is None:
        raise ValueError("CSV lacks a header")
    return [dict(row) for row in reader]


def as_bool(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)
    text = str(value).strip().casefold()
    if text in {"1", "true", "yes", "y"}:
        return True
    if text in {"0", "false", "no", "n", ""}:
        return False
    return None


def as_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def close_to(value: Any, expected: float, tolerance: float = FLOAT_TOL) -> bool:
    number = as_float(value)
    return number is not None and abs(number - expected) <= tolerance


def has_calculated_notice(value: Any) -> bool:
    text = str(value).casefold().replace("—", "-").replace("–", "-")
    return all(word in text for word in CALCULATED_NOTICE_WORDS)


def path_is_within(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except (OSError, ValueError):
        return False


def paths_equal(left: Path | str, right: Path | str) -> bool:
    try:
        return os.path.normcase(str(Path(left).resolve())) == os.path.normcase(str(Path(right).resolve()))
    except OSError:
        return False


def valid_hash(value: Any) -> bool:
    return bool(VALID_HASH.fullmatch(str(value).strip()))


def formula_elements(formula: str) -> Counter[str]:
    clean = re.sub(r"[+\-−·.]", "", str(formula))
    result: Counter[str] = Counter()
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", clean):
        result[element] += int(count or "1")
    return result


def file_has_orca_normal_marker(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    with path.open("rb") as handle:
        handle.seek(max(0, path.stat().st_size - 2 * 1024 * 1024))
        text = handle.read().decode("utf-8", errors="replace")
    return any(marker in text for marker in ORCA_NORMAL_MARKERS)


def resolve_referenced_coordinates(input_path: Path, *, goat: bool = False) -> list[tuple[str, Path]]:
    pattern = GOAT_XYZFILE_PATTERN if goat else XYZFILE_PATTERN
    text = input_path.read_text(encoding="utf-8", errors="replace")
    records: list[tuple[str, Path]] = []
    seen: set[str] = set()
    lines = text.splitlines() if goat else [text]
    for block in lines:
        matches = [pattern.match(block)] if goat else list(pattern.finditer(block))
        for match in matches:
            if match is None:
                continue
            token = next(value for value in match.groups() if value is not None).strip()
            candidate = Path(token)
            if not candidate.is_absolute():
                candidate = input_path.parent / candidate
            resolved = candidate.resolve()
            key = os.path.normcase(str(resolved))
            if key not in seen:
                records.append((token, resolved))
                seen.add(key)
    return records


def goat_bundle_hash(auditor: Auditor, job_directory: Path, input_name: str = "goat.inp") -> tuple[str, list[Path]]:
    root = job_directory.resolve(strict=True)
    input_path = (root / input_name).resolve(strict=True)
    if not path_is_within(input_path, root):
        raise ValueError("GOAT input escapes its job directory")
    referenced = resolve_referenced_coordinates(input_path, goat=True)
    files = [input_path]
    for _token, path in referenced:
        if not path.is_file() or not path_is_within(path, root):
            raise ValueError(f"Invalid GOAT coordinate reference: {path}")
        files.append(path)
    unique: dict[str, Path] = {}
    for path in files:
        unique.setdefault(os.path.normcase(str(path.resolve().relative_to(root))), path)
    files = list(unique.values())
    digest = hashlib.sha256()
    for path in sorted(files, key=lambda item: os.path.normcase(str(item.resolve()))):
        relative = path.resolve().relative_to(root).as_posix().encode("utf-8")
        digest.update(len(relative).to_bytes(4, "big"))
        digest.update(relative)
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest(), files


def orca_bundle(auditor: Auditor, input_path: Path) -> dict[str, Any]:
    resolved_input = input_path.resolve(strict=True)
    records: list[dict[str, Any]] = [
        {"role": "input", "logical_name": input_path.name, "path": resolved_input}
    ]
    for token, path in resolve_referenced_coordinates(resolved_input):
        records.append({"role": "coordinate", "logical_name": token.replace("\\", "/"), "path": path})
    identities: list[dict[str, Any]] = []
    for record in records:
        path = Path(record["path"])
        if not path.is_file():
            raise FileNotFoundError(path)
        identity = {
            "role": record["role"],
            "logical_name": record["logical_name"],
            "size_bytes": path.stat().st_size,
            "sha256": auditor.sha256(path),
        }
        identities.append(identity)
        record.update(identity)
        record["path"] = str(path)
    bundle_hash = hashlib.sha256(
        json.dumps(identities, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    ).hexdigest()
    return {"sha256": bundle_hash, "files": records}


def load_required(
    auditor: Auditor,
    path: Path,
    section: str,
    subject: str,
    kind: str,
    *,
    comments: bool = False,
) -> Any | None:
    if not path.is_file():
        auditor.incomplete(
            f"{section}.{subject}.missing",
            section,
            subject,
            f"Required {kind} is not available yet",
            paths=[path],
        )
        return None
    try:
        return read_csv(path, comments=comments) if kind == "CSV" else read_json(path)
    except (OSError, UnicodeError, csv.Error, json.JSONDecodeError, ValueError) as exc:
        auditor.invalid(
            f"{section}.{subject}.parse",
            section,
            subject,
            f"Required {kind} is malformed: {type(exc).__name__}: {exc}",
            paths=[path],
        )
        return None


def audit_id_set(
    auditor: Auditor,
    section: str,
    subject: str,
    observed: Iterable[str],
    *,
    missing_is_incomplete: bool = True,
) -> bool:
    values = set(observed)
    missing = sorted(EXPECTED_ID_SET - values)
    extra = sorted(values - EXPECTED_ID_SET)
    if extra:
        auditor.invalid(
            f"{section}.{subject}.unexpected_ids",
            section,
            subject,
            "Unexpected Figure 3 IDs are present; F3_005 is forbidden",
            expected=list(EXPECTED_IDS),
            observed={"extra": extra},
        )
    if missing:
        method = auditor.incomplete if missing_is_incomplete else auditor.invalid
        method(
            f"{section}.{subject}.missing_ids",
            section,
            subject,
            "The canonical 51-compound set is not complete",
            expected=list(EXPECTED_IDS),
            observed={"missing": missing, "count": len(values)},
        )
    if not missing and not extra and len(values) == 51:
        auditor.pass_(
            f"{section}.{subject}.ids",
            section,
            subject,
            "Exact 51-ID set verified and F3_005 is absent",
            expected=51,
            observed=51,
        )
        return True
    return False


def audit_goat(auditor: Auditor) -> tuple[bool, dict[str, dict[str, str]]]:
    section = "GOAT"
    project = auditor.project
    manifest_path = project / "data" / "GOAT_job_manifest.csv"
    rows = load_required(auditor, manifest_path, section, "manifest", "CSV")
    if rows is None:
        return False, {}
    if not isinstance(rows, list):
        auditor.invalid("GOAT.manifest.type", section, "manifest", "GOAT manifest is not tabular")
        return False, {}
    job_ids = [str(row.get("job_id", "")).strip() for row in rows]
    figure_ids = [str(row.get("figure_id", "")).strip() for row in rows]
    full_ids = audit_id_set(auditor, section, "manifest", figure_ids)
    if len(rows) < EXPECTED_GOAT_JOBS:
        auditor.incomplete(
            "GOAT.manifest.row_count",
            section,
            "manifest",
            "GOAT manifest has not reached the canonical 103 jobs",
            expected=EXPECTED_GOAT_JOBS,
            observed=len(rows),
        )
    elif len(rows) > EXPECTED_GOAT_JOBS:
        auditor.invalid(
            "GOAT.manifest.row_count",
            section,
            "manifest",
            "GOAT manifest contains unexpected additional jobs",
            expected=EXPECTED_GOAT_JOBS,
            observed=len(rows),
        )
    else:
        auditor.pass_(
            "GOAT.manifest.row_count",
            section,
            "manifest",
            "Canonical 103-job GOAT manifest row count verified",
            expected=EXPECTED_GOAT_JOBS,
            observed=len(rows),
        )
    if not all(job_ids) or len(job_ids) != len(set(job_ids)):
        auditor.invalid(
            "GOAT.manifest.unique_jobs",
            section,
            "manifest",
            "GOAT job IDs are empty or non-unique",
            observed={"rows": len(job_ids), "unique": len(set(job_ids)), "empty": job_ids.count("")},
        )
    else:
        auditor.pass_(
            "GOAT.manifest.unique_jobs",
            section,
            "manifest",
            "All GOAT job IDs are non-empty and unique",
            observed=len(job_ids),
        )
    pairs = [(str(row.get("figure_id", "")), str(row.get("replica", ""))) for row in rows]
    if len(pairs) != len(set(pairs)):
        auditor.invalid("GOAT.manifest.replica_pairs", section, "manifest", "Duplicate (figure_id, replica) rows found")
    counts = Counter(figure_ids)
    for figure_id, count in sorted(counts.items()):
        declared = {str(row.get("expected_replica_count", "")).strip() for row in rows if row.get("figure_id") == figure_id}
        if declared != {str(count)}:
            auditor.invalid(
                "GOAT.manifest.replica_count",
                section,
                figure_id,
                "Declared replica count does not equal manifest rows",
                expected=count,
                observed=sorted(declared),
            )

    current_hashes: dict[str, str] = {}
    manifest_by_job: dict[str, dict[str, str]] = {}
    input_valid = True
    for row in rows:
        job_id = str(row.get("job_id", "")).strip()
        figure_id = str(row.get("figure_id", "")).strip()
        if not job_id:
            continue
        manifest_by_job[job_id] = row
        job_dir = Path(str(row.get("job_directory", "")))
        if not path_is_within(job_dir, project / "01_goat"):
            auditor.invalid(
                "GOAT.input.job_directory",
                section,
                job_id,
                "GOAT job directory is outside the canonical 01_goat root",
                paths=[job_dir],
            )
            input_valid = False
            continue
        try:
            bundle_hash, files = goat_bundle_hash(auditor, job_dir)
        except (OSError, ValueError) as exc:
            auditor.incomplete(
                "GOAT.input.bundle",
                section,
                job_id,
                f"Current GOAT input bundle cannot be verified: {type(exc).__name__}: {exc}",
                paths=[job_dir],
            )
            input_valid = False
            continue
        current_hashes[job_id] = bundle_hash
        expected_input = str(row.get("goat_input_sha256", "")).strip().casefold()
        expected_xyz = str(row.get("start_xyz_sha256", "")).strip().casefold()
        input_path = job_dir / "goat.inp"
        xyz_path = job_dir / "start.xyz"
        if not valid_hash(expected_input) or not valid_hash(expected_xyz):
            auditor.incomplete(
                "GOAT.input.manifest_hash_fields",
                section,
                job_id,
                "GOAT manifest has not recorded both current input hashes",
                paths=[manifest_path],
            )
            input_valid = False
        else:
            actual_input = auditor.sha256(input_path)
            actual_xyz = auditor.sha256(xyz_path)
            if expected_input != actual_input or expected_xyz != actual_xyz:
                auditor.invalid(
                    "GOAT.input.manifest_hash_mismatch",
                    section,
                    job_id,
                    "GOAT manifest hash does not match a current source input",
                    expected={"goat.inp": expected_input, "start.xyz": expected_xyz},
                    observed={"goat.inp": actual_input, "start.xyz": actual_xyz},
                    paths=files,
                )
                input_valid = False

    state_root = project / "logs" / "goat_batch" / "jobs"
    state_files = list(state_root.glob("*.json")) if state_root.is_dir() else []
    state_ids = {path.stem for path in state_files}
    manifest_jobs = set(job_ids) - {""}
    missing_states = sorted(manifest_jobs - state_ids)
    extra_states = sorted(state_ids - manifest_jobs)
    if missing_states:
        auditor.incomplete(
            "GOAT.status.coverage",
            section,
            "job status files",
            "GOAT status files do not yet cover the manifest",
            expected=len(manifest_jobs),
            observed={"available": len(state_ids), "missing": missing_states},
            paths=[state_root],
        )
    if extra_states:
        auditor.invalid(
            "GOAT.status.extra",
            section,
            "job status files",
            "GOAT status root contains jobs outside the canonical manifest",
            observed=extra_states,
            paths=[state_root],
        )
    state_valid = not missing_states and not extra_states
    completed_current = 0
    for job_id in sorted(manifest_jobs & state_ids):
        path = state_root / f"{job_id}.json"
        try:
            state = read_json(path)
        except (OSError, json.JSONDecodeError, UnicodeError) as exc:
            auditor.invalid(
                "GOAT.status.parse",
                section,
                job_id,
                f"GOAT status JSON is malformed: {exc}",
                paths=[path],
            )
            state_valid = False
            continue
        if not isinstance(state, dict):
            auditor.invalid("GOAT.status.type", section, job_id, "GOAT status is not an object", paths=[path])
            state_valid = False
            continue
        current_hash = current_hashes.get(job_id, "")
        status = str(state.get("status", ""))
        if status != "completed":
            auditor.incomplete(
                "GOAT.status.not_complete",
                section,
                job_id,
                "GOAT job has not completed for the current source bundle",
                expected="completed",
                observed=status or "missing status",
                paths=[path],
            )
            state_valid = False
            continue
        attempts = state.get("attempts")
        latest = attempts[-1] if isinstance(attempts, list) and attempts and isinstance(attempts[-1], dict) else None
        errors: list[str] = []
        if str(state.get("job_id", job_id)) != job_id:
            errors.append("state job_id mismatch")
        if str(state.get("figure_id", manifest_by_job[job_id].get("figure_id", ""))) != manifest_by_job[job_id].get("figure_id", ""):
            errors.append("state figure_id mismatch")
        if not current_hash or str(state.get("input_bundle_sha256", "")).casefold() != current_hash.casefold():
            errors.append("state bundle hash differs from current source input")
        if latest is None:
            errors.append("latest attempt is absent")
        else:
            if latest.get("status") != "completed" or latest.get("return_code") != 0:
                errors.append("latest attempt is not a successful completion")
            if str(latest.get("input_bundle_sha256", "")).casefold() != current_hash.casefold():
                errors.append("latest attempt bundle hash differs from current source input")
            if latest.get("source_input_unchanged") is not True:
                errors.append("latest attempt did not attest unchanged source input")
            if str(latest.get("source_input_bundle_sha256_after", "")).casefold() != current_hash.casefold():
                errors.append("post-run source bundle hash differs from current input")
            output = Path(str(latest.get("output_file", "")))
            if not path_is_within(output, project / "02_goat_runs") or not file_has_orca_normal_marker(output):
                errors.append("GOAT output is missing, outside run root, or lacks normal termination")
            elif str(latest.get("output_sha256", "")).casefold() != auditor.sha256(output):
                errors.append("GOAT output hash mismatch")
            artifacts = latest.get("artifacts") if isinstance(latest.get("artifacts"), dict) else {}
            for name in ("globalminimum", "finalensemble"):
                artifact = Path(str(artifacts.get(name, "")))
                if not path_is_within(artifact, project / "02_goat_runs") or not artifact.is_file() or artifact.stat().st_size == 0:
                    errors.append(f"missing/noncanonical {name} artifact")
        if errors:
            auditor.invalid(
                "GOAT.status.current_hash_verification",
                section,
                job_id,
                "; ".join(errors),
                expected=current_hash,
                observed=state.get("input_bundle_sha256"),
                paths=[path],
            )
            state_valid = False
        else:
            completed_current += 1

    if completed_current == EXPECTED_GOAT_JOBS and len(rows) == EXPECTED_GOAT_JOBS and input_valid:
        auditor.pass_(
            "GOAT.status.current_hash_verification",
            section,
            "all jobs",
            "All 103 GOAT jobs completed against their current verified input hashes",
            expected=EXPECTED_GOAT_JOBS,
            observed=completed_current,
        )
    else:
        auditor.info(
            "GOAT.status.progress",
            section,
            "all jobs",
            "Current-hash-verified GOAT completion progress",
            expected=EXPECTED_GOAT_JOBS,
            observed=completed_current,
        )

    aggregate_path = project / "data" / "GOAT_batch_status.json"
    if completed_current == EXPECTED_GOAT_JOBS:
        aggregate = load_required(auditor, aggregate_path, section, "aggregate status", "JSON")
        if isinstance(aggregate, dict):
            jobs = aggregate.get("jobs")
            aggregate_ids = {str(item.get("job_id", "")) for item in jobs if isinstance(item, dict)} if isinstance(jobs, list) else set()
            if aggregate_ids != manifest_jobs or len(jobs or []) != EXPECTED_GOAT_JOBS:
                auditor.invalid(
                    "GOAT.aggregate.coverage",
                    section,
                    "aggregate status",
                    "Aggregate GOAT status does not exactly cover the canonical manifest",
                    expected=EXPECTED_GOAT_JOBS,
                    observed=len(aggregate_ids),
                    paths=[aggregate_path],
                )
            elif any(
                str(item.get("status")) != "completed"
                or str(item.get("input_bundle_sha256", "")).casefold() != current_hashes.get(str(item.get("job_id")), "").casefold()
                for item in jobs
            ):
                auditor.invalid(
                    "GOAT.aggregate.hash_status",
                    section,
                    "aggregate status",
                    "Aggregate GOAT status contains a non-complete or non-current row",
                    paths=[aggregate_path],
                )
            else:
                auditor.pass_("GOAT.aggregate.hash_status", section, "aggregate status", "Aggregate status exactly matches all current completed jobs", paths=[aggregate_path])

    complete = (
        full_ids
        and len(rows) == EXPECTED_GOAT_JOBS
        and len(job_ids) == len(set(job_ids))
        and input_valid
        and state_valid
        and completed_current == EXPECTED_GOAT_JOBS
    )
    return complete, manifest_by_job


def required_artifact_names(row: Mapping[str, str], stage: str) -> list[str]:
    raw = str(row.get("required_artifacts", row.get("required_artifact_files", ""))).strip()
    if raw:
        if raw.startswith("["):
            values = json.loads(raw)
            if not isinstance(values, list):
                raise ValueError("required_artifacts JSON is not a list")
            return [str(value).strip() for value in values if str(value).strip()]
        return [value.strip() for value in re.split(r"[;,|]", raw) if value.strip()]
    output_name = str(row.get("output_file", ""))
    stem = Path(output_name).stem
    context = " ".join((stage, str(row.get("job_id", "")), output_name)).casefold()
    if "dft_opt_freq" in context or ("opt" in context and "freq" in context):
        return [f"{stem}.gbw", f"{stem}.hess", f"{stem}.property.txt", f"{stem}.xyz"]
    if "shield" in context or "coupling" in context:
        return [f"{stem}.gbw", f"{stem}.property.txt"]
    return []


def audit_orca_stage(
    auditor: Auditor,
    *,
    section: str,
    manifest_path: Path,
    status_root: Path,
    expected_stage: str,
    job_root: Path,
    suffix: str | None,
) -> tuple[bool, list[dict[str, str]], dict[str, dict[str, Any]]]:
    rows = load_required(auditor, manifest_path, section, "manifest", "CSV")
    if rows is None:
        return False, [], {}
    if not isinstance(rows, list) or not rows:
        auditor.incomplete(f"{section}.manifest.empty", section, "manifest", "Manifest is empty", paths=[manifest_path])
        return False, [], {}
    job_ids = [str(row.get("job_id", "")).strip() for row in rows]
    figure_ids = [str(row.get("figure_id", "")).strip() for row in rows]
    full_ids = audit_id_set(auditor, section, "manifest", figure_ids)
    unique = all(job_ids) and len(job_ids) == len(set(job_ids))
    if not unique:
        auditor.invalid(
            f"{section}.manifest.unique_jobs",
            section,
            "manifest",
            "Manifest job IDs are empty or non-unique",
            observed={"rows": len(job_ids), "unique": len(set(job_ids))},
            paths=[manifest_path],
        )
    else:
        auditor.pass_(f"{section}.manifest.unique_jobs", section, "manifest", "Manifest job IDs are non-empty and unique", observed=len(job_ids))

    status_files = list(status_root.glob("*.json")) if status_root.is_dir() else []
    status_files = [path for path in status_files if path.name != "aggregate_status.json"]
    state_ids = {path.stem for path in status_files}
    manifest_ids = set(job_ids) - {""}
    missing = sorted(manifest_ids - state_ids)
    extra = sorted(state_ids - manifest_ids)
    if missing:
        auditor.incomplete(
            f"{section}.status.missing",
            section,
            "status coverage",
            "Status files do not yet cover every manifest job",
            expected=len(manifest_ids),
            observed={"available": len(state_ids), "missing": missing},
            paths=[status_root],
        )
    if extra:
        auditor.invalid(
            f"{section}.status.extra",
            section,
            "status coverage",
            "Status root contains jobs outside this canonical manifest",
            observed=extra,
            paths=[status_root],
        )
    if not missing and not extra:
        auditor.pass_(
            f"{section}.status.coverage",
            section,
            "status coverage",
            "Status-file set exactly covers the manifest",
            observed=len(state_ids),
        )

    status_by_job: dict[str, dict[str, Any]] = {}
    valid = full_ids and unique and not missing and not extra
    complete_count = 0
    for row in rows:
        job_id = str(row.get("job_id", "")).strip()
        figure_id = str(row.get("figure_id", "")).strip()
        if not job_id:
            continue
        if suffix and not job_id.endswith(suffix):
            auditor.invalid(f"{section}.manifest.job_suffix", section, job_id, f"Job ID does not end in {suffix!r}")
            valid = False
        job_dir = Path(str(row.get("job_directory", "")))
        if not path_is_within(job_dir, job_root):
            auditor.invalid(
                f"{section}.manifest.job_root",
                section,
                job_id,
                "Job directory is outside the canonical production root",
                expected=str(job_root),
                observed=str(job_dir),
            )
            valid = False
            continue
        input_name = str(row.get("input_file", "")).strip()
        output_name = str(row.get("output_file", "")).strip()
        input_path = job_dir / input_name
        output_path = job_dir / output_name
        try:
            bundle = orca_bundle(auditor, input_path)
        except (OSError, ValueError) as exc:
            auditor.incomplete(
                f"{section}.input.bundle",
                section,
                job_id,
                f"Input bundle cannot be verified yet: {type(exc).__name__}: {exc}",
                paths=[input_path],
            )
            valid = False
            continue
        expected_bundle = str(row.get("expected_input_bundle_sha256", row.get("input_bundle_sha256", ""))).strip().casefold()
        expected_input = str(row.get("expected_input_sha256", row.get("input_file_sha256", ""))).strip().casefold()
        coordinate_records = [item for item in bundle["files"] if item["role"] == "coordinate"]
        expected_coordinate = str(row.get("expected_coordinate_sha256", row.get("conformer_xyz_sha256", ""))).strip().casefold()
        hash_errors: list[str] = []
        if not valid_hash(expected_bundle):
            hash_errors.append("manifest lacks expected bundle hash")
        elif expected_bundle != bundle["sha256"]:
            hash_errors.append("manifest bundle hash mismatch")
        input_record = next(item for item in bundle["files"] if item["role"] == "input")
        if not valid_hash(expected_input):
            hash_errors.append("manifest lacks expected input hash")
        elif expected_input != input_record["sha256"]:
            hash_errors.append("manifest input hash mismatch")
        if expected_coordinate:
            if not valid_hash(expected_coordinate) or len(coordinate_records) != 1 or expected_coordinate != coordinate_records[0]["sha256"]:
                hash_errors.append("manifest coordinate hash mismatch")
        if hash_errors:
            # Missing hash fields indicate that a canonical manifest has not yet
            # been rebuilt; a recorded mismatch is an integrity failure.
            only_missing = all("lacks" in error for error in hash_errors)
            method = auditor.incomplete if only_missing else auditor.invalid
            method(
                f"{section}.manifest.current_hashes",
                section,
                job_id,
                "; ".join(hash_errors),
                expected={"bundle": expected_bundle, "input": expected_input, "coordinate": expected_coordinate},
                observed={
                    "bundle": bundle["sha256"],
                    "input": input_record["sha256"],
                    "coordinates": [record["sha256"] for record in coordinate_records],
                },
                paths=[input_path],
            )
            valid = False

        status_path = status_root / f"{job_id}.json"
        if not status_path.is_file():
            continue
        try:
            status = read_json(status_path)
        except (OSError, json.JSONDecodeError, UnicodeError) as exc:
            auditor.invalid(f"{section}.status.parse", section, job_id, f"Malformed status JSON: {exc}", paths=[status_path])
            valid = False
            continue
        if not isinstance(status, dict):
            auditor.invalid(f"{section}.status.type", section, job_id, "Status JSON is not an object", paths=[status_path])
            valid = False
            continue
        status_by_job[job_id] = status
        state = str(status.get("status", ""))
        if state == "manifest_hash_mismatch":
            auditor.invalid(f"{section}.status.manifest_hash_mismatch", section, job_id, "Runner detected a manifest hash mismatch", paths=[status_path])
            valid = False
            continue
        if state not in COMPLETE_ORCA_STATUSES:
            auditor.incomplete(
                f"{section}.status.not_complete",
                section,
                job_id,
                "ORCA job is not complete",
                expected=sorted(COMPLETE_ORCA_STATUSES),
                observed=state or "missing status",
                paths=[status_path],
            )
            valid = False
            continue
        errors: list[str] = []
        if str(status.get("job_id", "")) != job_id:
            errors.append("status job_id mismatch")
        if str(status.get("figure_id", "")) != figure_id:
            errors.append("status figure_id mismatch")
        if str(status.get("stage", "")) != expected_stage:
            errors.append("runner stage mismatch")
        if str(status.get("input_bundle_sha256", "")).casefold() != bundle["sha256"].casefold():
            errors.append("status input-bundle hash mismatch")
        nested = status.get("input_bundle")
        if not isinstance(nested, dict) or str(nested.get("sha256", "")).casefold() != bundle["sha256"].casefold():
            errors.append("status nested input-bundle hash mismatch")
        else:
            nested_files = {
                (str(item.get("role")), str(item.get("logical_name"))): str(item.get("sha256", "")).casefold()
                for item in nested.get("files", [])
                if isinstance(item, dict)
            }
            current_files = {(item["role"], item["logical_name"]): item["sha256"].casefold() for item in bundle["files"]}
            if nested_files != current_files:
                errors.append("status input-file hash records differ from current bundle")
        validation = status.get("manifest_hash_validation")
        if not isinstance(validation, dict) or validation.get("matches") is not True:
            errors.append("runner manifest-hash validation did not pass")
        elif not validation.get("checks"):
            errors.append("runner manifest-hash validation contains no explicit checks")
        elif any(not isinstance(check, dict) or check.get("matches") is not True for check in validation.get("checks", [])):
            errors.append("one or more runner manifest-hash checks failed")
        if status.get("normal_termination_marker") is not True or not file_has_orca_normal_marker(output_path):
            errors.append("normal termination marker is absent")
        if str(status.get("output", "")) and not paths_equal(str(status.get("output")), output_path):
            errors.append("status output path mismatch")
        if not output_path.is_file() or not valid_hash(status.get("output_sha256")) or str(status.get("output_sha256")).casefold() != auditor.sha256(output_path):
            errors.append("output file/hash mismatch")
        try:
            expected_artifacts = required_artifact_names(row, expected_stage)
        except (ValueError, json.JSONDecodeError) as exc:
            expected_artifacts = []
            errors.append(f"invalid required_artifacts field: {exc}")
        recorded_artifacts = {
            str(item.get("name", "")): item
            for item in status.get("required_artifacts", [])
            if isinstance(item, dict)
        }
        if set(recorded_artifacts) != set(expected_artifacts):
            errors.append("status artifact-name set differs from manifest")
        for name in expected_artifacts:
            artifact = job_dir / name
            record = recorded_artifacts.get(name, {})
            if not artifact.is_file() or artifact.stat().st_size == 0:
                errors.append(f"required artifact absent/empty: {name}")
            elif (
                not paths_equal(str(record.get("path", "")), artifact)
                or
                record.get("exists_nonempty") is not True
                or not valid_hash(record.get("sha256"))
                or str(record.get("sha256")).casefold() != auditor.sha256(artifact)
            ):
                errors.append(f"required artifact hash mismatch: {name}")
        if status.get("required_artifacts_complete") is not True:
            errors.append("required_artifacts_complete is not true")
        orca_path = Path(str(status.get("orca_executable", "")))
        if not orca_path.is_file() or not valid_hash(status.get("orca_executable_sha256")):
            errors.append("ORCA executable identity is incomplete")
        elif str(status.get("orca_executable_sha256")).casefold() != auditor.sha256(orca_path):
            errors.append("ORCA executable hash mismatch")
        if errors:
            auditor.invalid(
                f"{section}.status.provenance",
                section,
                job_id,
                "; ".join(errors),
                paths=[status_path, input_path, output_path],
            )
            valid = False
        else:
            complete_count += 1

    if complete_count == len(rows) and full_ids and unique:
        auditor.pass_(
            f"{section}.status.all_verified",
            section,
            "all jobs",
            "Every manifest job has exact status coverage and current verified hashes",
            expected=len(rows),
            observed=complete_count,
        )
    else:
        auditor.info(
            f"{section}.status.progress",
            section,
            "all jobs",
            "Fully verified job progress",
            expected=len(rows),
            observed=complete_count,
        )
    return valid and complete_count == len(rows), rows, status_by_job


def audit_dft_goat_provenance(
    auditor: Auditor,
    rows: list[dict[str, str]],
    goat_rows: Mapping[str, Mapping[str, str]],
) -> bool:
    section = "DFT"
    if not rows:
        return False
    goat_manifest = auditor.project / "data" / "GOAT_job_manifest.csv"
    if not goat_manifest.is_file():
        return False
    current_manifest_hash = auditor.sha256(goat_manifest)
    valid = True
    for row in rows:
        job_id = str(row.get("job_id", ""))
        source_job = str(row.get("source_goat_job", ""))
        source = goat_rows.get(source_job)
        errors: list[str] = []
        if source is None:
            errors.append("source GOAT job is absent from canonical manifest")
        source_manifest_path = Path(str(row.get("source_goat_manifest", "")))
        if not paths_equal(source_manifest_path, goat_manifest):
            errors.append("source GOAT manifest path is not canonical")
        if str(row.get("source_goat_verification_status", "")) != "VERIFIED_BATCH_STATUS_AND_INPUT_BUNDLE":
            errors.append("source GOAT verification is not strict/current-hash verified")
        if str(row.get("source_goat_manifest_sha256", "")).casefold() != current_manifest_hash:
            errors.append("source GOAT manifest hash mismatch")
        goat_state_path = auditor.project / "logs" / "goat_batch" / "jobs" / f"{source_job}.json"
        try:
            goat_state = read_json(goat_state_path)
        except (OSError, json.JSONDecodeError, UnicodeError):
            goat_state = {}
            errors.append("source GOAT completed-state JSON is absent or malformed")
        source_bundle_hash = str(row.get("source_goat_input_bundle_sha256", "")).casefold()
        if (
            not isinstance(goat_state, dict)
            or goat_state.get("status") != "completed"
            or not valid_hash(source_bundle_hash)
            or str(goat_state.get("input_bundle_sha256", "")).casefold() != source_bundle_hash
        ):
            errors.append("source GOAT input-bundle hash does not match completed state")
        source_output = Path(str(row.get("source_goat_output_file", "")))
        source_ensemble = Path(str(row.get("source_goat_finalensemble_file", "")))
        if not source_output.is_file() or str(row.get("source_goat_output_sha256", "")).casefold() != auditor.sha256(source_output):
            errors.append("source GOAT output hash mismatch")
        if not source_ensemble.is_file() or str(row.get("source_goat_finalensemble_sha256", "")).casefold() != auditor.sha256(source_ensemble):
            errors.append("source GOAT final ensemble hash mismatch")
        if errors:
            method = auditor.incomplete if any("absent" in error for error in errors) and not any("mismatch" in error for error in errors) else auditor.invalid
            method(
                "DFT.source_goat.provenance",
                section,
                job_id,
                "; ".join(errors),
                paths=[source_output, source_ensemble, goat_state_path],
            )
            valid = False
    if valid:
        auditor.pass_("DFT.source_goat.provenance", section, "all jobs", "Every DFT input has strict current-hash GOAT provenance", observed=len(rows))
    return valid


def audit_property_strictness(
    auditor: Auditor,
    dft_rows: list[dict[str, str]],
    shielding_rows: list[dict[str, str]],
    scalar_rows: list[dict[str, str]],
) -> bool:
    section = "PROPERTY_STRICTNESS"
    if not shielding_rows or not scalar_rows:
        return False
    shield_map = {str(row.get("conformer_id", "")): row for row in shielding_rows}
    scalar_map = {str(row.get("conformer_id", "")): row for row in scalar_rows}
    dft_map = {str(row.get("job_id", "")): row for row in dft_rows}
    if set(shield_map) != set(scalar_map) or len(shield_map) != len(shielding_rows) or len(scalar_map) != len(scalar_rows):
        auditor.invalid(
            "PROPERTY.manifest.pairing",
            section,
            "shielding/scalar manifests",
            "Shielding and scalar-coupling manifests do not have an exact one-to-one conformer set",
            observed={"shielding": len(shield_map), "scalar": len(scalar_map)},
        )
        return False
    valid = True
    for conformer_id in sorted(shield_map):
        dft_row = dft_map.get(conformer_id)
        if dft_row is None:
            auditor.invalid(
                "PROPERTY.source_dft.manifest",
                section,
                conformer_id,
                "Selected property conformer is absent from the canonical DFT manifest",
            )
            valid = False
            continue
        dft_job_dir = Path(str(dft_row.get("job_directory", "")))
        dft_output_name = str(dft_row.get("output_file", "dft_opt_freq.out"))
        dft_stem = Path(dft_output_name).stem
        source_paths = (
            (dft_job_dir / str(dft_row.get("input_file", "dft_opt_freq.inp")), "source_dft_input_sha256", "DFT input"),
            (dft_job_dir / "start.xyz", "source_dft_start_xyz_sha256", "DFT start geometry"),
            (dft_job_dir / dft_output_name, "source_dft_output_sha256", "DFT output"),
            (dft_job_dir / f"{dft_stem}.xyz", "source_dft_final_xyz_sha256", "DFT final geometry"),
            (dft_job_dir / f"{dft_stem}.property.txt", "source_dft_property_sha256", "DFT property file"),
        )
        for kind, row in (("shielding", shield_map[conformer_id]), ("scalar", scalar_map[conformer_id])):
            errors: list[str] = []
            if as_bool(row.get("canonical_strict_build")) is not True:
                errors.append("canonical_strict_build is not true")
            if str(row.get("property_build_tag", "")).casefold() != "canonical":
                errors.append("property_build_tag is not canonical")
            if as_bool(row.get("include_all_valid_minima")) is not True:
                errors.append("include_all_valid_minima is not true")
            if str(row.get("property_selection_policy", "")) != "all_valid_post_dft_unique_minima_strict":
                errors.append("property selection policy is not strict all-valid-minima")
            if not close_to(row.get("property_selected_population_of_all_valid_dft"), 1.0):
                errors.append("property selected population is below 100%")
            if as_bool(row.get("quasi_rrho")) is not True:
                errors.append("Quasi-RRHO is not true")
            if not close_to(row.get("thermochemistry_temperature_K"), 298.15, 0.011):
                errors.append("thermochemistry temperature is not 298.15 K")
            if not str(row.get("thermochemistry_validation", "")).startswith("passed_Quasi_RRHO_True"):
                errors.append("thermochemistry validation did not pass")
            if str(row.get("significant_imaginary_policy", "")) != "strict_pass":
                errors.append("significant-imaginary policy is not strict_pass")
            if str(row.get("provisional_override_audit", "")).strip():
                errors.append("a provisional override is recorded")
            if str(row.get("source_dft_runner_stage", "")) != "dft_production":
                errors.append("source DFT stage is not dft_production")
            if str(row.get("source_dft_runner_status", "")) not in COMPLETE_ORCA_STATUSES:
                errors.append("source DFT status is not complete")
            if not path_is_within(Path(str(row.get("job_directory", ""))), auditor.project / "03_nmr_properties"):
                errors.append("property job is outside canonical 03_nmr_properties")
            status_path = Path(str(row.get("source_dft_status_path", "")))
            if not path_is_within(status_path, auditor.project / "logs" / "dft_production"):
                errors.append("source DFT status is outside dft_production")
            elif not status_path.is_file() or str(row.get("source_dft_status_sha256", "")).casefold() != auditor.sha256(status_path):
                errors.append("source DFT status hash mismatch")
            for path, field_hash, label in source_paths:
                if not path.is_file() or not valid_hash(row.get(field_hash)) or str(row.get(field_hash)).casefold() != auditor.sha256(path):
                    errors.append(f"{label} hash mismatch")
            try:
                current_dft_bundle = orca_bundle(
                    auditor,
                    dft_job_dir / str(dft_row.get("input_file", "dft_opt_freq.inp")),
                )
            except (OSError, ValueError) as exc:
                errors.append(f"current DFT bundle cannot be verified: {exc}")
            else:
                if str(row.get("source_dft_input_bundle_sha256", "")).casefold() != current_dft_bundle["sha256"].casefold():
                    errors.append("source DFT input-bundle hash mismatch")
            conformer_path = Path(str(row.get("job_directory", ""))) / "conformer.xyz"
            if not conformer_path.is_file() or str(row.get("conformer_xyz_sha256", "")).casefold() != auditor.sha256(conformer_path):
                errors.append("copied conformer hash mismatch")
            elif str(row.get("conformer_xyz_sha256", "")).casefold() != str(
                row.get("source_dft_final_xyz_sha256", "")
            ).casefold():
                errors.append("copied conformer is not identical to the verified DFT final geometry")
            if errors:
                auditor.invalid(
                    "PROPERTY.strictness",
                    section,
                    f"{conformer_id}:{kind}",
                    "; ".join(errors),
                    paths=[status_path],
                )
                valid = False
        shared = (
            "figure_id",
            "boltzmann_weight",
            "conformer_xyz_sha256",
            "source_dft_status_sha256",
            "source_dft_output_sha256",
            "source_dft_final_xyz_sha256",
        )
        if any(str(shield_map[conformer_id].get(field, "")) != str(scalar_map[conformer_id].get(field, "")) for field in shared):
            auditor.invalid("PROPERTY.manifest.shared_provenance", section, conformer_id, "Shielding/scalar rows disagree on shared conformer provenance")
            valid = False
    if valid:
        auditor.pass_(
            "PROPERTY.strictness",
            section,
            "all selected conformers",
            "Canonical, strict, Quasi-RRHO property preparation and DFT provenance verified",
            observed=len(shield_map),
        )
    return valid


def audit_boltzmann(
    auditor: Auditor,
    dft_rows: list[dict[str, str]],
    shielding_rows: list[dict[str, str]],
    scalar_rows: list[dict[str, str]],
) -> bool:
    section = "BOLTZMANN"
    audit_path = auditor.project / "data" / "DFT_minimum_and_weight_audit.csv"
    audit_rows = load_required(auditor, audit_path, section, "DFT minimum/weight audit", "CSV")
    if audit_rows is None or not isinstance(audit_rows, list):
        return False
    dft_ids = {str(row.get("job_id", "")) for row in dft_rows}
    audit_ids = [str(row.get("job_id", "")) for row in audit_rows]
    if set(audit_ids) != dft_ids or len(audit_ids) != len(set(audit_ids)):
        auditor.invalid(
            "BOLTZMANN.audit.coverage",
            section,
            "DFT audit",
            "DFT minimum/weight audit does not exactly cover the DFT manifest",
            expected=len(dft_ids),
            observed={"rows": len(audit_ids), "unique": len(set(audit_ids)), "missing": sorted(dft_ids - set(audit_ids)), "extra": sorted(set(audit_ids) - dft_ids)},
            paths=[audit_path],
        )
        return False
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    valid = True
    for row in audit_rows:
        grouped[str(row.get("figure_id", ""))].append(row)
        errors: list[str] = []
        if str(row.get("dft_provenance_validation", "")) != "passed_current_inp_start_output_and_artifact_hashes":
            errors.append("DFT provenance validation did not pass")
        if str(row.get("parse_status", "")) != "parsed":
            errors.append("DFT output is not fully parsed")
        if as_bool(row.get("significant_imaginary")) is not False:
            errors.append("significant imaginary frequency is present/unknown")
        if as_bool(row.get("quasi_rrho")) is not True:
            errors.append("Quasi-RRHO is not true")
        if not close_to(row.get("thermochemistry_temperature_K"), 298.15, 0.011):
            errors.append("thermochemistry temperature is not 298.15 K")
        if as_bool(row.get("canonical_strict_build")) is not True:
            errors.append("canonical_strict_build is not true")
        if str(row.get("property_build_tag", "")).casefold() != "canonical":
            errors.append("property build tag is not canonical")
        if as_bool(row.get("include_all_valid_minima")) is not True:
            errors.append("include_all_valid_minima is not true")
        if str(row.get("property_selection_policy", "")) != "all_valid_post_dft_unique_minima_strict":
            errors.append("property selection policy is not strict all-valid-minima")
        if as_bool(row.get("significant_imaginary_override_enabled")) is not False:
            errors.append("imaginary-frequency override is enabled")
        if as_bool(row.get("incomplete_dft_override_enabled")) is not False:
            errors.append("incomplete-DFT override is enabled")
        if str(row.get("significant_imaginary_policy", "")) != "strict_pass":
            errors.append("imaginary-frequency policy is not strict_pass")
        if str(row.get("provisional_override_audit", "")).strip():
            errors.append("provisional override audit is non-empty")
        if errors:
            auditor.invalid("BOLTZMANN.strict_dft", section, str(row.get("job_id", "")), "; ".join(errors))
            valid = False
    selected_ids: set[str] = set()
    for figure_id in EXPECTED_IDS:
        rows = grouped.get(figure_id, [])
        if not rows:
            auditor.incomplete("BOLTZMANN.compound.missing", section, figure_id, "No DFT minimum/weight rows are available")
            valid = False
            continue
        unique_rows = [row for row in rows if as_bool(row.get("dft_unique_minimum")) is True]
        duplicate_rows = [row for row in rows if as_bool(row.get("dft_unique_minimum")) is False]
        if not unique_rows:
            auditor.invalid("BOLTZMANN.unique_minimum", section, figure_id, "No unique DFT minimum remains after deduplication")
            valid = False
            continue
        all_weights = [as_float(row.get("boltzmann_weight_all_valid")) for row in unique_rows]
        selected = [row for row in unique_rows if as_bool(row.get("selected_for_nmr")) is True]
        selected_weights = [as_float(row.get("boltzmann_weight_selected_normalized")) for row in selected]
        if any(value is None or value < 0.0 for value in all_weights) or not close_to(sum(value or 0.0 for value in all_weights), 1.0):
            auditor.invalid(
                "BOLTZMANN.partition.all_valid",
                section,
                figure_id,
                "All-valid unique-minimum Boltzmann weights do not sum to one",
                expected=1.0,
                observed=sum(value or 0.0 for value in all_weights),
            )
            valid = False
        if any(value is None or value < 0.0 for value in selected_weights) or not close_to(sum(value or 0.0 for value in selected_weights), 1.0):
            auditor.invalid(
                "BOLTZMANN.partition.selected",
                section,
                figure_id,
                "Selected-normalized Boltzmann weights do not sum to one",
                expected=1.0,
                observed=sum(value or 0.0 for value in selected_weights),
            )
            valid = False
        selected_population = sum(as_float(row.get("boltzmann_weight_all_valid")) or 0.0 for row in selected)
        # Quality-first final delivery requires property calculations for every
        # valid unique minimum, rather than any deliberately truncated partition.
        if len(selected) != len(unique_rows) or not close_to(selected_population, 1.0):
            auditor.invalid(
                "BOLTZMANN.population.full_coverage",
                section,
                figure_id,
                "Property selection does not cover 100% of the valid unique-minimum population",
                expected={"selected_population": 1.0, "selected_minima": len(unique_rows)},
                observed={"selected_population": selected_population, "selected_minima": len(selected)},
            )
            valid = False
        for row in unique_rows:
            if not close_to(row.get("property_selected_population_of_all_valid_dft"), 1.0):
                auditor.invalid(
                    "BOLTZMANN.audit.recorded_coverage",
                    section,
                    str(row.get("job_id", "")),
                    "Recorded property population coverage is not 1.0",
                    expected=1.0,
                    observed=row.get("property_selected_population_of_all_valid_dft"),
                )
                valid = False
        selected_ids.update(str(row.get("job_id", "")) for row in selected)
        for row in duplicate_rows:
            if as_bool(row.get("selected_for_nmr")) is True or str(row.get("boltzmann_weight_all_valid", "")).strip():
                auditor.invalid("BOLTZMANN.duplicate_exclusion", section, str(row.get("job_id", "")), "A post-DFT duplicate entered the partition/property selection")
                valid = False
    shield_conformers = {str(row.get("conformer_id", "")) for row in shielding_rows}
    scalar_conformers = {str(row.get("conformer_id", "")) for row in scalar_rows}
    if selected_ids != shield_conformers or selected_ids != scalar_conformers:
        auditor.invalid(
            "BOLTZMANN.property_manifest_coverage",
            section,
            "selected conformers",
            "Property manifests do not exactly cover the selected unique minima",
            expected=sorted(selected_ids),
            observed={"shielding_missing": sorted(selected_ids - shield_conformers), "scalar_missing": sorted(selected_ids - scalar_conformers), "shielding_extra": sorted(shield_conformers - selected_ids), "scalar_extra": sorted(scalar_conformers - selected_ids)},
        )
        valid = False
    else:
        audit_map = {str(row.get("job_id", "")): row for row in audit_rows}
        for row in shielding_rows + scalar_rows:
            conformer_id = str(row.get("conformer_id", ""))
            expected_weight = as_float(audit_map[conformer_id].get("boltzmann_weight_selected_normalized"))
            if expected_weight is None or not close_to(row.get("boltzmann_weight"), expected_weight):
                auditor.invalid("BOLTZMANN.property_weight", section, str(row.get("job_id", "")), "Property-manifest weight differs from DFT audit")
                valid = False
    if valid:
        auditor.pass_(
            "BOLTZMANN.full_coverage",
            section,
            "all compounds",
            "All valid unique minima have normalized 100% shielding/scalar property coverage",
            observed={"compounds": 51, "selected_conformers": len(selected_ids)},
        )
    return valid


def verify_source_record(auditor: Auditor, record: Mapping[str, Any]) -> list[str]:
    errors: list[str] = []
    path = Path(str(record.get("path", "")))
    if record.get("exists") is not True or not path.is_file():
        errors.append(f"source missing: {path}")
    elif not valid_hash(record.get("sha256")) or str(record.get("sha256")).casefold() != auditor.sha256(path):
        errors.append(f"source hash mismatch: {path}")
    elif int(record.get("size_bytes") or -1) != path.stat().st_size:
        errors.append(f"source size mismatch: {path}")
    return errors


def audit_nmr_predictions(auditor: Auditor) -> tuple[bool, dict[str, dict[str, Any]]]:
    section = "NMR_PREDICTIONS"
    root = auditor.project / "NMR_predictions_CALCULATED_NOT_EXPERIMENTAL"
    summary_path = root / "processing_summary__ALL.json"
    summary = load_required(auditor, summary_path, section, "ALL summary", "JSON")
    directory_ids = {path.name for path in root.glob("F3_*") if path.is_dir()} if root.is_dir() else set()
    extras = directory_ids - EXPECTED_ID_SET
    if extras:
        auditor.invalid("NMR_PREDICTIONS.directories.extra", section, "output directories", "Unexpected/F3_005 prediction directories are present", observed=sorted(extras), paths=[root])
    if summary is None:
        # Pilot directories are intentionally ignored until a canonical ALL
        # summary exists; otherwise an in-progress production run could be
        # misclassified because of stale pilot artifacts.
        auditor.info("NMR_PREDICTIONS.directory_progress", section, "output directories", "Canonical bundle directory progress", expected=51, observed=len(directory_ids & EXPECTED_ID_SET))
        return False, {}
    if not isinstance(summary, dict):
        auditor.invalid("NMR_PREDICTIONS.summary.type", section, "ALL summary", "Summary JSON is not an object", paths=[summary_path])
        return False, {}
    selected = [str(value) for value in summary.get("selected_figure_ids", [])]
    results = summary.get("results")
    if not isinstance(results, list):
        auditor.invalid("NMR_PREDICTIONS.summary.results", section, "ALL summary", "Summary results is not a list", paths=[summary_path])
        return False, {}
    result_ids = [str(row.get("figure_id", "")) for row in results if isinstance(row, dict)]
    full_ids = audit_id_set(auditor, section, "summary", selected)
    if set(result_ids) != EXPECTED_ID_SET or len(result_ids) != 51 or len(set(result_ids)) != 51:
        auditor.invalid("NMR_PREDICTIONS.summary.result_ids", section, "ALL summary", "Summary result rows are not the exact unique 51-ID set", observed=result_ids, paths=[summary_path])
        full_ids = False
    if set(directory_ids) != EXPECTED_ID_SET:
        missing = sorted(EXPECTED_ID_SET - directory_ids)
        if missing:
            auditor.incomplete("NMR_PREDICTIONS.directories.missing", section, "output directories", "Prediction bundle directories are incomplete", observed={"missing": missing})
        full_ids = False
    if not has_calculated_notice(summary.get("notice")):
        auditor.invalid("NMR_PREDICTIONS.summary.notice", section, "ALL summary", "Summary is not explicitly labelled calculated/not experimental")
        full_ids = False
    valid = full_ids
    for row in results:
        if not isinstance(row, dict):
            auditor.invalid("NMR_PREDICTIONS.summary.row", section, "ALL summary", "Non-object result row found")
            valid = False
            continue
        figure_id = str(row.get("figure_id", ""))
        if row.get("status") != "processed_calculated_prediction_not_experimental":
            auditor.invalid("NMR_PREDICTIONS.summary.status", section, figure_id, "Prediction result status is not a successful calculated prediction", observed=row.get("status"))
            valid = False
        if not has_calculated_notice(row.get("notice")):
            auditor.invalid("NMR_PREDICTIONS.summary.notice", section, figure_id, "Prediction result lacks calculated/not-experimental notice")
            valid = False
        for field in ("shielding_population_coverage", "coupling_population_coverage"):
            if not close_to(row.get(field), 1.0):
                auditor.invalid("NMR_PREDICTIONS.summary.population_coverage", section, figure_id, f"{field} is below 100%", expected=1.0, observed=row.get(field))
                valid = False

    bundles: dict[str, dict[str, Any]] = {}
    for figure_id in EXPECTED_IDS:
        bundle_path = root / figure_id / f"{figure_id}_prediction_bundle.json"
        bundle = load_required(auditor, bundle_path, section, figure_id, "JSON")
        if bundle is None:
            valid = False
            continue
        if not isinstance(bundle, dict):
            auditor.invalid("NMR_PREDICTIONS.bundle.type", section, figure_id, "Prediction bundle is not an object", paths=[bundle_path])
            valid = False
            continue
        bundles[figure_id] = bundle
        errors: list[str] = []
        if str(bundle.get("compound", {}).get("figure_id", "")) != figure_id:
            errors.append("bundle compound ID mismatch")
        if not has_calculated_notice(bundle.get("notice")):
            errors.append("calculated/not-experimental notice absent")
        weighting = bundle.get("weighting") if isinstance(bundle.get("weighting"), dict) else {}
        for field in (
            "shielding_selected_population_of_all_valid_dft",
            "shielding_completed_population_of_all_valid_dft",
            "shielding_completed_fraction_of_selected_population",
            "coupling_selected_population_of_all_valid_dft",
            "coupling_completed_population_of_all_valid_dft",
            "coupling_completed_fraction_of_selected_population",
        ):
            if not close_to(weighting.get(field), 1.0):
                errors.append(f"{field} is not 1.0")
        if not close_to(weighting.get("temperature_K"), 298.15, 0.011):
            errors.append("Boltzmann temperature is not 298.15 K")
        conformers = bundle.get("conformers") if isinstance(bundle.get("conformers"), list) else []
        for conformer in conformers:
            if not isinstance(conformer, dict):
                errors.append("non-object conformer audit row")
                continue
            if conformer.get("weighting_status") == "included":
                if as_bool(conformer.get("qrrho_verified")) is not True:
                    errors.append("included conformer lacks QRRHO verification")
                if str(conformer.get("dft_audit_validation_status", "")) != "validated_unique_minimum_and_weights":
                    errors.append("included conformer failed DFT audit validation")
                if conformer.get("processing_error") not in (None, ""):
                    errors.append("included conformer has a processing error")
        source_files = bundle.get("source_files")
        if not isinstance(source_files, list) or not source_files:
            errors.append("source_files hash ledger is absent")
        else:
            source_paths: set[str] = set()
            for record in source_files:
                if not isinstance(record, dict):
                    errors.append("non-object source-file record")
                else:
                    errors.extend(verify_source_record(auditor, record))
                    try:
                        source_paths.add(os.path.normcase(str(Path(str(record.get("path", ""))).resolve())))
                    except OSError:
                        pass
            mandatory_sources = (
                auditor.project / "data" / "DFT_job_manifest.csv",
                auditor.project / "data" / "DFT_minimum_and_weight_audit.csv",
                auditor.project / "data" / "NMR_shielding_job_manifest.csv",
                auditor.project / "data" / "NMR_scalar_coupling_job_manifest.csv",
            )
            for mandatory in mandatory_sources:
                if os.path.normcase(str(mandatory.resolve())) not in source_paths:
                    errors.append(f"canonical source manifest absent from hash ledger: {mandatory.name}")
            for record in source_files:
                if not isinstance(record, dict):
                    continue
                source_path = Path(str(record.get("path", "")))
                if "03_nmr_properties" in str(source_path).casefold() and not path_is_within(
                    source_path, auditor.project / "03_nmr_properties"
                ):
                    errors.append(f"noncanonical property source in hash ledger: {source_path}")
        atom_map = bundle.get("atom_map") if isinstance(bundle.get("atom_map"), list) else []
        averaged = bundle.get("averaged_nuclei") if isinstance(bundle.get("averaged_nuclei"), list) else []
        expected_elements = Counter(str(row.get("element", "")) for row in atom_map if str(row.get("element", "")) in {"H", "C", "F", "P"})
        observed_elements = Counter(str(row.get("element", "")) for row in averaged if isinstance(row, dict))
        if expected_elements != observed_elements:
            errors.append(f"averaged-nucleus coverage mismatch: expected {dict(expected_elements)}, got {dict(observed_elements)}")
        references = bundle.get("reference_shieldings") if isinstance(bundle.get("reference_shieldings"), dict) else {}
        for element in expected_elements:
            reference = references.get(element) if isinstance(references.get(element), dict) else {}
            if reference.get("status") != "available_calculated_same_level_reference":
                errors.append(f"{element} same-level raw reference is unavailable")
            reference_path = Path(str(reference.get("reference_output", "")))
            if not reference_path.is_file() or str(reference.get("reference_output_sha256", "")).casefold() != auditor.sha256(reference_path):
                errors.append(f"{element} reference output hash mismatch")
        for nucleus in averaged:
            if not isinstance(nucleus, dict):
                continue
            if nucleus.get("reference_status") != "available_calculated_same_level_reference":
                errors.append(f"{nucleus.get('atom_label_1based')} lacks same-level raw reference")
            if as_float(nucleus.get("raw_reference_shift_ppm")) is None:
                errors.append(f"{nucleus.get('atom_label_1based')} lacks primary raw-reference shift")
        regressions = bundle.get("provisional_external_regressions")
        if not isinstance(regressions, dict):
            errors.append("secondary external-regression metadata is absent")
        else:
            for element in ("H", "C"):
                regression = regressions.get(element)
                if not isinstance(regression, dict) or regression.get("status") != "provisional_external_regression":
                    errors.append(f"{element} external regression is not explicitly labelled provisional/secondary")
        assumptions = " ".join(str(value) for value in bundle.get("assumptions_and_limitations", [])).casefold()
        if "raw chemical shifts use" not in assumptions or "external h/c regression" not in assumptions or "provisional" not in assumptions:
            errors.append("bundle does not clearly distinguish primary raw shifts from provisional regression")
        expected_files = [
            root / figure_id / f"{figure_id}_conformer_QRRHO_weights.csv",
            root / figure_id / f"{figure_id}_calculated_shifts.csv",
            root / figure_id / f"{figure_id}_calculated_scalar_couplings.csv",
            root / figure_id / f"{figure_id}_README.txt",
        ]
        if any(not path.is_file() or path.stat().st_size == 0 for path in expected_files):
            errors.append("one or more mandatory per-compound prediction files are absent/empty")
        carbon_count = observed_elements.get("C", 0)
        formula_carbon = formula_elements(str(bundle.get("compound", {}).get("formula", ""))).get("C", 0)
        if figure_id in NO_CARBON_IDS and (carbon_count != 0 or formula_carbon != 0):
            errors.append("carbon-free compound unexpectedly contains/calculates carbon")
        if figure_id not in NO_CARBON_IDS and (carbon_count == 0 or formula_carbon == 0):
            errors.append("carbon-containing compound has no calculated/formula carbon")
        if errors:
            auditor.invalid("NMR_PREDICTIONS.bundle.integrity", section, figure_id, "; ".join(errors[:30]), paths=[bundle_path])
            valid = False
        else:
            auditor.pass_("NMR_PREDICTIONS.bundle.integrity", section, figure_id, "Bundle, source hashes, 100% population coverage, and primary same-level reference verified", paths=[bundle_path])
    if valid and len(bundles) == 51:
        auditor.pass_("NMR_PREDICTIONS.all_valid", section, "all compounds", "All 51 calculated NMR prediction bundles passed strict audit", observed=51)
    return valid and len(bundles) == 51, bundles


def audit_nmr_spectra(auditor: Auditor, bundles: Mapping[str, Mapping[str, Any]]) -> bool:
    section = "NMR_SPECTRA"
    root = auditor.project / "NMR_spectra_CALCULATED_NOT_EXPERIMENTAL"
    summary_path = root / "spectrum_generation_summary__ALL.json"
    summary = load_required(auditor, summary_path, section, "ALL summary", "JSON")
    directory_ids = {path.name for path in root.glob("F3_*") if path.is_dir()} if root.is_dir() else set()
    extras = directory_ids - EXPECTED_ID_SET
    if extras:
        auditor.invalid("NMR_SPECTRA.directories.extra", section, "output directories", "Unexpected/F3_005 spectrum directories are present", observed=sorted(extras), paths=[root])
    if summary is None:
        auditor.info("NMR_SPECTRA.directory_progress", section, "output directories", "Canonical spectrum directory progress", expected=51, observed=len(directory_ids & EXPECTED_ID_SET))
        return False
    if not isinstance(summary, dict) or not isinstance(summary.get("results"), list):
        auditor.invalid("NMR_SPECTRA.summary.type", section, "ALL summary", "Spectrum summary structure is invalid", paths=[summary_path])
        return False
    selected = [str(value) for value in summary.get("selected_figure_ids", [])]
    results = [row for row in summary["results"] if isinstance(row, dict)]
    result_ids = [str(row.get("figure_id", "")) for row in results]
    valid = audit_id_set(auditor, section, "summary", selected)
    if set(result_ids) != EXPECTED_ID_SET or len(result_ids) != 51 or len(set(result_ids)) != 51:
        auditor.invalid("NMR_SPECTRA.summary.result_ids", section, "ALL summary", "Spectrum summary results are not the exact unique 51-ID set", observed=result_ids)
        valid = False
    if set(directory_ids) != EXPECTED_ID_SET:
        missing = sorted(EXPECTED_ID_SET - directory_ids)
        if missing:
            auditor.incomplete("NMR_SPECTRA.directories.missing", section, "output directories", "Spectrum directories are incomplete", observed={"missing": missing})
        valid = False
    result_map = {str(row.get("figure_id", "")): row for row in results}
    for figure_id in EXPECTED_IDS:
        result = result_map.get(figure_id, {})
        expected_carbon_status = "not_applicable" if figure_id in NO_CARBON_IDS else "calculated_profile_generated"
        if result.get("status") != "calculated_spectra_generated_not_experimental" or result.get("carbon_status") != expected_carbon_status:
            auditor.invalid(
                "NMR_SPECTRA.summary.status",
                section,
                figure_id,
                "Spectrum summary status/carbon disposition is incorrect",
                expected={"status": "calculated_spectra_generated_not_experimental", "carbon_status": expected_carbon_status},
                observed={"status": result.get("status"), "carbon_status": result.get("carbon_status")},
            )
            valid = False
        metadata_path = root / figure_id / f"{figure_id}_spectrum_generation_metadata.json"
        metadata = load_required(auditor, metadata_path, section, figure_id, "JSON")
        if metadata is None:
            valid = False
            continue
        if not isinstance(metadata, dict):
            auditor.invalid("NMR_SPECTRA.metadata.type", section, figure_id, "Spectrum metadata is not an object")
            valid = False
            continue
        errors: list[str] = []
        if not has_calculated_notice(metadata.get("notice")):
            errors.append("calculated/not-experimental notice absent")
        if str(metadata.get("compound", {}).get("figure_id", "")) != figure_id:
            errors.append("metadata compound ID mismatch")
        source = metadata.get("input_prediction_bundle") if isinstance(metadata.get("input_prediction_bundle"), dict) else {}
        expected_bundle_path = auditor.project / "NMR_predictions_CALCULATED_NOT_EXPERIMENTAL" / figure_id / f"{figure_id}_prediction_bundle.json"
        if not paths_equal(str(source.get("path", "")), expected_bundle_path):
            errors.append("input prediction-bundle path mismatch")
        elif not expected_bundle_path.is_file() or str(source.get("sha256", "")).casefold() != auditor.sha256(expected_bundle_path):
            errors.append("input prediction-bundle hash mismatch")
        elif int(source.get("size_bytes") or -1) != expected_bundle_path.stat().st_size:
            errors.append("input prediction-bundle size mismatch")
        policy = metadata.get("data_policy") if isinstance(metadata.get("data_policy"), dict) else {}
        for field in ("new_random_or_synthetic_measurements", "noise_added", "baseline_added", "solvent_peak_added", "water_peak_added"):
            if as_bool(policy.get(field)) is not False:
                errors.append(f"data_policy.{field} is not false")
        variants = metadata.get("shift_variants") if isinstance(metadata.get("shift_variants"), dict) else {}
        keys = list(variants)
        raw = variants.get("raw_TMS") if isinstance(variants.get("raw_TMS"), dict) else {}
        provisional = variants.get("provisional_regression") if isinstance(variants.get("provisional_regression"), dict) else {}
        if not keys or keys[0] != "raw_TMS" or raw.get("shift_field") != "raw_tms_reference_shift_ppm":
            errors.append("raw same-level TMS view is not the primary/first shift variant")
        if provisional.get("shift_field") != "provisional_external_regression_shift_ppm" or "provisional" not in str(provisional.get("label", "")).casefold():
            errors.append("external regression is not explicitly secondary/provisional")
        weights = metadata.get("upstream_weighting") if isinstance(metadata.get("upstream_weighting"), dict) else {}
        for field in (
            "shielding_selected_population_of_all_valid_dft",
            "shielding_completed_population_of_all_valid_dft",
            "shielding_completed_fraction_of_selected_population",
            "coupling_selected_population_of_all_valid_dft",
            "coupling_completed_population_of_all_valid_dft",
            "coupling_completed_fraction_of_selected_population",
        ):
            if not close_to(weights.get(field), 1.0):
                errors.append(f"upstream_weighting.{field} is not 1.0")
        carbon = metadata.get("carbon_status") if isinstance(metadata.get("carbon_status"), dict) else {}
        if carbon.get("status") != expected_carbon_status:
            errors.append("metadata carbon_status mismatch")
        if figure_id in NO_CARBON_IDS and "carbon-free" not in str(carbon.get("reason", "")).casefold():
            errors.append("13C N/A reason does not explicitly identify the carbon-free reference state")
        outputs = metadata.get("output_files")
        if not isinstance(outputs, list) or not outputs:
            errors.append("output-file hash ledger is absent")
            outputs = []
        output_paths: set[str] = set()
        for record in outputs:
            if not isinstance(record, dict):
                errors.append("non-object output-file record")
                continue
            path = Path(str(record.get("path", "")))
            output_paths.add(path.name)
            if not path_is_within(path, root / figure_id):
                errors.append(f"output path outside compound directory: {path}")
            elif not path.is_file() or path.stat().st_size == 0:
                errors.append(f"output missing/empty: {path.name}")
            elif not valid_hash(record.get("sha256")) or str(record.get("sha256")).casefold() != auditor.sha256(path):
                errors.append(f"output hash mismatch: {path.name}")
            elif int(record.get("size_bytes") or -1) != path.stat().st_size:
                errors.append(f"output size mismatch: {path.name}")
        if len(output_paths) != len(outputs):
            errors.append("duplicate output-file records are present")
        mandatory_names = {
            f"{figure_id}_1H_raw_TMS_sticks.csv",
            f"{figure_id}_1H_raw_TMS_profile.txt",
            f"{figure_id}_1H_raw_TMS_calculated.png",
            f"{figure_id}_1H_raw_TMS_calculated.svg",
            f"{figure_id}_1H_raw_TMS_calculated.pdf",
            f"{figure_id}_1H_provisional_regression_sticks.csv",
            f"{figure_id}_1H_provisional_regression_profile.txt",
            f"{figure_id}_1H_provisional_regression_calculated.png",
            f"{figure_id}_1H_provisional_regression_calculated.svg",
            f"{figure_id}_1H_provisional_regression_calculated.pdf",
            f"{figure_id}_NMR_prediction_comparison.png",
            f"{figure_id}_NMR_prediction_comparison.svg",
            f"{figure_id}_NMR_prediction_comparison.pdf",
            f"{figure_id}_spin_component_audit.csv",
            f"{figure_id}_spectrum_generation_README.txt",
        }
        carbon_names = {
            f"{figure_id}_13C_{variant}_{suffix}"
            for variant in ("raw_TMS", "provisional_regression")
            for suffix in ("sticks.csv", "profile.txt", "calculated.png", "calculated.svg", "calculated.pdf")
        }
        if figure_id in NO_CARBON_IDS:
            mandatory_names.add(f"{figure_id}_13C_NOT_APPLICABLE.txt")
            if any(
                f"{figure_id}_13C_raw_TMS_" in name
                or f"{figure_id}_13C_provisional_regression_" in name
                for name in output_paths
            ):
                errors.append("calculated 13C profile exists for a carbon-free compound")
        else:
            mandatory_names.update(carbon_names)
        missing_names = mandatory_names - output_paths
        if missing_names:
            errors.append("mandatory output(s) absent from hash ledger: " + ", ".join(sorted(missing_names)))
        if errors:
            auditor.invalid("NMR_SPECTRA.metadata.integrity", section, figure_id, "; ".join(errors[:35]), paths=[metadata_path])
            valid = False
        else:
            auditor.pass_("NMR_SPECTRA.metadata.integrity", section, figure_id, "Plot source hash, output hashes, primary raw view, and carbon disposition verified", paths=[metadata_path])
    if valid:
        auditor.pass_("NMR_SPECTRA.all_valid", section, "all compounds", "All 51 calculated NMR spectrum sets passed strict audit", observed=51)
    return valid


def audit_hrms(auditor: Auditor) -> bool:
    section = "HRMS"
    root = auditor.project / "HRMS_formula_predictions_CALCULATED_NOT_EXPERIMENTAL"
    run_path = root / "RUN_METADATA.json"
    summary_path = root / "Figure3_formula_HRMS_prediction_summary.csv"
    run = load_required(auditor, run_path, section, "run metadata", "JSON")
    rows = load_required(auditor, summary_path, section, "summary", "CSV", comments=True)
    if run is None or rows is None or not isinstance(run, dict) or not isinstance(rows, list):
        return False
    valid = True
    processed = [str(value) for value in run.get("processed_figure_ids", [])]
    if not audit_id_set(auditor, section, "run metadata", processed):
        valid = False
    skipped = run.get("skipped") if isinstance(run.get("skipped"), list) else []
    skipped_ids = {str(item.get("figure_id", "")) for item in skipped if isinstance(item, dict)}
    if skipped_ids != {"F3_005"}:
        auditor.invalid("HRMS.run.skipped", section, "run metadata", "The skipped set must be exactly F3_005", expected=["F3_005"], observed=sorted(skipped_ids))
        valid = False
    scope = str(run.get("calculation_scope", "")).casefold()
    if not has_calculated_notice(run.get("prediction_status")) or "formula-level" not in scope or "precursor/adduct" not in scope or "no fragment" not in scope or run.get("fragment_ions_generated") is not False:
        auditor.invalid("HRMS.run.scope", section, "run metadata", "HRMS run is not strictly labelled formula-only precursor/adduct with fragment=false", paths=[run_path])
        valid = False
    summary_ids = [str(row.get("figure_id", "")) for row in rows]
    if not audit_id_set(auditor, section, "summary", summary_ids):
        valid = False
    if any(
        not has_calculated_notice(row.get("prediction_status"))
        or as_bool(row.get("fragment_ions_generated")) is not False
        for row in rows
    ):
        auditor.invalid("HRMS.summary.scope", section, "summary", "At least one HRMS summary row is not calculated/formula-only fragment=false", paths=[summary_path])
        valid = False
    directory_ids = {path.name for path in root.glob("F3_*") if path.is_dir()}
    if directory_ids != EXPECTED_ID_SET:
        extra = sorted(directory_ids - EXPECTED_ID_SET)
        missing = sorted(EXPECTED_ID_SET - directory_ids)
        if extra:
            auditor.invalid("HRMS.directories.extra", section, "output directories", "Unexpected/F3_005 HRMS directories are present", observed=extra)
        if missing:
            auditor.incomplete("HRMS.directories.missing", section, "output directories", "HRMS output directories are incomplete", observed=missing)
        valid = False
    summary_by_key = Counter((str(row.get("figure_id", "")), str(row.get("ion_label", ""))) for row in rows)
    ion_count = 0
    for figure_id in EXPECTED_IDS:
        metadata_path = root / figure_id / f"{figure_id}_calculation_metadata.json"
        metadata = load_required(auditor, metadata_path, section, figure_id, "JSON")
        if metadata is None or not isinstance(metadata, dict):
            valid = False
            continue
        errors: list[str] = []
        if str(metadata.get("figure_id", "")) != figure_id:
            errors.append("metadata Figure ID mismatch")
        scope = str(metadata.get("calculation_scope", "")).casefold()
        if not has_calculated_notice(metadata.get("prediction_status")) or "formula-level" not in scope or "precursor/adduct" not in scope or "no fragment" not in scope:
            errors.append("metadata is not explicitly calculated formula-only precursor/adduct")
        if metadata.get("fragment_ions_generated") is not False:
            errors.append("top-level fragment_ions_generated is not false")
        models = metadata.get("ion_models") if isinstance(metadata.get("ion_models"), list) else []
        if not models:
            errors.append("no formula-level ion model")
        ion_count += len(models)
        for model in models:
            if not isinstance(model, dict):
                errors.append("non-object ion model")
                continue
            label = str(model.get("label", ""))
            if model.get("fragment_ions_generated") is not False:
                errors.append(f"{label} fragment_ions_generated is not false")
            if summary_by_key[(figure_id, label)] != 1:
                errors.append(f"{label} lacks exactly one matching summary row")
        sticks = list((root / figure_id).glob("*_isotope_sticks.csv"))
        profiles = list((root / figure_id).glob("*_gaussian_profile_R60000.csv"))
        plots = list((root / figure_id).glob("*_gaussian_profile_R60000.png"))
        if not (len(sticks) == len(profiles) == len(plots) == len(models)):
            errors.append(f"ion-model/file count mismatch models={len(models)}, sticks={len(sticks)}, profiles={len(profiles)}, png={len(plots)}")
        if any(path.stat().st_size == 0 for path in sticks + profiles + plots):
            errors.append("empty HRMS output file")
        if errors:
            auditor.invalid("HRMS.compound.integrity", section, figure_id, "; ".join(errors), paths=[metadata_path])
            valid = False
        else:
            auditor.pass_("HRMS.compound.integrity", section, figure_id, "Formula-only precursor/adduct metadata, fragment=false, summary, and output files verified", paths=[metadata_path])
    if ion_count != len(rows):
        auditor.invalid("HRMS.summary.ion_count", section, "all compounds", "Total ion-model count differs from summary rows", expected=ion_count, observed=len(rows))
        valid = False
    if valid:
        auditor.pass_("HRMS.all_valid", section, "all compounds", "All 51 formula-only HRMS precursor/adduct envelopes passed strict audit", observed={"compounds": 51, "ion_models": ion_count})
    return valid


def audit_pipeline_status(auditor: Auditor) -> None:
    path = auditor.project / "logs" / "full_pipeline" / "pipeline_status.json"
    if not path.is_file():
        auditor.info("PIPELINE.status.absent", "PIPELINE", "driver status", "No full-pipeline driver status file is present", paths=[path])
        return
    try:
        payload = read_json(path)
    except (OSError, json.JSONDecodeError, UnicodeError) as exc:
        auditor.invalid("PIPELINE.status.parse", "PIPELINE", "driver status", f"Malformed pipeline status: {exc}", paths=[path])
        return
    if not isinstance(payload, dict):
        auditor.invalid("PIPELINE.status.type", "PIPELINE", "driver status", "Pipeline status is not an object", paths=[path])
        return
    status = str(payload.get("status", ""))
    current_stage = str(payload.get("current_stage", ""))
    auditor.info(
        "PIPELINE.status.snapshot",
        "PIPELINE",
        "driver status",
        "Pipeline-driver state is informational; artifact checks determine the audit verdict",
        observed={"status": status, "current_stage": current_stage, "updated_at": payload.get("updated_at")},
        paths=[path],
    )
    if status == "complete" and current_stage == "complete":
        auditor.pass_("PIPELINE.status.complete", "PIPELINE", "driver status", "Pipeline driver reports completion", paths=[path])


def run_audit(project: Path) -> Auditor:
    auditor = Auditor(project)
    if not auditor.project.is_dir():
        auditor.invalid("PROJECT.root", "PROJECT", "root", "Project root does not exist", paths=[auditor.project])
        return auditor
    auditor.pass_("PROJECT.root", "PROJECT", "root", "Project root is readable", paths=[auditor.project])
    audit_pipeline_status(auditor)
    goat_complete, goat_rows = audit_goat(auditor)
    if not goat_complete:
        auditor.info("PIPELINE.goat_gate", "PIPELINE", "GOAT gate", "GOAT production is not yet fully current-hash verified; downstream canonical outputs may legitimately be absent")

    dft_complete, dft_rows, _dft_status = audit_orca_stage(
        auditor,
        section="DFT",
        manifest_path=auditor.project / "data" / "DFT_job_manifest.csv",
        status_root=auditor.project / "logs" / "dft_production",
        expected_stage="dft_production",
        job_root=auditor.project / "02_dft_opt_freq",
        suffix=None,
    )
    dft_manifest_has_full_ids = {str(row.get("figure_id", "")) for row in dft_rows} == EXPECTED_ID_SET
    if dft_rows and dft_manifest_has_full_ids:
        audit_dft_goat_provenance(auditor, dft_rows, goat_rows)
    elif dft_rows:
        auditor.info(
            "DFT.source_goat.deferred",
            "DFT",
            "source GOAT provenance",
            "Deep source-provenance validation is deferred until the canonical DFT manifest covers all 51 IDs",
        )

    shield_complete, shielding_rows, _shield_status = audit_orca_stage(
        auditor,
        section="NMR_SHIELDING_JOBS",
        manifest_path=auditor.project / "data" / "NMR_shielding_job_manifest.csv",
        status_root=auditor.project / "logs" / "shielding_production",
        expected_stage="shielding_production",
        job_root=auditor.project / "03_nmr_properties",
        suffix="_shield",
    )
    scalar_complete, scalar_rows, _scalar_status = audit_orca_stage(
        auditor,
        section="NMR_SCALAR_JOBS",
        manifest_path=auditor.project / "data" / "NMR_scalar_coupling_job_manifest.csv",
        status_root=auditor.project / "logs" / "scalar_coupling_production",
        expected_stage="scalar_coupling_production",
        job_root=auditor.project / "03_nmr_properties",
        suffix="_scalar_coupling",
    )
    if shielding_rows and scalar_rows:
        audit_property_strictness(auditor, dft_rows, shielding_rows, scalar_rows)
    if dft_rows and shielding_rows and scalar_rows:
        audit_boltzmann(auditor, dft_rows, shielding_rows, scalar_rows)

    predictions_complete, bundles = audit_nmr_predictions(auditor)
    spectra_complete = audit_nmr_spectra(auditor, bundles)
    hrms_complete = audit_hrms(auditor)
    auditor.info(
        "PIPELINE.stage_gates",
        "PIPELINE",
        "strict stage gates",
        "Strict stage completion snapshot",
        observed={
            "GOAT": goat_complete,
            "DFT": dft_complete,
            "shielding": shield_complete,
            "scalar_coupling": scalar_complete,
            "NMR_predictions": predictions_complete,
            "NMR_spectra": spectra_complete,
            "HRMS_formula_only": hrms_complete,
        },
    )
    return auditor


def json_payload(auditor: Auditor) -> dict[str, Any]:
    summary = auditor.summary()
    return {
        "schema_version": "1.0.0",
        "audit_kind": "read_only_final_output_integrity_and_completeness",
        "generated_utc": now_utc(),
        "project": str(auditor.project),
        "expected_figure_ids": list(EXPECTED_IDS),
        "excluded_figure_ids": ["F3_005"],
        "expected_goat_jobs": EXPECTED_GOAT_JOBS,
        "strict_policy": {
            "calculated_not_experimental_required": True,
            "goat_current_input_hash_required": True,
            "exact_manifest_status_coverage_required": True,
            "all_valid_unique_minima_property_population_required": 1.0,
            "allow_partial_forbidden": True,
            "imaginary_frequency_override_forbidden": True,
            "incomplete_dft_override_forbidden": True,
            "provisional_build_tag_forbidden": True,
            "external_HC_regression_allowed_only_as_secondary": True,
            "raw_same_level_reference_is_primary": True,
            "hrms_formula_only": True,
            "hrms_fragment_ions_generated": False,
        },
        "summary": summary,
        "checks": [asdict(check) for check in auditor.checks],
    }


def markdown_report(payload: Mapping[str, Any]) -> str:
    summary = payload["summary"]
    lines = [
        "# Final output audit",
        "",
        f"- Generated (UTC): `{payload['generated_utc']}`",
        f"- Project: `{payload['project']}`",
        f"- Overall status: **{summary['overall_status']}**",
        f"- Exit code: `{summary['exit_code']}`",
        f"- Checks: PASS {summary['counts']['PASS']}; INFO {summary['counts']['INFO']}; INCOMPLETE {summary['counts']['INCOMPLETE']}; INVALID {summary['counts']['INVALID']}",
        "",
        "This is a read-only audit of calculated predictions, not experimental characterization.",
        "The external H/C regression is accepted only as a clearly labelled secondary view; same-level raw referencing is primary.",
        "",
        "## Section summary",
        "",
        "| Section | PASS | INFO | INCOMPLETE | INVALID |",
        "|---|---:|---:|---:|---:|",
    ]
    for section, counts in summary["sections"].items():
        lines.append(f"| {section} | {counts['PASS']} | {counts['INFO']} | {counts['INCOMPLETE']} | {counts['INVALID']} |")
    issues = [check for check in payload["checks"] if check["status"] in {"INCOMPLETE", "INVALID"}]
    lines.extend(["", "## Blocking findings", ""])
    if not issues:
        lines.append("No blocking findings.")
    else:
        lines.extend(["| Status | Section | Subject | Finding |", "|---|---|---|---|"])
        for check in issues:
            message = str(check["message"]).replace("|", "\\|").replace("\n", " ")
            subject = str(check["subject"]).replace("|", "\\|")
            lines.append(f"| {check['status']} | {check['section']} | {subject} | {message} |")
    lines.extend(
        [
            "",
            "## Exit-code contract",
            "",
            "- `0`: COMPLETE",
            "- `2`: INCOMPLETE (jobs/results still absent or unfinished)",
            "- `3`: INVALID (integrity, provenance, strict-method, or semantic failure)",
            "- `4`: AUDITOR_ERROR",
            "",
        ]
    )
    return "\n".join(lines)


def csv_text(payload: Mapping[str, Any]) -> str:
    fields = ("check_id", "section", "subject", "status", "message", "expected", "observed", "paths")
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    for check in payload["checks"]:
        row = dict(check)
        for field in ("expected", "observed", "paths"):
            row[field] = "" if row.get(field) is None else json.dumps(row[field], ensure_ascii=False, sort_keys=True)
        writer.writerow({field: row.get(field, "") for field in fields})
    return stream.getvalue()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_reports(payload: Mapping[str, Any], output_dir: Path) -> list[Path]:
    output_dir = output_dir.resolve()
    paths = [
        output_dir / "final_output_audit.json",
        output_dir / "final_output_audit.csv",
        output_dir / "final_output_audit.md",
    ]
    atomic_write(paths[0], json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    atomic_write(paths[1], csv_text(payload))
    atomic_write(paths[2], markdown_report(payload))
    return paths


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description="Read-only strict audit of all 51 final ORCA/NMR/formula-HRMS outputs")
    result.add_argument("--project-root", type=Path, default=DEFAULT_PROJECT)
    result.add_argument(
        "--output-dir",
        type=Path,
        help="Audit-report directory only (default: <project>/audit_final_outputs)",
    )
    result.add_argument("--stdout-only", action="store_true", help="Do not write audit report files; emit JSON to stdout")
    result.add_argument("--compact", action="store_true", help="Emit compact rather than pretty JSON to stdout")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    project = args.project_root.expanduser().resolve()
    try:
        auditor = run_audit(project)
        payload = json_payload(auditor)
        if args.stdout_only:
            print(json.dumps(payload, ensure_ascii=False, indent=None if args.compact else 2))
        else:
            output_dir = (args.output_dir or (project / "audit_final_outputs")).expanduser().resolve()
            paths = write_reports(payload, output_dir)
            print(json.dumps({"summary": payload["summary"], "reports": [str(path) for path in paths]}, ensure_ascii=False, indent=2))
        return int(payload["summary"]["exit_code"])
    except Exception as exc:  # last-resort machine-readable failure
        failure = {
            "schema_version": "1.0.0",
            "audit_kind": "read_only_final_output_integrity_and_completeness",
            "generated_utc": now_utc(),
            "project": str(project),
            "summary": {"overall_status": "AUDITOR_ERROR", "exit_code": 4},
            "error": {"type": type(exc).__name__, "message": str(exc)},
        }
        print(json.dumps(failure, ensure_ascii=False, indent=2), file=sys.stderr)
        return 4


if __name__ == "__main__":
    raise SystemExit(main())
