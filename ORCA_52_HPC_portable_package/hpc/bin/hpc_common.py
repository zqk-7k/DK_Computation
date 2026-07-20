from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping


NORMAL_MARKER = "ORCA TERMINATED NORMALLY"
TRUE_VALUES = {"1", "true", "yes", "y", "on", "enabled", "ready", "approved", "pass"}
FALSE_VALUES = {"0", "false", "no", "n", "off", "disabled"}
WINDOWS_ABSOLUTE = re.compile(r"^[A-Za-z]:[\\/]")
XYZFILE_PATTERN = re.compile(
    r"(?im)^\s*\*\s*(?:xyzfile|gzmtfile)\s+[-+]?\d+\s+\d+\s+"
    r'(?:"([^"]+)"|\'([^\']+)\'|(\S+))'
)


def utc_now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        dir=str(path.parent), prefix=f".{path.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def atomic_json(path: Path, payload: object) -> None:
    atomic_text(path, json.dumps(payload, ensure_ascii=False, indent=2) + "\n")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest().upper()


def canonical_row_hash(row: Mapping[str, Any]) -> str:
    value = json.dumps(dict(row), sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return sha256_text(value)


def read_manifest(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"Manifest has no header: {path}")
        rows = [{str(k): str(v or "").strip() for k, v in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def parse_bool(value: Any, default: bool = True) -> bool:
    token = str(value or "").strip().casefold()
    if not token:
        return default
    if token in TRUE_VALUES:
        return True
    if token in FALSE_VALUES:
        return False
    raise ValueError(f"Invalid boolean token: {value!r}")


def row_value(row: Mapping[str, str], *names: str) -> str:
    lowered = {str(key).casefold(): str(value).strip() for key, value in row.items()}
    for name in names:
        value = lowered.get(name.casefold(), "")
        if value:
            return value
    return ""


def safe_component(value: str, fallback: str = "job") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()).strip("._")
    return cleaned or fallback


def is_windows_absolute(value: str) -> bool:
    return bool(WINDOWS_ABSOLUTE.match(value.strip()))


def resolve_package_path(
    package_root: Path,
    raw: str,
    *,
    field: str,
    allow_missing: bool = False,
) -> Path:
    if not raw:
        raise ValueError(f"Missing required path field: {field}")
    if is_windows_absolute(raw):
        raise ValueError(
            f"{field} contains a Windows absolute path and is not portable: {raw}. "
            "Regenerate the HPC manifest with package-relative paths."
        )
    portable = raw.replace("\\", "/")
    candidate = Path(portable)
    if not candidate.is_absolute():
        candidate = package_root / candidate
    resolved = candidate.resolve()
    root = package_root.resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"{field} escapes the package root: {raw}") from exc
    if not allow_missing and not resolved.exists():
        raise FileNotFoundError(f"{field} does not exist: {resolved}")
    return resolved


def resolve_manifest_path(package_root: Path, raw: str | Path) -> Path:
    text = str(raw)
    candidate = Path(text)
    if not candidate.is_absolute():
        candidate = package_root / candidate
    result = candidate.resolve()
    try:
        result.relative_to(package_root.resolve())
    except ValueError as exc:
        raise ValueError(f"Manifest must be inside the package: {result}") from exc
    if not result.is_file():
        raise FileNotFoundError(result)
    return result


def input_coordinate_references(input_path: Path) -> list[Path]:
    text = input_path.read_text(encoding="utf-8", errors="replace")
    paths: list[Path] = []
    for match in XYZFILE_PATTERN.finditer(text):
        token = next(part for part in match.groups() if part is not None)
        candidate = Path(token)
        if not candidate.is_absolute():
            candidate = input_path.parent / candidate
        paths.append(candidate.resolve())
    return paths


def parse_orca_resources(input_path: Path) -> dict[str, int]:
    text = input_path.read_text(encoding="utf-8", errors="replace")
    nprocs_matches = re.findall(r"(?i)\bnprocs\s+(\d+)", text)
    worker_matches = re.findall(r"(?i)\bNWorkers\s+(\d+)", text)
    maxcore_matches = re.findall(r"(?i)^\s*%maxcore\s+(\d+)", text, flags=re.MULTILINE)
    nprocs = max((int(value) for value in nprocs_matches), default=1)
    nworkers = max((int(value) for value in worker_matches), default=1)
    maxcore = max((int(value) for value in maxcore_matches), default=0)
    return {
        "nprocs": nprocs,
        "nworkers": nworkers,
        "required_cpu_slots": max(nprocs, nworkers),
        "maxcore_mb_per_rank": maxcore,
    }


def enabled_row(row: Mapping[str, str]) -> bool:
    return parse_bool(row_value(row, "enabled", "run", "active"), default=True)


def approval_problem(row: Mapping[str, str]) -> str:
    readiness = row_value(row, "readiness")
    if readiness and readiness.casefold() not in TRUE_VALUES:
        return f"readiness={readiness!r}"
    approval = row_value(row, "approval_status")
    if approval and any(token in approval.casefold() for token in ("draft", "not_approved", "blocked", "reject")):
        return f"approval_status={approval!r}"
    return ""


def glob_nonempty(directory: Path, pattern: str) -> list[Path]:
    return sorted(path for path in directory.glob(pattern) if path.is_file() and path.stat().st_size > 0)


def write_csv(path: Path, fieldnames: Iterable[str], rows: Iterable[Mapping[str, Any]]) -> None:
    field_list = list(fieldnames)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        dir=str(path.parent), prefix=f".{path.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=field_list, extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise
