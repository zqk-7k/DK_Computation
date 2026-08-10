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
STAGE_APPROVALS: dict[str, frozenset[str]] = {
    "goat": frozenset({"APPROVED_FOR_GOAT_CALCULATION"}),
    "optfreq": frozenset({"APPROVED_FOR_OPT_FREQ"}),
    "reference_optfreq": frozenset({"APPROVED_FOR_REFERENCE_GEOMETRY_CALCULATION"}),
    "reference_shielding": frozenset({"APPROVED_FOR_REFERENCE_SHIELDING_CALCULATION"}),
    "shielding": frozenset({"APPROVED_FOR_NMR_PROPERTIES"}),
    "jcoupling": frozenset({"APPROVED_FOR_NMR_PROPERTIES"}),
    "j_method_benchmark": frozenset({"APPROVED_FOR_J_BENCHMARK_CALCULATION"}),
    "high_level_benchmark": frozenset({"APPROVED_FOR_HIGH_LEVEL_BENCHMARK"}),
    "vpt2_force_field": frozenset({"APPROVED_FOR_VPT2"}),
    "vpt2_nmr": frozenset({"APPROVED_FOR_VPT2"}),
}
DEFAULT_PARENT_STAGE: dict[str, str] = {
    "vpt2_nmr": "vpt2_force_field",
}
DEFAULT_DEPENDENCY_ARTIFACTS: dict[tuple[str, str], tuple[str, ...]] = {
    ("vpt2_nmr", "vpt2_force_field"): ("vpt2_ff.hess", "vpt2_ff.vpt2"),
}
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
    # Missing execution intent is never permission to spend HPC resources.
    return parse_bool(row_value(row, "enabled", "run", "active"), default=False)


def approval_problem(row: Mapping[str, str], stage: str | None = None) -> str:
    readiness = row_value(row, "readiness")
    if readiness.casefold() != "ready":
        return f"readiness={readiness!r}"
    approval = row_value(row, "approval_status")
    declared_stage = row_value(row, "stage")
    selected_stage = str(stage or declared_stage).strip()
    if not selected_stage:
        return "stage is missing; no approval allowlist can be selected"
    if declared_stage and declared_stage != selected_stage:
        return f"stage={declared_stage!r} does not match requested stage={selected_stage!r}"
    allowed = STAGE_APPROVALS.get(selected_stage)
    if not allowed:
        return f"stage={selected_stage!r} has no configured approval allowlist"
    if approval not in allowed:
        return f"approval_status={approval!r}; expected one of {sorted(allowed)!r}"
    return ""


def manifest_name_list(
    row: Mapping[str, str], field: str, *, allow_globs: bool = False
) -> tuple[str, ...]:
    raw = row_value(row, field)
    if not raw:
        return ()
    if raw.lstrip().startswith("["):
        parsed = json.loads(raw)
        if not isinstance(parsed, list):
            raise ValueError(f"{field} must contain a JSON list or delimited names")
        values = [str(item).strip() for item in parsed if str(item).strip()]
    else:
        values = [item.strip() for item in re.split(r"[;,|]", raw) if item.strip()]
    unique: list[str] = []
    for value in values:
        token = Path(value)
        if (
            token.is_absolute()
            or token.name != value
            or value in {".", ".."}
            or "/" in value
            or "\\" in value
        ):
            raise ValueError(f"Unsafe {field} entry: {value!r}")
        if not allow_globs and any(character in value for character in "*?[]"):
            raise ValueError(f"Globs are forbidden in {field}: {value!r}")
        if value not in unique:
            unique.append(value)
    return tuple(unique)


def required_artifact_patterns(
    stage_config: Mapping[str, Any], row: Mapping[str, str]
) -> tuple[str, ...]:
    patterns: list[str] = []
    for value in stage_config.get("required_artifact_globs", []):
        token = str(value).strip()
        if token and token not in patterns:
            patterns.append(token)
    for value in manifest_name_list(row, "required_artifacts"):
        if value not in patterns:
            patterns.append(value)
    for value in manifest_name_list(row, "required_output_globs", allow_globs=True):
        if value not in patterns:
            patterns.append(value)
    for value in patterns:
        token = Path(value)
        if (
            token.is_absolute()
            or token.name != value
            or value in {".", ".."}
            or "/" in value
            or "\\" in value
        ):
            raise ValueError(f"Unsafe required artifact pattern: {value!r}")
    return tuple(patterns)


def resolve_required_helpers(
    site_config: Mapping[str, Any], orca_executable: Path, row: Mapping[str, str]
) -> list[dict[str, str]]:
    configured = site_config.get("helper_executables", {})
    if configured is None:
        configured = {}
    if not isinstance(configured, Mapping):
        raise ValueError("site_config helper_executables must be an object")
    records: list[dict[str, str]] = []
    for name in manifest_name_list(row, "required_helpers"):
        raw = str(configured.get(name, "")).strip()
        candidate = Path(raw).expanduser() if raw else orca_executable.parent / name
        if not candidate.is_absolute() or not candidate.is_file() or not os.access(candidate, os.X_OK):
            raise ValueError(
                f"Required ORCA helper {name!r} is unavailable/executable; "
                "configure helper_executables or install it beside ORCA"
            )
        records.append(
            {
                "name": name,
                "path": str(candidate.resolve()),
                "sha256": sha256_file(candidate.resolve()),
            }
        )
    return records


def dependency_contract(
    row: Mapping[str, str], child_stage: str
) -> tuple[str, str, tuple[str, ...]]:
    """Return and validate the declared cross-stage dependency contract.

    Dependency artifact names are deliberately restricted to flat files.  A
    child stage can never authorize an arbitrary path through its manifest.
    """

    parent_id = row_value(row, "depends_on_job_id")
    parent_stage = row_value(row, "depends_on_stage") or DEFAULT_PARENT_STAGE.get(
        child_stage, ""
    )
    if not parent_id:
        if child_stage in DEFAULT_PARENT_STAGE:
            raise ValueError(f"{child_stage} row lacks depends_on_job_id")
        return "", "", ()
    if not parent_stage:
        raise ValueError(
            f"Dependency {parent_id!r} lacks depends_on_stage for child stage {child_stage!r}"
        )
    raw_artifacts = row_value(row, "depends_on_artifacts")
    artifacts: list[str]
    if raw_artifacts:
        if raw_artifacts.lstrip().startswith("["):
            parsed = json.loads(raw_artifacts)
            if not isinstance(parsed, list):
                raise ValueError("depends_on_artifacts must be a JSON list or delimited names")
            artifacts = [str(item).strip() for item in parsed if str(item).strip()]
        else:
            artifacts = [
                item.strip()
                for item in re.split(r"[;,|]", raw_artifacts)
                if item.strip()
            ]
    else:
        artifacts = list(
            DEFAULT_DEPENDENCY_ARTIFACTS.get((child_stage, parent_stage), ())
        )
    if parent_stage != child_stage and not artifacts:
        raise ValueError(
            f"Cross-stage dependency {parent_stage}/{parent_id} declares no required artifacts"
        )
    for name in artifacts:
        token = Path(name)
        if token.is_absolute() or token.name != name or name in {".", ".."}:
            raise ValueError(f"Unsafe dependency artifact name: {name!r}")
    return parent_id, parent_stage, tuple(artifacts)


def verify_materialized_dependency_receipt(
    package_root: Path,
    row: Mapping[str, str],
    job_dir: Path,
    child_stage: str,
) -> dict[str, Any]:
    """Verify a cross-stage materialization against its original attempt ledger."""

    parent_id, parent_stage, required_artifacts = dependency_contract(row, child_stage)
    if not parent_id or parent_stage == child_stage:
        return {}
    receipt_path = job_dir / (
        f".hpc_materialization__{safe_component(parent_stage)}__{safe_component(parent_id)}.json"
    )
    receipt = load_json(receipt_path)
    if (
        receipt.get("status") != "HASH_VERIFIED_MATERIALIZATION"
        or receipt.get("stage") != parent_stage
        or receipt.get("job_id") != parent_id
    ):
        raise ValueError(f"Dependency receipt identity/status mismatch: {receipt_path}")

    attempt_raw = str(receipt.get("attempt_path", "")).strip()
    attempt = resolve_package_path(
        package_root, attempt_raw, field="dependency attempt_path"
    )
    ledger_path = attempt / "FILE_MANIFEST_SHA256.csv"
    expected_ledger_hash = str(receipt.get("file_manifest_sha256", "")).strip().upper()
    if not expected_ledger_hash or not ledger_path.is_file():
        raise ValueError(f"Dependency receipt lacks its source ledger: {receipt_path}")
    if sha256_file(ledger_path) != expected_ledger_hash:
        raise ValueError(f"Dependency source ledger hash mismatch: {ledger_path}")
    with ledger_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not {"file", "size_bytes", "sha256"}.issubset(
            reader.fieldnames
        ):
            raise ValueError(f"Malformed dependency source ledger: {ledger_path}")
        ledger_rows = list(reader)
    ledger: dict[str, tuple[int, str]] = {}
    for ledger_row in ledger_rows:
        name = str(ledger_row.get("file", "")).strip()
        token = Path(name)
        if not name or token.is_absolute() or token.name != name or name in ledger:
            raise ValueError(f"Unsafe/duplicate dependency ledger entry: {name!r}")
        try:
            size = int(str(ledger_row.get("size_bytes", "")))
        except ValueError as exc:
            raise ValueError(f"Invalid dependency ledger size for {name}") from exc
        digest = str(ledger_row.get("sha256", "")).strip().upper()
        source = attempt / name
        if (
            not digest
            or not source.is_file()
            or source.stat().st_size != size
            or sha256_file(source) != digest
        ):
            raise ValueError(f"Dependency source artifact differs from ledger: {source}")
        ledger[name] = (size, digest)

    receipt_artifacts: dict[str, tuple[int, str]] = {}
    for item in receipt.get("artifacts", []):
        if not isinstance(item, Mapping):
            continue
        name = str(item.get("file", "")).strip()
        try:
            size = int(str(item.get("size_bytes", "")))
        except ValueError as exc:
            raise ValueError(f"Invalid dependency receipt size for {name}") from exc
        digest = str(item.get("sha256", "")).strip().upper()
        if not name or name in receipt_artifacts:
            raise ValueError(f"Unsafe/duplicate dependency receipt entry: {name!r}")
        receipt_artifacts[name] = (size, digest)
    for name in required_artifacts:
        if name not in ledger or receipt_artifacts.get(name) != ledger[name]:
            raise ValueError(f"Dependency artifact is not ledger-bound: {name}")
        materialized = job_dir / name
        size, digest = ledger[name]
        if (
            not materialized.is_file()
            or materialized.stat().st_size != size
            or sha256_file(materialized) != digest
        ):
            raise ValueError(f"Materialized dependency artifact is absent or altered: {materialized}")
    return {
        "parent_job_id": parent_id,
        "parent_stage": parent_stage,
        "receipt": str(receipt_path),
        "file_manifest_sha256": expected_ledger_hash,
        "artifacts": list(required_artifacts),
    }


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
