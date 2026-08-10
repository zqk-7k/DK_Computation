"""Shared, fail-closed helpers for the explicit-cluster GOAT -> DFT hand-off.

This module deliberately contains no subprocess or ORCA launch code.  It only
validates immutable upstream evidence, parses coordinates, writes atomic files,
and renders reviewed downstream inputs.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np


NOTICE = "CALCULATED MODEL - NOT EXPERIMENTAL DATA"
DEUTERIUM_MASS_U = 2.01410177812
HARTREE_TO_KCAL_MOL = 627.5094740631
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
XYZFILE_LINE = re.compile(
    r"^\s*\*\s*xyzfile\s+[+-]?\d+\s+\d+\s+"
    r"(?:\"(?P<double>[^\"]+)\"|'(?P<single>[^']+)'|(?P<bare>\S+))",
    flags=re.IGNORECASE,
)


class WorkflowValidationError(RuntimeError):
    """Raised whenever provenance or scientific hand-off validation fails."""


@dataclass(frozen=True)
class XYZFrame:
    symbols: tuple[str, ...]
    coordinates: tuple[tuple[float, float, float], ...]
    comment: str
    energy_hartree: float | None = None
    converged: bool | None = None
    source_frame_index_1based: int | None = None


def normalize_sha256(value: str, field: str = "sha256") -> str:
    value = str(value).strip().lower()
    if not SHA256_RE.fullmatch(value):
        raise WorkflowValidationError(f"{field} is not a 64-character SHA-256: {value!r}")
    return value


def sha256_file(path: Path) -> str:
    path = Path(path)
    if not path.is_file():
        raise WorkflowValidationError(f"Required file is missing: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def atomic_write_bytes(path: Path, data: bytes) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def atomic_write_text(path: Path, text: str) -> None:
    atomic_write_bytes(Path(path), text.encode("utf-8"))


def atomic_write_json(path: Path, payload: object) -> None:
    atomic_write_text(Path(path), json.dumps(payload, indent=2, ensure_ascii=False) + "\n")


def atomic_write_csv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row.get(key, "") for key in fields})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def read_json(path: Path) -> dict:
    path = Path(path)
    try:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError) as exc:
        raise WorkflowValidationError(f"Cannot read valid JSON from {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise WorkflowValidationError(f"JSON root must be an object: {path}")
    return payload


def read_text_auto(path: Path) -> str:
    """Read ordinary ORCA text plus PowerShell-created UTF-16 evidence safely."""

    data = Path(path).read_bytes()
    if data.startswith((b"\xff\xfe", b"\xfe\xff")):
        return data.decode("utf-16", errors="replace")
    return data.decode("utf-8", errors="replace")


def read_csv(path: Path) -> list[dict[str, str]]:
    path = Path(path)
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = [dict(row) for row in csv.DictReader(handle)]
    except OSError as exc:
        raise WorkflowValidationError(f"Cannot read CSV {path}: {exc}") from exc
    if not rows:
        raise WorkflowValidationError(f"CSV contains no data rows: {path}")
    return rows


def ensure_within(path: Path, root: Path, label: str) -> Path:
    resolved = Path(path).resolve(strict=False)
    root_resolved = Path(root).resolve(strict=False)
    try:
        resolved.relative_to(root_resolved)
    except ValueError as exc:
        raise WorkflowValidationError(f"{label} escapes project root: {resolved}") from exc
    return resolved


def source_input_files(job_directory: Path, input_name: str = "goat.inp") -> list[Path]:
    """Return exactly the source bundle recognized by ``run_goat_batch.py``."""

    root = Path(job_directory).resolve(strict=True)
    input_path = (root / input_name).resolve(strict=True)
    ensure_within(input_path, root, "GOAT input")
    referenced: list[Path] = []
    for line_number, line in enumerate(
        input_path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1
    ):
        match = XYZFILE_LINE.match(line)
        if match is None:
            continue
        token = next(value for value in match.group("double", "single", "bare") if value)
        raw = Path(token).expanduser()
        resolved = (raw if raw.is_absolute() else root / raw).resolve(strict=True)
        ensure_within(resolved, root, f"xyzfile reference at line {line_number}")
        if not resolved.is_file():
            raise WorkflowValidationError(f"xyzfile reference is not a file: {resolved}")
        referenced.append(resolved)
    unique: dict[str, Path] = {}
    for path in (input_path, *referenced):
        unique.setdefault(os.path.normcase(str(path.relative_to(root))), path)
    return list(unique.values())


def runner_input_bundle_sha256(files: Sequence[Path]) -> str:
    """Reproduce ``run_goat_batch.py``'s length-prefixed bundle hash."""

    if not files:
        raise WorkflowValidationError("Cannot fingerprint an empty input bundle")
    root = Path(files[0]).resolve().parent
    digest = hashlib.sha256()
    for path in sorted((Path(item).resolve() for item in files), key=lambda p: os.path.normcase(str(p))):
        relative = path.relative_to(root).as_posix().encode("utf-8")
        digest.update(len(relative).to_bytes(4, "big"))
        digest.update(relative)
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest()


def _parse_comment(comment: str) -> tuple[float | None, bool | None]:
    fields = comment.strip().split()
    energy: float | None = None
    if fields:
        try:
            energy = float(fields[0])
        except ValueError:
            match = re.search(r"(?:energy|E)\s*=\s*([-+0-9.eE]+)", comment)
            if match:
                energy = float(match.group(1))
    lower = comment.lower()
    converged: bool | None = None
    if "converged=true" in lower or "converged = true" in lower:
        converged = True
    elif "converged=false" in lower or "converged = false" in lower:
        converged = False
    return energy, converged


def parse_xyz_frames(path: Path) -> list[XYZFrame]:
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    frames: list[XYZFrame] = []
    cursor = 0
    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue
        try:
            atom_count = int(lines[cursor].strip())
        except ValueError as exc:
            raise WorkflowValidationError(
                f"Expected XYZ atom count at {path}:{cursor + 1}, got {lines[cursor]!r}"
            ) from exc
        if atom_count <= 0 or cursor + 2 + atom_count > len(lines):
            raise WorkflowValidationError(f"Truncated or invalid XYZ frame in {path} at line {cursor + 1}")
        comment = lines[cursor + 1].strip()
        symbols: list[str] = []
        coordinates: list[tuple[float, float, float]] = []
        for offset in range(atom_count):
            fields = lines[cursor + 2 + offset].split()
            if len(fields) < 4:
                raise WorkflowValidationError(
                    f"XYZ atom line has fewer than four fields in {path}: {fields!r}"
                )
            try:
                xyz = (float(fields[1]), float(fields[2]), float(fields[3]))
            except ValueError as exc:
                raise WorkflowValidationError(f"Invalid XYZ coordinate in {path}: {fields!r}") from exc
            if not all(math.isfinite(value) for value in xyz):
                raise WorkflowValidationError(f"Non-finite XYZ coordinate in {path}: {fields!r}")
            symbols.append(fields[0])
            coordinates.append(xyz)
        energy, converged = _parse_comment(comment)
        frames.append(
            XYZFrame(
                tuple(symbols),
                tuple(coordinates),
                comment,
                energy,
                converged,
                len(frames) + 1,
            )
        )
        cursor += atom_count + 2
    if not frames:
        raise WorkflowValidationError(f"No XYZ frames found: {path}")
    return frames


def xyz_text(frame: XYZFrame, comment: str | None = None, decimals: int = 12) -> str:
    lines = [str(len(frame.symbols)), comment if comment is not None else frame.comment]
    for symbol, (x, y, z) in zip(frame.symbols, frame.coordinates, strict=True):
        lines.append(f"{symbol:<2s} {x: .{decimals}f} {y: .{decimals}f} {z: .{decimals}f}")
    return "\n".join(lines) + "\n"


def fixed_order_kabsch_rmsd(frame_a: XYZFrame, frame_b: XYZFrame) -> float:
    """All-atom RMSD without atom permutation or reflection."""

    if frame_a.symbols != frame_b.symbols:
        raise WorkflowValidationError("Cannot compare frames with different ordered element lists")
    a = np.asarray(frame_a.coordinates, dtype=float)
    b = np.asarray(frame_b.coordinates, dtype=float)
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    u, _s, vt = np.linalg.svd(a.T @ b)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    delta = a @ rotation - b
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def frames_numerically_equal(a: XYZFrame, b: XYZFrame, tolerance_A: float = 2e-8) -> bool:
    if a.symbols != b.symbols or len(a.coordinates) != len(b.coordinates):
        return False
    return bool(
        np.allclose(
            np.asarray(a.coordinates),
            np.asarray(b.coordinates),
            atol=tolerance_A,
            rtol=0.0,
        )
    )


def validate_manifest_row_source(row: Mapping[str, str], project: Path) -> dict[str, object]:
    required = (
        "job_id",
        "figure_id",
        "scenario_id",
        "job_directory",
        "input_bundle_sha256",
        "charge",
        "multiplicity",
        "continuum_keyword",
    )
    missing = [field for field in required if not str(row.get(field, "")).strip()]
    if missing:
        raise WorkflowValidationError(f"Manifest row is missing {missing}: {row}")
    job_dir = ensure_within(Path(row["job_directory"]), project, "source GOAT directory")
    bundle_hash = runner_input_bundle_sha256(source_input_files(job_dir))
    expected = normalize_sha256(row["input_bundle_sha256"], "manifest input_bundle_sha256")
    if bundle_hash != expected:
        raise WorkflowValidationError(
            f"Current source bundle differs from manifest for {row['job_id']}: {bundle_hash} != {expected}"
        )
    return {"job_directory": job_dir, "input_bundle_sha256": bundle_hash}


def validate_completed_goat_state(
    row: Mapping[str, str], state_path: Path, project: Path
) -> dict[str, object]:
    """Validate a completed runner state against current source inputs and artifacts."""

    source = validate_manifest_row_source(row, project)
    state_path = ensure_within(Path(state_path), project, "GOAT state")
    state = read_json(state_path)
    if state.get("job_id") != row["job_id"] or state.get("figure_id") != row["figure_id"]:
        raise WorkflowValidationError(f"State identity mismatch: {state_path}")
    if state.get("status") != "completed":
        raise WorkflowValidationError(f"GOAT state is not completed: {state_path}")
    attempts = state.get("attempts")
    if not isinstance(attempts, list) or not attempts:
        raise WorkflowValidationError(f"Completed GOAT state has no attempts: {state_path}")
    attempt = attempts[-1]
    if not isinstance(attempt, dict) or attempt.get("status") != "completed":
        raise WorkflowValidationError(f"Latest GOAT attempt is not completed: {state_path}")
    try:
        return_code = int(attempt.get("return_code"))
    except (TypeError, ValueError) as exc:
        raise WorkflowValidationError(f"Invalid GOAT return code in {state_path}") from exc
    if return_code != 0 or attempt.get("source_input_unchanged") is not True:
        raise WorkflowValidationError(f"GOAT attempt failed or source input changed: {state_path}")
    current_hash = str(source["input_bundle_sha256"])
    for label, value in (
        ("state input bundle", state.get("input_bundle_sha256", "")),
        ("attempt input bundle", attempt.get("input_bundle_sha256", "")),
        ("attempt post-run input bundle", attempt.get("source_input_bundle_sha256_after", "")),
    ):
        if normalize_sha256(str(value), label) != current_hash:
            raise WorkflowValidationError(f"{label} does not match current source for {row['job_id']}")
    output_path = ensure_within(Path(str(attempt.get("output_file", ""))), project, "GOAT output")
    output_text = output_path.read_text(encoding="utf-8", errors="replace")
    if "ORCA TERMINATED NORMALLY" not in output_text:
        raise WorkflowValidationError(f"GOAT output lacks normal-termination marker: {output_path}")
    recorded_output_hash = normalize_sha256(str(attempt.get("output_sha256", "")), "output_sha256")
    if sha256_file(output_path) != recorded_output_hash:
        raise WorkflowValidationError(f"GOAT output hash mismatch: {output_path}")
    artifacts = attempt.get("artifacts")
    if not isinstance(artifacts, dict):
        raise WorkflowValidationError(f"GOAT state lacks artifact mapping: {state_path}")
    ensemble_path = ensure_within(Path(str(artifacts.get("finalensemble", ""))), project, "GOAT ensemble")
    globalminimum_value = str(artifacts.get("globalminimum", "")).strip()
    globalminimum_path = (
        ensure_within(Path(globalminimum_value), project, "GOAT global minimum")
        if globalminimum_value
        else None
    )
    if not ensemble_path.is_file() or ensemble_path.stat().st_size == 0:
        raise WorkflowValidationError(f"GOAT final ensemble is missing or empty: {ensemble_path}")
    if globalminimum_path is not None and (
        not globalminimum_path.is_file() or globalminimum_path.stat().st_size == 0
    ):
        raise WorkflowValidationError(f"GOAT global minimum is missing or empty: {globalminimum_path}")
    return {
        "state": state,
        "state_path": state_path,
        "state_sha256": sha256_file(state_path),
        "latest_attempt": attempt,
        "output_path": output_path,
        "output_sha256": recorded_output_hash,
        "ensemble_path": ensemble_path,
        "ensemble_sha256": sha256_file(ensemble_path),
        "globalminimum_path": globalminimum_path,
        "globalminimum_sha256": sha256_file(globalminimum_path) if globalminimum_path else "",
        **source,
    }


def validate_cluster_maps(
    atom_map_path: Path, isotope_map_path: Path, frame: XYZFrame
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[int, float]]:
    atom_rows = read_csv(atom_map_path)
    isotope_rows = read_csv(isotope_map_path) if Path(isotope_map_path).is_file() else []
    if len(atom_rows) != len(frame.symbols):
        raise WorkflowValidationError(
            f"atom_map length {len(atom_rows)} != coordinate atom count {len(frame.symbols)}"
        )
    for expected_index, (row, symbol) in enumerate(zip(atom_rows, frame.symbols, strict=True), start=1):
        try:
            index = int(row["atom_index_1based"])
        except (KeyError, TypeError, ValueError) as exc:
            raise WorkflowValidationError(f"Invalid atom_map index in {atom_map_path}: {row}") from exc
        if index != expected_index or row.get("element", "").strip().upper() != symbol.upper():
            raise WorkflowValidationError(
                f"atom_map order/element mismatch at atom {expected_index}: {row} versus {symbol}"
            )
        if not row.get("component", "").strip():
            raise WorkflowValidationError(f"atom_map component is blank at atom {expected_index}")
    isotope_masses: dict[int, float] = {}
    for row in isotope_rows:
        try:
            index = int(row["atom_index_1based"])
            mass = float(row["mass_u"])
        except (KeyError, TypeError, ValueError) as exc:
            raise WorkflowValidationError(f"Invalid isotope_map row in {isotope_map_path}: {row}") from exc
        if index < 1 or index > len(frame.symbols) or index in isotope_masses:
            raise WorkflowValidationError(f"Invalid or duplicate isotope atom index {index}")
        if frame.symbols[index - 1].upper() != row.get("element", "").strip().upper():
            raise WorkflowValidationError(f"isotope_map element mismatch at atom {index}")
        if frame.symbols[index - 1].upper() != "H" or not math.isclose(
            mass, DEUTERIUM_MASS_U, rel_tol=0.0, abs_tol=5e-12
        ):
            raise WorkflowValidationError(
                f"Only exact deuterium annotations are accepted; atom {index} has mass {mass}"
            )
        isotope_masses[index] = DEUTERIUM_MASS_U
    return atom_rows, isotope_rows, isotope_masses


def downstream_atom_map_rows(
    source_rows: Sequence[Mapping[str, str]],
    isotope_masses: Mapping[int, float],
    figure_id: str,
    scenario_id: str,
    species_smiles: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in source_rows:
        index = int(row["atom_index_1based"])
        include = row.get("include_in_analyte_nmr", "").strip().lower() in {"1", "true", "yes"}
        environment = row.get("environment_atom", "").strip().lower() in {"1", "true", "yes"}
        if include:
            nmr_scope = "ANALYTE"
        elif environment:
            nmr_scope = "ENVIRONMENT_EXCLUDED_FROM_ANALYTE_NMR"
        else:
            nmr_scope = "ANALYTE_EXCLUDED_FROM_ORDINARY_NMR"
        rows.append(
            {
                **row,
                "figure_id": figure_id,
                "scenario_id": scenario_id,
                "scenario_species_smiles": species_smiles,
                "nmr_scope": nmr_scope,
                "effective_isotope_mass_u": (
                    f"{isotope_masses[index]:.11f}" if index in isotope_masses else "NATURAL_ABUNDANCE"
                ),
                "source_atom_map_immutable": "true",
                "notice": NOTICE,
            }
        )
    return rows


def continuum_token(value: str) -> str:
    token = str(value).strip().upper()
    aliases = {"DMSO-D6": "DMSO", "DMSO": "DMSO", "D2O": "WATER", "H2O": "WATER", "WATER": "WATER"}
    if token not in aliases:
        raise WorkflowValidationError(f"Unsupported or unsafe CPCM solvent token: {value!r}")
    return aliases[token]


def render_opt_freq_input(
    frame: XYZFrame,
    charge: int,
    multiplicity: int,
    continuum: str,
    isotope_masses: Mapping[int, float],
    nprocs: int = 16,
    maxcore_mb: int = 6000,
    temperature_K: float | Sequence[float] = 298.15,
) -> str:
    if nprocs < 1 or nprocs > 32:
        raise WorkflowValidationError(f"nprocs must be between 1 and 32, got {nprocs}")
    if maxcore_mb < 256:
        raise WorkflowValidationError(f"maxcore_mb is implausibly small: {maxcore_mb}")
    solvent = continuum_token(continuum)
    raw_temperatures = (
        [float(temperature_K)]
        if isinstance(temperature_K, (int, float))
        else [float(value) for value in temperature_K]
    )
    temperatures = sorted(set(raw_temperatures))
    if not temperatures or len(temperatures) > 16 or any(value <= 0.0 for value in temperatures):
        raise WorkflowValidationError(
            "Thermochemistry requires 1-16 unique positive temperatures"
        )
    temperature_text = ", ".join(
        f"{value:.8f}".rstrip("0").rstrip(".") for value in temperatures
    )
    lines = [
        f"# {NOTICE}",
        "# Explicit cluster retained exactly from the approved GOAT conformer.",
        f"! PBE0 D4 def2-TZVPP def2/J RIJCOSX CPCM({solvent}) VeryTightOpt Freq VeryTightSCF DefGrid3",
        f"%pal nprocs {nprocs} end",
        f"%maxcore {maxcore_mb}",
        "%geom",
        "  MaxIter 500",
        "end",
        "%scf",
        "  MaxIter 500",
        "end",
        "%freq",
        f"  Temp {temperature_text}",
        "end",
        f"* xyz {charge} {multiplicity}",
    ]
    for index, (symbol, (x, y, z)) in enumerate(
        zip(frame.symbols, frame.coordinates, strict=True), start=1
    ):
        suffix = f" M={isotope_masses[index]:.11f}" if index in isotope_masses else ""
        lines.append(f"  {symbol:<2s} {x: .12f} {y: .12f} {z: .12f}{suffix}")
    lines.extend(("*", ""))
    text = "\n".join(lines)
    validate_rendered_mass_annotations(text, isotope_masses)
    return text


def validate_rendered_mass_annotations(text: str, isotope_masses: Mapping[int, float]) -> None:
    expected = len(isotope_masses)
    exact = text.count(f"M={DEUTERIUM_MASS_U:.11f}")
    all_annotations = len(re.findall(r"\bM\s*=", text))
    if exact != expected or all_annotations != expected:
        raise WorkflowValidationError(
            f"Rendered isotope annotation count mismatch: exact={exact}, all={all_annotations}, expected={expected}"
        )


def scientific_dependency_sha256(records: Iterable[tuple[str, str]]) -> str:
    digest = hashlib.sha256()
    for logical_name, value in sorted(records):
        name = logical_name.encode("utf-8")
        content = value.encode("utf-8")
        digest.update(len(name).to_bytes(4, "big"))
        digest.update(name)
        digest.update(len(content).to_bytes(8, "big"))
        digest.update(content)
    return digest.hexdigest()


def inline_orca_input_bundle_sha256(input_name: str, input_bytes: bytes) -> str:
    """Hash an inline-coordinate ORCA input as ``run_orca_manifest.py`` does."""

    identity = [
        {
            "role": "input",
            "logical_name": input_name,
            "size_bytes": len(input_bytes),
            "sha256": sha256_bytes(input_bytes),
        }
    ]
    return hashlib.sha256(
        json.dumps(identity, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode(
            "utf-8"
        )
    ).hexdigest()


def orca_input_bundle_sha256(
    input_name: str,
    input_bytes: bytes,
    coordinate_files: Sequence[tuple[str, bytes]] = (),
) -> str:
    """Canonical bundle hash used by ``run_orca_manifest.py``.

    ``coordinate_files`` contains the exact logical token written after
    ``xyzfile`` and its bytes.  Order follows occurrence in the input.
    """

    identity = [
        {
            "role": "input",
            "logical_name": input_name,
            "size_bytes": len(input_bytes),
            "sha256": sha256_bytes(input_bytes),
        }
    ]
    seen: set[str] = set()
    for logical_name, data in coordinate_files:
        key = logical_name.replace("\\", "/")
        if key in seen:
            continue
        seen.add(key)
        identity.append(
            {
                "role": "coordinate",
                "logical_name": key,
                "size_bytes": len(data),
                "sha256": sha256_bytes(data),
            }
        )
    return hashlib.sha256(
        json.dumps(identity, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode(
            "utf-8"
        )
    ).hexdigest()


def validate_completed_opt_freq_state(
    row: Mapping[str, str], state_path: Path, project: Path
) -> dict[str, object]:
    """Validate one completed ``run_opt_freq_serial.py`` state and artifacts."""

    project = Path(project).resolve()
    job_id = str(row.get("job_id", "")).strip()
    figure_id = str(row.get("figure_id", "")).strip()
    scenario_id = str(row.get("scenario_id", "")).strip()
    if not all((job_id, figure_id, scenario_id)):
        raise WorkflowValidationError(f"Incomplete opt/freq manifest identity: {row}")
    job_dir = ensure_within(Path(str(row.get("job_directory", ""))), project, "opt/freq job")
    input_path = ensure_within(job_dir / str(row.get("input_file", "")), project, "opt/freq input")
    input_bytes = input_path.read_bytes()
    input_hash = sha256_file(input_path)
    if input_hash != normalize_sha256(str(row.get("input_sha256", "")), "opt/freq input hash"):
        raise WorkflowValidationError(f"Current opt/freq input hash differs for {job_id}")
    bundle_hash = inline_orca_input_bundle_sha256(input_path.name, input_bytes)
    if bundle_hash != normalize_sha256(
        str(row.get("expected_input_bundle_sha256", "")), "opt/freq bundle hash"
    ):
        raise WorkflowValidationError(f"Current opt/freq bundle hash differs for {job_id}")
    provenance_path = job_dir / "provenance.json"
    provenance_hash = sha256_file(provenance_path)
    provenance = read_json(provenance_path)
    if provenance.get("job_id") != job_id or provenance.get("scenario_id") != scenario_id:
        raise WorkflowValidationError(f"Opt/freq provenance identity differs for {job_id}")
    if normalize_sha256(str(provenance.get("dft_input_sha256", ""))) != input_hash:
        raise WorkflowValidationError(f"Opt/freq provenance input hash differs for {job_id}")
    state_path = ensure_within(Path(state_path), project, "opt/freq state")
    state = read_json(state_path)
    if (
        state.get("status") != "completed"
        or state.get("job_id") != job_id
        or state.get("figure_id") != figure_id
        or state.get("scenario_id") != scenario_id
    ):
        raise WorkflowValidationError(f"Opt/freq state is not completed or identity-bound: {state_path}")
    if normalize_sha256(str(state.get("input_bundle_sha256", ""))) != bundle_hash:
        raise WorkflowValidationError(f"Opt/freq state bundle hash differs for {job_id}")
    if normalize_sha256(str(state.get("provenance_sha256", ""))) != provenance_hash:
        raise WorkflowValidationError(f"Opt/freq state provenance hash differs for {job_id}")
    attempts = state.get("attempts")
    if not isinstance(attempts, list) or not attempts or not isinstance(attempts[-1], dict):
        raise WorkflowValidationError(f"Opt/freq state has no final attempt: {state_path}")
    attempt = attempts[-1]
    if (
        attempt.get("status") != "completed"
        or int(attempt.get("return_code", 1)) != 0
        or attempt.get("normal_termination_marker") is not True
        or attempt.get("source_input_unchanged") is not True
    ):
        raise WorkflowValidationError(f"Opt/freq final attempt is not valid: {state_path}")
    if normalize_sha256(str(attempt.get("input_sha256", ""))) != input_hash:
        raise WorkflowValidationError(f"Opt/freq attempt input hash differs for {job_id}")
    if normalize_sha256(str(attempt.get("input_bundle_sha256", ""))) != bundle_hash:
        raise WorkflowValidationError(f"Opt/freq attempt bundle hash differs for {job_id}")
    resource_policy_path = project / "config" / "resource_policy.json"
    if resource_policy_path.is_file():
        resource_policy = read_json(resource_policy_path)
        locked_orca_hash = str(resource_policy.get("machine", {}).get("orca_sha256", "")).strip()
        if locked_orca_hash:
            if normalize_sha256(str(state.get("resource_policy_sha256", ""))) != sha256_file(
                resource_policy_path
            ):
                raise WorkflowValidationError(f"Opt/freq resource-policy hash differs for {job_id}")
            if normalize_sha256(str(attempt.get("orca_executable_sha256", ""))) != normalize_sha256(
                locked_orca_hash, "locked ORCA executable hash"
            ):
                raise WorkflowValidationError(f"Opt/freq ORCA executable hash differs for {job_id}")
    output_path = ensure_within(Path(str(attempt.get("output_file", ""))), project, "opt/freq output")
    if "ORCA TERMINATED NORMALLY" not in read_text_auto(output_path):
        raise WorkflowValidationError(f"Opt/freq output lacks normal termination: {output_path}")
    output_hash = sha256_file(output_path)
    if output_hash != normalize_sha256(str(attempt.get("output_sha256", ""))):
        raise WorkflowValidationError(f"Opt/freq output hash differs for {job_id}")
    artifact_payload = attempt.get("artifacts")
    if not isinstance(artifact_payload, dict):
        raise WorkflowValidationError(f"Opt/freq state lacks artifacts: {state_path}")
    artifacts: dict[str, dict[str, object]] = {}
    for role in ("gbw", "hess", "property", "optimized_xyz"):
        record = artifact_payload.get(role)
        if not isinstance(record, dict) or record.get("exists_nonempty") is not True:
            raise WorkflowValidationError(f"Opt/freq state lacks complete {role} artifact: {job_id}")
        artifact_path = ensure_within(Path(str(record.get("path", ""))), project, f"opt/freq {role}")
        artifact_hash = sha256_file(artifact_path)
        if artifact_hash != normalize_sha256(str(record.get("sha256", "")), f"{role} sha256"):
            raise WorkflowValidationError(f"Opt/freq {role} artifact hash differs for {job_id}")
        artifacts[role] = {"path": artifact_path, "sha256": artifact_hash}
    return {
        "row": dict(row),
        "state": state,
        "state_path": state_path,
        "state_sha256": sha256_file(state_path),
        "attempt": attempt,
        "job_directory": job_dir,
        "input_path": input_path,
        "input_sha256": input_hash,
        "input_bundle_sha256": bundle_hash,
        "provenance_path": provenance_path,
        "provenance_sha256": provenance_hash,
        "provenance": provenance,
        "output_path": output_path,
        "output_sha256": output_hash,
        "artifacts": artifacts,
    }
