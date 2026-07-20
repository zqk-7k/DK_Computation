"""Strict readers for ORCA NMR outputs and project atom maps.

The parsers in this module never infer experimental observations.  They only
extract values that ORCA wrote to a ``.property.txt`` or ``.out`` file and
normalize atom-map tables used by the Figure 3 workflow.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


class OrcaParseError(RuntimeError):
    """Raised when a required ORCA section is absent or internally inconsistent."""


@dataclass(frozen=True)
class AtomRecord:
    index0: int
    element: str
    label: str
    component: str
    include_in_analyte_nmr: bool
    exchangeable: bool


@dataclass(frozen=True)
class CouplingRecord:
    atom_i0: int
    atom_j0: int
    j_hz: float


def _read(path: Path) -> str:
    if not path.is_file():
        raise FileNotFoundError(path)
    return path.read_text(encoding="utf-8", errors="replace")


def _truth(value: str | None, default: bool) -> bool:
    if value is None or value.strip() == "":
        return default
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "y"}:
        return True
    if normalized in {"0", "false", "no", "n"}:
        return False
    raise OrcaParseError(f"Invalid boolean value in atom map: {value!r}")


def read_atom_map(path: Path) -> list[AtomRecord]:
    """Read both legacy and explicit-solvent atom-map schemas.

    Supported index columns are ``orca_index_0based`` and
    ``atom_index_1based``.  Explicit environment atoms default to exclusion.
    """

    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise OrcaParseError(f"Atom map is empty: {path}")

    records: list[AtomRecord] = []
    for position, row in enumerate(rows):
        if row.get("orca_index_0based", "").strip():
            index0 = int(row["orca_index_0based"])
        elif row.get("atom_index_1based", "").strip():
            index0 = int(row["atom_index_1based"]) - 1
        else:
            index0 = position
        element = (row.get("element") or "").strip()
        if not re.fullmatch(r"[A-Z][a-z]?", element):
            raise OrcaParseError(f"Invalid element at atom {index0}: {element!r}")
        component = (row.get("component") or "analyte").strip() or "analyte"
        default_include = component.lower() == "analyte" and not _truth(
            row.get("environment_atom"), False
        )
        label = (
            row.get("label")
            or row.get("atom_label")
            or f"{element}{index0 + 1}"
        ).strip()
        records.append(
            AtomRecord(
                index0=index0,
                element=element,
                label=label,
                component=component,
                include_in_analyte_nmr=_truth(
                    row.get("include_in_analyte_nmr"), default_include
                ),
                exchangeable=_truth(row.get("exchangeable"), False),
            )
        )
    records.sort(key=lambda item: item.index0)
    expected = list(range(len(records)))
    observed = [item.index0 for item in records]
    if observed != expected:
        raise OrcaParseError(
            f"Atom map indices must be contiguous and zero-based after normalization: {observed[:12]}"
        )
    return records


def _shieldings_from_property(text: str) -> list[float]:
    return [
        float(value)
        for value in re.findall(
            rf'&siso\s+\[&Type\s+"Double"\]\s+({FLOAT})', text
        )
    ]


def _shieldings_from_output(text: str) -> list[tuple[int, str, float]]:
    marker = "CHEMICAL SHIELDING SUMMARY (ppm)"
    if marker not in text:
        return []
    block = text.rsplit(marker, 1)[1]
    values: list[tuple[int, str, float]] = []
    started = False
    for line in block.splitlines():
        match = re.match(
            rf"^\s*(\d+)\s+([A-Z][a-z]?)\s+({FLOAT})\s+({FLOAT})\s*$",
            line,
        )
        if match:
            started = True
            values.append((int(match.group(1)), match.group(2), float(match.group(3))))
        elif started and not line.strip():
            break
    return values


def parse_shieldings(path: Path, atom_map: Iterable[AtomRecord] | None = None) -> list[float]:
    """Parse isotropic shieldings in ORCA atom order.

    ``.property.txt`` is preferred because it preserves full precision.  The
    human-readable ``.out`` summary is accepted as a rounded fallback.
    """

    text = _read(path)
    values = _shieldings_from_property(text)
    expected_atoms = list(atom_map) if atom_map is not None else None
    if not values:
        summary = _shieldings_from_output(text)
        if not summary:
            raise OrcaParseError(f"No ORCA isotropic shielding section found: {path}")
        indices = [index for index, _, _ in summary]
        if indices != list(range(len(indices))):
            raise OrcaParseError(f"Non-contiguous shielding atom indices in {path}")
        if expected_atoms is not None:
            for (index, element, _), atom in zip(summary, expected_atoms):
                if index != atom.index0 or element != atom.element:
                    raise OrcaParseError(
                        f"Shielding/atom-map mismatch at {index}: {element} != {atom.element}"
                    )
        values = [shielding for _, _, shielding in summary]
    if expected_atoms is not None and len(values) != len(expected_atoms):
        raise OrcaParseError(
            f"Shielding count {len(values)} != atom-map count {len(expected_atoms)} for {path}"
        )
    return values


def _section(text: str, start: str, end: str) -> str:
    try:
        start_index = text.index(start)
        end_index = text.index(end, start_index + len(start))
    except ValueError as exc:
        raise OrcaParseError(f"Missing ORCA property section {start!r} or {end!r}") from exc
    return text[start_index:end_index]


def _couplings_from_property(text: str) -> list[CouplingRecord]:
    if "&pairsInfo" not in text or "&pairsTotalSSCIso" not in text:
        return []
    pairs_block = _section(text, "&pairsInfo", "&pairsDistances")
    pairs: list[tuple[int, int]] = []
    for line in pairs_block.splitlines():
        match = re.match(r"^\s*\d+\s+(\d+)\s+\d+\s+(\d+)\s+\d+\s*$", line)
        if match:
            pairs.append((int(match.group(1)), int(match.group(2))))
    values_block = _section(text, "&pairsTotalSSCIso", "&Method")
    values: list[float] = []
    for line in values_block.splitlines():
        match = re.match(rf"^\s*\d+\s+({FLOAT})\s*$", line)
        if match:
            values.append(float(match.group(1)))
    if len(pairs) != len(values):
        raise OrcaParseError(
            f"ORCA coupling pair/value mismatch: {len(pairs)} pairs, {len(values)} values"
        )
    return [CouplingRecord(i, j, value) for (i, j), value in zip(pairs, values)]


def _couplings_from_output(text: str) -> list[CouplingRecord]:
    marker = "SUMMARY OF ISOTROPIC COUPLING CONSTANTS J (Hz)"
    if marker not in text:
        return []
    block = text.rsplit(marker, 1)[1]
    current_columns: list[int] = []
    values: dict[tuple[int, int], float] = {}
    header_pattern = re.compile(r"(\d+)\s+([A-Z][a-z]?)")
    row_pattern = re.compile(rf"^\s*(\d+)\s+([A-Z][a-z]?)\s+((?:{FLOAT}\s*)+)$")
    for line in block.splitlines():
        row_match = row_pattern.match(line)
        if row_match and current_columns:
            row_index = int(row_match.group(1))
            row_values = [float(item) for item in re.findall(FLOAT, row_match.group(3))]
            for column, value in zip(current_columns, row_values):
                if row_index != column:
                    key = tuple(sorted((row_index, column)))
                    previous = values.get(key)
                    if previous is None or abs(value) > abs(previous):
                        values[key] = value
            continue
        header = header_pattern.findall(line)
        if len(header) >= 2 and not re.match(r"^\s*\d+\s+[A-Z][a-z]?", line):
            current_columns = [int(index) for index, _ in header]
    return [
        CouplingRecord(i, j, value)
        for (i, j), value in sorted(values.items())
    ]


def parse_j_couplings(path: Path) -> list[CouplingRecord]:
    """Parse total isotropic J values in hertz.

    Full-precision property data are preferred; the printed summary matrix is
    a rounded fallback and is marked as such by the caller's metadata.
    """

    text = _read(path)
    records = _couplings_from_property(text)
    if not records:
        records = _couplings_from_output(text)
    if not records:
        raise OrcaParseError(f"No ORCA isotropic J-coupling section found: {path}")
    return records


def coupling_matrix(atom_indices0: list[int], records: Iterable[CouplingRecord]) -> list[list[float]]:
    lookup = {atom: local for local, atom in enumerate(atom_indices0)}
    matrix = [[0.0 for _ in atom_indices0] for _ in atom_indices0]
    for record in records:
        if record.atom_i0 not in lookup or record.atom_j0 not in lookup:
            continue
        i = lookup[record.atom_i0]
        j = lookup[record.atom_j0]
        matrix[i][j] = matrix[j][i] = record.j_hz
    return matrix
