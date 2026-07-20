from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
import uuid
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence


PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data"
DEFAULT_DFT_MANIFEST = DATA / "DFT_job_manifest.csv"
DEFAULT_DFT_AUDIT = DATA / "DFT_minimum_and_weight_audit.csv"
DEFAULT_SHIELDING_MANIFEST = DATA / "NMR_shielding_job_manifest.csv"
DEFAULT_COUPLING_MANIFEST = DATA / "NMR_scalar_coupling_job_manifest.csv"
LEGACY_HH_COUPLING_MANIFEST = DATA / "NMR_HH_coupling_job_manifest.csv"
DEFAULT_TMS_OUTPUT = (
    PROJECT / "reference" / "TMS_PBE0_pcSseg2_CPCM_DMSO" / "tms_shield.out"
)
DEFAULT_CFCL3_OUTPUT = (
    PROJECT
    / "reference"
    / "CFCl3_r2SCAN3c_PBE0_pcSseg2_CPCM_DMSO"
    / "shielding.out"
)
DEFAULT_H3PO4_OUTPUT = (
    PROJECT
    / "reference"
    / "H3PO4_r2SCAN3c_PBE0_pcSseg2_CPCM_DMSO"
    / "shielding.out"
)
DEFAULT_OUTPUT_ROOT = PROJECT / "NMR_predictions_CALCULATED_NOT_EXPERIMENTAL"

NOTICE = "CALCULATED PREDICTION — NOT EXPERIMENTAL"
SCHEMA_VERSION = "1.2.0"
HARTREE_TO_KCAL_MOL = 627.5094740631
R_KCAL_MOL_K = 0.00198720425864083
DEFAULT_SIGNIFICANT_IMAGINARY_CM = -20.0
WEIGHT_TOLERANCE = 5.0e-6

# External PBE0/pcS-2 regressions from Rusakov, Semenov & Rusakova,
# Int. J. Mol. Sci. 24, 14623 (2023), DOI 10.3390/ijms241914623, Table 2.
# That paper fitted 713 1H and 767 13C shifts across 23 natural products.  The
# present calculations use pcSseg-2 rather than pcS-2 and a different geometry/
# solvation protocol, so transfer of these coefficients is deliberately marked
# provisional and never replaces the raw same-level TMS-referenced value.
PROVISIONAL_REGRESSION: dict[str, dict[str, Any]] = {
    "H": {
        "slope": -0.90872,
        "intercept_ppm": 28.87741,
        "r_squared_in_source": 0.96807,
        "source_n_shifts": 713,
        "status": "provisional_external_regression",
        "source_basis": "pcS-2",
        "target_basis": "pcSseg-2",
        "source_doi": "10.3390/ijms241914623",
        "provenance": "Rusakov, Semenov & Rusakova, Int. J. Mol. Sci. 24, 14623 (2023), Table 2",
    },
    "C": {
        "slope": -0.93833,
        "intercept_ppm": 173.0559,
        "r_squared_in_source": 0.99674,
        "source_n_shifts": 767,
        "status": "provisional_external_regression",
        "source_basis": "pcS-2",
        "target_basis": "pcSseg-2",
        "source_doi": "10.3390/ijms241914623",
        "provenance": "Rusakov, Semenov & Rusakova, Int. J. Mol. Sci. 24, 14623 (2023), Table 2",
    },
}

ATOMIC_NUMBER = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "I": 53,
}
ELEMENT_FROM_ATOMIC_NUMBER = {number: element for element, number in ATOMIC_NUMBER.items()}
EXCHANGEABLE_H_PARENT_ELEMENTS = {"N", "O", "S"}


class ProcessingError(RuntimeError):
    """Raised when an input cannot support an auditable prediction."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_text(path: Path, text: str, encoding: str = "utf-8") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{uuid.uuid4().hex}")
    try:
        temporary.write_text(text, encoding=encoding, newline="\n")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_write_json(path: Path, value: Any) -> None:
    payload = json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    atomic_write_text(path, payload)


def atomic_write_csv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{uuid.uuid4().hex}")
    try:
        with temporary.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def read_manifest(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        missing = set(required) - fields
        if missing:
            raise ProcessingError(f"{path.name} lacks required columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ProcessingError(f"Manifest is empty: {path}")
    job_ids: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        job_id = row.get("job_id", "").strip()
        if not job_id:
            raise ProcessingError(f"{path.name}:{line_number} has no job_id")
        if job_id in job_ids:
            raise ProcessingError(f"Duplicate job_id in {path.name}: {job_id}")
        job_ids.add(job_id)
    return rows


def _manifest_suffix(path: Path, stem_prefix: str) -> str | None:
    stem = path.stem
    return stem[len(stem_prefix) :] if stem.startswith(stem_prefix) else None


def infer_dft_audit_path(
    explicit: Path | None,
    dft_manifest: Path,
    shielding_manifest: Path,
    coupling_manifest: Path,
) -> Path:
    """Find the audit paired with full or subset manifests, preferring subset suffixes."""
    if explicit is not None:
        resolved = explicit.expanduser().resolve()
        if not resolved.is_file():
            raise FileNotFoundError(resolved)
        return resolved

    suffix_sources = (
        (dft_manifest, "DFT_job_manifest"),
        (shielding_manifest, "NMR_shielding_job_manifest"),
        (coupling_manifest, "NMR_scalar_coupling_job_manifest"),
        (coupling_manifest, "NMR_HH_coupling_job_manifest"),
    )
    suffixes: list[str] = []
    for path, prefix in suffix_sources:
        suffix = _manifest_suffix(path, prefix)
        if suffix and suffix not in suffixes:
            suffixes.append(suffix)
    suffixes.append("")
    search_directories: list[Path] = []
    for directory in (
        dft_manifest.parent,
        shielding_manifest.parent,
        coupling_manifest.parent,
        DATA,
    ):
        resolved = directory.resolve()
        if resolved not in search_directories:
            search_directories.append(resolved)
    candidates = [
        directory / f"DFT_minimum_and_weight_audit{suffix}.csv"
        for suffix in suffixes
        for directory in search_directories
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError(
        "Could not infer a DFT minimum/weight audit. Tried: "
        + ", ".join(str(path) for path in candidates)
    )


def parse_audit_flag(row: dict[str, str], field: str, default: bool | None = None) -> bool:
    text = row.get(field, "").strip().casefold()
    if not text and default is not None:
        return default
    if text in {"1", "true", "yes"}:
        return True
    if text in {"0", "false", "no"}:
        return False
    raise ProcessingError(f"Invalid or missing {field} flag for audit job {row.get('job_id')}: {text!r}")


def normal_termination(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    with path.open("rb") as handle:
        handle.seek(max(0, path.stat().st_size - 131072))
        tail = handle.read().decode("utf-8", errors="replace")
    return "****ORCA TERMINATED NORMALLY****" in tail


FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def as_float(value: str) -> float:
    result = float(value.replace("D", "E").replace("d", "e"))
    if not math.isfinite(result):
        raise ValueError(f"Non-finite numeric value: {value!r}")
    return result


def last_float(text: str, pattern: str, label: str) -> float:
    matches = re.findall(pattern, text, flags=re.I | re.M)
    if not matches:
        raise ProcessingError(f"Could not parse {label}")
    return as_float(matches[-1])


def last_frequency_block(text: str) -> list[float]:
    if "VIBRATIONAL FREQUENCIES" not in text:
        return []
    block = text.rsplit("VIBRATIONAL FREQUENCIES", 1)[1]
    if "NORMAL MODES" in block:
        block = block.split("NORMAL MODES", 1)[0]
    return [
        as_float(value)
        for value in re.findall(
            rf"^\s*\d+:\s+({FLOAT_PATTERN})\s+cm\*\*-1", block, flags=re.M
        )
    ]


def parse_dft_qrrho(path: Path, significant_imaginary_cm: float) -> dict[str, Any]:
    if not normal_termination(path):
        raise ProcessingError(f"DFT job did not terminate normally: {path}")
    text = path.read_text(encoding="utf-8", errors="replace")
    qrrho_values = re.findall(r"Quasi\s+RRHO\s+\.\.\.\s+(True|False)", text, flags=re.I)
    if not qrrho_values:
        raise ProcessingError(f"No Quasi RRHO flag in {path}")
    qrrho = qrrho_values[-1].casefold() == "true"
    if not qrrho:
        raise ProcessingError(f"Quasi RRHO is not enabled in {path}")
    temperature = last_float(
        text,
        rf"THERMOCHEMISTRY\s+AT\s+({FLOAT_PATTERN})\s*K",
        "thermochemistry temperature",
    )
    gibbs = last_float(
        text,
        rf"Final\s+Gibbs\s+free\s+energy\s+\.\.\.\s+({FLOAT_PATTERN})\s+Eh",
        "QRRHO Gibbs free energy",
    )
    electronic_matches = re.findall(
        rf"FINAL\s+SINGLE\s+POINT\s+ENERGY\s+({FLOAT_PATTERN})", text, flags=re.I
    )
    electronic = as_float(electronic_matches[-1]) if electronic_matches else None
    frequencies = [value for value in last_frequency_block(text) if abs(value) > 0.5]
    negative = [value for value in frequencies if value < 0.0]
    significant = [value for value in negative if value < significant_imaginary_cm]
    return {
        "qrrho_verified": True,
        "thermochemistry_temperature_K": temperature,
        "gibbs_free_energy_hartree": gibbs,
        "electronic_energy_hartree": electronic,
        "n_vibrational_frequencies": len(frequencies),
        "n_negative_frequencies": len(negative),
        "most_negative_frequency_cm-1": min(negative) if negative else 0.0,
        "significant_imaginary": bool(significant),
        "significant_imaginary_threshold_cm-1": significant_imaginary_cm,
    }


def output_path_from_manifest(row: dict[str, str]) -> Path:
    return Path(row["job_directory"]).expanduser().resolve() / row["output_file"]


def property_path_from_output(output: Path) -> Path:
    return output.with_suffix(".property.txt")


def parse_shielding_summary(path: Path) -> list[dict[str, Any]]:
    if not normal_termination(path):
        raise ProcessingError(f"Shielding job did not terminate normally: {path}")
    text = path.read_text(encoding="utf-8", errors="replace")
    marker = "CHEMICAL SHIELDING SUMMARY (ppm)"
    position = text.rfind(marker)
    if position < 0:
        raise ProcessingError(f"No chemical shielding summary in {path}")
    summary = text[position + len(marker) :]
    row_pattern = re.compile(
        rf"^\s*(\d+)\s+([A-Z][a-z]?)\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s*$",
        flags=re.M,
    )
    rows = [
        {
            "nucleus_index_0based": int(match.group(1)),
            "element": match.group(2),
            "isotropic_shielding_ppm": as_float(match.group(3)),
            "anisotropy_ppm": as_float(match.group(4)),
        }
        for match in row_pattern.finditer(summary)
    ]
    if not rows:
        raise ProcessingError(f"Empty chemical shielding summary in {path}")
    indices = [row["nucleus_index_0based"] for row in rows]
    if len(indices) != len(set(indices)):
        raise ProcessingError(f"Duplicate nucleus indices in shielding summary: {path}")
    return rows


def extract_property_section(text: str, section_name: str) -> str:
    marker = f"${section_name}"
    start = text.rfind(marker)
    if start < 0:
        raise ProcessingError(f"Property section {marker} not found")
    end = text.find("$End", start)
    if end < 0:
        raise ProcessingError(f"Property section {marker} has no $End")
    return text[start : end + len("$End")]


def parse_property_array(
    section: str,
    name: str,
    converter: Callable[[str], int | float],
) -> list[list[int | float]]:
    declaration = re.search(
        rf"&{re.escape(name)}\s+\[[^\n]*?&Dim\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)[^\n]*\][^\n]*\n",
        section,
        flags=re.I,
    )
    if not declaration:
        raise ProcessingError(f"Property array {name} not found")
    n_rows = int(declaration.group(1))
    n_columns = int(declaration.group(2))
    rows: dict[int, list[int | float]] = {}
    for line in section[declaration.end() :].splitlines():
        stripped = line.strip()
        if stripped.startswith("&") or stripped.startswith("$End"):
            break
        tokens = stripped.split()
        if len(tokens) < n_columns + 1 or not re.fullmatch(r"\d+", tokens[0]):
            continue
        try:
            row_index = int(tokens[0])
            values = [converter(value) for value in tokens[1 : n_columns + 1]]
        except (TypeError, ValueError):
            continue
        if 0 <= row_index < n_rows:
            rows[row_index] = values
    if len(rows) != n_rows:
        raise ProcessingError(
            f"Property array {name} expected {n_rows} rows, parsed {len(rows)}"
        )
    return [rows[index] for index in range(n_rows)]


def element_from_atomic_number(atomic_number: int) -> str:
    """Return a stable element label without silently discarding unknown nuclei."""
    return ELEMENT_FROM_ATOMIC_NUMBER.get(atomic_number, f"Z{atomic_number}")


def scalar_pair_class(element_a: str, element_b: str) -> str:
    """Canonical, index-independent pair class (for example H-F, C-P or P-P)."""
    ordered = sorted(
        (element_a, element_b),
        key=lambda element: (ATOMIC_NUMBER.get(element, 10_000), element),
    )
    return "-".join(ordered)


def parse_j_pairs(property_path: Path, output_path: Path) -> list[dict[str, Any]]:
    if not normal_termination(output_path):
        raise ProcessingError(f"Coupling job did not terminate normally: {output_path}")
    if not property_path.is_file():
        raise ProcessingError(f"Missing ORCA property file: {property_path}")
    text = property_path.read_text(encoding="utf-8", errors="replace")
    section = extract_property_section(text, "SCF_Spin_Spin_Coupling")
    pairs_info = parse_property_array(section, "pairsInfo", int)
    distances = parse_property_array(section, "pairsDistances", as_float)
    couplings = parse_property_array(section, "pairsTotalSSCIso", as_float)
    if not (len(pairs_info) == len(distances) == len(couplings)):
        raise ProcessingError(f"Inconsistent J-pair arrays in {property_path}")
    result: list[dict[str, Any]] = []
    for info, distance, coupling in zip(pairs_info, distances, couplings):
        if len(info) != 4 or len(distance) != 1 or len(coupling) != 1:
            raise ProcessingError(f"Unexpected J-pair array shape in {property_path}")
        index_a, atomic_a, index_b, atomic_b = (int(value) for value in info)
        if index_b < index_a:
            index_a, index_b = index_b, index_a
            atomic_a, atomic_b = atomic_b, atomic_a
        element_a = element_from_atomic_number(atomic_a)
        element_b = element_from_atomic_number(atomic_b)
        result.append(
            {
                "nucleus_A_index_0based": index_a,
                "nucleus_A_atomic_number": atomic_a,
                "nucleus_A_element_from_property": element_a,
                "nucleus_B_index_0based": index_b,
                "nucleus_B_atomic_number": atomic_b,
                "nucleus_B_element_from_property": element_b,
                "scalar_coupling_pair_class": scalar_pair_class(element_a, element_b),
                "is_HH_pair": int(element_a == "H" and element_b == "H"),
                "distance_angstrom_as_reported_by_orca": float(distance[0]),
                "J_isotropic_Hz": float(coupling[0]),
            }
        )
    keys = [
        (row["nucleus_A_index_0based"], row["nucleus_B_index_0based"])
        for row in result
    ]
    if len(keys) != len(set(keys)):
        raise ProcessingError(f"Duplicate J pairs in {property_path}")
    return result


def _split_semicolon_values(value: str) -> list[str]:
    return [item.strip() for item in str(value).split(";") if item.strip()]


def annotate_exchangeable_hydrogens(
    atom_map: dict[int, dict[str, Any]], source: Path | str
) -> dict[int, dict[str, Any]]:
    """Derive exchangeable-H identity from explicit connectivity and cross-check flags.

    The structural assignment (H directly bound to N/O/S) is deterministic.  It does
    not imply that a peak will be observable, sharp or quantitatively integrable in a
    real solution spectrum; those observation-level properties remain provisional.
    """
    for index, atom in atom_map.items():
        attached_indices = [
            int(value)
            for value in _split_semicolon_values(atom["attached_atom_indices_0based"])
        ]
        attached_elements = _split_semicolon_values(atom["attached_elements"])
        if len(attached_indices) != len(attached_elements):
            raise ProcessingError(
                f"Connectivity field-length mismatch at atom {index} in {source}: "
                f"{attached_indices} vs {attached_elements}"
            )
        atom["attached_atom_indices_0based_parsed"] = attached_indices
        atom["attached_elements_parsed"] = attached_elements
        atom["exchangeable_H_parent_index_0based"] = None
        atom["exchangeable_H_parent_index_1based"] = None
        atom["exchangeable_H_parent_element"] = None
        atom["exchangeable_H_class"] = "not_applicable"
        atom["exchangeability_assignment_method"] = "not_applicable_non_hydrogen"
        atom["exchangeability_structural_confidence"] = "not_applicable"
        atom["exchangeable_proton_shift_status"] = "not_applicable"
        atom["exchangeable_proton_linewidth_status"] = "not_applicable"
        atom["exchangeable_proton_integral_status"] = "not_applicable"
        if atom["element"] != "H":
            if int(atom["exchangeable_H_flag"]) != 0:
                raise ProcessingError(
                    f"Non-hydrogen atom {index} carries exchangeable_H_flag=1 in {source}"
                )
            continue
        if len(attached_indices) != 1:
            raise ProcessingError(
                f"Hydrogen atom {index} must have exactly one bonded parent in {source}; "
                f"found {attached_indices}"
            )
        parent_index = attached_indices[0]
        if parent_index not in atom_map:
            raise ProcessingError(
                f"Hydrogen atom {index} references absent parent {parent_index} in {source}"
            )
        parent_element = str(atom_map[parent_index]["element"])
        if attached_elements != [parent_element]:
            raise ProcessingError(
                f"Hydrogen atom {index} parent-element mismatch in {source}: "
                f"connectivity says {attached_elements}, atom map says {parent_element}"
            )
        expected_exchangeable = int(parent_element in EXCHANGEABLE_H_PARENT_ELEMENTS)
        declared_exchangeable = int(atom["exchangeable_H_flag"])
        if declared_exchangeable != expected_exchangeable:
            raise ProcessingError(
                f"Hydrogen atom {index} exchangeable flag conflicts with connectivity in {source}: "
                f"flag={declared_exchangeable}, parent={parent_element}{parent_index + 1}"
            )
        atom["exchangeable_H_parent_index_0based"] = parent_index
        atom["exchangeable_H_parent_index_1based"] = parent_index + 1
        atom["exchangeable_H_parent_element"] = parent_element
        atom["exchangeable_H_class"] = (
            f"{parent_element}-H" if expected_exchangeable else "not_exchangeable"
        )
        atom["exchangeability_assignment_method"] = (
            "atom_map_connectivity_H_bound_to_N_O_or_S;upstream_flag_cross_checked"
        )
        atom["exchangeability_structural_confidence"] = "high"
        if expected_exchangeable:
            atom["exchangeable_proton_shift_status"] = (
                "provisional_calculated_shift_solution_exchange_not_modelled"
            )
            atom["exchangeable_proton_linewidth_status"] = (
                "provisional_display_broadening_not_predicted"
            )
            atom["exchangeable_proton_integral_status"] = (
                "provisional_one_explicit_H_nominal_area_not_quantitative"
            )
        else:
            atom["exchangeable_proton_shift_status"] = "not_exchangeable"
            atom["exchangeable_proton_linewidth_status"] = "not_exchangeable"
            atom["exchangeable_proton_integral_status"] = "not_exchangeable"
    return atom_map


def load_atom_map(figure_id: str) -> tuple[Path, dict[int, dict[str, Any]]]:
    candidates = [
        PROJECT / "02_dft_opt_freq" / figure_id / "atom_map.csv",
        PROJECT / "01_goat" / figure_id / "atom_map.csv",
    ]
    path = next((candidate for candidate in candidates if candidate.is_file()), None)
    if path is None:
        raise FileNotFoundError(f"No atom_map.csv found for {figure_id}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ProcessingError(f"Empty atom map: {path}")
    result: dict[int, dict[str, Any]] = {}
    for row in rows:
        index = int(row["orca_index_0based"])
        if index in result:
            raise ProcessingError(f"Duplicate atom index {index} in {path}")
        one_based = int(row["orca_index_1based"])
        if one_based != index + 1:
            raise ProcessingError(f"Inconsistent one-based atom index in {path}: {row}")
        result[index] = {
            "nucleus_index_0based": index,
            "nucleus_index_1based": one_based,
            "element": row["element"],
            "formal_charge": int(row["formal_charge"]),
            "canonical_symmetry_rank": int(row["canonical_symmetry_rank"]),
            "attached_atom_indices_0based": row["attached_atom_indices_0based"],
            "attached_elements": row["attached_elements"],
            "exchangeable_H_flag": int(row["exchangeable_H_flag"]),
        }
    if sorted(result) != list(range(len(result))):
        raise ProcessingError(f"Atom map indices are not contiguous in {path}")
    return path, annotate_exchangeable_hydrogens(result, path)


def xyz_elements(path: Path) -> list[str]:
    if not path.is_file():
        raise FileNotFoundError(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        raise ProcessingError(f"Malformed XYZ: {path}")
    try:
        count = int(lines[0].strip())
    except ValueError as exc:
        raise ProcessingError(f"Malformed XYZ atom count: {path}") from exc
    atoms = [line.split()[0] for line in lines[2 : 2 + count] if line.split()]
    if len(atoms) != count:
        raise ProcessingError(f"XYZ expected {count} atoms, read {len(atoms)}: {path}")
    return atoms


def validate_atom_order(xyz_path: Path, atom_map: dict[int, dict[str, Any]]) -> None:
    elements = xyz_elements(xyz_path)
    expected = [atom_map[index]["element"] for index in sorted(atom_map)]
    if elements != expected:
        raise ProcessingError(
            f"XYZ/atom-map element order mismatch: {xyz_path}; got {elements}, expected {expected}"
        )


def group_manifest_by_figure_and_conformer(
    rows: Sequence[dict[str, str]],
) -> dict[str, dict[str, dict[str, str]]]:
    grouped: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        figure_id = row["figure_id"]
        conformer_id = row["conformer_id"]
        if conformer_id in grouped[figure_id]:
            raise ProcessingError(f"Duplicate {figure_id}/{conformer_id} in property manifest")
        grouped[figure_id][conformer_id] = row
    return dict(grouped)


def prepare_dft_records(
    rows: Sequence[dict[str, str]],
    dft_audit_rows: Sequence[dict[str, str]],
    selected_figures: set[str],
    significant_imaginary_cm: float,
) -> tuple[dict[str, dict[str, dict[str, Any]]], list[dict[str, Any]]]:
    audit_by_job: dict[str, dict[str, str]] = {}
    for row in dft_audit_rows:
        if row.get("figure_id") not in selected_figures:
            continue
        job_id = row.get("job_id", "").strip()
        if not job_id:
            raise ProcessingError("DFT audit contains a selected-figure row without job_id")
        if job_id in audit_by_job:
            raise ProcessingError(f"Duplicate job_id in DFT audit: {job_id}")
        audit_by_job[job_id] = row

    grouped_raw: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["figure_id"] in selected_figures:
            grouped_raw[row["figure_id"]].append(row)
    missing = selected_figures - set(grouped_raw)
    if missing:
        raise ProcessingError(f"Selected Figure3 IDs absent from DFT manifest: {sorted(missing)}")

    grouped: dict[str, dict[str, dict[str, Any]]] = {}
    audit: list[dict[str, Any]] = []
    for figure_id, figure_rows in sorted(grouped_raw.items()):
        if len(figure_rows) > 1 and not all(
            "dft_unique_minimum" in audit_by_job.get(row["job_id"], {})
            for row in figure_rows
        ):
            raise ProcessingError(
                f"{figure_id}: multiple DFT rows require an audit with dft_unique_minimum; "
                "regenerate the post-DFT minimum/weight audit before processing"
            )
        parsed_valid: list[dict[str, Any]] = []
        all_records: dict[str, dict[str, Any]] = {}
        for row in figure_rows:
            conformer_id = row["job_id"]
            audit_row = audit_by_job.get(conformer_id)
            if audit_row is None:
                raise ProcessingError(
                    f"{figure_id}/{conformer_id} is absent from the inferred DFT minimum/weight audit"
                )
            if audit_row.get("figure_id") != figure_id:
                raise ProcessingError(
                    f"DFT audit figure mismatch for {conformer_id}: {audit_row.get('figure_id')} vs {figure_id}"
                )
            output = output_path_from_manifest(row)
            audit_has_unique_flag = "dft_unique_minimum" in audit_row
            unique_minimum = parse_audit_flag(
                audit_row, "dft_unique_minimum", default=True
            )
            selected_for_nmr = parse_audit_flag(
                audit_row, "selected_for_nmr", default=False
            )
            record: dict[str, Any] = {
                **row,
                "conformer_id": conformer_id,
                "dft_output": str(output),
                "dft_output_sha256": sha256_file(output) if output.is_file() else None,
                "dft_audit_parse_status": audit_row.get("parse_status"),
                "dft_audit_has_unique_minimum_field": int(audit_has_unique_flag),
                "dft_unique_minimum": int(unique_minimum),
                "dft_duplicate_of": audit_row.get("dft_duplicate_of", ""),
                "dft_duplicate_rmsd_angstrom": audit_row.get(
                    "dft_duplicate_rmsd_angstrom", ""
                ),
                "dft_duplicate_delta_g_kcal_mol": audit_row.get(
                    "dft_duplicate_delta_g_kcal_mol", ""
                ),
                "selected_for_nmr": int(selected_for_nmr),
                "audit_boltzmann_weight_all_valid": audit_row.get(
                    "boltzmann_weight_all_valid", ""
                ),
                "audit_boltzmann_weight_selected_normalized": audit_row.get(
                    "boltzmann_weight_selected_normalized", ""
                ),
            }
            try:
                parsed = parse_dft_qrrho(output, significant_imaginary_cm)
                record.update(parsed)
                if audit_row.get("parse_status", "").strip().casefold() != "parsed":
                    record["weighting_status"] = "excluded_by_dft_audit_parse_status"
                elif not unique_minimum:
                    record["weighting_status"] = "excluded_post_dft_duplicate"
                elif parsed["significant_imaginary"]:
                    record["weighting_status"] = "excluded_significant_imaginary_frequency"
                else:
                    audit_significant_text = audit_row.get("significant_imaginary", "").strip()
                    if audit_significant_text:
                        audit_significant = parse_audit_flag(
                            audit_row, "significant_imaginary"
                        )
                        if audit_significant != bool(parsed["significant_imaginary"]):
                            raise ProcessingError(
                                f"{figure_id}/{conformer_id}: DFT audit significant-imaginary flag "
                                "does not match the parsed ORCA output"
                            )
                    audit_weight_text = audit_row.get(
                        "boltzmann_weight_all_valid", ""
                    ).strip()
                    if not audit_weight_text:
                        raise ProcessingError(
                            f"{figure_id}/{conformer_id}: unique valid minimum lacks "
                            "boltzmann_weight_all_valid in the DFT audit"
                        )
                    audit_weight = as_float(audit_weight_text)
                    if audit_weight < 0.0:
                        raise ProcessingError(
                            f"{figure_id}/{conformer_id}: negative audit Boltzmann weight"
                        )
                    record["audit_boltzmann_weight_all_valid"] = audit_weight
                    if selected_for_nmr:
                        selected_text = audit_row.get(
                            "boltzmann_weight_selected_normalized", ""
                        ).strip()
                        if not selected_text:
                            raise ProcessingError(
                                f"{figure_id}/{conformer_id}: selected_for_nmr=1 but selected-normalized "
                                "weight is absent from the DFT audit"
                            )
                        record["audit_boltzmann_weight_selected_normalized"] = as_float(
                            selected_text
                        )
                    record["weighting_status"] = "included"
                    parsed_valid.append(record)
            except Exception as exc:
                record["weighting_status"] = "excluded_parse_or_completion_failure"
                record["processing_error"] = repr(exc)
            if conformer_id in all_records:
                raise ProcessingError(f"Duplicate DFT conformer ID: {conformer_id}")
            all_records[conformer_id] = record
            audit.append(record)
        if not parsed_valid:
            raise ProcessingError(f"No valid QRRHO minima for {figure_id}")
        temperatures = [float(record["thermochemistry_temperature_K"]) for record in parsed_valid]
        if max(temperatures) - min(temperatures) > 1.0e-3:
            raise ProcessingError(f"Inconsistent thermochemistry temperatures for {figure_id}")
        temperature = statistics.fmean(temperatures)
        minimum = min(float(record["gibbs_free_energy_hartree"]) for record in parsed_valid)
        terms: list[float] = []
        for record in parsed_valid:
            relative = (
                float(record["gibbs_free_energy_hartree"]) - minimum
            ) * HARTREE_TO_KCAL_MOL
            record["relative_qrrho_gibbs_kcal_mol"] = relative
            terms.append(math.exp(-relative / (R_KCAL_MOL_K * temperature)))
        partition = sum(terms)
        for record, term in zip(parsed_valid, terms):
            record["boltzmann_weight_all_valid_qrrho"] = term / partition
            audit_weight = float(record["audit_boltzmann_weight_all_valid"])
            if abs(record["boltzmann_weight_all_valid_qrrho"] - audit_weight) > WEIGHT_TOLERANCE:
                raise ProcessingError(
                    f"{figure_id}/{record['conformer_id']}: recomputed all-minimum QRRHO weight "
                    f"{record['boltzmann_weight_all_valid_qrrho']:.9g} differs from DFT audit "
                    f"{audit_weight:.9g}"
                )
        audit_weight_sum = sum(
            float(record["audit_boltzmann_weight_all_valid"]) for record in parsed_valid
        )
        if abs(audit_weight_sum - 1.0) > WEIGHT_TOLERANCE:
            raise ProcessingError(
                f"{figure_id}: unique-valid DFT audit weights sum to {audit_weight_sum:.9g}, not 1"
            )
        selected_valid = [record for record in parsed_valid if record["selected_for_nmr"]]
        if not selected_valid:
            raise ProcessingError(f"DFT audit selects no unique minimum for NMR: {figure_id}")
        selected_population = sum(
            float(record["boltzmann_weight_all_valid_qrrho"])
            for record in selected_valid
        )
        for record in selected_valid:
            expected_selected = (
                float(record["boltzmann_weight_all_valid_qrrho"]) / selected_population
            )
            supplied_selected = float(
                record["audit_boltzmann_weight_selected_normalized"]
            )
            if abs(expected_selected - supplied_selected) > WEIGHT_TOLERANCE:
                raise ProcessingError(
                    f"{figure_id}/{record['conformer_id']}: DFT audit selected-normalized weight "
                    f"{supplied_selected:.9g} differs from recomputed value {expected_selected:.9g}"
                )
            record["dft_audit_validation_status"] = (
                "validated_unique_minimum_and_weights"
            )
        grouped[figure_id] = all_records
    return grouped, audit


def validate_property_manifest_weights(
    figure_id: str,
    manifest_rows: dict[str, dict[str, str]],
    dft_records: dict[str, dict[str, Any]],
    warnings: list[str],
) -> float:
    selected_population = 0.0
    for conformer_id, row in manifest_rows.items():
        if conformer_id not in dft_records:
            raise ProcessingError(f"{figure_id}/{conformer_id} absent from DFT manifest")
        dft = dft_records[conformer_id]
        if dft.get("weighting_status") != "included":
            raise ProcessingError(
                f"{figure_id}/{conformer_id} is in a property manifest but lacks valid QRRHO Gibbs energy"
            )
        if int(dft.get("dft_unique_minimum", 0)) != 1:
            raise ProcessingError(
                f"{figure_id}/{conformer_id} is in a property manifest but the DFT audit marks it as a post-DFT duplicate"
            )
        if int(dft.get("selected_for_nmr", 0)) != 1:
            raise ProcessingError(
                f"{figure_id}/{conformer_id} is in a property manifest but selected_for_nmr is not 1 in the DFT audit"
            )
        selected_population += float(dft["boltzmann_weight_all_valid_qrrho"])
    if selected_population <= 0.0:
        raise ProcessingError(f"No positive selected population for {figure_id}")
    for conformer_id, row in manifest_rows.items():
        expected = float(dft_records[conformer_id]["boltzmann_weight_all_valid_qrrho"]) / selected_population
        supplied_text = row.get("boltzmann_weight", "").strip()
        if supplied_text:
            supplied = as_float(supplied_text)
            if abs(expected - supplied) > WEIGHT_TOLERANCE:
                raise ProcessingError(
                    f"{figure_id}/{conformer_id}: manifest Boltzmann weight {supplied:.9g} "
                    f"differs from recomputed QRRHO-selected weight {expected:.9g}"
                )
            audit_selected = float(
                dft["audit_boltzmann_weight_selected_normalized"]
            )
            if abs(supplied - audit_selected) > WEIGHT_TOLERANCE:
                raise ProcessingError(
                    f"{figure_id}/{conformer_id}: manifest Boltzmann weight {supplied:.9g} "
                    f"differs from DFT audit selected-normalized weight {audit_selected:.9g}"
                )
    return selected_population


def completed_weight_context(
    dft_records: dict[str, dict[str, Any]],
    selected_ids: set[str],
    completed_ids: set[str],
) -> dict[str, Any]:
    selected_population = sum(
        float(dft_records[conformer_id]["boltzmann_weight_all_valid_qrrho"])
        for conformer_id in selected_ids
    )
    completed_population = sum(
        float(dft_records[conformer_id]["boltzmann_weight_all_valid_qrrho"])
        for conformer_id in completed_ids
    )
    if completed_population <= 0.0:
        raise ProcessingError("Completed property jobs have zero QRRHO population")
    normalized = {
        conformer_id: float(dft_records[conformer_id]["boltzmann_weight_all_valid_qrrho"])
        / completed_population
        for conformer_id in completed_ids
    }
    return {
        "selected_population_of_all_valid_dft": selected_population,
        "completed_population_of_all_valid_dft": completed_population,
        "completed_fraction_of_selected_population": completed_population / selected_population,
        "normalized_completed_weights": normalized,
    }


def weighted_mean_sd(values: Sequence[float], weights: Sequence[float]) -> tuple[float, float]:
    if not values or len(values) != len(weights):
        raise ProcessingError("Weighted statistics received inconsistent arrays")
    total = sum(weights)
    if total <= 0.0:
        raise ProcessingError("Weighted statistics received non-positive total weight")
    normalized = [weight / total for weight in weights]
    mean = sum(value * weight for value, weight in zip(values, normalized))
    variance = sum(weight * (value - mean) ** 2 for value, weight in zip(values, normalized))
    return mean, math.sqrt(max(0.0, variance))


def _reference_record(
    output: Path | None,
    element: str,
    compound: str,
    convention: str,
    limitation: str,
    required: bool,
) -> dict[str, Any]:
    """Parse one same-level reference, or retain an explicit auditable pending state."""
    base: dict[str, Any] = {
        "notice": NOTICE,
        "element": element,
        "reference_compound": compound,
        "reference_output": str(output) if output is not None else None,
        "reference_output_sha256": (
            sha256_file(output) if output is not None and output.is_file() else None
        ),
        "reference_shift_convention_ppm": 0.0,
        "reference_convention": convention,
        "reference_limitation": limitation,
        "raw_shift_equation": (
            "delta_raw_ppm = mean(reference shielding) - sample shielding "
            "+ reference_shift_convention_ppm"
        ),
    }
    if output is None or not output.is_file():
        if required:
            raise ProcessingError(f"Required {element} reference output is missing: {output}")
        return {
            **base,
            "status": "pending_reference_calculation",
            "mean_isotropic_shielding_ppm": None,
            "equivalent_atom_population_sd_ppm": None,
            "minimum_isotropic_shielding_ppm": None,
            "maximum_isotropic_shielding_ppm": None,
            "n_reference_nuclei": 0,
            "nucleus_indices_0based": [],
            "processing_error": "reference output does not exist",
        }
    try:
        rows = parse_shielding_summary(output)
        matching = [row for row in rows if row["element"] == element]
        values = [float(row["isotropic_shielding_ppm"]) for row in matching]
        if not values:
            raise ProcessingError(f"Reference lacks {element} shielding rows")
    except (OSError, ValueError, ProcessingError) as exc:
        if required:
            raise ProcessingError(f"Invalid required {element} reference {output}: {exc}") from exc
        return {
            **base,
            "status": "invalid_reference_output",
            "mean_isotropic_shielding_ppm": None,
            "equivalent_atom_population_sd_ppm": None,
            "minimum_isotropic_shielding_ppm": None,
            "maximum_isotropic_shielding_ppm": None,
            "n_reference_nuclei": 0,
            "nucleus_indices_0based": [],
            "processing_error": repr(exc),
        }
    return {
        **base,
        "status": "available_calculated_same_level_reference",
        "mean_isotropic_shielding_ppm": statistics.fmean(values),
        "equivalent_atom_population_sd_ppm": statistics.pstdev(values),
        "minimum_isotropic_shielding_ppm": min(values),
        "maximum_isotropic_shielding_ppm": max(values),
        "n_reference_nuclei": len(values),
        "nucleus_indices_0based": [row["nucleus_index_0based"] for row in matching],
        "processing_error": None,
    }


def reference_summary(
    tms_output: Path,
    cfcl3_output: Path | None = None,
    h3po4_output: Path | None = None,
) -> dict[str, dict[str, Any]]:
    return {
        "H": _reference_record(
            tms_output,
            "H",
            "tetramethylsilane (TMS)",
            "TMS is assigned delta = 0 ppm for 1H",
            "Continuum DMSO omits explicit solute-solvent and concentration effects.",
            required=True,
        ),
        "C": _reference_record(
            tms_output,
            "C",
            "tetramethylsilane (TMS)",
            "TMS is assigned delta = 0 ppm for 13C",
            "Continuum DMSO omits explicit solute-solvent and concentration effects.",
            required=True,
        ),
        "F": _reference_record(
            cfcl3_output,
            "F",
            "trichlorofluoromethane (CFCl3)",
            "CFCl3 is assigned delta = 0 ppm for 19F",
            "The isolated CFCl3/CPCM(DMSO) model does not reproduce a physical internal or external reference environment exactly.",
            required=False,
        ),
        "P": _reference_record(
            h3po4_output,
            "P",
            "orthophosphoric acid (H3PO4)",
            "The conventional 31P zero is an external 85% aqueous H3PO4 reference",
            "An isolated neutral H3PO4 molecule in CPCM(DMSO) is not the experimental 85% aqueous external standard; hydrogen bonding, ionization, water content, concentration and external-susceptibility effects are not represented.",
            required=False,
        ),
    }


def method_line(path: Path) -> str | None:
    if not path.is_file():
        return None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.lstrip().startswith("!"):
            return line.strip()
    return None


def file_record(path: Path) -> dict[str, Any]:
    return {
        "path": str(path),
        "exists": path.is_file(),
        "sha256": sha256_file(path) if path.is_file() else None,
        "size_bytes": path.stat().st_size if path.is_file() else None,
    }


def format_optional_number(value: Any, precision: int = 6) -> str:
    return "PENDING" if value is None else f"{float(value):.{precision}f}"


def build_human_report(bundle: dict[str, Any]) -> str:
    compound = bundle["compound"]
    weighting = bundle["weighting"]
    lines = [
        NOTICE,
        "=" * len(NOTICE),
        "",
        f"Compound: {compound['figure_id']} — {compound['display_name']}",
        f"Formula: {compound['formula']}",
        f"InChIKey: {compound['inchikey']}",
        "",
        "This file contains calculated ORCA predictions. It is not an experimental NMR report,",
        "does not establish identity or purity, and must not be presented as measured characterization.",
        "",
        "QRRHO conformer weighting",
        "--------------------------",
        f"Temperature: {weighting['temperature_K']:.2f} K",
        f"Unique valid DFT minima: {weighting['n_valid_dft_conformers']}",
        f"Post-DFT duplicates excluded: {weighting['n_post_dft_duplicates_excluded']}",
        f"Shielding-manifest population coverage: {weighting['shielding_selected_population_of_all_valid_dft']:.8f}",
        f"Completed shielding population coverage: {weighting['shielding_completed_population_of_all_valid_dft']:.8f}",
        "",
        "Calculated chemical shifts",
        "---------------------------",
        "Raw values use the same-level nucleus-specific reference (TMS for H/C, CFCl3 for F, H3PO4 for P).",
        "F/P values remain PENDING when their reference output is absent or invalid; H/C regression is never applied to F/P.",
        "H/C regression values use an external PBE0/pcS-2 fit (DOI 10.3390/ijms241914623, Table 2)",
        "and are provisional here because this workflow uses pcSseg-2 and a different protocol.",
        "The reported conformer SD is population spread, not predictive error.",
        "",
        "Atom       sigma_avg / ppm   delta_reference_raw / ppm   delta_HC_regression_provisional / ppm",
    ]
    for row in bundle["averaged_nuclei"]:
        exchange_note = (
            " [exchangeable; shift/linewidth/integral provisional]"
            if int(row.get("exchangeable_H_flag", 0))
            else ""
        )
        lines.append(
            f"{row['atom_label_1based']:<10s} {row['boltzmann_averaged_shielding_ppm']:>15.6f} "
            f"{format_optional_number(row['raw_reference_shift_ppm']):>27s} "
            f"{format_optional_number(row['provisional_external_regression_shift_ppm']):>39s}"
            f"{exchange_note}"
        )
    lines.extend(
        [
            "",
            "Calculated scalar couplings",
            "---------------------------",
            "These are Boltzmann-averaged calculated values, not measured splittings.",
            "The separate legacy JHH tables are strict H-H filtered compatibility views.",
            "",
            "Pair                 class      J_avg / Hz   conformer_SD / Hz",
        ]
    )
    for row in bundle["averaged_scalar_couplings"]:
        lines.append(
            f"{row['pair_label_1based']:<20s} {row['scalar_coupling_pair_class']:<8s} "
            f"{row['boltzmann_averaged_J_Hz']:>10.6f} "
            f"{row['conformer_population_sd_J_Hz']:>19.6f}"
        )
    lines.extend(["", "Assumptions and limitations", "---------------------------"])
    lines.extend(f"- {item}" for item in bundle["assumptions_and_limitations"])
    if bundle["warnings"]:
        lines.extend(["", "Warnings", "--------"])
        lines.extend(f"- {item}" for item in bundle["warnings"])
    return "\n".join(lines) + "\n"


def process_figure(
    figure_id: str,
    dft_records: dict[str, dict[str, Any]],
    shielding_manifest: dict[str, dict[str, str]],
    coupling_manifest: dict[str, dict[str, str]],
    references: dict[str, dict[str, Any]],
    tms_output: Path,
    cfcl3_output: Path,
    h3po4_output: Path,
    output_root: Path,
    allow_partial: bool,
    overwrite: bool,
    source_manifest_paths: dict[str, Path],
) -> dict[str, Any]:
    output_directory = output_root / figure_id
    if output_directory.exists() and any(output_directory.iterdir()) and not overwrite:
        raise ProcessingError(
            f"Output directory already contains files: {output_directory}; pass --overwrite to replace this script's files"
        )

    warnings: list[str] = []
    atom_map_path, atom_map = load_atom_map(figure_id)
    if not shielding_manifest:
        raise ProcessingError(f"No shielding jobs in manifest for {figure_id}")
    if not coupling_manifest:
        raise ProcessingError(f"No coupling jobs in manifest for {figure_id}")

    shielding_selected_population = validate_property_manifest_weights(
        figure_id, shielding_manifest, dft_records, warnings
    )
    coupling_selected_population = validate_property_manifest_weights(
        figure_id, coupling_manifest, dft_records, warnings
    )

    parsed_shieldings: dict[str, list[dict[str, Any]]] = {}
    parsed_couplings: dict[str, list[dict[str, Any]]] = {}
    source_files: list[dict[str, Any]] = [
        file_record(atom_map_path),
        file_record(tms_output),
        file_record(cfcl3_output),
        file_record(h3po4_output),
    ]

    for conformer_id, row in sorted(shielding_manifest.items()):
        output = output_path_from_manifest(row)
        xyz = Path(row["job_directory"]) / "conformer.xyz"
        try:
            validate_atom_order(xyz, atom_map)
            parsed = parse_shielding_summary(output)
            for shielding in parsed:
                index = int(shielding["nucleus_index_0based"])
                if index not in atom_map:
                    raise ProcessingError(f"Shielding nucleus {index} absent from atom map")
                if shielding["element"] != atom_map[index]["element"]:
                    raise ProcessingError(
                        f"Shielding element mismatch at atom {index}: {shielding['element']} vs {atom_map[index]['element']}"
                    )
            parsed_shieldings[conformer_id] = parsed
            source_files.extend([file_record(output), file_record(xyz)])
        except Exception as exc:
            message = f"{figure_id}/{conformer_id} shielding unavailable: {exc!r}"
            warnings.append(message)
            if not allow_partial:
                raise ProcessingError(message) from exc

    for conformer_id, row in sorted(coupling_manifest.items()):
        output = output_path_from_manifest(row)
        property_path = property_path_from_output(output)
        xyz = Path(row["job_directory"]) / "conformer.xyz"
        try:
            validate_atom_order(xyz, atom_map)
            parsed = parse_j_pairs(property_path, output)
            for pair in parsed:
                for side in ("A", "B"):
                    index = int(pair[f"nucleus_{side}_index_0based"])
                    atomic_number = int(pair[f"nucleus_{side}_atomic_number"])
                    if index not in atom_map:
                        raise ProcessingError(f"J-pair nucleus {index} absent from atom map")
                    expected_number = ATOMIC_NUMBER.get(atom_map[index]["element"])
                    if expected_number != atomic_number:
                        raise ProcessingError(
                            f"J-pair atomic-number mismatch at atom {index}: {atomic_number} vs {expected_number}"
                        )
            parsed_couplings[conformer_id] = parsed
            source_files.extend([file_record(output), file_record(property_path), file_record(xyz)])
        except Exception as exc:
            message = f"{figure_id}/{conformer_id} coupling unavailable: {exc!r}"
            warnings.append(message)
            if not allow_partial:
                raise ProcessingError(message) from exc

    if not parsed_shieldings:
        raise ProcessingError(f"No completed shielding results for {figure_id}")
    if not parsed_couplings:
        raise ProcessingError(f"No completed coupling results for {figure_id}")

    shield_context = completed_weight_context(
        dft_records, set(shielding_manifest), set(parsed_shieldings)
    )
    coupling_context = completed_weight_context(
        dft_records, set(coupling_manifest), set(parsed_couplings)
    )
    if shield_context["completed_fraction_of_selected_population"] < 1.0 - 1.0e-10:
        warnings.append("Shieldings were renormalized over an incomplete selected-conformer subset")
    if coupling_context["completed_fraction_of_selected_population"] < 1.0 - 1.0e-10:
        warnings.append("J couplings were renormalized over an incomplete selected-conformer subset")

    shielding_raw_rows: list[dict[str, Any]] = []
    shielding_by_nucleus: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for conformer_id, rows in sorted(parsed_shieldings.items()):
        dft = dft_records[conformer_id]
        normalized_weight = shield_context["normalized_completed_weights"][conformer_id]
        for row in rows:
            index = int(row["nucleus_index_0based"])
            atom = atom_map[index]
            combined = {
                "notice": NOTICE,
                "figure_id": figure_id,
                "conformer_id": conformer_id,
                "relative_qrrho_gibbs_kcal_mol": dft["relative_qrrho_gibbs_kcal_mol"],
                "boltzmann_weight_all_valid_qrrho": dft["boltzmann_weight_all_valid_qrrho"],
                "boltzmann_weight_completed_shielding_subset": normalized_weight,
                **atom,
                "atom_label_1based": f"{atom['element']}{atom['nucleus_index_1based']}",
                "isotropic_shielding_ppm": row["isotropic_shielding_ppm"],
                "anisotropy_ppm": row["anisotropy_ppm"],
                "source_output": str(output_path_from_manifest(shielding_manifest[conformer_id])),
            }
            shielding_raw_rows.append(combined)
            shielding_by_nucleus[index].append(combined)

    averaged_nuclei: list[dict[str, Any]] = []
    for index, rows in sorted(shielding_by_nucleus.items()):
        atom = atom_map[index]
        element = atom["element"]
        if element not in references:
            warnings.append(
                f"No reference model configured for {element} at atom {index}; shielding retained only in the by-conformer table"
            )
            continue
        values = [float(row["isotropic_shielding_ppm"]) for row in rows]
        weights = [
            float(row["boltzmann_weight_completed_shielding_subset"]) for row in rows
        ]
        mean, spread = weighted_mean_sd(values, weights)
        anisotropy_mean, anisotropy_spread = weighted_mean_sd(
            [float(row["anisotropy_ppm"]) for row in rows], weights
        )
        reference = references[element]
        regression = PROVISIONAL_REGRESSION.get(element)
        reference_available = (
            reference.get("status") == "available_calculated_same_level_reference"
            and reference.get("mean_isotropic_shielding_ppm") is not None
        )
        raw_shift = (
            float(reference["mean_isotropic_shielding_ppm"])
            - mean
            + float(reference.get("reference_shift_convention_ppm", 0.0))
            if reference_available
            else None
        )
        provisional_shift = (
            float(regression["slope"]) * mean + float(regression["intercept_ppm"])
            if regression is not None
            else None
        )
        flags = [NOTICE]
        if regression is not None:
            flags.append("PROVISIONAL_EXTERNAL_PCS2_REGRESSION_APPLIED_TO_PCSSEG2")
        else:
            flags.append("NO_HC_REGRESSION_APPLIED")
        if not reference_available:
            flags.append(f"RAW_{element}_SHIFT_PENDING_OR_INVALID_REFERENCE")
            warnings.append(
                f"{element}{atom['nucleus_index_1based']} retains calculated shielding but has no raw shift: "
                f"{reference.get('status')} ({reference.get('reference_output')})"
            )
        if element == "P":
            flags.append("H3PO4_EXTERNAL_REFERENCE_AND_SOLVENT_MODEL_LIMITATION")
        if int(atom["exchangeable_H_flag"]):
            flags.append("EXCHANGEABLE_H_SOLVENT_CONCENTRATION_AND_EXCHANGE_SENSITIVE")
        if shield_context["completed_population_of_all_valid_dft"] < 0.9999:
            flags.append("TRUNCATED_OR_INCOMPLETE_CONFORMER_POPULATION")
        averaged_nuclei.append(
            {
                "notice": NOTICE,
                "figure_id": figure_id,
                **atom,
                "atom_label_1based": f"{element}{atom['nucleus_index_1based']}",
                "n_conformers_with_shielding": len(rows),
                "completed_population_of_all_valid_dft": shield_context[
                    "completed_population_of_all_valid_dft"
                ],
                "boltzmann_averaged_shielding_ppm": mean,
                "conformer_population_sd_shielding_ppm": spread,
                "minimum_conformer_shielding_ppm": min(values),
                "maximum_conformer_shielding_ppm": max(values),
                "boltzmann_averaged_anisotropy_ppm": anisotropy_mean,
                "conformer_population_sd_anisotropy_ppm": anisotropy_spread,
                "reference_status": reference["status"],
                "reference_compound": reference["reference_compound"],
                "reference_output": reference["reference_output"],
                "reference_output_sha256": reference["reference_output_sha256"],
                "reference_convention": reference["reference_convention"],
                "reference_limitation": reference["reference_limitation"],
                "reference_mean_shielding_ppm": reference[
                    "mean_isotropic_shielding_ppm"
                ],
                "reference_equivalent_atom_sd_ppm": reference[
                    "equivalent_atom_population_sd_ppm"
                ],
                "raw_reference_shift_ppm": raw_shift,
                "raw_reference_shift_conformer_population_sd_ppm": (
                    spread if reference_available else None
                ),
                "raw_shift_equation": reference["raw_shift_equation"],
                "tms_reference_mean_shielding_ppm": (
                    reference["mean_isotropic_shielding_ppm"]
                    if element in {"H", "C"}
                    else None
                ),
                "tms_reference_equivalent_atom_sd_ppm": (
                    reference["equivalent_atom_population_sd_ppm"]
                    if element in {"H", "C"}
                    else None
                ),
                "raw_tms_reference_shift_ppm": (
                    raw_shift if element in {"H", "C"} else None
                ),
                "raw_tms_shift_conformer_population_sd_ppm": (
                    spread if element in {"H", "C"} and reference_available else None
                ),
                "provisional_regression_slope": (
                    regression["slope"] if regression is not None else None
                ),
                "provisional_regression_intercept_ppm": (
                    regression["intercept_ppm"] if regression is not None else None
                ),
                "provisional_regression_source_r_squared": (
                    regression["r_squared_in_source"] if regression is not None else None
                ),
                "provisional_regression_source_n_shifts": (
                    regression["source_n_shifts"] if regression is not None else None
                ),
                "provisional_regression_source_basis": (
                    regression["source_basis"] if regression is not None else None
                ),
                "provisional_regression_target_basis": (
                    regression["target_basis"] if regression is not None else None
                ),
                "provisional_regression_source_doi": (
                    regression["source_doi"] if regression is not None else None
                ),
                "provisional_external_regression_shift_ppm": provisional_shift,
                "provisional_regression_conformer_population_sd_ppm": (
                    abs(float(regression["slope"])) * spread
                    if regression is not None
                    else None
                ),
                "regression_provenance": (
                    regression["provenance"] if regression is not None else None
                ),
                "predictive_uncertainty_status": (
                    "not quantified for this workflow; any reported conformer SD is population spread only. "
                    "H/C regression transfer and F/P reference-model errors are not validated uncertainty estimates"
                ),
                "flags": ";".join(flags),
            }
        )

    coupling_raw_rows: list[dict[str, Any]] = []
    coupling_by_pair: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for conformer_id, rows in sorted(parsed_couplings.items()):
        dft = dft_records[conformer_id]
        normalized_weight = coupling_context["normalized_completed_weights"][conformer_id]
        for row in rows:
            index_a = int(row["nucleus_A_index_0based"])
            index_b = int(row["nucleus_B_index_0based"])
            atom_a = atom_map[index_a]
            atom_b = atom_map[index_b]
            pair_class = scalar_pair_class(atom_a["element"], atom_b["element"])
            touches_exchangeable_h = int(
                int(atom_a["exchangeable_H_flag"]) == 1
                or int(atom_b["exchangeable_H_flag"]) == 1
            )
            combined = {
                "notice": NOTICE,
                "figure_id": figure_id,
                "conformer_id": conformer_id,
                "relative_qrrho_gibbs_kcal_mol": dft["relative_qrrho_gibbs_kcal_mol"],
                "boltzmann_weight_all_valid_qrrho": dft["boltzmann_weight_all_valid_qrrho"],
                "boltzmann_weight_completed_coupling_subset": normalized_weight,
                "nucleus_A_index_0based": index_a,
                "nucleus_A_index_1based": index_a + 1,
                "nucleus_A_element": atom_a["element"],
                "nucleus_B_index_0based": index_b,
                "nucleus_B_index_1based": index_b + 1,
                "nucleus_B_element": atom_b["element"],
                "pair_label_1based": (
                    f"{atom_a['element']}{index_a + 1}-{atom_b['element']}{index_b + 1}"
                ),
                "scalar_coupling_pair_class": pair_class,
                "is_HH_pair": int(pair_class == "H-H"),
                "nucleus_A_exchangeable_H_flag": int(atom_a["exchangeable_H_flag"]),
                "nucleus_B_exchangeable_H_flag": int(atom_b["exchangeable_H_flag"]),
                "touches_exchangeable_H": touches_exchangeable_h,
                "default_spectrum_model_disposition": (
                    "retained_as_calculated_for_audit_but_excluded_from_coherent_fine_structure"
                    if touches_exchangeable_h
                    else "eligible_for_spectrum_model_subject_to_threshold_and_decoupling"
                ),
                "distance_angstrom_as_reported_by_orca": row[
                    "distance_angstrom_as_reported_by_orca"
                ],
                "J_isotropic_Hz": row["J_isotropic_Hz"],
                "source_property_file": str(
                    property_path_from_output(
                        output_path_from_manifest(coupling_manifest[conformer_id])
                    )
                ),
            }
            coupling_raw_rows.append(combined)
            coupling_by_pair[(index_a, index_b)].append(combined)

    averaged_couplings: list[dict[str, Any]] = []
    expected_coupling_ids = set(parsed_couplings)
    for (index_a, index_b), rows in sorted(coupling_by_pair.items()):
        available_ids = {str(row["conformer_id"]) for row in rows}
        pair_population = sum(
            float(dft_records[conformer_id]["boltzmann_weight_all_valid_qrrho"])
            for conformer_id in available_ids
        )
        if available_ids != expected_coupling_ids:
            message = (
                f"J pair {index_a}-{index_b} is absent from "
                f"{sorted(expected_coupling_ids - available_ids)}; absent pairs were not treated as zero"
            )
            warnings.append(message)
            if not allow_partial:
                raise ProcessingError(message)
        values = [float(row["J_isotropic_Hz"]) for row in rows]
        full_weights = [
            float(dft_records[str(row["conformer_id"])]["boltzmann_weight_all_valid_qrrho"])
            for row in rows
        ]
        mean, spread = weighted_mean_sd(values, full_weights)
        distance_mean, distance_spread = weighted_mean_sd(
            [float(row["distance_angstrom_as_reported_by_orca"]) for row in rows],
            full_weights,
        )
        atom_a = atom_map[index_a]
        atom_b = atom_map[index_b]
        pair_class = scalar_pair_class(atom_a["element"], atom_b["element"])
        touches_exchangeable_h = int(
            int(atom_a["exchangeable_H_flag"]) == 1
            or int(atom_b["exchangeable_H_flag"]) == 1
        )
        averaged_couplings.append(
            {
                "notice": NOTICE,
                "figure_id": figure_id,
                "nucleus_A_index_0based": index_a,
                "nucleus_A_index_1based": index_a + 1,
                "nucleus_A_element": atom_a["element"],
                "nucleus_B_index_0based": index_b,
                "nucleus_B_index_1based": index_b + 1,
                "nucleus_B_element": atom_b["element"],
                "pair_label_1based": (
                    f"{atom_a['element']}{index_a + 1}-{atom_b['element']}{index_b + 1}"
                ),
                "scalar_coupling_pair_class": pair_class,
                "is_HH_pair": int(pair_class == "H-H"),
                "nucleus_A_exchangeable_H_flag": int(atom_a["exchangeable_H_flag"]),
                "nucleus_B_exchangeable_H_flag": int(atom_b["exchangeable_H_flag"]),
                "touches_exchangeable_H": touches_exchangeable_h,
                "default_spectrum_model_disposition": (
                    "retained_as_calculated_for_audit_but_excluded_from_coherent_fine_structure"
                    if touches_exchangeable_h
                    else "eligible_for_spectrum_model_subject_to_threshold_and_decoupling"
                ),
                "n_conformers_with_pair": len(rows),
                "pair_population_of_all_valid_dft": pair_population,
                "boltzmann_averaged_J_Hz": mean,
                "conformer_population_sd_J_Hz": spread,
                "minimum_conformer_J_Hz": min(values),
                "maximum_conformer_J_Hz": max(values),
                "boltzmann_averaged_distance_angstrom_as_reported_by_orca": distance_mean,
                "conformer_population_sd_distance_angstrom": distance_spread,
                "predictive_uncertainty_status": (
                    "not quantified: conformer SD is population spread only; no empirical J-coupling calibration error supplied"
                ),
                "flags": (
                    NOTICE
                    + (
                        ";EXCHANGEABLE_H_COUPLING_RETAINED_FOR_AUDIT_NOT_DEFAULT_FINE_STRUCTURE"
                        if touches_exchangeable_h
                        else ""
                    )
                ),
            }
        )

    # Generic scalar-coupling outputs are authoritative.  The legacy JHH keys and
    # CSVs remain as strict H-H filtered compatibility views, never aliases for
    # H-F, C-F, F-F, H-P, C-P, P-P or any other heteronuclear pair.
    jhh_raw_rows = [row for row in coupling_raw_rows if int(row["is_HH_pair"]) == 1]
    averaged_jhh = [row for row in averaged_couplings if int(row["is_HH_pair"]) == 1]

    valid_dft = [
        row for row in dft_records.values() if row.get("weighting_status") == "included"
    ]
    n_post_dft_duplicates_excluded = sum(
        row.get("weighting_status") == "excluded_post_dft_duplicate"
        for row in dft_records.values()
    )
    valid_dft.sort(key=lambda row: float(row["relative_qrrho_gibbs_kcal_mol"]))
    temperatures = [float(row["thermochemistry_temperature_K"]) for row in valid_dft]
    temperature = statistics.fmean(temperatures)
    first_dft = valid_dft[0]
    first_shield = next(iter(shielding_manifest.values()))
    first_coupling = next(iter(coupling_manifest.values()))
    compound = {
        "figure_id": figure_id,
        "compound_number": int(first_dft["compound_number"]),
        "display_name": first_dft["display_name"],
        "formula": first_dft["formula"],
        "isomeric_smiles": first_dft["isomeric_smiles"],
        "inchikey": first_dft["inchikey"],
        "charge": int(first_dft["charge"]),
        "multiplicity": int(first_dft["multiplicity"]),
    }
    conformer_rows: list[dict[str, Any]] = []
    for conformer_id, row in sorted(
        dft_records.items(),
        key=lambda pair: (
            0 if pair[1].get("weighting_status") == "included" else 1,
            float(pair[1].get("relative_qrrho_gibbs_kcal_mol", math.inf)),
            pair[0],
        ),
    ):
        conformer_rows.append(
            {
                "notice": NOTICE,
                "figure_id": figure_id,
                "conformer_id": conformer_id,
                "weighting_status": row.get("weighting_status"),
                "dft_audit_parse_status": row.get("dft_audit_parse_status"),
                "dft_unique_minimum": row.get("dft_unique_minimum"),
                "dft_duplicate_of": row.get("dft_duplicate_of"),
                "dft_duplicate_rmsd_angstrom": row.get(
                    "dft_duplicate_rmsd_angstrom"
                ),
                "dft_duplicate_delta_g_kcal_mol": row.get(
                    "dft_duplicate_delta_g_kcal_mol"
                ),
                "selected_for_nmr": row.get("selected_for_nmr"),
                "audit_boltzmann_weight_all_valid": row.get(
                    "audit_boltzmann_weight_all_valid"
                ),
                "audit_boltzmann_weight_selected_normalized": row.get(
                    "audit_boltzmann_weight_selected_normalized"
                ),
                "dft_audit_validation_status": row.get(
                    "dft_audit_validation_status"
                ),
                "qrrho_verified": row.get("qrrho_verified"),
                "thermochemistry_temperature_K": row.get("thermochemistry_temperature_K"),
                "gibbs_free_energy_hartree": row.get("gibbs_free_energy_hartree"),
                "relative_qrrho_gibbs_kcal_mol": row.get("relative_qrrho_gibbs_kcal_mol"),
                "boltzmann_weight_all_valid_qrrho": row.get(
                    "boltzmann_weight_all_valid_qrrho"
                ),
                "n_negative_frequencies": row.get("n_negative_frequencies"),
                "most_negative_frequency_cm-1": row.get(
                    "most_negative_frequency_cm-1"
                ),
                "in_shielding_manifest": int(conformer_id in shielding_manifest),
                "completed_shielding": int(conformer_id in parsed_shieldings),
                "in_coupling_manifest": int(conformer_id in coupling_manifest),
                "completed_coupling": int(conformer_id in parsed_couplings),
                "dft_output": row.get("dft_output"),
                "dft_output_sha256": row.get("dft_output_sha256"),
                "processing_error": row.get("processing_error"),
            }
        )

    assumptions = [
        NOTICE,
        "All reported values are calculations and cannot replace measured NMR characterization, identity or purity evidence.",
        "The paired DFT minimum/weight audit is authoritative for post-DFT deduplication and NMR selection. Rows with dft_unique_minimum=0 never enter the QRRHO partition function or property coverage.",
        "Conformer populations use ORCA 6.1.1 final QRRHO Gibbs free energies at the reported thermochemistry temperature.",
        "Property averages are renormalized over completed conformers listed in the corresponding property manifest; population coverage is reported explicitly.",
        "Raw chemical shifts use delta = mean same-level nucleus-specific reference shielding minus Boltzmann-averaged sample shielding, with the stated reference convention.",
        "The external H/C regression is from DOI 10.3390/ijms241914623 (PBE0/pcS-2, 23 natural products; 713 1H and 767 13C shifts). It is provisional here because the present basis is pcSseg-2 and the protocol/domain differ; its source R-squared is not a validated error bar for these compounds.",
        "No H/C regression is applied to 19F or 31P. If the corresponding CFCl3 or H3PO4 reference output is absent or invalid, the calculated shielding is retained and the raw shift is explicitly pending.",
        "The conventional 31P zero is an external 85% aqueous H3PO4 standard. An isolated neutral H3PO4/CPCM(DMSO) calculation does not reproduce that experimental reference environment; hydrogen bonding, ionization, concentration, water and susceptibility effects remain limitations.",
        "Conformer-population standard deviations quantify conformational spread only and are not experimental or predictive error bars.",
        "Exchangeable NH/OH shifts are especially sensitive to solvent, water, concentration, temperature and proton exchange; continuum DMSO does not model these effects explicitly.",
        "Hydrogens directly bound to N, O or S are structurally classified as exchangeable from the atom-map connectivity and cross-checked against the upstream flag. Their calculated shieldings and all scalar couplings are retained, but shift interpretation, linewidth and nominal integral are provisional solution-observation quantities.",
        "Scalar couplings touching a structurally exchangeable H are retained in the authoritative coupling tables for traceability but are excluded by default from the coherent fine-structure network used for display spectra.",
        "Scalar-coupling pairs absent because of the ORCA distance threshold are treated as missing, never silently as zero. Generic tables retain all parsed nucleus pairs; legacy JHH tables contain H-H pairs only.",
        "No linewidth, solvent peak, water peak, noise, baseline or experimental-looking spectrum was generated by this script.",
        "This processor reports configured H, C, F and P shieldings and scalar couplings. It does not infer 33S shifts without a nucleus-specific calculation and reference.",
    ]
    bundle = {
        "schema_version": SCHEMA_VERSION,
        "notice": NOTICE,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "generated_by": "scripts/process_nmr_predictions.py",
        "compound": compound,
        "methods": {
            "dft_opt_freq_input_keywords": method_line(
                Path(first_dft["job_directory"]) / first_dft["input_file"]
            ),
            "shielding_input_keywords": method_line(
                Path(first_shield["job_directory"]) / first_shield["input_file"]
            ),
            "coupling_input_keywords": method_line(
                Path(first_coupling["job_directory"]) / first_coupling["input_file"]
            ),
            "tms_reference_output": str(tms_output),
            "cfcl3_19f_reference_output": str(cfcl3_output),
            "h3po4_31p_reference_output": str(h3po4_output),
        },
        "weighting": {
            "energy_type": "ORCA final QRRHO Gibbs free energy",
            "temperature_K": temperature,
            "gas_constant_kcal_mol_K": R_KCAL_MOL_K,
            "hartree_to_kcal_mol": HARTREE_TO_KCAL_MOL,
            "n_valid_dft_conformers": len(valid_dft),
            "n_post_dft_duplicates_excluded": n_post_dft_duplicates_excluded,
            "shielding_selected_population_of_all_valid_dft": shielding_selected_population,
            "shielding_completed_population_of_all_valid_dft": shield_context[
                "completed_population_of_all_valid_dft"
            ],
            "shielding_completed_fraction_of_selected_population": shield_context[
                "completed_fraction_of_selected_population"
            ],
            "coupling_selected_population_of_all_valid_dft": coupling_selected_population,
            "coupling_completed_population_of_all_valid_dft": coupling_context[
                "completed_population_of_all_valid_dft"
            ],
            "coupling_completed_fraction_of_selected_population": coupling_context[
                "completed_fraction_of_selected_population"
            ],
        },
        "reference_shieldings": references,
        "provisional_external_regressions": PROVISIONAL_REGRESSION,
        "exchangeable_proton_policy": {
            "structural_definition": "explicit H directly bonded to N, O or S in atom_map connectivity",
            "assignment_validation": "upstream exchangeable_H_flag must agree with connectivity",
            "calculated_shieldings_retained": True,
            "calculated_scalar_couplings_retained": True,
            "default_coherent_fine_structure": "excluded for every coupling touching exchangeable H",
            "shift_status": "provisional because explicit hydrogen bonding and proton exchange are not modelled",
            "linewidth_status": "provisional display parameter, not a calculated relaxation/exchange rate",
            "integral_status": "provisional nominal area of one explicit H, not a quantitative experimental integral",
        },
        "atom_map": [atom_map[index] for index in sorted(atom_map)],
        "conformers": conformer_rows,
        "shieldings_by_conformer": shielding_raw_rows,
        "averaged_nuclei": averaged_nuclei,
        "scalar_couplings_by_conformer": coupling_raw_rows,
        "averaged_scalar_couplings": averaged_couplings,
        "JHH_by_conformer": jhh_raw_rows,
        "averaged_JHH": averaged_jhh,
        "assumptions_and_limitations": assumptions,
        "warnings": warnings,
        "source_files": source_files
        + [file_record(path) for path in source_manifest_paths.values()],
    }

    output_directory.mkdir(parents=True, exist_ok=True)
    prefix = figure_id
    atomic_write_csv(
        output_directory / f"{prefix}_conformer_QRRHO_weights.csv",
        conformer_rows,
        list(conformer_rows[0]),
    )
    atomic_write_csv(
        output_directory / f"{prefix}_shieldings_by_conformer.csv",
        shielding_raw_rows,
        list(shielding_raw_rows[0]),
    )
    atomic_write_csv(
        output_directory / f"{prefix}_calculated_shifts.csv",
        averaged_nuclei,
        list(averaged_nuclei[0]),
    )
    atomic_write_csv(
        output_directory / f"{prefix}_scalar_couplings_by_conformer.csv",
        coupling_raw_rows,
        list(coupling_raw_rows[0]),
    )
    atomic_write_csv(
        output_directory / f"{prefix}_calculated_scalar_couplings.csv",
        averaged_couplings,
        list(averaged_couplings[0]),
    )
    atomic_write_json(
        output_directory / f"{prefix}_scalar_couplings_by_conformer.json",
        {
            "schema_version": SCHEMA_VERSION,
            "notice": NOTICE,
            "figure_id": figure_id,
            "rows": coupling_raw_rows,
        },
    )
    atomic_write_json(
        output_directory / f"{prefix}_calculated_scalar_couplings.json",
        {
            "schema_version": SCHEMA_VERSION,
            "notice": NOTICE,
            "figure_id": figure_id,
            "rows": averaged_couplings,
        },
    )
    atomic_write_csv(
        output_directory / f"{prefix}_JHH_by_conformer.csv",
        jhh_raw_rows,
        list(coupling_raw_rows[0]),
    )
    atomic_write_csv(
        output_directory / f"{prefix}_calculated_JHH.csv",
        averaged_jhh,
        list(averaged_couplings[0]),
    )
    atomic_write_json(output_directory / f"{prefix}_prediction_bundle.json", bundle)
    atomic_write_text(output_directory / f"{prefix}_README.txt", build_human_report(bundle))
    return {
        "notice": NOTICE,
        "figure_id": figure_id,
        "display_name": compound["display_name"],
        "status": "processed_calculated_prediction_not_experimental",
        "n_valid_dft_conformers": len(valid_dft),
        "n_post_dft_duplicates_excluded": n_post_dft_duplicates_excluded,
        "n_averaged_nuclei": len(averaged_nuclei),
        "n_exchangeable_H": sum(
            int(row.get("exchangeable_H_flag", 0)) for row in averaged_nuclei
        ),
        "n_averaged_scalar_coupling_pairs": len(averaged_couplings),
        "n_averaged_JHH_pairs": len(averaged_jhh),
        "shielding_population_coverage": shield_context[
            "completed_population_of_all_valid_dft"
        ],
        "coupling_population_coverage": coupling_context[
            "completed_population_of_all_valid_dft"
        ],
        "n_warnings": len(warnings),
        "output_directory": str(output_directory),
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description=(
            "Parse ORCA 6.1.1 shielding/J outputs and calculate QRRHO-Boltzmann-averaged "
            "NMR predictions. Outputs are explicitly labelled calculated, not experimental."
        )
    )
    result.add_argument("--dft-manifest", type=Path, default=DEFAULT_DFT_MANIFEST)
    result.add_argument(
        "--dft-audit",
        type=Path,
        default=None,
        help=(
            "Optional DFT_minimum_and_weight_audit CSV. If omitted, a matching full/subset "
            "audit is inferred from the DFT and property manifest suffixes."
        ),
    )
    result.add_argument("--shielding-manifest", type=Path, default=DEFAULT_SHIELDING_MANIFEST)
    result.add_argument(
        "--coupling-manifest",
        type=Path,
        default=DEFAULT_COUPLING_MANIFEST,
        help=(
            "Scalar-coupling manifest (default: NMR_scalar_coupling_job_manifest.csv). "
            f"A legacy H-H-only manifest remains available only when passed explicitly, e.g. {LEGACY_HH_COUPLING_MANIFEST.name}."
        ),
    )
    result.add_argument("--tms-output", type=Path, default=DEFAULT_TMS_OUTPUT)
    result.add_argument(
        "--cfcl3-output",
        type=Path,
        default=DEFAULT_CFCL3_OUTPUT,
        help=(
            "Optional same-level CFCl3 shielding output for raw 19F referencing. "
            "If absent, 19F shieldings are retained and shifts are marked pending."
        ),
    )
    result.add_argument(
        "--h3po4-output",
        type=Path,
        default=DEFAULT_H3PO4_OUTPUT,
        help=(
            "Optional same-level H3PO4 shielding output for raw 31P referencing. "
            "If absent, 31P shieldings are retained and shifts are marked pending."
        ),
    )
    result.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    result.add_argument("--only", nargs="+", help="Optional Figure3 IDs, e.g. --only F3_001")
    result.add_argument(
        "--significant-imaginary-cm",
        type=float,
        default=DEFAULT_SIGNIFICANT_IMAGINARY_CM,
        help="Conformers with a frequency below this value are excluded (default: -20 cm^-1)",
    )
    result.add_argument(
        "--allow-partial",
        action="store_true",
        help="Write clearly flagged partial averages when selected property jobs are missing",
    )
    result.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace files previously written by this processor in selected compound directories",
    )
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    manifest_paths = {
        "dft_manifest": args.dft_manifest.expanduser().resolve(),
        "shielding_manifest": args.shielding_manifest.expanduser().resolve(),
        "coupling_manifest": args.coupling_manifest.expanduser().resolve(),
    }
    tms_output = args.tms_output.expanduser().resolve()
    cfcl3_output = args.cfcl3_output.expanduser().resolve()
    h3po4_output = args.h3po4_output.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()
    try:
        dft_audit_path = infer_dft_audit_path(
            args.dft_audit,
            manifest_paths["dft_manifest"],
            manifest_paths["shielding_manifest"],
            manifest_paths["coupling_manifest"],
        )
        manifest_paths["dft_audit"] = dft_audit_path
        dft_rows = read_manifest(
            manifest_paths["dft_manifest"],
            {
                "job_id",
                "figure_id",
                "compound_number",
                "display_name",
                "formula",
                "isomeric_smiles",
                "inchikey",
                "charge",
                "multiplicity",
                "job_directory",
                "input_file",
                "output_file",
            },
        )
        dft_minimum_audit_rows = read_manifest(
            dft_audit_path,
            {
                "job_id",
                "figure_id",
                "parse_status",
                "significant_imaginary",
                "boltzmann_weight_all_valid",
                "selected_for_nmr",
                "boltzmann_weight_selected_normalized",
            },
        )
        shielding_rows = read_manifest(
            manifest_paths["shielding_manifest"],
            {
                "job_id",
                "figure_id",
                "conformer_id",
                "job_directory",
                "input_file",
                "output_file",
            },
        )
        coupling_rows = read_manifest(
            manifest_paths["coupling_manifest"],
            {
                "job_id",
                "figure_id",
                "conformer_id",
                "job_directory",
                "input_file",
                "output_file",
            },
        )
        shielding_grouped = group_manifest_by_figure_and_conformer(shielding_rows)
        coupling_grouped = group_manifest_by_figure_and_conformer(coupling_rows)
        available = set(shielding_grouped) & set(coupling_grouped)
        if args.only:
            selected = set(args.only)
            missing = selected - available
            if missing:
                raise ProcessingError(
                    f"Selected IDs lack both shielding and coupling manifest rows: {sorted(missing)}"
                )
        else:
            selected = available
        if not selected:
            raise ProcessingError("No compounds selected")
        dft_grouped, dft_audit = prepare_dft_records(
            dft_rows,
            dft_minimum_audit_rows,
            selected,
            args.significant_imaginary_cm,
        )
        references = reference_summary(tms_output, cfcl3_output, h3po4_output)
        summaries: list[dict[str, Any]] = []
        for figure_id in sorted(selected):
            summaries.append(
                process_figure(
                    figure_id=figure_id,
                    dft_records=dft_grouped[figure_id],
                    shielding_manifest=shielding_grouped[figure_id],
                    coupling_manifest=coupling_grouped[figure_id],
                    references=references,
                    tms_output=tms_output,
                    cfcl3_output=cfcl3_output,
                    h3po4_output=h3po4_output,
                    output_root=output_root,
                    allow_partial=args.allow_partial,
                    overwrite=args.overwrite,
                    source_manifest_paths=manifest_paths,
                )
            )
    except (OSError, ValueError, ProcessingError) as exc:
        print(f"Processing failed: {exc}")
        return 2

    label = "ALL" if not args.only else "_".join(sorted(selected))
    summary_payload = {
        "schema_version": SCHEMA_VERSION,
        "notice": NOTICE,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "selected_figure_ids": sorted(selected),
        "results": summaries,
        "dft_parse_audit": [
            {
                "figure_id": row.get("figure_id"),
                "conformer_id": row.get("conformer_id"),
                "weighting_status": row.get("weighting_status"),
                "dft_unique_minimum": row.get("dft_unique_minimum"),
                "dft_duplicate_of": row.get("dft_duplicate_of"),
                "selected_for_nmr": row.get("selected_for_nmr"),
                "processing_error": row.get("processing_error"),
            }
            for row in dft_audit
        ],
    }
    output_root.mkdir(parents=True, exist_ok=True)
    atomic_write_json(output_root / f"processing_summary__{label}.json", summary_payload)
    atomic_write_csv(
        output_root / f"processing_summary__{label}.csv",
        summaries,
        list(summaries[0]),
    )
    print(NOTICE)
    print(f"Processed {len(summaries)} compound(s): {', '.join(sorted(selected))}")
    print(output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
