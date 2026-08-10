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

from nmr_method_validation import MethodValidationError, validate_nmr_method_evidence
from runner_provenance import RunnerProvenanceError, validate_runner_output


PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data"
DEFAULT_DFT_MANIFEST = DATA / "DFT_explicit_job_manifest.csv"
DEFAULT_WEIGHT_TABLE = DATA / "NMR_conformer_weights.csv"
DEFAULT_SHIELDING_MANIFEST = DATA / "NMR_shielding_job_manifest.csv"
DEFAULT_COUPLING_MANIFEST = DATA / "NMR_j_coupling_job_manifest.csv"
LEGACY_HH_COUPLING_MANIFEST = DATA / "NMR_HH_coupling_job_manifest.csv"
DEFAULT_REFERENCE_MANIFEST = DATA / "reference_job_manifest.csv"
DEFAULT_REFERENCE_REGISTRY = PROJECT / "config" / "reference_job_registry.json"
DEFAULT_NUCLEUS_CONDITIONS = PROJECT / "config" / "nucleus_conditions.csv"
DEFAULT_COMPOUND_REGISTRY = (
    PROJECT / "registry" / "Figure3_52_ORCA_reference_state_registry.csv"
)
DEFAULT_SHIELDING_METHOD = "production_revdsd_pcsseg3"
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
SCHEMA_VERSION = "1.5.0"
HARTREE_TO_KCAL_MOL = 627.5094740631
DEFAULT_SIGNIFICANT_IMAGINARY_CM = -20.0
WEIGHT_TOLERANCE = 5.0e-6

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
NUCLEUS_FROM_ELEMENT = {"H": "1H", "C": "13C", "F": "19F", "P": "31P"}


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


def read_csv_table(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
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
        raise ProcessingError(f"CSV table is empty: {path}")
    return rows


def select_shielding_method_rows(
    rows: Sequence[dict[str, str]], shielding_method: str
) -> tuple[list[dict[str, str]], str]:
    """Select exactly one shielding branch and prove that it has one fingerprint.

    The property manifest intentionally contains independent production and
    cross-check calculations for each conformer.  Averaging those rows together
    would silently combine different electronic-structure models.  Selection is
    therefore performed before grouping by conformer and is fail-closed.
    """
    method = str(shielding_method).strip()
    if not method:
        raise ProcessingError("--shielding-method must not be empty")
    required = {"property_branch", "method_fingerprint"}
    available_fields = set().union(*(row.keys() for row in rows)) if rows else set()
    missing = required - available_fields
    if missing:
        raise ProcessingError(
            "Shielding manifest cannot isolate electronic-structure methods; "
            f"missing columns: {sorted(missing)}"
        )
    selected = [
        row for row in rows if row.get("property_branch", "").strip() == method
    ]
    if not selected:
        available = sorted(
            {
                row.get("property_branch", "").strip()
                for row in rows
                if row.get("property_branch", "").strip()
            }
        )
        raise ProcessingError(
            f"No shielding rows match method {method!r}; available methods: {available}"
        )
    missing_fingerprint_jobs = [
        row.get("job_id", "<unknown>")
        for row in selected
        if not row.get("method_fingerprint", "").strip()
    ]
    if missing_fingerprint_jobs:
        raise ProcessingError(
            f"Selected shielding rows lack method_fingerprint: {missing_fingerprint_jobs}"
        )
    fingerprints = {
        row["method_fingerprint"].strip()
        for row in selected
        if row.get("method_fingerprint", "").strip()
    }
    if len(fingerprints) != 1:
        raise ProcessingError(
            f"Shielding method {method!r} has multiple fingerprints: {sorted(fingerprints)}"
        )
    return selected, next(iter(fingerprints))


def resolve_same_method_tms_output(
    reference_manifest: Path,
    shielding_method: str,
    shielding_method_fingerprint: str,
    explicit_tms_output: Path | None = None,
) -> tuple[Path, dict[str, str]]:
    """Resolve TMS from the reference row locked to the selected sample method.

    An explicitly relocated output is allowed, but it does not bypass the
    reference-manifest fingerprint check.
    """
    rows = read_manifest(
        reference_manifest,
        {
            "job_id",
            "stage",
            "method_id",
            "method_fingerprint",
            "job_directory",
            "output_file",
        },
    )
    candidates = [
        row
        for row in rows
        if row.get("stage", "").strip() == "reference_shielding"
        and row.get("method_id", "").strip() == shielding_method
    ]
    if len(candidates) != 1:
        raise ProcessingError(
            f"Expected exactly one TMS reference row for {shielding_method!r}, "
            f"found {len(candidates)} in {reference_manifest}"
        )
    reference_row = candidates[0]
    reference_fingerprint = reference_row.get("method_fingerprint", "").strip()
    if reference_fingerprint != shielding_method_fingerprint:
        raise ProcessingError(
            "Sample/TMS shielding method fingerprint mismatch: "
            f"sample={shielding_method_fingerprint!r}, TMS={reference_fingerprint!r}"
        )
    if explicit_tms_output is not None:
        output = explicit_tms_output.expanduser().resolve()
    else:
        directory = Path(reference_row["job_directory"]).expanduser()
        if not directory.is_absolute():
            directory = PROJECT / directory
        output = (directory / reference_row["output_file"]).resolve()
    return output, reference_row


def validate_approved_reference_registry(
    registry_path: Path,
    reference_manifest_row: dict[str, str],
    shielding_method: str,
    shielding_method_fingerprint: str,
    requested_tms_output: Path | None,
    reference_state_root: Path | None = None,
) -> tuple[Path, dict[str, Any]]:
    if not registry_path.is_file():
        raise FileNotFoundError(registry_path)
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    if registry.get("status") != "APPROVED_FOR_NMR_PROPERTIES":
        raise ProcessingError(f"Reference registry is not approved: {registry_path}")
    if not str(registry.get("reviewer", "")).strip() or not str(
        registry.get("approved_at_utc", "")
    ).strip():
        raise ProcessingError("Approved reference registry lacks reviewer/time provenance")
    entries = registry.get("reference_jobs")
    if not isinstance(entries, list):
        raise ProcessingError("Reference registry lacks reference_jobs")
    matched = [
        entry
        for entry in entries
        if isinstance(entry, dict)
        and entry.get("method_id") == shielding_method
        and entry.get("method_fingerprint") == shielding_method_fingerprint
        and entry.get("approval_status") == "APPROVED_FOR_NMR_PROPERTIES"
    ]
    by_nucleus = {entry.get("nucleus"): entry for entry in matched}
    if set(by_nucleus) != {"1H", "13C"} or len(matched) != 2:
        raise ProcessingError(
            f"Approved registry must contain exactly one 1H and one 13C TMS entry for {shielding_method}"
        )
    output_hashes = {str(entry.get("output_sha256", "")).upper() for entry in matched}
    output_paths = {str(entry.get("output_file", "")) for entry in matched}
    if len(output_hashes) != 1 or "" in output_hashes or len(output_paths) != 1:
        raise ProcessingError("1H/13C TMS registry entries do not bind one identical output")
    runner = validate_runner_output(
        reference_manifest_row,
        "reference_shielding",
        PROJECT,
        reference_state_root,
    )
    runner_output = Path(runner["output_path"])
    if runner["output_sha256"].upper() != next(iter(output_hashes)):
        raise ProcessingError("Reference registry output hash differs from runner provenance")
    registered_output = Path(next(iter(output_paths))).expanduser().resolve()
    if not registered_output.is_file() or sha256_file(registered_output).upper() != next(iter(output_hashes)):
        raise ProcessingError("Reference registry output path/hash is invalid")
    if registered_output != runner_output.resolve():
        raise ProcessingError("Reference registry output path differs from trusted runner output")
    reference_input = Path(str(matched[0].get("input_file", ""))).expanduser().resolve()
    method_validation = validate_nmr_method_evidence(
        reference_input,
        runner_output,
        shielding_method,
        shielding_method_fingerprint,
    )
    for entry in matched:
        if str(entry.get("reviewer", "")).strip() != str(registry["reviewer"]).strip():
            raise ProcessingError("Reference entry reviewer differs from registry approver")
        registered_state = Path(str(entry.get("runner_state_file", ""))).expanduser().resolve()
        if registered_state != Path(runner["runner_state"]).resolve():
            raise ProcessingError("Reference registry runner-state path differs from trusted state")
        if str(entry.get("runner_state_sha256", "")).upper() != str(
            runner["runner_state_sha256"]
        ).upper():
            raise ProcessingError("Reference registry runner-state hash is invalid")
        if str(entry.get("runner_kind", "")) != str(runner["runner_kind"]):
            raise ProcessingError("Reference registry runner kind differs from trusted state")
        if entry.get("vibrational_correction_policy") != (
            "DIAGNOSTIC_ONLY_NOT_APPLIED_TO_MAIN_DELTA"
        ):
            raise ProcessingError(
                "Reference vibrational/VPT2 correction is not approved as diagnostic-only"
            )
        forbidden_main_correction = any(
            entry.get(field) not in (None, "", False, 0, 0.0)
            for field in (
                "vpt2_correction_applied_to_main_delta",
                "vibrational_correction_applied_to_main_delta",
                "main_delta_correction_ppm",
            )
        )
        if forbidden_main_correction:
            raise ProcessingError(
                "VPT2/TMS vibrational corrections may not enter the main chemical-shift delta"
            )
        for path_field, hash_field in (
            ("input_file", "input_sha256"),
            ("coordinate_file", "coordinate_sha256"),
        ):
            artifact = Path(str(entry.get(path_field, ""))).expanduser().resolve()
            if not artifact.is_file() or sha256_file(artifact).upper() != str(
                entry.get(hash_field, "")
            ).upper():
                raise ProcessingError(f"Reference registry {path_field} hash is invalid")
    selected = runner_output
    if requested_tms_output is not None:
        explicit = requested_tms_output.expanduser().resolve()
        if not explicit.is_file() or sha256_file(explicit).upper() != runner["output_sha256"].upper():
            raise ProcessingError("Explicit --tms-output is not byte-identical to approved runner output")
        selected = explicit
    return selected, {
        "registry_file": str(registry_path.resolve()),
        "registry_sha256": sha256_file(registry_path),
        "registry_status": registry["status"],
        "registry_reviewer": registry["reviewer"],
        "approved_at_utc": registry["approved_at_utc"],
        "entries": matched,
        "runner_provenance": {
            **runner,
            "output_path": str(runner["output_path"]),
        },
        "method_validation": method_validation,
    }


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
    return [item for item in re.split(r"[;,\s]+", str(value).strip()) if item]


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
        if not attached_indices:
            declared_exchangeable = int(atom["exchangeable_H_flag"])
            atom["exchangeable_H_class"] = (
                "exchangeable_parent_unspecified"
                if declared_exchangeable
                else "not_exchangeable"
            )
            atom["exchangeability_assignment_method"] = (
                "upstream_nmr_atom_map_exchangeable_flag;adjacency_not_available"
            )
            atom["exchangeability_structural_confidence"] = "medium"
            if declared_exchangeable:
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


def _bool_from_map(value: Any, field: str, default: bool | None = None) -> bool:
    text = str(value or "").strip().casefold()
    if not text and default is not None:
        return default
    if text in {"1", "true", "yes", "y"}:
        return True
    if text in {"0", "false", "no", "n"}:
        return False
    raise ProcessingError(f"Invalid atom-map Boolean {field}={value!r}")


def _index_values(row: dict[str, str]) -> list[int]:
    zero_based_fields = (
        "attached_atom_indices_0based",
        "adjacent_atom_indices_0based",
        "adjacency_atom_indices_0based",
        "bonded_atom_indices_0based",
    )
    one_based_fields = (
        "attached_atom_indices_1based",
        "adjacent_atom_indices_1based",
        "adjacency_atom_indices_1based",
        "bonded_atom_indices_1based",
    )
    for field in zero_based_fields:
        if str(row.get(field, "")).strip():
            return [int(value) for value in re.findall(r"\d+", row[field])]
    for field in one_based_fields:
        if str(row.get(field, "")).strip():
            return [int(value) - 1 for value in re.findall(r"\d+", row[field])]
    return []


def parse_atom_map(path: Path) -> dict[int, dict[str, Any]]:
    """Parse legacy ORCA-index maps and the portable nmr_atom_map_v2 schema."""
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ProcessingError(f"Empty atom map: {path}")
    result: dict[int, dict[str, Any]] = {}
    for row in rows:
        if str(row.get("orca_index_0based", "")).strip():
            index = int(row["orca_index_0based"])
        elif str(row.get("orca_index_1based", "")).strip():
            index = int(row["orca_index_1based"]) - 1
        elif str(row.get("atom_index_1based", "")).strip():
            index = int(row["atom_index_1based"]) - 1
        else:
            raise ProcessingError(f"Atom map has no supported index column: {path}")
        if index in result:
            raise ProcessingError(f"Duplicate atom index {index} in {path}")
        one_based = int(row.get("orca_index_1based") or row.get("atom_index_1based") or index + 1)
        if one_based != index + 1:
            raise ProcessingError(f"Inconsistent one-based atom index in {path}: {row}")
        element = str(row.get("element", "")).strip()
        if not element:
            raise ProcessingError(f"Blank element at atom {one_based} in {path}")
        environment = _bool_from_map(row.get("environment_atom"), "environment_atom", False)
        if str(row.get("include_in_analyte_nmr", "")).strip():
            include_in_analyte_nmr = _bool_from_map(
                row.get("include_in_analyte_nmr"), "include_in_analyte_nmr"
            )
        else:
            include_in_analyte_nmr = not environment
        nmr_scope = str(row.get("nmr_scope", "")).strip()
        if not nmr_scope:
            if include_in_analyte_nmr:
                nmr_scope = "ANALYTE"
            elif environment:
                nmr_scope = "ENVIRONMENT_EXCLUDED_FROM_ANALYTE_NMR"
            else:
                nmr_scope = "ANALYTE_EXCLUDED_FROM_ORDINARY_NMR"
        declared_exchangeable = _bool_from_map(
            row.get("exchangeable_H_flag", row.get("exchangeable")),
            "exchangeable",
            False,
        )
        attached = _index_values(row)
        result[index] = {
            "nucleus_index_0based": index,
            "nucleus_index_1based": one_based,
            "element": element,
            "formal_charge": int(row.get("formal_charge") or 0),
            "canonical_symmetry_rank": int(
                row.get("canonical_symmetry_rank")
                or row.get("source_atom_index_1based")
                or one_based
            ),
            "component": str(row.get("component", "analyte")).strip() or "analyte",
            "environment_atom": int(environment),
            "include_in_analyte_nmr": int(include_in_analyte_nmr),
            "nmr_scope": nmr_scope,
            "attached_atom_indices_0based": ";".join(str(value) for value in attached),
            "attached_elements": str(row.get("attached_elements", "")).strip(),
            "exchangeable_H_flag": int(declared_exchangeable),
        }
    if sorted(result) != list(range(len(result))):
        raise ProcessingError(f"Atom map indices are not contiguous in {path}")
    for index, atom in result.items():
        attached = [
            int(value)
            for value in _split_semicolon_values(atom["attached_atom_indices_0based"])
        ]
        if attached and not atom["attached_elements"]:
            try:
                atom["attached_elements"] = ";".join(result[value]["element"] for value in attached)
            except KeyError as exc:
                raise ProcessingError(
                    f"Atom {index} has out-of-range adjacency in {path}: {attached}"
                ) from exc
    return annotate_exchangeable_hydrogens(result, path)


def load_scenario_atom_map(
    figure_id: str,
    scenario_id: str,
    *property_manifests: dict[str, dict[str, str]],
) -> tuple[Path, dict[int, dict[str, Any]], str, list[Path]]:
    """Load job-local maps and require byte-identical maps within one scenario."""
    paths: list[Path] = []
    for manifest in property_manifests:
        for row in manifest.values():
            path = Path(row["job_directory"]).expanduser().resolve() / "nmr_atom_map_v2.csv"
            if path not in paths:
                paths.append(path)
    if not paths:
        raise ProcessingError(f"No property atom maps for {figure_id}/{scenario_id}")
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            f"Missing job-local nmr_atom_map_v2.csv for {figure_id}/{scenario_id}: {missing[0]}"
        )
    hashes = {sha256_file(path) for path in paths}
    if len(hashes) != 1:
        raise ProcessingError(
            f"Atom-map hash mismatch inside scenario {figure_id}/{scenario_id}: "
            f"{sorted((str(path), sha256_file(path)) for path in paths)}"
        )
    actual_hash = next(iter(hashes)).upper()
    for manifest in property_manifests:
        for row in manifest.values():
            expected = str(row.get("nmr_atom_map_sha256", "")).strip().upper()
            if expected and expected != actual_hash:
                raise ProcessingError(
                    f"Manifest atom-map hash mismatch for {row.get('job_id')}: "
                    f"manifest={expected}, actual={actual_hash}"
                )
    path = paths[0]
    return path, parse_atom_map(path), next(iter(hashes)), paths


def validate_scenario_isotope_maps(
    figure_id: str,
    scenario_id: str,
    atom_map: dict[int, dict[str, Any]],
    *property_manifests: dict[str, dict[str, str]],
) -> dict[str, Any]:
    paths: list[Path] = []
    for manifest in property_manifests:
        for row in manifest.values():
            directory = Path(row["job_directory"]).expanduser().resolve()
            path = directory / "source_isotope_map.csv"
            if path not in paths:
                paths.append(path)
            expected = str(row.get("source_isotope_map_sha256", "")).strip().upper()
            if not expected:
                raise ProcessingError(
                    f"Property manifest lacks source_isotope_map_sha256: {row.get('job_id')}"
                )
            if not path.is_file() or sha256_file(path).upper() != expected:
                raise ProcessingError(f"Property isotope-map hash mismatch: {path}")
    hashes = {sha256_file(path).upper() for path in paths}
    if len(hashes) != 1:
        raise ProcessingError(
            f"Isotope-map hash mismatch inside scenario {figure_id}/{scenario_id}"
        )
    with paths[0].open("r", encoding="utf-8-sig", newline="") as handle:
        isotope_rows = list(csv.DictReader(handle))
    deuterium: dict[int, float] = {}
    for row in isotope_rows:
        try:
            index = int(row["atom_index_1based"]) - 1
            mass = as_float(row["mass_u"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ProcessingError(f"Malformed isotope-map row: {row}") from exc
        if index not in atom_map or atom_map[index]["element"].upper() != "H":
            raise ProcessingError(f"Isotope map points to a non-H/absent atom: {row}")
        if index in deuterium or abs(mass - 2.01410177812) > 5.0e-11:
            raise ProcessingError(f"Invalid/duplicate exact deuterium mass: {row}")
        deuterium[index] = mass
    expected_dmso_d = {
        index
        for index, atom in atom_map.items()
        if atom["element"].upper() == "H"
        and "dmso" in str(atom.get("component", "")).strip().casefold()
        and int(atom.get("environment_atom", 0)) == 1
    }
    if set(deuterium) != expected_dmso_d:
        raise ProcessingError(
            f"DMSO-d6 isotope mapping was lost for {figure_id}/{scenario_id}: "
            f"expected={sorted(expected_dmso_d)}, mapped={sorted(deuterium)}"
        )
    return {
        "source_isotope_map_sha256": next(iter(hashes)),
        "source_isotope_map_copies_checked": len(paths),
        "deuterium_atom_indices_1based": [index + 1 for index in sorted(deuterium)],
        "deuterium_mass_u": 2.01410177812,
        "status": "PASS_EXPLICIT_DMSO_D6_MASS_MAPPING_PRESERVED",
        "paths": paths,
    }


def property_coordinate_path(row: dict[str, str]) -> Path:
    directory = Path(row["job_directory"]).expanduser().resolve()
    candidates: list[Path] = []
    configured = str(row.get("coordinate_file", "")).strip()
    if configured:
        candidates.append(directory / configured)
    candidates.extend((directory / "optimized.xyz", directory / "conformer.xyz"))
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(
        f"No property coordinate file found in {directory}; tried "
        + ", ".join(path.name for path in candidates)
    )


def validate_property_provenance(row: dict[str, str]) -> dict[str, Any]:
    path = Path(row["job_directory"]).expanduser().resolve() / "provenance.json"
    if not path.is_file():
        raise ProcessingError(f"Missing property provenance: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    for field in ("job_id", "figure_id", "scenario_id", "method_fingerprint"):
        if str(payload.get(field, "")) != str(row.get(field, "")):
            raise ProcessingError(f"Property provenance {field} mismatch: {path}")
    hash_bindings = (
        ("nmr_atom_map_sha256", "nmr_atom_map_sha256"),
        ("source_isotope_map_sha256", "source_isotope_map_sha256"),
        ("source_opt_freq_state_sha256", "source_opt_freq_state_sha256"),
        ("approved_minimum_gate_sha256", "approved_minimum_gate_sha256"),
        ("conformer_xyz_sha256", "optimized_coordinate_sha256"),
    )
    for manifest_field, provenance_field in hash_bindings:
        expected = str(row.get(manifest_field, "")).strip().upper()
        actual = str(payload.get(provenance_field, "")).strip().upper()
        if not expected or not re.fullmatch(r"[0-9A-F]{64}", expected):
            raise ProcessingError(
                f"Property manifest lacks valid {manifest_field}: {row.get('job_id')}"
            )
        if actual != expected:
            raise ProcessingError(
                f"Property provenance {provenance_field} hash mismatch: {path}"
            )
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "payload": payload,
    }


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


def group_manifest_by_state_and_conformer(
    rows: Sequence[dict[str, str]],
) -> dict[tuple[str, str], dict[str, dict[str, str]]]:
    """Partition property jobs by immutable (compound, solvent scenario) state."""
    grouped: dict[tuple[str, str], dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        figure_id = row.get("figure_id", "").strip()
        scenario_id = row.get("scenario_id", "").strip()
        conformer_id = row.get("conformer_id", "").strip()
        if not figure_id or not scenario_id or not conformer_id:
            raise ProcessingError(
                "Property manifest rows require non-empty figure_id, scenario_id and conformer_id"
            )
        key = (figure_id, scenario_id)
        if conformer_id in grouped[key]:
            raise ProcessingError(
                f"Duplicate {figure_id}/{scenario_id}/{conformer_id} in property manifest"
            )
        grouped[key][conformer_id] = row
    return dict(grouped)


def select_state_keys(
    available: set[tuple[str, str]],
    only_figures: Sequence[str] | None = None,
    only_scenarios: Sequence[str] | None = None,
) -> set[tuple[str, str]]:
    selected = set(available)
    if only_figures:
        requested_figures = {str(value).strip() for value in only_figures}
        present_figures = {figure_id for figure_id, _scenario_id in available}
        missing = requested_figures - present_figures
        if missing:
            raise ProcessingError(
                f"Selected IDs lack both shielding and coupling manifest rows: {sorted(missing)}"
            )
        selected = {key for key in selected if key[0] in requested_figures}
    if only_scenarios:
        requested_scenarios = {str(value).strip() for value in only_scenarios}
        present_scenarios = {scenario_id for _figure_id, scenario_id in selected}
        missing = requested_scenarios - present_scenarios
        if missing:
            raise ProcessingError(
                f"Selected scenario IDs lack paired property rows: {sorted(missing)}"
            )
        selected = {key for key in selected if key[1] in requested_scenarios}
    if not selected:
        raise ProcessingError("No compound/scenario states selected")
    return selected


def load_locked_nucleus_conditions(path: Path) -> dict[tuple[str, str], dict[str, Any]]:
    rows = read_csv_table(path, {"figure_id", "nucleus", "temperature_K", "condition_status"})
    result: dict[tuple[str, str], dict[str, Any]] = {}
    for row in rows:
        key = (row["figure_id"].strip(), row["nucleus"].strip())
        if key in result:
            raise ProcessingError(f"Duplicate nucleus condition: {key}")
        temperature = as_float(row["temperature_K"])
        if temperature <= 0.0 or not row["condition_status"].strip():
            raise ProcessingError(f"Invalid locked nucleus condition: {row}")
        result[key] = {**row, "temperature_K": temperature}
    return result


def load_nucleus_weight_sets(
    rows: Sequence[dict[str, str]],
    selected_states: set[tuple[str, str]],
    locked_conditions: dict[tuple[str, str], dict[str, Any]],
) -> dict[tuple[str, str], dict[str, dict[str, Any]]]:
    """Validate authoritative scenario/nucleus/temperature-specific conformer weights."""
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        state = (row["figure_id"].strip(), row["scenario_id"].strip())
        if state in selected_states:
            grouped[(state[0], state[1], row["nucleus"].strip())].append(row)
    result: dict[tuple[str, str], dict[str, dict[str, Any]]] = defaultdict(dict)
    for (figure_id, scenario_id, nucleus), values in sorted(grouped.items()):
        condition = locked_conditions.get((figure_id, nucleus))
        if condition is None:
            raise ProcessingError(f"No locked acquisition condition for {figure_id}/{nucleus}")
        temperatures = {as_float(row["temperature_K"]) for row in values}
        if len(temperatures) != 1:
            raise ProcessingError(
                f"Weight table mixes temperatures for {figure_id}/{scenario_id}/{nucleus}: {temperatures}"
            )
        temperature = next(iter(temperatures))
        if abs(temperature - float(condition["temperature_K"])) > 1.0e-4:
            raise ProcessingError(
                f"Weight temperature {temperature:.6f} K does not match locked "
                f"{figure_id}/{nucleus} temperature {condition['temperature_K']:.6f} K"
            )
        by_conformer: dict[str, dict[str, Any]] = {}
        for row in values:
            conformer_id = row["conformer_id"].strip()
            if not conformer_id or conformer_id in by_conformer:
                raise ProcessingError(
                    f"Blank/duplicate conformer weight for {figure_id}/{scenario_id}/{nucleus}"
                )
            weight = as_float(row["weight"])
            energy = as_float(row["qrrho_gibbs_free_energy_hartree"])
            if weight < 0.0:
                raise ProcessingError(f"Negative conformer weight: {row}")
            by_conformer[conformer_id] = {
                **row,
                "weight": weight,
                "qrrho_gibbs_free_energy_hartree": energy,
            }
        total = math.fsum(record["weight"] for record in by_conformer.values())
        if abs(total - 1.0) > WEIGHT_TOLERANCE:
            raise ProcessingError(
                f"Weights sum to {total:.12g}, not 1 for {figure_id}/{scenario_id}/{nucleus}"
            )
        minimum = min(
            record["qrrho_gibbs_free_energy_hartree"] for record in by_conformer.values()
        )
        for record in by_conformer.values():
            record["relative_qrrho_gibbs_kcal_mol"] = (
                record["qrrho_gibbs_free_energy_hartree"] - minimum
            ) * HARTREE_TO_KCAL_MOL
        gate_hashes = {
            str(row.get("approved_gate_sha256", "")).strip().upper()
            for row in values
        }
        if len(gate_hashes) != 1 or not re.fullmatch(
            r"[0-9A-F]{64}", next(iter(gate_hashes), "")
        ):
            raise ProcessingError(
                f"Weight rows require one valid approved-gate SHA-256 for "
                f"{figure_id}/{scenario_id}/{nucleus}: {sorted(gate_hashes)}"
            )
        result[(figure_id, scenario_id)][nucleus] = {
            "figure_id": figure_id,
            "scenario_id": scenario_id,
            "nucleus": nucleus,
            "temperature_K": temperature,
            "condition_status": condition["condition_status"],
            "condition_source_note": condition.get("source_note", ""),
            "conformers": by_conformer,
            "approved_gate_sha256": next(iter(gate_hashes)),
        }
    for state in selected_states:
        missing = {"1H", "13C"} - set(result.get(state, {}))
        if missing:
            raise ProcessingError(
                f"NMR_conformer_weights.csv lacks required nucleus sets for {state}: {sorted(missing)}"
            )
    return dict(result)


def validate_weight_gate_binding(
    figure_id: str,
    scenario_id: str,
    weight_sets: dict[str, dict[str, Any]],
    *property_manifests: dict[str, dict[str, str]],
) -> dict[str, Any]:
    """Bind every property row to the exact approved gate that emitted the weights."""
    weight_hashes = {
        str(weight_set.get("approved_gate_sha256", "")).strip().upper()
        for weight_set in weight_sets.values()
    }
    if len(weight_hashes) != 1 or not re.fullmatch(
        r"[0-9A-F]{64}", next(iter(weight_hashes), "")
    ):
        raise ProcessingError(
            f"Nucleus weight sets do not share one approved gate for "
            f"{figure_id}/{scenario_id}: {sorted(weight_hashes)}"
        )
    expected = next(iter(weight_hashes))
    checked_jobs: list[str] = []
    gate_paths: set[Path] = set()
    for manifest in property_manifests:
        for row in manifest.values():
            actual = str(row.get("approved_minimum_gate_sha256", "")).strip().upper()
            if actual != expected:
                raise ProcessingError(
                    f"Property/weight approved-gate mismatch for {row.get('job_id')}: "
                    f"property={actual or '<missing>'}, weights={expected}"
                )
            provenance_path = (
                Path(row["job_directory"]).expanduser().resolve() / "provenance.json"
            )
            if not provenance_path.is_file():
                raise ProcessingError(f"Missing property provenance: {provenance_path}")
            provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
            gate_text = str(provenance.get("approved_minimum_gate_file", "")).strip()
            if not gate_text:
                raise ProcessingError(
                    f"Property provenance lacks approved gate path: {row.get('job_id')}"
                )
            gate_path = Path(gate_text).expanduser().resolve()
            if not gate_path.is_file() or sha256_file(gate_path).upper() != expected:
                raise ProcessingError(
                    f"Approved minimum gate path/hash is invalid for {row.get('job_id')}: {gate_path}"
                )
            gate_paths.add(gate_path)
            checked_jobs.append(str(row.get("job_id", "")))
    if not checked_jobs:
        raise ProcessingError(f"No property jobs available for gate binding: {figure_id}/{scenario_id}")
    if len(gate_paths) != 1:
        raise ProcessingError(
            f"Property jobs do not share one approved gate for {figure_id}/{scenario_id}: "
            f"{sorted(str(path) for path in gate_paths)}"
        )
    gate_path = next(iter(gate_paths))
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    if (
        gate.get("approval_status") != "APPROVED_FOR_NMR"
        or gate.get("figure_id") != figure_id
        or not str(gate.get("reviewer", "")).strip()
        or not str(gate.get("reviewed_at", "")).strip()
    ):
        raise ProcessingError(f"Minimum/Boltzmann gate is not fully approved: {gate_path}")
    gate_weight_sets = gate.get("boltzmann_weight_sets")
    if not isinstance(gate_weight_sets, list):
        raise ProcessingError(f"Approved gate lacks Boltzmann weight sets: {gate_path}")
    upstream_states = gate.get("upstream_state_hashes")
    if not isinstance(upstream_states, list):
        raise ProcessingError(f"Approved gate lacks upstream runner-state hashes: {gate_path}")
    upstream_by_job = {
        str(item.get("job_id", "")): str(item.get("sha256", "")).strip().upper()
        for item in upstream_states
        if isinstance(item, dict)
    }
    for manifest in property_manifests:
        for conformer_id, row in manifest.items():
            source_state_hash = str(
                row.get("source_opt_freq_state_sha256", "")
            ).strip().upper()
            if (
                not re.fullmatch(r"[0-9A-F]{64}", source_state_hash)
                or upstream_by_job.get(conformer_id) != source_state_hash
            ):
                raise ProcessingError(
                    f"Property/gate upstream DFT-state mismatch for {conformer_id}"
                )
    for nucleus, weight_set in weight_sets.items():
        candidates = [
            item
            for item in gate_weight_sets
            if isinstance(item, dict)
            and item.get("condition_scenario_id") == scenario_id
            and item.get("nucleus") == nucleus
        ]
        if len(candidates) != 1 or candidates[0].get("approval_status") != "APPROVED_FOR_NMR":
            raise ProcessingError(
                f"Gate lacks one approved {scenario_id}/{nucleus} weight set: {gate_path}"
            )
        gate_set = candidates[0]
        if abs(float(gate_set.get("temperature_K")) - float(weight_set["temperature_K"])) > 1.0e-4:
            raise ProcessingError(f"Gate/table temperature mismatch for {scenario_id}/{nucleus}")
        supplied = gate_set.get("weights")
        if not isinstance(supplied, list):
            raise ProcessingError(f"Malformed gate weights for {scenario_id}/{nucleus}")
        gate_rows = {
            str(item.get("conformer_id", "")): item
            for item in supplied
            if isinstance(item, dict)
        }
        if set(gate_rows) != set(weight_set["conformers"]):
            raise ProcessingError(f"Gate/table conformer mismatch for {scenario_id}/{nucleus}")
        for conformer_id, table_row in weight_set["conformers"].items():
            gate_row = gate_rows[conformer_id]
            if (
                abs(float(gate_row.get("weight")) - float(table_row["weight"]))
                > WEIGHT_TOLERANCE
                or abs(
                    float(gate_row.get("qrrho_gibbs_free_energy_hartree"))
                    - float(table_row["qrrho_gibbs_free_energy_hartree"])
                )
                > 1.0e-10
            ):
                raise ProcessingError(
                    f"Gate/table value mismatch for {scenario_id}/{nucleus}/{conformer_id}"
                )
    return {
        "approved_minimum_gate_sha256": expected,
        "approved_minimum_gate_file": str(gate_path),
        "reviewer": gate["reviewer"],
        "reviewed_at": gate["reviewed_at"],
        "property_rows_checked": len(checked_jobs),
        "status": "PASS_PROPERTY_AND_WEIGHT_TABLE_SHARE_APPROVED_GATE",
    }


def validate_property_dft_state_binding(
    dft_records: dict[str, dict[str, Any]],
    *property_manifests: dict[str, dict[str, str]],
) -> dict[str, Any]:
    checked = 0
    for manifest in property_manifests:
        for conformer_id, row in manifest.items():
            if conformer_id not in dft_records:
                raise ProcessingError(f"Property conformer lacks trusted DFT state: {conformer_id}")
            expected = str(
                dft_records[conformer_id]
                .get("runner_provenance", {})
                .get("runner_state_sha256", "")
            ).strip().upper()
            supplied = str(row.get("source_opt_freq_state_sha256", "")).strip().upper()
            if not expected or supplied != expected:
                raise ProcessingError(
                    f"Property/current DFT runner-state mismatch for {conformer_id}: "
                    f"property={supplied or '<missing>'}, DFT={expected or '<missing>'}"
                )
            checked += 1
    return {
        "property_rows_checked": checked,
        "status": "PASS_PROPERTY_JOBS_BIND_CURRENT_TRUSTED_DFT_RUNNER_STATES",
    }


def weight_context(
    weight_set: dict[str, Any],
    selected_ids: set[str],
    completed_ids: set[str],
) -> dict[str, Any]:
    records = weight_set["conformers"]
    if set(records) != selected_ids:
        raise ProcessingError(
            f"Property conformers differ from {weight_set['nucleus']} weight set: "
            f"weights_only={sorted(set(records) - selected_ids)}, "
            f"properties_only={sorted(selected_ids - set(records))}"
        )
    missing_completed = completed_ids - set(records)
    if missing_completed:
        raise ProcessingError(f"Completed conformers absent from weight table: {sorted(missing_completed)}")
    completed_population = math.fsum(records[key]["weight"] for key in completed_ids)
    if completed_population <= 0.0:
        raise ProcessingError("Completed property jobs have zero authoritative population")
    return {
        "nucleus": weight_set["nucleus"],
        "temperature_K": weight_set["temperature_K"],
        "selected_population": 1.0,
        "completed_population": completed_population,
        "completed_fraction": completed_population,
        "normalized_completed_weights": {
            key: records[key]["weight"] / completed_population for key in completed_ids
        },
    }


def prepare_dft_records_from_weight_sets(
    rows: Sequence[dict[str, str]],
    selected_states: set[tuple[str, str]],
    weight_sets: dict[tuple[str, str], dict[str, dict[str, Any]]],
    significant_imaginary_cm: float,
    dft_state_root: Path | None = None,
) -> tuple[dict[tuple[str, str], dict[str, dict[str, Any]]], list[dict[str, Any]]]:
    """Bind selected conformers to trusted opt/freq outputs without reweighting them."""
    rows_by_state_job: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in rows:
        key = (
            row.get("figure_id", "").strip(),
            row.get("scenario_id", "").strip(),
            row.get("job_id", "").strip(),
        )
        if key in rows_by_state_job:
            raise ProcessingError(f"Duplicate DFT state/job row: {key}")
        rows_by_state_job[key] = row
    grouped: dict[tuple[str, str], dict[str, dict[str, Any]]] = {}
    audit: list[dict[str, Any]] = []
    for state in sorted(selected_states):
        nucleus_sets = weight_sets[state]
        conformer_sets = {
            nucleus: set(weight_set["conformers"])
            for nucleus, weight_set in nucleus_sets.items()
        }
        first_set = next(iter(conformer_sets.values()))
        for nucleus, conformers in conformer_sets.items():
            if conformers != first_set:
                raise ProcessingError(
                    f"Nucleus-specific conformer membership differs for {state}: "
                    f"{nucleus}={sorted(conformers)}, expected={sorted(first_set)}"
                )
        records: dict[str, dict[str, Any]] = {}
        for conformer_id in sorted(first_set):
            manifest_row = rows_by_state_job.get((state[0], state[1], conformer_id))
            if manifest_row is None:
                raise ProcessingError(f"Weight conformer absent from DFT manifest: {state}/{conformer_id}")
            try:
                runner = validate_runner_output(
                    manifest_row, "optfreq", PROJECT, dft_state_root
                )
                parsed = parse_dft_qrrho(
                    Path(runner["output_path"]), significant_imaginary_cm
                )
                if parsed["significant_imaginary"]:
                    raise ProcessingError(
                        f"Selected conformer has a significant imaginary mode: {state}/{conformer_id}"
                    )
                record = {
                    **manifest_row,
                    **parsed,
                    "conformer_id": conformer_id,
                    "weighting_status": "included_from_authoritative_nucleus_weight_table",
                    "dft_output": str(runner["output_path"]),
                    "dft_output_sha256": runner["output_sha256"],
                    "runner_provenance": {
                        **runner,
                        "output_path": str(runner["output_path"]),
                    },
                    "weights_by_nucleus": {
                        nucleus: weight_set["conformers"][conformer_id]["weight"]
                        for nucleus, weight_set in nucleus_sets.items()
                    },
                    "qrrho_gibbs_by_nucleus_hartree": {
                        nucleus: weight_set["conformers"][conformer_id][
                            "qrrho_gibbs_free_energy_hartree"
                        ]
                        for nucleus, weight_set in nucleus_sets.items()
                    },
                    "relative_qrrho_gibbs_by_nucleus_kcal_mol": {
                        nucleus: weight_set["conformers"][conformer_id][
                            "relative_qrrho_gibbs_kcal_mol"
                        ]
                        for nucleus, weight_set in nucleus_sets.items()
                    },
                }
            except (OSError, ValueError, ProcessingError, RunnerProvenanceError) as exc:
                raise ProcessingError(
                    f"DFT runner/provenance validation failed for {state}/{conformer_id}: {exc}"
                ) from exc
            records[conformer_id] = record
            audit.append(record)
        grouped[state] = records
    return grouped, audit


def prepare_dft_records(
    rows: Sequence[dict[str, str]],
    dft_audit_rows: Sequence[dict[str, str]],
    selected_states: set[tuple[str, str]],
    significant_imaginary_cm: float,
) -> tuple[
    dict[tuple[str, str], dict[str, dict[str, Any]]], list[dict[str, Any]]
]:
    selected_job_ids = {
        row["job_id"]
        for row in rows
        if (row.get("figure_id", "").strip(), row.get("scenario_id", "").strip())
        in selected_states
    }
    audit_by_job: dict[str, dict[str, str]] = {}
    for row in dft_audit_rows:
        job_id = row.get("job_id", "").strip()
        if job_id not in selected_job_ids:
            continue
        if not job_id:
            raise ProcessingError("DFT audit contains a selected-figure row without job_id")
        if job_id in audit_by_job:
            raise ProcessingError(f"Duplicate job_id in DFT audit: {job_id}")
        audit_by_job[job_id] = row

    grouped_raw: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (row.get("figure_id", "").strip(), row.get("scenario_id", "").strip())
        if key in selected_states:
            grouped_raw[key].append(row)
    missing = selected_states - set(grouped_raw)
    if missing:
        raise ProcessingError(
            f"Selected compound/scenario states absent from DFT manifest: {sorted(missing)}"
        )

    grouped: dict[tuple[str, str], dict[str, dict[str, Any]]] = {}
    audit: list[dict[str, Any]] = []
    for state_key, figure_rows in sorted(grouped_raw.items()):
        figure_id, scenario_id = state_key
        state_label = f"{figure_id}/{scenario_id}"
        if len(figure_rows) > 1 and not all(
            "dft_unique_minimum" in audit_by_job.get(row["job_id"], {})
            for row in figure_rows
        ):
            raise ProcessingError(
                f"{state_label}: multiple DFT rows require an audit with dft_unique_minimum; "
                "regenerate the post-DFT minimum/weight audit before processing"
            )
        parsed_valid: list[dict[str, Any]] = []
        all_records: dict[str, dict[str, Any]] = {}
        for row in figure_rows:
            conformer_id = row["job_id"]
            audit_row = audit_by_job.get(conformer_id)
            if audit_row is None:
                raise ProcessingError(
                    f"{state_label}/{conformer_id} is absent from the inferred DFT minimum/weight audit"
                )
            if audit_row.get("figure_id") != figure_id:
                raise ProcessingError(
                    f"DFT audit figure mismatch for {conformer_id}: {audit_row.get('figure_id')} vs {figure_id}"
                )
            audit_scenario = audit_row.get("scenario_id", "").strip()
            if audit_scenario and audit_scenario != scenario_id:
                raise ProcessingError(
                    f"DFT audit scenario mismatch for {conformer_id}: "
                    f"{audit_scenario} vs {scenario_id}"
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
                                f"{state_label}/{conformer_id}: DFT audit significant-imaginary flag "
                                "does not match the parsed ORCA output"
                            )
                    audit_weight_text = audit_row.get(
                        "boltzmann_weight_all_valid", ""
                    ).strip()
                    if not audit_weight_text:
                        raise ProcessingError(
                            f"{state_label}/{conformer_id}: unique valid minimum lacks "
                            "boltzmann_weight_all_valid in the DFT audit"
                        )
                    audit_weight = as_float(audit_weight_text)
                    if audit_weight < 0.0:
                        raise ProcessingError(
                            f"{state_label}/{conformer_id}: negative audit Boltzmann weight"
                        )
                    record["audit_boltzmann_weight_all_valid"] = audit_weight
                    if selected_for_nmr:
                        selected_text = audit_row.get(
                            "boltzmann_weight_selected_normalized", ""
                        ).strip()
                        if not selected_text:
                            raise ProcessingError(
                                f"{state_label}/{conformer_id}: selected_for_nmr=1 but selected-normalized "
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
            raise ProcessingError(f"No valid QRRHO minima for {state_label}")
        temperatures = [float(record["thermochemistry_temperature_K"]) for record in parsed_valid]
        if max(temperatures) - min(temperatures) > 1.0e-3:
            raise ProcessingError(f"Inconsistent thermochemistry temperatures for {state_label}")
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
                    f"{state_label}/{record['conformer_id']}: recomputed all-minimum QRRHO weight "
                    f"{record['boltzmann_weight_all_valid_qrrho']:.9g} differs from DFT audit "
                    f"{audit_weight:.9g}"
                )
        audit_weight_sum = sum(
            float(record["audit_boltzmann_weight_all_valid"]) for record in parsed_valid
        )
        if abs(audit_weight_sum - 1.0) > WEIGHT_TOLERANCE:
            raise ProcessingError(
                f"{state_label}: unique-valid DFT audit weights sum to {audit_weight_sum:.9g}, not 1"
            )
        selected_valid = [record for record in parsed_valid if record["selected_for_nmr"]]
        if not selected_valid:
            raise ProcessingError(f"DFT audit selects no unique minimum for NMR: {state_label}")
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
                    f"{state_label}/{record['conformer_id']}: DFT audit selected-normalized weight "
                    f"{supplied_selected:.9g} differs from recomputed value {expected_selected:.9g}"
                )
            record["dft_audit_validation_status"] = (
                "validated_unique_minimum_and_weights"
            )
        grouped[state_key] = all_records
    return grouped, audit


def validate_property_manifest_weights(
    figure_id: str,
    manifest_rows: dict[str, dict[str, str]],
    dft_records: dict[str, dict[str, Any]],
    weight_sets: dict[str, dict[str, Any]],
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
        f"Explicit-solvent scenario: {compound['scenario_id']}",
        f"Formula: {compound['formula']}",
        f"InChIKey: {compound['inchikey']}",
        "",
        "This file contains calculated ORCA predictions. It is not an experimental NMR report,",
        "does not establish identity or purity, and must not be presented as measured characterization.",
        "",
        "QRRHO conformer weighting",
        "--------------------------",
        f"Unique valid DFT minima: {weighting['n_valid_dft_conformers']}",
        f"Post-DFT duplicates excluded: {weighting['n_post_dft_duplicates_excluded']}",
    ]
    for nucleus, weight_set in sorted(weighting["nucleus_weight_sets"].items()):
        coverage = weight_set.get("completed_population")
        lines.append(
            f"{nucleus}: T={weight_set['temperature_K']:.4f} K; completed population="
            + ("not requested" if coverage is None else f"{float(coverage):.8f}")
        )
    lines.extend(
        [
        "",
        "Calculated chemical shifts",
        "---------------------------",
        "Raw values use the same-level nucleus-specific reference (TMS for H/C, CFCl3 for F, H3PO4 for P).",
        "No empirical H/C regression is applied. A regression fitted for a different method/basis is not transferable.",
        "F/P values remain PENDING when their reference output is absent or invalid.",
        "The reported conformer SD is population spread, not predictive error.",
        "",
        "Atom       sigma_avg / ppm   delta_same_method_reference / ppm",
        ]
    )
    for row in bundle["averaged_nuclei"]:
        exchange_note = (
            " [exchangeable; shift/linewidth/integral provisional]"
            if int(row.get("exchangeable_H_flag", 0))
            else ""
        )
        lines.append(
            f"{row['atom_label_1based']:<10s} {row['boltzmann_averaged_shielding_ppm']:>15.6f} "
            f"{format_optional_number(row['raw_reference_shift_ppm']):>35s}"
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
    scenario_id: str,
    dft_records: dict[str, dict[str, Any]],
    weight_sets: dict[str, dict[str, Any]],
    shielding_manifest: dict[str, dict[str, str]],
    coupling_manifest: dict[str, dict[str, str]],
    shielding_method: str,
    shielding_method_fingerprint: str,
    references: dict[str, dict[str, Any]],
    tms_output: Path,
    cfcl3_output: Path,
    h3po4_output: Path,
    output_root: Path,
    allow_partial: bool,
    overwrite: bool,
    source_manifest_paths: dict[str, Path],
    shielding_state_root: Path | None = None,
    coupling_state_root: Path | None = None,
    reference_approval: dict[str, Any] | None = None,
    compound_metadata: dict[str, str] | None = None,
) -> dict[str, Any]:
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", scenario_id):
        raise ProcessingError(f"Unsafe scenario_id for output path: {scenario_id!r}")
    output_directory = output_root / figure_id / scenario_id
    if output_directory.exists() and any(output_directory.iterdir()) and not overwrite:
        raise ProcessingError(
            f"Output directory already contains files: {output_directory}; pass --overwrite to replace this script's files"
        )

    warnings: list[str] = []
    if not shielding_manifest:
        raise ProcessingError(f"No shielding jobs in manifest for {figure_id}/{scenario_id}")
    if not coupling_manifest:
        raise ProcessingError(f"No coupling jobs in manifest for {figure_id}/{scenario_id}")
    atom_map_path, atom_map, atom_map_sha256, atom_map_paths = load_scenario_atom_map(
        figure_id, scenario_id, shielding_manifest, coupling_manifest
    )
    state_label = f"{figure_id}/{scenario_id}"
    weight_gate_provenance = validate_weight_gate_binding(
        figure_id,
        scenario_id,
        weight_sets,
        shielding_manifest,
        coupling_manifest,
    )
    dft_state_binding = validate_property_dft_state_binding(
        dft_records,
        shielding_manifest,
        coupling_manifest,
    )
    isotope_provenance = validate_scenario_isotope_maps(
        figure_id, scenario_id, atom_map, shielding_manifest, coupling_manifest
    )

    parsed_shieldings: dict[str, list[dict[str, Any]]] = {}
    parsed_couplings: dict[str, list[dict[str, Any]]] = {}
    shielding_runner_records: dict[str, dict[str, Any]] = {}
    coupling_runner_records: dict[str, dict[str, Any]] = {}
    shielding_method_validations: dict[str, dict[str, str]] = {}
    source_files: list[dict[str, Any]] = [
        *[file_record(path) for path in atom_map_paths],
        *[file_record(path) for path in isotope_provenance["paths"]],
        file_record(tms_output),
        file_record(cfcl3_output),
        file_record(h3po4_output),
    ]

    for conformer_id, row in sorted(shielding_manifest.items()):
        runner = validate_runner_output(row, "shielding", PROJECT, shielding_state_root)
        output = Path(runner["output_path"])
        xyz = property_coordinate_path(row)
        try:
            property_provenance = validate_property_provenance(row)
            method_validation = validate_nmr_method_evidence(
                Path(row["job_directory"]) / row["input_file"],
                output,
                row["property_branch"],
                row["method_fingerprint"],
            )
            validate_atom_order(xyz, atom_map)
            parsed_all = parse_shielding_summary(output)
            parsed: list[dict[str, Any]] = []
            for shielding in parsed_all:
                index = int(shielding["nucleus_index_0based"])
                if index not in atom_map:
                    raise ProcessingError(f"Shielding nucleus {index} absent from atom map")
                if shielding["element"] != atom_map[index]["element"]:
                    raise ProcessingError(
                        f"Shielding element mismatch at atom {index}: {shielding['element']} vs {atom_map[index]['element']}"
                    )
                atom = atom_map[index]
                if (
                    int(atom.get("environment_atom", 0)) == 1
                    or int(atom.get("include_in_analyte_nmr", 0)) != 1
                    or atom.get("nmr_scope") != "ANALYTE"
                ):
                    warnings.append(
                        f"Excluded non-analyte shielding atom {index} from {state_label}/{conformer_id}"
                    )
                    continue
                parsed.append(shielding)
            if not parsed:
                raise ProcessingError("Shielding output contains no analyte nuclei")
            parsed_shieldings[conformer_id] = parsed
            shielding_runner_records[conformer_id] = runner
            shielding_method_validations[conformer_id] = method_validation
            source_files.extend(
                [
                    file_record(output),
                    file_record(xyz),
                    file_record(Path(property_provenance["path"])),
                    file_record(Path(runner["runner_state"])),
                ]
            )
        except Exception as exc:
            message = (
                f"{figure_id}/{scenario_id}/{conformer_id} shielding unavailable: {exc!r}"
            )
            warnings.append(message)
            if not allow_partial:
                raise ProcessingError(message) from exc

    for conformer_id, row in sorted(coupling_manifest.items()):
        runner = validate_runner_output(row, "jcoupling", PROJECT, coupling_state_root)
        output = Path(runner["output_path"])
        property_path = property_path_from_output(output)
        xyz = property_coordinate_path(row)
        try:
            property_provenance = validate_property_provenance(row)
            validate_atom_order(xyz, atom_map)
            parsed_all = parse_j_pairs(property_path, output)
            parsed: list[dict[str, Any]] = []
            for pair in parsed_all:
                pair_is_analyte = True
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
                    atom = atom_map[index]
                    if (
                        int(atom.get("environment_atom", 0)) == 1
                        or int(atom.get("include_in_analyte_nmr", 0)) != 1
                        or atom.get("nmr_scope") != "ANALYTE"
                    ):
                        pair_is_analyte = False
                if pair_is_analyte:
                    parsed.append(pair)
                else:
                    warnings.append(
                        f"Excluded coupling touching a non-analyte atom from {state_label}/{conformer_id}"
                    )
            if not parsed:
                raise ProcessingError("Coupling output contains no analyte-analyte pairs")
            parsed_couplings[conformer_id] = parsed
            coupling_runner_records[conformer_id] = runner
            source_files.extend(
                [
                    file_record(output),
                    file_record(property_path),
                    file_record(xyz),
                    file_record(Path(property_provenance["path"])),
                    file_record(Path(runner["runner_state"])),
                ]
            )
        except Exception as exc:
            message = (
                f"{figure_id}/{scenario_id}/{conformer_id} coupling unavailable: {exc!r}"
            )
            warnings.append(message)
            if not allow_partial:
                raise ProcessingError(message) from exc

    if not parsed_shieldings:
        raise ProcessingError(f"No completed shielding results for {figure_id}/{scenario_id}")
    if not parsed_couplings:
        raise ProcessingError(f"No completed coupling results for {figure_id}/{scenario_id}")

    shielding_nuclei = {
        NUCLEUS_FROM_ELEMENT[row["element"]]
        for rows in parsed_shieldings.values()
        for row in rows
        if row["element"] in NUCLEUS_FROM_ELEMENT
    }
    shield_contexts = {
        nucleus: weight_context(
            weight_sets[nucleus], set(shielding_manifest), set(parsed_shieldings)
        )
        for nucleus in shielding_nuclei
    }
    coupling_weight_nucleus = "1H" if "1H" in weight_sets else "13C"
    coupling_context = weight_context(
        weight_sets[coupling_weight_nucleus],
        set(coupling_manifest),
        set(parsed_couplings),
    )
    for nucleus, context in shield_contexts.items():
        if context["completed_fraction"] < 1.0 - 1.0e-10:
            warnings.append(
                f"{nucleus} shieldings were renormalized over an incomplete selected-conformer subset"
            )
    if coupling_context["completed_fraction"] < 1.0 - 1.0e-10:
        warnings.append("J couplings were renormalized over an incomplete selected-conformer subset")

    shielding_raw_rows: list[dict[str, Any]] = []
    shielding_by_nucleus: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for conformer_id, rows in sorted(parsed_shieldings.items()):
        for row in rows:
            index = int(row["nucleus_index_0based"])
            atom = atom_map[index]
            nucleus = NUCLEUS_FROM_ELEMENT[atom["element"]]
            context = shield_contexts[nucleus]
            weight_record = weight_sets[nucleus]["conformers"][conformer_id]
            normalized_weight = context["normalized_completed_weights"][conformer_id]
            combined = {
                "notice": NOTICE,
                "figure_id": figure_id,
                "scenario_id": scenario_id,
                "conformer_id": conformer_id,
                "weight_nucleus": nucleus,
                "weight_temperature_K": weight_sets[nucleus]["temperature_K"],
                "relative_qrrho_gibbs_kcal_mol": weight_record[
                    "relative_qrrho_gibbs_kcal_mol"
                ],
                "boltzmann_weight_authoritative": weight_record["weight"],
                "boltzmann_weight_completed_shielding_subset": normalized_weight,
                **atom,
                "atom_label_1based": f"{atom['element']}{atom['nucleus_index_1based']}",
                "isotropic_shielding_ppm": row["isotropic_shielding_ppm"],
                "anisotropy_ppm": row["anisotropy_ppm"],
                "source_output": str(shielding_runner_records[conformer_id]["output_path"]),
            }
            shielding_raw_rows.append(combined)
            shielding_by_nucleus[index].append(combined)

    averaged_nuclei: list[dict[str, Any]] = []
    for index, rows in sorted(shielding_by_nucleus.items()):
        atom = atom_map[index]
        element = atom["element"]
        nucleus = NUCLEUS_FROM_ELEMENT.get(element)
        if nucleus is None or nucleus not in shield_contexts:
            warnings.append(f"No authoritative nucleus weight set for {element}{index + 1}")
            continue
        shield_context = shield_contexts[nucleus]
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
        flags = [NOTICE, "NO_EMPIRICAL_REGRESSION_APPLIED"]
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
        if shield_context["completed_population"] < 0.9999:
            flags.append("TRUNCATED_OR_INCOMPLETE_CONFORMER_POPULATION")
        averaged_nuclei.append(
            {
                "notice": NOTICE,
                "figure_id": figure_id,
                "scenario_id": scenario_id,
                **atom,
                "atom_label_1based": f"{element}{atom['nucleus_index_1based']}",
                "n_conformers_with_shielding": len(rows),
                "weight_nucleus": nucleus,
                "weight_temperature_K": shield_context["temperature_K"],
                "completed_population_of_authoritative_weight_set": shield_context[
                    "completed_population"
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
                "predictive_uncertainty_status": (
                    "not quantified for this workflow; any reported conformer SD is population spread only. "
                    "Reference-model errors are not validated uncertainty estimates"
                ),
                "flags": ";".join(flags),
            }
        )

    coupling_raw_rows: list[dict[str, Any]] = []
    coupling_by_pair: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for conformer_id, rows in sorted(parsed_couplings.items()):
        normalized_weight = coupling_context["normalized_completed_weights"][conformer_id]
        weight_record = weight_sets[coupling_weight_nucleus]["conformers"][conformer_id]
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
                "scenario_id": scenario_id,
                "conformer_id": conformer_id,
                "weight_nucleus": coupling_weight_nucleus,
                "weight_temperature_K": weight_sets[coupling_weight_nucleus]["temperature_K"],
                "relative_qrrho_gibbs_kcal_mol": weight_record[
                    "relative_qrrho_gibbs_kcal_mol"
                ],
                "boltzmann_weight_authoritative": weight_record["weight"],
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
                        Path(coupling_runner_records[conformer_id]["output_path"])
                    )
                ),
            }
            coupling_raw_rows.append(combined)
            coupling_by_pair[(index_a, index_b)].append(combined)

    averaged_couplings: list[dict[str, Any]] = []
    expected_coupling_ids = set(parsed_couplings)
    for (index_a, index_b), rows in sorted(coupling_by_pair.items()):
        available_ids = {str(row["conformer_id"]) for row in rows}
        pair_population = math.fsum(
            weight_sets[coupling_weight_nucleus]["conformers"][conformer_id]["weight"]
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
        full_weights = [float(row["boltzmann_weight_completed_coupling_subset"]) for row in rows]
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
                "scenario_id": scenario_id,
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
                "weight_nucleus": coupling_weight_nucleus,
                "weight_temperature_K": coupling_context["temperature_K"],
                "pair_population_of_authoritative_weight_set": pair_population,
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

    valid_dft = list(dft_records.values())
    n_post_dft_duplicates_excluded = 0
    valid_dft.sort(
        key=lambda row: min(row["relative_qrrho_gibbs_by_nucleus_kcal_mol"].values())
    )
    first_dft = valid_dft[0]
    first_shield = next(iter(shielding_manifest.values()))
    first_coupling = next(iter(coupling_manifest.values()))
    metadata = compound_metadata or {}
    compound = {
        "figure_id": figure_id,
        "scenario_id": scenario_id,
        "compound_number": int(first_dft["compound_number"]),
        "display_name": first_dft["display_name"],
        "formula": first_dft.get("formula") or metadata.get("formula", ""),
        "isomeric_smiles": first_dft.get("isomeric_smiles") or metadata.get(
            "isomeric_smiles", ""
        ),
        "inchikey": first_dft.get("inchikey") or metadata.get("inchikey", ""),
        "charge": int(first_dft["charge"]),
        "multiplicity": int(first_dft["multiplicity"]),
    }
    conformer_rows: list[dict[str, Any]] = []
    for conformer_id, row in sorted(
        dft_records.items(),
        key=lambda pair: (
            min(pair[1]["relative_qrrho_gibbs_by_nucleus_kcal_mol"].values()),
            pair[0],
        ),
    ):
        conformer_rows.append(
            {
                "notice": NOTICE,
                "figure_id": figure_id,
                "scenario_id": scenario_id,
                "conformer_id": conformer_id,
                "weighting_status": row.get("weighting_status"),
                "weights_by_nucleus": json.dumps(
                    row["weights_by_nucleus"], sort_keys=True, separators=(",", ":")
                ),
                "qrrho_gibbs_by_nucleus_hartree": json.dumps(
                    row["qrrho_gibbs_by_nucleus_hartree"],
                    sort_keys=True,
                    separators=(",", ":"),
                ),
                "relative_qrrho_gibbs_by_nucleus_kcal_mol": json.dumps(
                    row["relative_qrrho_gibbs_by_nucleus_kcal_mol"],
                    sort_keys=True,
                    separators=(",", ":"),
                ),
                "qrrho_verified": row.get("qrrho_verified"),
                "thermochemistry_temperature_K": row.get("thermochemistry_temperature_K"),
                "gibbs_free_energy_hartree": row.get("gibbs_free_energy_hartree"),
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
                "runner_state": row.get("runner_provenance", {}).get("runner_state"),
                "runner_state_sha256": row.get("runner_provenance", {}).get(
                    "runner_state_sha256"
                ),
            }
        )

    assumptions = [
        NOTICE,
        "All reported values are calculations and cannot replace measured NMR characterization, identity or purity evidence.",
        "NMR_conformer_weights.csv is the authoritative population source. Weight sets are separately hash-traced by explicit-solvent scenario, nucleus and locked acquisition temperature; this processor never recomputes them from a legacy single-temperature audit.",
        "For benzamide, 1H weights are locked to 312.609 K and 13C weights to 304.0651 K. F3_015 uses its explicitly locked provisional condition until raw acquisition metadata are available.",
        "Property averages are renormalized over completed conformers listed in the corresponding property manifest; population coverage is reported explicitly.",
        "Boltzmann populations and every reported average are restricted to one explicit-solvent scenario. Different solvent stoichiometries are never combined.",
        "Raw chemical shifts use delta = mean same-level nucleus-specific reference shielding minus Boltzmann-averaged sample shielding, with the stated reference convention.",
        "No empirical regression is applied. In particular, the former external PBE0/pcS-2 fit is disabled because it is not method-matched to revDSD-PBEP86-D4/2021/pcSseg-3 or PBE0/pcSseg-3.",
        "VPT2 and any vibrational shielding correction are diagnostic branches only. Unless a separately approved same-level sample/reference correction exists, no VPT2 term is allowed in the main chemical-shift delta.",
        "If the corresponding nucleus-specific reference output is absent or invalid, the calculated shielding is retained and the raw shift is explicitly pending.",
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
            "shielding_method": shielding_method,
            "shielding_method_fingerprint": shielding_method_fingerprint,
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
            "method_validation_by_conformer": shielding_method_validations,
            "vpt2_main_delta_policy": "FORBIDDEN_DIAGNOSTIC_ONLY",
        },
        "condition_state": {
            "figure_id": figure_id,
            "scenario_id": scenario_id,
            "state_isolation_status": "PASS_NO_CROSS_SCENARIO_AVERAGING",
            "job_local_atom_map_sha256": atom_map_sha256,
            "job_local_atom_map_copies_checked": len(atom_map_paths),
            "isotope_mass_mapping": {
                key: value
                for key, value in isotope_provenance.items()
                if key != "paths"
            },
            "weight_gate_binding": weight_gate_provenance,
            "property_to_dft_runner_state_binding": dft_state_binding,
        },
        "weighting": {
            "source": "data/NMR_conformer_weights.csv",
            "policy": "authoritative scenario/nucleus/temperature-specific weights; no processor-side reweighting",
            "n_valid_dft_conformers": len(valid_dft),
            "n_post_dft_duplicates_excluded": n_post_dft_duplicates_excluded,
            "nucleus_weight_sets": {
                nucleus: {
                    "temperature_K": weight_set["temperature_K"],
                    "condition_status": weight_set["condition_status"],
                    "condition_source_note": weight_set["condition_source_note"],
                    "approved_gate_sha256": weight_set["approved_gate_sha256"],
                    "completed_population": shield_contexts.get(
                        nucleus, {"completed_population": None}
                    )["completed_population"],
                }
                for nucleus, weight_set in weight_sets.items()
            },
            "coupling_weight_nucleus": coupling_weight_nucleus,
            "coupling_completed_population": coupling_context["completed_population"],
        },
        "reference_shieldings": references,
        "reference_registry_approval": reference_approval,
        "empirical_regression_policy": {
            "status": "disabled",
            "reason": "No independently validated, exact-method regression is bundled.",
        },
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
    prefix = f"{figure_id}__{scenario_id}"
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
            "scenario_id": scenario_id,
            "rows": coupling_raw_rows,
        },
    )
    atomic_write_json(
        output_directory / f"{prefix}_calculated_scalar_couplings.json",
        {
            "schema_version": SCHEMA_VERSION,
            "notice": NOTICE,
            "figure_id": figure_id,
            "scenario_id": scenario_id,
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
        "scenario_id": scenario_id,
        "display_name": compound["display_name"],
        "status": "processed_calculated_prediction_not_experimental",
        "shielding_method": shielding_method,
        "shielding_method_fingerprint": shielding_method_fingerprint,
        "n_valid_dft_conformers": len(valid_dft),
        "n_post_dft_duplicates_excluded": n_post_dft_duplicates_excluded,
        "n_averaged_nuclei": len(averaged_nuclei),
        "n_exchangeable_H": sum(
            int(row.get("exchangeable_H_flag", 0)) for row in averaged_nuclei
        ),
        "n_averaged_scalar_coupling_pairs": len(averaged_couplings),
        "n_averaged_JHH_pairs": len(averaged_jhh),
        "shielding_population_coverage_by_nucleus": json.dumps(
            {
                nucleus: context["completed_population"]
                for nucleus, context in shield_contexts.items()
            },
            sort_keys=True,
            separators=(",", ":"),
        ),
        "coupling_population_coverage": coupling_context["completed_population"],
        "n_warnings": len(warnings),
        "output_directory": str(output_directory),
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description=(
            "Parse ORCA 6.1.1 shielding/J outputs and apply approved nucleus-specific "
            "QRRHO-Boltzmann weights to "
            "NMR predictions. Outputs are explicitly labelled calculated, not experimental."
        )
    )
    result.add_argument("--dft-manifest", type=Path, default=DEFAULT_DFT_MANIFEST)
    result.add_argument(
        "--weights",
        type=Path,
        default=DEFAULT_WEIGHT_TABLE,
        help="Authoritative scenario/nucleus/temperature-specific NMR_conformer_weights.csv.",
    )
    result.add_argument(
        "--nucleus-conditions", type=Path, default=DEFAULT_NUCLEUS_CONDITIONS
    )
    result.add_argument(
        "--compound-registry", type=Path, default=DEFAULT_COMPOUND_REGISTRY
    )
    result.add_argument("--shielding-manifest", type=Path, default=DEFAULT_SHIELDING_MANIFEST)
    result.add_argument(
        "--shielding-method",
        default=DEFAULT_SHIELDING_METHOD,
        help=(
            "property_branch to process. The default is the production revDSD branch; "
            "different shielding branches are never averaged together."
        ),
    )
    result.add_argument(
        "--coupling-manifest",
        type=Path,
        default=DEFAULT_COUPLING_MANIFEST,
        help=(
            "Scalar-coupling manifest (default: NMR_j_coupling_job_manifest.csv). "
            f"A legacy H-H-only manifest remains available only when passed explicitly, e.g. {LEGACY_HH_COUPLING_MANIFEST.name}."
        ),
    )
    result.add_argument(
        "--reference-manifest",
        type=Path,
        default=DEFAULT_REFERENCE_MANIFEST,
        help="Reference-job manifest used to lock TMS to the selected shielding fingerprint.",
    )
    result.add_argument(
        "--reference-registry",
        type=Path,
        default=DEFAULT_REFERENCE_REGISTRY,
        help="Named-reviewer-approved reference registry with output hashes.",
    )
    result.add_argument("--dft-state-root", type=Path)
    result.add_argument("--shielding-state-root", type=Path)
    result.add_argument("--coupling-state-root", type=Path)
    result.add_argument("--reference-state-root", type=Path)
    result.add_argument(
        "--tms-output",
        type=Path,
        default=None,
        help=(
            "Optional relocated TMS output. If omitted, the output is resolved from the "
            "same-method row in --reference-manifest. This override never bypasses the "
            "sample/reference fingerprint check."
        ),
    )
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
        "--scenario-id",
        nargs="+",
        help=(
            "Optional explicit-solvent scenario IDs. Each selected scenario is still "
            "processed and written independently; cross-scenario averaging is forbidden."
        ),
    )
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
        "weight_table": args.weights.expanduser().resolve(),
        "nucleus_conditions": args.nucleus_conditions.expanduser().resolve(),
        "compound_registry": args.compound_registry.expanduser().resolve(),
        "shielding_manifest": args.shielding_manifest.expanduser().resolve(),
        "coupling_manifest": args.coupling_manifest.expanduser().resolve(),
        "reference_manifest": args.reference_manifest.expanduser().resolve(),
        "reference_registry": args.reference_registry.expanduser().resolve(),
    }
    cfcl3_output = args.cfcl3_output.expanduser().resolve()
    h3po4_output = args.h3po4_output.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()
    try:
        dft_rows = read_manifest(
            manifest_paths["dft_manifest"],
            {
                "job_id",
                "figure_id",
                "scenario_id",
                "compound_number",
                "display_name",
                "charge",
                "multiplicity",
                "job_directory",
                "input_file",
                "output_file",
                "input_sha256",
            },
        )
        weight_rows = read_csv_table(
            manifest_paths["weight_table"],
            {
                "figure_id",
                "scenario_id",
                "nucleus",
                "temperature_K",
                "conformer_id",
                "qrrho_gibbs_free_energy_hartree",
                "weight",
                "approved_gate_sha256",
            },
        )
        locked_conditions = load_locked_nucleus_conditions(
            manifest_paths["nucleus_conditions"]
        )
        compound_rows = read_csv_table(
            manifest_paths["compound_registry"],
            {"figure_id", "formula", "isomeric_smiles", "inchikey"},
        )
        compounds = {row["figure_id"]: row for row in compound_rows}
        shielding_rows = read_manifest(
            manifest_paths["shielding_manifest"],
            {
                "job_id",
                "figure_id",
                "scenario_id",
                "conformer_id",
                "job_directory",
                "input_file",
                "output_file",
                "property_branch",
                "method_fingerprint",
                "input_sha256",
                "conformer_xyz_sha256",
                "nmr_atom_map_sha256",
                "source_isotope_map_sha256",
                "approved_minimum_gate_sha256",
                "source_opt_freq_state_sha256",
            },
        )
        coupling_rows = read_manifest(
            manifest_paths["coupling_manifest"],
            {
                "job_id",
                "figure_id",
                "scenario_id",
                "conformer_id",
                "job_directory",
                "input_file",
                "output_file",
                "input_sha256",
                "conformer_xyz_sha256",
                "nmr_atom_map_sha256",
                "source_isotope_map_sha256",
                "approved_minimum_gate_sha256",
                "source_opt_freq_state_sha256",
            },
        )
        shielding_rows, shielding_method_fingerprint = select_shielding_method_rows(
            shielding_rows, args.shielding_method
        )
        _unapproved_default_tms, tms_reference_row = resolve_same_method_tms_output(
            manifest_paths["reference_manifest"],
            args.shielding_method,
            shielding_method_fingerprint,
            None,
        )
        tms_output, reference_approval = validate_approved_reference_registry(
            manifest_paths["reference_registry"],
            tms_reference_row,
            args.shielding_method,
            shielding_method_fingerprint,
            args.tms_output,
            args.reference_state_root,
        )
        shielding_grouped = group_manifest_by_state_and_conformer(shielding_rows)
        coupling_grouped = group_manifest_by_state_and_conformer(coupling_rows)
        shielding_states = set(shielding_grouped)
        coupling_states = set(coupling_grouped)
        if shielding_states != coupling_states:
            raise ProcessingError(
                "Shielding/coupling scenario coverage differs; "
                f"shielding_only={sorted(shielding_states - coupling_states)}, "
                f"coupling_only={sorted(coupling_states - shielding_states)}"
            )
        available = shielding_states
        selected = select_state_keys(available, args.only, args.scenario_id)
        weight_sets = load_nucleus_weight_sets(
            weight_rows, selected, locked_conditions
        )
        dft_grouped, dft_audit = prepare_dft_records_from_weight_sets(
            dft_rows,
            selected,
            weight_sets,
            args.significant_imaginary_cm,
            args.dft_state_root,
        )
        references = reference_summary(tms_output, cfcl3_output, h3po4_output)
        summaries: list[dict[str, Any]] = []
        for figure_id, scenario_id in sorted(selected):
            state_key = (figure_id, scenario_id)
            summaries.append(
                process_figure(
                    figure_id=figure_id,
                    scenario_id=scenario_id,
                    dft_records=dft_grouped[state_key],
                    weight_sets=weight_sets[state_key],
                    shielding_manifest=shielding_grouped[state_key],
                    coupling_manifest=coupling_grouped[state_key],
                    shielding_method=args.shielding_method,
                    shielding_method_fingerprint=shielding_method_fingerprint,
                    references=references,
                    tms_output=tms_output,
                    cfcl3_output=cfcl3_output,
                    h3po4_output=h3po4_output,
                    output_root=output_root,
                    allow_partial=args.allow_partial,
                    overwrite=args.overwrite,
                    source_manifest_paths=manifest_paths,
                    shielding_state_root=args.shielding_state_root,
                    coupling_state_root=args.coupling_state_root,
                    reference_approval=reference_approval,
                    compound_metadata=compounds.get(figure_id),
                )
            )
    except (
        OSError,
        ValueError,
        ProcessingError,
        MethodValidationError,
        RunnerProvenanceError,
    ) as exc:
        print(f"Processing failed: {exc}")
        return 2

    label = (
        "ALL_STATES"
        if not args.only and not args.scenario_id
        else "__".join(f"{figure_id}_{scenario_id}" for figure_id, scenario_id in sorted(selected))
    )
    summary_payload = {
        "schema_version": SCHEMA_VERSION,
        "notice": NOTICE,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "selected_states": [
            {"figure_id": figure_id, "scenario_id": scenario_id}
            for figure_id, scenario_id in sorted(selected)
        ],
        "shielding_method": args.shielding_method,
        "shielding_method_fingerprint": shielding_method_fingerprint,
        "tms_reference_job_id": tms_reference_row["job_id"],
        "tms_reference_output": str(tms_output),
        "results": summaries,
        "dft_runner_and_weight_binding_audit": [
            {
                "figure_id": row.get("figure_id"),
                "scenario_id": row.get("scenario_id"),
                "conformer_id": row.get("conformer_id"),
                "weighting_status": row.get("weighting_status"),
                "weights_by_nucleus": row.get("weights_by_nucleus"),
                "dft_output_sha256": row.get("dft_output_sha256"),
                "runner_state_sha256": row.get("runner_provenance", {}).get(
                    "runner_state_sha256"
                ),
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
    print(
        f"Processed {len(summaries)} independently isolated compound/scenario state(s): "
        + ", ".join(f"{figure_id}/{scenario_id}" for figure_id, scenario_id in sorted(selected))
    )
    print(output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
