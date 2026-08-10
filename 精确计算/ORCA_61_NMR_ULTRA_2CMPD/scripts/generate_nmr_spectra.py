from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import uuid
from collections import deque
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_ROOT = PROJECT / "NMR_predictions_CALCULATED_NOT_EXPERIMENTAL"
DEFAULT_OUTPUT_ROOT = PROJECT / "NMR_spectra_CALCULATED_NOT_EXPERIMENTAL"
DEFAULT_SOLUTION_EQUIVALENCE = PROJECT / "config" / "solution_equivalence_groups.json"
DEFAULT_APPEARANCE_CONFIG = PROJECT / "config" / "appearance_profiles.json"

NOTICE = "CALCULATED PREDICTION — NOT EXPERIMENTAL"
SCHEMA_VERSION = "1.3.0"
NO_CARBON_FIGURE_IDS = {"F3_035", "F3_043", "F3_044"}
EXCHANGEABLE_H_PARENT_ELEMENTS = {"N", "O", "S"}
H1_MATH_LABEL = r"$^{1}\mathrm{H}$"
C13_MATH_LABEL = r"$^{13}\mathrm{C}$"
C13_H1_DECOUPLED_MATH_LABEL = r"$^{13}\mathrm{C}\{^{1}\mathrm{H}\}$"

VARIANTS: dict[str, dict[str, str]] = {
    "raw_TMS": {
        "shift_field": "raw_tms_reference_shift_ppm",
        "label": "raw TMS-referenced calculated shifts",
        "short_label": "Raw TMS reference",
        "color": "#2C6E9F",
    },
    "provisional_regression": {
        "shift_field": "provisional_external_regression_shift_ppm",
        "label": "provisional external-regression calculated shifts",
        "short_label": "Provisional regression",
        "color": "#C46D2D",
    },
}


class SpectrumError(RuntimeError):
    """Raised when a prediction bundle cannot support an auditable spectrum."""


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
    atomic_write_text(
        path,
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
    )


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


def field_union(rows: Sequence[dict[str, Any]]) -> list[str]:
    """Preserve existing first-row order while retaining conditional audit fields."""
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for field in row:
            if field not in seen:
                fields.append(field)
                seen.add(field)
    return fields


def safe_token(value: str) -> str:
    token = "".join(character if character.isalnum() or character in "_.-" else "_" for character in value)
    token = token.strip("_.")
    if not token or token in {".", ".."}:
        raise SpectrumError(f"Unsafe output identifier: {value!r}")
    return token


def available_variant_keys(bundle: dict[str, Any]) -> list[str]:
    """Return only shift variants actually present for every analyte H/C row."""

    analyte_rows = [
        row for row in bundle.get("averaged_nuclei", []) if row.get("element") in {"H", "C"}
    ]
    if not analyte_rows:
        raise SpectrumError("Prediction bundle contains no analyte 1H/13C rows")
    available: list[str] = []
    for variant_key, variant in VARIANTS.items():
        field = variant["shift_field"]
        try:
            complete = all(
                row.get(field) is not None and math.isfinite(float(row[field]))
                for row in analyte_rows
            )
        except (TypeError, ValueError):
            complete = False
        if complete:
            available.append(variant_key)
    if not available:
        raise SpectrumError("No complete calculated shift variant is present in the prediction bundle")
    return available


def atomic_write_profile(
    path: Path,
    x_ppm: np.ndarray,
    intensity: np.ndarray,
    metadata: dict[str, Any],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{uuid.uuid4().hex}")
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(f"# {NOTICE}\n")
            for key, value in metadata.items():
                handle.write(f"# {key}: {value}\n")
            handle.write("chemical_shift_ppm\trelative_intensity\n")
            for shift, value in zip(x_ppm, intensity):
                handle.write(f"{shift:.8f}\t{value:.10g}\n")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def configure_matplotlib() -> str:
    try:
        arial_path = font_manager.findfont("Arial", fallback_to_default=False)
    except ValueError as exc:
        raise SpectrumError("Arial is required but was not found by matplotlib") from exc
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7.0,
            "axes.linewidth": 0.75,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )
    return arial_path


def atomic_save_figure(fig: plt.Figure, stem: Path, dpi: int) -> list[Path]:
    outputs: list[Path] = []
    for suffix, format_name in ((".png", "png"), (".svg", "svg"), (".pdf", "pdf")):
        path = stem.with_suffix(suffix)
        temporary = path.with_name(f".{path.stem}.tmp.{uuid.uuid4().hex}{suffix}")
        try:
            fig.savefig(
                temporary,
                format=format_name,
                dpi=dpi if format_name == "png" else None,
                bbox_inches="tight",
                facecolor="white",
            )
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)
        outputs.append(path)
    return outputs


def load_bundle(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)
    bundle = json.loads(path.read_text(encoding="utf-8"))
    if bundle.get("notice") != NOTICE:
        raise SpectrumError(f"Prediction bundle lacks the required non-experimental notice: {path}")
    for key in ("compound", "averaged_nuclei", "averaged_JHH", "weighting"):
        if key not in bundle:
            raise SpectrumError(f"Prediction bundle lacks {key}: {path}")
    if not isinstance(bundle["averaged_nuclei"], list) or not isinstance(
        bundle["averaged_JHH"], list
    ):
        raise SpectrumError(f"Prediction bundle has malformed NMR arrays: {path}")
    if "averaged_scalar_couplings" in bundle and not isinstance(
        bundle["averaged_scalar_couplings"], list
    ):
        raise SpectrumError(
            f"Prediction bundle has malformed averaged_scalar_couplings: {path}"
        )
    return bundle


def apply_declared_solution_equivalence(
    source_bundle: dict[str, Any], config_path: Path
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Apply only pre-declared, identity-locked solution symmetry operations.

    Explicit solvent clusters may break an analyte's solution-averaged symmetry.
    The raw atom-resolved bundle remains unchanged on disk; this function works
    on a copy and records every averaged shift and J orbit for audit.
    """
    config = json.loads(config_path.read_text(encoding="utf-8"))
    figure_id = str(source_bundle.get("compound", {}).get("figure_id", ""))
    compound_rules = config.get("compounds", {}).get(figure_id)
    if compound_rules is None:
        return copy.deepcopy(source_bundle), []
    if not isinstance(compound_rules, dict):
        raise SpectrumError(f"Malformed solution-equivalence entry for {figure_id}")
    expected_inchikey = str(compound_rules.get("inchikey", ""))
    observed_inchikey = str(source_bundle.get("compound", {}).get("inchikey", ""))
    if observed_inchikey and expected_inchikey and observed_inchikey != expected_inchikey:
        raise SpectrumError(f"Solution-equivalence identity mismatch for {figure_id}")
    bundle = copy.deepcopy(source_bundle)
    nuclei = {
        int(row["nucleus_index_0based"]): row for row in bundle["averaged_nuclei"]
    }
    audit: list[dict[str, Any]] = []
    for group in compound_rules.get("groups", []):
        if not isinstance(group, dict):
            raise SpectrumError(f"Malformed solution-equivalence group for {figure_id}")
        members = [int(value) for value in group.get("atom_indices_0based", [])]
        element = str(group.get("element", ""))
        if len(members) < 2 or len(set(members)) != len(members):
            raise SpectrumError(f"Invalid member list for {group.get('group_id')}")
        if any(index not in nuclei or nuclei[index].get("element") != element for index in members):
            raise SpectrumError(f"Equivalence group element/index mismatch: {group.get('group_id')}")
        for shift_field in {value["shift_field"] for value in VARIANTS.values()}:
            if any(nuclei[index].get(shift_field) is None for index in members):
                continue
            raw = [float(nuclei[index][shift_field]) for index in members]
            mean = float(sum(raw) / len(raw))
            for index in members:
                nuclei[index][shift_field] = mean
            audit.append(
                {
                    "notice": NOTICE,
                    "figure_id": figure_id,
                    "operation": "solution_equivalence_shift_average",
                    "group_id": group["group_id"],
                    "element": element,
                    "atom_indices_0based": ";".join(map(str, members)),
                    "atom_indices_1based": ";".join(str(value + 1) for value in members),
                    "shift_field": shift_field,
                    "raw_values": ";".join(f"{value:.12g}" for value in raw),
                    "averaged_value": mean,
                    "symmetry_source": group.get("symmetry_source", ""),
                    "config_sha256": sha256_file(config_path),
                    "raw_source_retained": True,
                }
            )

    # Symmetrize only complete J orbits under the declared molecular operation.
    j_by_pair = {
        tuple(
            sorted(
                (
                    int(row["nucleus_A_index_0based"]),
                    int(row["nucleus_B_index_0based"]),
                )
            )
        ): row
        for row in bundle["averaged_JHH"]
    }
    visited: set[tuple[int, int]] = set()
    for raw_permutation in compound_rules.get("proton_symmetry_permutations_0based", []):
        permutation = {int(key): int(value) for key, value in raw_permutation.items()}
        for pair in sorted(j_by_pair):
            if pair in visited:
                continue
            transformed = tuple(sorted((permutation.get(pair[0], pair[0]), permutation.get(pair[1], pair[1]))))
            orbit = sorted({pair, transformed})
            if len(orbit) == 1:
                visited.add(pair)
                continue
            if any(item not in j_by_pair for item in orbit):
                raise SpectrumError(
                    f"Incomplete declared J symmetry orbit for {figure_id}: {orbit}"
                )
            raw = [float(j_by_pair[item]["boltzmann_averaged_J_Hz"]) for item in orbit]
            mean = float(sum(raw) / len(raw))
            for item in orbit:
                j_by_pair[item]["boltzmann_averaged_J_Hz"] = mean
                visited.add(item)
            audit.append(
                {
                    "notice": NOTICE,
                    "figure_id": figure_id,
                    "operation": "solution_equivalence_J_orbit_average",
                    "group_id": "DECLARED_PROTON_SYMMETRY_OPERATION",
                    "element": "H-H",
                    "atom_indices_0based": ";".join(f"{a}-{b}" for a, b in orbit),
                    "atom_indices_1based": ";".join(f"{a + 1}-{b + 1}" for a, b in orbit),
                    "shift_field": "not_applicable",
                    "raw_values": ";".join(f"{value:.12g}" for value in raw),
                    "averaged_value": mean,
                    "symmetry_source": compound_rules.get("symmetry_source", ""),
                    "config_sha256": sha256_file(config_path),
                    "raw_source_retained": True,
                }
            )
    return bundle, audit


def discover_bundles(input_root: Path, only: set[str]) -> list[Path]:
    candidates = sorted(input_root.rglob("F3_*_prediction_bundle.json"))
    selected: list[Path] = []
    identities: dict[tuple[str, str], Path] = {}
    found: set[str] = set()
    for path in candidates:
        bundle = load_bundle(path)
        compound = bundle["compound"]
        figure_id = str(compound.get("figure_id", ""))
        scenario_id = str(compound.get("scenario_id", ""))
        if not figure_id or not scenario_id:
            raise SpectrumError(f"Prediction bundle lacks figure/scenario identity: {path}")
        if only and figure_id not in only:
            continue
        key = (figure_id, scenario_id)
        if key in identities:
            raise SpectrumError(
                f"Duplicate prediction bundles for {figure_id}/{scenario_id}: "
                f"{identities[key]} and {path}"
            )
        identities[key] = path
        found.add(figure_id)
        selected.append(path)
    if only:
        missing = only - found
        if missing:
            raise SpectrumError(f"Prediction bundles not found for: {sorted(missing)}")
    if not selected:
        raise SpectrumError(f"No prediction bundles found under {input_root}")
    return selected


def appearance_line_rows(
    sticks: Sequence[dict[str, Any]],
    figure_id: str,
    scenario_id: str,
    nucleus: str,
    variant_key: str,
) -> list[dict[str, Any]]:
    """Export already-resolved calculated transitions for the appearance layer.

    No linewidth is copied: the acquisition profile owns broadening, including
    exchangeable-H broadening.  Multiplicity is fixed to singlet because J
    splitting is already represented by the individual transition rows.
    """

    result: list[dict[str, Any]] = []
    for index, row in enumerate(sticks, start=1):
        assignment = str(
            row.get("assignment_atom_label_1based")
            or row.get("component_atom_labels_1based")
            or row.get("atom_label_1based")
            or row.get("component_id")
            or f"transition_{index}"
        )
        result.append(
            {
                "figure_id": figure_id,
                "scenario_id": scenario_id,
                "nucleus": nucleus,
                "shift_variant": variant_key,
                "representation": "resolved_transitions",
                "ppm": float(row["shift_ppm"]),
                "intensity": float(row["intensity"]),
                "linewidth_hz": "",
                "exchangeable_H_flag": int(row.get("is_exchangeable_H", 0)) if nucleus == "1H" else 0,
                "multiplicity": "s",
                "J_hz": "",
                "assignment": assignment,
                "source_transition_id": str(row.get("component_id") or f"{nucleus}_{index:05d}"),
                "source_model": str(row.get("model", "")),
                "notice": NOTICE,
            }
        )
    return result


def _split_semicolon_values(value: Any) -> list[str]:
    return [item.strip() for item in str(value or "").split(";") if item.strip()]


def exchangeable_proton_map(bundle: dict[str, Any]) -> dict[int, dict[str, Any]]:
    """Return structurally exchangeable H atoms, validating bundle connectivity.

    New bundles contain explicit parent fields.  For an older audited bundle, the
    parent is deterministically recovered from its atom-map connectivity.  A flag
    that conflicts with H-to-N/O/S connectivity is rejected instead of trusted.
    """
    atom_rows = bundle.get("atom_map")
    if not isinstance(atom_rows, list) or not atom_rows:
        raise SpectrumError("Prediction bundle lacks a non-empty atom_map")
    by_index = {int(row["nucleus_index_0based"]): row for row in atom_rows}
    if len(by_index) != len(atom_rows):
        raise SpectrumError("Prediction bundle atom_map contains duplicate indices")
    result: dict[int, dict[str, Any]] = {}
    for index, row in sorted(by_index.items()):
        if str(row.get("element")) != "H":
            continue
        attached_indices = [
            int(value)
            for value in _split_semicolon_values(row.get("attached_atom_indices_0based"))
        ]
        attached_elements = _split_semicolon_values(row.get("attached_elements"))
        if len(attached_indices) != 1 or len(attached_elements) != 1:
            raise SpectrumError(
                f"Hydrogen atom {index} lacks one unambiguous bonded parent in atom_map"
            )
        parent_index = attached_indices[0]
        if parent_index not in by_index:
            raise SpectrumError(
                f"Hydrogen atom {index} references absent parent {parent_index}"
            )
        parent_element = str(by_index[parent_index].get("element"))
        if attached_elements[0] != parent_element:
            raise SpectrumError(
                f"Hydrogen atom {index} parent-element mismatch: "
                f"{attached_elements[0]} vs {parent_element}"
            )
        expected = int(parent_element in EXCHANGEABLE_H_PARENT_ELEMENTS)
        declared = int(row.get("exchangeable_H_flag", expected))
        if declared != expected:
            raise SpectrumError(
                f"Hydrogen atom {index} exchangeable flag conflicts with parent "
                f"{parent_element}{parent_index + 1}"
            )
        explicit_parent = row.get("exchangeable_H_parent_index_0based")
        if explicit_parent not in (None, "") and int(explicit_parent) != parent_index:
            raise SpectrumError(
                f"Hydrogen atom {index} explicit exchangeable parent conflicts with connectivity"
            )
        if expected:
            result[index] = {
                "nucleus_index_0based": index,
                "atom_label_1based": str(row.get("atom_label_1based", f"H{index + 1}")),
                "parent_index_0based": parent_index,
                "parent_index_1based": parent_index + 1,
                "parent_element": parent_element,
                "exchangeable_H_class": f"{parent_element}-H",
                "assignment_method": (
                    "atom_map_connectivity_H_bound_to_N_O_or_S;flag_cross_checked"
                ),
                "structural_confidence": "high",
                "solution_observation_status": "provisional_not_experimentally_validated",
            }
    averaged_h = {
        int(row["nucleus_index_0based"]): int(row.get("exchangeable_H_flag", 0))
        for row in bundle["averaged_nuclei"]
        if str(row.get("element")) == "H"
    }
    for index in averaged_h:
        if bool(averaged_h[index]) != (index in result):
            raise SpectrumError(
                f"Averaged-nucleus exchangeable flag conflicts with atom_map at H{index + 1}"
            )
    return result


def spin_operators(n_spins: int) -> tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray]]:
    identity = np.eye(2, dtype=complex)
    ix_single = 0.5 * np.array([[0, 1], [1, 0]], dtype=complex)
    iy_single = 0.5 * np.array([[0, -1j], [1j, 0]], dtype=complex)
    iz_single = 0.5 * np.array([[1, 0], [0, -1]], dtype=complex)

    def embedded(single: np.ndarray, position: int) -> np.ndarray:
        result = np.array([[1.0 + 0.0j]])
        for index in range(n_spins):
            result = np.kron(result, single if index == position else identity)
        return result

    return (
        [embedded(ix_single, index) for index in range(n_spins)],
        [embedded(iy_single, index) for index in range(n_spins)],
        [embedded(iz_single, index) for index in range(n_spins)],
    )


def exact_spin_sticks(
    shifts_ppm: np.ndarray,
    j_hz: np.ndarray,
    frequency_mhz: float,
    intensity_cutoff: float,
) -> list[dict[str, float]]:
    n_spins = len(shifts_ppm)
    if n_spins == 0:
        return []
    ix, iy, iz = spin_operators(n_spins)
    hamiltonian = np.zeros_like(ix[0])
    for index, shift in enumerate(shifts_ppm):
        hamiltonian += float(shift) * frequency_mhz * iz[index]
    for first in range(n_spins):
        for second in range(first + 1, n_spins):
            coupling = float(j_hz[first, second])
            if coupling:
                hamiltonian += coupling * (
                    ix[first] @ ix[second]
                    + iy[first] @ iy[second]
                    + iz[first] @ iz[second]
                )
    energies, vectors = np.linalg.eigh(hamiltonian)
    mx = sum(ix)
    mx_eigen = vectors.conj().T @ mx @ vectors
    rows: list[dict[str, float]] = []
    for low in range(len(energies)):
        for high in range(low + 1, len(energies)):
            intensity = float(abs(mx_eigen[low, high]) ** 2)
            if intensity < intensity_cutoff:
                continue
            frequency_hz = float(energies[high] - energies[low])
            if frequency_hz <= 0.0:
                continue
            rows.append(
                {
                    "shift_ppm": frequency_hz / frequency_mhz,
                    "frequency_Hz": frequency_hz,
                    "intensity": intensity,
                }
            )
    total = sum(row["intensity"] for row in rows)
    if total <= 0.0:
        raise SpectrumError("Exact spin Hamiltonian produced no observable transitions")
    for row in rows:
        row["intensity"] *= n_spins / total
    return rows


def merge_sticks(
    rows: Sequence[dict[str, Any]], tolerance_hz: float, frequency_mhz: float
) -> list[dict[str, Any]]:
    merged: list[dict[str, Any]] = []
    for row in sorted(rows, key=lambda item: float(item["shift_ppm"])):
        if (
            merged
            and row.get("component_id") == merged[-1].get("component_id")
            and row.get("model") == merged[-1].get("model")
            and abs(float(row["frequency_Hz"]) - float(merged[-1]["frequency_Hz"]))
            <= tolerance_hz
        ):
            old_intensity = float(merged[-1]["intensity"])
            new_intensity = old_intensity + float(row["intensity"])
            merged[-1]["shift_ppm"] = (
                float(merged[-1]["shift_ppm"]) * old_intensity
                + float(row["shift_ppm"]) * float(row["intensity"])
            ) / new_intensity
            merged[-1]["frequency_Hz"] = merged[-1]["shift_ppm"] * frequency_mhz
            merged[-1]["intensity"] = new_intensity
            merged[-1]["merged_transition_count"] = int(
                merged[-1].get("merged_transition_count", 1)
            ) + int(row.get("merged_transition_count", 1))
        else:
            copied = dict(row)
            copied.setdefault("merged_transition_count", 1)
            merged.append(copied)
    return merged


def first_order_component_sticks(
    atom_indices: Sequence[int],
    shifts: dict[int, float],
    j_lookup: dict[tuple[int, int], float],
    frequency_mhz: float,
    edge_threshold_hz: float,
    maximum_neighbors: int,
) -> tuple[list[dict[str, Any]], dict[int, list[int]]]:
    rows: list[dict[str, Any]] = []
    omitted_by_atom: dict[int, list[int]] = {}
    component_set = set(atom_indices)
    for atom_index in atom_indices:
        neighbors = [
            other
            for other in component_set
            if other != atom_index
            and abs(j_lookup.get(tuple(sorted((atom_index, other))), 0.0))
            >= edge_threshold_hz
        ]
        neighbors.sort(
            key=lambda other: abs(j_lookup[tuple(sorted((atom_index, other)))]),
            reverse=True,
        )
        used = neighbors[:maximum_neighbors]
        omitted = neighbors[maximum_neighbors:]
        omitted_by_atom[atom_index] = omitted
        lines = [(0.0, 1.0)]
        for neighbor in used:
            coupling = j_lookup[tuple(sorted((atom_index, neighbor)))]
            split = coupling / (2.0 * frequency_mhz)
            next_lines: list[tuple[float, float]] = []
            for offset, intensity in lines:
                next_lines.append((offset - split, 0.5 * intensity))
                next_lines.append((offset + split, 0.5 * intensity))
            lines = next_lines
        for offset, intensity in lines:
            shift = shifts[atom_index] + offset
            rows.append(
                {
                    "shift_ppm": shift,
                    "frequency_Hz": shift * frequency_mhz,
                    "intensity": intensity,
                    "assignment_atom_index_0based": atom_index,
                    "first_order_neighbors_used": ";".join(str(value) for value in used),
                    "first_order_neighbors_omitted": ";".join(
                        str(value) for value in omitted
                    ),
                }
            )
    return rows, omitted_by_atom


def connected_components(nodes: Sequence[int], edges: Sequence[tuple[int, int]]) -> list[list[int]]:
    adjacency: dict[int, set[int]] = {node: set() for node in nodes}
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    unseen = set(nodes)
    components: list[list[int]] = []
    while unseen:
        start = min(unseen)
        queue = deque([start])
        unseen.remove(start)
        component: list[int] = []
        while queue:
            node = queue.popleft()
            component.append(node)
            for neighbor in sorted(adjacency[node]):
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    queue.append(neighbor)
        components.append(sorted(component))
    return components


def collect_heteronuclear_couplings(
    bundle: dict[str, Any],
    observed_element: str,
    decouple_f: bool,
    decouple_p: bool,
    threshold_hz: float,
) -> tuple[dict[int, list[dict[str, Any]]], list[dict[str, Any]], list[str]]:
    """Collect H/C-to-F/P pairs without requiring an F/P chemical shift."""
    exchangeable_indices = (
        set(exchangeable_proton_map(bundle)) if observed_element == "H" else set()
    )
    rows = bundle.get("averaged_scalar_couplings", bundle["averaged_JHH"])
    active: dict[int, list[dict[str, Any]]] = {}
    audit: list[dict[str, Any]] = []
    warnings: list[str] = []
    if "averaged_scalar_couplings" not in bundle:
        structure_elements = {
            str(row.get("element", ""))
            for row in bundle.get("atom_map", bundle["averaged_nuclei"])
        }
        if structure_elements & {"F", "P"}:
            warnings.append(
                "The bundle contains F/P but lacks averaged_scalar_couplings; no heteronuclear "
                "splitting could be applied from the legacy averaged_JHH-only schema"
            )
    seen: set[tuple[int, str, int]] = set()
    observed_shift_indices = {
        int(row["nucleus_index_0based"])
        for row in bundle["averaged_nuclei"]
        if row.get("element") == observed_element
    }
    for row in rows:
        element_a = str(row.get("nucleus_A_element", ""))
        element_b = str(row.get("nucleus_B_element", ""))
        index_a = int(row["nucleus_A_index_0based"])
        index_b = int(row["nucleus_B_index_0based"])
        if element_a == observed_element and element_b in {"F", "P"}:
            observed_index, hetero_index, hetero_element = index_a, index_b, element_b
        elif element_b == observed_element and element_a in {"F", "P"}:
            observed_index, hetero_index, hetero_element = index_b, index_a, element_a
        else:
            continue
        key = (observed_index, hetero_element, hetero_index)
        if key in seen:
            raise SpectrumError(
                f"Duplicate averaged {observed_element}-{hetero_element} scalar coupling: {key}"
            )
        seen.add(key)
        coupling = float(row["boltzmann_averaged_J_Hz"])
        if not math.isfinite(coupling):
            raise SpectrumError(f"Non-finite scalar coupling for {key}: {coupling}")
        decoupled = (hetero_element == "F" and decouple_f) or (
            hetero_element == "P" and decouple_p
        )
        touches_exchangeable_h = observed_index in exchangeable_indices
        below_threshold = abs(coupling) < threshold_hz
        missing_observed_shift = observed_index not in observed_shift_indices
        if missing_observed_shift:
            disposition = "omitted_observed_nucleus_has_no_shift"
            warnings.append(
                f"{observed_element}{observed_index + 1}-{hetero_element}{hetero_index + 1} "
                f"was not applied because {observed_element}{observed_index + 1} has no shift"
            )
        elif touches_exchangeable_h:
            disposition = "omitted_exchangeable_proton_default_observation_model"
        elif decoupled:
            disposition = f"omitted_by_{hetero_element.lower()}_decoupling"
        elif below_threshold:
            disposition = "omitted_below_coupling_threshold"
            warnings.append(
                f"{observed_element}{observed_index + 1}-{hetero_element}{hetero_index + 1}: "
                f"J={coupling:.6g} Hz was below the {threshold_hz:.6g} Hz display threshold"
            )
        else:
            disposition = "applied_high_field_secular_weak_heteronuclear_splitting"
        record = {
            "notice": NOTICE,
            "observed_element": observed_element,
            "observed_index_0based": observed_index,
            "heteronucleus_element": hetero_element,
            "heteronucleus_index_0based": hetero_index,
            "pair_label_1based": (
                f"{observed_element}{observed_index + 1}-"
                f"{hetero_element}{hetero_index + 1}"
            ),
            "scalar_coupling_pair_class": f"{observed_element}-{hetero_element}",
            "J_Hz": coupling,
            "coupling_threshold_Hz": threshold_hz,
            "decoupling_requested": int(decoupled),
            "touches_exchangeable_H": int(touches_exchangeable_h),
            "disposition": disposition,
            "heteronuclear_spin_quantum_number": 0.5,
            "fp_chemical_shift_required_for_this_model": False,
            "heteronucleus_shift_reference_status": (
                bundle.get("reference_shieldings", {})
                .get(hetero_element, {})
                .get("status", "not_reported")
            ),
            "model_limitation": (
                "High-field secular/weak-heteronuclear treatment: heteronuclear flip-flop, "
                "second-order, relaxation and line-shape effects are omitted"
            ),
        }
        audit.append(record)
        if (
            not missing_observed_shift
            and not touches_exchangeable_h
            and not decoupled
            and not below_threshold
        ):
            active.setdefault(observed_index, []).append(record)
    for observed_index in active:
        active[observed_index].sort(
            key=lambda item: (
                str(item["heteronucleus_element"]),
                int(item["heteronucleus_index_0based"]),
            )
        )
    return active, audit, warnings


def heteronuclear_shift_states(
    observed_indices: Sequence[int],
    base_shifts: dict[int, float],
    coupling_map: dict[int, list[dict[str, Any]]],
    frequency_mhz: float,
) -> tuple[list[tuple[str, int]], list[dict[int, float]]]:
    """Return all equally populated F/P m=±1/2 secular shift states."""
    heteronuclei = sorted(
        {
            (
                str(record["heteronucleus_element"]),
                int(record["heteronucleus_index_0based"]),
            )
            for observed_index in observed_indices
            for record in coupling_map.get(observed_index, [])
        }
    )
    if not heteronuclei:
        return [], [dict(base_shifts)]
    states: list[dict[int, float]] = []
    for projections in product((-0.5, 0.5), repeat=len(heteronuclei)):
        projection_by_nucleus = dict(zip(heteronuclei, projections))
        shifted = dict(base_shifts)
        for observed_index in observed_indices:
            offset_hz = sum(
                float(record["J_Hz"])
                * projection_by_nucleus[
                    (
                        str(record["heteronucleus_element"]),
                        int(record["heteronucleus_index_0based"]),
                    )
                ]
                for record in coupling_map.get(observed_index, [])
            )
            shifted[observed_index] += offset_hz / frequency_mhz
        states.append(shifted)
    return heteronuclei, states


def proton_dynamic_topology(
    bundle: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Identify only strict methyl rotors; flag CH2/NH2 without merging them."""
    atom_rows = bundle.get("atom_map")
    if not isinstance(atom_rows, list) or not atom_rows:
        raise SpectrumError("Prediction bundle lacks atom_map for dynamic-equivalence audit")
    by_index = {int(row["nucleus_index_0based"]): row for row in atom_rows}
    methyl_groups: list[dict[str, Any]] = []
    unmerged: list[dict[str, Any]] = []
    closed_shell = int(bundle.get("compound", {}).get("multiplicity", 1)) == 1
    for parent_index, parent in sorted(by_index.items()):
        attached = [
            int(value)
            for value in _split_semicolon_values(parent.get("attached_atom_indices_0based"))
        ]
        hydrogen_indices = sorted(
            index for index in attached if str(by_index[index].get("element")) == "H"
        )
        parent_element = str(parent.get("element"))
        if (
            parent_element == "C"
            and len(attached) == 4
            and len(hydrogen_indices) == 3
            and int(parent.get("formal_charge", 0)) == 0
            and closed_shell
        ):
            heavy_neighbors = [index for index in attached if index not in hydrogen_indices]
            if len(heavy_neighbors) != 1:
                raise SpectrumError(
                    f"Methyl topology at C{parent_index + 1} lacks exactly one heavy neighbor"
                )
            for hydrogen_index in hydrogen_indices:
                h_attached = [
                    int(value)
                    for value in _split_semicolon_values(
                        by_index[hydrogen_index].get("attached_atom_indices_0based")
                    )
                ]
                if h_attached != [parent_index]:
                    raise SpectrumError(
                        f"Methyl H{hydrogen_index + 1} connectivity is not reciprocal"
                    )
            methyl_groups.append(
                {
                    "group_id": f"METHYL_C{parent_index + 1}",
                    "parent_index_0based": parent_index,
                    "parent_label_1based": f"C{parent_index + 1}",
                    "heavy_neighbor_index_0based": heavy_neighbors[0],
                    "hydrogen_indices_0based": hydrogen_indices,
                    "assignment_rule": (
                        "neutral closed-shell carbon with exactly three H and one heavy-atom neighbor"
                    ),
                }
            )
        elif parent_element in {"C", "N"} and len(hydrogen_indices) == 2:
            unmerged.append(
                {
                    "group_id": f"STATIC_{parent_element}{parent_index + 1}_H2",
                    "parent_index_0based": parent_index,
                    "parent_label_1based": f"{parent_element}{parent_index + 1}",
                    "hydrogen_indices_0based": hydrogen_indices,
                    "status": "not_merged",
                    "reason": (
                        "CH2 protons may be enantiotopic/diastereotopic and are not automatically equivalent"
                        if parent_element == "C"
                        else "NH2 protons may be exchange-sensitive or inequivalent under restricted rotation"
                    ),
                }
            )
    return methyl_groups, unmerged


def apply_methyl_rotation_model(
    shifts: dict[int, float],
    j_lookup: dict[tuple[int, int], float],
    heteronuclear_couplings: dict[int, list[dict[str, Any]]],
    methyl_groups: Sequence[dict[str, Any]],
    unmerged_groups: Sequence[dict[str, Any]],
    all_proton_indices: Sequence[int],
    variant_key: str,
) -> tuple[
    dict[int, float],
    dict[tuple[int, int], float],
    dict[int, list[dict[str, Any]]],
    list[dict[str, Any]],
]:
    """Apply a conservative fast-methyl-rotation observation model.

    Three methyl shifts are averaged. Internal methyl H-H J values are omitted as
    equivalent-spin couplings. External H-H or H-F/P J values are symmetrized only
    when all three member values exist; incomplete triplets are excluded from the
    display model but remain untouched in the prediction bundle/audit.
    """
    averaged_shifts = dict(shifts)
    display_j = dict(j_lookup)
    display_hetero = {
        index: [dict(record) for record in records]
        for index, records in heteronuclear_couplings.items()
    }
    audit: list[dict[str, Any]] = []
    all_indices = set(all_proton_indices)
    for group in methyl_groups:
        members = list(group["hydrogen_indices_0based"])
        member_set = set(members)
        original_shifts = [float(shifts[index]) for index in members]
        mean_shift = float(sum(original_shifts) / 3.0)
        for index in members:
            averaged_shifts[index] = mean_shift

        internal_removed: list[str] = []
        for first_position, first in enumerate(members):
            for second in members[first_position + 1 :]:
                key = tuple(sorted((first, second)))
                if key in display_j:
                    internal_removed.append(f"H{first + 1}-H{second + 1}")
                    display_j.pop(key)

        complete_external_hh = 0
        incomplete_external_hh = 0
        external_hh_actions: list[str] = []
        for outside in sorted(all_indices - member_set):
            keys = [tuple(sorted((member, outside))) for member in members]
            present = [key for key in keys if key in display_j]
            if len(present) == 3:
                mean_j = float(sum(display_j[key] for key in present) / 3.0)
                for key in keys:
                    display_j[key] = mean_j
                complete_external_hh += 1
                external_hh_actions.append(
                    f"H{outside + 1}:complete_triplet_symmetrized_to_{mean_j:.12g}_Hz"
                )
            elif present:
                for key in present:
                    display_j.pop(key)
                incomplete_external_hh += 1
                external_hh_actions.append(
                    f"H{outside + 1}:incomplete_{len(present)}_of_3_excluded_from_display"
                )

        heteronucleus_keys = sorted(
            {
                (
                    str(record["heteronucleus_element"]),
                    int(record["heteronucleus_index_0based"]),
                )
                for member in members
                for record in display_hetero.get(member, [])
            }
        )
        complete_hetero = 0
        incomplete_hetero = 0
        heteronuclear_actions: list[str] = []
        for hetero_key in heteronucleus_keys:
            found: list[tuple[int, dict[str, Any]]] = []
            for member in members:
                for record in display_hetero.get(member, []):
                    key = (
                        str(record["heteronucleus_element"]),
                        int(record["heteronucleus_index_0based"]),
                    )
                    if key == hetero_key:
                        found.append((member, record))
            if len(found) == 3:
                mean_j = float(sum(float(record["J_Hz"]) for _, record in found) / 3.0)
                for _, record in found:
                    record["J_Hz"] = mean_j
                    record["methyl_rotation_symmetrized"] = 1
                complete_hetero += 1
                heteronuclear_actions.append(
                    f"{hetero_key[0]}{hetero_key[1] + 1}:complete_triplet_"
                    f"symmetrized_to_{mean_j:.12g}_Hz"
                )
            elif found:
                for member, record in found:
                    display_hetero[member].remove(record)
                incomplete_hetero += 1
                heteronuclear_actions.append(
                    f"{hetero_key[0]}{hetero_key[1] + 1}:incomplete_{len(found)}_"
                    "of_3_excluded_from_display"
                )

        audit.append(
            {
                "notice": NOTICE,
                "variant": variant_key,
                **group,
                "status": "applied_topology_based_fast_methyl_rotation_model",
                "hydrogen_labels_1based": ";".join(f"H{index + 1}" for index in members),
                "pre_average_shifts_ppm": ";".join(
                    f"{value:.12g}" for value in original_shifts
                ),
                "post_average_shift_ppm_each": mean_shift,
                "sum_shifts_before_ppm": sum(original_shifts),
                "sum_shifts_after_ppm": 3.0 * mean_shift,
                "nominal_integral_before_H": 3.0,
                "nominal_integral_after_H": 3.0,
                "internal_HH_pairs_removed_from_display": ";".join(internal_removed),
                "complete_external_HH_triplets_symmetrized": complete_external_hh,
                "incomplete_external_HH_triplets_excluded": incomplete_external_hh,
                "external_HH_triplet_actions": ";".join(external_hh_actions),
                "complete_HF_HP_triplets_symmetrized": complete_hetero,
                "incomplete_HF_HP_triplets_excluded": incomplete_hetero,
                "HF_HP_triplet_actions": ";".join(heteronuclear_actions),
                "assumption": (
                    "methyl rotation is fast on the solution-NMR timescale; this may fail "
                    "at very low temperature or for an unusually high rotational barrier"
                ),
            }
        )
    for group in unmerged_groups:
        audit.append(
            {
                "notice": NOTICE,
                "variant": variant_key,
                **group,
                "hydrogen_labels_1based": ";".join(
                    f"H{index + 1}" for index in group["hydrogen_indices_0based"]
                ),
                "nominal_integral_before_H": 2.0,
                "nominal_integral_after_H": 2.0,
            }
        )
    return averaged_shifts, display_j, display_hetero, audit


def build_proton_sticks(
    bundle: dict[str, Any],
    variant_key: str,
    heteronuclear_couplings: dict[int, list[dict[str, Any]]],
    frequency_mhz: float,
    edge_threshold_hz: float,
    exact_max_spins: int,
    first_order_max_neighbors: int,
    merge_tolerance_hz: float,
    transition_intensity_cutoff: float,
    regular_linewidth_hz: float,
    exchangeable_linewidth_hz: float,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[str],
]:
    variant = VARIANTS[variant_key]
    proton_rows = [row for row in bundle["averaged_nuclei"] if row["element"] == "H"]
    if not proton_rows:
        raise SpectrumError("Prediction bundle contains no calculated proton shifts")
    indices = sorted(int(row["nucleus_index_0based"]) for row in proton_rows)
    if len(indices) != len(set(indices)):
        raise SpectrumError("Prediction bundle contains duplicate proton indices")
    by_index = {int(row["nucleus_index_0based"]): row for row in proton_rows}
    shifts = {index: float(by_index[index][variant["shift_field"]]) for index in indices}
    labels = {index: str(by_index[index]["atom_label_1based"]) for index in indices}
    exchangeable = exchangeable_proton_map(bundle)
    exchangeable_indices = set(exchangeable)
    nonexchangeable_indices = [index for index in indices if index not in exchangeable_indices]

    j_lookup: dict[tuple[int, int], float] = {}
    ignored_pairs: list[str] = []
    exchangeable_coupling_audit: list[dict[str, Any]] = []
    for row in bundle["averaged_JHH"]:
        first = int(row["nucleus_A_index_0based"])
        second = int(row["nucleus_B_index_0based"])
        if first not in by_index or second not in by_index:
            ignored_pairs.append(f"{first}-{second}: pair references a nucleus without a proton shift")
            continue
        key = tuple(sorted((first, second)))
        if key in j_lookup:
            raise SpectrumError(f"Duplicate averaged J pair: {key}")
        j_lookup[key] = float(row["boltzmann_averaged_J_Hz"])

    methyl_groups, unmerged_dynamic_groups = proton_dynamic_topology(bundle)
    (
        shifts,
        display_j_lookup,
        display_heteronuclear_couplings,
        dynamic_equivalence_audit,
    ) = apply_methyl_rotation_model(
        shifts,
        j_lookup,
        heteronuclear_couplings,
        methyl_groups,
        unmerged_dynamic_groups,
        indices,
        variant_key,
    )

    nonexchangeable_j_lookup: dict[tuple[int, int], float] = {}
    for (first, second), coupling in sorted(j_lookup.items()):
        touches_exchangeable = first in exchangeable_indices or second in exchangeable_indices
        excluded_by_dynamic_model = (
            not touches_exchangeable and (first, second) not in display_j_lookup
        )
        record = {
            "notice": NOTICE,
            "variant": variant_key,
            "nucleus_A_index_0based": first,
            "nucleus_A_label_1based": labels[first],
            "nucleus_B_index_0based": second,
            "nucleus_B_label_1based": labels[second],
            "pair_label_1based": f"{labels[first]}-{labels[second]}",
            "J_Hz": coupling,
            "J_Hz_used_in_display": (
                display_j_lookup.get((first, second))
                if not touches_exchangeable
                else None
            ),
            "touches_exchangeable_H": int(touches_exchangeable),
            "disposition": (
                "retained_in_prediction_bundle_but_excluded_from_coherent_fine_structure"
                if touches_exchangeable
                else (
                    "retained_in_prediction_bundle_but_excluded_by_dynamic_equivalence_model"
                    if excluded_by_dynamic_model
                    else "eligible_for_nonexchangeable_HH_model_subject_to_threshold"
                )
            ),
            "reason": (
                "rapid/condition-dependent proton exchange and relaxation are not represented "
                "by a static isotropic spin Hamiltonian"
                if touches_exchangeable
                else (
                    "methyl equivalent-spin/internal or incomplete-triplet coupling handling"
                    if excluded_by_dynamic_model
                    else "neither endpoint is structurally exchangeable"
                )
            ),
        }
        if touches_exchangeable:
            exchangeable_coupling_audit.append(record)
        if not touches_exchangeable and not excluded_by_dynamic_model:
            nonexchangeable_j_lookup[(first, second)] = float(
                display_j_lookup[(first, second)]
            )

    active_edges = [
        pair
        for pair, value in sorted(nonexchangeable_j_lookup.items())
        if abs(value) >= edge_threshold_hz
    ]
    components = connected_components(nonexchangeable_indices, active_edges)
    sticks: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    for component_number, component in enumerate(components, start=1):
        component_id = f"HCOMP_{component_number:02d}"
        component_labels = [labels[index] for index in component]
        component_methyl_groups = [
            str(group["group_id"])
            for group in methyl_groups
            if set(group["hydrogen_indices_0based"]) & set(component)
        ]
        contains_exchangeable = False
        dimension = 2 ** len(component)
        heteronuclei, heteronuclear_shifted_states = heteronuclear_shift_states(
            component,
            shifts,
            display_heteronuclear_couplings,
            frequency_mhz,
        )
        heteronuclear_records = [
            record
            for atom_index in component
            for record in display_heteronuclear_couplings.get(atom_index, [])
        ]
        has_heteronuclear_splitting = bool(heteronuclei)
        heteronuclear_state_count = len(heteronuclear_shifted_states)
        if len(component) <= exact_max_spins:
            j_matrix = np.zeros((len(component), len(component)), dtype=float)
            local = {atom_index: index for index, atom_index in enumerate(component)}
            retained_edges = 0
            for (first, second), value in nonexchangeable_j_lookup.items():
                if (
                    first in local
                    and second in local
                    and abs(value) >= edge_threshold_hz
                ):
                    j_matrix[local[first], local[second]] = value
                    j_matrix[local[second], local[first]] = value
                    retained_edges += 1
            if not has_heteronuclear_splitting:
                component_rows = exact_spin_sticks(
                    np.asarray([shifts[index] for index in component], dtype=float),
                    j_matrix,
                    frequency_mhz,
                    transition_intensity_cutoff,
                )
                model = "exact_isotropic_spin_hamiltonian"
                reason = (
                    f"component_size={len(component)} <= exact_max_spins={exact_max_spins}"
                )
                approximation_note = "none"
            else:
                component_rows = []
                for shifted in heteronuclear_shifted_states:
                    state_rows = exact_spin_sticks(
                        np.asarray([shifted[index] for index in component], dtype=float),
                        j_matrix,
                        frequency_mhz,
                        transition_intensity_cutoff,
                    )
                    for row in state_rows:
                        row["intensity"] = (
                            float(row["intensity"]) / heteronuclear_state_count
                        )
                    component_rows.extend(state_rows)
                model = "exact_HH_hamiltonian_plus_secular_weak_heteronuclear_splitting"
                reason = (
                    f"component_size={len(component)} <= exact_max_spins={exact_max_spins}; "
                    f"enumerated {heteronuclear_state_count} equally populated F/P m-state combinations"
                )
                approximation_note = (
                    "H-H state mixing was treated by the exact isotropic Hamiltonian within each "
                    "static F/P m-state. H-F/H-P couplings used the high-field secular weak-"
                    "heteronuclear limit; heteronuclear flip-flop, second-order, relaxation and "
                    "line-shape effects were omitted."
                )
            omitted_neighbors = ""
        else:
            if not has_heteronuclear_splitting:
                component_rows, omitted = first_order_component_sticks(
                    component,
                    shifts,
                    nonexchangeable_j_lookup,
                    frequency_mhz,
                    edge_threshold_hz,
                    first_order_max_neighbors,
                )
            else:
                component_rows = []
                omitted = {}
                for shifted in heteronuclear_shifted_states:
                    state_rows, state_omitted = first_order_component_sticks(
                        component,
                        shifted,
                        nonexchangeable_j_lookup,
                        frequency_mhz,
                        edge_threshold_hz,
                        first_order_max_neighbors,
                    )
                    omitted = state_omitted
                    for row in state_rows:
                        row["intensity"] = (
                            float(row["intensity"]) / heteronuclear_state_count
                        )
                    component_rows.extend(state_rows)
            retained_edges = sum(
                1
                for (first, second), value in nonexchangeable_j_lookup.items()
                if first in component
                and second in component
                and abs(value) >= edge_threshold_hz
            )
            model = (
                "first_order_HH_plus_secular_weak_heteronuclear_splitting"
                if has_heteronuclear_splitting
                else "first_order_weak_coupling_approximation"
            )
            reason = f"component_size={len(component)} > exact_max_spins={exact_max_spins}"
            if has_heteronuclear_splitting:
                reason += (
                    f"; enumerated {heteronuclear_state_count} equally populated F/P m-state combinations"
                )
            omitted_flat = sorted({value for values in omitted.values() for value in values})
            omitted_neighbors = ";".join(str(value) for value in omitted_flat)
            approximation_note = (
                "Each spin was split independently by retained |J| values; strong-coupling state mixing "
                "was not represented. Neighbor enumeration was capped at "
                f"{first_order_max_neighbors} per spin."
            )
            if has_heteronuclear_splitting:
                approximation_note += (
                    " H-F/H-P couplings used enumerated spin-1/2 m states in the high-field "
                    "secular weak-heteronuclear limit; flip-flop, second-order, relaxation and "
                    "line-shape effects were omitted."
                )
        for row in component_rows:
            row.update(
                {
                    "notice": NOTICE,
                    "variant": variant_key,
                    "shift_source": variant["label"],
                    "component_id": component_id,
                    "component_size": len(component),
                    "component_atom_indices_0based": ";".join(
                        str(value) for value in component
                    ),
                    "component_atom_labels_1based": ";".join(component_labels),
                    "model": model,
                    "model_reason": reason,
                    "contains_exchangeable_H": int(contains_exchangeable),
                    "is_exchangeable_H": 0,
                    "exchangeable_H_class": "not_exchangeable",
                    "shift_observation_status": "calculated_prediction_nonexchangeable",
                    "linewidth_observation_status": "display_parameter_not_predicted",
                    "integral_observation_status": "nominal_one_per_explicit_H_not_quantitative",
                    "display_linewidth_Hz": regular_linewidth_hz,
                    "dynamic_equivalence_model": (
                        "fast_methyl_rotation_symmetrization"
                        if component_methyl_groups
                        else "none"
                    ),
                    "dynamic_equivalence_group_ids": ";".join(
                        component_methyl_groups
                    ),
                }
            )
            if has_heteronuclear_splitting:
                row.update(
                    {
                        "heteronuclear_secular_splitting_applied": 1,
                        "heteronuclear_spin_state_count": heteronuclear_state_count,
                        "heteronuclear_nuclei": ";".join(
                            f"{element}{index + 1}" for element, index in heteronuclei
                        ),
                        "heteronuclear_pair_labels_1based": ";".join(
                            str(record["pair_label_1based"])
                            for record in heteronuclear_records
                        ),
                        "heteronuclear_J_Hz": ";".join(
                            f"{float(record['J_Hz']):.12g}"
                            for record in heteronuclear_records
                        ),
                    }
                )
        component_rows = merge_sticks(component_rows, merge_tolerance_hz, frequency_mhz)
        sticks.extend(component_rows)
        audit_row = {
            "notice": NOTICE,
            "variant": variant_key,
            "component_id": component_id,
            "component_size": len(component),
            "hilbert_dimension": dimension,
            "atom_indices_0based": ";".join(str(value) for value in component),
            "atom_labels_1based": ";".join(component_labels),
            "contains_exchangeable_H": int(contains_exchangeable),
            "retained_coupling_edges": retained_edges,
            "coupling_edge_threshold_Hz": edge_threshold_hz,
            "model": model,
            "model_reason": reason,
            "approximation_note": approximation_note,
            "omitted_neighbor_indices_due_to_cap": omitted_neighbors,
            "n_output_sticks_after_merge": len(component_rows),
            "merge_tolerance_Hz": merge_tolerance_hz,
            "dynamic_equivalence_model": (
                "fast_methyl_rotation_symmetrization"
                if component_methyl_groups
                else "none"
            ),
            "dynamic_equivalence_group_ids": ";".join(component_methyl_groups),
        }
        if has_heteronuclear_splitting:
            audit_row.update(
                {
                    "heteronuclear_model": "high_field_secular_weak_heteronuclear_m_state_enumeration",
                    "heteronuclear_spin_state_count": heteronuclear_state_count,
                    "HH_hilbert_dimension_times_heteronuclear_states": (
                        dimension * heteronuclear_state_count
                    ),
                    "heteronuclear_nuclei": ";".join(
                        f"{element}{index + 1}" for element, index in heteronuclei
                    ),
                    "heteronuclear_pair_labels_1based": ";".join(
                        str(record["pair_label_1based"])
                        for record in heteronuclear_records
                    ),
                    "fp_chemical_shift_required_for_splitting": False,
                }
            )
        audit.append(audit_row)

    # Exchangeable protons remain atom-resolved calculated shifts, but are not
    # allowed to create or receive coherent fine-structure in the default display.
    # Each is rendered as a separate provisional broad line with nominal unit area.
    for exchange_number, index in enumerate(sorted(exchangeable_indices), start=1):
        info = exchangeable[index]
        component_id = f"HEXCH_{exchange_number:02d}"
        shift = shifts[index]
        excluded_pairs = [
            row
            for row in exchangeable_coupling_audit
            if int(row["touches_exchangeable_H"]) == 1
            and index
            in {
                int(row["nucleus_A_index_0based"]),
                int(row["nucleus_B_index_0based"]),
            }
        ]
        sticks.append(
            {
                "notice": NOTICE,
                "variant": variant_key,
                "shift_source": variant["label"],
                "component_id": component_id,
                "component_size": 1,
                "component_atom_indices_0based": str(index),
                "component_atom_labels_1based": labels[index],
                "assignment_atom_index_0based": index,
                "assignment_atom_label_1based": labels[index],
                "shift_ppm": shift,
                "frequency_Hz": shift * frequency_mhz,
                "intensity": 1.0,
                "model": "provisional_exchangeable_broad_singlet",
                "model_reason": (
                    f"{info['exchangeable_H_class']} identified from explicit connectivity; "
                    "static calculated J values were retained for audit but excluded from "
                    "the default coherent fine-structure model"
                ),
                "contains_exchangeable_H": 1,
                "is_exchangeable_H": 1,
                "exchangeable_H_class": info["exchangeable_H_class"],
                "exchangeable_H_parent_index_0based": info["parent_index_0based"],
                "exchangeable_H_parent_index_1based": info["parent_index_1based"],
                "exchangeable_H_parent_element": info["parent_element"],
                "exchangeability_assignment_method": info["assignment_method"],
                "exchangeability_structural_confidence": info["structural_confidence"],
                "shift_observation_status": (
                    "provisional_calculated_shift_explicit_H_bonding_and_exchange_not_modelled"
                ),
                "linewidth_observation_status": (
                    "provisional_user_selected_display_broadening_not_predicted"
                ),
                "integral_observation_status": (
                    "provisional_nominal_one_explicit_H_not_quantitative"
                ),
                "display_linewidth_Hz": exchangeable_linewidth_hz,
                "nominal_integral_H": 1.0,
                "excluded_calculated_J_pair_count": len(excluded_pairs),
                "excluded_calculated_J_pair_labels_1based": ";".join(
                    str(row["pair_label_1based"]) for row in excluded_pairs
                ),
            }
        )
        audit.append(
            {
                "notice": NOTICE,
                "variant": variant_key,
                "component_id": component_id,
                "component_size": 1,
                "hilbert_dimension": 0,
                "atom_indices_0based": str(index),
                "atom_labels_1based": labels[index],
                "contains_exchangeable_H": 1,
                "retained_coupling_edges": 0,
                "excluded_calculated_J_pairs": len(excluded_pairs),
                "excluded_calculated_J_pair_labels_1based": ";".join(
                    str(row["pair_label_1based"]) for row in excluded_pairs
                ),
                "coupling_edge_threshold_Hz": edge_threshold_hz,
                "model": "provisional_exchangeable_broad_singlet",
                "model_reason": (
                    f"{info['exchangeable_H_class']}; coherent fine structure disabled by default"
                ),
                "approximation_note": (
                    "Shift, linewidth and nominal integral are provisional; no relaxation, "
                    "explicit hydrogen bonding, water content or proton-exchange kinetics were modelled."
                ),
                "omitted_neighbor_indices_due_to_cap": "",
                "n_output_sticks_after_merge": 1,
                "merge_tolerance_Hz": merge_tolerance_hz,
                "display_linewidth_Hz": exchangeable_linewidth_hz,
            }
        )
    sticks.sort(key=lambda row: float(row["shift_ppm"]))
    dropped_below_threshold = [
        f"{first}-{second}: J={value:.6g} Hz"
        for (first, second), value in sorted(nonexchangeable_j_lookup.items())
        if abs(value) < edge_threshold_hz
    ]
    warnings = ignored_pairs
    if dropped_below_threshold:
        warnings.append(
            "Couplings below the graph/model threshold were treated as zero for display: "
            + "; ".join(dropped_below_threshold)
        )
    warnings.append(
        "H-H pairs absent from the prediction bundle were treated as uncoupled for spectrum construction; "
        "the property workflow used a finite inter-nuclear distance threshold."
    )
    if exchangeable_indices:
        warnings.append(
            "Structurally exchangeable N-H/O-H/S-H protons were retained as atom-resolved "
            "provisional broad lines, while every calculated J pair touching them was excluded "
            "from the default coherent fine-structure model."
        )
    if methyl_groups:
        warnings.append(
            "Strict topological methyl groups used a fast-rotation shift/J symmetrization "
            "model; the original atom-resolved calculated values remain in the prediction bundle."
        )
    if unmerged_dynamic_groups:
        warnings.append(
            "CH2/NH2 hydrogens were deliberately not auto-merged; possible dynamic or "
            "stereochemical equivalence requires molecule- and condition-specific evidence."
        )
    return (
        sticks,
        audit,
        exchangeable_coupling_audit,
        dynamic_equivalence_audit,
        warnings,
    )


def build_carbon_sticks(
    bundle: dict[str, Any],
    variant_key: str,
    heteronuclear_couplings: dict[int, list[dict[str, Any]]],
    frequency_mhz: float,
    merge_tolerance_hz: float,
) -> list[dict[str, Any]]:
    variant = VARIANTS[variant_key]
    rows: list[dict[str, Any]] = []
    for nucleus in bundle["averaged_nuclei"]:
        if nucleus["element"] != "C":
            continue
        nucleus_index = int(nucleus["nucleus_index_0based"])
        shift = float(nucleus[variant["shift_field"]])
        heteronuclear_records = heteronuclear_couplings.get(nucleus_index, [])
        if not heteronuclear_records:
            rows.append(
                {
                    "notice": NOTICE,
                    "variant": variant_key,
                    "shift_source": variant["label"],
                    "nucleus_index_0based": nucleus_index,
                    "nucleus_index_1based": int(nucleus["nucleus_index_1based"]),
                    "atom_label_1based": nucleus["atom_label_1based"],
                    "shift_ppm": shift,
                    "intensity": 1.0,
                    "model": "proton_decoupled_13C_stick",
                    "model_reason": "one calculated stick per atom-resolved carbon; proton couplings omitted",
                }
            )
            continue
        lines = [(0.0, 1.0)]
        for record in heteronuclear_records:
            half_coupling_hz = float(record["J_Hz"]) / 2.0
            next_lines: list[tuple[float, float]] = []
            for offset_hz, intensity in lines:
                next_lines.append((offset_hz - half_coupling_hz, 0.5 * intensity))
                next_lines.append((offset_hz + half_coupling_hz, 0.5 * intensity))
            lines = next_lines
        carbon_rows: list[dict[str, Any]] = []
        for offset_hz, intensity in lines:
            split_shift = shift + offset_hz / frequency_mhz
            carbon_rows.append(
                {
                    "notice": NOTICE,
                    "variant": variant_key,
                    "shift_source": variant["label"],
                    "nucleus_index_0based": nucleus_index,
                    "nucleus_index_1based": int(nucleus["nucleus_index_1based"]),
                    "atom_label_1based": nucleus["atom_label_1based"],
                    "shift_ppm": split_shift,
                    "frequency_Hz": split_shift * frequency_mhz,
                    "intensity": intensity,
                    "component_id": f"C{nucleus_index + 1}",
                    "model": "proton_decoupled_13C_first_order_heteronuclear_split",
                    "model_reason": (
                        "13C{1H} stick split by retained C-F/C-P J values using equally populated "
                        "spin-1/2 m states in the high-field first-order limit"
                    ),
                    "heteronuclear_secular_splitting_applied": 1,
                    "heteronuclear_spin_state_count": len(lines),
                    "heteronuclear_nuclei": ";".join(
                        f"{record['heteronucleus_element']}"
                        f"{int(record['heteronucleus_index_0based']) + 1}"
                        for record in heteronuclear_records
                    ),
                    "heteronuclear_pair_labels_1based": ";".join(
                        str(record["pair_label_1based"])
                        for record in heteronuclear_records
                    ),
                    "heteronuclear_J_Hz": ";".join(
                        f"{float(record['J_Hz']):.12g}"
                        for record in heteronuclear_records
                    ),
                    "fp_chemical_shift_required_for_splitting": False,
                }
            )
        rows.extend(merge_sticks(carbon_rows, merge_tolerance_hz, frequency_mhz))
    return sorted(rows, key=lambda row: float(row["shift_ppm"]))


def lorentzian_profile(
    sticks: Sequence[dict[str, Any]],
    minimum_ppm: float,
    maximum_ppm: float,
    points: int,
    linewidth_hz: float,
    frequency_mhz: float,
) -> tuple[np.ndarray, np.ndarray]:
    if not sticks:
        raise SpectrumError("Cannot build a profile without sticks")
    if maximum_ppm <= minimum_ppm:
        raise SpectrumError("Spectrum range must have maximum > minimum")
    if points < 1000:
        raise SpectrumError("At least 1000 profile points are required")
    if linewidth_hz <= 0.0 or frequency_mhz <= 0.0:
        raise SpectrumError("Frequency and Lorentzian linewidth must be positive")
    x = np.linspace(minimum_ppm, maximum_ppm, points)
    y = np.zeros_like(x)
    for row in sticks:
        shift = float(row["shift_ppm"])
        intensity = float(row["intensity"])
        row_linewidth_hz = float(row.get("display_linewidth_Hz", linewidth_hz))
        if not math.isfinite(row_linewidth_hz) or row_linewidth_hz <= 0.0:
            raise SpectrumError(
                f"Invalid stick-specific display linewidth: {row_linewidth_hz}"
            )
        fwhm_ppm = row_linewidth_hz / frequency_mhz
        half_width = 0.5 * fwhm_ppm
        y += intensity * (half_width / math.pi) / ((x - shift) ** 2 + half_width**2)
    maximum = float(y.max())
    if maximum <= 0.0:
        raise SpectrumError("Lorentzian profile has no positive intensity")
    return x, y / maximum


def style_axis(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    sticks: Sequence[dict[str, Any]],
    color: str,
    nucleus_label: str,
    subtitle: str,
    frequency_mhz: float,
    linewidth_hz: float,
) -> None:
    ax.plot(x, y, color=color, linewidth=0.95, solid_joinstyle="round")
    maximum_stick = max(float(row["intensity"]) for row in sticks)
    for row in sticks:
        height = 0.085 * float(row["intensity"]) / maximum_stick
        is_exchangeable = int(row.get("is_exchangeable_H", 0)) == 1
        ax.vlines(
            float(row["shift_ppm"]),
            0.0,
            height,
            color="#9A5B3D" if is_exchangeable else "#59636E",
            linewidth=0.55 if is_exchangeable else 0.45,
            linestyles="dashed" if is_exchangeable else "solid",
        )
    ax.axhline(0.0, color="#9AA3AA", linewidth=0.45)
    ax.set_xlim(float(x.max()), float(x.min()))
    ax.set_ylim(-0.018, 1.07)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#4D555C")
    ax.tick_params(axis="x", direction="out", length=2.8, width=0.7, colors="#343A40")
    ax.set_xlabel(r"$\delta$ (ppm)", fontsize=7.5)
    ax.set_title(subtitle, loc="left", fontsize=8.0, fontweight="bold", pad=7)
    exchangeable_widths = sorted(
        {
            float(row["display_linewidth_Hz"])
            for row in sticks
            if int(row.get("is_exchangeable_H", 0)) == 1
            and row.get("display_linewidth_Hz") not in (None, "")
        }
    )
    linewidth_text = f"display FWHM {linewidth_hz:.2f} Hz"
    if exchangeable_widths:
        linewidth_text += (
            f"; exch. provisional {exchangeable_widths[0]:.2f} Hz"
            if len(exchangeable_widths) == 1
            else "; exch. provisional "
            + "/".join(f"{value:.2f}" for value in exchangeable_widths)
            + " Hz"
        )
    ax.text(
        0.995,
        0.965,
        f"{nucleus_label}; {frequency_mhz:.2f} MHz; {linewidth_text}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.2,
        color="#59636E",
    )


def individual_figure(
    figure_id: str,
    display_name: str,
    nucleus_label: str,
    variant_key: str,
    x: np.ndarray,
    y: np.ndarray,
    sticks: Sequence[dict[str, Any]],
    frequency_mhz: float,
    linewidth_hz: float,
) -> plt.Figure:
    variant = VARIANTS[variant_key]
    fig, ax = plt.subplots(figsize=(7.20, 2.55))
    fig.subplots_adjust(left=0.065, right=0.99, bottom=0.20, top=0.72)
    fig.suptitle(NOTICE, fontsize=9.0, fontweight="bold", y=0.97, color="#202428")
    subtitle = f"{figure_id} | {display_name} | {nucleus_label} | {variant['short_label']}"
    style_axis(
        ax,
        x,
        y,
        sticks,
        variant["color"],
        nucleus_label,
        subtitle,
        frequency_mhz,
        linewidth_hz,
    )
    fig.text(
        0.99,
        0.035,
        "Display broadening only; no noise, baseline, solvent peak or water peak",
        ha="right",
        va="bottom",
        fontsize=6.1,
        color="#6A737B",
    )
    return fig


def comparison_figure(
    figure_id: str,
    display_name: str,
    spectra: dict[tuple[str, str], tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]],
    h_frequency_mhz: float,
    c_frequency_mhz: float,
    h_linewidth_hz: float,
    c_linewidth_hz: float,
    carbon_na_reason: str | None,
) -> plt.Figure:
    variant_keys = [
        variant_key
        for variant_key in VARIANTS
        if any((nucleus, variant_key) in spectra for nucleus in ("1H", "13C"))
    ]
    if not variant_keys:
        raise SpectrumError("No calculated spectra are available for the comparison figure")
    fig, axes = plt.subplots(
        2,
        len(variant_keys),
        figsize=(max(4.2, 3.60 * len(variant_keys)), 5.30),
        squeeze=False,
    )
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.10, top=0.84, hspace=0.56, wspace=0.16)
    fig.suptitle(NOTICE, fontsize=9.2, fontweight="bold", y=0.97, color="#202428")
    fig.text(
        0.07,
        0.905,
        f"{figure_id} | {display_name} | calculated NMR display",
        fontsize=8.1,
        fontweight="bold",
        color="#343A40",
    )
    panel_specs = []
    panel_index = 0
    for row, (nucleus, nucleus_label, frequency, linewidth) in enumerate(
        (
            ("1H", H1_MATH_LABEL, h_frequency_mhz, h_linewidth_hz),
            ("13C", C13_MATH_LABEL, c_frequency_mhz, c_linewidth_hz),
        )
    ):
        for column, variant_key in enumerate(variant_keys):
            panel_specs.append(
                (
                    row,
                    column,
                    nucleus,
                    nucleus_label,
                    variant_key,
                    chr(ord("a") + panel_index),
                    frequency,
                    linewidth,
                )
            )
            panel_index += 1
    for row, column, nucleus, nucleus_label, variant_key, panel, frequency, linewidth in panel_specs:
        ax = axes[row, column]
        key = (nucleus, variant_key)
        if key in spectra:
            x, y, sticks = spectra[key]
            style_axis(
                ax,
                x,
                y,
                sticks,
                VARIANTS[variant_key]["color"],
                nucleus_label,
                f"{panel}  {nucleus_label} — {VARIANTS[variant_key]['short_label']}",
                frequency,
                linewidth,
            )
        else:
            ax.axis("off")
            ax.text(
                0.02,
                0.90,
                f"{panel}  {nucleus_label} — {VARIANTS[variant_key]['short_label']}",
                transform=ax.transAxes,
                fontsize=8.0,
                fontweight="bold",
                va="top",
            )
            not_available = (
                C13_MATH_LABEL + ": NOT APPLICABLE\n(no carbon in the reference structure)"
                if nucleus == "13C" and carbon_na_reason
                else f"{nucleus_label}: shift variant not available"
            )
            ax.text(
                0.5,
                0.48,
                not_available,
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=8.0,
                color="#59636E",
            )
            if carbon_na_reason:
                ax.text(
                    0.5,
                    0.27,
                    carbon_na_reason,
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=6.0,
                    color="#7A838B",
                    wrap=True,
                )
    fig.text(
        0.99,
        0.025,
        "Profiles use configurable Lorentzian display broadening only; no noise, baseline, solvent or water peaks",
        ha="right",
        va="bottom",
        fontsize=6.0,
        color="#6A737B",
    )
    return fig


def output_file_record(path: Path) -> dict[str, Any]:
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def build_readme(metadata: dict[str, Any]) -> str:
    compound = metadata["compound"]
    display = metadata["display_parameters"]
    lines = [
        NOTICE,
        "=" * len(NOTICE),
        "",
        f"Compound: {compound['figure_id']} — {compound['display_name']}",
        f"Source prediction bundle: {metadata['input_prediction_bundle']['path']}",
        "",
        "Scope",
        "-----",
        "These files visualize values already present in the calculated prediction bundle.",
        "The script did not generate new chemical shifts, J couplings, experimental peaks, noise,",
        "baseline, solvent peaks or water peaks.",
        "",
        "Display parameters (not experimental observations)",
        "--------------------------------------------------",
        f"1H frequency: {display['h_frequency_MHz']:.4f} MHz",
        f"13C frequency: {display['c_frequency_MHz']:.4f} MHz",
        f"1H Lorentzian FWHM: {display['h_lorentzian_fwhm_Hz']:.4f} Hz",
        f"Exchangeable 1H provisional display FWHM: {display['exchangeable_h_lorentzian_fwhm_Hz_provisional']:.4f} Hz",
        f"13C Lorentzian FWHM: {display['c_lorentzian_fwhm_Hz']:.4f} Hz",
        f"Coupling graph threshold: |J| >= {display['coupling_edge_threshold_Hz']:.6g} Hz",
        f"Exact Hamiltonian component limit: {display['exact_max_spins']} spin-1/2 nuclei",
        f"First-order neighbor cap: {display['first_order_max_neighbors']} per spin",
        f"Transition merge tolerance: {display['transition_merge_tolerance_Hz']:.6g} Hz",
        f"19F decoupling: {'ON' if display['decouple_19F'] else 'OFF'}",
        f"31P decoupling: {'ON' if display['decouple_31P'] else 'OFF'}",
        "13C display is proton-decoupled by default; F/P couplings remain unless explicitly decoupled.",
        "",
        "Spin-model audit",
        "----------------",
    ]
    for row in metadata["spin_component_audit"]:
        lines.append(
            f"- {row['variant']} {row['component_id']}: {row['atom_labels_1based']}; "
            f"{row['model']} ({row['model_reason']}); {row['n_output_sticks_after_merge']} sticks."
        )
        if row["approximation_note"] != "none":
            lines.append(f"  Approximation: {row['approximation_note']}")
    lines.extend(["", "Heteronuclear-coupling audit", "----------------------------"])
    lines.append(
        "The disposition below records pre-dynamic eligibility; for methyl H, the later "
        "dynamic-equivalence audit is authoritative for final symmetrization/exclusion."
    )
    if metadata["heteronuclear_coupling_audit"]:
        for row in metadata["heteronuclear_coupling_audit"]:
            lines.append(
                f"- {row['pair_label_1based']}: J={float(row['J_Hz']):.6g} Hz; "
                f"{row['disposition']}; F/P shift reference status="
                f"{row['heteronucleus_shift_reference_status']} (not required for splitting)."
            )
    else:
        lines.append("- No H-F, H-P, C-F or C-P pair was present in the source bundle.")
    lines.extend(["", "Exchangeable-proton policy", "----------------------------"])
    exchange_policy = metadata["exchangeable_proton_policy"]
    if exchange_policy["exchangeable_atom_records"]:
        for row in exchange_policy["exchangeable_atom_records"]:
            lines.append(
                f"- {row['atom_label_1based']} ({row['exchangeable_H_class']}, parent "
                f"{row['parent_element']}{int(row['parent_index_1based'])}): calculated shift retained; "
                "shift/linewidth/integral provisional; coherent fine J excluded."
            )
    else:
        lines.append("- No structurally exchangeable N-H/O-H/S-H proton was identified.")
    lines.extend(["", "Dynamic-equivalence audit", "-------------------------"])
    if metadata["dynamic_equivalence_audit"]:
        for row in metadata["dynamic_equivalence_audit"]:
            lines.append(
                f"- {row['variant']} {row['group_id']}: {row['hydrogen_labels_1based']}; "
                f"{row['status']} ({row['reason'] if 'reason' in row else row['assumption']})."
            )
    else:
        lines.append("- No strict methyl, CH2 or NH2 topology was identified.")
    if metadata["carbon_status"]["status"] == "not_applicable":
        lines.extend(
            [
                "",
                "13C status",
                "-----------",
                "13C: NOT APPLICABLE because the reference structure contains no carbon.",
            ]
        )
    lines.extend(
        [
            "",
            "Interpretive limits",
            "-------------------",
            "- Raw TMS and provisional regression spectra are separate calibration views of the same calculated shieldings.",
            "- The provisional regression source dataset and predictive RMSE were not supplied.",
            "- Exact spin Hamiltonians are used only for connected components at or below the recorded size threshold.",
            "- Larger components use a recorded first-order weak-coupling approximation; strong-coupling state mixing is omitted.",
            "- H-H pairs absent from the bundle are assumed uncoupled for display; the upstream calculation used a finite distance cutoff.",
            "- Every calculated scalar coupling touching a structurally exchangeable N-H/O-H/S-H proton is retained in the source/audit tables but excluded from the default coherent fine-structure model.",
            "- Exchangeable-proton shifts, broadening and nominal one-H areas are provisional observation-model values; they are not predictions of exchange kinetics, relaxation or quantitative experimental integrals.",
            "- Strict topological methyl groups use a fast-rotation observation model with conserved three-H area and shift centroid; the assumption can fail at low temperature or for an unusually high rotational barrier.",
            "- CH2 and NH2 hydrogens remain atom-resolved and are never automatically merged, because diastereotopy, exchange and restricted rotation require molecule- and condition-specific evidence.",
            "- H-F/H-P splitting uses high-field secular spin-1/2 F/P m-state enumeration on top of the established H-H exact/first-order model.",
            "- C-F/C-P splitting in 13C{1H} uses a first-order spin-1/2 stick construction; F/P are not automatically decoupled.",
            "- F/P chemical shifts are not needed for these J-based splittings, so a pending F/P reference shift does not suppress them.",
            "- Heteronuclear flip-flop, second-order, relaxation and detailed line-shape effects are not represented.",
            "- Lorentzian linewidth, magnetic-field frequency and plotted chemical-shift window are display parameters, not predictions.",
            "- Exchangeable-proton positions and line shapes remain sensitive to solvent, water, concentration and exchange.",
            "- These calculated spectra cannot replace experimental characterization, identity or purity data.",
        ]
    )
    if metadata["warnings"]:
        lines.extend(["", "Warnings", "--------"])
        lines.extend(f"- {value}" for value in metadata["warnings"])
    return "\n".join(lines) + "\n"


def process_bundle(path: Path, args: argparse.Namespace, arial_path: str) -> dict[str, Any]:
    raw_bundle = load_bundle(path)
    bundle, solution_equivalence_audit = apply_declared_solution_equivalence(
        raw_bundle, args.solution_equivalence_config
    )
    compound = bundle["compound"]
    figure_id = str(compound["figure_id"])
    scenario_id = str(compound.get("scenario_id", ""))
    if not scenario_id:
        raise SpectrumError(f"Prediction bundle has no scenario_id: {path}")
    display_name = str(compound["display_name"])
    variant_keys = available_variant_keys(bundle)
    output_directory = args.output_root / safe_token(figure_id) / safe_token(scenario_id)
    if output_directory.exists() and any(output_directory.iterdir()) and not args.overwrite:
        raise SpectrumError(
            f"Output directory already contains files: {output_directory}; pass --overwrite to replace this script's files"
        )
    output_directory.mkdir(parents=True, exist_ok=True)

    spectra: dict[tuple[str, str], tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]] = {}
    component_audit: list[dict[str, Any]] = []
    exchangeable_proton_coupling_audit: list[dict[str, Any]] = []
    dynamic_equivalence_audit: list[dict[str, Any]] = []
    warnings: list[str] = []
    output_files: list[Path] = []
    proton_heteronuclear, proton_heteronuclear_audit, proton_heteronuclear_warnings = (
        collect_heteronuclear_couplings(
            bundle,
            "H",
            args.decouple_f,
            args.decouple_p,
            args.coupling_edge_threshold_hz,
        )
    )
    carbon_heteronuclear, carbon_heteronuclear_audit, carbon_heteronuclear_warnings = (
        collect_heteronuclear_couplings(
            bundle,
            "C",
            args.decouple_f,
            args.decouple_p,
            args.coupling_edge_threshold_hz,
        )
    )
    heteronuclear_coupling_audit = (
        proton_heteronuclear_audit + carbon_heteronuclear_audit
    )
    warnings.extend(proton_heteronuclear_warnings)
    warnings.extend(carbon_heteronuclear_warnings)

    for variant_key in variant_keys:
        variant = VARIANTS[variant_key]
        (
            h_sticks,
            audit,
            exchange_audit,
            dynamic_audit,
            h_warnings,
        ) = build_proton_sticks(
            bundle=bundle,
            variant_key=variant_key,
            heteronuclear_couplings=proton_heteronuclear,
            frequency_mhz=args.h_frequency_mhz,
            edge_threshold_hz=args.coupling_edge_threshold_hz,
            exact_max_spins=args.exact_max_spins,
            first_order_max_neighbors=args.first_order_max_neighbors,
            merge_tolerance_hz=args.transition_merge_tolerance_hz,
            transition_intensity_cutoff=args.transition_intensity_cutoff,
            regular_linewidth_hz=args.h_linewidth_hz,
            exchangeable_linewidth_hz=args.exchangeable_h_linewidth_hz,
        )
        component_audit.extend(audit)
        exchangeable_proton_coupling_audit.extend(exchange_audit)
        dynamic_equivalence_audit.extend(dynamic_audit)
        warnings.extend(h_warnings)
        x_h, y_h = lorentzian_profile(
            h_sticks,
            args.h_range[0],
            args.h_range[1],
            args.profile_points,
            args.h_linewidth_hz,
            args.h_frequency_mhz,
        )
        spectra[("1H", variant_key)] = (x_h, y_h, h_sticks)
        h_sticks_path = output_directory / f"{figure_id}_1H_{variant_key}_sticks.csv"
        atomic_write_csv(h_sticks_path, h_sticks, field_union(h_sticks))
        output_files.append(h_sticks_path)
        h_appearance_rows = appearance_line_rows(
            h_sticks, figure_id, scenario_id, "1H", variant_key
        )
        h_appearance_path = (
            output_directory / f"{figure_id}_{safe_token(scenario_id)}_1H_{variant_key}_appearance_lines.csv"
        )
        atomic_write_csv(
            h_appearance_path, h_appearance_rows, field_union(h_appearance_rows)
        )
        output_files.append(h_appearance_path)
        if args.render_appearance:
            from render_experimental_appearance import render as render_appearance

            rendered = render_appearance(
                h_appearance_path,
                args.appearance_config,
                output_directory / "SYNTHETIC_ACQUISITION_APPEARANCE" / "1H" / variant_key,
                args.appearance_profile_1h or args.appearance_profile,
                args.appearance_seed,
            )
            output_files.extend(rendered.values())
        h_profile_path = output_directory / f"{figure_id}_1H_{variant_key}_profile.txt"
        h_profile_metadata: dict[str, Any] = {
            "figure_id": figure_id,
            "nucleus": "1H",
            "shift_variant": variant["label"],
            "frequency_MHz_display_parameter": args.h_frequency_mhz,
            "lorentzian_FWHM_Hz_display_parameter": args.h_linewidth_hz,
            "exchangeable_H_lorentzian_FWHM_Hz_provisional_display_parameter": (
                args.exchangeable_h_linewidth_hz
            ),
            "exchangeable_H_shift_status": "provisional",
            "exchangeable_H_linewidth_status": "provisional_not_predicted",
            "exchangeable_H_integral_status": "provisional_nominal_explicit_H_count",
            "exchangeable_H_coherent_fine_J_network": "excluded_by_default",
            "noise": "none",
            "baseline": "none",
            "solvent_peak": "none",
            "water_peak": "none",
        }
        if proton_heteronuclear_audit or args.decouple_f or args.decouple_p:
            h_profile_metadata.update(
                {
                    "19F_decoupled": args.decouple_f,
                    "31P_decoupled": args.decouple_p,
                    "heteronuclear_model": (
                        "high-field secular weak-heteronuclear m-state enumeration"
                    ),
                    "F_or_P_chemical_shift_required_for_splitting": False,
                }
            )
        atomic_write_profile(
            h_profile_path,
            x_h,
            y_h,
            h_profile_metadata,
        )
        output_files.append(h_profile_path)
        fig_h = individual_figure(
            figure_id,
            display_name,
            H1_MATH_LABEL,
            variant_key,
            x_h,
            y_h,
            h_sticks,
            args.h_frequency_mhz,
            args.h_linewidth_hz,
        )
        output_files.extend(
            atomic_save_figure(
                fig_h,
                output_directory / f"{figure_id}_1H_{variant_key}_calculated",
                args.dpi,
            )
        )
        plt.close(fig_h)

    carbon_rows = [row for row in bundle["averaged_nuclei"] if row["element"] == "C"]
    carbon_na = figure_id in NO_CARBON_FIGURE_IDS or not carbon_rows
    carbon_na_reason: str | None = None
    if carbon_na:
        carbon_na_reason = (
            f"{figure_id} is an explicitly carbon-free reference state."
            if figure_id in NO_CARBON_FIGURE_IDS
            else "No carbon shift was present in the prediction bundle."
        )
        na_path = output_directory / f"{figure_id}_13C_NOT_APPLICABLE.txt"
        atomic_write_text(
            na_path,
            f"{NOTICE}\n\n13C: NOT APPLICABLE\n{carbon_na_reason}\n"
            "No 13C sticks, profile or experimental-looking spectrum was generated.\n",
        )
        output_files.append(na_path)
    else:
        for variant_key in variant_keys:
            variant = VARIANTS[variant_key]
            c_sticks = build_carbon_sticks(
                bundle,
                variant_key,
                carbon_heteronuclear,
                args.c_frequency_mhz,
                args.transition_merge_tolerance_hz,
            )
            x_c, y_c = lorentzian_profile(
                c_sticks,
                args.c_range[0],
                args.c_range[1],
                args.profile_points,
                args.c_linewidth_hz,
                args.c_frequency_mhz,
            )
            spectra[("13C", variant_key)] = (x_c, y_c, c_sticks)
            c_sticks_path = output_directory / f"{figure_id}_13C_{variant_key}_sticks.csv"
            atomic_write_csv(c_sticks_path, c_sticks, field_union(c_sticks))
            output_files.append(c_sticks_path)
            c_appearance_rows = appearance_line_rows(
                c_sticks, figure_id, scenario_id, "13C", variant_key
            )
            c_appearance_path = (
                output_directory
                / f"{figure_id}_{safe_token(scenario_id)}_13C_{variant_key}_appearance_lines.csv"
            )
            atomic_write_csv(
                c_appearance_path, c_appearance_rows, field_union(c_appearance_rows)
            )
            output_files.append(c_appearance_path)
            if args.render_appearance:
                from render_experimental_appearance import render as render_appearance

                rendered = render_appearance(
                    c_appearance_path,
                    args.appearance_config,
                    output_directory / "SYNTHETIC_ACQUISITION_APPEARANCE" / "13C" / variant_key,
                    args.appearance_profile_13c or args.appearance_profile,
                    args.appearance_seed,
                )
                output_files.extend(rendered.values())
            c_profile_path = output_directory / f"{figure_id}_13C_{variant_key}_profile.txt"
            c_profile_metadata: dict[str, Any] = {
                "figure_id": figure_id,
                "nucleus": "13C{1H} proton-decoupled display",
                "shift_variant": variant["label"],
                "frequency_MHz_display_parameter": args.c_frequency_mhz,
                "lorentzian_FWHM_Hz_display_parameter": args.c_linewidth_hz,
                "noise": "none",
                "baseline": "none",
                "solvent_peak": "none",
                "water_peak": "none",
            }
            if carbon_heteronuclear_audit or args.decouple_f or args.decouple_p:
                c_profile_metadata.update(
                    {
                        "19F_decoupled": args.decouple_f,
                        "31P_decoupled": args.decouple_p,
                        "heteronuclear_model": "high-field first-order spin-1/2 splitting",
                        "F_or_P_chemical_shift_required_for_splitting": False,
                    }
                )
            atomic_write_profile(
                c_profile_path,
                x_c,
                y_c,
                c_profile_metadata,
            )
            output_files.append(c_profile_path)
            fig_c = individual_figure(
                figure_id,
                display_name,
                C13_H1_DECOUPLED_MATH_LABEL,
                variant_key,
                x_c,
                y_c,
                c_sticks,
                args.c_frequency_mhz,
                args.c_linewidth_hz,
            )
            output_files.extend(
                atomic_save_figure(
                    fig_c,
                    output_directory / f"{figure_id}_13C_{variant_key}_calculated",
                    args.dpi,
                )
            )
            plt.close(fig_c)

    comparison = comparison_figure(
        figure_id,
        display_name,
        spectra,
        args.h_frequency_mhz,
        args.c_frequency_mhz,
        args.h_linewidth_hz,
        args.c_linewidth_hz,
        carbon_na_reason,
    )
    output_files.extend(
        atomic_save_figure(
            comparison,
            output_directory / f"{figure_id}_NMR_prediction_comparison",
            args.dpi,
        )
    )
    plt.close(comparison)

    audit_path = output_directory / f"{figure_id}_spin_component_audit.csv"
    atomic_write_csv(audit_path, component_audit, field_union(component_audit))
    output_files.append(audit_path)
    exchange_audit_path = (
        output_directory / f"{figure_id}_exchangeable_proton_coupling_audit.csv"
    )
    atomic_write_csv(
        exchange_audit_path,
        exchangeable_proton_coupling_audit,
        field_union(exchangeable_proton_coupling_audit)
        or ["notice", "variant", "pair_label_1based", "J_Hz", "disposition"],
    )
    output_files.append(exchange_audit_path)
    dynamic_audit_path = output_directory / f"{figure_id}_dynamic_equivalence_audit.csv"
    atomic_write_csv(
        dynamic_audit_path,
        dynamic_equivalence_audit,
        field_union(dynamic_equivalence_audit)
        or ["notice", "variant", "group_id", "status", "reason"],
    )
    output_files.append(dynamic_audit_path)
    solution_equivalence_path = (
        output_directory / f"{figure_id}_declared_solution_equivalence_audit.csv"
    )
    atomic_write_csv(
        solution_equivalence_path,
        solution_equivalence_audit,
        field_union(solution_equivalence_audit)
        or ["notice", "figure_id", "operation", "group_id", "raw_source_retained"],
    )
    output_files.append(solution_equivalence_path)

    metadata: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "notice": NOTICE,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "generated_by": "scripts/generate_nmr_spectra.py",
        "compound": compound,
        "scenario_id": scenario_id,
        "input_prediction_bundle": {
            "path": str(path),
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        },
        "data_policy": {
            "shift_source": "averaged_nuclei values in the input prediction bundle",
            "coupling_source": (
                "averaged_JHH for the established H-H model; averaged_scalar_couplings "
                "for H-F, H-P, C-F and C-P splitting"
            ),
            "averaged_scalar_couplings_key_present": (
                "averaged_scalar_couplings" in bundle
            ),
            "fp_chemical_shift_required_for_heteronuclear_splitting": False,
            "new_random_or_synthetic_measurements": False,
            "noise_added": False,
            "baseline_added": False,
            "solvent_peak_added": False,
            "water_peak_added": False,
            "separate_synthetic_appearance_rendered": bool(args.render_appearance),
        },
        "display_parameters": {
            "h_frequency_MHz": args.h_frequency_mhz,
            "c_frequency_MHz": args.c_frequency_mhz,
            "h_lorentzian_fwhm_Hz": args.h_linewidth_hz,
            "exchangeable_h_lorentzian_fwhm_Hz_provisional": (
                args.exchangeable_h_linewidth_hz
            ),
            "c_lorentzian_fwhm_Hz": args.c_linewidth_hz,
            "h_range_ppm": list(args.h_range),
            "c_range_ppm": list(args.c_range),
            "profile_points": args.profile_points,
            "coupling_edge_threshold_Hz": args.coupling_edge_threshold_hz,
            "exact_max_spins": args.exact_max_spins,
            "first_order_max_neighbors": args.first_order_max_neighbors,
            "transition_merge_tolerance_Hz": args.transition_merge_tolerance_hz,
            "transition_intensity_cutoff": args.transition_intensity_cutoff,
            "decouple_19F": args.decouple_f,
            "decouple_31P": args.decouple_p,
            "default_13C_decoupling": (
                "1H decoupled; 19F and 31P coupled unless their explicit CLI flags are set"
            ),
            "png_dpi": args.dpi,
            "font": "Arial",
            "font_file_resolved_for_rendering": arial_path,
        },
        "shift_variants": {variant_key: VARIANTS[variant_key] for variant_key in variant_keys},
        "unavailable_shift_variants": [
            variant_key for variant_key in VARIANTS if variant_key not in variant_keys
        ],
        "spin_component_audit": component_audit,
        "exchangeable_proton_policy": {
            "structural_definition": "explicit H directly bonded to N, O or S",
            "exchangeable_atom_records": list(exchangeable_proton_map(bundle).values()),
            "calculated_shieldings_retained": True,
            "calculated_scalar_couplings_retained": True,
            "coherent_fine_structure_default": "excluded for all J pairs touching exchangeable H",
            "shift_status": "provisional",
            "linewidth_status": "provisional display parameter, not predicted",
            "integral_status": "provisional nominal explicit-H count, not quantitative",
        },
        "exchangeable_proton_coupling_audit": exchangeable_proton_coupling_audit,
        "dynamic_equivalence_policy": {
            "methyl_rule": (
                "neutral closed-shell carbon with exactly three explicit H and one heavy-atom "
                "neighbor: average three shifts; remove internal equivalent-H J; symmetrize "
                "complete external J triplets; exclude incomplete external triplets from display"
            ),
            "methyl_assumption": (
                "fast methyl rotation on the solution-NMR timescale; may fail at very low "
                "temperature or for unusually high rotational barriers"
            ),
            "CH2_policy": "never merged automatically because enantiotopic/diastereotopic status is context-dependent",
            "NH2_policy": "never merged automatically because exchange and restricted rotation are context-dependent",
            "calculated_atom_resolved_source_values_retained_in_prediction_bundle": True,
        },
        "dynamic_equivalence_audit": dynamic_equivalence_audit,
        "declared_solution_equivalence": {
            "config_file": str(args.solution_equivalence_config),
            "config_sha256": sha256_file(args.solution_equivalence_config),
            "operations": solution_equivalence_audit,
            "raw_atom_resolved_prediction_bundle_unchanged": True,
        },
        "heteronuclear_coupling_audit": heteronuclear_coupling_audit,
        "heteronuclear_coupling_audit_scope": (
            "pre-dynamic eligibility/decoupling audit; for methyl H the dynamic_equivalence_audit "
            "is authoritative for final J symmetrization or incomplete-triplet exclusion"
        ),
        "heteronuclear_model": {
            "proton": (
                "existing H-H exact/first-order model evaluated for every equally populated "
                "active F/P spin-1/2 m-state combination in the high-field secular weak-"
                "heteronuclear limit"
            ),
            "carbon": (
                "13C{1H} sticks split sequentially by active C-F/C-P J values in the "
                "high-field first-order spin-1/2 limit"
            ),
            "limitations": (
                "heteronuclear flip-flop, second-order, relaxation and line-shape effects are omitted"
            ),
        },
        "carbon_status": {
            "status": "not_applicable" if carbon_na else "calculated_profile_generated",
            "reason": carbon_na_reason,
        },
        "warnings": sorted(set(warnings)),
        "upstream_weighting": bundle["weighting"],
    }
    readme_path = output_directory / f"{figure_id}_spectrum_generation_README.txt"
    atomic_write_text(readme_path, build_readme(metadata))
    output_files.append(readme_path)
    metadata["output_files"] = [output_file_record(file) for file in output_files]
    metadata_path = output_directory / f"{figure_id}_spectrum_generation_metadata.json"
    atomic_write_json(metadata_path, metadata)
    output_files.append(metadata_path)
    return {
        "notice": NOTICE,
        "figure_id": figure_id,
        "display_name": display_name,
        "status": "calculated_spectra_generated_not_experimental",
        "carbon_status": metadata["carbon_status"]["status"],
        "n_spin_components_per_variant": len(component_audit) // len(variant_keys),
        "n_output_files": len(output_files),
        "output_directory": str(output_directory),
        "scenario_id": scenario_id,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Render explicitly labelled calculated 1H/13C spectra from prediction-bundle JSON. "
            "The clean layer is always retained; an explicitly labelled synthetic acquisition "
            "appearance can be rendered as a separate output layer."
        )
    )
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--solution-equivalence-config",
        type=Path,
        default=DEFAULT_SOLUTION_EQUIVALENCE,
        help="Identity-locked solution-equivalence groups; raw bundle values remain unchanged.",
    )
    parser.add_argument("--only", nargs="+", help="Optional Figure3 IDs, e.g. --only F3_001")
    parser.add_argument("--h-frequency-mhz", type=float, default=400.0)
    parser.add_argument("--c-frequency-mhz", type=float, default=100.6)
    parser.add_argument("--h-linewidth-hz", type=float, default=1.0)
    parser.add_argument(
        "--exchangeable-h-linewidth-hz",
        type=float,
        default=12.0,
        help=(
            "Provisional display-only Lorentzian FWHM for structurally exchangeable "
            "N-H/O-H/S-H signals (default: 12 Hz; not a relaxation/exchange prediction)."
        ),
    )
    parser.add_argument("--c-linewidth-hz", type=float, default=1.5)
    parser.add_argument("--h-range", type=float, nargs=2, default=(0.0, 12.0), metavar=("MIN", "MAX"))
    parser.add_argument("--c-range", type=float, nargs=2, default=(0.0, 220.0), metavar=("MIN", "MAX"))
    parser.add_argument("--profile-points", type=int, default=50000)
    parser.add_argument("--coupling-edge-threshold-hz", type=float, default=0.05)
    parser.add_argument("--exact-max-spins", type=int, default=7)
    parser.add_argument("--first-order-max-neighbors", type=int, default=12)
    parser.add_argument("--transition-merge-tolerance-hz", type=float, default=0.02)
    parser.add_argument("--transition-intensity-cutoff", type=float, default=1.0e-10)
    parser.add_argument(
        "--decouple-f",
        action="store_true",
        help="Omit calculated H-F/C-F splitting as an explicit 19F-decoupled display (default: coupled).",
    )
    parser.add_argument(
        "--decouple-p",
        action="store_true",
        help="Omit calculated H-P/C-P splitting as an explicit 31P-decoupled display (default: coupled).",
    )
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument(
        "--render-appearance",
        action="store_true",
        help=(
            "Also render separately labelled synthetic acquisition appearance outputs with "
            "noise, baseline, solvent/water components and empirical linewidths."
        ),
    )
    parser.add_argument(
        "--appearance-config", type=Path, default=DEFAULT_APPEARANCE_CONFIG
    )
    parser.add_argument(
        "--appearance-profile",
        help=(
            "Optional shared appearance profile override. Prefer the nucleus-specific options "
            "when both 1H and 13C are generated."
        ),
    )
    parser.add_argument("--appearance-profile-1h", help="Optional explicit 1H appearance profile.")
    parser.add_argument("--appearance-profile-13c", help="Optional explicit 13C appearance profile.")
    parser.add_argument("--appearance-seed", type=int)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace files previously written by this script in selected compound directories",
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.solution_equivalence_config.is_file():
        raise SpectrumError(
            f"Solution-equivalence config is missing: {args.solution_equivalence_config}"
        )
    if args.render_appearance and not args.appearance_config.is_file():
        raise SpectrumError(f"Appearance config is missing: {args.appearance_config}")
    if args.h_frequency_mhz <= 0.0 or args.c_frequency_mhz <= 0.0:
        raise SpectrumError("Frequencies must be positive")
    if (
        args.h_linewidth_hz <= 0.0
        or args.exchangeable_h_linewidth_hz <= 0.0
        or args.c_linewidth_hz <= 0.0
    ):
        raise SpectrumError("Lorentzian linewidths must be positive")
    if args.exchangeable_h_linewidth_hz < args.h_linewidth_hz:
        raise SpectrumError(
            "--exchangeable-h-linewidth-hz must be at least --h-linewidth-hz"
        )
    if args.h_range[1] <= args.h_range[0] or args.c_range[1] <= args.c_range[0]:
        raise SpectrumError("Spectrum range maximum must exceed minimum")
    if args.profile_points < 1000:
        raise SpectrumError("--profile-points must be at least 1000")
    if args.coupling_edge_threshold_hz < 0.0:
        raise SpectrumError("--coupling-edge-threshold-hz cannot be negative")
    if args.exact_max_spins < 1 or args.exact_max_spins > 10:
        raise SpectrumError("--exact-max-spins must be between 1 and 10")
    if args.first_order_max_neighbors < 0 or args.first_order_max_neighbors > 20:
        raise SpectrumError("--first-order-max-neighbors must be between 0 and 20")
    if args.transition_merge_tolerance_hz < 0.0:
        raise SpectrumError("--transition-merge-tolerance-hz cannot be negative")
    if args.transition_intensity_cutoff <= 0.0:
        raise SpectrumError("--transition-intensity-cutoff must be positive")
    if args.dpi < 150:
        raise SpectrumError("--dpi must be at least 150")


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.input_root = args.input_root.expanduser().resolve()
    args.output_root = args.output_root.expanduser().resolve()
    args.solution_equivalence_config = args.solution_equivalence_config.expanduser().resolve()
    args.appearance_config = args.appearance_config.expanduser().resolve()
    try:
        validate_args(args)
        arial_path = configure_matplotlib()
        bundle_paths = discover_bundles(args.input_root, set(args.only or []))
        summaries = [process_bundle(path, args, arial_path) for path in bundle_paths]
    except (OSError, ValueError, SpectrumError, np.linalg.LinAlgError) as exc:
        print(f"Spectrum generation failed: {exc}")
        return 2

    label = "ALL" if not args.only else "_".join(sorted(set(args.only)))
    args.output_root.mkdir(parents=True, exist_ok=True)
    summary = {
        "schema_version": SCHEMA_VERSION,
        "notice": NOTICE,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "selected_figure_ids": sorted({row["figure_id"] for row in summaries}),
        "selected_figure_scenarios": [
            {"figure_id": row["figure_id"], "scenario_id": row["scenario_id"]}
            for row in summaries
        ],
        "results": summaries,
    }
    atomic_write_json(args.output_root / f"spectrum_generation_summary__{label}.json", summary)
    atomic_write_csv(
        args.output_root / f"spectrum_generation_summary__{label}.csv",
        summaries,
        list(summaries[0]),
    )
    print(NOTICE)
    print(f"Generated calculated spectra for {len(summaries)} compound(s)")
    print(args.output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
