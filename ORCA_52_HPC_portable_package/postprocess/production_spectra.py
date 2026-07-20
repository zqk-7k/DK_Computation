"""Reproducible post-processing for calculated NMR and formula-level HRMS.

This module deliberately separates three NMR layers:

P0  ORCA-derived calculated transitions (sticks)
P1  clean broadened calculated analyte spectrum
P2  synthetic acquisition appearance, never acquired experimentally

The HRMS branch derives exact masses and natural-abundance isotope envelopes
from an explicitly declared ion formula.  ORCA is not used for this branch and
no LC-MS fragmentation, retention time, abundance, or measured ``found`` value
is predicted.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Mandatory editable-text publication settings.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["legend.frameon"] = False

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
SCRIPTS = PROJECT / "scripts"
CALIBRATION = PROJECT / "data" / "real_nmr_calibration_audit"
if not CALIBRATION.is_dir():
    # Portable HPC package layout. Keep this relative to the unpacked package
    # root so the archive works under any Linux filesystem path.
    CALIBRATION = PROJECT / "real_spectrum_calibration"

if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from orca_io import (  # noqa: E402
    AtomRecord,
    CouplingRecord,
    OrcaParseError,
    coupling_matrix,
    parse_j_couplings,
    parse_shieldings,
    read_atom_map,
)


CALCULATED = "CALCULATED_NOT_EXPERIMENTAL"
SYNTHETIC = "SYNTHETIC_APPEARANCE_NOT_ACQUIRED"
TRANSITION_STATUS = "CALCULATED"
HRMS_SCOPE = (
    "FORMULA-DERIVED PRECURSOR/ADDUCT EXACT MASS AND NATURAL-ABUNDANCE "
    "ISOTOPE ENVELOPE ONLY; NO LC-MS FRAGMENTATION OR RETENTION TIME"
)
ELEMENT_BY_NUCLEUS = {"1H": "H", "13C": "C", "19F": "F", "31P": "P"}
FREQUENCY_DEFAULTS = {"1H": 400.13, "13C": 100.61, "19F": 376.50, "31P": 162.00}
LINEWIDTH_DEFAULTS = {"1H": 3.32, "13C": 2.40, "19F": 2.00, "31P": 2.00}
EXCHANGEABLE_H_LINEWIDTH_HZ = 18.16
PALETTE = {
    "analyte": "#0F4D92",
    "sticks": "#272727",
    "accent": "#B64342",
    "neutral": "#767676",
    "hrms": "#42949E",
}


class PostprocessError(RuntimeError):
    pass


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{hashlib.sha256(text.encode()).hexdigest()[:12]}.tmp")
    temporary.write_text(text, encoding="utf-8", newline="\n")
    temporary.replace(path)


def write_json(path: Path, payload: Any) -> None:
    atomic_write(path, json.dumps(payload, ensure_ascii=False, indent=2) + "\n")


def write_rows(path: Path, rows: Sequence[Mapping[str, Any]], delimiter: str = ",") -> None:
    if not rows:
        raise PostprocessError(f"Refusing to write an empty table: {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    import io

    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter=delimiter, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
    atomic_write(path, "\ufeff" + stream.getvalue())


def resolve_path(config_path: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = config_path.parent / path
    return path.resolve()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise PostprocessError(f"Cannot import module from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def spin_operators(n_spins: int) -> tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray]]:
    identity = np.eye(2, dtype=complex)
    ix0 = 0.5 * np.array([[0, 1], [1, 0]], dtype=complex)
    iy0 = 0.5 * np.array([[0, -1j], [1j, 0]], dtype=complex)
    iz0 = 0.5 * np.array([[1, 0], [0, -1]], dtype=complex)

    def embed(single: np.ndarray, position: int) -> np.ndarray:
        result = np.array([[1.0 + 0.0j]])
        for index in range(n_spins):
            result = np.kron(result, single if index == position else identity)
        return result

    return (
        [embed(ix0, index) for index in range(n_spins)],
        [embed(iy0, index) for index in range(n_spins)],
        [embed(iz0, index) for index in range(n_spins)],
    )


def exact_spin_sticks(
    shifts_ppm: np.ndarray,
    j_hz: np.ndarray,
    frequency_mhz: float,
    cutoff: float = 1.0e-10,
) -> list[dict[str, float]]:
    n_spins = len(shifts_ppm)
    ix, iy, iz = spin_operators(n_spins)
    hamiltonian = np.zeros_like(ix[0])
    for index, shift in enumerate(shifts_ppm):
        hamiltonian += shift * frequency_mhz * iz[index]
    for i in range(n_spins):
        for j in range(i + 1, n_spins):
            hamiltonian += j_hz[i, j] * (
                ix[i] @ ix[j] + iy[i] @ iy[j] + iz[i] @ iz[j]
            )
    energies, vectors = np.linalg.eigh(hamiltonian)
    detector = vectors.conj().T @ sum(ix) @ vectors
    rows: list[dict[str, float]] = []
    for low in range(len(energies)):
        for high in range(low + 1, len(energies)):
            intensity = float(abs(detector[low, high]) ** 2)
            if intensity <= cutoff:
                continue
            frequency_hz = float(energies[high] - energies[low])
            ppm = frequency_hz / frequency_mhz
            if -1000.0 < ppm < 1000.0:
                rows.append(
                    {
                        "shift_ppm": ppm,
                        "frequency_Hz": frequency_hz,
                        "intensity": intensity,
                    }
                )
    total = sum(row["intensity"] for row in rows)
    if total <= 0:
        raise PostprocessError("Exact-spin Hamiltonian yielded no allowed transitions")
    for row in rows:
        row["intensity"] *= n_spins / total
    return sorted(rows, key=lambda row: row["shift_ppm"])


def merge_sticks(
    rows: Sequence[dict[str, float]], tolerance_hz: float, frequency_mhz: float
) -> list[dict[str, float]]:
    tolerance_ppm = tolerance_hz / frequency_mhz
    merged: list[dict[str, float]] = []
    for row in sorted(rows, key=lambda item: item["shift_ppm"]):
        if merged and abs(row["shift_ppm"] - merged[-1]["shift_ppm"]) <= tolerance_ppm:
            old = merged[-1]["intensity"]
            total = old + row["intensity"]
            merged[-1]["shift_ppm"] = (
                merged[-1]["shift_ppm"] * old + row["shift_ppm"] * row["intensity"]
            ) / total
            merged[-1]["intensity"] = total
            merged[-1]["frequency_Hz"] = merged[-1]["shift_ppm"] * frequency_mhz
        else:
            merged.append(dict(row))
    return merged


def first_order_sticks(
    shifts_ppm: Sequence[float],
    j_hz: np.ndarray,
    frequency_mhz: float,
    threshold_hz: float,
    max_neighbors: int,
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for spin, center in enumerate(shifts_ppm):
        neighbors = sorted(
            (abs(j_hz[spin, other]), j_hz[spin, other], other)
            for other in range(len(shifts_ppm))
            if other != spin and abs(j_hz[spin, other]) >= threshold_hz
        )[-max_neighbors:]
        local = [(float(center), 1.0)]
        for _, coupling, _ in neighbors:
            delta = coupling / (2.0 * frequency_mhz)
            local = [
                item
                for shift, intensity in local
                for item in (
                    (shift - delta, intensity / 2.0),
                    (shift + delta, intensity / 2.0),
                )
            ]
        for shift, intensity in local:
            rows.append(
                {
                    "shift_ppm": shift,
                    "frequency_Hz": shift * frequency_mhz,
                    "intensity": intensity,
                    "spin_local_index": str(spin),
                }
            )
    return rows


def pseudo_voigt_area(
    axis: np.ndarray, center: float, fwhm_ppm: float, area: float, eta: float = 0.65
) -> np.ndarray:
    if fwhm_ppm <= 0:
        raise PostprocessError("Line width must be positive")
    gamma = 0.5 * fwhm_ppm
    lorentz = gamma / math.pi / ((axis - center) ** 2 + gamma**2)
    sigma = fwhm_ppm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    gaussian = np.exp(-0.5 * ((axis - center) / sigma) ** 2) / (
        sigma * math.sqrt(2.0 * math.pi)
    )
    return area * (eta * lorentz + (1.0 - eta) * gaussian)


def clean_profile(
    sticks: Sequence[Mapping[str, Any]],
    frequency_mhz: float,
    low: float,
    high: float,
    points: int,
) -> tuple[np.ndarray, np.ndarray]:
    axis = np.linspace(low, high, points, dtype=float)
    signal = np.zeros_like(axis)
    for row in sticks:
        signal += pseudo_voigt_area(
            axis,
            float(row["shift_ppm"]),
            float(row["linewidth_hz"]) / frequency_mhz,
            float(row["intensity"]),
        )
    if signal.max() > 0:
        signal /= signal.max()
    return axis, signal


def save_nmr_figure(
    stem: Path,
    compound_id: str,
    display_name: str,
    nucleus: str,
    solvent: str,
    sticks: Sequence[Mapping[str, Any]],
    axis: np.ndarray,
    signal: np.ndarray,
    dpi: int,
) -> list[Path]:
    fig, ax = plt.subplots(figsize=(7.20, 2.15))
    ax.plot(axis, signal, color=PALETTE["analyte"], lw=0.85, zorder=2)
    baseline = -0.055
    scale = 0.05 / max(float(row["intensity"]) for row in sticks)
    for row in sticks:
        ax.vlines(
            float(row["shift_ppm"]),
            baseline,
            baseline + float(row["intensity"]) * scale,
            color=PALETTE["sticks"],
            lw=0.45,
            zorder=1,
        )
    ax.axhline(0, color=PALETTE["neutral"], lw=0.45)
    ax.set_xlim(float(axis.max()), float(axis.min()))
    ax.set_ylim(-0.075, 1.10)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(r"$\delta$ (ppm)")
    ax.set_title(f"{compound_id} | {display_name} | calculated {nucleus} NMR", loc="left", weight="bold")
    ax.text(
        0.995,
        0.965,
        f"P0 sticks + P1 clean model; {solvent}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        color=PALETTE["neutral"],
        fontsize=6,
    )
    ax.text(
        0.995,
        0.06,
        "CALCULATED MODEL — NOT EXPERIMENTAL",
        transform=ax.transAxes,
        ha="right",
        color=PALETTE["accent"],
        fontsize=6,
        weight="bold",
    )
    fig.tight_layout(pad=1.2)
    outputs: list[Path] = []
    for suffix, kwargs in ((".svg", {}), (".png", {"dpi": dpi})):
        path = stem.with_suffix(suffix)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, bbox_inches="tight", **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def display_range(config: Mapping[str, Any], nucleus: str, shifts: Sequence[float]) -> tuple[float, float]:
    ranges = config.get("display_ranges_ppm", {})
    if nucleus in ranges:
        low, high = (float(item) for item in ranges[nucleus])
        return min(low, high), max(low, high)
    if nucleus == "1H":
        return 0.0, 12.0
    if nucleus == "13C":
        return 0.0, 220.0
    minimum, maximum = min(shifts), max(shifts)
    margin = max(5.0, 0.15 * max(1.0, maximum - minimum))
    return minimum - margin, maximum + margin


def make_standard_transitions(
    nucleus: str,
    atoms: Sequence[AtomRecord],
    shifts: Mapping[int, float],
    j_records: Sequence[CouplingRecord],
    frequency_mhz: float,
    solvent: str,
    exact_max_spins: int,
    transition_merge_hz: float,
) -> tuple[list[dict[str, Any]], str]:
    target_element = ELEMENT_BY_NUCLEUS[nucleus]
    selected = [
        atom
        for atom in atoms
        if atom.element == target_element and atom.include_in_analyte_nmr
    ]
    if nucleus == "1H" and solvent.upper().replace("₂", "2") == "D2O":
        selected = [atom for atom in selected if not atom.exchangeable]
    if not selected:
        raise PostprocessError(f"No observable analyte atoms for {nucleus} in {solvent}")

    default_width = LINEWIDTH_DEFAULTS[nucleus]
    rows: list[dict[str, Any]] = []
    method = "direct calculated singlets"
    if nucleus == "1H":
        stable = [atom for atom in selected if not atom.exchangeable]
        exchangeable = [atom for atom in selected if atom.exchangeable]
        if stable:
            indices = [atom.index0 for atom in stable]
            matrix = np.asarray(coupling_matrix(indices, j_records), dtype=float)
            stable_shifts = np.asarray([shifts[atom.index0] for atom in stable], dtype=float)
            if len(stable) <= exact_max_spins:
                raw = exact_spin_sticks(stable_shifts, matrix, frequency_mhz)
                raw = merge_sticks(raw, transition_merge_hz, frequency_mhz)
                method = f"exact {len(stable)}-spin Hamiltonian for non-exchangeable 1H"
                for sequence, item in enumerate(raw, 1):
                    rows.append(
                        {
                            "transition_id": f"H_EXACT_{sequence:05d}",
                            "component_id": "analyte_nonexchangeable_1H_exact_system",
                            "nucleus": nucleus,
                            "shift_ppm": f"{item['shift_ppm']:.9f}",
                            "intensity": f"{item['intensity']:.12g}",
                            "linewidth_hz": f"{default_width:.6f}",
                            "assignment": ";".join(atom.label for atom in stable),
                            "data_status": TRANSITION_STATUS,
                        }
                    )
            else:
                raw = first_order_sticks(
                    stable_shifts,
                    matrix,
                    frequency_mhz,
                    threshold_hz=max(0.20, 0.2 * default_width),
                    max_neighbors=12,
                )
                method = f"first-order fallback for {len(stable)} non-exchangeable 1H"
                for sequence, item in enumerate(raw, 1):
                    local = int(str(item["spin_local_index"]))
                    rows.append(
                        {
                            "transition_id": f"H_FIRST_{sequence:05d}",
                            "component_id": f"analyte_{stable[local].label}",
                            "nucleus": nucleus,
                            "shift_ppm": f"{float(item['shift_ppm']):.9f}",
                            "intensity": f"{float(item['intensity']):.12g}",
                            "linewidth_hz": f"{default_width:.6f}",
                            "assignment": stable[local].label,
                            "data_status": TRANSITION_STATUS,
                        }
                    )
        for atom in exchangeable:
            rows.append(
                {
                    "transition_id": f"H_EXCHANGE_{atom.index0 + 1:04d}",
                    "component_id": f"analyte_{atom.label}",
                    "nucleus": nucleus,
                    "shift_ppm": f"{shifts[atom.index0]:.9f}",
                    "intensity": "1",
                    "linewidth_hz": f"{EXCHANGEABLE_H_LINEWIDTH_HZ:.6f}",
                    "assignment": f"{atom.label}; empirical exchange broadening, not ORCA-predicted",
                    "data_status": TRANSITION_STATUS,
                }
            )
    else:
        for atom in selected:
            rows.append(
                {
                    "transition_id": f"{nucleus}_{atom.index0 + 1:04d}",
                    "component_id": f"analyte_{atom.label}",
                    "nucleus": nucleus,
                    "shift_ppm": f"{shifts[atom.index0]:.9f}",
                    "intensity": "1",
                    "linewidth_hz": f"{default_width:.6f}",
                    "assignment": atom.label,
                    "data_status": TRANSITION_STATUS,
                }
            )
    if not rows:
        raise PostprocessError(f"No transitions produced for {nucleus}")
    return rows, method


def invoke_appearance(
    transitions: Path,
    compound_id: str,
    display_name: str,
    nucleus: str,
    solvent: str,
    output_root: Path,
    seed: int,
    dpi: int,
) -> Path:
    module = load_module(
        "nmr_acquisition_appearance_postprocess",
        SCRIPTS / "nmr_acquisition_appearance.py",
    )
    args = module.parser().parse_args(
        [
            "--input",
            str(transitions),
            "--spectrum-id",
            compound_id,
            "--display-name",
            display_name,
            "--nucleus",
            nucleus,
            "--solvent",
            solvent,
            "--calibration-dir",
            str(CALIBRATION),
            "--output-root",
            str(output_root),
            "--seed",
            str(seed),
            "--dpi",
            str(dpi),
            "--overwrite",
        ]
    )
    return module.run_generation(args)


def process_nmr(config_path: Path, config: Mapping[str, Any], output_root: Path) -> dict[str, Any]:
    compound_id = str(config["compound_id"])
    display_name = str(config["display_name"])
    solvent = str(config["solvent"])
    nuclei = [str(item) for item in config.get("nuclei", ["1H", "13C"])]
    unsupported = sorted(set(nuclei) - set(ELEMENT_BY_NUCLEUS))
    if unsupported:
        raise PostprocessError(f"Unsupported nuclei: {unsupported}")
    reference = {str(key): float(value) for key, value in config["reference_shielding_ppm"].items()}
    missing_reference = sorted(set(nuclei) - set(reference))
    if missing_reference:
        raise PostprocessError(f"Missing calculated reference shieldings: {missing_reference}")

    conformers = list(config.get("conformers", []))
    if not conformers:
        conformers = [
            {
                "conformer_id": "conf_01",
                "weight": 1.0,
                "shielding_output": config["shielding_output"],
                "j_output": config.get("j_output"),
                "atom_map": config["atom_map"],
            }
        ]
    weights = np.asarray([float(row["weight"]) for row in conformers], dtype=float)
    if np.any(weights < 0) or float(weights.sum()) <= 0:
        raise PostprocessError("Conformer weights must be non-negative and sum to a positive value")
    supplied_weight_sum = float(weights.sum())
    weights /= weights.sum()

    atom_map: list[AtomRecord] | None = None
    shielding_arrays: list[np.ndarray] = []
    coupling_sets: list[list[CouplingRecord]] = []
    source_files: list[Path] = []
    for conformer in conformers:
        map_path = resolve_path(config_path, conformer["atom_map"])
        current_map = read_atom_map(map_path)
        if atom_map is None:
            atom_map = current_map
        elif [(row.element, row.label) for row in current_map] != [
            (row.element, row.label) for row in atom_map
        ]:
            raise PostprocessError("Conformer atom maps differ in element or atom-label order")
        shielding_path = resolve_path(config_path, conformer["shielding_output"])
        shielding_arrays.append(np.asarray(parse_shieldings(shielding_path, current_map), dtype=float))
        source_files.extend([map_path, shielding_path])
        j_value = conformer.get("j_output")
        if j_value:
            j_path = resolve_path(config_path, j_value)
            coupling_sets.append(parse_j_couplings(j_path))
            source_files.append(j_path)
        else:
            coupling_sets.append([])
    assert atom_map is not None
    shielding_stack = np.vstack(shielding_arrays)
    shielding_mean = np.average(shielding_stack, axis=0, weights=weights)
    shielding_sd = np.sqrt(np.average((shielding_stack - shielding_mean) ** 2, axis=0, weights=weights))

    averaged_j: list[CouplingRecord] = []
    all_pairs = sorted(
        {
            tuple(sorted((row.atom_i0, row.atom_j0)))
            for records in coupling_sets
            for row in records
        }
    )
    for atom_i0, atom_j0 in all_pairs:
        values = []
        for records in coupling_sets:
            lookup = {
                tuple(sorted((row.atom_i0, row.atom_j0))): row.j_hz for row in records
            }
            values.append(lookup.get((atom_i0, atom_j0), 0.0))
        averaged_j.append(
            CouplingRecord(atom_i0, atom_j0, float(np.average(values, weights=weights)))
        )

    compound_root = output_root / compound_id
    compound_root.mkdir(parents=True, exist_ok=True)
    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "generated_utc": utc_now(),
        "compound_id": compound_id,
        "display_name": display_name,
        "solvent": solvent,
        "status": CALCULATED,
        "conformer_count": len(conformers),
        "supplied_weight_sum": supplied_weight_sum,
        "normalized_weights": weights.tolist(),
        "nuclei": {},
        "source_files": [
            {"file": str(path), "sha256": sha256_file(path)} for path in sorted(set(source_files))
        ],
    }
    workbook_rows: list[dict[str, Any]] = []
    transition_workbook_rows: list[dict[str, Any]] = []
    exact_max_spins = int(config.get("exact_max_spins", 10))
    merge_hz = float(config.get("transition_merge_tolerance_hz", 0.02))
    profile_points = int(config.get("profile_points", 50000))
    dpi = int(config.get("dpi", 300))
    seed = int(config.get("appearance_seed", 20260720))
    make_appearance = bool(config.get("make_synthetic_appearance", True))

    for nucleus in nuclei:
        target = ELEMENT_BY_NUCLEUS[nucleus]
        shifts = {
            atom.index0: reference[nucleus] - float(shielding_mean[atom.index0])
            for atom in atom_map
            if atom.element == target and atom.include_in_analyte_nmr
        }
        transitions, spin_method = make_standard_transitions(
            nucleus,
            atom_map,
            shifts,
            averaged_j,
            float(config.get("frequencies_MHz", {}).get(nucleus, FREQUENCY_DEFAULTS[nucleus])),
            solvent,
            exact_max_spins,
            merge_hz,
        )
        nucleus_root = compound_root / "NMR" / nucleus
        nucleus_root.mkdir(parents=True, exist_ok=True)
        transition_csv = nucleus_root / "P0_calculated_sticks.csv"
        write_rows(transition_csv, transitions)
        write_rows(nucleus_root / "P0_calculated_sticks.txt", transitions, delimiter="\t")

        atom_rows: list[dict[str, Any]] = []
        for atom in atom_map:
            if atom.index0 not in shifts:
                continue
            row = {
                "compound_id": compound_id,
                "display_name": display_name,
                "solvent": solvent,
                "nucleus": nucleus,
                "atom_index_0based": atom.index0,
                "atom_index_1based": atom.index0 + 1,
                "element": atom.element,
                "atom_label": atom.label,
                "exchangeable": atom.exchangeable,
                "reference_shielding_ppm": f"{reference[nucleus]:.12g}",
                "weighted_shielding_ppm": f"{shielding_mean[atom.index0]:.12g}",
                "conformer_shielding_sd_ppm": f"{shielding_sd[atom.index0]:.12g}",
                "calculated_shift_ppm": f"{shifts[atom.index0]:.12g}",
                "data_status": CALCULATED,
            }
            atom_rows.append(row)
            workbook_rows.append(row)
        write_rows(nucleus_root / "atom_shifts.csv", atom_rows)
        write_rows(nucleus_root / "atom_shifts.txt", atom_rows, delimiter="\t")

        if nucleus == "1H" and averaged_j:
            j_rows = [
                {
                    "compound_id": compound_id,
                    "atom_i_0based": row.atom_i0,
                    "atom_j_0based": row.atom_j0,
                    "J_Hz": f"{row.j_hz:.12g}",
                    "data_status": CALCULATED,
                }
                for row in averaged_j
            ]
            write_rows(nucleus_root / "J_couplings_weighted.csv", j_rows)
            write_rows(nucleus_root / "J_couplings_weighted.txt", j_rows, delimiter="\t")

        frequency = float(config.get("frequencies_MHz", {}).get(nucleus, FREQUENCY_DEFAULTS[nucleus]))
        shift_values = [float(row["shift_ppm"]) for row in transitions]
        low, high = display_range(config, nucleus, shift_values)
        axis, signal = clean_profile(transitions, frequency, low, high, profile_points)
        profile_rows = [
            {
                "shift_ppm": f"{ppm:.9f}",
                "analyte_relative_intensity": f"{value:.12g}",
                "data_status": CALCULATED,
            }
            for ppm, value in zip(axis, signal)
        ]
        write_rows(nucleus_root / "P1_clean_calculated_profile.csv", profile_rows)
        write_rows(nucleus_root / "P1_clean_calculated_profile.txt", profile_rows, delimiter="\t")
        figure_paths = save_nmr_figure(
            nucleus_root / "P0_P1_calculated_model",
            compound_id,
            display_name,
            nucleus,
            solvent,
            transitions,
            axis,
            signal,
            dpi,
        )
        appearance_path = None
        if make_appearance:
            appearance_root_value = config.get("appearance_output_root")
            appearance_root = (
                resolve_path(config_path, appearance_root_value)
                if appearance_root_value
                else output_root.parent / "_P2"
            )
            appearance_path = invoke_appearance(
                transition_csv,
                compound_id,
                display_name,
                nucleus,
                solvent,
                appearance_root,
                seed,
                dpi,
            )
        for row in transitions:
            transition_workbook_rows.append(
                {
                    "compound_id": compound_id,
                    "display_name": display_name,
                    "solvent": solvent,
                    **row,
                }
            )
        nucleus_metadata = {
            "nucleus": nucleus,
            "frequency_MHz": frequency,
            "reference_shielding_ppm": reference[nucleus],
            "reference_policy": "calculated reference shielding supplied by config; no experimental peak fitting",
            "spin_simulation": spin_method,
            "observable_atom_count": len(atom_rows),
            "transition_count": len(transitions),
            "P0_notice": CALCULATED,
            "P1_notice": CALCULATED,
            "P2_notice": SYNTHETIC if appearance_path else "NOT_REQUESTED",
            "13C_height_warning": (
                "direct/deconvoluted display intensity is non-quantitative and is not a carbon count"
                if nucleus == "13C"
                else None
            ),
            "files": {
                "P0_csv": transition_csv.name,
                "P1_csv": "P1_clean_calculated_profile.csv",
                "figures": [path.name for path in figure_paths],
                "P2_directory": str(appearance_path) if appearance_path else None,
            },
        }
        write_json(nucleus_root / "nmr_postprocess_metadata.json", nucleus_metadata)
        summary["nuclei"][nucleus] = nucleus_metadata

    write_rows(compound_root / "workbook_ready_nmr_atom_shifts.csv", workbook_rows)
    write_rows(compound_root / "workbook_ready_nmr_transitions.csv", transition_workbook_rows)
    write_json(compound_root / "compound_nmr_summary.json", summary)
    return summary


def load_hrms_engine():
    return load_module("formula_hrms_postprocess_engine", SCRIPTS / "generate_formula_hrms.py")


def save_hrms_figure(
    stem: Path,
    compound_id: str,
    display_name: str,
    ion_label: str,
    ion_formula: str,
    monoisotopic_mz: float,
    sticks: Sequence[Any],
    axis: np.ndarray,
    signal: np.ndarray,
    resolving_power: float,
    dpi: int,
) -> list[Path]:
    fig, ax = plt.subplots(figsize=(7.20, 2.15))
    ax.plot(axis, signal, color=PALETTE["hrms"], lw=0.85)
    for stick in sticks:
        ax.vlines(stick.mz, 0, stick.relative_intensity, color=PALETTE["sticks"], lw=0.45)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Rel. intensity (%)")
    ax.set_ylim(0, 108)
    ax.set_title(f"{compound_id} | {display_name} | {ion_label}", loc="left", weight="bold")
    ax.text(
        0.995,
        0.95,
        f"{ion_formula}; monoisotopic m/z {monoisotopic_mz:.6f}; R = {resolving_power:g}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6,
        color=PALETTE["neutral"],
    )
    ax.text(
        0.995,
        0.06,
        "FORMULA-DERIVED ONLY — NO FRAGMENTATION PREDICTED",
        transform=ax.transAxes,
        ha="right",
        fontsize=6,
        color=PALETTE["accent"],
        weight="bold",
    )
    fig.tight_layout(pad=1.2)
    outputs: list[Path] = []
    for suffix, kwargs in ((".svg", {}), (".png", {"dpi": dpi})):
        path = stem.with_suffix(suffix)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, bbox_inches="tight", **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def process_hrms(config: Mapping[str, Any], output_root: Path) -> dict[str, Any] | None:
    ions = list(config.get("hrms_ions", []))
    if not ions:
        return None
    engine = load_hrms_engine()
    compound_id = str(config["compound_id"])
    display_name = str(config["display_name"])
    compound_root = output_root / compound_id
    hrms_root = compound_root / "HRMS"
    hrms_root.mkdir(parents=True, exist_ok=True)
    resolving_power = float(config.get("hrms_resolving_power", 60000.0))
    min_relative = float(config.get("hrms_min_relative_percent", 0.01))
    max_states = int(config.get("hrms_max_isotope_states", 20000))
    dpi = int(config.get("dpi", 300))
    summary_rows: list[dict[str, Any]] = []
    metadata: dict[str, Any] = {
        "schema_version": "1.0",
        "generated_utc": utc_now(),
        "compound_id": compound_id,
        "display_name": display_name,
        "scope": HRMS_SCOPE,
        "orca_used_for_hrms": False,
        "fragment_ions_generated": False,
        "ions": [],
    }
    for sequence, ion in enumerate(ions, 1):
        label = str(ion["label"])
        formula = str(ion["ion_formula"])
        charge = int(ion["charge"])
        if charge == 0:
            raise PostprocessError(f"HRMS ion {label} has zero charge")
        counts = engine.parse_formula(formula)
        envelope = engine.calculate_envelope(
            counts,
            charge=charge,
            max_states=max_states,
            min_relative_percent=min_relative,
        )
        axis, signal, step = engine.gaussian_profile(
            envelope.sticks,
            resolving_power=resolving_power,
            points_per_fwhm=12,
            max_points=250000,
        )
        slug = f"ion_{sequence:02d}"
        ion_root = hrms_root / slug
        ion_root.mkdir(parents=True, exist_ok=True)
        stick_rows = [
            {
                "compound_id": compound_id,
                "display_name": display_name,
                "ion_label": label,
                "ion_formula": engine.charged_formula(counts, charge),
                "charge": charge,
                "isotopologue_index": index,
                "calculated_mz": f"{stick.mz:.9f}",
                "natural_abundance_probability": f"{stick.probability:.12g}",
                "relative_intensity_percent": f"{stick.relative_intensity:.9f}",
                "mass_shift_from_monoisotopic_mz_u": f"{stick.mass_shift_u:.9f}",
                "data_status": CALCULATED,
            }
            for index, stick in enumerate(envelope.sticks, 1)
        ]
        profile_rows = [
            {
                "calculated_mz": f"{mz:.9f}",
                "calculated_relative_intensity_percent": f"{value:.9f}",
                "data_status": CALCULATED,
            }
            for mz, value in zip(axis, signal)
        ]
        write_rows(ion_root / "isotope_sticks.csv", stick_rows)
        write_rows(ion_root / "isotope_sticks.txt", stick_rows, delimiter="\t")
        write_rows(ion_root / "gaussian_profile.csv", profile_rows)
        write_rows(ion_root / "gaussian_profile.txt", profile_rows, delimiter="\t")
        figure_paths = save_hrms_figure(
            ion_root / "formula_derived_isotope_envelope",
            compound_id,
            display_name,
            label,
            engine.charged_formula(counts, charge),
            envelope.monoisotopic_mz,
            envelope.sticks,
            axis,
            signal,
            resolving_power,
            dpi,
        )
        summary_row = {
            "compound_id": compound_id,
            "display_name": display_name,
            "ion_label": label,
            "ion_formula": engine.charged_formula(counts, charge),
            "charge": charge,
            "monoisotopic_mz": f"{envelope.monoisotopic_mz:.9f}",
            "resolving_power_fwhm": resolving_power,
            "isotope_stick_count": len(envelope.sticks),
            "retained_probability_before_display_filter": f"{envelope.retained_probability_before_display_filter:.12g}",
            "retained_probability_after_display_filter": f"{envelope.retained_probability_after_display_filter:.12g}",
            "fragment_ions_generated": False,
            "retention_time_predicted": False,
            "orca_used_for_hrms": False,
            "scope": HRMS_SCOPE,
            "data_status": CALCULATED,
        }
        summary_rows.append(summary_row)
        metadata["ions"].append(
            {
                **summary_row,
                "profile_step_mz": step,
                "figures": [path.name for path in figure_paths],
            }
        )
        write_json(ion_root / "hrms_metadata.json", metadata["ions"][-1])
    write_rows(compound_root / "workbook_ready_hrms_rows.csv", summary_rows)
    write_json(hrms_root / "compound_hrms_summary.json", metadata)
    return metadata


def aggregate(output_root: Path) -> dict[str, Any]:
    nmr_rows: list[dict[str, str]] = []
    transition_rows: list[dict[str, str]] = []
    hrms_rows: list[dict[str, str]] = []
    for path in sorted(output_root.glob("*/workbook_ready_nmr_atom_shifts.csv")):
        with path.open(encoding="utf-8-sig", newline="") as handle:
            nmr_rows.extend(csv.DictReader(handle))
    for path in sorted(output_root.glob("*/workbook_ready_nmr_transitions.csv")):
        with path.open(encoding="utf-8-sig", newline="") as handle:
            transition_rows.extend(csv.DictReader(handle))
    for path in sorted(output_root.glob("*/workbook_ready_hrms_rows.csv")):
        with path.open(encoding="utf-8-sig", newline="") as handle:
            hrms_rows.extend(csv.DictReader(handle))
    aggregate_root = output_root / "_aggregate_workbook_ready"
    aggregate_root.mkdir(parents=True, exist_ok=True)
    if nmr_rows:
        write_rows(aggregate_root / "NMR_atom_shifts.csv", nmr_rows)
        write_rows(aggregate_root / "NMR_atom_shifts.txt", nmr_rows, delimiter="\t")
    if transition_rows:
        write_rows(aggregate_root / "NMR_transitions.csv", transition_rows)
        write_rows(aggregate_root / "NMR_transitions.txt", transition_rows, delimiter="\t")
    if hrms_rows:
        write_rows(aggregate_root / "HRMS_formula_ions.csv", hrms_rows)
        write_rows(aggregate_root / "HRMS_formula_ions.txt", hrms_rows, delimiter="\t")
    manifest = {
        "generated_utc": utc_now(),
        "status": CALCULATED,
        "nmr_atom_shift_rows": len(nmr_rows),
        "nmr_transition_rows": len(transition_rows),
        "hrms_ion_rows": len(hrms_rows),
        "workbook_policy": (
            "These normalized CSV/TXT tables are inputs for an eventual workbook; "
            "they are not experimental characterization."
        ),
        "hrms_scope": HRMS_SCOPE,
    }
    write_json(aggregate_root / "aggregate_manifest.json", manifest)
    return manifest


def run_config(config_path: Path, output_root: Path | None = None) -> dict[str, Any]:
    config_path = config_path.expanduser().resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if output_root is None:
        output_value = config.get("output_root", "results")
        output_root = resolve_path(config_path, output_value)
    else:
        output_root = output_root.expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    nmr_summary = process_nmr(config_path, config, output_root)
    hrms_summary = process_hrms(config, output_root)
    aggregate_summary = aggregate(output_root)
    run_summary = {
        "generated_utc": utc_now(),
        "config": str(config_path),
        "config_sha256": sha256_file(config_path),
        "output_root": str(output_root),
        "nmr": nmr_summary,
        "hrms": hrms_summary,
        "aggregate": aggregate_summary,
        "global_warning": (
            "All NMR P0/P1 values are calculated; P2 is synthetic and not acquired. "
            "HRMS is formula-derived and contains no LC-MS fragmentation or retention time."
        ),
    }
    write_json(output_root / str(config["compound_id"]) / "postprocess_run_summary.json", run_summary)
    return run_summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Post-process completed ORCA NMR properties and explicit formula HRMS ions. "
            "This program never launches ORCA and never creates experimental claims."
        )
    )
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Run NMR and optional formula-HRMS from one JSON config")
    run.add_argument("--config", type=Path, required=True)
    run.add_argument("--output-root", type=Path, default=None)
    agg = sub.add_parser("aggregate", help="Refresh workbook-ready aggregate CSV/TXT tables")
    agg.add_argument("--output-root", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "run":
        summary = run_config(args.config, args.output_root)
    else:
        summary = aggregate(args.output_root.expanduser().resolve())
    print(json.dumps(summary.get("aggregate", summary), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
