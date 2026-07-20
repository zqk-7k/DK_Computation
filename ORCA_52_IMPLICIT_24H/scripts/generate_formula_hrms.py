#!/usr/bin/env python3
"""Generate formula-level HRMS isotope predictions for the Figure 3 registry.

This program deliberately does *not* simulate collision-induced dissociation,
in-source fragmentation, contaminants, background ions, detector noise, or an
experimental baseline.  It converts the registered molecular formula/reference
state into conservative ESI ion formulas, convolves tabulated natural isotope
abundances, and applies a Gaussian resolving-power model.

Every generated artefact is labelled:

    CALCULATED PREDICTION — NOT EXPERIMENTAL

The default input is the ordered, machine-readable Figure 3 SDF and its matching
registry.  All 52 compounds, including the benzamide validation example, are
included by default.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.inchi import MolToInchiKey


PREDICTION_LABEL = "CALCULATED HRMS PREDICTION - NOT EXPERIMENTAL"
NO_FRAGMENT_NOTICE = (
    "Formula-level precursor/adduct isotope envelope only; no fragment ions, "
    "noise, contaminants, or experimental peaks were generated."
)
SKIP_IDS = frozenset()
ELECTRON_MASS_U = 0.000548579909065


@dataclass(frozen=True)
class Isotope:
    mass_number: int
    exact_mass: float
    abundance: float


# Exact isotope masses and representative terrestrial natural abundances.
# Values are embedded so the calculation does not depend on an external isotope
# package.  The elements cover every formula in the current 52-structure set and
# Na for the explicitly declared [M+Na]+ adduct.
ISOTOPES: dict[str, tuple[Isotope, ...]] = {
    "H": (
        Isotope(1, 1.00782503223, 0.999885),
        Isotope(2, 2.01410177812, 0.000115),
    ),
    "C": (
        Isotope(12, 12.00000000000, 0.9893),
        Isotope(13, 13.00335483507, 0.0107),
    ),
    "N": (
        Isotope(14, 14.00307400443, 0.99636),
        Isotope(15, 15.00010889888, 0.00364),
    ),
    "O": (
        Isotope(16, 15.99491461957, 0.99757),
        Isotope(17, 16.99913175650, 0.00038),
        Isotope(18, 17.99915961286, 0.00205),
    ),
    "F": (Isotope(19, 18.99840316273, 1.0),),
    "Na": (Isotope(23, 22.98976928200, 1.0),),
    "P": (Isotope(31, 30.97376199842, 1.0),),
    "S": (
        Isotope(32, 31.97207117440, 0.9499),
        Isotope(33, 32.97145890980, 0.0075),
        Isotope(34, 33.96786700400, 0.0425),
        Isotope(36, 35.96708071000, 0.0001),
    ),
    "Cl": (
        Isotope(35, 34.96885268200, 0.7576),
        Isotope(37, 36.96590260200, 0.2424),
    ),
    "I": (Isotope(127, 126.90447300000, 1.0),),
}


@dataclass(frozen=True)
class IonModel:
    label: str
    mode: str
    formula_counts: Mapping[str, int]
    charge: int
    rationale: str
    salt_handling: str


@dataclass(frozen=True)
class Stick:
    mz: float
    probability: float
    relative_intensity: float
    mass_shift_u: float


@dataclass(frozen=True)
class EnvelopeResult:
    sticks: tuple[Stick, ...]
    monoisotopic_mz: float
    retained_probability_before_display_filter: float
    retained_probability_after_display_filter: float


FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula: str) -> Counter[str]:
    """Parse a neutral Hill-style formula (charge is stored separately)."""

    compact = formula.strip().replace(" ", "")
    counts: Counter[str] = Counter()
    cursor = 0
    for match in FORMULA_TOKEN.finditer(compact):
        if match.start() != cursor:
            raise ValueError(f"Unsupported formula syntax: {formula!r}")
        element, count_text = match.groups()
        counts[element] += int(count_text or "1")
        cursor = match.end()
    if cursor != len(compact) or not counts:
        raise ValueError(f"Unsupported formula syntax: {formula!r}")
    unsupported = sorted(set(counts) - set(ISOTOPES))
    if unsupported:
        raise ValueError(
            f"No isotope table for {', '.join(unsupported)} in {formula!r}"
        )
    return counts


def hill_formula(counts: Mapping[str, int]) -> str:
    positive = {element: int(n) for element, n in counts.items() if int(n) > 0}
    if "C" in positive:
        order = ["C"]
        if "H" in positive:
            order.append("H")
        order.extend(sorted(set(positive) - set(order)))
    else:
        order = sorted(positive)
    return "".join(
        element + (str(positive[element]) if positive[element] != 1 else "")
        for element in order
    )


def charged_formula(counts: Mapping[str, int], charge: int) -> str:
    base = hill_formula(counts)
    if charge == 0:
        return base
    magnitude = "" if abs(charge) == 1 else str(abs(charge))
    sign = "+" if charge > 0 else "-"
    return f"{base}{magnitude}{sign}"


def formula_counts_from_mol(mol: Chem.Mol) -> Counter[str]:
    """Count all atoms after RDKit has materialized implicit hydrogens."""

    with_h = Chem.AddHs(Chem.Mol(mol), addCoords=False)
    counts = Counter(atom.GetSymbol() for atom in with_h.GetAtoms())
    unsupported = sorted(set(counts) - set(ISOTOPES))
    if unsupported:
        raise ValueError(f"Unsupported element(s) in molecule: {unsupported}")
    return counts


def add_to_formula(
    counts: Mapping[str, int], element: str, delta: int
) -> Counter[str]:
    updated = Counter({key: int(value) for key, value in counts.items()})
    updated[element] += delta
    if updated[element] < 0:
        raise ValueError(
            f"Cannot apply {delta:+d} {element} to {hill_formula(counts)}"
        )
    if updated[element] == 0:
        del updated[element]
    return updated


def monoisotopic_ionic_mass(counts: Mapping[str, int], charge: int) -> float:
    neutral_atom_sum = sum(
        int(number) * ISOTOPES[element][0].exact_mass
        for element, number in counts.items()
    )
    return neutral_atom_sum - charge * ELECTRON_MASS_U


def monoisotopic_mz(counts: Mapping[str, int], charge: int) -> float:
    if charge == 0:
        return monoisotopic_ionic_mass(counts, 0)
    return monoisotopic_ionic_mass(counts, charge) / abs(charge)


def _compositions(total: int, parts: int) -> Iterator[tuple[int, ...]]:
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for remainder in _compositions(total - first, parts - 1):
            yield (first,) + remainder


@lru_cache(maxsize=None)
def element_isotope_states(element: str, atom_count: int) -> tuple[tuple[float, float], ...]:
    """Return exact-mass/probability states for one element via multinomials."""

    if atom_count == 0:
        return ((0.0, 1.0),)
    isotopes = ISOTOPES[element]
    factorial_n = math.factorial(atom_count)
    states: list[tuple[float, float]] = []
    for allocation in _compositions(atom_count, len(isotopes)):
        coefficient = factorial_n
        probability = 1.0
        mass = 0.0
        for isotope_count, isotope in zip(allocation, isotopes):
            coefficient //= math.factorial(isotope_count)
            probability *= isotope.abundance**isotope_count
            mass += isotope.exact_mass * isotope_count
        probability *= coefficient
        if probability > 0.0:
            states.append((mass, probability))
    return tuple(states)


def _merge_near_identical_states(
    states: Iterable[tuple[float, float]], digits: int = 8
) -> list[tuple[float, float]]:
    """Merge mathematically equivalent paths without losing their mean mass."""

    probability_by_bin: defaultdict[float, float] = defaultdict(float)
    weighted_mass_by_bin: defaultdict[float, float] = defaultdict(float)
    for mass, probability in states:
        key = round(mass, digits)
        probability_by_bin[key] += probability
        weighted_mass_by_bin[key] += mass * probability
    return [
        (weighted_mass_by_bin[key] / probability, probability)
        for key, probability in probability_by_bin.items()
        if probability > 0.0
    ]


def isotope_distribution(
    counts: Mapping[str, int], max_states: int
) -> tuple[list[tuple[float, float]], float]:
    """Self-contained fine-isotope convolution with probability pruning."""

    component_patterns = [
        element_isotope_states(element, int(number))
        for element, number in counts.items()
        if int(number) > 0
    ]
    component_patterns.sort(key=len)
    states: list[tuple[float, float]] = [(0.0, 1.0)]
    for component in component_patterns:
        expanded = (
            (mass_a + mass_b, probability_a * probability_b)
            for mass_a, probability_a in states
            for mass_b, probability_b in component
        )
        states = _merge_near_identical_states(expanded)
        if len(states) > max_states:
            states.sort(key=lambda item: item[1], reverse=True)
            states = states[:max_states]
    states.sort(key=lambda item: item[0])
    return states, sum(probability for _, probability in states)


def calculate_envelope(
    counts: Mapping[str, int],
    charge: int,
    max_states: int,
    min_relative_percent: float,
) -> EnvelopeResult:
    if charge == 0:
        raise ValueError("An HRMS ion envelope requires a non-zero charge")
    states, retained_before = isotope_distribution(counts, max_states=max_states)
    if not states:
        raise RuntimeError("Isotope convolution returned no states")
    maximum_probability = max(probability for _, probability in states)
    mono_neutral_atom_mass = sum(
        int(number) * ISOTOPES[element][0].exact_mass
        for element, number in counts.items()
    )
    sticks: list[Stick] = []
    for neutral_atom_mass, probability in states:
        relative = 100.0 * probability / maximum_probability
        if relative + 1e-15 < min_relative_percent:
            continue
        mz = (neutral_atom_mass - charge * ELECTRON_MASS_U) / abs(charge)
        sticks.append(
            Stick(
                mz=mz,
                probability=probability,
                relative_intensity=relative,
                mass_shift_u=(neutral_atom_mass - mono_neutral_atom_mass)
                / abs(charge),
            )
        )
    sticks.sort(key=lambda item: item.mz)
    return EnvelopeResult(
        sticks=tuple(sticks),
        monoisotopic_mz=monoisotopic_mz(counts, charge),
        retained_probability_before_display_filter=retained_before,
        retained_probability_after_display_filter=sum(
            item.probability for item in sticks
        ),
    )


def gaussian_profile(
    sticks: Sequence[Stick],
    resolving_power: float,
    points_per_fwhm: int,
    max_points: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Apply an area-normalized Gaussian R=m/delta-m response."""

    if not sticks:
        raise ValueError("Cannot build a profile without isotope sticks")
    sigma_factor = 2.0 * math.sqrt(2.0 * math.log(2.0))
    sigmas = np.array(
        [stick.mz / resolving_power / sigma_factor for stick in sticks],
        dtype=float,
    )
    minimum_sigma = float(sigmas.min())
    maximum_sigma = float(sigmas.max())
    start = min(stick.mz for stick in sticks) - 6.0 * maximum_sigma
    stop = max(stick.mz for stick in sticks) + 6.0 * maximum_sigma
    target_step = minimum_sigma * sigma_factor / points_per_fwhm
    point_count = int(math.ceil((stop - start) / target_step)) + 1
    point_count = max(3, min(point_count, max_points))
    x = np.linspace(start, stop, point_count, dtype=float)
    actual_step = float(x[1] - x[0])
    y = np.zeros_like(x)
    root_two_pi = math.sqrt(2.0 * math.pi)
    for stick, sigma in zip(sticks, sigmas):
        lower = int(np.searchsorted(x, stick.mz - 6.0 * sigma, side="left"))
        upper = int(np.searchsorted(x, stick.mz + 6.0 * sigma, side="right"))
        local_x = x[lower:upper]
        if local_x.size == 0:
            continue
        y[lower:upper] += (
            stick.probability
            / (sigma * root_two_pi)
            * np.exp(-0.5 * ((local_x - stick.mz) / sigma) ** 2)
        )
    if float(y.max()) > 0.0:
        y = 100.0 * y / float(y.max())
    return x, y, actual_step


def _fragment_is_counterion(fragment: Chem.Mol) -> bool:
    atoms = list(fragment.GetAtoms())
    if not atoms or len(atoms) > 3:
        return False
    allowed = {"F", "Cl", "Br", "I", "Li", "Na", "K"}
    return all(atom.GetSymbol() in allowed for atom in atoms) and (
        Chem.GetFormalCharge(fragment) != 0
    )


def _combine_fragments(fragments: Sequence[Chem.Mol]) -> Chem.Mol:
    if not fragments:
        raise ValueError("No analyte fragment remains after counterion removal")
    combined = Chem.Mol(fragments[0])
    for fragment in fragments[1:]:
        combined = Chem.CombineMols(combined, fragment)
    Chem.SanitizeMol(combined)
    return combined


ACIDIC_PATTERNS = tuple(
    Chem.MolFromSmarts(smarts)
    for smarts in (
        "[O;H1]-[C,S,P](=[O,S])",
        "[O;H1]-[c]",
        "[N;H1,H2]-[S](=O)(=O)",
        "[N;H1,H2]-[P](=O)",
        "[nH]",
    )
)


def _has_acidic_site(mol: Chem.Mol) -> bool:
    return any(
        pattern is not None and mol.HasSubstructMatch(pattern)
        for pattern in ACIDIC_PATTERNS
    )


def _charge_suffix(charge: int) -> str:
    if charge == 0:
        return ""
    magnitude = "" if abs(charge) == 1 else str(abs(charge))
    return magnitude + ("+" if charge > 0 else "-")


def _counterion_loss_text(counterion_counts: Mapping[str, int]) -> str:
    if len(counterion_counts) == 1:
        element, number = next(iter(counterion_counts.items()))
        return (str(number) if number != 1 else "") + element
    return hill_formula(counterion_counts)


def select_ion_models(mol: Chem.Mol) -> list[IonModel]:
    """Select conservative precursor/adduct ions without fragmentation."""

    fragments = list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True))
    counterions = [fragment for fragment in fragments if _fragment_is_counterion(fragment)]
    analyte_fragments = [
        fragment for fragment in fragments if not _fragment_is_counterion(fragment)
    ]
    if counterions:
        analyte = _combine_fragments(analyte_fragments)
        analyte_charge = Chem.GetFormalCharge(analyte)
        analyte_counts = formula_counts_from_mol(analyte)
        counterion_counts: Counter[str] = Counter()
        for counterion in counterions:
            counterion_counts.update(formula_counts_from_mol(counterion))
        if analyte_charge != 0:
            loss_text = _counterion_loss_text(counterion_counts)
            label = f"[M-{loss_text}]{_charge_suffix(analyte_charge)}"
            return [
                IonModel(
                    label=label,
                    mode="ESI+" if analyte_charge > 0 else "ESI-",
                    formula_counts=analyte_counts,
                    charge=analyte_charge,
                    rationale=(
                        "Pre-existing analyte ion after excluding explicit "
                        "counterion(s); this is a salt/desolvation model, not a "
                        "fragmentation prediction."
                    ),
                    salt_handling=(
                        f"Explicit counterion formula {hill_formula(counterion_counts)} "
                        "excluded from the measured analyte ion."
                    ),
                )
            ]
        mol = analyte

    base_counts = formula_counts_from_mol(mol)
    net_charge = Chem.GetFormalCharge(mol)
    if net_charge != 0:
        return [
            IonModel(
                label=f"[M]{_charge_suffix(net_charge)}",
                mode="ESI+" if net_charge > 0 else "ESI-",
                formula_counts=base_counts,
                charge=net_charge,
                rationale="Pre-existing registered molecular ion.",
                salt_handling="No explicit counterion was included in this ion formula.",
            )
        ]

    has_positive_atom = any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms())
    has_negative_atom = any(atom.GetFormalCharge() < 0 for atom in mol.GetAtoms())
    zwitterion = has_positive_atom and has_negative_atom
    models = [
        IonModel(
            label="[M+H]+",
            mode="ESI+",
            formula_counts=add_to_formula(base_counts, "H", +1),
            charge=+1,
            rationale="Protonated intact precursor/adduct ion.",
            salt_handling="No explicit counterion in the registered reference state.",
        )
    ]
    if not zwitterion and (base_counts.get("O", 0) + base_counts.get("N", 0) > 0):
        models.append(
            IonModel(
                label="[M+Na]+",
                mode="ESI+",
                formula_counts=add_to_formula(base_counts, "Na", +1),
                charge=+1,
                rationale="Sodiated intact adduct ion; no bond fragmentation applied.",
                salt_handling="No explicit counterion in the registered reference state.",
            )
        )
    if zwitterion or _has_acidic_site(mol):
        models.append(
            IonModel(
                label="[M-H]-",
                mode="ESI-",
                formula_counts=add_to_formula(base_counts, "H", -1),
                charge=-1,
                rationale="Deprotonated intact precursor/adduct ion.",
                salt_handling="No explicit counterion in the registered reference state.",
            )
        )
    return models


def safe_slug(text: str) -> str:
    replacements = {
        "+": "_plus_",
        "-": "_minus_",
        "[": "",
        "]": "",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    text = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return text or "ion"


def metadata_comment_lines(metadata: Mapping[str, object]) -> list[str]:
    lines = [
        f"# {PREDICTION_LABEL}",
        f"# {NO_FRAGMENT_NOTICE}",
    ]
    for key, value in metadata.items():
        lines.append(f"# {key}: {value}")
    return lines


def write_csv(
    path: Path,
    metadata: Mapping[str, object],
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        for line in metadata_comment_lines(metadata):
            handle.write(line + "\n")
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_profile(
    path: Path,
    figure_id: str,
    display_name: str,
    ion: IonModel,
    envelope: EnvelopeResult,
    x: np.ndarray,
    y: np.ndarray,
    resolving_power: float,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
        }
    )
    figure, axis = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    axis.plot(x, y, color="#202020", linewidth=1.0, label=f"Gaussian profile, R={resolving_power:g}")
    stick_scale = np.array([stick.relative_intensity for stick in envelope.sticks])
    axis.vlines(
        [stick.mz for stick in envelope.sticks],
        0,
        stick_scale,
        color="#2574A9",
        linewidth=0.75,
        alpha=0.80,
        label="Natural-abundance isotope sticks",
    )
    axis.set_ylim(0, 108)
    axis.set_xlabel("Calculated m/z")
    axis.set_ylabel("Calculated relative intensity (%)")
    axis.set_title(
        f"{PREDICTION_LABEL}\n"
        f"{figure_id} | {display_name} | {ion.label} | "
        f"{charged_formula(ion.formula_counts, ion.charge)}",
        fontsize=9.5,
        fontweight="bold",
    )
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.legend(frameon=False, fontsize=7.5, loc="upper right")
    axis.text(
        0.0,
        -0.22,
        "Formula isotope envelope only. No fragment ions or experimental peaks.",
        transform=axis.transAxes,
        fontsize=7.2,
        color="#555555",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        path,
        dpi=300,
        metadata={
            "Title": PREDICTION_LABEL,
            "Description": NO_FRAGMENT_NOTICE,
            "Author": "Formula-level Python calculation",
        },
    )
    plt.close(figure)


def load_registry(path: Path) -> tuple[list[dict[str, str]], dict[str, dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Registry is empty: {path}")
    by_id = {row["figure_id"]: row for row in rows}
    if len(by_id) != len(rows):
        raise ValueError("Duplicate figure_id values in registry")
    return rows, by_id


def load_ordered_sdf(path: Path) -> list[tuple[str, Chem.Mol]]:
    supplier = Chem.SDMolSupplier(
        str(path), removeHs=False, sanitize=True, strictParsing=False
    )
    records: list[tuple[str, Chem.Mol]] = []
    for index, mol in enumerate(supplier, start=1):
        if mol is None:
            raise ValueError(f"SDF record {index} failed RDKit parsing")
        if mol.HasProp("Figure3_ID"):
            figure_id = mol.GetProp("Figure3_ID")
        else:
            figure_id = f"F3_{index:03d}"
        records.append((figure_id, mol))
    if not records:
        raise ValueError(f"SDF is empty: {path}")
    return records


def validate_record(
    figure_id: str, mol: Chem.Mol, registry_row: Mapping[str, str]
) -> Counter[str]:
    if registry_row["figure_id"] != figure_id:
        raise ValueError(f"SDF/registry ID mismatch at {figure_id}")
    calculated_key = MolToInchiKey(mol)
    if calculated_key != registry_row["inchikey"]:
        raise ValueError(
            f"{figure_id}: InChIKey mismatch: {calculated_key} != "
            f"{registry_row['inchikey']}"
        )
    counts = formula_counts_from_mol(mol)
    calculated_formula = hill_formula(counts)
    if calculated_formula != registry_row["formula"]:
        raise ValueError(
            f"{figure_id}: formula mismatch: {calculated_formula} != "
            f"{registry_row['formula']}"
        )
    calculated_charge = Chem.GetFormalCharge(mol)
    if calculated_charge != int(registry_row["net_charge_as_drawn"]):
        raise ValueError(
            f"{figure_id}: charge mismatch: {calculated_charge} != "
            f"{registry_row['net_charge_as_drawn']}"
        )
    return counts


def process_compound(
    output_root: Path,
    figure_id: str,
    mol: Chem.Mol,
    registry_row: Mapping[str, str],
    resolving_power: float,
    min_relative_percent: float,
    max_states: int,
    points_per_fwhm: int,
    max_profile_points: int,
) -> list[dict[str, object]]:
    reference_counts = validate_record(figure_id, mol, registry_row)
    display_name = registry_row["display_name"]
    reference_charge = int(registry_row["net_charge_as_drawn"])
    reference_exact_mass = monoisotopic_ionic_mass(
        reference_counts, reference_charge
    )
    # The project path is already long on Windows.  Keep the on-disk compound
    # folder keyed only by the stable Figure ID; the full display name remains
    # in every CSV/JSON/PNG metadata record.
    folder = output_root / figure_id
    folder.mkdir(parents=True, exist_ok=True)
    ion_models = select_ion_models(mol)
    compound_metadata: dict[str, object] = {
        "prediction_status": PREDICTION_LABEL,
        "figure_id": figure_id,
        "compound_number": registry_row["compound_number"],
        "display_name": display_name,
        "registered_reference_formula": registry_row["formula"],
        "registered_net_charge": reference_charge,
        "registered_inchikey": registry_row["inchikey"],
        "registered_reference_exact_mass_u": reference_exact_mass,
        "calculation_scope": NO_FRAGMENT_NOTICE,
        "fragment_ions_generated": False,
        "natural_abundance_model": "Embedded exact-isotope multinomial convolution",
        "resolving_power_fwhm": resolving_power,
        "display_stick_threshold_percent_of_base_peak": min_relative_percent,
        "ion_models": [],
    }
    summary_rows: list[dict[str, object]] = []
    for ion in ion_models:
        envelope = calculate_envelope(
            ion.formula_counts,
            charge=ion.charge,
            max_states=max_states,
            min_relative_percent=min_relative_percent,
        )
        x, y, profile_step = gaussian_profile(
            envelope.sticks,
            resolving_power=resolving_power,
            points_per_fwhm=points_per_fwhm,
            max_points=max_profile_points,
        )
        ion_slug = safe_slug(ion.label)
        common_metadata = {
            "prediction_status": PREDICTION_LABEL,
            "figure_id": figure_id,
            "display_name": display_name,
            "ion_label": ion.label,
            "ion_mode": ion.mode,
            "ion_formula": charged_formula(ion.formula_counts, ion.charge),
            "ion_charge": ion.charge,
            "monoisotopic_mz": f"{envelope.monoisotopic_mz:.8f}",
            "resolving_power_fwhm": resolving_power,
            "fragment_ions_generated": False,
        }
        stick_rows = [
            {
                "isotopologue_index": index,
                "calculated_mz": f"{stick.mz:.8f}",
                "natural_abundance_probability": f"{stick.probability:.12g}",
                "relative_intensity_percent": f"{stick.relative_intensity:.8f}",
                "mass_shift_from_monoisotopic_mz_u": f"{stick.mass_shift_u:.8f}",
            }
            for index, stick in enumerate(envelope.sticks, start=1)
        ]
        write_csv(
            folder / f"{figure_id}_{ion_slug}_isotope_sticks.csv",
            common_metadata,
            (
                "isotopologue_index",
                "calculated_mz",
                "natural_abundance_probability",
                "relative_intensity_percent",
                "mass_shift_from_monoisotopic_mz_u",
            ),
            stick_rows,
        )
        profile_rows = (
            {
                "calculated_mz": f"{mz:.8f}",
                "calculated_relative_intensity_percent": f"{intensity:.8f}",
            }
            for mz, intensity in zip(x, y)
        )
        write_csv(
            folder
            / f"{figure_id}_{ion_slug}_gaussian_profile_R{resolving_power:g}.csv",
            {**common_metadata, "profile_grid_step_u": f"{profile_step:.10g}"},
            ("calculated_mz", "calculated_relative_intensity_percent"),
            profile_rows,
        )
        plot_profile(
            folder
            / f"{figure_id}_{ion_slug}_gaussian_profile_R{resolving_power:g}.png",
            figure_id,
            display_name,
            ion,
            envelope,
            x,
            y,
            resolving_power,
        )
        ion_metadata = {
            **asdict(ion),
            "formula_counts": dict(ion.formula_counts),
            "ion_formula": charged_formula(ion.formula_counts, ion.charge),
            "monoisotopic_mz": envelope.monoisotopic_mz,
            "isotope_stick_count": len(envelope.sticks),
            "retained_probability_before_display_filter": (
                envelope.retained_probability_before_display_filter
            ),
            "retained_probability_after_display_filter": (
                envelope.retained_probability_after_display_filter
            ),
            "profile_grid_step_u": profile_step,
            "fragment_ions_generated": False,
        }
        compound_metadata["ion_models"].append(ion_metadata)
        summary_rows.append(
            {
                "prediction_status": PREDICTION_LABEL,
                "figure_id": figure_id,
                "compound_number": registry_row["compound_number"],
                "display_name": display_name,
                "reference_formula": registry_row["formula"],
                "reference_net_charge": reference_charge,
                "reference_exact_mass_u": f"{reference_exact_mass:.8f}",
                "ion_label": ion.label,
                "ion_mode": ion.mode,
                "ion_formula": charged_formula(ion.formula_counts, ion.charge),
                "ion_charge": ion.charge,
                "monoisotopic_mz": f"{envelope.monoisotopic_mz:.8f}",
                "isotope_stick_count": len(envelope.sticks),
                "resolving_power_fwhm": resolving_power,
                "salt_handling": ion.salt_handling,
                "fragment_ions_generated": False,
            }
        )
    with (folder / f"{figure_id}_calculation_metadata.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(compound_metadata, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return summary_rows


def build_parser(project_root: Path) -> argparse.ArgumentParser:
    default_base = project_root.parent
    parser = argparse.ArgumentParser(
        description=(
            "Generate formula-level natural-abundance HRMS predictions. "
            + PREDICTION_LABEL
        )
    )
    parser.add_argument(
        "--sdf",
        type=Path,
        default=default_base / "Figure3_52_ORCA_ordered_reference_states.sdf",
        help="Ordered Figure 3 SDF.",
    )
    parser.add_argument(
        "--registry",
        type=Path,
        default=default_base / "Figure3_52_ORCA_reference_state_registry.csv",
        help="Matching structure registry CSV.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=project_root / "HRMS_formula_predictions_CALCULATED_NOT_EXPERIMENTAL",
    )
    parser.add_argument(
        "--only",
        default="",
        help="Optional comma-separated Figure3 IDs for a limited run.",
    )
    parser.add_argument("--resolution", type=float, default=60000.0)
    parser.add_argument(
        "--min-relative-percent",
        type=float,
        default=0.01,
        help="Minimum isotope-stick intensity as percent of the base peak.",
    )
    parser.add_argument("--max-isotope-states", type=int, default=20000)
    parser.add_argument("--points-per-fwhm", type=int, default=12)
    parser.add_argument("--max-profile-points", type=int, default=250000)
    return parser


SUMMARY_FIELDS = (
    "prediction_status",
    "figure_id",
    "compound_number",
    "display_name",
    "reference_formula",
    "reference_net_charge",
    "reference_exact_mass_u",
    "ion_label",
    "ion_mode",
    "ion_formula",
    "ion_charge",
    "monoisotopic_mz",
    "isotope_stick_count",
    "resolving_power_fwhm",
    "salt_handling",
    "fragment_ions_generated",
)


def run(args: argparse.Namespace) -> int:
    if args.resolution <= 0:
        raise ValueError("--resolution must be positive")
    if not 0 <= args.min_relative_percent <= 100:
        raise ValueError("--min-relative-percent must be between 0 and 100")
    if args.max_isotope_states < 10:
        raise ValueError("--max-isotope-states must be at least 10")
    if args.points_per_fwhm < 3:
        raise ValueError("--points-per-fwhm must be at least 3")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    registry_rows, registry_by_id = load_registry(args.registry)
    sdf_records = load_ordered_sdf(args.sdf)
    registry_order = [row["figure_id"] for row in registry_rows]
    sdf_order = [figure_id for figure_id, _ in sdf_records]
    if sdf_order != registry_order:
        raise ValueError(
            "Ordered SDF and registry do not have the same Figure3_ID sequence"
        )

    only_ids = {
        value.strip() for value in args.only.split(",") if value.strip()
    }
    unknown_only = sorted(only_ids - set(registry_by_id))
    if unknown_only:
        raise ValueError(f"Unknown --only ID(s): {', '.join(unknown_only)}")

    selected = []
    skipped: list[dict[str, str]] = []
    for figure_id, mol in sdf_records:
        if figure_id in SKIP_IDS:
            skipped.append(
                {
                    "figure_id": figure_id,
                    "reason": "Always skipped: separately handled benzamide example.",
                }
            )
            continue
        if only_ids and figure_id not in only_ids:
            continue
        selected.append((figure_id, mol))
    if not selected:
        raise ValueError("No compounds selected after applying skip/only filters")

    summary_rows: list[dict[str, object]] = []
    for figure_id, mol in selected:
        print(f"{PREDICTION_LABEL}: processing {figure_id}")
        summary_rows.extend(
            process_compound(
                output_root=args.output_dir,
                figure_id=figure_id,
                mol=mol,
                registry_row=registry_by_id[figure_id],
                resolving_power=args.resolution,
                min_relative_percent=args.min_relative_percent,
                max_states=args.max_isotope_states,
                points_per_fwhm=args.points_per_fwhm,
                max_profile_points=args.max_profile_points,
            )
        )

    write_csv(
        args.output_dir / "Figure3_formula_HRMS_prediction_summary.csv",
        {
            "prediction_status": PREDICTION_LABEL,
            "calculation_scope": NO_FRAGMENT_NOTICE,
            "ordered_sdf": str(args.sdf.resolve()),
            "registry": str(args.registry.resolve()),
            "resolution_fwhm": args.resolution,
            "always_skipped": ",".join(sorted(SKIP_IDS)),
        },
        SUMMARY_FIELDS,
        summary_rows,
    )
    run_metadata = {
        "prediction_status": PREDICTION_LABEL,
        "calculation_scope": NO_FRAGMENT_NOTICE,
        "fragment_ions_generated": False,
        "input_sdf": str(args.sdf.resolve()),
        "input_registry": str(args.registry.resolve()),
        "output_directory": str(args.output_dir.resolve()),
        "processed_figure_ids": [figure_id for figure_id, _ in selected],
        "skipped": skipped,
        "resolving_power_fwhm": args.resolution,
        "minimum_relative_stick_percent": args.min_relative_percent,
        "max_isotope_states": args.max_isotope_states,
        "points_per_fwhm": args.points_per_fwhm,
        "isotope_table": {
            element: [asdict(isotope) for isotope in isotopes]
            for element, isotopes in ISOTOPES.items()
        },
    }
    with (args.output_dir / "RUN_METADATA.json").open("w", encoding="utf-8") as handle:
        json.dump(run_metadata, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(
        f"{PREDICTION_LABEL}: wrote {len(summary_rows)} ion envelope(s) "
        f"for {len(selected)} compound(s) to {args.output_dir}"
    )
    print(NO_FRAGMENT_NOTICE)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    RDLogger.DisableLog("rdApp.warning")
    project_root = Path(__file__).resolve().parents[1]
    parser = build_parser(project_root)
    args = parser.parse_args(argv)
    try:
        return run(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
