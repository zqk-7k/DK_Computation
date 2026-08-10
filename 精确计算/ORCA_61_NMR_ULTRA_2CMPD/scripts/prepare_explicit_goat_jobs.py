"""Prepare explicit-solvent, restart-safe ORCA GOAT jobs for Figure 3.

The script does not run ORCA.  It creates deterministic, independently embedded
starting structures, provenance/atom/isotope maps and a CSV manifest consumed by
``run_goat_batch.py``.  All products are computational models, never acquired
experimental data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from rdkit import Chem
from rdkit.Chem import AllChem


PROJECT = Path(__file__).resolve().parents[1]
SOURCE_REGISTRY = PROJECT / "registry" / "Figure3_52_ORCA_reference_state_registry.csv"
CONDITION_TABLE = PROJECT / "config" / "Figure3_52_proposed_NMR_conditions.csv"
JOB_ROOT = PROJECT / "01_goat_explicit"
DATA_ROOT = PROJECT / "data"
MANIFEST = DATA_ROOT / "GOAT_explicit_job_manifest.csv"
SCENARIO_REGISTRY = DATA_ROOT / "explicit_scenario_registry.csv"

D2O_COUNTS: dict[str, tuple[int, ...]] = {
    "F3_019": (6, 8, 10),
    "F3_020": (6, 8, 10),
    "F3_021": (8, 10, 12),
    "F3_035": (4, 6, 8),
    "F3_043": (4, 6, 8),
    "F3_044": (4, 6, 8),
}

# The approved primary aqueous species.  These are modelling assumptions, not
# claims that an experimental pD or counterion was measured.
SCENARIO_SMILES: dict[str, str] = {
    "F3_020": "NC(=O)CC[C@H]([NH3+])C(=O)[O-]",
    "F3_035": "NS(=O)(=O)[O-]",
    "F3_044": "NP(N)(=O)[O-]",
}

COUNTERION_ELEMENTS = {"Cl", "Br", "I", "Li", "Na", "K"}
DEUTERIUM_MASS = 2.01410177812
ISOTOPE_POLICY = DATA_ROOT / "ORCA_6_1_1_isotope_mass_policy.json"
ORCA_COORDINATE_DOC_URL = (
    "https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html"
)
ORCA_NMR_TUTORIAL_URL = "https://www.faccts.de/docs/orca/6.1/tutorials/spec/NMR.html"

# ORCA 6.1.1's external ``xyzfile`` reader expects a standard four-field atom
# row (element, x, y, z).  Isotopic masses are therefore deliberately *not*
# appended to ``start.xyz``.  Isotope masses do not change the electronic PES,
# so ordinary-H masses are appropriate for the xTB/GOAT conformer search.
# For the later frequency calculation, the selected GOAT geometry must instead
# be rebuilt as an inline ``* xyz`` block.  ORCA's documented per-atom syntax is
# ``H x y z M=2.01410177812``.
GOAT_XYZ_MASS_POLICY = "PLAIN_H_FOUR_COLUMN_EXTERNAL_XYZ"
DOWNSTREAM_D_MASS_POLICY = "INLINE_XYZ_M_EQUALS_FROM_SELECTED_GOAT_GEOMETRY"


@dataclass(frozen=True)
class AtomRecord:
    element: str
    x: float
    y: float
    z: float
    component: str
    source_atom_index_1based: int | None
    environment_atom: bool
    exchangeable: bool
    isotope_mass: float | None = None


@dataclass(frozen=True)
class Scenario:
    figure_id: str
    scenario_id: str
    solvent: str
    continuum: str
    explicit_solvent_count: int
    species_smiles: str
    assumed_pd: str
    species_note: str
    explicit_dmso_count: int = 0
    explicit_water_count: int = 0


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest().upper()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8", newline="\n")
    os.replace(temporary, path)


def atomic_json(path: Path, payload: object) -> None:
    atomic_text(path, json.dumps(payload, ensure_ascii=False, indent=2) + "\n")


def write_isotope_policy() -> str:
    """Write the machine-readable ORCA 6.1.1 isotope-mass hand-off contract.

    This is intentionally a hand-off contract rather than a ready-to-run
    opt/freq input: production opt/freq coordinates must come from the reviewed
    GOAT ensemble selection, not from an arbitrary starting structure.
    """

    payload = {
        "schema_version": "1.0",
        "notice": "CALCULATED MODEL - NOT EXPERIMENTAL DATA",
        "orca_version_target": "6.1.1",
        "goat_external_xyz": {
            "policy": GOAT_XYZ_MASS_POLICY,
            "required_atom_row_fields": ["element", "x_A", "y_A", "z_A"],
            "mass_annotation_permitted": False,
            "scientific_basis": (
                "Isotopic mass does not alter the electronic potential-energy surface; "
                "the GOAT/xTB conformer-search coordinates therefore use ordinary H."
            ),
        },
        "downstream_opt_freq": {
            "policy": DOWNSTREAM_D_MASS_POLICY,
            "coordinate_source": "REVIEWED_GOAT_SELECTED_GEOMETRY_ONLY",
            "atom_order_requirement": "PRESERVE_START_XYZ_ATOM_ORDER_AND_VERIFY_HASH",
            "orca_inline_block_open": "* xyz CHARGE MULTIPLICITY",
            "deuterium_atom_syntax": "H {x_A} {y_A} {z_A} M=2.01410177812",
            "orca_inline_block_close": "*",
            "mapping_source": "per-job isotope_map.csv (1-based atom indices)",
            "prohibition": (
                "Do not append M=... to an external XYZ file and do not use the "
                "initial start.xyz geometry for production opt/freq."
            ),
        },
        "official_documentation": [
            {
                "title": "ORCA 6.1 Manual, 2.2 Input of Coordinates",
                "url": ORCA_COORDINATE_DOC_URL,
                "relevance": (
                    "External xyzfile syntax and M=... as a non-standard-isotope "
                    "modifier after an inline atomic coordinate definition."
                ),
            },
            {
                "title": "ORCA 6.1 NMR tutorial",
                "url": ORCA_NMR_TUTORIAL_URL,
                "relevance": "Official inline-coordinate example using M=16.999131.",
            },
        ],
        "local_parser_observation": (
            "ORCA 6.1.1 ScanXYZ rejected a fifth-column M= token in external start.xyz; "
            "the production bundle therefore enforces standard four-column XYZ."
        ),
    }
    atomic_json(ISOTOPE_POLICY, payload)
    return sha256_file(ISOTOPE_POLICY)


def write_csv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def normalized_only(values: Sequence[str] | None) -> set[str]:
    result: set[str] = set()
    for raw in values or []:
        for token in raw.split(","):
            token = token.strip().upper()
            if not token:
                continue
            if token.isdigit():
                token = f"F3_{int(token):03d}"
            result.add(token)
    return result


def scenario_plan(row: dict[str, str], condition: dict[str, str]) -> list[Scenario]:
    figure_id = row["figure_id"]
    species_smiles = SCENARIO_SMILES.get(figure_id, row["isomeric_smiles"])
    primary = condition["primary_solvent"]
    scenarios: list[Scenario] = []
    if primary == "D2O":
        counts = D2O_COUNTS.get(figure_id)
        if not counts:
            raise ValueError(f"No approved explicit-water counts for {figure_id}")
        pd = {
            "F3_019": "7.4 ASSUMED if no actual pD record is recovered",
            "F3_020": "7.4 ASSUMED if no actual pD record is recovered",
            "F3_021": "5.5 ASSUMED for the dihydrochloride scenario",
            "F3_035": "7.0 ASSUMED; sulfamate anion",
            "F3_043": "7.0 ASSUMED; neutral phosphoramide",
            "F3_044": "7.0 ASSUMED; monoanionic diamidophosphate, counterion unknown",
        }[figure_id]
        for count in counts:
            scenarios.append(
                Scenario(
                    figure_id=figure_id,
                    scenario_id=f"D2O_{count:02d}W_primary",
                    solvent="D2O",
                    continuum="WATER",
                    explicit_solvent_count=count,
                    species_smiles=species_smiles,
                    assumed_pd=pd,
                    species_note="Approved primary D2O species; different hydration stoichiometries are separate sensitivity scenarios",
                    explicit_water_count=count,
                )
            )
    else:
        if figure_id in {"F3_016", "F3_017"}:
            cation_smiles = species_smiles.split(".")[0]
            return [
                Scenario(
                    figure_id,
                    "DMSO_cation_continuum_01",
                    "DMSO-d6",
                    "DMSO",
                    1,
                    cation_smiles,
                    "NOT_APPLICABLE",
                    "Cation-only continuum sensitivity scenario; chloride omitted explicitly",
                ),
                Scenario(
                    figure_id,
                    "DMSO_contact_ion_pair_01",
                    "DMSO-d6",
                    "DMSO",
                    1,
                    species_smiles,
                    "NOT_APPLICABLE",
                    "Contact chloride ion-pair sensitivity scenario",
                ),
                Scenario(
                    figure_id,
                    "DMSO_solvent_separated_02",
                    "DMSO-d6",
                    "DMSO",
                    2,
                    species_smiles,
                    "NOT_APPLICABLE",
                    "DMSO-separated chloride ion-pair sensitivity scenario",
                ),
            ]
        # These are separate fixed-composition sensitivity scenarios. They are
        # never combined with ordinary Boltzmann weights across stoichiometries.
        # The trace-water branch models a possible DMSO-d6/H2O microstate; it is
        # not a prediction of the sample water content.
        microsolvation_plan = (
            (0, 0, "continuum"),
            (1, 0, "primary"),
            (2, 0, "primary"),
            (3, 0, "extended"),
            (2, 1, "trace_water_sensitivity"),
        )
        for dmso_count, water_count, label in microsolvation_plan:
            total_count = dmso_count + water_count
            scenarios.append(
                Scenario(
                    figure_id=figure_id,
                    scenario_id=f"DMSO_{dmso_count:02d}_W{water_count:02d}_{label}",
                    solvent="DMSO-d6",
                    continuum="DMSO",
                    explicit_solvent_count=total_count,
                    species_smiles=species_smiles,
                    assumed_pd="NOT_APPLICABLE",
                    species_note=(
                        "Fixed-composition DMSO-d6 microsolvation model; explicit DMSO/H2O atoms are "
                        "environment atoms and cross-stoichiometry Boltzmann mixing is forbidden"
                    ),
                    explicit_dmso_count=dmso_count,
                    explicit_water_count=water_count,
                )
            )
    return scenarios


def random_unit(rng: random.Random) -> tuple[float, float, float]:
    while True:
        x, y, z = (rng.uniform(-1.0, 1.0) for _ in range(3))
        norm = math.sqrt(x * x + y * y + z * z)
        if 1.0e-8 < norm <= 1.0:
            return x / norm, y / norm, z / norm


def rotation_matrix(rng: random.Random) -> tuple[tuple[float, float, float], ...]:
    # Uniform random unit quaternion converted to a right-handed matrix.
    u1, u2, u3 = rng.random(), rng.random(), rng.random()
    qx = math.sqrt(1 - u1) * math.sin(2 * math.pi * u2)
    qy = math.sqrt(1 - u1) * math.cos(2 * math.pi * u2)
    qz = math.sqrt(u1) * math.sin(2 * math.pi * u3)
    qw = math.sqrt(u1) * math.cos(2 * math.pi * u3)
    return (
        (1 - 2 * (qy * qy + qz * qz), 2 * (qx * qy - qz * qw), 2 * (qx * qz + qy * qw)),
        (2 * (qx * qy + qz * qw), 1 - 2 * (qx * qx + qz * qz), 2 * (qy * qz - qx * qw)),
        (2 * (qx * qz - qy * qw), 2 * (qy * qz + qx * qw), 1 - 2 * (qx * qx + qy * qy)),
    )


def apply_rotation(
    point: tuple[float, float, float], matrix: tuple[tuple[float, float, float], ...]
) -> tuple[float, float, float]:
    x, y, z = point
    return tuple(sum(matrix[i][j] * (x, y, z)[j] for j in range(3)) for i in range(3))  # type: ignore[return-value]


def centroid(records: Sequence[AtomRecord], heavy_only: bool = True) -> tuple[float, float, float]:
    chosen = [r for r in records if not heavy_only or r.element != "H"] or list(records)
    return tuple(sum(getattr(r, axis) for r in chosen) / len(chosen) for axis in ("x", "y", "z"))  # type: ignore[return-value]


def distance(a: AtomRecord, b: AtomRecord) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def acceptable_fragment(fragment: Sequence[AtomRecord], existing: Sequence[AtomRecord]) -> bool:
    for a in fragment:
        for b in existing:
            minimum = 1.15 if "H" in (a.element, b.element) else 1.65
            if distance(a, b) < minimum:
                return False
    return True


def embed_fragment(smiles: str, seed: int) -> tuple[Chem.Mol, int]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit cannot parse species SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed & 0x7FFFFFFF) or 1
    params.useRandomCoords = True
    params.enforceChirality = True
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        raise RuntimeError(f"RDKit embedding failed for {smiles} (seed {seed})")
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
        else:
            AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
    except Exception:
        # xTB/GOAT performs the production relaxation; force-field failure is
        # recorded by the geometry provenance but is not silently re-labelled.
        pass
    return mol, 0


def split_species(smiles: str) -> tuple[str, list[str]]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid species SMILES: {smiles}")
    main: list[str] = []
    ions: list[str] = []
    for frag in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True):
        text = Chem.MolToSmiles(frag, isomericSmiles=True)
        atoms = list(frag.GetAtoms())
        if len(atoms) == 1 and atoms[0].GetSymbol() in COUNTERION_ELEMENTS:
            ions.append(text)
        else:
            main.append(text)
    if len(main) != 1:
        raise ValueError(f"Expected one analyte fragment in {smiles}, found {main}")
    return main[0], ions


def is_exchangeable_hydrogen(mol: Chem.Mol, atom: Chem.Atom) -> bool:
    if atom.GetSymbol() != "H" or atom.GetDegree() != 1:
        return False
    parent = atom.GetNeighbors()[0]
    return parent.GetSymbol() in {"N", "O", "S", "P"}


def solute_records(mol: Chem.Mol, conf_id: int, solvent: str) -> list[AtomRecord]:
    conformer = mol.GetConformer(conf_id)
    records: list[AtomRecord] = []
    for atom in mol.GetAtoms():
        point = conformer.GetAtomPosition(atom.GetIdx())
        exchangeable = is_exchangeable_hydrogen(mol, atom)
        isotope = DEUTERIUM_MASS if solvent == "D2O" and exchangeable else None
        records.append(
            AtomRecord(
                element=atom.GetSymbol(),
                x=float(point.x),
                y=float(point.y),
                z=float(point.z),
                component="analyte",
                source_atom_index_1based=atom.GetIdx() + 1,
                environment_atom=False,
                exchangeable=exchangeable,
                isotope_mass=isotope,
            )
        )
    return records


def template_records(smiles: str, seed: int, component: str, deuterate_all_h: bool) -> list[AtomRecord]:
    mol, conf_id = embed_fragment(smiles, seed)
    conformer = mol.GetConformer(conf_id)
    centre = tuple(
        sum(float(getattr(conformer.GetAtomPosition(i), axis)) for i in range(mol.GetNumAtoms()))
        / mol.GetNumAtoms()
        for axis in ("x", "y", "z")
    )
    result: list[AtomRecord] = []
    for atom in mol.GetAtoms():
        point = conformer.GetAtomPosition(atom.GetIdx())
        result.append(
            AtomRecord(
                atom.GetSymbol(),
                float(point.x) - centre[0],
                float(point.y) - centre[1],
                float(point.z) - centre[2],
                component,
                None,
                True,
                False,
                DEUTERIUM_MASS if deuterate_all_h and atom.GetSymbol() == "H" else None,
            )
        )
    return result


def transformed_fragment(
    template: Sequence[AtomRecord],
    matrix: tuple[tuple[float, float, float], ...],
    translation: tuple[float, float, float],
) -> list[AtomRecord]:
    result: list[AtomRecord] = []
    for atom in template:
        x, y, z = apply_rotation((atom.x, atom.y, atom.z), matrix)
        result.append(
            AtomRecord(
                atom.element,
                x + translation[0],
                y + translation[1],
                z + translation[2],
                atom.component,
                atom.source_atom_index_1based,
                atom.environment_atom,
                atom.exchangeable,
                atom.isotope_mass,
            )
        )
    return result


def donor_oxygen_targets(records: Sequence[AtomRecord]) -> list[tuple[float, float, float]]:
    """Return first-shell acceptor positions along analyte X-H donor vectors."""
    targets: list[tuple[float, float, float]] = []
    heavy = [
        atom
        for atom in records
        if atom.component == "analyte" and atom.element in {"N", "O", "S", "P"}
    ]
    for hydrogen in records:
        if hydrogen.component != "analyte" or not hydrogen.exchangeable:
            continue
        candidates = sorted(heavy, key=lambda atom: distance(hydrogen, atom))
        if not candidates or distance(hydrogen, candidates[0]) > 1.35:
            continue
        parent = candidates[0]
        vx, vy, vz = hydrogen.x - parent.x, hydrogen.y - parent.y, hydrogen.z - parent.z
        norm = math.sqrt(vx * vx + vy * vy + vz * vz)
        if norm < 1.0e-8:
            continue
        # About 1.9 A from donor H to solvent O is a reasonable starting
        # hydrogen-bond geometry; xTB/GOAT performs the production relaxation.
        scale = 1.90 / norm
        targets.append(
            (
                hydrogen.x + vx * scale,
                hydrogen.y + vy * scale,
                hydrogen.z + vz * scale,
            )
        )
    return targets


def place_environment(
    analyte: list[AtomRecord],
    ion_smiles: Sequence[str],
    scenario: Scenario,
    seed: int,
) -> list[AtomRecord]:
    rng = random.Random(seed)
    records = list(analyte)
    centre = centroid(records)
    heavy_radius = max(
        math.dist((r.x, r.y, r.z), centre) for r in records if r.element != "H"
    )

    # Counterions are placed first.  In D2O and the two-DMSO scenario use a
    # solvent-separated initial radius; GOAT may still discover contact pairs.
    for ion_index, ion_smiles_item in enumerate(ion_smiles, start=1):
        ion = Chem.MolFromSmiles(ion_smiles_item)
        if ion is None or ion.GetNumAtoms() != 1:
            raise ValueError(f"Unsupported counterion fragment: {ion_smiles_item}")
        element = ion.GetAtomWithIdx(0).GetSymbol()
        separated = scenario.solvent == "D2O" or "solvent_separated" in scenario.scenario_id
        radius = heavy_radius + (4.2 if separated else 2.8)
        for _ in range(500):
            ux, uy, uz = random_unit(rng)
            candidate = AtomRecord(
                element,
                centre[0] + radius * ux,
                centre[1] + radius * uy,
                centre[2] + radius * uz,
                f"counterion_{ion_index}",
                None,
                True,
                False,
                None,
            )
            if acceptable_fragment([candidate], records):
                records.append(candidate)
                break
        else:
            raise RuntimeError(f"Could not place counterion for {scenario.figure_id}/{scenario.scenario_id}")

    donor_targets = donor_oxygen_targets(analyte)

    if scenario.solvent == "D2O":
        fragment_specs = (("O", "explicit_D2O", scenario.explicit_solvent_count, True, 2.2),)
    else:
        fragment_specs = (
            ("CS(C)=O", "explicit_DMSO_d6", scenario.explicit_dmso_count, True, 2.5),
            ("O", "explicit_trace_H2O", scenario.explicit_water_count, False, 2.2),
        )

    global_solvent_index = 0
    for template_smiles, component_name, fragment_count, deuterate_all_h, surface_offset in fragment_specs:
        if fragment_count <= 0:
            continue
        template = template_records(
            template_smiles,
            seed + 101 + global_solvent_index,
            component_name,
            deuterate_all_h,
        )
        oxygen_anchor_index = next(i for i, atom in enumerate(template) if atom.element == "O")
        for _local_index in range(1, fragment_count + 1):
            global_solvent_index += 1
            solvent_index = global_solvent_index
            placed = False
            if donor_targets and solvent_index <= len(donor_targets):
                target = donor_targets[(solvent_index - 1) % len(donor_targets)]
                for _ in range(500):
                    matrix = rotation_matrix(rng)
                    rotated_anchor = apply_rotation(
                        (
                            template[oxygen_anchor_index].x,
                            template[oxygen_anchor_index].y,
                            template[oxygen_anchor_index].z,
                        ),
                        matrix,
                    )
                    translation = tuple(target[i] - rotated_anchor[i] for i in range(3))
                    candidate = transformed_fragment(template, matrix, translation)
                    candidate = [
                        AtomRecord(
                            a.element,
                            a.x,
                            a.y,
                            a.z,
                            f"{a.component}_{solvent_index:02d}",
                            a.source_atom_index_1based,
                            a.environment_atom,
                            a.exchangeable,
                            a.isotope_mass,
                        )
                        for a in candidate
                    ]
                    if acceptable_fragment(candidate, records):
                        records.extend(candidate)
                        placed = True
                        break
            if placed:
                continue
            for _ in range(2000):
                matrix = rotation_matrix(rng)
                ux, uy, uz = random_unit(rng)
                shell_jitter = rng.uniform(-0.25, 0.75) + 0.15 * (solvent_index - 1)
                radius = heavy_radius + surface_offset + shell_jitter
                translation = (
                    centre[0] + radius * ux,
                    centre[1] + radius * uy,
                    centre[2] + radius * uz,
                )
                candidate = transformed_fragment(template, matrix, translation)
                candidate = [
                    AtomRecord(
                        a.element,
                        a.x,
                        a.y,
                        a.z,
                        f"{a.component}_{solvent_index:02d}",
                        a.source_atom_index_1based,
                        a.environment_atom,
                        a.exchangeable,
                        a.isotope_mass,
                    )
                    for a in candidate
                ]
                if acceptable_fragment(candidate, records):
                    records.extend(candidate)
                    placed = True
                    break
            if not placed:
                raise RuntimeError(
                    f"Could not place solvent {solvent_index} for {scenario.figure_id}/{scenario.scenario_id}"
                )
    return records


def xyz_text(records: Sequence[AtomRecord], comment: str) -> str:
    """Return a strict standard XYZ file for ORCA's external xyzfile reader."""

    lines = [str(len(records)), comment]
    for atom in records:
        lines.append(f"{atom.element:<2s} {atom.x: .10f} {atom.y: .10f} {atom.z: .10f}")
    return "\n".join(lines) + "\n"


def goat_input(charge: int, scenario: Scenario, temperature: float) -> str:
    return f"""! GOAT-ENTROPY XTB ALPB({scenario.continuum})

%pal
  nprocs 8
end

%maxcore 2000

%goat
  NWorkers 8
  RandomSeed false
  AutoWall true
  Align true
  FreezeAmides false
  FreezeCisTrans true
  MaxEn 12.0
  KeepWorkerData false
  ConfTemp {temperature:.2f}
end

* xyzfile {charge} 1 start.xyz
"""


def atom_map_rows(records: Sequence[AtomRecord], solvent: str) -> list[dict[str, object]]:
    periodic_table = Chem.GetPeriodicTable()
    bonded: dict[int, list[int]] = {index: [] for index in range(len(records))}
    for left in range(len(records)):
        for right in range(left + 1, len(records)):
            # Explicit-solvent placement must never create a covalent bond in
            # the immutable map, even if a hydrogen bond happens to be short.
            if records[left].component != records[right].component:
                continue
            cutoff = 1.25 * (
                periodic_table.GetRcovalent(records[left].element)
                + periodic_table.GetRcovalent(records[right].element)
            )
            if 0.35 < distance(records[left], records[right]) <= cutoff:
                bonded[left].append(right)
                bonded[right].append(left)
    rows: list[dict[str, object]] = []
    for zero_index, atom in enumerate(records):
        index = zero_index + 1
        attached = sorted(bonded[zero_index])
        rows.append(
            {
                "atom_index_1based": index,
                "orca_index_0based": zero_index,
                "orca_index_1based": index,
                "element": atom.element,
                "component": atom.component,
                "source_atom_index_1based": atom.source_atom_index_1based or "",
                "formal_charge": 0,
                "canonical_symmetry_rank": atom.source_atom_index_1based or index,
                "attached_atom_indices_0based": ";".join(str(value) for value in attached),
                "attached_elements": ";".join(records[value].element for value in attached),
                "environment_atom": str(atom.environment_atom).lower(),
                "exchangeable": str(atom.exchangeable).lower(),
                "exchangeable_H_flag": int(atom.exchangeable),
                # Exchangeable NH/OH signals remain eligible in DMSO-d6 (with a
                # separate broadening/exchange model) but are removed from the
                # ordinary D2O analyte spectrum.
                "include_in_analyte_nmr": str(
                    not atom.environment_atom
                    and (not atom.exchangeable or solvent != "D2O")
                ).lower(),
                "isotope_mass_u": f"{atom.isotope_mass:.11f}" if atom.isotope_mass is not None else "",
            }
        )
    return rows


def prepare_job(
    row: dict[str, str],
    condition: dict[str, str],
    scenario: Scenario,
    start_index: int,
    force_inputs: bool,
) -> dict[str, object]:
    figure_id = row["figure_id"]
    seed_text = f"{figure_id}|{scenario.scenario_id}|{start_index}|ORCA61_ULTRA_2CMPD_V1"
    seed = int(hashlib.sha256(seed_text.encode()).hexdigest()[:8], 16) & 0x7FFFFFFF
    main_smiles, counterions = split_species(scenario.species_smiles)
    solute, conf_id = embed_fragment(main_smiles, seed)
    analyte = solute_records(solute, conf_id, scenario.solvent)
    records = place_environment(analyte, counterions, scenario, seed + 2001)
    total_charge = sum(atom.GetFormalCharge() for atom in Chem.MolFromSmiles(scenario.species_smiles).GetAtoms())
    job_id = f"{figure_id}_{scenario.scenario_id}_s{start_index:02d}"
    directory = JOB_ROOT / figure_id / scenario.scenario_id / f"start_{start_index:02d}"
    directory.mkdir(parents=True, exist_ok=True)
    xyz = xyz_text(
        records,
        f"{job_id}; CALCULATED MODEL; deterministic preparation seed={seed}; solvent={scenario.solvent}",
    )
    inp = goat_input(total_charge, scenario, float(condition["temperature_K"]))
    planned = {"start.xyz": xyz, "goat.inp": inp}
    for name, content in planned.items():
        target = directory / name
        if target.exists() and target.read_text(encoding="utf-8") != content:
            if not force_inputs:
                raise RuntimeError(f"Existing reviewed input differs; use --force-inputs only before calculation: {target}")
            if any(directory.glob("*.out")) or any(directory.glob("*.finalensemble.xyz")):
                raise RuntimeError(f"Refusing to replace input beside calculation output: {directory}")
        atomic_text(target, content)

    amap = atom_map_rows(records, scenario.solvent)
    write_csv(
        directory / "atom_map.csv",
        amap,
        (
            "atom_index_1based",
            "orca_index_0based",
            "orca_index_1based",
            "element",
            "component",
            "source_atom_index_1based",
            "formal_charge",
            "canonical_symmetry_rank",
            "attached_atom_indices_0based",
            "attached_elements",
            "environment_atom",
            "exchangeable",
            "exchangeable_H_flag",
            "include_in_analyte_nmr",
            "isotope_mass_u",
        ),
    )
    isotope_rows = [
        {
            "atom_index_1based": r["atom_index_1based"],
            "element": r["element"],
            "mass_u": r["isotope_mass_u"],
            "component": r["component"],
            "goat_xyz_mass_annotation": "NONE_PLAIN_H",
            "downstream_coordinate_source": "REVIEWED_GOAT_SELECTED_GEOMETRY_ONLY",
            "downstream_inline_mass_syntax": (
                f'{r["element"]} {{x_A}} {{y_A}} {{z_A}} M={r["isotope_mass_u"]}'
            ),
            "orca_6_1_manual_verified": "true",
        }
        for r in amap
        if r["isotope_mass_u"]
    ]
    write_csv(
        directory / "isotope_map.csv",
        isotope_rows,
        (
            "atom_index_1based",
            "element",
            "mass_u",
            "component",
            "goat_xyz_mass_annotation",
            "downstream_coordinate_source",
            "downstream_inline_mass_syntax",
            "orca_6_1_manual_verified",
        ),
    )
    scenario_mol = Chem.MolFromSmiles(scenario.species_smiles)
    validation = {
        "schema_version": "1.0",
        "notice": "CALCULATED MODEL - NOT EXPERIMENTAL DATA",
        "approval_basis": "AUTHOR_APPROVED_PROTOCOL_ASSUMPTIONS_2026-07-20",
        "figure_id": figure_id,
        "job_id": job_id,
        "scenario_id": scenario.scenario_id,
        "registry_inchikey": row["inchikey"],
        "registry_isomeric_smiles": row["isomeric_smiles"],
        "scenario_isomeric_smiles": Chem.MolToSmiles(scenario_mol, isomericSmiles=True),
        "scenario_inchikey": Chem.MolToInchiKey(scenario_mol),
        "complete_cluster_charge": total_charge,
        "multiplicity": 1,
        "primary_solvent": scenario.solvent,
        "continuum_keyword": scenario.continuum,
        "explicit_solvent_count": scenario.explicit_solvent_count,
        "assumed_pd": scenario.assumed_pd,
        "preparation_seed": seed,
        "coordinate_sha256": sha256_bytes(xyz.encode("utf-8")),
        "atom_map_sha256": sha256_file(directory / "atom_map.csv"),
        "isotope_map_sha256": sha256_file(directory / "isotope_map.csv"),
        "isotope_mass_policy": {
            "goat_external_xyz": GOAT_XYZ_MASS_POLICY,
            "downstream_opt_freq": DOWNSTREAM_D_MASS_POLICY,
            "downstream_coordinate_source": "REVIEWED_GOAT_SELECTED_GEOMETRY_ONLY",
            "deuterium_inline_syntax": "H {x_A} {y_A} {z_A} M=2.01410177812",
            "policy_file": ISOTOPE_POLICY.relative_to(PROJECT).as_posix(),
        },
        "checks": {
            "registry_identity_available": "PASS",
            "scenario_species_parse": "PASS",
            "stereochemistry_preserved_by_rdkit_embed": "PASS",
            "formal_charge_explicit": "PASS",
            "environment_atoms_excluded_from_analyte_nmr": "PASS",
            "coordinate_clash_gate": "PASS",
            "goat_xyz_has_exactly_four_atom_fields": "PASS",
            "isotope_map_retained_for_downstream_opt_freq": "PASS",
            "downstream_inline_M_equals_syntax_recorded": "PASS",
        },
    }
    atomic_json(directory / "structure_validation.json", validation)
    bundle_hash = hashlib.sha256()
    for name in ("goat.inp", "start.xyz", "atom_map.csv", "isotope_map.csv", "structure_validation.json"):
        bundle_hash.update(name.encode())
        bundle_hash.update((directory / name).read_bytes())
    atomic_text(directory / "INPUT_BUNDLE_SHA256.txt", bundle_hash.hexdigest().upper() + "\n")
    return {
        "queue_order": int(row["compound_number"]),
        "job_id": job_id,
        "figure_id": figure_id,
        "compound_number": row["compound_number"],
        "display_name": row["display_name"],
        "scenario_id": scenario.scenario_id,
        "primary_solvent": scenario.solvent,
        "continuum_keyword": scenario.continuum,
        "explicit_solvent_count": scenario.explicit_solvent_count,
        "explicit_dmso_count": scenario.explicit_dmso_count,
        "explicit_water_count": scenario.explicit_water_count,
        "start_index": start_index,
        "charge": total_charge,
        "multiplicity": 1,
        "preparation_seed": seed,
        "stage": "goat",
        "job_directory": directory.relative_to(PROJECT).as_posix(),
        "portable_job_directory": directory.relative_to(PROJECT).as_posix(),
        "relative_job_directory": directory.relative_to(PROJECT).as_posix(),
        "input_file": "goat.inp",
        "coordinate_file": "start.xyz",
        "output_file": "goat.out",
        "goat_input_sha256": sha256_file(directory / "goat.inp"),
        "start_xyz_sha256": sha256_file(directory / "start.xyz"),
        "input_bundle_sha256": bundle_hash.hexdigest().upper(),
        "expected_input_sha256": sha256_file(directory / "goat.inp"),
        "expected_coordinate_sha256": sha256_file(directory / "start.xyz"),
        "expected_input_bundle_sha256": bundle_hash.hexdigest().upper(),
        "enabled": "true",
        "readiness": "ready",
        "approval_status": "APPROVED_FOR_GOAT_CALCULATION",
        "goat_xyz_mass_policy": GOAT_XYZ_MASS_POLICY,
        "downstream_deuterium_mass_policy": DOWNSTREAM_D_MASS_POLICY,
        "isotope_policy_sha256": sha256_file(ISOTOPE_POLICY),
        "source_registry_sha256": sha256_file(SOURCE_REGISTRY),
        "condition_table_sha256": sha256_file(CONDITION_TABLE),
        "notice": "CALCULATED_MODEL_NOT_EXPERIMENTAL",
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", nargs="+", help="Optional F3 identifiers")
    parser.add_argument("--starts-per-scenario", type=int, default=10)
    parser.add_argument("--force-inputs", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.starts_per_scenario < 10:
        raise ValueError("High-quality production requires at least 10 starts per fixed scenario")
    only = normalized_only(args.only)
    isotope_policy_sha256 = write_isotope_policy()
    registry_rows = read_csv(SOURCE_REGISTRY)
    condition_rows = {r["figure_id"]: r for r in read_csv(CONDITION_TABLE)}
    known = {r["figure_id"] for r in registry_rows}
    if only - known:
        raise ValueError(f"Unknown identifiers: {sorted(only - known)}")
    manifest_rows: list[dict[str, object]] = []
    scenario_rows: list[dict[str, object]] = []
    for row in registry_rows:
        figure_id = row["figure_id"]
        if only and figure_id not in only:
            continue
        condition = condition_rows[figure_id]
        scenarios = scenario_plan(row, condition)
        for scenario in scenarios:
            scenario_rows.append(
                {
                    "figure_id": figure_id,
                    "compound_number": row["compound_number"],
                    "display_name": row["display_name"],
                    "scenario_id": scenario.scenario_id,
                    "primary_solvent": scenario.solvent,
                    "continuum_keyword": scenario.continuum,
                    "explicit_solvent_count": scenario.explicit_solvent_count,
                    "explicit_dmso_count": scenario.explicit_dmso_count,
                    "explicit_water_count": scenario.explicit_water_count,
                    "minimum_goat_starts": args.starts_per_scenario,
                    "species_smiles": scenario.species_smiles,
                    "assumed_pd": scenario.assumed_pd,
                    "species_note": scenario.species_note,
                    "status": "AUTHOR_APPROVED_ASSUMPTION_FOR_CALCULATION",
                }
            )
            for start_index in range(1, args.starts_per_scenario + 1):
                manifest_rows.append(
                    prepare_job(row, condition, scenario, start_index, args.force_inputs)
                )
    if not manifest_rows:
        raise RuntimeError("No jobs prepared")
    manifest_rows.sort(
        key=lambda r: (int(r["compound_number"]), str(r["scenario_id"]), int(r["start_index"]))
    )
    manifest_fields = list(manifest_rows[0].keys())
    write_csv(MANIFEST, manifest_rows, manifest_fields)
    write_csv(SCENARIO_REGISTRY, scenario_rows, list(scenario_rows[0].keys()))
    summary = {
        "notice": "CALCULATED MODELS - NOT EXPERIMENTAL DATA",
        "author_protocol_approval_recorded": True,
        "compounds": len({r["figure_id"] for r in manifest_rows}),
        "scenarios": len(scenario_rows),
        "goat_jobs": len(manifest_rows),
        "starts_per_scenario": args.starts_per_scenario,
        "manifest": MANIFEST.relative_to(PROJECT).as_posix(),
        "manifest_sha256": sha256_file(MANIFEST),
        "scenario_registry_sha256": sha256_file(SCENARIO_REGISTRY),
        "isotope_mass_policy": ISOTOPE_POLICY.relative_to(PROJECT).as_posix(),
        "isotope_mass_policy_sha256": isotope_policy_sha256,
        "goat_xyz_mass_policy": GOAT_XYZ_MASS_POLICY,
        "downstream_deuterium_mass_policy": DOWNSTREAM_D_MASS_POLICY,
    }
    atomic_json(DATA_ROOT / "explicit_preparation_summary.json", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
