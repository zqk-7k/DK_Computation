#!/usr/bin/env python3
"""Build the frozen, solute-only structures used by the 24 h workflow.

The approved 3D coordinates are taken from the already-audited explicit-solvent
preparation package.  Every atom marked as an environment atom (solvent or
counterion) is removed.  The resulting structures therefore contain *no*
explicit DMSO or water molecules.  Salt identities remain in the registry and
HRMS policy, while the electronic NMR model uses the approved continuum species.

This script is a reproducibility helper.  The generated structures are committed
to the portable package, so it is not required on the HPC server unless the
source preparation package is also present.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT.parent / "ORCA_52_HPC_portable_package"

SCENARIO_SELECTION = {
    "F3_016": "DMSO_cation_continuum_01",
    "F3_017": "DMSO_cation_continuum_01",
    "F3_019": "D2O_06W_primary",
    "F3_020": "D2O_06W_primary",
    "F3_021": "D2O_08W_primary",
    "F3_035": "D2O_04W_primary",
    "F3_043": "D2O_04W_primary",
    "F3_044": "D2O_04W_primary",
}

MODEL_SMILES_OVERRIDES = {
    "F3_016": "NC(=O)[C@@H]([NH3+])CO",
    "F3_017": "CC(C)C[C@H]([NH3+])C(N)=O",
    "F3_019": "NC(=O)C[C@H]([NH3+])C(=O)[O-]",
    "F3_020": "NC(=O)CC[C@H]([NH3+])C(=O)[O-]",
    "F3_021": "NC(=[NH2+])NCCC[C@H]([NH3+])C(N)=O",
    "F3_023": "[NH3+][C@@H](Cc1ccc(O)c(NC(=O)c2ccccc2O)c1)C(=O)[O-]",
    "F3_035": "NS(=O)(=O)[O-]",
    "F3_043": "NP(N)(N)=O",
    "F3_044": "NP(N)(=O)[O-]",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def truth(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_xyz(path: Path) -> tuple[str, list[tuple[str, float, float, float]]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        raise ValueError(f"Malformed XYZ: {path}")
    count = int(lines[0].strip())
    atoms: list[tuple[str, float, float, float]] = []
    for line in lines[2 : 2 + count]:
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"Malformed XYZ atom line in {path}: {line!r}")
        atoms.append((fields[0], float(fields[1]), float(fields[2]), float(fields[3])))
    if len(atoms) != count:
        raise ValueError(f"XYZ atom count mismatch in {path}: {len(atoms)} != {count}")
    return lines[1], atoms


def write_xyz(path: Path, atoms: list[tuple[str, float, float, float]], comment: str) -> None:
    lines = [str(len(atoms)), comment]
    lines.extend(f"{el:<2s} {x: .10f} {y: .10f} {z: .10f}" for el, x, y, z in atoms)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def formula_from_elements(elements: list[str]) -> str:
    counts = Counter(elements)
    order: list[str] = []
    if "C" in counts:
        order.append("C")
        if "H" in counts:
            order.append("H")
    order.extend(sorted(set(counts) - set(order)))
    return "".join(el + (str(counts[el]) if counts[el] != 1 else "") for el in order)


def formal_charge(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit cannot parse approved model SMILES: {smiles}")
    return int(Chem.GetFormalCharge(mol))


def scenario_for(figure_id: str, rows: list[dict[str, str]]) -> dict[str, str]:
    wanted = SCENARIO_SELECTION.get(figure_id, "DMSO_01_primary")
    matches = [row for row in rows if row["figure_id"] == figure_id and row["scenario_id"] == wanted]
    if len(matches) != 1:
        raise ValueError(f"Expected one scenario for {figure_id}/{wanted}; found {len(matches)}")
    return matches[0]


def build_compounds(source: Path) -> list[dict[str, Any]]:
    registry_path = ROOT / "registry" / "Figure3_52_ORCA_reference_state_registry.csv"
    scenario_path = ROOT / "registry" / "source_explicit_scenario_registry.csv"
    registry = read_csv(registry_path)
    scenarios = read_csv(scenario_path)
    if [row["figure_id"] for row in registry] != [f"F3_{i:03d}" for i in range(1, 53)]:
        raise ValueError("The Figure 3 registry is not the expected ordered F3_001..F3_052 set")

    summary: list[dict[str, Any]] = []
    for reg in registry:
        figure_id = reg["figure_id"]
        scenario = scenario_for(figure_id, scenarios)
        scenario_id = scenario["scenario_id"]
        start_dir = source / "01_goat_explicit" / figure_id / scenario_id / "start_01"
        atom_map_source = start_dir / "atom_map.csv"
        xyz_source = start_dir / "start.xyz"
        if not atom_map_source.is_file() or not xyz_source.is_file():
            raise FileNotFoundError(f"Approved 3D source missing for {figure_id}: {start_dir}")

        source_map = read_csv(atom_map_source)
        _, source_atoms = read_xyz(xyz_source)
        if len(source_map) != len(source_atoms):
            raise ValueError(f"Atom-map/XYZ size mismatch for {figure_id}")

        kept_rows: list[dict[str, Any]] = []
        kept_atoms: list[tuple[str, float, float, float]] = []
        for source_index, (row, atom) in enumerate(zip(source_map, source_atoms), start=1):
            if row["element"] != atom[0]:
                raise ValueError(f"Element mismatch for {figure_id} source atom {source_index}")
            if truth(row.get("environment_atom", "false")):
                continue
            exchangeable = truth(row.get("exchangeable", "false"))
            include = truth(row.get("include_in_analyte_nmr", "true"))
            if scenario["primary_solvent"] == "D2O" and exchangeable:
                include = False
            new_index = len(kept_atoms) + 1
            kept_atoms.append(atom)
            kept_rows.append(
                {
                    "atom_index_1based": new_index,
                    "element": atom[0],
                    "component": "analyte",
                    "source_atom_index_1based": source_index,
                    "environment_atom": "false",
                    "exchangeable": str(exchangeable).lower(),
                    "include_in_analyte_nmr": str(include).lower(),
                    "reference_signal": "false",
                }
            )

        model_smiles = MODEL_SMILES_OVERRIDES.get(figure_id, scenario["species_smiles"])
        model_mol = Chem.AddHs(Chem.MolFromSmiles(model_smiles))
        model_element_order = [atom.GetSymbol() for atom in model_mol.GetAtoms()]
        expected_elements = Counter(model_element_order)
        observed_elements = Counter(atom[0] for atom in kept_atoms)
        if expected_elements != observed_elements:
            raise ValueError(
                f"Approved species/3D composition mismatch for {figure_id}: "
                f"SMILES={dict(expected_elements)}, XYZ={dict(observed_elements)}"
            )
        observed_element_order = [atom[0] for atom in kept_atoms]
        if model_element_order != observed_element_order:
            raise ValueError(
                f"Approved RDKit/3D atom order mismatch for {figure_id}; "
                "topology could not be transferred without an atom mapping"
            )

        out_dir = ROOT / "structures" / figure_id
        out_dir.mkdir(parents=True, exist_ok=True)
        xyz_out = out_dir / "start.xyz"
        map_out = out_dir / "atom_map.csv"
        bonds_out = out_dir / "bonds.csv"
        write_xyz(
            xyz_out,
            kept_atoms,
            f"{figure_id}; solute-only approved reference state; implicit {scenario['continuum_keyword']}; CALCULATED MODEL",
        )
        write_csv(
            map_out,
            kept_rows,
            [
                "atom_index_1based",
                "element",
                "component",
                "source_atom_index_1based",
                "environment_atom",
                "exchangeable",
                "include_in_analyte_nmr",
                "reference_signal",
            ],
        )
        bond_rows: list[dict[str, Any]] = []
        for bond in model_mol.GetBonds():
            bond_rows.append(
                {
                    "atom_i_1based": bond.GetBeginAtomIdx() + 1,
                    "atom_j_1based": bond.GetEndAtomIdx() + 1,
                    "bond_order": float(bond.GetBondTypeAsDouble()),
                    "is_aromatic": str(bond.GetIsAromatic()).lower(),
                }
            )
        write_csv(
            bonds_out,
            bond_rows,
            ["atom_i_1based", "atom_j_1based", "bond_order", "is_aromatic"],
        )
        charge = formal_charge(model_smiles)
        observed_h = sum(row["element"] == "H" and row["include_in_analyte_nmr"] == "true" for row in kept_rows)
        observed_c = sum(row["element"] == "C" for row in kept_rows)
        metadata = {
            "figure_id": figure_id,
            "compound_number": int(reg["compound_number"]),
            "display_name": reg["display_name"],
            "model_smiles": model_smiles,
            "inchikey_as_drawn": reg["inchikey"],
            "model_formula": formula_from_elements([atom[0] for atom in kept_atoms]),
            "charge": charge,
            "multiplicity": 1,
            "experimental_solvent_label": scenario["primary_solvent"],
            "orca_continuum_keyword": scenario["continuum_keyword"],
            "assumed_pd": scenario["assumed_pd"],
            "selected_source_scenario": scenario_id,
            "explicit_environment_atoms_in_model": 0,
            "observable_1H_count_after_exchange_policy": observed_h,
            "observable_13C_count": observed_c,
            "source_xyz_sha256": sha256(xyz_source),
            "source_atom_map_sha256": sha256(atom_map_source),
            "start_xyz_sha256": sha256(xyz_out),
            "atom_map_sha256": sha256(map_out),
            "bonds_sha256": sha256(bonds_out),
            "species_note": scenario["species_note"],
            "status": "FROZEN_SOLUTE_ONLY_IMPLICIT_MODEL",
        }
        (out_dir / "structure.json").write_text(
            json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        summary.append(metadata)

    fields = list(summary[0].keys())
    write_csv(ROOT / "config" / "primary_species.csv", summary, fields)
    (ROOT / "config" / "primary_species.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def embed_reference(name: str, smiles: str, charge: int, solvent: str) -> dict[str, Any]:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = 20260720 if name == "TMS_DMSO" else 20260721
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError(f"ETKDG failed for {name}")
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
    conformer = mol.GetConformer()
    atoms: list[tuple[str, float, float, float]] = []
    reference_indices: set[int] = set()
    silicon = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == "Si"]
    if len(silicon) != 1:
        raise ValueError(f"Expected exactly one Si in {name}")
    si_atom = silicon[0]
    methyl_carbons = [
        nbr for nbr in si_atom.GetNeighbors()
        if nbr.GetSymbol() == "C" and sum(n.GetSymbol() == "H" for n in nbr.GetNeighbors()) == 3
    ]
    for carbon in methyl_carbons:
        reference_indices.add(carbon.GetIdx())
        reference_indices.update(n.GetIdx() for n in carbon.GetNeighbors() if n.GetSymbol() == "H")
    expected_methyl_count = 4 if name == "TMS_DMSO" else 3
    if len(methyl_carbons) != expected_methyl_count:
        raise ValueError(
            f"Expected {expected_methyl_count} Si-methyl reference carbons in {name}; "
            f"found {len(methyl_carbons)}"
        )

    rows: list[dict[str, Any]] = []
    for index, atom in enumerate(mol.GetAtoms()):
        position = conformer.GetAtomPosition(index)
        atoms.append((atom.GetSymbol(), position.x, position.y, position.z))
        rows.append(
            {
                "atom_index_1based": index + 1,
                "element": atom.GetSymbol(),
                "component": "reference",
                "source_atom_index_1based": index + 1,
                "environment_atom": "false",
                "exchangeable": "false",
                "include_in_analyte_nmr": "true",
                "reference_signal": str(index in reference_indices).lower(),
            }
        )
    out_dir = ROOT / "references" / name
    out_dir.mkdir(parents=True, exist_ok=True)
    xyz = out_dir / "start.xyz"
    amap = out_dir / "atom_map.csv"
    write_xyz(xyz, atoms, f"{name}; computational zero-shift reference; implicit {solvent}")
    write_csv(amap, rows, list(rows[0].keys()))
    return {
        "reference_id": name,
        "display_name": "tetramethylsilane" if name == "TMS_DMSO" else "3-(trimethylsilyl)propionate anion",
        "smiles": smiles,
        "charge": charge,
        "multiplicity": 1,
        "orca_continuum_keyword": solvent,
        "reference_atom_indices_1based": ";".join(str(index + 1) for index in sorted(reference_indices)),
        "reference_nucleus_policy": (
            "mean of four Si-methyl carbons / twelve Si-methyl protons"
            if name == "TMS_DMSO"
            else "mean of three Si-methyl carbons / nine Si-methyl protons"
        ),
        "xyz_sha256": sha256(xyz),
        "atom_map_sha256": sha256(amap),
    }


def build_references() -> list[dict[str, Any]]:
    rows = [
        embed_reference("TMS_DMSO", "C[Si](C)(C)C", 0, "DMSO"),
        embed_reference("TSP_MINUS_WATER", "[O-]C(=O)CC[Si](C)(C)C", -1, "WATER"),
    ]
    write_csv(ROOT / "references" / "reference_registry.csv", rows, list(rows[0].keys()))
    return rows


def write_manifest() -> None:
    rows: list[dict[str, Any]] = []
    for folder in (ROOT / "structures", ROOT / "references", ROOT / "registry"):
        for path in sorted(folder.rglob("*")):
            if path.is_file():
                rows.append(
                    {
                        "path": path.relative_to(ROOT).as_posix(),
                        "size_bytes": path.stat().st_size,
                        "sha256": sha256(path),
                    }
                )
    write_csv(ROOT / "config" / "input_manifest_sha256.csv", rows, ["path", "size_bytes", "sha256"])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-package", type=Path, default=DEFAULT_SOURCE)
    args = parser.parse_args()
    if not args.source_package.is_dir():
        raise FileNotFoundError(f"Source preparation package not found: {args.source_package}")
    compounds = build_compounds(args.source_package)
    references = build_references()
    write_manifest()
    print(f"Built {len(compounds)} solute-only structures and {len(references)} references in {ROOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
