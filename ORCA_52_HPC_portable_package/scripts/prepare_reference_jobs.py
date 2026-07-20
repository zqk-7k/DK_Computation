"""Prepare hash-locked, DRAFT method-matched NMR reference jobs.

No ORCA process is launched.  Raw calculated shieldings, molecular/cluster
models and experimental processing conventions are kept as separate records.
In particular, an isolated or microsolvated H3PO4 model is never labelled as
equivalent to an 85 wt% phosphoric-acid external standard.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from rdkit import Chem
from rdkit.Chem import AllChem

from workflow_schema import (
    DEUTERIUM_MASS_U,
    NOTICE,
    WorkflowValidationError,
    atomic_write_csv,
    atomic_write_json,
    atomic_write_text,
    sha256_bytes,
    sha256_file,
)


PROJECT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Component:
    name: str
    smiles: str
    count: int = 1
    environment: bool = False


@dataclass(frozen=True)
class Model:
    model_id: str
    label: str
    solvent_domain: str
    continuum: str
    nuclei: tuple[str, ...]
    components: tuple[Component, ...]
    target_rule: str
    raw_reference_role: str
    processing_convention: str
    limitation: str
    charge: int = 0
    multiplicity: int = 1


MODELS = (
    Model(
        "REF_DMSO_TMS",
        "Tetramethylsilane raw shielding model in DMSO continuum",
        "DMSO-d6",
        "DMSO",
        ("1H", "13C"),
        (Component("TMS", "C[Si](C)(C)C"),),
        "TMS_METHYL_H_AND_C",
        "Raw sigma(1H/13C) for TMS calculated at the target shielding method",
        "Experimental DMSO-d6 spectra may be processed against residual DMSO; the measured residual-solvent value and spectrometer convention remain separate metadata",
        "TMS raw shielding is a calculated zero-scale intermediate, not evidence that TMS was added to the author's sample",
    ),
    Model(
        "REF_DMSO_CFCL3",
        "CFCl3 raw 19F shielding model in DMSO continuum",
        "DMSO-d6",
        "DMSO",
        ("19F",),
        (Component("CFCl3", "FC(Cl)(Cl)Cl"),),
        "ALL_F",
        "Raw sigma(19F) for CFCl3 calculated at the target shielding method",
        "The experimental 19F convention (internal, external or indirect CFCl3 scale) must be copied from the acquisition/processing record",
        "A CPCM(DMSO) CFCl3 molecule is a declared calculation model and does not prove physical internal CFCl3 was present",
    ),
    Model(
        "REF_DMSO_31P_H3PO4_2DMSO_MODEL",
        "Microsolvated H3PO4 + 2 DMSO raw 31P shielding model",
        "DMSO-d6",
        "DMSO",
        ("31P",),
        (
            Component("H3PO4_model", "OP(=O)(O)O"),
            Component("DMSO-d6_environment", "CS(C)=O", count=2, environment=True),
        ),
        "ALL_P",
        "Raw sigma(31P) for a declared H3PO4·(DMSO)2 cluster model",
        "If the experiment used 85 wt% H3PO4 externally, record that external convention and its assigned delta=0 separately",
        "MODEL_ONLY_NOT_EQUIVALENT_TO_85WT_PERCENT_H3PO4; concentration, water activity, proton exchange and intermolecular structure of the external standard are not reproduced",
    ),
    Model(
        "REF_D2O_31P_H3PO4_4H2O_MODEL",
        "Microsolvated H3PO4 + 4 D2O raw 31P shielding model",
        "D2O",
        "WATER",
        ("31P",),
        (
            Component("H3PO4_model", "OP(=O)(O)O"),
            Component("D2O_environment", "O", count=4, environment=True),
        ),
        "ALL_P",
        "Raw sigma(31P) for a declared H3PO4·(H2O)4 electronic cluster with D isotope masses in the frequency input",
        "The actual D2O 31P reference may be external 85 wt% H3PO4 or an indirect scale; record the laboratory convention independently",
        "MODEL_ONLY_NOT_EQUIVALENT_TO_85WT_PERCENT_H3PO4; this dilute microsolvated model cannot reproduce concentrated acid speciation or bulk susceptibility",
    ),
    Model(
        "REF_D2O_DSS_NA_6D2O",
        "DSS sodium-salt + 6 D2O candidate raw 1H/13C model",
        "D2O",
        "WATER",
        ("1H", "13C"),
        (
            Component("DSS_anion", "[O-]S(=O)(=O)CCC[Si](C)(C)C"),
            Component("Na_counterion", "[Na+]", environment=True),
            Component("D2O_environment", "O", count=6, environment=True),
        ),
        "SILYL_METHYL_H_AND_C",
        "Raw sigma(1H/13C) for the trimethylsilyl resonance of a declared DSS-/Na+/6D2O cluster",
        "Use only if the actual acquisition used DSS with a compatible internal/external convention; concentration, pD and temperature must be recorded",
        "DSS ion pairing and hydration are model-dependent; this candidate cannot be substituted for TSP or an indirect reference without a new calibration decision",
    ),
    Model(
        "REF_D2O_TSP_NA_6D2O",
        "TSP sodium-salt + 6 D2O candidate raw 1H/13C model",
        "D2O",
        "WATER",
        ("1H", "13C"),
        (
            Component("TSP_anion", "[O-]C(=O)CC[Si](C)(C)C"),
            Component("Na_counterion", "[Na+]", environment=True),
            Component("D2O_environment", "O", count=6, environment=True),
        ),
        "SILYL_METHYL_H_AND_C",
        "Raw sigma(1H/13C) for the trimethylsilyl resonance of a declared TSP-/Na+/6D2O cluster",
        "Use only if the actual acquisition used the matching TSP/TSP-d4 convention; indirect referencing remains a separate non-molecular processing route",
        "Exact isotopologue, pD, sodium ion pairing and hydration must be confirmed; non-deuterated and TSP-d4 models are not silently interchangeable",
    ),
)


def method_fingerprint(model: Model) -> str:
    return (
        f"PBE0|pcSseg-2|AutoAux|RIJCOSX|GIAO|CPCM({model.continuum})|"
        "TightSCF|DefGrid3|ORCA_6.1.1"
    )


def _embed(smiles: str, seed: int) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise WorkflowValidationError(f"Cannot parse reference SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    if mol.GetNumAtoms() == 1:
        conformer = Chem.Conformer(1)
        conformer.SetAtomPosition(0, (0.0, 0.0, 0.0))
        mol.AddConformer(conformer, assignId=True)
        return mol
    params = AllChem.ETKDGv3()
    params.randomSeed = seed & 0x7FFFFFFF
    params.useRandomCoords = True
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise WorkflowValidationError(f"RDKit embedding failed for {smiles}")
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
        elif AllChem.UFFHasAllMoleculeParams(mol):
            AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
    except Exception as exc:
        raise WorkflowValidationError(f"Reference fragment force-field cleanup failed: {smiles}") from exc
    return mol


def _directions(count: int) -> list[tuple[float, float, float]]:
    if count <= 0:
        return []
    golden = math.pi * (3.0 - math.sqrt(5.0))
    result = []
    for index in range(count):
        z = 1.0 - 2.0 * (index + 0.5) / count
        radius = math.sqrt(max(0.0, 1.0 - z * z))
        angle = index * golden
        result.append((radius * math.cos(angle), radius * math.sin(angle), z))
    return result


def _target_selector(rule: str, mol: Chem.Mol) -> Callable[[int], bool]:
    if rule == "ALL_F":
        return lambda index: mol.GetAtomWithIdx(index).GetSymbol() == "F"
    if rule == "ALL_P":
        return lambda index: mol.GetAtomWithIdx(index).GetSymbol() == "P"
    if rule in {"TMS_METHYL_H_AND_C", "SILYL_METHYL_H_AND_C"}:
        methyl_carbons = {
            atom.GetIdx()
            for atom in mol.GetAtoms()
            if atom.GetSymbol() == "C"
            and any(neighbor.GetSymbol() == "Si" for neighbor in atom.GetNeighbors())
            and sum(neighbor.GetSymbol() == "H" for neighbor in atom.GetNeighbors()) == 3
        }
        target = set(methyl_carbons)
        for carbon in methyl_carbons:
            target.update(
                neighbor.GetIdx()
                for neighbor in mol.GetAtomWithIdx(carbon).GetNeighbors()
                if neighbor.GetSymbol() == "H"
            )
        return lambda index: index in target
    raise WorkflowValidationError(f"Unknown reference-site rule: {rule}")


def _minimum_allowed_intercomponent_distance(element_a: str, element_b: str) -> float:
    """Conservative clash threshold that still permits ordinary H bonds.

    The threshold is intentionally a little below the sum of van der Waals
    radii.  It is a collision screen, not a claim that the initial cluster is a
    converged solvation structure.
    """

    table = Chem.GetPeriodicTable()
    radius_a = float(table.GetRvdw(table.GetAtomicNumber(element_a)))
    radius_b = float(table.GetRvdw(table.GetAtomicNumber(element_b)))
    return max(1.60, 0.62 * (radius_a + radius_b))


def intercomponent_contact_audit(records: list[dict[str, object]]) -> dict[str, object]:
    closest_distance = math.inf
    closest_pair: tuple[int, int] | None = None
    minimum_clearance_ratio = math.inf
    limiting_pair: tuple[int, int] | None = None
    for left_position, left in enumerate(records):
        for right_position in range(left_position + 1, len(records)):
            right = records[right_position]
            if left["component"] == right["component"]:
                continue
            distance = math.sqrt(
                sum(
                    (float(left[key]) - float(right[key])) ** 2
                    for key in ("x_A", "y_A", "z_A")
                )
            )
            required = _minimum_allowed_intercomponent_distance(
                str(left["element"]), str(right["element"])
            )
            ratio = distance / required
            if distance < closest_distance:
                closest_distance = distance
                closest_pair = (left_position + 1, right_position + 1)
            if ratio < minimum_clearance_ratio:
                minimum_clearance_ratio = ratio
                limiting_pair = (left_position + 1, right_position + 1)
    if closest_pair is None or limiting_pair is None:
        return {
            "cross_component_pair_count": 0,
            "minimum_cross_component_distance_A": None,
            "minimum_nonbonded_clearance_ratio": None,
            "closest_cross_component_pair_1based": None,
            "limiting_clearance_pair_1based": None,
            "automatic_intercomponent_collision_check": "NOT_APPLICABLE_SINGLE_COMPONENT",
        }
    return {
        "minimum_cross_component_distance_A": round(closest_distance, 6),
        "minimum_nonbonded_clearance_ratio": round(minimum_clearance_ratio, 6),
        "closest_cross_component_pair_1based": list(closest_pair),
        "limiting_clearance_pair_1based": list(limiting_pair),
        "automatic_intercomponent_collision_check": (
            "PASS" if minimum_clearance_ratio >= 1.0 - 1e-10 else "FAIL"
        ),
    }


def build_model_records(model: Model) -> list[dict[str, object]]:
    instances: list[tuple[Component, int, Chem.Mol]] = []
    seed_base = int(hashlib.sha256(model.model_id.encode()).hexdigest()[:8], 16)
    for component in model.components:
        for copy_index in range(1, component.count + 1):
            instances.append((component, copy_index, _embed(component.smiles, seed_base + len(instances) * 101)))
    directions = _directions(max(0, len(instances) - 1))
    records: list[dict[str, object]] = []
    for instance_index, (component, copy_index, mol) in enumerate(instances):
        conformer = mol.GetConformer()
        points = [conformer.GetAtomPosition(index) for index in range(mol.GetNumAtoms())]
        center = (
            sum(point.x for point in points) / len(points),
            sum(point.y for point in points) / len(points),
            sum(point.z for point in points) / len(points),
        )
        centered = [
            (point.x - center[0], point.y - center[1], point.z - center[2]) for point in points
        ]
        if instance_index == 0:
            translation = (0.0, 0.0, 0.0)
        else:
            direction = directions[instance_index - 1]
            translation = None
            for shell_step in range(0, 45):
                shell = 3.0 + 0.20 * shell_step
                candidate_translation = tuple(shell * value for value in direction)
                collision = False
                for local_index, atom in enumerate(mol.GetAtoms()):
                    candidate = tuple(
                        centered[local_index][axis] + candidate_translation[axis]
                        for axis in range(3)
                    )
                    for prior in records:
                        distance = math.sqrt(
                            sum(
                                (candidate[axis] - float(prior[key])) ** 2
                                for axis, key in enumerate(("x_A", "y_A", "z_A"))
                            )
                        )
                        required = _minimum_allowed_intercomponent_distance(
                            atom.GetSymbol(), str(prior["element"])
                        ) * 1.08
                        if distance < required:
                            collision = True
                            break
                    if collision:
                        break
                if not collision:
                    translation = candidate_translation
                    break
            if translation is None:
                raise WorkflowValidationError(
                    f"Could not place a collision-free component in reference model {model.model_id}"
                )
        selector = _target_selector(model.target_rule, mol) if instance_index == 0 else lambda _index: False
        for local_index, atom in enumerate(mol.GetAtoms()):
            x = centered[local_index][0] + translation[0]
            y = centered[local_index][1] + translation[1]
            z = centered[local_index][2] + translation[2]
            element = atom.GetSymbol()
            target = selector(local_index)
            mass = ""
            if element == "H" and model.solvent_domain == "D2O":
                if component.name == "D2O_environment":
                    mass = f"{DEUTERIUM_MASS_U:.11f}"
                elif instance_index == 0 and any(
                    neighbor.GetSymbol() in {"O", "N", "S", "P"} for neighbor in atom.GetNeighbors()
                ):
                    mass = f"{DEUTERIUM_MASS_U:.11f}"
            if element == "H" and component.name == "DMSO-d6_environment":
                mass = f"{DEUTERIUM_MASS_U:.11f}"
            records.append(
                {
                    "atom_index_1based": len(records) + 1,
                    "element": element,
                    "x_A": x,
                    "y_A": y,
                    "z_A": z,
                    "component": f"{component.name}_{copy_index:02d}",
                    "component_smiles": component.smiles,
                    "environment_atom": str(component.environment).lower(),
                    "reference_site": str(target).lower(),
                    "requested_nucleus": (
                        {"H": "1H", "C": "13C", "F": "19F", "P": "31P"}.get(element, "")
                        if target
                        else ""
                    ),
                    "isotope_mass_u": mass,
                }
            )
    return records


def xyz_text(records: list[dict[str, object]], model: Model) -> str:
    lines = [str(len(records)), f"{NOTICE}; {model.model_id}; DRAFT reference model"]
    for row in records:
        lines.append(
            f"{row['element']:<2s} {float(row['x_A']): .12f} {float(row['y_A']): .12f} {float(row['z_A']): .12f}"
        )
    return "\n".join(lines) + "\n"


def opt_input(records: list[dict[str, object]], model: Model) -> str:
    lines = [
        f"# {NOTICE}",
        f"# {model.model_id}: declared raw-reference model; not an experimental spectrum.",
        f"! r2SCAN-3c CPCM({model.continuum}) Opt TightOpt Freq TightSCF",
        "%pal nprocs 6 end",
        "%maxcore 4000",
        "%geom MaxIter 500 end",
        "%freq Temp 298.15 end",
        f"* xyz {model.charge} {model.multiplicity}",
    ]
    for row in records:
        mass = f" M={row['isotope_mass_u']}" if row["isotope_mass_u"] else ""
        lines.append(
            f"  {row['element']:<2s} {float(row['x_A']): .12f} {float(row['y_A']): .12f} {float(row['z_A']): .12f}{mass}"
        )
    lines.extend(("*", ""))
    return "\n".join(lines)


def shielding_template(records: list[dict[str, object]], model: Model) -> str:
    by_element: dict[str, list[int]] = {}
    for row in records:
        if row["reference_site"] == "true":
            by_element.setdefault(str(row["element"]), []).append(int(row["atom_index_1based"]))
    lines = [
        f"# {NOTICE}",
        f"# Raw shielding only for the declared reference resonance of {model.model_id}.",
        f"! PBE0 pcSseg-2 AutoAux RIJCOSX TightSCF DefGrid3 CPCM({model.continuum})",
        "%pal nprocs 6 end",
        "%maxcore 4000",
        f"* xyzfile {model.charge} {model.multiplicity} optimized.xyz",
        "%eprnmr",
        "  Ori = GIAO",
    ]
    isotope = {"H": 1, "C": 13, "F": 19, "P": 31}
    for element in ("H", "C", "F", "P"):
        if element in by_element:
            lines.append(
                f"  Nuclei = {','.join(str(value) for value in by_element[element])} "
                f"{{ shift, ist = {isotope[element]} }};"
            )
    lines.extend(("end", ""))
    return "\n".join(lines)


def prepare_reference_jobs(project: Path, replace_draft: bool = False) -> list[dict[str, object]]:
    project = Path(project).resolve()
    root = project / "reference_jobs"
    manifest_rows: list[dict[str, object]] = []
    registry_entries: list[dict[str, object]] = []
    gate_decisions: list[dict[str, object]] = []
    for model in MODELS:
        directory = root / model.model_id
        if directory.exists():
            calculated = list(directory.glob("*.out")) + list(directory.glob("*.gbw"))
            if calculated:
                raise WorkflowValidationError(f"Refusing to replace a started reference job: {calculated[0]}")
            if not replace_draft:
                raise WorkflowValidationError(f"Reference DRAFT directory exists: {directory}")
        records = build_model_records(model)
        contact_audit = intercomponent_contact_audit(records)
        if contact_audit["automatic_intercomponent_collision_check"] == "FAIL":
            raise WorkflowValidationError(
                f"Reference cluster contains an intercomponent collision: {model.model_id}"
            )
        target_nuclei = {
            str(row["requested_nucleus"]) for row in records if row["requested_nucleus"]
        }
        if target_nuclei != set(model.nuclei):
            raise WorkflowValidationError(
                f"Reference target map mismatch for {model.model_id}: {target_nuclei} != {model.nuclei}"
            )
        xyz = xyz_text(records, model)
        opt = opt_input(records, model)
        shielding = shielding_template(records, model)
        map_fields = (
            "atom_index_1based",
            "element",
            "x_A",
            "y_A",
            "z_A",
            "component",
            "component_smiles",
            "environment_atom",
            "reference_site",
            "requested_nucleus",
            "isotope_mass_u",
        )
        atomic_write_text(directory / "start.xyz", xyz)
        atomic_write_text(directory / "opt_freq.inp", opt)
        atomic_write_text(directory / "shielding.inp.template", shielding)
        atomic_write_csv(directory / "reference_atom_map.csv", records, map_fields)
        declaration = {
            "schema_version": "2.0",
            "notice": NOTICE,
            "approval_status": "DRAFT_NOT_APPROVED",
            "model_id": model.model_id,
            "label": model.label,
            "solvent_domain": model.solvent_domain,
            "continuum_keyword": model.continuum,
            "requested_nuclei": model.nuclei,
            "raw_calculation_reference": model.raw_reference_role,
            "experimental_processing_convention": model.processing_convention,
            "model_limitation": model.limitation,
            "raw_shielding_and_processing_convention_separated": True,
            "raw_shielding_aggregation_policy": (
                "PRESERVE_EACH_REQUESTED_ATOM_VALUE; within each nucleus and this single approved "
                "reference model, report n, arithmetic mean, population standard deviation, minimum "
                "and maximum over the declared chemically equivalent reference-site atoms; never "
                "average across different models, isotopologues, solvent domains, conformers or "
                "experimental processing conventions"
            ),
            "equivalent_to_85wt_percent_H3PO4_external_standard": False,
            "method_fingerprint": method_fingerprint(model),
            "method_fingerprint_sha256": sha256_bytes(method_fingerprint(model).encode()),
            "structure_sha256": sha256_bytes(xyz.encode()),
            "opt_input_sha256": sha256_bytes(opt.encode()),
            "shielding_template_sha256": sha256_bytes(shielding.encode()),
            "initial_cluster_contact_audit": contact_audit,
            "initial_cluster_geometry_role": (
                "DRAFT_COLLISION_SCREENED_START_ONLY; requires reviewed optimization/frequency "
                "minimum before shielding and is not a converged solvation ensemble"
            ),
        }
        atomic_write_json(directory / "model_declaration.json", declaration)
        manifest_rows.append(
            {
                "queue_order": len(manifest_rows) + 1,
                "model_id": model.model_id,
                "label": model.label,
                "solvent_domain": model.solvent_domain,
                "continuum_keyword": model.continuum,
                "requested_nuclei": ";".join(model.nuclei),
                "charge": model.charge,
                "multiplicity": model.multiplicity,
                "job_directory": str(directory.resolve()),
                "start_xyz_sha256": sha256_file(directory / "start.xyz"),
                "reference_atom_map_sha256": sha256_file(directory / "reference_atom_map.csv"),
                "model_declaration_sha256": sha256_file(directory / "model_declaration.json"),
                "opt_input_file": "opt_freq.inp",
                "opt_input_sha256": sha256_file(directory / "opt_freq.inp"),
                "shielding_template_file": "shielding.inp.template",
                "shielding_template_sha256": sha256_file(directory / "shielding.inp.template"),
                "shielding_method_fingerprint": method_fingerprint(model),
                "shielding_method_fingerprint_sha256": sha256_bytes(
                    method_fingerprint(model).encode()
                ),
                "initial_collision_check": contact_audit[
                    "automatic_intercomponent_collision_check"
                ],
                "minimum_cross_component_distance_A": contact_audit[
                    "minimum_cross_component_distance_A"
                ],
                "minimum_nonbonded_clearance_ratio": contact_audit[
                    "minimum_nonbonded_clearance_ratio"
                ],
                "approval_status": "DRAFT_NOT_APPROVED",
                "notice": NOTICE,
            }
        )
        gate_decisions.append(
            {
                "model_id": model.model_id,
                "approval_status": "DRAFT_NOT_APPROVED",
                "automated_intercomponent_collision_screen": contact_audit[
                    "automatic_intercomponent_collision_check"
                ],
                "structure_and_speciation_review": "REVIEW_REQUIRED",
                "reference_site_index_review": "REVIEW_REQUIRED",
                "method_match_review": "REVIEW_REQUIRED",
                "experimental_convention_review": "REVIEW_REQUIRED",
                "model_declaration_file": str((directory / "model_declaration.json").resolve()),
                "model_declaration_sha256": sha256_file(directory / "model_declaration.json"),
                "decision_reason": "",
            }
        )
        for nucleus in model.nuclei:
            registry_entries.append(
                {
                    "reference_id": f"{model.model_id}__{nucleus}",
                    "approval_status": "DRAFT_NOT_APPROVED",
                    "model_id": model.model_id,
                    "nucleus": nucleus,
                    "solvent_domain": model.solvent_domain,
                    "raw_calculation_reference": model.raw_reference_role,
                    "experimental_reference_convention": model.processing_convention,
                    "model_limitation": model.limitation,
                    "method_fingerprint": method_fingerprint(model),
                    "input_file": "",
                    "output_file": "",
                    "input_sha256": "",
                    "output_sha256": "",
                    "raw_shielding_value_ppm": None,
                    "processing_offset_or_regression": "TO_LOCK_FROM_ACTUAL_ACQUISITION_AND_REFERENCE_DATA",
                }
            )
    manifest_path = project / "data" / "reference_job_manifest_DRAFT.csv"
    atomic_write_csv(manifest_path, manifest_rows, list(manifest_rows[0]))
    manifest_hash = sha256_file(manifest_path)
    gate = {
        "schema_version": "2.0",
        "notice": NOTICE,
        "approval_status": "DRAFT_NOT_APPROVED",
        "stage": "reference_model_protocol",
        "reference_manifest_file": str(manifest_path.resolve()),
        "reference_manifest_sha256": manifest_hash,
        "raw_shielding_reference_check": "REVIEW_REQUIRED",
        "experimental_processing_convention_check": "REVIEW_REQUIRED",
        "no_false_external_standard_equivalence_check": "REVIEW_REQUIRED",
        "d2o_1h_13c_reference_route_decision": {
            "status": "UNRESOLVED_REVIEW_REQUIRED",
            "allowed_routes": [
                "REF_D2O_DSS_NA_6D2O",
                "REF_D2O_TSP_NA_6D2O",
                "INDIRECT_REFERENCE_NO_ORCA_JOB",
            ],
            "selected_route": "",
            "actual_standard_or_indirect_metadata": "",
            "decision_reason": "",
        },
        "model_decisions": gate_decisions,
        "reviewer": "",
        "reviewed_at": "",
        "decision_notes": [
            "Approve a raw reference model only after structure, reference-site indices and method matching are reviewed.",
            "DSS and TSP are alternative D2O candidates; retain only the model matching the actual acquisition, or use an explicitly documented indirect-reference route.",
            "When the protocol is promoted, mark every model either APPROVED_FOR_REFERENCE_CALCULATION or REJECTED_REFERENCE_MODEL; a rejected model requires a reason.",
            "A microsolvated H3PO4 calculation is not equivalent to an 85 wt% H3PO4 external standard. The latter remains an experimental processing convention.",
            "After reference calculations, bind the actual shielding input/output hashes and raw shielding values in the reference registry; do not enter fitted experimental offsets as raw quantum-chemical shieldings.",
        ],
    }
    atomic_write_json(project / "jobs" / "reference_model_protocol_DRAFT.json", gate)
    registry = {
        "schema_version": "2.0",
        "status": "DRAFT_NOT_APPROVED",
        "notice": NOTICE,
        "policy": (
            "Raw calculated shieldings, physical reference models and experimental processing "
            "conventions are distinct. A reference is reusable only after nucleus, method, solvent "
            "domain, model, input/output hashes and processing convention are all locked."
        ),
        "reference_model_protocol_gate": str(
            (project / "jobs" / "reference_model_protocol_DRAFT.json").resolve()
        ),
        "reference_model_protocol_gate_sha256": sha256_file(
            project / "jobs" / "reference_model_protocol_DRAFT.json"
        ),
        "indirect_reference_routes": [
            {
                "solvent_domain": "D2O",
                "nuclei": ["1H", "13C"],
                "approval_status": "DRAFT_NOT_APPROVED",
                "route": "instrument/vendor or IUPAC-compatible indirect referencing",
                "required_metadata": [
                    "spectrometer reference frequency",
                    "reference nucleus and ratio table/version",
                    "processing software and exact command/parameter",
                    "temperature and lock solvent",
                ],
                "limitation": "This is a processing convention and has no molecular ORCA reference job.",
            }
        ],
        "reference_jobs": registry_entries,
    }
    atomic_write_json(project / "config" / "reference_job_registry.json", registry)
    return manifest_rows


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--project", type=Path, default=PROJECT)
    value.add_argument("--replace-draft", action="store_true")
    return value


def main() -> None:
    args = parser().parse_args()
    rows = prepare_reference_jobs(args.project, replace_draft=args.replace_draft)
    print(f"Prepared {len(rows)} DRAFT reference models; ORCA was not started.")


if __name__ == "__main__":
    main()
