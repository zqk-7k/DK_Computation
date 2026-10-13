"""Independent, read-only static audit of every prepared explicit-solvent job.

Does not import the generator or launcher, execute QM, or write any package file.
JSON evidence is printed to stdout for a reviewer to save after inspection.
This validates prepared starting models, not sampling or prediction accuracy.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[1]
SCENARIOS = {
    "DMSO_00_W00_continuum": (0, 0, 3),
    "DMSO_02_W00_primary": (2, 0, 4),
    "DMSO_02_W01_trace_water_sensitivity": (2, 1, 7),
}
EXACT_ISOTOPE_MASS = {("O", 18): 17.99915961286, ("H", 2): 2.01410177812}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def readj(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def mol(path):
    value = Chem.MolFromMolBlock(Path(path).read_text(encoding='utf-8'), removeHs=False)
    if value is None:
        raise ValueError("Unreadable MOL: " + str(path))
    return value


def canon(value):
    return Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(value)), isomericSmiles=True)


def atom_signature(value):
    return [(a.GetAtomicNum(), a.GetIsotope(), a.GetFormalCharge()) for a in value.GetAtoms()]


def bond_signature(value):
    return sorted((min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                   max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                   str(b.GetBondType()), b.GetIsAromatic()) for b in value.GetBonds())


def analyse(root):
    errors, warnings = [], []
    checks = Counter()
    job_results = []

    def check(condition, name, job="GLOBAL", detail=None):
        checks[name] += 1
        if not condition:
            errors.append({"scope": job, "check": name, "detail": detail})

    provenance = readj(root / "registry/source_provenance.json")
    source = Path(provenance["source_v1"])
    jobs = readj(root / "registry/jobs.json")
    species = readj(root / "registry/species.json")
    products = readj(root / "registry/products.json")
    species_by_id = {s["species_id"]: s for s in species}
    source_species = readj(source / "registry/species.json")
    source_active = [s for s in source_species if s["solvent"] == "DMSO" and s.get("default_run_enabled")]
    check(species == source_active, "species_registry_exact_source_subset")
    check(sha(source / "registry/species.json") == provenance["source_species_sha256"], "source_species_hash")
    check(sha(source / "registry/products.json") == provenance["source_products_sha256"], "source_products_hash")
    check(products == readj(source / "registry/products.json"), "product_registry_preserved")
    check(len(products) == 64 and len(species) == 67, "expected_products_and_models")
    plan_path=root/'registry/sampling_plan_v4.json'
    planned=readj(plan_path) if plan_path.is_file() else None
    expected_counts={(r['species_id'],r['scenario_id']):r['n_starts'] for r in planned['species']} if planned else {(s['species_id'],scenario):starts for s in species for scenario,(_,_,starts) in SCENARIOS.items()}
    check(len(jobs) == sum(expected_counts.values()), "expected_initial_job_count", detail=len(jobs))
    check(len({j["job_id"] for j in jobs}) == len(jobs), "unique_job_ids")
    check(len({j["job_directory"] for j in jobs}) == len(jobs), "unique_job_directories")
    templates = {s["species_id"]: mol(source / "structures" / s["species_id"] / "model.mol") for s in species}
    for sid in species_by_id:
        check(sha(source / "structures" / sid / "model.mol") == sha(root / "solute_models" / sid / "model.mol"),
              "copied_source_model_hash", sid)
    for item in provenance["reference_files"]:
        check(sha(root / "reference_51" / item["path"]) == item["sha256"], "reference_snapshot_hash", item["path"])

    groups = defaultdict(list)
    initial_methods = Counter()
    coordinate_hashes = defaultdict(list)
    analyte_coordinates = defaultdict(dict)
    isotope_totals = Counter()
    cluster_mass_range = [float("inf"), 0.0]
    min_interfragment_distance = float("inf")
    minimum_clearance_margin = float("inf")
    max_xyz_mol_delta = 0.0
    max_mass_difference = 0.0
    max_mass_adjusted_difference = 0.0
    pt = Chem.GetPeriodicTable()
    hashes = {
        "goat.inp": "input_sha256", "start.xyz": "xyz_sha256",
        "atom_map.json": "atom_map_sha256", "cluster.mol": "cluster_mol_sha256",
        "provenance.json": "provenance_sha256",
    }
    for job in jobs:
        jid, sid = job["job_id"], job["species_id"]
        starting_error_count = len(errors)
        try:
            path = (root / job["job_directory"]).resolve()
            check(path.is_relative_to(root.resolve()) and path != root.resolve(), "job_path_confined", jid)
            if not path.is_relative_to(root.resolve()) or path == root.resolve():
                continue
            for name, field in hashes.items():
                check(sha(path / name) == job[field], "manifest_file_hash", jid, name)
            amap = readj(path / "atom_map.json")
            imap = readj(path / "isotope_map.json")
            prov = readj(path / "provenance.json")
            cluster = mol(path / "cluster.mol")
            template = templates[sid]
            solute_n, total_n = template.GetNumAtoms(), cluster.GetNumAtoms()
            nd, nw, starts = SCENARIOS[job["scenario_id"]]
            groups[(sid, job["scenario_id"])].append(job["start_index"])
            check(job["DMSO_d6_count"] == nd and job["H2O_count"] == nw, "manifest_composition", jid)
            check(total_n == solute_n + 10 * nd + 3 * nw == job["atom_count"], "cluster_atom_count", jid)
            check(job["solute_atom_count"] == solute_n and len(amap) == total_n, "solute_and_map_atom_count", jid)
            check(job["status"] == "PREPARED_NOT_RUN" and job["orca_runtime_validated"] is False,
                  "no_runtime_success_claim", jid)
            check(prov["cluster_optimization"] == "NOT_RUN" and prov["orca_calculation"] == "NOT_RUN",
                  "no_qm_success_claim", jid)
            check(prov["analyte_start"]["global_sampling_converged"] is False, "no_sampling_convergence_claim", jid)
            check(prov["source_species"] == species_by_id[sid], "provenance_species_consistency", jid)
            check(prov["source_model_sha256"] == sha(source / "structures" / sid / "model.mol"),
                  "provenance_source_hash", jid)
            check(prov["composition"] == {"DMSO_d6": nd, "H2O": nw}, "provenance_composition", jid)
            initial_methods[prov["analyte_start"]["method"]] += 1

            lines = (path / "start.xyz").read_text(encoding="utf-8").splitlines()
            parsed = [line.split() for line in lines[2:]]
            check(int(lines[0]) == total_n and len(parsed) == total_n, "xyz_atom_count", jid)
            check(all(len(row) == 4 for row in parsed), "xyz_plain_coordinate_columns", jid)
            xyz = np.array([[float(x) for x in row[1:]] for row in parsed])
            check(bool(np.all(np.isfinite(xyz))), "finite_xyz", jid)
            check([r[0] for r in parsed] == [a.GetSymbol() for a in cluster.GetAtoms()], "xyz_mol_atom_order", jid)
            delta = float(np.max(np.abs(xyz - cluster.GetConformer().GetPositions())))
            max_xyz_mol_delta = max(max_xyz_mol_delta, delta)
            check(delta <= 0.000051, "xyz_mol_coordinate_rounding", jid, delta)
            coordinate_digest = hashlib.sha256(np.round(xyz, 8).astype("<f8").tobytes()).hexdigest()
            coordinate_hashes[coordinate_digest].append(jid)
            analyte_coordinates[sid][job["start_index"]] = xyz[:solute_n]

            check(atom_signature(cluster)[:solute_n] == atom_signature(template), "analyte_atom_order_isotopes_charges", jid)
            fragment_indices = Chem.GetMolFrags(cluster)
            fragment_molecules = Chem.GetMolFrags(cluster, asMols=True, sanitizeFrags=True)
            check(len(fragment_indices) == 1 + nd + nw, "fragment_count", jid)
            check(tuple(fragment_indices[0]) == tuple(range(solute_n)), "analyte_is_first_intact_fragment", jid)
            analyte = fragment_molecules[0]
            check(bond_signature(analyte) == bond_signature(template), "analyte_bond_graph_order", jid)
            check(canon(analyte) == canon(template), "stored_analyte_stereochemistry", jid)
            q = Chem.Mol(analyte)
            conf = q.GetConformer()
            for i in range(solute_n):
                conf.SetAtomPosition(i, tuple(float(x) for x in xyz[i]))
            for bond in q.GetBonds():
                bond.SetBondDir(Chem.BondDir.NONE)
            Chem.AssignStereochemistryFrom3D(q, confId=0, replaceExistingTags=True)
            check(canon(q) == canon(template), "xyz_derived_analyte_stereochemistry", jid)

            expected_18 = [a.GetIdx() for a in template.GetAtoms() if a.GetSymbol() == "O" and a.GetIsotope() == 18]
            actual_18 = [a.GetIdx() for a in cluster.GetAtoms() if a.GetSymbol() == "O" and a.GetIsotope() == 18]
            expected_label_count = 0 if species_by_id[sid]["compound_id"] == "TMS" else species_by_id[sid]["labelled_18O_count"]
            check(actual_18 == expected_18 and len(expected_18) == expected_label_count,
                  "18O_count_and_original_indices", jid)
            actual_D = [a.GetIdx() for a in cluster.GetAtoms() if a.GetSymbol() == "H" and a.GetIsotope() == 2]
            check(len(actual_D) == 6 * nd and all(i >= solute_n for i in actual_D), "12D_in_two_DMSO_only", jid)
            isotope_totals["18O_atoms_across_jobs"] += len(actual_18)
            isotope_totals["D_atoms_across_jobs"] += len(actual_D)
            component_ids = [r["component"] for r in amap]
            check(component_ids[:solute_n] == ["analyte"] * solute_n, "analyte_component_map", jid)
            for idx, (indices, fragment) in enumerate(zip(fragment_indices[1:], fragment_molecules[1:]), start=1):
                expected_fragment_counts = Counter({("C", 0): 2, ("H", 2): 6, ("O", 0): 1, ("S", 0): 1}) if idx <= nd else Counter({("H", 0): 2, ("O", 0): 1})
                observed = Counter((a.GetSymbol(), a.GetIsotope()) for a in fragment.GetAtoms())
                check(observed == expected_fragment_counts, "environment_fragment_isotopic_formula", jid, idx)
                check(len({component_ids[i] for i in indices}) == 1, "fragment_component_label_consistency", jid, idx)
                plain = Chem.Mol(fragment)
                for atom in plain.GetAtoms():
                    atom.SetIsotope(0)
                check(canon(plain) == ("CS(C)=O" if idx <= nd else "O"), "environment_covalent_topology", jid, idx)

            check(imap == [r for r in amap if r["isotope"]], "isotope_sidecar_exact_subset", jid)
            custom_mass = 0.0
            isotope_table_adjustment = 0.0
            for i, (atom, row) in enumerate(zip(cluster.GetAtoms(), amap)):
                symbol, isotope = atom.GetSymbol(), atom.GetIsotope()
                env = i >= solute_n
                check(row["index_zero_based"] == i and row["index_one_based"] == i + 1 and row["element"] == symbol
                      and row["isotope"] == isotope and row["formal_charge"] == atom.GetFormalCharge(), "map_atom_identity", jid, i)
                check(sorted(row["neighbours_zero_based"]) == sorted(a.GetIdx() for a in atom.GetNeighbors()), "map_bond_neighbours", jid, i)
                check(row["environment_atom"] == env and row["source_atom_index_zero_based"] == (None if env else i),
                      "map_source_and_environment_indices", jid, i)
                check(row["include_in_analyte_nmr"] == (not env and symbol in ("H", "C")), "solvent_excluded_from_analyte_nmr", jid, i)
                expected_mass = EXACT_ISOTOPE_MASS.get((symbol, isotope))
                check(row["isotope_mass_u"] == expected_mass, "exact_isotope_mass_constants", jid, i)
                check(not isotope or expected_mass is not None, "all_isotopes_supported", jid, i)
                custom_mass += expected_mass if expected_mass is not None else pt.GetMostCommonIsotopeMass(atom.GetAtomicNum())
                if expected_mass is not None:
                    isotope_table_adjustment += expected_mass - pt.GetMassForIsotope(atom.GetAtomicNum(), isotope)
            rdkit_mass = Descriptors.ExactMolWt(cluster)
            mass_delta = abs(custom_mass - rdkit_mass)
            # RDKit 2026.03.3 uses 17.999161 for 18O rather than the explicitly
            # declared 17.99915961286. Account exactly for its table convention,
            # rather than mistaking the known table offset for a composition error.
            max_mass_difference = max(max_mass_difference, mass_delta)
            adjusted_delta = abs(custom_mass - (rdkit_mass + isotope_table_adjustment))
            max_mass_adjusted_difference = max(max_mass_adjusted_difference, adjusted_delta)
            check(adjusted_delta < 1e-8, "cluster_exact_mass_independent_crosscheck", jid, adjusted_delta)
            cluster_mass_range = [min(cluster_mass_range[0], custom_mass), max(cluster_mass_range[1], custom_mass)]
            check(Chem.GetFormalCharge(cluster) == job["charge"] == species_by_id[sid]["charge"], "cluster_total_charge", jid)

            radii = np.array([pt.GetRvdw(a.GetAtomicNum()) for a in cluster.GetAtoms()])
            atomic_numbers = np.array([a.GetAtomicNum() for a in cluster.GetAtoms()])
            covalent_radii = np.array([pt.GetRcovalent(a.GetAtomicNum()) for a in cluster.GetAtoms()])
            for a, b in np.array([(a, b) for a in range(total_n) for b in range(a)], dtype=int):
                distance = float(np.linalg.norm(xyz[a] - xyz[b]))
                covalent_sum = covalent_radii[a] + covalent_radii[b]
                bonded = cluster.GetBondBetweenAtoms(int(a), int(b)) is not None
                if bonded:
                    check(0.57 * covalent_sum < distance < 1.65 * covalent_sum, "bond_length_sanity", jid, [int(a), int(b), distance])
                elif component_ids[a] == component_ids[b]:
                    check(distance >= 0.55 * covalent_sum, "intrafragment_nonbonded_sanity", jid, [int(a), int(b), distance])
                else:
                    threshold = max(0.64 * (radii[a] + radii[b]), 1.15 if 1 in (atomic_numbers[a], atomic_numbers[b]) else 1.65)
                    min_interfragment_distance = min(min_interfragment_distance, distance)
                    minimum_clearance_margin = min(minimum_clearance_margin, distance - threshold)
                    check(distance >= threshold - 1e-8, "interfragment_clearance", jid, [int(a), int(b), distance, threshold])

            inp = (path / "goat.inp").read_text(encoding="utf-8")
            check(re.search(r"(?m)^\* xyzfile " + str(job["charge"]) + " " + str(job["multiplicity"]) + r" start\.xyz\s*$", inp) is not None,
                  "input_charge_multiplicity_local_xyz", jid)
            check("GOAT-ENTROPY XTB ALPB(DMSO)" in inp, "explicit_goat_alpb_input", jid)
            check(not any(p.suffix.lower() in (".out", ".gbw", ".hess") for p in path.iterdir()), "no_qm_output_in_start_directory", jid)
            job_results.append({"job_id": jid, "passed": len(errors) == starting_error_count,
                                "xyz_sha256": sha(path / "start.xyz"), "exact_cluster_mass_u": round(custom_mass, 9)})
        except Exception as exc:
            errors.append({"scope": jid, "check": "audit_exception", "detail": type(exc).__name__ + ": " + str(exc)})

    for (sid,scenario),starts in expected_counts.items():
        check(sorted(groups[(sid, scenario)]) == list(range(1, starts + 1)), "complete_initial_start_indices", sid, scenario)
    check(set(groups)==set(expected_counts),'planned_scenarios_only')
    duplicate_groups = [values for values in coordinate_hashes.values() if len(values) > 1]
    check(not duplicate_groups, "no_exact_duplicate_cluster_coordinates", detail=duplicate_groups)
    aligned_pair_distances = []
    near_duplicates = []
    for sid, entries in analyte_coordinates.items():
        heavy = [a.GetIdx() for a in templates[sid].GetAtoms() if a.GetAtomicNum() > 1]
        if len(heavy) < 2:
            continue
        pair_matrix = []
        for index, coordinates in sorted(entries.items()):
            xyz = coordinates[heavy]
            dist = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
            pair_matrix.append((index, dist[np.triu_indices(len(heavy), 1)]))
        for i, (first_index, first) in enumerate(pair_matrix):
            for second_index, second in pair_matrix[:i]:
                rms = float(np.sqrt(np.mean((first-second)**2)))
                aligned_pair_distances.append(rms)
                if rms < 0.01:
                    near_duplicates.append([sid, first_index, second_index, round(rms, 8)])
    warnings.append({"code": "NO_QM_OR_ACCURACY_VALIDATION", "detail": "All tests concern input models, atom maps and simple geometry. They do not establish stable solvation structures, equilibrium populations, sampling convergence, or NMR accuracy."})
    warnings.append({"code": "ISOTOPES_SIDECAR_ONLY_FOR_GOAT", "detail": "Standard XYZ and GOAT input do not apply the preserved 18O/D masses to ORCA. This initial-input audit verifies exact masses in maps only. The v3 downstream renderer restores isotope masses in inline coordinates; implementation tests and separate small engineering calculations do not constitute validation of target thermochemistry or NMR accuracy."})
    warnings.append({"code": "START_UNIQUENESS_NOT_CONFORMER_INDEPENDENCE", "detail": "Exact coordinate uniqueness is not independent basin discovery. Near-identical heavy-atom pair-distance starts, including rigid molecules and fallback perturbations, are counted but are not converged conformers.", "near_identical_pairs_below_0.01_A_RMS": len(near_duplicates), "examples": near_duplicates[:10]})
    report = {
        "audit": "INDEPENDENT_ALL_PREPARED_INITIAL_CLUSTERS_STATIC_QA",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": "STATIC_INPUT_CHECKS_PASS_WITH_LIMITATIONS" if not errors else "FAIL",
        "package_root": str(root), "original_source_root": str(source),
        "reference_commit": provenance["reference_commit"],
        "auditor_script_sha256": sha(Path(__file__)), "rdkit_version": rdBase.rdkitVersion,
        "scope": {"products": len(products), "source_active_models_including_TMS": len(species),
                  "fixed_composition_scenarios": len(groups), "manifest_jobs": len(jobs),
                  "jobs_fully_audited": len(job_results), "jobs_passing": sum(j["passed"] for j in job_results),
                  "all_jobs_audited": len(job_results) == len(jobs), "QM_jobs_executed_by_auditor": 0},
        "inputs": {"jobs_manifest_sha256": sha(root / "registry/jobs.json"), "source_provenance_sha256": sha(root / "registry/source_provenance.json")},
        "test_counts": dict(sorted(checks.items())), "error_count": len(errors), "errors": errors,
        "initial_geometry_methods": dict(initial_methods), "isotope_inventory": dict(isotope_totals),
        "geometry": {"minimum_interfragment_distance_A": min_interfragment_distance,
                     "minimum_vdw_clearance_margin_A": minimum_clearance_margin,
                     "max_xyz_mol_component_rounding_difference_A": max_xyz_mol_delta,
                     "exact_duplicate_coordinate_groups": duplicate_groups,
                     "analyte_near_identical_pair_count": len(near_duplicates)},
        "mass": {"explicit_isotope_constants_u": {"18O": EXACT_ISOTOPE_MASS[("O", 18)], "2H": EXACT_ISOTOPE_MASS[("H", 2)]},
                 "primary_mass_sources_verified_2026_10_13": [
                     "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=O",
                     "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=H"],
                 "minimum_cluster_exact_mass_u": cluster_mass_range[0], "maximum_cluster_exact_mass_u": cluster_mass_range[1],
                 "maximum_RDKit_comparison_difference_u": max_mass_difference,
                 "RDKit_isotope_table_u": {"18O": pt.GetMassForIsotope(8, 18), "2H": pt.GetMassForIsotope(1, 2)},
                 "maximum_RDKit_isotope_table_adjusted_difference_u": max_mass_adjusted_difference,
                 "interpretation": "Neutral sum of isotope masses; not m/z or an experimental mass spectrum."},
        "all_job_result_digest_sha256": hashlib.sha256(json.dumps(job_results, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
        "warnings": warnings,
        "scientific_acceptance": "NOT_EVALUATED_NO_QM_NO_SAMPLING_CONVERGENCE_NO_PREDICTED_SPECTRA",
    }
    if len(job_results) != len(jobs):
        report["status"] = "FAIL"
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    report = analyse(args.root)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["status"] == "STATIC_INPUT_CHECKS_PASS_WITH_LIMITATIONS" else 1


if __name__ == "__main__":
    sys.exit(main())
