"""Independent scientific failure-mode regression tests; never execute QM.

Uses existing initial coordinates only as test fixtures, never as accepted
GOAT/DFT results. Mutations are in-memory or in temporary test directories.
"""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def readj(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def initial_fixture(species="F2_E01_rep_R__DMSO", scenario="DMSO_02_W01_trace_water_sensitivity"):
    folder = ROOT / "01_goat_explicit" / species / scenario / "start_01"
    # RDKit's native narrow Windows path API cannot reliably open the required
    # Chinese package name; Python reads Unicode paths before MOL parsing.
    value = Chem.MolFromMolBlock((folder / "cluster.mol").read_text(encoding="utf-8"), removeHs=False)
    if value is None:
        raise RuntimeError("Independent fixture MOL is unreadable")
    return value, readj(folder / "atom_map.json")


class IndependentGeometryAndInputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from clusterflow import geometry, inputs
        cls.geometry, cls.inputs = geometry, inputs
        cls.protocol = readj(ROOT / "config/protocol.json")

    def setUp(self):
        self.template, self.amap = initial_fixture()

    def issues(self, candidate, amap=None):
        return self.geometry.validate_geometry(self.template, candidate, self.amap if amap is None else amap, self.protocol)

    def test_initial_fixture_is_only_valid_geometry_not_completed_qm(self):
        self.assertFalse(self.issues(self.template))
        self.assertFalse(self.geometry.validate_atom_map(self.template, self.amap))

    def test_solvent_dissociation_is_rejected(self):
        changed = Chem.Mol(self.template)
        conf = changed.GetConformer()
        component = next(r["component"] for r in self.amap if r["environment_atom"])
        for row in self.amap:
            if row["component"] == component:
                idx = row["index_zero_based"]
                conf.SetAtomPosition(idx, tuple(np.array(conf.GetAtomPosition(idx)) + [30, 0, 0]))
        self.assertTrue(self.issues(changed), "Detached DMSO must not pass a fixed-cluster geometry gate")

    def test_water_proton_transfer_is_rejected(self):
        changed = Chem.Mol(self.template)
        conf = changed.GetConformer()
        oxygen = next(r["index_zero_based"] for r in self.amap if not r["environment_atom"] and r["element"] == "O")
        water_h = next(r["index_zero_based"] for r in self.amap if r["component"].startswith("H2O") and r["element"] == "H")
        water_o = next(r["index_zero_based"] for r in self.amap if r["component"].startswith("H2O") and r["element"] == "O")
        target = np.array(conf.GetAtomPosition(oxygen))
        direction = target - np.array(conf.GetAtomPosition(water_o))
        direction /= np.linalg.norm(direction)
        conf.SetAtomPosition(water_h, tuple(target + direction))
        self.assertTrue(self.issues(changed), "Water H transferred to analyte must not retain old topology acceptance")

    def test_reflected_enantiomer_is_rejected(self):
        changed = Chem.Mol(self.template)
        conf = changed.GetConformer()
        for i in range(changed.GetNumAtoms()):
            x, y, z = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (-x, y, z))
        self.assertTrue(self.issues(changed), "A global mirror reflection must not preserve the chosen epoxide stereoisomer")

    def test_18O_isotope_loss_is_rejected(self):
        changed = Chem.Mol(self.template)
        oxygen = next(a for a in changed.GetAtoms() if a.GetIsotope() == 18)
        oxygen.SetIsotope(0)
        self.assertTrue(self.issues(changed))

    def test_D_isotope_loss_is_rejected(self):
        changed = Chem.Mol(self.template)
        deuterium = next(a for a in changed.GetAtoms() if a.GetIsotope() == 2)
        deuterium.SetIsotope(0)
        self.assertTrue(self.issues(changed))

    def test_solvent_cannot_forge_analyte_mask(self):
        mask = copy.deepcopy(self.amap)
        row = next(r for r in mask if r["environment_atom"] and r["element"] == "C")
        row.update(environment_atom=False, include_in_analyte_nmr=True, source_atom_index_zero_based=0)
        self.assertTrue(self.geometry.validate_atom_map(self.template, mask))

    def test_atom_map_isotope_mismatch_is_rejected(self):
        mask = copy.deepcopy(self.amap)
        row = next(r for r in mask if r["isotope"] == 18)
        row["isotope"] = 0
        self.assertTrue(self.geometry.validate_atom_map(self.template, mask))

    def test_duplicate_index_in_atom_map_is_rejected(self):
        mask = copy.deepcopy(self.amap)
        mask[1]["index_zero_based"] = 0
        self.assertTrue(self.geometry.validate_atom_map(self.template, mask))

    def render(self, stage, config=None, value=None, amap=None):
        return self.inputs.render_input(value or self.template, stage, config or self.protocol,
                                        "DMSO", 0, 1, atom_map=amap or self.amap)

    def test_shielding_whitelist_is_one_based_and_analyte_only(self):
        text, _ = self.render("shielding")
        lines = re.findall(r"(?im)^\s*Nuclei\s*=\s*([^\{]+)\{[^\}]*shift[^\}]*\}", text)
        self.assertTrue(lines, "Explicit property-index whitelist is required")
        observed = set()
        for line in lines:
            self.assertNotIn("all", line.lower())
            observed.update(int(x) for x in re.findall(r"\d+", line))
        expected = {r["index_zero_based"] + 1 for r in self.amap if not r["environment_atom"] and r["element"] in ("H", "C")}
        self.assertEqual(expected, observed)
        self.assertGreater(text.lower().index("%eprnmr"), text.lower().index("* xyz"))
        self.assertNotIn("revdsd", text.lower())

    def test_frequency_input_carries_all_exact_isotope_masses(self):
        text, _ = self.render("optfreq")
        masses = [float(v) for v in re.findall(r"\bM\s*=\s*([0-9.]+)", text)]
        self.assertEqual(sum(abs(v - 17.99915961286) < 1e-9 for v in masses), 1)
        self.assertEqual(sum(abs(v - 2.01410177812) < 1e-9 for v in masses), 12)
        self.assertNotIn("revdsd", text.lower())

    def test_J_requested_pairs_exclude_all_environment_atoms(self):
        pairs = self.inputs.selected_j_pairs(self.template, self.amap, cutoff=1000.0)
        allowed = {r["index_zero_based"] for r in self.amap if not r["environment_atom"]}
        observed = set()
        for pair in pairs:
            if isinstance(pair, dict):
                i, j = pair["atom_i_zero_based"], pair["atom_j_zero_based"]
            else:
                i, j = pair[:2]
            self.assertTrue(i in allowed and j in allowed)
            self.assertNotEqual(i, j)
            observed.add(tuple(sorted((i, j))))
        hydrogens = [r["index_zero_based"] for r in self.amap if not r["environment_atom"] and r["element"] == "H"]
        expected = {(i, j) for i in hydrogens for j in hydrogens if i < j}
        self.assertEqual(expected, observed)

    def test_Br_J_requires_explicit_mixed_basis_configuration(self):
        value, amap = initial_fixture("F3_N15_source__DMSO")
        settings = copy.deepcopy(self.protocol)
        settings["j_PBE_pcJ3"].pop("bromine_basis", None)
        settings["j"].pop("bromine_basis", None)
        with self.assertRaises(ValueError):
            self.render("j_PBE_pcJ3", settings, value, amap)

    def test_fluorine_is_a_J_endpoint_even_without_shift_observation_flag(self):
        value, amap = initial_fixture("F4_D04_rep_R__DMSO")
        fluorines = {r["index_zero_based"] for r in amap if not r["environment_atom"] and r["element"] == "F"}
        self.assertTrue(fluorines)
        self.assertTrue(all(not amap[i]["include_in_analyte_nmr"] for i in fluorines))
        pairs = self.inputs.selected_j_pairs(value, amap, cutoff=1000.0)
        endpoints = {i for pair in pairs for i in pair}
        self.assertTrue(fluorines.issubset(endpoints), "19F coupling endpoints are distinct from the H/C shift mask")
        expected = {(min(i, j), max(i, j)) for i in fluorines for j, row in enumerate(amap)
                    if not row["environment_atom"] and row["element"] in ("H", "C")}
        self.assertTrue(expected.issubset({tuple(p) for p in pairs}))

    def test_J_atom_pairs_are_zero_based_and_S_basis_is_supported(self):
        text, pairs = self.render("j_PBE_pcJ3")
        block = next(line for line in text.splitlines() if "SpinSpinAtomPairs" in line)
        printed = [list(map(int, p)) for p in re.findall(r"\{\s*(\d+)\s+(\d+)\s*\}", block)]
        self.assertEqual(printed, pairs)
        self.assertTrue(any(a.GetSymbol() == "S" for a in self.template.GetAtoms()))
        self.assertNotIn("NewGTO S", text)

    def test_Br_mixed_basis_fingerprint_cannot_reuse_light_atom_J_benchmark(self):
        light = self.inputs.method_fingerprint(self.protocol, "j_PBE_pcJ3", "DMSO", has_br=False)
        bromine = self.inputs.method_fingerprint(self.protocol, "j_PBE_pcJ3", "DMSO", has_br=True)
        self.assertNotEqual(light, bromine)

    def test_prohibited_revDSD_is_never_rendered(self):
        settings = copy.deepcopy(self.protocol)
        settings["shielding"]["keywords"] = "revDSD-PBEP86-D4 pcSseg-3"
        with self.assertRaises(ValueError):
            self.render("shielding", settings)


class IndependentCollection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from clusterflow import collect
        cls.collect = collect

    def test_weighted_average_preserves_given_complete_weights(self):
        result = self.collect.weighted_average([{0: 10.0}, {0: 30.0}], [0.25, 0.75])
        self.assertEqual(result, {0: 25.0})

    def test_unormalized_negative_or_nonfinite_weights_rejected(self):
        for weights in ([0.2, 0.3], [-0.1, 1.1], [float("nan"), 1.0], [True, 0.0]):
            with self.subTest(weights=weights), self.assertRaises(ValueError):
                self.collect.weighted_average([{0: 1.0}, {0: 2.0}], weights)

    def test_missing_property_never_becomes_zero(self):
        with self.assertRaises(ValueError):
            self.collect.weighted_average([{0: 1.0, 1: 4.0}, {0: 2.0}], [0.5, 0.5])

    def test_nonfinite_property_rejected_even_in_zero_weight_conformer(self):
        with self.assertRaises(ValueError):
            self.collect.weighted_average([{0: 1.0}, {0: float("nan")}], [1.0, 0.0])

    def test_collector_rejects_environment_scope_forgery(self):
        _, amap = initial_fixture()
        row = next(r for r in amap if r["environment_atom"] and r["element"] == "C")
        row.update(environment_atom=False, include_in_analyte_nmr=True)
        with self.assertRaises(ValueError):
            self.collect.validate_atom_map(amap)

    def test_extra_solvent_shielding_rejected_not_silently_filtered(self):
        _, amap = initial_fixture()
        atoms = self.collect.validate_atom_map(amap)
        rows = [{"atom_index_zero_based": i, "element": r["element"], "isotropic_ppm": 1.0}
                for i, r in atoms.items() if r["include_in_analyte_nmr"]]
        solvent_c = next(i for i, r in atoms.items() if r["environment_atom"] and r["element"] == "C")
        rows.append({"atom_index_zero_based": solvent_c, "element": "C", "isotropic_ppm": 1.0})
        with self.assertRaises(ValueError):
            self.collect.extract_shieldings({"shielding": rows}, atoms)

    def test_missing_J_pair_rejected_not_filled_with_zero(self):
        _, amap = initial_fixture()
        atoms = self.collect.validate_atom_map(amap)
        pairs = sorted(self.collect.required_pairs(atoms))
        rows = [{"atom_i_zero_based": i, "atom_j_zero_based": j, "J_Hz": 0.0} for i, j in pairs[:-1]]
        with self.assertRaises(ValueError):
            self.collect.extract_couplings({"j_pairs": rows}, {"expected_j_pairs_zero_based": pairs}, atoms)

    def test_HF_and_CF_expected_by_collection(self):
        _, amap = initial_fixture("F4_D04_rep_R__DMSO")
        atoms = self.collect.validate_atom_map(amap)
        pairs = self.collect.required_pairs(atoms)
        kinds = {tuple(sorted((atoms[i]["element"], atoms[j]["element"]))) for i, j in pairs}
        self.assertIn(("F", "H"), kinds)
        self.assertIn(("C", "F"), kinds)

    def test_methyl_averaging_never_merges_CH2_hydrogens(self):
        _, amap = initial_fixture("F2_E09_rep_R__DMSO")
        atoms = self.collect.validate_atom_map(amap)
        sigma = {i: float(i*i) for i, a in atoms.items() if a["include_in_analyte_nmr"]}
        js = {pair: float(sum(pair)) for pair in self.collect.required_pairs(atoms)}
        averaged, _, groups = self.collect.methyl_average(atoms, sigma, js)
        self.assertTrue(groups)
        methyl_h = {i for group in groups for i in group}
        ch2_h = set()
        for i, row in atoms.items():
            if row["environment_atom"] or row["element"] != "C":
                continue
            hs = [j for j in row["neighbours_zero_based"] if atoms[j]["element"] == "H"]
            if len(hs) == 2:
                ch2_h.update(hs)
        self.assertTrue(ch2_h)
        self.assertFalse(ch2_h & methyl_h)
        for i in ch2_h:
            self.assertEqual(averaged[i], sigma[i])

    def mocked_collection(self, production=False, wrong_reference=False):
        # A unit test of arithmetic/control flow only. No synthetic values are
        # written into a real ensemble or asserted to be calculated properties.
        c = self.collect
        _, amap = initial_fixture()
        _, tmap = initial_fixture("TMS_DMSO")
        atoms, tatoms = c.validate_atom_map(amap), c.validate_atom_map(tmap)
        scenario, ph = "DMSO_02_W01_trace_water_sensitivity", "test_protocol"
        target = {"species_id": "TEST_ANALYTE", "compound_id": "TEST_COMPOUND", "scenario_id": scenario,
                  "protocol_hash": ph, "solvent": "DMSO", "temperature_K": 298.15,
                  "sampling_convergence_established": False,
                  "conformers": [{"id": "t0", "shielding_job": "t0"}, {"id": "t1", "shielding_job": "t1"}]}
        if wrong_reference:
            target["reference_ensemble"] = f"ensembles/{ph}/TMS_DMSO/DMSO_00_W00_continuum/ensemble.json"
        ref = {**target, "species_id": "TMS_DMSO", "compound_id": "TMS",
               "conformers": [{"id": "r0", "shielding_job": "r0"}, {"id": "r1", "shielding_job": "r1"}]}
        candidate = {"ensemble": target, "atoms": atoms, "weights": [0.2, 0.8], "source": {"atom_map_sha256": "tmap"},
                     "raw_couplings": {}, "couplings": {}, "method_fingerprint": "j_test", "provenance": []}
        ref_evidence = {"atom_map_sha256": "rmap"}

        def load(root, path, ensemble, row, stage, settings, map_hash):
            choices = {"t0": (20.0, 100.0), "t1": (24.0, 120.0), "r0": (30.0, 180.0), "r1": (40.0, 200.0)}
            h, carbon = choices[row["id"]]
            selected = tatoms if ensemble["species_id"] == "TMS_DMSO" else atoms
            values = [{"atom_index_zero_based": i, "element": a["element"], "isotropic_ppm": h if a["element"] == "H" else carbon}
                      for i, a in selected.items() if a["include_in_analyte_nmr"]]
            return {"method_fingerprint": "same_shielding_test"}, {"shielding": values}, {"test_fixture_only": True}

        with mock.patch.object(c, "_context", return_value=({"temperature_K": 298.15}, ph)), \
             mock.patch("clusterflow.policy.compound_for_species", return_value="TEST_COMPOUND"), \
             mock.patch.object(c, "collect_candidate_couplings", return_value=candidate), \
             mock.patch.object(c, "_ensemble", return_value=(ref, tatoms, [0.8, 0.2], ref_evidence)), \
             mock.patch.object(c, "_load_job", side_effect=load):
            return c.collect_species(ROOT, "TEST_ANALYTE", scenario, exploratory=not production,
                                     benchmark_gate="not_a_real_gate.json" if production else None, write=False)

    def test_TMS_uses_its_own_complete_ensemble_weights(self):
        result = self.mocked_collection()
        self.assertEqual(result["metadata"]["reference_weights"], [0.8, 0.2])
        self.assertAlmostEqual(result["metadata"]["reference_shieldings_ppm"]["H"], 32.0)
        self.assertAlmostEqual(result["metadata"]["reference_shieldings_ppm"]["C"], 184.0)
        for row in result["atom_rows"]:
            self.assertAlmostEqual(row["shift_ppm"], 8.8 if row["element"] == "H" else 68.0)
        self.assertEqual(result["metadata"]["status"], "UNVALIDATED_EXPLORATORY_CANDIDATE")

    def test_production_requires_convergence_even_after_J_approval(self):
        with mock.patch("clusterflow.benchmark.validate_benchmark_gate", return_value={"unit_test_stub": True}):
            with self.assertRaisesRegex(ValueError, "sampling|convergence"):
                self.mocked_collection(production=True)

    def test_cross_composition_TMS_reference_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "cross-composition|matching TMS"):
            self.mocked_collection(wrong_reference=True)


class IndependentFrequencyGate(unittest.TestCase):
    def test_incomplete_frequency_dimension_is_rejected(self):
        from clusterflow import gates, common, geometry
        settings = readj(ROOT / "config/protocol.json")
        value, amap = initial_fixture()
        with tempfile.TemporaryDirectory(prefix="synthetic_gate_unit_test_") as td:
            path = Path(td)
            (path / "job.inp").write_text("SYNTHETIC TEST FIXTURE; NOT ORCA INPUT\n", encoding="utf-8")
            (path / "atom_map.json").write_text(json.dumps(amap), encoding="utf-8")
            (path / "input_structure.mol").write_text(Chem.MolToMolBlock(value), encoding="utf-8")
            (path / "job.xyz").write_text(geometry.xyz_text(value, "SYNTHETIC_TEST_NOT_QM_RESULT"), encoding="utf-8")
            meta = {"stage": "optfreq", "protocol_hash": common.hash_protocol(settings),
                    "input_sha256": common.sha(path / "job.inp"), "atom_map_sha256": common.sha(path / "atom_map.json")}
            (path / "job.json").write_text(json.dumps(meta), encoding="utf-8")
            (path / "job.out").write_text("Program Version 6.1.0\nSYNTHETIC TEST OUTPUT\n", encoding="utf-8")
            (path / "receipt.json").write_text(json.dumps({"status": "completed", "synthetic_test_fixture": True}), encoding="utf-8")
            parsed = {"normal_termination": True, "optimization_converged": True, "frequencies_cm1": [0.0, 0.0, 100.0],
                      "gibbs_free_energy_Eh": -100.0, "temperature_K": settings["temperature_K"],
                      "qrrho_enabled": True, "qrrho_reference_frequency_cm1": settings["qrrho_reference_frequency_cm1"]}
            with mock.patch.object(gates, "receipt_errors", return_value=[]), mock.patch.object(gates, "parse_output", return_value=parsed):
                result = gates.inspect_job(path, settings)
            self.assertFalse(result["accepted"], "A 3-number frequency block for a multiatom cluster must never pass")
            self.assertTrue(any("frequen" in x.lower() for x in result["errors"]))


class IndependentPipelineAndDuplicates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from clusterflow import pipeline
        cls.pipeline = pipeline

    def test_initial_seed_without_GOAT_receipt_cannot_be_reviewed_as_result(self):
        jobs = readj(ROOT / "registry/jobs.json")
        job = next(j for j in jobs if j["species_id"] == "F2_E01_rep_R__DMSO" and j["scenario_id"] == "DMSO_02_W01_trace_water_sensitivity" and j["start_index"] == 7)
        with self.assertRaises((ValueError, OSError)):
            self.pipeline.draft_goat_review(ROOT, job["job_id"])

    def test_prescreen_without_reviewed_import_cannot_promote_initial_cluster(self):
        with self.assertRaises((ValueError, OSError)):
            self.pipeline.prepare_prescreen(ROOT, "F2_E01_rep_R__DMSO", "DMSO_02_W01_trace_water_sensitivity")

    def test_two_identical_DMSO_component_exchange_is_duplicate(self):
        value, amap = initial_fixture()
        groups = {}
        for row in amap:
            if row["component"].startswith("DMSO_d6"):
                groups.setdefault(row["component"], []).append(row["index_zero_based"])
        first, second = list(groups.values())
        changed = Chem.Mol(value)
        conf, original = changed.GetConformer(), value.GetConformer().GetPositions()
        for i, j in zip(first, second):
            conf.SetAtomPosition(i, tuple(original[j]))
            conf.SetAtomPosition(j, tuple(original[i]))
        self.assertLess(self.pipeline.cluster_rms(value, changed, amap), 1e-5)

    def test_cluster_rigid_rotation_and_translation_is_duplicate(self):
        value, amap = initial_fixture()
        changed = Chem.Mol(value)
        conf = changed.GetConformer()
        for i, (x, y, z) in enumerate(value.GetConformer().GetPositions()):
            conf.SetAtomPosition(i, (-y + 10, x - 20, z + 5))
        self.assertLess(self.pipeline.cluster_rms(value, changed, amap), 1e-5)

    def test_distinct_solvent_motif_not_merged_by_solute_only_RMSD(self):
        value, amap = initial_fixture()
        changed = Chem.Mol(value)
        conf = changed.GetConformer()
        component = next(r["component"] for r in amap if r["environment_atom"])
        for row in amap:
            if row["component"] == component:
                i = row["index_zero_based"]
                conf.SetAtomPosition(i, tuple(np.array(conf.GetAtomPosition(i)) + [1.5, 0.0, 0.0]))
        self.assertGreater(self.pipeline.cluster_rms(value, changed, amap), 0.1)

    def test_dedup_does_not_double_count_swapped_solvent_basins(self):
        value, amap = initial_fixture()
        settings = readj(ROOT / "config/protocol.json")
        entries = [{"mol": value, "job": {"conformer_id": "first"}, "parsed": {"gibbs_free_energy_Eh": -100.0}},
                   {"mol": Chem.Mol(value), "job": {"conformer_id": "duplicate"}, "parsed": {"gibbs_free_energy_Eh": -100.0}}]
        unique, duplicates = self.pipeline.deduplicate(entries, settings, amap)
        self.assertEqual(len(unique), 1)
        self.assertEqual(len(duplicates), 1)


class IndependentBenchmarkAndConvergence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from clusterflow import benchmark, convergence
        cls.benchmark, cls.convergence = benchmark, convergence

    def test_software_cannot_supply_experimental_panel_approval(self):
        for name in ("codex", "agent", "automatic", "system"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.benchmark._review({"decision": "APPROVED_INDEPENDENT_EXPERIMENTAL_PANEL", "reviewer": name,
                                        "reviewed_at": "2026-10-13T00:00:00Z"}, "APPROVED_INDEPENDENT_EXPERIMENTAL_PANEL")

    def test_known_J_sign_is_not_replaced_by_magnitude(self):
        policy = {"criterion": "lowest_MAE_Hz", "max_MAE_Hz": 2, "max_abs_error_Hz": 2,
                  "min_sign_accuracy": 1, "min_compounds": 2, "min_pairs": 2}
        rows = [{"compound_id": "A", "experimental_J_Hz": -10, "sign_known": True,
                 "predicted_J_Hz": {"PBE_pcJ3": -9, "PBE0_pcJ2": 10}},
                {"compound_id": "B", "experimental_J_Hz": 5, "sign_known": True,
                 "predicted_J_Hz": {"PBE_pcJ3": 5, "PBE0_pcJ2": -5}}]
        comparison, recommendation = self.benchmark.evaluate_comparison(rows, policy)
        self.assertEqual(recommendation, ["PBE_pcJ3"])
        self.assertEqual(comparison["PBE0_pcJ2"]["sign_accuracy"], 0)

    def test_unsigned_experimental_J_does_not_fabricate_sign_accuracy(self):
        policy = {"criterion": "lowest_MAE_Hz", "max_MAE_Hz": 1, "max_abs_error_Hz": 1,
                  "min_sign_accuracy": 0, "min_compounds": 2, "min_pairs": 2}
        rows = [{"compound_id": compound, "experimental_J_Hz": 10, "sign_known": False,
                 "predicted_J_Hz": {"PBE_pcJ3": -10, "PBE0_pcJ2": 10}} for compound in ("A", "B")]
        comparison, recommendation = self.benchmark.evaluate_comparison(rows, policy)
        self.assertEqual(set(recommendation), {"PBE_pcJ3", "PBE0_pcJ2"})
        self.assertIsNone(comparison["PBE_pcJ3"]["sign_accuracy"])

    def test_synthetic_artifact_cannot_be_relabelled_experimental(self):
        with tempfile.TemporaryDirectory(prefix="synthetic_benchmark_unit_test_") as td:
            root = Path(td)
            path = root / "claimed_experiment.json"
            path.write_text(json.dumps({"data_origin": "CALCULATED_SYNTHETIC_TEST_FIXTURE"}), encoding="utf-8")
            evidence = [{"path": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                         "source_type": "acquired_experimental", "source_description": "unit test only"}]
            with self.assertRaises(ValueError):
                self.benchmark._evidence(root, evidence)

    def test_simulated_panel_rejected_before_any_calculated_read(self):
        with tempfile.TemporaryDirectory(prefix="synthetic_panel_unit_test_") as td:
            root = Path(td)
            panel = root / "panel.json"
            panel.write_text(json.dumps({"data_origin": "experimental", "simulated": True}), encoding="utf-8")
            with mock.patch.object(self.benchmark, "collect_candidate_couplings") as calculations:
                with self.assertRaises(ValueError):
                    self.benchmark.run_benchmark(root, "panel.json", ["TARGET"], write=False)
                calculations.assert_not_called()

    def test_convergence_review_boolean_without_evidence_is_rejected(self):
        with tempfile.TemporaryDirectory(prefix="synthetic_convergence_unit_test_") as td:
            root = Path(td)
            gate = {"status": "REVIEWED_CONVERGENCE", "reviewer": "Test Human", "reviewed_utc": "2026-10-13T00:00:00Z",
                    "reason": "synthetic test fixture", "checks": {k: True for k in
                        ("basin_rediscovery", "no_new_important_basin", "motif_coverage", "peak_order_stable")},
                    "sampling_converged": True, "review_evidence": [], "waves": []}
            (root / "gate.json").write_text(json.dumps(gate), encoding="utf-8")
            with self.assertRaises(ValueError):
                self.convergence.validate_gate(root, "gate.json", "TARGET", "DMSO_02_W01_trace_water_sensitivity", "PBE0_pcJ2")

    def test_duplicated_sampling_IDs_rejected(self):
        jobs = readj(ROOT / "registry/jobs.json")
        job = next(j for j in jobs if j["species_id"] == "F2_E01_rep_R__DMSO")
        with self.assertRaisesRegex(ValueError, "Repeated|duplicate"):
            self.convergence._completed_sampling(ROOT, [job["job_id"], job["job_id"]], job["species_id"], job["scenario_id"])

    def convergence_fixture(self, root, duplicate_result=False, detached_sampling=False):
        def save(name, obj):
            path = root / name
            path.write_text(json.dumps(obj, sort_keys=True), encoding="utf-8")
            return {"path": name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}

        evidence = save("review_evidence.json", {"synthetic_unit_test": True})
        results, waves = [], []
        for number in (1, 2):
            target_ids = [f"t{i}" for i in range(1, number + 1)]
            reference_ids = [f"r{i}" for i in range(1, number + 1)]
            snapshot = {}
            for role, ids in (("target", target_ids), ("reference", reference_ids)):
                snapshot[role] = {"import_index_snapshot": {"sources": [{"goat_job_id": j} for j in ids],
                                                              "conformers": [{"goat_job_id": j} for j in ids]}}
            record = save(f"result{number}.json", {"metadata": {"ensemble_snapshot": snapshot}, "synthetic_unit_test_wave": number})
            results.append(record)
            waves.append({"result": record, "released_goat_jobs": target_ids + reference_ids})
        if duplicate_result:
            waves[1]["result"] = waves[0]["result"]
        if detached_sampling:
            # Completion of unrelated, unimported starts cannot be attached to
            # an unchanged calculated result as evidence of an expanded wave.
            waves[1]["released_goat_jobs"].append("t_unrepresented")
        gate = {"status": "REVIEWED_CONVERGENCE", "reviewer": "Test Human", "reviewed_utc": "2026-10-13T00:00:00Z",
                "reason": "synthetic unit test only", "checks": {k: True for k in
                    ("basin_rediscovery", "no_new_important_basin", "motif_coverage", "peak_order_stable")},
                "review_evidence": [evidence], "waves": waves}
        save("gate.json", gate)
        return gate

    def test_same_prediction_twice_cannot_establish_convergence(self):
        with tempfile.TemporaryDirectory(prefix="duplicate_snapshot_unit_test_") as td:
            root = Path(td)
            self.convergence_fixture(root, duplicate_result=True)
            with mock.patch.object(self.convergence, "_completed_sampling", return_value={"TARGET": 5, "TMS_DMSO": 5}), \
                 mock.patch.object(self.convergence, "_recompute_snapshot", return_value=({0: 1.0}, {0: {"element": "H"}})):
                with self.assertRaisesRegex(ValueError, "Repeated prediction|duplicate|same.*prediction"):
                    self.convergence.validate_gate(root, "gate.json", "TARGET", "DMSO_02_W01_trace_water_sensitivity", "PBE0_pcJ2")

    def test_completed_start_not_in_snapshot_is_not_convergence_evidence(self):
        with tempfile.TemporaryDirectory(prefix="unbound_wave_unit_test_") as td:
            root = Path(td)
            self.convergence_fixture(root, detached_sampling=True)
            with mock.patch.object(self.convergence, "_completed_sampling", return_value={"TARGET": 5, "TMS_DMSO": 5}), \
                 mock.patch.object(self.convergence, "_recompute_snapshot", return_value=({0: 1.0}, {0: {"element": "H"}})):
                with self.assertRaisesRegex(ValueError, "not exactly bound|source|import selection"):
                    self.convergence.validate_gate(root, "gate.json", "TARGET", "DMSO_02_W01_trace_water_sensitivity", "PBE0_pcJ2")

    def test_embedded_ensemble_cannot_differ_from_frozen_hash_bound_source(self):
        from clusterflow.common import hash_protocol
        with tempfile.TemporaryDirectory(prefix="modified_snapshot_unit_test_") as td:
            root = Path(td)
            (root / "config").mkdir()
            settings = readj(ROOT / "config/protocol.json")
            (root / "config/protocol.json").write_text(json.dumps(settings), encoding="utf-8")
            frozen = root / "frozen.json"
            frozen.write_text(json.dumps({"conformers": ["full_original"]}), encoding="utf-8")
            result = {"metadata": {"species_id": "TARGET", "scenario_id": "DMSO_02_W01_trace_water_sensitivity",
                        "j_method": "PBE0_pcJ2", "protocol_hash": hash_protocol(settings),
                        "ensemble_snapshot": {"target": {"conformers": ["arbitrarily_reduced"]}},
                        "frozen_ensembles": {"target": {"path": frozen.name, "sha256": hashlib.sha256(frozen.read_bytes()).hexdigest()}}}}
            with self.assertRaisesRegex(ValueError, "Embedded ensemble|frozen"):
                self.convergence._recompute_snapshot(root, result, "TARGET", "DMSO_02_W01_trace_water_sensitivity", "PBE0_pcJ2")

    def test_software_cannot_approve_convergence(self):
        with tempfile.TemporaryDirectory(prefix="software_review_unit_test_") as td:
            root = Path(td)
            gate = self.convergence_fixture(root)
            gate["reviewer"] = "codex"
            (root / "gate.json").write_text(json.dumps(gate), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Software|software"):
                self.convergence.validate_gate(root, "gate.json", "TARGET", "DMSO_02_W01_trace_water_sensitivity", "PBE0_pcJ2")


class IndependentOutputParsing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from clusterflow import orcaio
        cls.orcaio = orcaio

    def test_empty_output_has_no_success_or_fabricated_property(self):
        parsed = self.orcaio.parse_output("")
        self.assertFalse(parsed["normal_termination"])
        self.assertIsNone(parsed["final_energy_Eh"])
        self.assertIsNone(parsed["gibbs_free_energy_Eh"])
        self.assertEqual(parsed["shielding"], [])

    def test_appended_failed_run_cannot_inherit_old_success(self):
        text = "* O R C A *\nORCA TERMINATED NORMALLY\n* O R C A *\nSCF failed\n"
        self.assertFalse(self.orcaio.parse_output(text)["normal_termination"])

    def test_later_truncated_energy_does_not_use_earlier_energy(self):
        text = "FINAL SINGLE POINT ENERGY -40.0\nFINAL SINGLE POINT ENERGY\n"
        self.assertIsNone(self.orcaio.parse_output(text)["final_energy_Eh"])

    def test_duplicate_shielding_index_is_withheld(self):
        text = "CHEMICAL SHIELDING SUMMARY (ppm)\n 0 C 100.0 1.0\n 0 C 110.0 1.0\n\n"
        self.assertEqual(self.orcaio.parse_output(text)["shielding"], [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
