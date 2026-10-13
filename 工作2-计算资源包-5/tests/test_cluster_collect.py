"""Synthetic UNIT TEST fixtures only; never QM predictions or experimental data.

Temporary fake ORCA text exercises parsing and provenance. Scientific geometry
gates are explicitly mocked in collection integration tests; no test approval
or fixture is written into the package's production registry.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from clusterflow import benchmark, collect
from clusterflow.common import hash_protocol


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True), encoding="utf-8")


def atom(index, element, neighbours, environment=False, isotope=0):
    return {"index_zero_based": index, "element": element, "neighbours_zero_based": neighbours,
            "component": "environment_UNIT_TEST" if environment else "analyte", "environment_atom": environment,
            "isotope": isotope, "exchangeable_H": False,
            "include_in_analyte_nmr": not environment and element in {"H", "C"}}


def target_atoms():
    # CH2FO test graph; 3/4 are F/O and must not enter a proton/carbon shielding request.
    return [atom(0, "C", [1, 2, 3, 4]), atom(1, "H", [0]), atom(2, "H", [0]),
            atom(3, "F", [0]), atom(4, "O", [0], isotope=18), atom(5, "O", [6], True), atom(6, "H", [5], True, 2)]


def reference_atoms():
    rows = [atom(0, "Si", [1, 2, 3, 4])]
    rows += [atom(i, "C", [0] + list(range(5 + (i-1)*3, 8 + (i-1)*3))) for i in range(1, 5)]
    rows += [atom(i, "H", [1 + (i-5)//3]) for i in range(5, 17)]
    return rows + [atom(17, "O", [18], True), atom(18, "H", [17], True, 2)]


def shield_text(rows, offset=0, reference=False):
    lines = ["Program Version 6.1.1", "CHEMICAL SHIELDING SUMMARY (ppm)", "Nucleus Element Isotropic Anisotropy"]
    for row in rows:
        if not row["include_in_analyte_nmr"]:
            continue
        base = (200 if row["element"] == "C" else 30) if reference else (100 if row["element"] == "C" else 10 + row["index_zero_based"])
        lines.append(f"{row['index_zero_based']} {row['element']} {base + offset} 0.0")
    return "\n".join(lines) + "\n\n****ORCA TERMINATED NORMALLY****\n"


def j_property(pairs, offset=0):
    text = f'$SCF_Spin_Spin_Coupling\n&pairsInfo [&Type "ArrayOfIntegers", &Dim ({len(pairs)},4)]\n'
    text += "\n".join(f"{n} {i} 1 {j} 1" for n, (i, j) in enumerate(pairs))
    text += f'\n&pairsTotalSSCIso [&Type "ArrayOfDoubles", &Dim ({len(pairs)},1)]\n'
    text += "\n".join(f"{n} {n+1+offset}" for n in range(len(pairs)))
    return text + "\n$End\n"


class ClusterCollectionIntegration(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.settings = {"temperature_K": 298.15, "orca_version": "6.1.1"}
        self.ph = hash_protocol(self.settings)
        save(self.root / "config/protocol.json", self.settings)
        self.scenario = "DMSO_02_W00_primary"
        self.species = "UNIT_TEST_TARGET"
        self.target = self.make_ensemble(self.species, target_atoms(), [0.25, 0.75], [0, 4])
        self.reference = self.make_ensemble("TMS_DMSO", reference_atoms(), [0.8, 0.2], [0, 20])
        p1 = mock.patch("clusterflow.pipeline.validate_ensemble_sources", return_value={"accepted": True})
        p2 = mock.patch("clusterflow.gates.inspect_job", return_value={"accepted": True})
        p1.start(); p2.start()
        self.addCleanup(p1.stop); self.addCleanup(p2.stop)

    def seal(self, folder):
        names = ("job.json", "job.inp", "atom_map.json")
        save(folder / "receipt.json", {"status": "completed", "notice": "SYNTHETIC UNIT TEST RECEIPT NOT SCIENTIFIC EVIDENCE",
            "input_sha256": {n: collect.sha256_file(folder / n) for n in names},
            "output_sha256": {p.name: collect.sha256_file(p) for p in folder.glob("job.*") if p.name not in names}})

    def make_ensemble(self, species, atoms, weights, offsets):
        relative = f"maps/{species}/{self.scenario}/atom_map.json"
        save(self.root / relative, atoms)
        digest = collect.sha256_file(self.root / relative)
        e = {"species_id": species, "compound_id": "TMS" if species == "TMS_DMSO" else "UNIT_TEST_COMPOUND",
             "scenario_id": self.scenario, "solvent": "DMSO", "protocol_hash": self.ph, "temperature_K": 298.15,
             "atom_map": relative, "conformers": [], "sampling_convergence_established": False}
        reference = species == "TMS_DMSO"
        pairs = sorted(collect.required_pairs(collect.validate_atom_map(atoms)))
        for n, (weight, offset) in enumerate(zip(weights, offsets)):
            row = {"id": f"c{n}", "weight": weight, "G_Eh": -100+n, "j_jobs": {}}
            stages = ["shielding"] if reference else ["shielding", "j_PBE_pcJ3", "j_PBE0_pcJ2"]
            for stage in stages:
                folder = self.root / f"jobs/{self.ph}/{species}/{self.scenario}/{stage}/c{n}"
                meta = {"species_id": species, "scenario_id": self.scenario, "conformer_id": f"c{n}", "solvent": "DMSO",
                    "protocol_hash": self.ph, "stage": stage, "method_fingerprint": stage + "_UNIT_TEST",
                    "atom_map_sha256": digest, "expected_j_pairs_zero_based": pairs}
                save(folder / "job.json", meta)
                save(folder / "atom_map.json", atoms)
                (folder / "job.inp").write_text("SYNTHETIC UNIT TEST; NOT RUN", encoding="utf-8")
                (folder / "job.out").write_text(shield_text(atoms, offset, reference) if stage == "shielding" else
                    "Program Version 6.1.1\n****ORCA TERMINATED NORMALLY****\n", encoding="utf-8")
                if stage != "shielding":
                    (folder / "job.property.txt").write_text(j_property(pairs, offset), encoding="utf-8")
                    row["j_jobs"][stage[2:]] = folder.relative_to(self.root).as_posix()
                else:
                    row["shielding_job"] = folder.relative_to(self.root).as_posix()
                self.seal(folder)
            e["conformers"].append(row)
        save(self.root / f"ensembles/{self.ph}/{species}/{self.scenario}/ensemble.json", e)
        return e

    def run_collect(self, **kwargs):
        return collect.collect_species(self.root, self.species, self.scenario, **{"exploratory": True, "write": False, **kwargs})

    def test_complete_reference_has_own_weights_and_CH2_is_not_averaged(self):
        result = self.run_collect()
        values = {r["atom_index_zero_based"]: r for r in result["atom_rows"]}
        self.assertAlmostEqual(result["metadata"]["reference_shieldings_ppm"]["C"], 204)
        self.assertAlmostEqual(values[0]["shift_ppm"], 101)
        self.assertAlmostEqual(values[1]["shift_ppm"], 20)
        self.assertAlmostEqual(values[2]["shift_ppm"], 19)
        self.assertEqual({r["atom_index_zero_based"] for r in result["spectrum_atom_rows"]}, {0, 1, 2, 3})
        self.assertTrue(all(r["atom_i_zero_based"] < 5 and r["atom_j_zero_based"] < 5 for r in result["j_pairs"]))
        self.assertEqual(result["metadata"]["status"], "UNVALIDATED_EXPLORATORY_CANDIDATE")

    def test_nonoverwrite_timestamp_outputs_and_two_distinct_methods(self):
        a = self.run_collect(write=True)
        b = self.run_collect(write=True, j_method="PBE_pcJ3")
        c = self.run_collect(write=True)
        self.assertEqual(len({a["result_path"], b["result_path"], c["result_path"]}), 3)
        self.assertTrue(all((self.root / x["result_path"]).is_file() for x in (a, b, c)))
        self.assertIn("input_sha256", a["metadata"]["provenance"][0])
        for label, frozen in a["metadata"]["frozen_ensembles"].items():
            source = a["metadata"]["source" if label == "target" else "reference_source"]
            self.assertEqual((self.root / frozen["path"]).read_bytes(), (self.root / source["ensemble_path"]).read_bytes())
            self.assertEqual(collect.sha256_file(self.root / frozen["path"]), frozen["sha256"])

    def test_production_requires_real_gates(self):
        with self.assertRaisesRegex(ValueError, "benchmark gate"):
            self.run_collect(exploratory=False)
        with self.assertRaisesRegex(ValueError, "convergence gate"):
            self.run_collect(exploratory=False, benchmark_gate="not_a_real_gate.json")

    def test_changed_output_receipt_rejected(self):
        p = self.root / self.target["conformers"][0]["shielding_job"] / "job.out"
        p.write_text(p.read_text() + "changed", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "artifact changed"):
            self.run_collect()

    def test_no_missing_J_filled_with_zero(self):
        folder = self.root / self.target["conformers"][0]["j_jobs"]["PBE0_pcJ2"]
        pairs = sorted(collect.required_pairs(collect.validate_atom_map(target_atoms())))
        (folder / "job.property.txt").write_text(j_property(pairs[:-1]), encoding="utf-8")
        self.seal(folder)
        with self.assertRaisesRegex(ValueError, "Missing J pairs"):
            self.run_collect()

    def test_cross_composition_reference_rejected(self):
        self.target["reference_ensemble"] = f"ensembles/{self.ph}/TMS_DMSO/DMSO_00_W00_continuum/ensemble.json"
        save(self.root / f"ensembles/{self.ph}/{self.species}/{self.scenario}/ensemble.json", self.target)
        with self.assertRaisesRegex(ValueError, "cross-composition"):
            self.run_collect()

    def test_TMS_different_fingerprint_rejected(self):
        for row in self.reference["conformers"]:
            folder = self.root / row["shielding_job"]
            meta = collect.read_json(folder / "job.json")
            meta["method_fingerprint"] = "OTHER_METHOD"
            save(folder / "job.json", meta)
            self.seal(folder)
        with self.assertRaisesRegex(ValueError, "TMS method"):
            self.run_collect()

    def test_weights_must_not_be_silently_renormalized(self):
        self.reference["conformers"][0]["weight"] = 0.7
        save(self.root / f"ensembles/{self.ph}/TMS_DMSO/{self.scenario}/ensemble.json", self.reference)
        with self.assertRaisesRegex(ValueError, "normalized"):
            self.run_collect()


class ExperimentalBenchmarkUnitTests(unittest.TestCase):
    """Mocked numerical comparison only; fixtures are NOT genuine measurements."""
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.settings = {"temperature_K": 298.15, "orca_version": "6.1.1"}
        self.ph = hash_protocol(self.settings)
        save(self.root / "config/protocol.json", self.settings)
        self.policy = {"criterion": "lowest_MAE_Hz", "max_MAE_Hz": 2.0, "max_abs_error_Hz": 3.0,
                       "min_sign_accuracy": 1.0, "min_compounds": 2, "min_pairs": 2}
        review = lambda decision: {"reviewer": "UNIT TEST REVIEWER (FICTIONAL)", "reviewed_at": "2026-01-01", "decision": decision}
        evidence = self.root / "experimental_fixture/raw.txt"
        evidence.parent.mkdir(parents=True)
        evidence.write_text("SYNTHETIC ENGINEERING UNIT TEST; NOT EXPERIMENTAL DATA", encoding="utf-8")
        self.panel = {"data_origin": "experimental", "protocol_hash": self.ph, "temperature_K": 298.15,
            "solvent": "DMSO", "scenario_id": "test_scenario", "chemical_domain": "UNIT TEST ONLY",
            "independence_statement": "MOCKED ENGINEERING TEST: NOT A SCIENTIFIC CLAIM", "predeclared_selection": self.policy,
            "review": review("APPROVED_INDEPENDENT_EXPERIMENTAL_PANEL"), "cases": []}
        for k in (1, 2):
            self.panel["cases"].append({"case_id": f"case{k}", "compound_id": f"independent{k}", "species_id": f"independent{k}",
                "scenario_id": "test_scenario", "temperature_K": 298.15,
                "experimental_evidence": [{"path": evidence.relative_to(self.root).as_posix(), "sha256": collect.sha256_file(evidence),
                    "source_type": "acquired_experimental", "source_description": "FICTIONAL TEST ATTESTATION"}],
                "assignment_review": review("APPROVED_ATOM_PAIR_ASSIGNMENT"),
                "measured_pairs": [{"atoms_zero_based": [1, 2], "J_Hz": 8.0, "sign_known": True, "evidence_index": 0, "source_locator": "TEST row 1"}]})
        self.save_panel()
        self.prediction = mock.patch("clusterflow.benchmark.collect_candidate_couplings", side_effect=self.candidate)
        self.prediction.start(); self.addCleanup(self.prediction.stop)
        # Nominal comparison tests model an authentic source acceptance; these
        # explicitly synthetic artifacts must never pass the real source gate.
        self.real_evidence = benchmark._evidence
        self.source_gate = mock.patch("clusterflow.benchmark._evidence", side_effect=self.mock_evidence)
        self.source_gate.start(); self.addCleanup(self.source_gate.stop)

    def mock_evidence(self, root, records):
        for entry in records:
            if collect.sha256_file(root / entry["path"]) != entry["sha256"]:
                raise ValueError("Experimental evidence artifact is missing or changed")
        return copy.deepcopy(records)

    def save_panel(self):
        save(self.root / "experimental_fixture/panel.json", self.panel)

    def candidate(self, root, species, scenario, method):
        return {"ensemble": {"compound_id": species, "protocol_hash": self.ph, "temperature_K": 298.15,
                             "conformers": [{"id": "UNIT_TEST_CONFORMER"}]},
                "couplings": {(1, 2): 8.2 if method == "PBE_pcJ3" else 9.5},
                "atoms": collect.validate_atom_map(target_atoms()), "source": {"test": "MOCKED_UNIT_TEST_ONLY"},
                "provenance": [{"notice": "NO_REAL_QUANTUM_OUTPUT"}], "weights": [1.0], "method_fingerprint": method}

    def run_comparison(self, targets=("future_target",), **kwargs):
        return benchmark.run_benchmark(self.root, "experimental_fixture/panel.json", targets, **{"write": False, **kwargs})

    def test_comparison_never_grants_approval(self):
        report = self.run_comparison()
        self.assertFalse(report["approval_granted"])
        self.assertEqual(report["status"], "EVALUATED_REQUIRES_SELECTION_REVIEW")
        self.assertEqual(report["recommended_methods"], ["PBE_pcJ3"])
        self.assertAlmostEqual(report["method_comparison"]["PBE_pcJ3"]["MAE_Hz"], 0.2)
        self.assertEqual(len(report["measurements"]), 2)

    def test_simulated_declared_panel_rejected(self):
        self.panel["data_origin"] = "simulated"
        self.save_panel()
        with self.assertRaisesRegex(ValueError, "experimental"):
            self.run_comparison()

    def test_actual_evidence_gate_rejects_our_synthetic_unit_fixture(self):
        with self.assertRaisesRegex(ValueError, "synthetic"):
            self.real_evidence(self.root, self.panel["cases"][0]["experimental_evidence"])

    def test_target_self_validation_rejected(self):
        with self.assertRaisesRegex(ValueError, "independent"):
            self.run_comparison(targets=["independent1"])

    def test_absent_review_rejected(self):
        self.panel.pop("review")
        self.save_panel()
        with self.assertRaisesRegex(ValueError, "human review"):
            self.run_comparison()

    def test_different_temperature_rejected(self):
        self.panel["cases"][0]["temperature_K"] = 310.0
        self.save_panel()
        with self.assertRaisesRegex(ValueError, "temperature mismatch"):
            self.run_comparison()

    def test_changed_experimental_artifact_rejected(self):
        (self.root / "experimental_fixture/raw.txt").write_text("modified", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "artifact is missing or changed"):
            self.run_comparison()

    def test_missing_pair_cannot_count_as_a_comparison(self):
        self.panel["cases"][0]["measured_pairs"][0]["atoms_zero_based"] = [5, 6]
        self.save_panel()
        with self.assertRaisesRegex(ValueError, "BOTH methods"):
            self.run_comparison()

    def test_unsigned_experiment_does_not_prove_sign_accuracy(self):
        measurements = [{"compound_id": str(k), "experimental_J_Hz": 8, "predicted_J_Hz": {m: -8 for m in collect.J_METHODS}, "sign_known": False} for k in (1, 2)]
        metrics, choices = benchmark.evaluate_comparison(measurements, self.policy)
        self.assertEqual(metrics["PBE_pcJ3"]["MAE_Hz"], 0)
        self.assertIsNone(metrics["PBE_pcJ3"]["sign_accuracy"])
        self.assertFalse(choices)

    def test_report_itself_is_not_a_production_gate(self):
        report = self.run_comparison(write=True)
        with self.assertRaisesRegex(ValueError, "human review"):
            benchmark.validate_benchmark_gate(self.root, report["result_path"], "future_target", "PBE_pcJ3", "test_scenario", 298.15, self.ph)

    def make_gate(self, report):
        gate = {"decision": "APPROVED_FOR_STATED_SCOPE", "reviewer": "FICTIONAL UNIT TEST REVIEWER",
                "reviewed_at": "2026-01-02", "selection_rationale": "MOCKED UNIT TEST ONLY",
                "report_path": report["result_path"], "report_sha256": collect.sha256_file(self.root / report["result_path"]),
                "selected_j_method": "PBE_pcJ3", "scope": report["scope"]}
        save(self.root / "experimental_fixture/gate.json", gate)
        return "experimental_fixture/gate.json"

    def test_explicit_mocked_review_gate_rechecks_both_methods(self):
        report = self.run_comparison(write=True)
        path = self.make_gate(report)
        gate = benchmark.validate_benchmark_gate(self.root, path, "future_target", "PBE_pcJ3", "test_scenario", 298.15,
                                                self.ph, target_pair_classes={"HH"}, method_fingerprint="PBE_pcJ3")
        self.assertEqual(gate["status"], "REVIEWED_J_METHOD_FOR_STATED_SCOPE_ONLY")
        with self.assertRaisesRegex(ValueError, "coupling classes"):
            benchmark.validate_benchmark_gate(self.root, path, "future_target", "PBE_pcJ3", "test_scenario", 298.15,
                                              self.ph, target_pair_classes={"CF"})
        with self.assertRaisesRegex(ValueError, "method/basis"):
            benchmark.validate_benchmark_gate(self.root, path, "future_target", "PBE_pcJ3", "test_scenario", 298.15,
                                              self.ph, method_fingerprint="Br_MIXED_BASIS")

    def test_forged_report_metrics_rejected_even_with_new_file_hash(self):
        report = self.run_comparison(write=True)
        report["method_comparison"]["PBE_pcJ3"]["MAE_Hz"] = 0.0
        save(self.root / report["result_path"], report)
        path = self.make_gate(report)
        with self.assertRaisesRegex(ValueError, "scientific content digest"):
            benchmark.validate_benchmark_gate(self.root, path, "future_target", "PBE_pcJ3", "test_scenario", 298.15, self.ph)


if __name__ == "__main__":
    unittest.main()
