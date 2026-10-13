"""SYNTHETIC engineering fixtures only: no experimental panel/QM produced.

Collection integration reuses temporary synthetic ORCA outputs from the
collector tests, whose geometry/provenance gates are explicitly mocked there.
These tests establish software behavior, never experimental NMR accuracy.
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).parent))
from clusterflow import collect, multiplet
from clusterflow.spectra import make_spectrum
import test_cluster_collect as collection_fixtures
from test_cluster_collect import save


class MultipletNumericTests(unittest.TestCase):
    def test_uniform_shift_preserves_pattern_not_absolute_agreement(self):
        a = [{"ppm": 1.0, "area": 1}, {"ppm": 1.02, "area": 1}]
        b = [{"ppm": 1.1, "area": 1}, {"ppm": 1.12, "area": 1}]
        result = multiplet.compare_lines(a, b, 400)
        self.assertAlmostEqual(result["absolute_position_rmse_ppm"], 0.1)
        self.assertAlmostEqual(result["centered_position_rmse_Hz"], 0.0)
        self.assertAlmostEqual(result["adjacent_spacing_rmse_Hz"], 0.0)

    def test_spacing_error_in_hz(self):
        a = [{"ppm": 1, "area": 1}, {"ppm": 1.02, "area": 1}]
        b = [{"ppm": 1, "area": 1}, {"ppm": 1.025, "area": 1}]
        result = multiplet.compare_lines(a, b, 400)
        self.assertAlmostEqual(result["adjacent_spacing_rmse_Hz"], 2)

    def test_integrated_area_ratios(self):
        a = [{"ppm": 1, "area": 2}, {"ppm": 2, "area": 4}, {"ppm": 3, "area": 2}]
        b = [{"ppm": 1, "area": 1}, {"ppm": 2, "area": 2}, {"ppm": 3, "area": 1}]
        self.assertEqual(multiplet.compare_lines(a, b, 400)["relative_line_area_mae"], 0)

    def test_missing_lines_never_silent_match(self):
        a = [{"ppm": 1, "area": 1}, {"ppm": 2, "area": 1}]
        b = [{"ppm": 1, "area": 1}]
        result = multiplet.compare_lines(a, b, 400)
        self.assertEqual(result["line_count_difference"], 1)
        self.assertIsNone(result["relative_line_area_mae"])
        self.assertIsNone(result["centered_position_rmse_Hz"])

    def test_empty_unresolved_no_false_pass(self):
        result = multiplet.compare_lines([], [{"ppm": 1, "area": 1}], 400)
        self.assertIsNone(result["absolute_position_rmse_ppm"])
        self.assertEqual(result["experimental_line_count"], 0)

    def test_resolution_group_centroid_and_area(self):
        result = multiplet.resolve_transitions([{"ppm": 1, "intensity": 1},
                {"ppm": 1.001, "intensity": 3}, {"ppm": 2, "intensity": 1}], 400, 1)
        self.assertEqual(len(result), 2)
        self.assertAlmostEqual(result[0]["ppm"], 1.00075)
        self.assertEqual(result[0]["area"], 4)

    def test_resolution_is_not_transitive_chain(self):
        lines = [{"ppm": 1+i*0.002, "intensity": 1} for i in range(3)]
        self.assertEqual(len(multiplet.resolve_transitions(lines, 400, 1)), 2)

    def test_intensity_floor_explicit(self):
        lines = [{"ppm": 1, "intensity": 1}, {"ppm": 2, "intensity": 0.001}]
        self.assertEqual(len(multiplet.resolve_transitions(lines, 400, 1, 0.01)), 1)

    def test_negative_area_invalid(self):
        with self.assertRaises(ValueError):
            multiplet.resolve_transitions([{"ppm": 1, "intensity": -1}], 400, 1)

    def test_invalid_floor_invalid(self):
        with self.assertRaises(ValueError):
            multiplet.resolve_transitions([], 400, 1, 1)

    def test_parameter_digest_binds_field(self):
        self.assertNotEqual(multiplet.parameter_digest([], [], "1H", 400),
                            multiplet.parameter_digest([], [], "1H", 600))

    def test_exact_AB_second_order_roofing_not_first_order_doublets(self):
        atoms = [{"atom_index_zero_based": 0, "element": "H", "shift_ppm": 4.0},
                 {"atom_index_zero_based": 1, "element": "H", "shift_ppm": 4.01}]
        model = make_spectrum(atoms, [{"atom_i_zero_based": 0, "atom_j_zero_based": 1, "J_Hz": 12}],
                              "1H", {"field_MHz_H": 400})
        lines = multiplet.resolve_transitions(model["sticks"], 400, 0.01, 0)
        splitting = math.sqrt(4**2 + 12**2)
        expected = [4.005-(splitting+12)/800, 4.005-(splitting-12)/800,
                    4.005+(splitting-12)/800, 4.005+(splitting+12)/800]
        self.assertEqual(len(lines), 4)
        for line, ppm in zip(lines, expected):
            self.assertAlmostEqual(line["ppm"], ppm, places=8)
        self.assertGreater(lines[1]["area"], 30*lines[0]["area"])


class MultipletWorkflowTests(unittest.TestCase):
    def setUp(self):
        # This helper activates explicitly mocked scientific gates, TEMP ONLY.
        self.fixture = collection_fixtures.ClusterCollectionIntegration()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        self.prediction = collect.collect_species(self.root, self.fixture.species,
                        self.fixture.scenario, "PBE_pcJ3", exploratory=True, write=True)
        self.prediction_path = self.root / self.prediction["result_path"]
        source = self.root / "SYNTHETIC_TEST_ONLY_NOT_EXPERIMENT.txt"
        source.write_text("SYNTHETIC ENGINEERING FIXTURE; NOT EXPERIMENTAL NMR", encoding="utf-8")
        model = make_spectrum(self.prediction["spectrum_atom_rows"], self.prediction["j_pairs"], "1H",
                              {"field_MHz_H": 400, "linewidth_Hz_H": 1})
        lines = multiplet.resolve_transitions(model["sticks"], 400, 0.1, 0.001)
        low, high = min(r["ppm"] for r in lines)-1, max(r["ppm"] for r in lines)+1
        self.panel = {"schema_version": "1.0", "data_kind": "EXPERIMENTAL_NMR",
            "engineering_test_notice": "THIS TEMP PANEL IS SYNTHETIC, NOT AN EXPERIMENTAL DATASET",
            "prediction": {"path": self.prediction["result_path"], "sha256": collect.sha256_file(self.prediction_path)},
            "experimental_artifacts": [{"path": source.name, "sha256": collect.sha256_file(source)}],
            "assignment_review": {"reviewer": "Fictional fixture reviewer", "date": "2000-01-01", "notes": "Synthetic unit-test schema only"},
            "experiment": {"sample_identity": self.fixture.species, "solvent": "DMSO-d6", "reference": "TMS",
               "processing": "Synthetic numerical fixture, no acquired data", "decoupling": "none", "field_MHz": 400,
               "temperature_K": 298.15, "nucleus": "1H"},
            "regions": [{"id": "unit_test_group", "assigned_atom_indices_zero_based": [1, 2],
                "experimental_window_ppm": [low, high], "overlap": False, "second_order_observed": "unknown",
                "multiplicity_description": "Synthetic numerical fixture", "lines": lines,
                "line_area_kind": "fitted_or_assigned_component_area", "line_area_quantitative": True,
                "resolution_Hz": 0.1, "relative_area_floor": 0.001,
                "integral": {"value": 2, "normalization_per_nucleus": 1, "quantitative": True,
                             "normalization_description": "Synthetic known area"}}]}
        self.panel_path = self.root / "UNIT_TEST_PANEL.json"

    def run_panel(self, write=False):
        save(self.panel_path, self.panel)
        return multiplet.validate_multiplets(self.root, self.panel_path, write=write)

    def test_report_is_not_automatic_approval(self):
        result = self.run_panel()
        self.assertEqual(result["status"], "METRICS_COMPLETE_AWAITING_HUMAN_REVIEW")
        self.assertEqual(result["regions"][0]["metrics"]["line_count_difference"], 0)
        self.assertNotIn("approved", result)
        self.assertTrue(any("first order" in x for x in result["regions"][0]["warnings"]))

    def test_write_new_report_and_no_overwrite(self):
        result = self.run_panel(write=True)
        self.assertTrue((self.root / result["report_path"]).is_file())
        with self.assertRaises(FileExistsError):
            multiplet.validate_multiplets(self.root, self.panel_path, output_path=result["report_path"])

    def test_missing_experimental_artifact_blocks(self):
        self.panel["experimental_artifacts"][0]["sha256"] = "0"*64
        with self.assertRaisesRegex(ValueError, "changed source"):
            self.run_panel()

    def test_changed_prediction_hash_blocks(self):
        self.panel["prediction"]["sha256"] = "0"*64
        with self.assertRaisesRegex(ValueError, "changed source"):
            self.run_panel()

    def test_modified_prediction_numbers_cannot_be_resealed(self):
        value = collect.read_json(self.prediction_path)
        value["atom_rows"][0]["shift_ppm"] += 0.2
        save(self.prediction_path, value)
        self.panel["prediction"]["sha256"] = collect.sha256_file(self.prediction_path)
        with self.assertRaisesRegex(ValueError, "re-collected"):
            self.run_panel()

    def test_changed_real_source_blocks(self):
        source = self.prediction["metadata"]["provenance"][0]["directory"]
        output = self.root / source / "job.out"
        output.write_text(output.read_text(encoding="utf-8")+"changed", encoding="utf-8")
        with self.assertRaises(ValueError):
            self.run_panel()

    def test_frozen_ensemble_changed_blocks(self):
        snap = self.root / self.prediction["metadata"]["frozen_ensembles"]["target"]["path"]
        snap.write_text("{}", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "changed source"):
            self.run_panel()

    def test_incomplete_assignments_block(self):
        self.panel["regions"][0]["assigned_atom_indices_zero_based"] = [1]
        with self.assertRaisesRegex(ValueError, "Every observed nucleus"):
            self.run_panel()

    def test_explicit_exclusion_reported(self):
        self.panel["regions"][0]["assigned_atom_indices_zero_based"] = [1]
        self.panel["excluded_atoms"] = [{"atom_index_zero_based": 2, "reason": "Unresolved in synthetic engineering fixture"}]
        self.assertEqual(self.run_panel()["coverage"]["excluded_nuclei"], 1)

    def test_peak_height_not_component_area(self):
        self.panel["regions"][0]["line_area_kind"] = "peak_height"
        with self.assertRaisesRegex(ValueError, "not peak heights"):
            self.run_panel()

    def test_unresolved_is_reported_not_passed(self):
        self.panel["regions"][0]["lines"] = []
        self.panel["regions"][0]["unresolved_reason"] = "Synthetic unresolved test"
        result = self.run_panel()
        self.assertTrue(any("remains incomplete" in x for x in result["regions"][0]["warnings"]))

    def test_large_error_is_recommendation_not_exception(self):
        self.panel["regions"][0]["integral"]["value"] = 10
        result = self.run_panel()
        self.assertIn("relative_region_integral_error", result["regions"][0]["metric_flags"])

    def test_nonquantitative_integral_not_scored(self):
        self.panel["regions"][0]["integral"]["quantitative"] = False
        self.assertIsNone(self.run_panel()["regions"][0]["metrics"]["relative_region_integral_error"])

    def test_nonquantitative_line_area_not_scored(self):
        self.panel["regions"][0]["line_area_quantitative"] = False
        self.assertIsNone(self.run_panel()["regions"][0]["metrics"]["relative_line_area_mae"])

    def test_software_cannot_assign_review(self):
        self.panel["assignment_review"]["reviewer"] = "Codex"
        with self.assertRaisesRegex(ValueError, "Software"):
            self.run_panel()

    def test_conditions_mismatch_is_advisory(self):
        self.panel["experiment"]["temperature_K"] = 320
        self.assertTrue(any("differs" in x for x in self.run_panel()["warnings"]))

    def test_field_invalid_blocks(self):
        self.panel["experiment"]["field_MHz"] = 0
        with self.assertRaises(ValueError):
            self.run_panel()

    def test_nested_settings_cannot_override_experimental_field(self):
        self.panel["spectrum_settings"] = {"field_MHz_H": 600, "spectrum": {"field_MHz_H": 800}}
        self.assertEqual(self.run_panel()["model"]["field_MHz"], 400)

    def test_nested_settings_cannot_inject_noise(self):
        self.panel["spectrum_settings"] = {"spectrum": {"noise": True}}
        with self.assertRaisesRegex(ValueError, "noise"):
            self.run_panel()

    def test_sample_identity_mismatch_blocks(self):
        self.panel["experiment"]["sample_identity"] = "ANOTHER_COMPOUND"
        with self.assertRaisesRegex(ValueError, "identity"):
            self.run_panel()

    def test_tolerance_unknown_key_blocks(self):
        self.panel["tolerances"] = {"misleading_tolerance": 100}
        with self.assertRaisesRegex(ValueError, "Unknown"):
            self.run_panel()

    def test_external_parameter_mismatch_blocks(self):
        external = self.root / "UNIT_TEST_EXTERNAL.json"
        save(external, {"parameter_digest": "bad", "review": self.panel["assignment_review"]})
        self.panel["external_transitions"] = {"path": external.name, "sha256": collect.sha256_file(external)}
        with self.assertRaisesRegex(ValueError, "not bound"):
            self.run_panel()

    def test_external_source_bound_output_is_separately_identified(self):
        model = make_spectrum(self.prediction["spectrum_atom_rows"], self.prediction["j_pairs"], "1H",
                              {"field_MHz_H": 400})
        value = {"parameter_digest": multiplet.parameter_digest(self.prediction["spectrum_atom_rows"],
                       self.prediction["j_pairs"], "1H", 400),
                 "prediction_sha256": self.panel["prediction"]["sha256"],
                 "review": self.panel["assignment_review"], "solver": "SYNTHETIC_UNIT_TEST_SOLVER",
                 "solver_version": "TEST_ONLY", "hamiltonian_description": "Numerical engineering fixture",
                 "spin_model_scope": "DECLARED_APPROXIMATION",
                 "limitations": "Not an actual external scientific calculation",
                 "solver_artifacts": self.panel["experimental_artifacts"],
                 "observed_atom_indices_zero_based": [1, 2], "sticks": model["sticks"]}
        external = self.root / "UNIT_TEST_EXTERNAL.json"
        save(external, value)
        self.panel["external_transitions"] = {"path": external.name, "sha256": collect.sha256_file(external)}
        result = self.run_panel()
        self.assertEqual(result["external_solver"]["solver"], "SYNTHETIC_UNIT_TEST_SOLVER")
        self.assertEqual(result["status"], "METRICS_COMPLETE_AWAITING_HUMAN_REVIEW")

    def test_full_F_model_needs_additional_parameter_sources(self):
        value = {"parameter_digest": multiplet.parameter_digest(self.prediction["spectrum_atom_rows"],
                       self.prediction["j_pairs"], "1H", 400),
                 "prediction_sha256": self.panel["prediction"]["sha256"],
                 "review": self.panel["assignment_review"], "solver": "SYNTHETIC_UNIT_TEST_SOLVER",
                 "solver_version": "TEST_ONLY", "hamiltonian_description": "Numerical engineering fixture",
                 "limitations": "Unit test only", "spin_model_scope": "REVIEWED_FULL_SPIN_MODEL"}
        external = self.root / "UNIT_TEST_EXTERNAL.json"
        save(external, value)
        self.panel["external_transitions"] = {"path": external.name, "sha256": collect.sha256_file(external)}
        with self.assertRaisesRegex(ValueError, "additional 19F"):
            self.run_panel()

    def test_input_source_cannot_escape_package(self):
        self.panel["experimental_artifacts"][0]["path"] = "../outside.txt"
        with self.assertRaisesRegex(ValueError, "escapes"):
            self.run_panel()


if __name__ == "__main__":
    unittest.main()
