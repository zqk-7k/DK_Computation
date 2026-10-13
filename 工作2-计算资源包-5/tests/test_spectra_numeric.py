"""Pure numerical synthetic engineering fixtures; NOT molecular predictions.

Uses the current protocol's spectrum settings. No image export, file generation,
experimental claim, ORCA invocation, or other quantum calculation is performed.
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from clusterflow.spectra import make_spectrum


def atom(index, element, shift=None, exchangeable=False):
    return {"atom_index_zero_based": index, "element": element, "shift_ppm": shift,
            "exchangeable_H": exchangeable, "source": "SYNTHETIC_NUMERICAL_UNIT_TEST_NOT_EXPERIMENTAL"}


def pair(left, right, value):
    return {"atom_i_zero_based": left, "atom_j_zero_based": right, "J_Hz": value}


def area(result):
    x, y = np.asarray(result["x_ppm"]), np.asarray(result["y"])
    return float(np.sum((y[:-1] + y[1:]) * np.diff(x) / 2))


class NumericalSpectrumTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.protocol = json.loads((ROOT / "config/protocol.json").read_text(encoding="utf-8-sig"))

    def test_AB_exact_positions_roofing_and_area_with_current_field(self):
        field = self.protocol["spectrum"]["proton_frequency_MHz"]
        delta_hz, coupling, first = 10.0, 8.0, 1.0
        second = first + delta_hz / field
        result = make_spectrum([atom(3, "H", first), atom(9, "H", second)], [pair(3, 9, coupling)], "1H", self.protocol)
        centre = (first + second) / 2
        expected = sorted(centre + sign * (math.hypot(delta_hz, coupling) + adjustment) / (2 * field)
                          for sign in (-1, 1) for adjustment in (-coupling, coupling))
        np.testing.assert_allclose([r["ppm"] for r in result["sticks"]], expected, atol=1e-9, rtol=0)
        self.assertEqual(result["method"], "exact H-H spin-1/2 block")
        self.assertGreater(result["sticks"][1]["intensity"], result["sticks"][0]["intensity"])
        self.assertGreater(result["sticks"][2]["intensity"], result["sticks"][3]["intensity"])
        self.assertAlmostEqual(sum(r["intensity"] for r in result["sticks"]), 2, places=10)
        self.assertAlmostEqual(area(result), 2, places=10)
        self.assertEqual(result["block_report"][0]["atom_indices_zero_based"], [3, 9])

    def test_C_F_doublet_spacing_and_area_at_current_carbon_field(self):
        field = self.protocol["spectrum"]["carbon_frequency_MHz"]
        result = make_spectrum([atom(2, "C", 110), atom(7, "F")], [pair(2, 7, 200)], "13C", self.protocol)
        np.testing.assert_allclose([r["ppm"] for r in result["sticks"]], [110-100/field, 110+100/field], atol=1e-9, rtol=0)
        np.testing.assert_allclose([r["intensity"] for r in result["sticks"]], [.5, .5], atol=1e-12)
        self.assertAlmostEqual(area(result), 1, places=10)
        self.assertIn("first-order C-F only", result["method"])
        self.assertTrue(any("C-F splitting is first order" in warning for warning in result["warnings"]))

    def test_equal_two_F_triplet_and_proton_decoupled_carbon(self):
        field = self.protocol["spectrum"]["carbon_frequency_MHz"]
        atoms = [atom(0, "C", 100), atom(1, "F"), atom(2, "F"), atom(3, "H", 4), atom(4, "C", 20)]
        couplings = [pair(0, 1, 200), pair(0, 2, 200), pair(0, 3, 150), pair(0, 4, 50)]
        result = make_spectrum(atoms, couplings, "13C", self.protocol)
        np.testing.assert_allclose([r["ppm"] for r in result["sticks"]], [20, 100-200/field, 100, 100+200/field], atol=1e-9, rtol=0)
        np.testing.assert_allclose([r["intensity"] for r in result["sticks"]], [1, .25, .5, .25], atol=1e-12)
        self.assertAlmostEqual(area(result), 2, places=10)
        self.assertEqual(result["n_observed_nuclei"], 2)

    def test_current_protocol_metadata_and_source_indices_preserved(self):
        atoms, pairs = [atom(11, "H", 1), atom(23, "H", 1.03)], [pair(11, 23, 7)]
        result = make_spectrum(atoms, pairs, "1H", self.protocol)
        self.assertEqual(result["data_status"], "CALCULATED_MODEL_NOT_EXPERIMENTAL")
        self.assertEqual(result["nucleus"], "1H")
        self.assertEqual(result["field_MHz"], self.protocol["spectrum"]["proton_frequency_MHz"])
        self.assertEqual(result["linewidth_FWHM_Hz"], self.protocol["spectrum"]["proton_linewidth_Hz"])
        self.assertEqual(result["threshold_Hz"], self.protocol["spectrum"]["J_graph_threshold_Hz"])
        self.assertEqual(result["source_atom_rows"], atoms)
        self.assertEqual(result["source_j_pairs"], pairs)
        self.assertIn("no noise", result["broadening"])
        self.assertTrue(np.isfinite(result["y"]).all())
        self.assertGreaterEqual(min(result["y"]), 0)

    def test_noise_free_output_is_deterministic(self):
        self.assertIs(self.protocol["spectrum"]["noise"], False)
        atoms, pairs = [atom(0, "H", 1), atom(1, "H", 1.025)], [pair(0, 1, 8)]
        first = make_spectrum(atoms, pairs, "H", self.protocol)
        second = make_spectrum(atoms, pairs, "H", self.protocol)
        self.assertEqual(first, second)

    def test_noise_request_is_explicitly_rejected(self):
        changed = copy.deepcopy(self.protocol)
        changed["spectrum"]["noise"] = True
        with self.assertRaisesRegex(ValueError, "noise"):
            make_spectrum([atom(0, "H", 1)], [], "H", changed)

    def test_missing_J_endpoint_rejected_in_both_nuclei(self):
        for nucleus, atoms in (("H", [atom(0, "H", 1)]), ("13C", [atom(0, "C", 100)])):
            with self.subTest(nucleus=nucleus), self.assertRaisesRegex(ValueError, "endpoint"):
                make_spectrum(atoms, [pair(0, 99, 8)], nucleus, self.protocol)

    def test_equivalent_protons_do_not_self_split(self):
        result = make_spectrum([atom(0, "H", 2), atom(1, "H", 2)], [pair(0, 1, 12)], "H", self.protocol)
        self.assertEqual(len(result["sticks"]), 1)
        self.assertEqual(result["sticks"][0]["ppm"], 2)
        self.assertAlmostEqual(result["sticks"][0]["intensity"], 2)
        self.assertAlmostEqual(area(result), 2)

    def test_exchangeable_peak_is_retained_without_coherent_J(self):
        result = make_spectrum([atom(0, "H", 1), atom(1, "H", 5, True)], [pair(0, 1, 8)], "H", self.protocol)
        self.assertEqual([r["ppm"] for r in result["sticks"]], [1, 5])
        self.assertTrue(any("exchangeable H" in message for message in result["warnings"]))
        self.assertEqual(len(result["source_j_pairs"]), 1)
        self.assertAlmostEqual(area(result), 2)


if __name__ == "__main__":
    unittest.main()
