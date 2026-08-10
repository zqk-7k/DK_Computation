from __future__ import annotations

import sys
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "scripts"))

from generate_nmr_spectra import apply_declared_solution_equivalence  # noqa: E402


class BenzamideSolutionEquivalenceTests(unittest.TestCase):
    def test_declared_pairs_are_averaged_without_mutating_raw_bundle(self) -> None:
        rows = []
        for index, element, shift in (
            (11, "H", 8.0),
            (15, "H", 7.8),
            (12, "H", 7.6),
            (14, "H", 7.4),
            (4, "C", 129.0),
            (8, "C", 127.0),
            (5, "C", 128.5),
            (7, "C", 127.5),
        ):
            rows.append(
                {
                    "nucleus_index_0based": index,
                    "element": element,
                    "raw_tms_reference_shift_ppm": shift,
                    "provisional_external_regression_shift_ppm": shift + 0.1,
                }
            )
        bundle = {
            "compound": {"figure_id": "F3_005", "inchikey": "KXDAEFPNCMNJSK-UHFFFAOYSA-N"},
            "averaged_nuclei": rows,
            "averaged_JHH": [
                {"nucleus_A_index_0based": 11, "nucleus_B_index_0based": 12, "boltzmann_averaged_J_Hz": 7.0},
                {"nucleus_A_index_0based": 14, "nucleus_B_index_0based": 15, "boltzmann_averaged_J_Hz": 9.0},
                {"nucleus_A_index_0based": 11, "nucleus_B_index_0based": 14, "boltzmann_averaged_J_Hz": 1.0},
                {"nucleus_A_index_0based": 12, "nucleus_B_index_0based": 15, "boltzmann_averaged_J_Hz": 3.0},
            ],
        }
        adjusted, audit = apply_declared_solution_equivalence(
            bundle, PROJECT / "config" / "solution_equivalence_groups.json"
        )
        adjusted_nuclei = {
            int(row["nucleus_index_0based"]): row for row in adjusted["averaged_nuclei"]
        }
        self.assertEqual(adjusted_nuclei[11]["raw_tms_reference_shift_ppm"], 7.9)
        self.assertEqual(adjusted_nuclei[15]["raw_tms_reference_shift_ppm"], 7.9)
        self.assertEqual(adjusted_nuclei[4]["raw_tms_reference_shift_ppm"], 128.0)
        self.assertEqual(bundle["averaged_nuclei"][0]["raw_tms_reference_shift_ppm"], 8.0)
        values = sorted(float(row["boltzmann_averaged_J_Hz"]) for row in adjusted["averaged_JHH"])
        self.assertEqual(values, [2.0, 2.0, 8.0, 8.0])
        self.assertTrue(audit)
        self.assertTrue(all(row["raw_source_retained"] for row in audit))


if __name__ == "__main__":
    unittest.main()
