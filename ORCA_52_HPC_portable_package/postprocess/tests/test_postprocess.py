from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


POSTPROCESS = Path(__file__).resolve().parents[1]
PROJECT_ROOT = POSTPROCESS.parents[3]
if str(POSTPROCESS) not in sys.path:
    sys.path.insert(0, str(POSTPROCESS))

from orca_io import parse_j_couplings, parse_shieldings, read_atom_map


BENZAMIDE = PROJECT_ROOT / "ORCA_benzamide_NMR_HRMS"


class PostprocessCompatibilityTests(unittest.TestCase):
    def test_benzamide_property_parsers(self):
        atom_map = read_atom_map(BENZAMIDE / "data" / "benzamide_atom_map.csv")
        shieldings = parse_shieldings(
            BENZAMIDE
            / "outputs"
            / "benzamide_shield"
            / "benzamide_shield.property.txt",
            atom_map,
        )
        couplings = parse_j_couplings(
            BENZAMIDE / "outputs" / "benzamide_j" / "benzamide_j.property.txt"
        )
        self.assertEqual(len(atom_map), 16)
        self.assertEqual(len(shieldings), 16)
        self.assertEqual(len(couplings), 16)
        self.assertAlmostEqual(shieldings[0], 144.09055868013186, places=9)

    def test_benzamide_qa_exports_are_explicitly_calculated(self):
        output = POSTPROCESS / "qa_outputs_CALCULATED_ONLY"
        manifest = json.loads(
            (output / "_aggregate_workbook_ready" / "aggregate_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(manifest["status"], "CALCULATED_NOT_EXPERIMENTAL")
        self.assertEqual(manifest["nmr_atom_shift_rows"], 14)
        self.assertEqual(manifest["hrms_ion_rows"], 1)
        self.assertIn("NO LC-MS FRAGMENTATION", manifest["hrms_scope"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
