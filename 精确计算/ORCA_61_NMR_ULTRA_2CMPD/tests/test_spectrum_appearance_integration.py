from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT / "scripts" / "generate_nmr_spectra.py"
SPEC = importlib.util.spec_from_file_location("generate_nmr_spectra_integration", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_bundle(path: Path, figure_id: str, scenario_id: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "notice": MODULE.NOTICE,
                "compound": {"figure_id": figure_id, "scenario_id": scenario_id},
                "averaged_nuclei": [],
                "averaged_JHH": [],
                "weighting": {},
            }
        ),
        encoding="utf-8",
    )


class SpectrumAppearanceIntegrationTests(unittest.TestCase):
    def test_recursive_discovery_retains_dry_and_wet_scenarios(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_bundle(
                root
                / "F3_005"
                / "DMSO_02_W00_primary"
                / "F3_005__DMSO_02_W00_primary_prediction_bundle.json",
                "F3_005",
                "DMSO_02_W00_primary",
            )
            write_bundle(
                root
                / "F3_005"
                / "DMSO_02_W01_trace_water_sensitivity"
                / "F3_005__DMSO_02_W01_trace_water_sensitivity_prediction_bundle.json",
                "F3_005",
                "DMSO_02_W01_trace_water_sensitivity",
            )
            found = MODULE.discover_bundles(root, {"F3_005"})
            self.assertEqual(len(found), 2)

    def test_appearance_export_preserves_resolved_ppm_and_area_without_linewidth(self):
        sticks = [
            {
                "shift_ppm": 7.443,
                "intensity": 0.25,
                "component_id": "HCOMP_01",
                "component_atom_labels_1based": "H3;H4",
                "model": "exact_second_order_spin_hamiltonian",
                "is_exchangeable_H": 0,
            },
            {
                "shift_ppm": 7.932,
                "intensity": 1.0,
                "component_id": "HEXCH_01",
                "assignment_atom_label_1based": "H8",
                "model": "provisional_exchangeable_broad_singlet",
                "is_exchangeable_H": 1,
            },
        ]
        rows = MODULE.appearance_line_rows(
            sticks,
            "F3_005",
            "DMSO_02_W01_trace_water_sensitivity",
            "1H",
            "raw_TMS",
        )
        self.assertEqual([row["ppm"] for row in rows], [7.443, 7.932])
        self.assertAlmostEqual(sum(row["intensity"] for row in rows), 1.25)
        self.assertTrue(all(row["representation"] == "resolved_transitions" for row in rows))
        self.assertTrue(all(row["multiplicity"] == "s" and row["J_hz"] == "" for row in rows))
        self.assertTrue(all(row["linewidth_hz"] == "" for row in rows))
        self.assertEqual(rows[1]["exchangeable_H_flag"], 1)

    def test_scenario_specific_output_paths_do_not_collide(self):
        class Args:
            output_root = Path("X")

        dry = Args.output_root / MODULE.safe_token("F3_005") / MODULE.safe_token("DMSO_02_W00_primary")
        wet = (
            Args.output_root
            / MODULE.safe_token("F3_005")
            / MODULE.safe_token("DMSO_02_W01_trace_water_sensitivity")
        )
        self.assertNotEqual(dry, wet)

    def test_missing_provisional_regression_does_not_block_raw_tms_spectrum(self):
        bundle = {
            "averaged_nuclei": [
                {
                    "element": "H",
                    "raw_tms_reference_shift_ppm": 7.5,
                },
                {
                    "element": "C",
                    "raw_tms_reference_shift_ppm": 128.0,
                },
            ]
        }
        self.assertEqual(MODULE.available_variant_keys(bundle), ["raw_TMS"])

    def test_partial_variant_is_never_silently_used(self):
        bundle = {
            "averaged_nuclei": [
                {
                    "element": "H",
                    "raw_tms_reference_shift_ppm": 7.5,
                    "provisional_external_regression_shift_ppm": 7.4,
                },
                {
                    "element": "C",
                    "raw_tms_reference_shift_ppm": 128.0,
                },
            ]
        }
        self.assertEqual(MODULE.available_variant_keys(bundle), ["raw_TMS"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
