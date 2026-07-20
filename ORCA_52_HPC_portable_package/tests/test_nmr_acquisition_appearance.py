from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT / "scripts" / "nmr_acquisition_appearance.py"
FIXTURE = PROJECT / "tests" / "fixtures" / "nmr_standard_transitions_fixture.csv"
CALIBRATION = PROJECT / "data" / "real_nmr_calibration_audit"

SPEC = importlib.util.spec_from_file_location("nmr_acquisition_appearance", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def arguments(output_root: Path, nucleus: str, solvent: str, spectrum_id: str):
    return MODULE.parser().parse_args(
        [
            "--input",
            str(FIXTURE),
            "--spectrum-id",
            spectrum_id,
            "--display-name",
            "benzamide test fixture",
            "--nucleus",
            nucleus,
            "--solvent",
            solvent,
            "--calibration-dir",
            str(CALIBRATION),
            "--output-root",
            str(output_root),
            "--seed",
            "123456",
            "--acquired-points",
            "1024",
            "--processed-points",
            "2048",
            "--dpi",
            "120",
        ]
    )


def load_column(path: Path, field: str) -> np.ndarray:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return np.asarray([float(row[field]) for row in csv.DictReader(handle)])


class NmrAppearanceTests(unittest.TestCase):
    def test_dmso_1h_is_deterministic_named_and_visibly_labelled(self):
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            first_dir = MODULE.run_generation(arguments(Path(first), "1H", "DMSO-d6", "FIX_1H"))
            second_dir = MODULE.run_generation(arguments(Path(second), "1H", "DMSO-d6", "FIX_1H"))
            first_values = load_column(
                first_dir / "layer_B_synthetic_acquisition.csv", "final_real_arbitrary_units"
            )
            second_values = load_column(
                second_dir / "layer_B_synthetic_acquisition.csv", "final_real_arbitrary_units"
            )
            np.testing.assert_array_equal(first_values, second_values)

            metadata = json.loads((first_dir / "spectrum_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(
                metadata["acquisition_profile"]["calibration_status"],
                "LOCAL_EMPIRICAL_DMSO_DISTRIBUTION",
            )
            self.assertFalse(metadata["scientific_scope"]["anonymous_impurity_peaks_added"])
            self.assertLess(
                abs(
                    metadata["processing"]["noise_achieved_frequency_domain_lag1"]
                    - metadata["processing"]["noise_target_frequency_domain_lag1_from_calibration"]
                ),
                0.005,
            )
            self.assertTrue((first_dir / "FIX_1H_1H_NMR_layers.png").is_file())
            self.assertTrue((first_dir / "FIX_1H_1H_NMR_layers.svg").is_file())
            svg = (first_dir / "FIX_1H_1H_NMR_layers.svg").read_text(encoding="utf-8")
            self.assertIn("CALCULATED MODEL - NOT EXPERIMENTAL", svg)
            self.assertIn("SYNTHETIC APPEARANCE - NOT ACQUIRED", svg)

            with (first_dir / "component_ledger.csv").open(
                encoding="utf-8-sig", newline=""
            ) as handle:
                ledger = list(csv.DictReader(handle))
            ids = {row["component_id"] for row in ledger}
            self.assertIn("residual_dmso_d6_1h", ids)
            self.assertIn("water_in_dmso_d6_1h", ids)
            self.assertIn("complex_fid_noise", ids)
            self.assertIn("slow_baseline", ids)
            self.assertTrue(all(row["component_name"].strip() for row in ledger))
            self.assertFalse(any("anonymous" in row["component_class"].lower() for row in ledger))

    def test_d2o_status_is_explicitly_assumed(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = MODULE.run_generation(
                arguments(Path(temporary), "1H", "D2O", "FIX_D2O")
            )
            metadata = json.loads((output / "spectrum_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(
                metadata["acquisition_profile"]["calibration_status"],
                "ASSUMED_NOT_LOCALLY_CALIBRATED",
            )
            components = {row["component_id"]: row for row in metadata["named_components"]}
            self.assertEqual(
                components["hdo_in_d2o_1h"]["provenance_status"],
                "ASSUMED_NOT_LOCALLY_CALIBRATED",
            )

    def test_13c_and_direct_heteronuclear_interfaces(self):
        cases = [
            ("13C", "DMSO-d6", "FIX_13C"),
            ("19F", "DMSO-d6", "FIX_19F"),
            ("31P", "D2O", "FIX_31P"),
        ]
        with tempfile.TemporaryDirectory() as temporary:
            for nucleus, solvent, spectrum_id in cases:
                output = MODULE.run_generation(
                    arguments(Path(temporary), nucleus, solvent, spectrum_id)
                )
                self.assertTrue((output / "layer_A_calculated_model.csv").is_file())
                self.assertTrue((output / "layer_B_synthetic_acquisition.csv").is_file())
                metadata = json.loads(
                    (output / "spectrum_metadata.json").read_text(encoding="utf-8")
                )
                self.assertEqual(metadata["nucleus"], nucleus)
                if nucleus == "31P":
                    self.assertEqual(
                        metadata["acquisition_profile"]["calibration_status"],
                        "ASSUMED_NOT_LOCALLY_CALIBRATED",
                    )

    def test_hash_manifest_is_valid(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = MODULE.run_generation(
                arguments(Path(temporary), "1H", "DMSO-d6", "FIX_HASH")
            )
            with (output / "output_file_hashes.csv").open(
                encoding="utf-8-sig", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))
            self.assertGreaterEqual(len(rows), 10)
            for row in rows:
                path = output / row["file"]
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                self.assertEqual(digest, row["sha256"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
