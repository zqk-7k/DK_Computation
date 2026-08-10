from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT / "scripts" / "render_experimental_appearance.py"
DEFAULT_CONFIG = PROJECT / "config" / "appearance_profiles.json"

SPEC = importlib.util.spec_from_file_location("render_experimental_appearance", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_lines(path: Path, nucleus: str = "1H") -> None:
    rows = [
        {
            "figure_id": "TEST_001",
            "nucleus": nucleus,
            "ppm": "7.25" if nucleus == "1H" else "128.4",
            "intensity": "1.0",
            "linewidth_hz": "",
            "exchangeable_H_flag": "false",
            "multiplicity": "d" if nucleus == "1H" else "s",
            "J_hz": "7.8" if nucleus == "1H" else "",
        },
        {
            "figure_id": "TEST_001",
            "nucleus": nucleus,
            "ppm": "8.10" if nucleus == "1H" else "167.8",
            "intensity": "1.0",
            "linewidth_hz": "",
            "exchangeable_H_flag": "true" if nucleus == "1H" else "false",
            "multiplicity": "s",
            "J_hz": "",
        },
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_column(path: Path, name: str) -> np.ndarray:
    with path.open(encoding="utf-8", newline="") as handle:
        return np.asarray([float(row[name]) for row in csv.DictReader(handle)])


def small_config(destination: Path) -> Path:
    payload = json.loads(DEFAULT_CONFIG.read_text(encoding="utf-8"))
    for profile in payload["profiles"].values():
        profile["axis"]["points"] = 4096
        profile["plot_dpi"] = 72
    destination.write_text(json.dumps(payload), encoding="utf-8")
    return destination


class AppearanceLayerTests(unittest.TestCase):
    def test_seed_reproducibility_and_two_distinct_layers(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "lines.csv"
            config = small_config(root / "profiles.json")
            write_lines(source)
            first = MODULE.render(source, config, root / "first", seed_override=12345)
            second = MODULE.render(source, config, root / "second", seed_override=12345)
            third = MODULE.render(source, config, root / "third", seed_override=54321)

            first_d = load_column(first["deterministic_csv"], "deterministic_total")
            second_d = load_column(second["deterministic_csv"], "deterministic_total")
            third_d = load_column(third["deterministic_csv"], "deterministic_total")
            np.testing.assert_array_equal(first_d, second_d)
            np.testing.assert_array_equal(first_d, third_d)

            first_s = load_column(first["stochastic_csv"], "stochastic_total")
            second_s = load_column(second["stochastic_csv"], "stochastic_total")
            third_s = load_column(third["stochastic_csv"], "stochastic_total")
            np.testing.assert_array_equal(first_s, second_s)
            self.assertFalse(np.array_equal(first_s, third_s))
            self.assertFalse(np.array_equal(first_d, first_s))
            for key in ("deterministic_png", "stochastic_png"):
                self.assertGreater(first[key].stat().st_size, 1000)

    def test_dmso_water_and_provisional_scope_are_explicit(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "lines.csv"
            config = small_config(root / "profiles.json")
            write_lines(source)
            outputs = MODULE.render(source, config, root / "rendered")
            metadata = json.loads(outputs["metadata"].read_text(encoding="utf-8"))
            self.assertEqual(metadata["profile_status"], "PROVISIONAL_NOT_FITTED")
            self.assertFalse(metadata["scientific_scope"]["calls_orca"])
            self.assertFalse(metadata["scientific_scope"]["first_principles_prediction"])
            self.assertFalse(metadata["scientific_scope"]["experimental_fit_performed"])
            identities = {item["identity"] for item in metadata["named_components"]}
            self.assertIn("residual protonated DMSO in DMSO-d6", identities)
            self.assertIn("water in DMSO-d6", identities)

            ppm = load_column(outputs["deterministic_csv"], "ppm")
            named = load_column(outputs["deterministic_csv"], "named_solvent_water_signal")
            dmso_window = named[np.abs(ppm - 2.50) < 0.02]
            water_window = named[np.abs(ppm - 3.33) < 0.04]
            self.assertGreater(float(dmso_window.max()), 0.0)
            self.assertGreater(float(water_window.max()), 0.0)

    def test_exchangeable_h_is_broader_and_lower_for_equal_area(self):
        profile = json.loads(DEFAULT_CONFIG.read_text(encoding="utf-8"))["profiles"][
            "dmso_d6_1H_400MHz_provisional"
        ]
        profile["axis"]["points"] = 20001
        profile["named_components"] = []
        profile["baseline"]["polynomial_coefficients"] = [0.0]
        profile["baseline"]["roll"]["amplitude"] = 0.0
        profile["noise"]["sigma_relative_to_signal_max"] = 0.0
        lines = [
            MODULE.CleanLine("X", "1H", 6.0, 1.0, None, False, "s", ()),
            MODULE.CleanLine("X", "1H", 9.0, 1.0, None, True, "s", ()),
        ]
        result = MODULE.synthesize(lines, profile, seed=1)
        ppm = result["ppm"]
        narrow_height = float(result["analyte"][np.argmin(np.abs(ppm - 6.0))])
        exchange_height = float(result["analyte"][np.argmin(np.abs(ppm - 9.0))])
        self.assertGreater(narrow_height, 3.0 * exchange_height)

    def test_all_three_lineshapes_and_13c_dmso_component(self):
        axis = np.linspace(-1.0, 1.0, 10001)
        for shape in ("lorentzian", "gaussian", "voigt"):
            values = MODULE.line_profile(axis, 0.0, 0.02, shape, 0.6)
            self.assertTrue(np.all(np.isfinite(values)))
            self.assertGreater(float(values.max()), 0.0)
            integrate = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
            self.assertAlmostEqual(float(integrate(values, axis)), 1.0, delta=0.02)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "c_lines.csv"
            config = small_config(root / "profiles.json")
            write_lines(source, "13C")
            outputs = MODULE.render(source, config, root / "carbon")
            ppm = load_column(outputs["deterministic_csv"], "ppm")
            named = load_column(outputs["deterministic_csv"], "named_solvent_water_signal")
            self.assertGreater(float(named[np.abs(ppm - 39.52) < 1.0].max()), 0.0)

    def test_mixed_spectra_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "mixed.csv"
            source.write_text(
                "figure_id,nucleus,ppm,intensity\nA,1H,7.0,1\nB,1H,8.0,1\n",
                encoding="utf-8",
            )
            with self.assertRaises(MODULE.AppearanceError):
                MODULE.read_clean_lines(source)


if __name__ == "__main__":
    unittest.main(verbosity=2)
