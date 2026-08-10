from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT / "scripts" / "compare_trace_water_sensitivity.py"
SPEC = importlib.util.spec_from_file_location("compare_trace_water_sensitivity_test", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_bundle(
    path: Path,
    scenario_id: str,
    shifts: dict[int, float],
    method_fingerprint: str = "revDSD-test-fingerprint",
    dft_keywords: str = "PBE0 D4 def2-TZVPP VeryTightOpt Freq",
) -> None:
    qrrho_protocol = {
        "electronic_structure_keywords": dft_keywords,
        "thermochemistry_engine": "ORCA_QRRHO_FROM_THE_SAME_ISOTOPE_AWARE_OPT_FREQ_HESSIAN",
        "temperature_policy": "NUCLEUS_SPECIFIC_LOCKED_ACQUISITION_TEMPERATURES",
    }
    qrrho_fingerprint = hashlib.sha256(
        json.dumps(qrrho_protocol, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest().upper()
    nuclei = [
        {
            "element": "H",
            "nucleus_index_0based": 7,
            "atom_label_1based": "H8",
            "exchangeable_H_flag": 1,
            "raw_tms_reference_shift_ppm": shifts[7],
        },
        {
            "element": "H",
            "nucleus_index_0based": 8,
            "atom_label_1based": "H9",
            "exchangeable_H_flag": 1,
            "raw_tms_reference_shift_ppm": shifts[8],
        },
        {
            "element": "H",
            "nucleus_index_0based": 9,
            "atom_label_1based": "H10",
            "exchangeable_H_flag": 0,
            "raw_tms_reference_shift_ppm": shifts[9],
        },
        {
            "element": "C",
            "nucleus_index_0based": 0,
            "atom_label_1based": "C1",
            "exchangeable_H_flag": 0,
            "raw_tms_reference_shift_ppm": shifts[0],
        },
    ]
    path.write_text(
        json.dumps(
            {
                "compound": {
                    "figure_id": "F3_015",
                    "scenario_id": scenario_id,
                    "inchikey": "VLJNHYLEOZPXFW-BYPYZUCNSA-N",
                },
                "methods": {
                    "shielding_method": "production_revdsd_pcsseg3",
                    "shielding_method_fingerprint": method_fingerprint,
                    "dft_opt_freq_input_keywords": dft_keywords,
                    "opt_freq_qrrho_protocol": qrrho_protocol,
                    "opt_freq_qrrho_protocol_fingerprint": qrrho_fingerprint,
                    "shielding_input_keywords": "revDSD-PBEP86-D4/2021 pcSseg-3 NMR",
                    "tms_reference_output": "reference_jobs/TMS/revdsd/tms_shield.out",
                },
                "weighting": {
                    "source": "data/NMR_conformer_weights.csv",
                    "policy": "authoritative scenario/nucleus/temperature-specific weights; no processor-side reweighting",
                    "nucleus_weight_sets": {
                        "1H": {
                            "temperature_K": 298.15,
                            "condition_status": "LOCKED_TEST",
                            "condition_source_note": "test",
                        },
                        "13C": {
                            "temperature_K": 298.15,
                            "condition_status": "LOCKED_TEST",
                            "condition_source_note": "test",
                        },
                    },
                },
                "condition_state": {
                    "weight_gate_binding": {
                        "approved_minimum_gate_sha256": "F" * 64,
                    }
                },
                "source_manifest_hashes": {
                    key: {"path": f"data/{key}.csv", "sha256": character * 64}
                    for key, character in {
                        "dft_manifest": "1",
                        "weight_table": "2",
                        "nucleus_conditions": "3",
                        "shielding_manifest": "4",
                        "coupling_manifest": "5",
                        "reference_manifest": "6",
                        "reference_registry": "7",
                    }.items()
                },
                "reference_shieldings": {
                    "H": {"reference_output_sha256": "A" * 64},
                    "C": {"reference_output_sha256": "B" * 64},
                },
                "reference_registry_approval": {
                    "registry_sha256": "C" * 64,
                    "status": "APPROVED",
                },
                "empirical_regression_policy": {"status": "disabled"},
                "averaged_nuclei": nuclei,
            }
        ),
        encoding="utf-8",
    )


class TraceWaterSensitivityTests(unittest.TestCase):
    def test_matched_comparison_flags_deltas_and_exchangeable_order_change(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dry = root / "dry.json"
            wet = root / "wet.json"
            write_bundle(dry, "DMSO_02_W00_primary", {7: 7.8, 8: 7.4, 9: 7.5, 0: 128.0})
            write_bundle(
                wet,
                "DMSO_02_W01_trace_water_sensitivity",
                {7: 7.2, 8: 8.1, 9: 7.63, 0: 129.2},
            )
            result = MODULE.compare(
                dry,
                wet,
                PROJECT / "config" / "workflow_profiles.json",
                root / "out",
            )
            report = result["report"]
            self.assertEqual(report["status"], "ENVIRONMENT_SENSITIVE")
            self.assertTrue(report["exchangeable_H_order_change"])
            self.assertEqual(report["n_environment_sensitive"], 4)
            self.assertTrue(result["csv"].is_file() and result["json"].is_file())
            self.assertEqual(report["cross_stoichiometry_boltzmann_mixing"], "FORBIDDEN")

    def test_unmatched_dry_scenario_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dry = root / "dry.json"
            wet = root / "wet.json"
            shifts = {7: 7.8, 8: 7.4, 9: 7.5, 0: 128.0}
            write_bundle(dry, "DMSO_03_W00_extended", shifts)
            write_bundle(wet, "DMSO_02_W01_trace_water_sensitivity", shifts)
            with self.assertRaises(MODULE.WaterSensitivityError):
                MODULE.compare(
                    dry,
                    wet,
                    PROJECT / "config" / "workflow_profiles.json",
                    root / "out",
                )

    def test_mixed_method_is_rejected_as_not_a_water_effect(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dry = root / "dry.json"
            wet = root / "wet.json"
            shifts = {7: 7.8, 8: 7.4, 9: 7.5, 0: 128.0}
            write_bundle(dry, "DMSO_02_W00_primary", shifts)
            write_bundle(
                wet,
                "DMSO_02_W01_trace_water_sensitivity",
                shifts,
                method_fingerprint="different-method",
            )
            with self.assertRaises(MODULE.WaterSensitivityError):
                MODULE.compare(
                    dry,
                    wet,
                    PROJECT / "config" / "workflow_profiles.json",
                    root / "out",
                )

    def test_missing_opt_freq_method_provenance_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dry = root / "dry.json"
            wet = root / "wet.json"
            shifts = {7: 7.8, 8: 7.4, 9: 7.5, 0: 128.0}
            write_bundle(dry, "DMSO_02_W00_primary", shifts)
            write_bundle(wet, "DMSO_02_W01_trace_water_sensitivity", shifts)
            payload = json.loads(wet.read_text(encoding="utf-8"))
            del payload["methods"]["dft_opt_freq_input_keywords"]
            wet.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaises(MODULE.WaterSensitivityError):
                MODULE.compare(
                    dry,
                    wet,
                    PROJECT / "config" / "workflow_profiles.json",
                    root / "out",
                )

    def test_different_opt_freq_qrrho_protocol_is_not_called_a_water_effect(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dry = root / "dry.json"
            wet = root / "wet.json"
            shifts = {7: 7.8, 8: 7.4, 9: 7.5, 0: 128.0}
            write_bundle(dry, "DMSO_02_W00_primary", shifts)
            write_bundle(
                wet,
                "DMSO_02_W01_trace_water_sensitivity",
                shifts,
                dft_keywords="B3LYP D3 def2-SVP Opt Freq",
            )
            with self.assertRaises(MODULE.WaterSensitivityError):
                MODULE.compare(
                    dry,
                    wet,
                    PROJECT / "config" / "workflow_profiles.json",
                    root / "out",
                )

    def test_different_approved_weight_gate_is_not_called_a_water_effect(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dry = root / "dry.json"
            wet = root / "wet.json"
            shifts = {7: 7.8, 8: 7.4, 9: 7.5, 0: 128.0}
            write_bundle(dry, "DMSO_02_W00_primary", shifts)
            write_bundle(wet, "DMSO_02_W01_trace_water_sensitivity", shifts)
            payload = json.loads(wet.read_text(encoding="utf-8"))
            payload["condition_state"]["weight_gate_binding"][
                "approved_minimum_gate_sha256"
            ] = "E" * 64
            wet.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaises(MODULE.WaterSensitivityError):
                MODULE.compare(
                    dry,
                    wet,
                    PROJECT / "config" / "workflow_profiles.json",
                    root / "out",
                )

    def test_different_weight_table_hash_is_not_called_a_water_effect(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dry = root / "dry.json"
            wet = root / "wet.json"
            shifts = {7: 7.8, 8: 7.4, 9: 7.5, 0: 128.0}
            write_bundle(dry, "DMSO_02_W00_primary", shifts)
            write_bundle(wet, "DMSO_02_W01_trace_water_sensitivity", shifts)
            payload = json.loads(wet.read_text(encoding="utf-8"))
            payload["source_manifest_hashes"]["weight_table"]["sha256"] = "8" * 64
            wet.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaises(MODULE.WaterSensitivityError):
                MODULE.compare(
                    dry,
                    wet,
                    PROJECT / "config" / "workflow_profiles.json",
                    root / "out",
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
