from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT / "scripts" / "build_workflow_profile.py"
SPEC = importlib.util.spec_from_file_location("build_workflow_profile", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)

PROPERTY_SCRIPT = PROJECT / "scripts" / "build_property_jobs.py"
PROPERTY_SPEC = importlib.util.spec_from_file_location(
    "build_property_jobs_profile_test", PROPERTY_SCRIPT
)
assert PROPERTY_SPEC is not None and PROPERTY_SPEC.loader is not None
PROPERTY_MODULE = importlib.util.module_from_spec(PROPERTY_SPEC)
sys.modules[PROPERTY_SPEC.name] = PROPERTY_MODULE
sys.path.insert(0, str(PROJECT / "scripts"))
PROPERTY_SPEC.loader.exec_module(PROPERTY_MODULE)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


class WorkflowProfileTests(unittest.TestCase):
    def test_default_practical_profile_is_small_and_water_is_mandatory(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "profiles"
            summary = MODULE.materialize(
                PROJECT / "config" / "workflow_profiles.json",
                PROJECT / "data" / "GOAT_explicit_job_manifest.csv",
                PROJECT / "data" / "explicit_scenario_registry.csv",
                output,
            )
            self.assertEqual(summary["workflow_profile"], "PRACTICAL_HIGH_ACCURACY")
            self.assertTrue(summary["default_profile"])
            self.assertEqual(summary["compounds"], 2)
            self.assertEqual(summary["fixed_composition_scenarios"], 7)
            self.assertEqual(summary["goat_jobs"], 70)
            self.assertEqual(summary["mandatory_trace_water_scenarios"], 2)

            manifest = rows(output / "PRACTICAL_HIGH_ACCURACY" / "GOAT_explicit_job_manifest.csv")
            self.assertEqual(len(manifest), 70)
            self.assertTrue(all(row["cross_stoichiometry_weight"] == "NONE" for row in manifest))
            water = [row for row in manifest if row["scenario_role"] == "MANDATORY_TRACE_WATER_SENSITIVITY"]
            self.assertEqual(len(water), 20)
            self.assertTrue(all(row["explicit_water_count"] == "1" for row in water))
            self.assertTrue(all(row["main_prediction_eligibility"] == "false" for row in water))
            matched_dry = [
                row
                for row in manifest
                if row["figure_id"] == "F3_015"
                and row["scenario_role"] == "MATCHED_DRY_CONTROL_FOR_TRACE_WATER"
            ]
            self.assertEqual(len(matched_dry), 10)
            self.assertTrue(all(row["scenario_id"] == "DMSO_02_W00_primary" for row in matched_dry))

    def test_ultra_profile_preserves_all_200_source_jobs(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "profiles"
            summary = MODULE.materialize(
                PROJECT / "config" / "workflow_profiles.json",
                PROJECT / "data" / "GOAT_explicit_job_manifest.csv",
                PROJECT / "data" / "explicit_scenario_registry.csv",
                output,
                "ULTRA_DIAGNOSTIC",
            )
            self.assertEqual(summary["goat_jobs"], 200)
            self.assertEqual(summary["fixed_composition_scenarios"], 10)
            self.assertFalse(summary["default_profile"])

    def test_profile_contract_declares_separate_clean_and_synthetic_deliverables(self):
        payload = json.loads(
            (PROJECT / "config" / "workflow_profiles.json").read_text(encoding="utf-8")
        )
        practical = payload["profiles"]["PRACTICAL_HIGH_ACCURACY"]
        self.assertEqual(
            practical["fixed_composition_policy"]["cross_stoichiometry_boltzmann_mixing"],
            "FORBIDDEN",
        )
        self.assertEqual(practical["deliverables"]["clean_calculated_spectrum"], "required")
        self.assertEqual(
            practical["deliverables"]["synthetic_acquisition_appearance"],
            "required_and_separately_labelled",
        )
        self.assertEqual(
            practical["fixed_composition_policy"]["matched_water_delta_pairs"]["F3_015"],
            {
                "dry": "DMSO_02_W00_primary",
                "wet": "DMSO_02_W01_trace_water_sensitivity",
            },
        )

    def test_adaptive_property_scope_uses_union_across_nucleus_temperatures(self):
        weights = {
            ("WET", "1H"): {
                "weights": [
                    {"conformer_id": "A", "weight": 0.70},
                    {"conformer_id": "B", "weight": 0.25},
                    {"conformer_id": "C", "weight": 0.05},
                ]
            },
            ("WET", "13C"): {
                "weights": [
                    {"conformer_id": "B", "weight": 0.60},
                    {"conformer_id": "C", "weight": 0.35},
                    {"conformer_id": "A", "weight": 0.05},
                ]
            },
        }
        selected = PROPERTY_MODULE.cumulative_population_union(weights, "WET", 0.95)
        self.assertEqual(selected, {"A", "B", "C"})

    def test_property_stage_requires_exact_materialized_profile_provenance(self):
        (
            profile_name,
            _scope,
            config_hash,
            manifest_hash,
            scenario_registry_hash,
            allowed_pairs,
            allowed_source_job_pairs,
        ) = PROPERTY_MODULE.load_workflow_property_scope(
            PROJECT / "config" / "workflow_profiles.json",
            "PRACTICAL_HIGH_ACCURACY",
            PROJECT,
        )
        source_rows = rows(
            PROJECT
            / "data"
            / "workflow_profiles"
            / profile_name
            / "GOAT_explicit_job_manifest.csv"
        )
        first_by_pair = {}
        for row in source_rows:
            first_by_pair.setdefault((row["figure_id"], row["scenario_id"]), row)
        dft_rows = [
            {
                "job_id": f"{figure_id}__{scenario_id}__c001",
                "figure_id": figure_id,
                "scenario_id": scenario_id,
                "source_goat_job_id": first_by_pair[(figure_id, scenario_id)]["job_id"],
                "workflow_profile": profile_name,
                "workflow_profile_config_sha256": config_hash,
                "source_goat_manifest_sha256": manifest_hash,
                "source_scenario_registry_sha256": scenario_registry_hash,
            }
            for figure_id, scenario_id in sorted(allowed_pairs)
        ]
        PROPERTY_MODULE.validate_workflow_profile_binding(
            dft_rows,
            profile_name,
            config_hash,
            manifest_hash,
            scenario_registry_hash,
            allowed_pairs,
            allowed_source_job_pairs,
        )
        tampered = [dict(row) for row in dft_rows]
        tampered[0]["source_goat_manifest_sha256"] = "0" * 64
        with self.assertRaises(PROPERTY_MODULE.WorkflowValidationError):
            PROPERTY_MODULE.validate_workflow_profile_binding(
                tampered,
                profile_name,
                config_hash,
                manifest_hash,
                scenario_registry_hash,
                allowed_pairs,
                allowed_source_job_pairs,
            )
        cross_scenario = [dict(row) for row in dft_rows]
        foreign_source = next(
            job_id
            for job_id, pair in allowed_source_job_pairs.items()
            if pair != (cross_scenario[0]["figure_id"], cross_scenario[0]["scenario_id"])
        )
        cross_scenario[0]["source_goat_job_id"] = foreign_source
        with self.assertRaises(PROPERTY_MODULE.WorkflowValidationError):
            PROPERTY_MODULE.validate_workflow_profile_binding(
                cross_scenario,
                profile_name,
                config_hash,
                manifest_hash,
                scenario_registry_hash,
                allowed_pairs,
                allowed_source_job_pairs,
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
