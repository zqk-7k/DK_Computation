"""Planning tests only: no quantum executable, no sampling, no convergence claims."""
from __future__ import annotations

import importlib.util
import json
import shutil
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("plan_sampling_v4", ROOT / "scripts/plan_sampling_v4.py")
planner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(planner)


class SamplingPlanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plan, cls.solvents = planner.plan(ROOT)
        cls.species = {r["species_id"]: r for r in cls.plan["species"]}
        cls.products = {r["compound_id"]: r for r in cls.plan["products"]}

    def test_64_products_66_models_and_tms(self):
        self.assertEqual(len(self.products), 64)
        self.assertEqual(len(self.species), 67)
        self.assertEqual(self.plan["summary"]["n_product_species"], 66)
        self.assertEqual(self.plan["tms_reference"]["compound_id"], "TMS")

    def test_counts_reconcile(self):
        summary = self.plan["summary"]
        self.assertEqual(sum(p["total_starts"] for p in self.products.values()), 1012)
        self.assertEqual(sum(s["n_starts"] for s in self.species.values()), 1020)
        self.assertEqual(summary["product_starts"] + summary["TMS_starts"], summary["total_starts"])
        self.assertEqual(summary["by_figure"], {
            "2": {"n_products": 30, "n_starts": 480},
            "3": {"n_products": 28, "n_starts": 372},
            "4": {"n_products": 6, "n_starts": 160},
        })

    def test_each_start_count_from_reviewed_allowed_tier(self):
        for row in self.species.values():
            self.assertIn(row["n_starts"], (8, 12, 16, 24, 32))
            self.assertEqual(row["n_starts"], row["n_starts_per_species"])
            self.assertGreater(row["suggested_unique_max"], row["suggested_unique_min"])
            self.assertTrue(row["advisory_no_quota"])
            self.assertTrue(row["n_starts_is_not_n_final_conformers"])
            self.assertFalse(row["sampling_converged"])

    def test_same_2_plus_1_environment_with_continuum(self):
        for row in self.species.values():
            self.assertEqual(row["scenario_id"], planner.SCENARIO)
            self.assertEqual(row["main_solvent"], "DMSO-d6")
            self.assertEqual((row["explicit_DMSO_d6"], row["explicit_H2O"], row["explicit_D2O"]), (2, 1, 0))
            self.assertEqual(row["outer_sampling_model"], "ALPB(DMSO)")
            self.assertEqual(row["outer_DFT_model"], "CPCM(DMSO)")
            self.assertEqual(row["temperature_K"], 298.15)
            self.assertIn("PROVISIONAL", row["temperature_status"])
            self.assertIn("NOT_MEASURED", row["solvent_composition_status"])

    def test_branch_totals_unambiguous(self):
        for cid, count, total in (("F2_E11", 8, 16), ("F2_E22", 12, 24)):
            row = self.products[cid]
            self.assertEqual(row["n_species"], 2)
            self.assertEqual(len(row["n_starts_per_species"]), 2)
            self.assertEqual(set(row["n_starts_per_species"].values()), {count})
            self.assertEqual(row["total_starts"], total)
            self.assertIn("separate_diastereomer_branches_do_not_pool_without_experimental_ratio", row["review_flags"])
        self.assertTrue(all(p["n_species"] == 1 for cid, p in self.products.items() if cid not in {"F2_E11", "F2_E22"}))

    def test_drug_j3_and_different_scope_not_silently_forced(self):
        drugs = [p for p in self.products.values() if p["figure"] == 4]
        self.assertEqual({p["compound_id"] for p in drugs}, {"F4_D01", "F4_D02", "F4_D03", "F4_D04", "F4_D05", "F4_D07"})
        self.assertEqual({p["J_method"] for p in drugs}, {"PBE_pcJ3"})
        for p in self.products.values():
            if p["figure"] != 4:
                self.assertEqual(p["J_method"], "BENCHMARK_SELECT_PBE_pcJ3_OR_PBE0_pcJ2")

    def test_water_and_identity_not_claimed_experimentally_confirmed(self):
        for p in self.products.values():
            self.assertIn("experimental_solvent_temperature_and_water_content_not_confirmed", p["review_flags"])
            self.assertIn("computational_structure_not_experimental_identity_or_stereo_assignment", p["review_flags"])
        self.assertIn("neutral_DMSO_model_provisional_HCl_salt_and_CD3OD_branches_held", self.products["F4_D07"]["review_flags"])
        self.assertIn("sample_state_protonation_or_tautomer_review_required", self.products["F2_S28"]["review_flags"])

    def test_not_blind_rotatable_bond_rule(self):
        flexible_ring = next(s for s in self.species.values() if s["compound_id"] == "F3_B03")
        rigid_ring = next(s for s in self.species.values() if s["compound_id"] == "F3_N01")
        self.assertEqual(flexible_ring["descriptors"]["rotatable_bonds_strict"], 0)
        self.assertEqual(rigid_ring["descriptors"]["rotatable_bonds_strict"], 0)
        self.assertEqual(flexible_ring["n_starts"], 16)
        self.assertEqual(rigid_ring["n_starts"], 8)
        self.assertEqual(self.products["F2_S28"]["total_starts"], 32)

    def test_optional_controls_not_active(self):
        controls = [p for p in self.products.values() if p["optional_control_recommended"]]
        self.assertEqual(len(controls), 8)
        self.assertTrue(all(not p["optional_controls_enabled"] for p in self.products.values()))
        self.assertEqual(self.plan["summary"]["n_enabled_control_starts"], 0)

    def test_aggregated_solvent_table_is_same_64_product_plan(self):
        self.assertEqual(self.solvents["products"], self.plan["products"])
        for p in self.products.values():
            self.assertEqual(p["total_starts"], sum(p["n_starts_per_species"].values()))
            self.assertEqual(set(p["species_ids"]), set(p["n_starts_per_species"]))

    def test_generation_reproducible_and_source_bound(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "registry").mkdir()
            (root / "docs").mkdir()
            for name in ("products.json", "species.json"):
                shutil.copy2(ROOT / "registry" / name, root / "registry" / name)
            planner.write_plan(root)
            first = (root / "registry/sampling_plan_v4.json").read_bytes()
            planner.write_plan(root)
            self.assertEqual(first, (root / "registry/sampling_plan_v4.json").read_bytes())
            generated = json.loads(first)
            self.assertEqual(generated, self.plan)
            self.assertEqual(generated["source_sha256"]["registry/species.json"], planner.sha(root / "registry/species.json"))

    def test_changed_product_set_requires_review(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "registry").mkdir()
            shutil.copy2(ROOT / "registry/species.json", root / "registry/species.json")
            products = json.loads((ROOT / "registry/products.json").read_text(encoding="utf-8"))
            products.pop()
            (root / "registry/products.json").write_text(json.dumps(products), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "64-product structure set changed"):
                planner.plan(root)


if __name__ == "__main__":
    unittest.main()
