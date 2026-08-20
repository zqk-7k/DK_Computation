from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from adaptive_goat_controller import evaluate_gate, plan_group
from classify_failure import classify_text
from render_slurm import render
from run_orca_job import validate_profile
from sync_scheduler_states import scheduler_state


def job(tmp: Path, name: str, wave: int, state: str, accepted: bool = False):
    return {"compound_id": "F3_001", "scenario_id": "S1", "job_id": name, "wave": wave,
            "target_success": 2, "job_dir": tmp / name, "state": state, "accepted": accepted, "submit_script": "x.sbatch"}


class FailureClassifierTests(unittest.TestCase):
    def test_memory_is_not_automatic_retry(self):
        result = classify_text("ORCA MP2 RESPONSE: cannot allocate memory")
        self.assertEqual(result["primary_category"], "MEMORY_EXHAUSTION")
        self.assertFalse(result["automatic_retry_permitted"])

    def test_filesystem_precedes_generic_orca_error(self):
        result = classify_text("Input/output error: Read-only file system\nORCA finished by error termination")
        self.assertEqual(result["primary_category"], "FILESYSTEM_READ_ONLY")


class AdaptiveControllerTests(unittest.TestCase):
    def test_initial_wave_releases_held_only(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td); jobs = [job(tmp, "j1", 1, "HELD"), job(tmp, "j2", 1, "HELD"), job(tmp, "j3", 2, "HELD")]
            result = plan_group(jobs, tmp / "gates")
            self.assertEqual(result["decision"], "RELEASE_WAVE")
            self.assertEqual(set(result["release_job_ids"]), {"j1", "j2"})

    def test_running_job_forces_wait(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td); jobs = [job(tmp, "j1", 1, "RUNNING"), job(tmp, "j2", 2, "HELD")]
            result = plan_group(jobs, tmp / "gates")
            self.assertEqual(result["decision"], "WAIT_RELEASED_TASKS_TERMINAL")
            self.assertEqual(result["release_job_ids"], [])

    def test_complete_gate_stops_without_cancellation(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td); gate_root = tmp / "gates"; path = gate_root / "F3_001" / "S1.json"; path.parent.mkdir(parents=True)
            jobs = [job(tmp, "j1", 1, "COMPLETED", True), job(tmp, "j2", 1, "COMPLETED", True), job(tmp, "j3", 2, "HELD")]
            gate = {"evaluated_job_ids": ["j1", "j2"], "accepted_success_count": 2,
                    "identity_gate_pass": True, "composition_gate_pass": True, "stereochemistry_gate_pass": True,
                    "low_energy_basin_rediscovered": True, "motif_coverage_pass": True,
                    "shift_stability": {"evaluated": True, "h_nonexchangeable_max_delta_ppm": 0.03,
                                        "c_max_delta_ppm": 0.3, "nh_max_delta_ppm": 0.1, "peak_order_stable": True}}
            path.write_text(json.dumps(gate), encoding="utf-8")
            result = plan_group(jobs, gate_root)
            self.assertEqual(result["decision"], "STOP_TARGET_AND_COVERAGE_PASSED")
            self.assertEqual(result["release_job_ids"], [])

    def test_stale_assessment_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            jobs = [job(tmp, "j1", 1, "COMPLETED", True), job(tmp, "j2", 1, "FAILED"), job(tmp, "j3", 2, "HELD")]
            gate = {"evaluated_job_ids": ["j1"], "accepted_success_count": 1,
                    "identity_gate_pass": True, "composition_gate_pass": True, "stereochemistry_gate_pass": True,
                    "low_energy_basin_rediscovered": False, "motif_coverage_pass": False, "shift_stability": {"evaluated": False}}
            path = tmp / "gates" / "F3_001" / "S1.json"; path.parent.mkdir(parents=True); path.write_text(json.dumps(gate), encoding="utf-8")
            result = plan_group(jobs, tmp / "gates")
            self.assertEqual(result["decision"], "WAIT_GATE_ASSESSMENT")
            self.assertIn("assessment does not cover every released job", result["reasons"])


class ResourceAndRenderingTests(unittest.TestCase):
    def test_all_profiles_obey_memory_contract(self):
        data = json.loads((ROOT / "hpc" / "resource_profiles.json").read_text(encoding="utf-8"))
        limit = data["policy"]["memory_fraction_limit"]
        for name, profile in data["profiles"].items():
            with self.subTest(name=name):
                self.assertLessEqual(profile["orca_nprocs"] * profile["maxcore_mb"], limit * profile["requested_memory_mb"])

    def test_orca_input_matches_profile(self):
        with tempfile.TemporaryDirectory() as td:
            inp = Path(td) / "x.inp"; inp.write_text("! PBE0\n%pal nprocs 8 end\n%maxcore 7000\n* xyz 0 1\nH 0 0 0\n*\n", encoding="utf-8")
            validate_profile(inp, {"orca_nprocs": 8, "maxcore_mb": 7000, "requested_memory_mb": 128000}, 0.70)

    def test_unresolved_slurm_token_fails(self):
        with self.assertRaises(ValueError): render("{{A}} {{B}}", {"A": "x"})

    def test_scheduler_state_normalisation(self):
        self.assertEqual(scheduler_state(" OUT_OF_MEMORY+\n"), "FAILED")
        self.assertEqual(scheduler_state(" COMPLETED\n"), "COMPLETED")


if __name__ == "__main__":
    unittest.main()
