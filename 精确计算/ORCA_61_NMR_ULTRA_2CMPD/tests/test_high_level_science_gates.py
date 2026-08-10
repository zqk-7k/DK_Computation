from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "scripts"))

from build_high_level_benchmarks import (  # noqa: E402
    BRANCHES,
    VPT2_METHOD_ID,
    activation_state,
    manifest_activation_fields,
    render_nmr,
    render_vpt2_ff,
    render_vpt2_nmr,
    validated_vpt2_reference_record,
)
from workflow_schema import WorkflowValidationError  # noqa: E402
from build_j_benchmark_jobs import j_benchmark_execution_state  # noqa: E402


class HighLevelScienceGateTests(unittest.TestCase):
    def test_current_high_level_gate_is_fail_closed(self) -> None:
        gate = PROJECT / "jobs" / "high_level_manual_opt_in_gate.json"
        for branch in [str(item["id"]) for item in BRANCHES] + [VPT2_METHOD_ID]:
            enabled, status = activation_state(PROJECT, gate, branch)
            self.assertFalse(enabled)
            self.assertEqual(status, "DRAFT_NOT_APPROVED")

    def test_autoci_is_named_diagnostic_and_collects_response_outputs(self) -> None:
        autoci = next(item for item in BRANCHES if str(item["id"]).startswith("diagnostic_autoci"))
        self.assertNotIn("gold", str(autoci["id"]).lower())
        self.assertIn("orca_autoci", autoci["required_helpers"])
        self.assertIn("benchmark.autociresp_*", autoci["required_output_globs"])

    def test_high_level_response_tolerances_are_explicit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            xyz = Path(temporary) / "optimized.xyz"
            xyz.write_text("1\nfixture\nH 0 0 0\n", encoding="utf-8")
            text = render_nmr(
                "AUTOCI-CCSD(T) NMR pcSseg-4",
                "%autoci\n  STol 1e-10\n  Density Relaxed\nend",
                0,
                1,
                {"H": [1]},
                2,
                1000,
                xyz,
                {},
            )
            self.assertIn("%method Z_Tol 1e-12 end", text)
            self.assertIn("Tol = 1e-12", text)

    def test_machine_approved_gate_maps_to_exact_hpc_stage_tokens(self) -> None:
        high = manifest_activation_fields(
            "high_level_benchmark", True, "APPROVED_FOR_DIAGNOSTIC_EXECUTION"
        )
        vpt2 = manifest_activation_fields(
            "vpt2_nmr", True, "APPROVED_FOR_DIAGNOSTIC_EXECUTION"
        )
        self.assertEqual(high["readiness"], "ready")
        self.assertEqual(high["approval_status"], "APPROVED_FOR_HIGH_LEVEL_BENCHMARK")
        self.assertEqual(vpt2["readiness"], "ready")
        self.assertEqual(vpt2["approval_status"], "APPROVED_FOR_VPT2")

    def test_vpt2_inputs_preserve_explicit_deuterium_mass_in_both_steps(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            xyz = Path(temporary) / "optimized.xyz"
            xyz.write_text(
                "3\nfixture\nO 0 0 0\nH 0 0 1\nH 0 1 0\n",
                encoding="utf-8",
            )
            mass = {2: 2.01410177812}
            ff = render_vpt2_ff(0, 1, "DMSO", 2, 1000, xyz, mass)
            nmr = render_vpt2_nmr(0, 1, "DMSO", {"H": [2, 3]}, 2, 1000, xyz, mass)
            for text in (ff, nmr):
                self.assertEqual(text.count("M=2.01410177812"), 1)
            self.assertIn('InHessName "vpt2_ff.hess"', nmr)
            self.assertIn("Requires both vpt2_ff.hess and vpt2_ff.vpt2", nmr)
            self.assertIn("PropDisp 0.05", nmr)

    def test_fake_nonempty_vpt2_reference_hashes_cannot_unlock_composite(self) -> None:
        with self.assertRaises(WorkflowValidationError):
            validated_vpt2_reference_record(
                PROJECT,
                {
                    "approval_status": "APPROVED_FOR_DIAGNOSTIC_EXECUTION",
                    "reference_id": "fake",
                    "method_fingerprint": "wrong",
                    "reference_geometry_sha256": "x",
                },
            )

    def test_j_benchmark_input_generation_does_not_authorize_execution(self) -> None:
        enabled, status = j_benchmark_execution_state(
            PROJECT,
            PROJECT / "jobs" / "j_benchmark_execution_gate.json",
            ["A" * 64, "B" * 64],
            "C" * 64,
        )
        self.assertFalse(enabled)
        self.assertEqual(status, "DRAFT_NOT_APPROVED")


if __name__ == "__main__":
    unittest.main()
