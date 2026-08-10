from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
EVIDENCE = PROJECT / "validation/local_orca_6_1_1_syntax_smoke.json"


class LocalSyntaxSmokeEvidenceTests(unittest.TestCase):
    def test_evidence_is_hash_bound_and_not_scientific_validation(self) -> None:
        payload = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        fixtures = PROJECT / payload["fixture_directory"]
        for name, expected in payload["fixture_sha256"].items():
            path = fixtures / name
            self.assertTrue(path.is_file(), path)
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest().upper(), expected)

        release = payload["release_gate"]
        self.assertFalse(release["production_target_calculations_run"])
        self.assertFalse(release["scientifically_validated_against_experiment"])
        self.assertFalse(release["autoci_branch_enabled"])
        self.assertFalse(release["vpt2_branch_enabled"])

    def test_blocked_autoci_is_not_reported_as_passed(self) -> None:
        payload = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        records = {item["id"]: item for item in payload["runs"]}
        autoci = records["autoci_ccsdt_giao"]
        self.assertEqual(
            autoci["status"],
            "EXECUTION_BLOCKED_LOCAL_DISTRIBUTION_HELPER_MISSING",
        )
        self.assertFalse(autoci["normal_termination"])
        self.assertTrue(records["vpt2_two_step_nmr_property_derivative"]["property_normal_termination"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
