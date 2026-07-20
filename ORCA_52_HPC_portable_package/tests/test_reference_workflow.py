from __future__ import annotations

import csv
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "scripts"))

from prepare_reference_jobs import prepare_reference_jobs  # noqa: E402
from run_reference_serial import (  # noqa: E402
    _parse_nuclei_requests,
    dry_run_validate,
    run_reference_serial,
)
from workflow_schema import WorkflowValidationError, read_csv, sha256_file  # noqa: E402


class ReferenceWorkflowTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.project = Path(self.temp.name) / "project"
        self.project.mkdir()
        (self.project / "config").mkdir()
        (self.project / "data").mkdir()
        (self.project / "jobs").mkdir()
        prepare_reference_jobs(self.project)
        self.manifest = self.project / "data" / "reference_job_manifest_DRAFT.csv"
        self.gate = self.project / "jobs" / "reference_model_protocol_DRAFT.json"

    def tearDown(self) -> None:
        self.temp.cleanup()

    def _manifest_rows(self) -> list[dict[str, str]]:
        return read_csv(self.manifest)

    def test_exact_models_hashes_reference_sites_and_collision_screen(self) -> None:
        rows = self._manifest_rows()
        self.assertEqual(len(rows), 6)
        self.assertEqual(
            {row["model_id"] for row in rows},
            {
                "REF_DMSO_TMS",
                "REF_DMSO_CFCL3",
                "REF_DMSO_31P_H3PO4_2DMSO_MODEL",
                "REF_D2O_31P_H3PO4_4H2O_MODEL",
                "REF_D2O_DSS_NA_6D2O",
                "REF_D2O_TSP_NA_6D2O",
            },
        )
        for row in rows:
            directory = Path(row["job_directory"])
            self.assertEqual(sha256_file(directory / "start.xyz"), row["start_xyz_sha256"])
            self.assertEqual(
                sha256_file(directory / "reference_atom_map.csv"),
                row["reference_atom_map_sha256"],
            )
            self.assertEqual(sha256_file(directory / "opt_freq.inp"), row["opt_input_sha256"])
            self.assertEqual(
                sha256_file(directory / "shielding.inp.template"),
                row["shielding_template_sha256"],
            )
            self.assertIn("AutoAux|RIJCOSX|GIAO", row["shielding_method_fingerprint"])
            self.assertIn(
                row["initial_collision_check"],
                {"PASS", "NOT_APPLICABLE_SINGLE_COMPONENT"},
            )
            if row["initial_collision_check"] == "PASS":
                self.assertGreaterEqual(float(row["minimum_nonbonded_clearance_ratio"]), 1.0)

        for model_id in ("REF_D2O_DSS_NA_6D2O", "REF_D2O_TSP_NA_6D2O"):
            directory = self.project / "reference_jobs" / model_id
            atom_rows = [
                row for row in read_csv(directory / "reference_atom_map.csv")
                if row["reference_site"] == "true"
            ]
            self.assertEqual(sum(row["element"] == "H" for row in atom_rows), 9)
            self.assertEqual(sum(row["element"] == "C" for row in atom_rows), 3)
            template = (directory / "shielding.inp.template").read_text(encoding="utf-8")
            requests = _parse_nuclei_requests(template)
            self.assertEqual(len(requests["H"]), 9)
            self.assertEqual(len(requests["C"]), 3)
            environment = {
                int(row["atom_index_1based"])
                for row in read_csv(directory / "reference_atom_map.csv")
                if row["environment_atom"] == "true"
            }
            self.assertTrue(set(requests["H"] + requests["C"]).isdisjoint(environment))

    def test_deuterium_counts_match_declared_models(self) -> None:
        expected = {
            "REF_DMSO_31P_H3PO4_2DMSO_MODEL": 12,
            "REF_D2O_31P_H3PO4_4H2O_MODEL": 11,
            "REF_D2O_DSS_NA_6D2O": 12,
            "REF_D2O_TSP_NA_6D2O": 12,
        }
        for model_id, count in expected.items():
            directory = self.project / "reference_jobs" / model_id
            mapped = sum(
                bool(row["isotope_mass_u"])
                for row in read_csv(directory / "reference_atom_map.csv")
            )
            self.assertEqual(mapped, count)
            self.assertEqual(
                (directory / "opt_freq.inp").read_text(encoding="utf-8").count(
                    "M=2.01410177812"
                ),
                count,
            )

    def test_raw_reference_processing_and_external_standard_are_separate(self) -> None:
        gate = json.loads(self.gate.read_text(encoding="utf-8"))
        route = gate["d2o_1h_13c_reference_route_decision"]
        self.assertEqual(route["status"], "UNRESOLVED_REVIEW_REQUIRED")
        self.assertEqual(
            set(route["allowed_routes"]),
            {
                "REF_D2O_DSS_NA_6D2O",
                "REF_D2O_TSP_NA_6D2O",
                "INDIRECT_REFERENCE_NO_ORCA_JOB",
            },
        )
        registry = json.loads(
            (self.project / "config" / "reference_job_registry.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(registry["status"], "DRAFT_NOT_APPROVED")
        self.assertEqual(len(registry["indirect_reference_routes"]), 1)
        self.assertIn(
            "no molecular ORCA reference job",
            registry["indirect_reference_routes"][0]["limitation"],
        )
        d2o_candidates: dict[str, set[str]] = {}
        for entry in registry["reference_jobs"]:
            self.assertNotEqual(
                entry["raw_calculation_reference"], entry["experimental_reference_convention"]
            )
            if entry["solvent_domain"] == "D2O":
                d2o_candidates.setdefault(entry["nucleus"], set()).add(entry["model_id"])
        self.assertEqual(
            d2o_candidates["1H"], {"REF_D2O_DSS_NA_6D2O", "REF_D2O_TSP_NA_6D2O"}
        )
        for model_id in (
            "REF_DMSO_31P_H3PO4_2DMSO_MODEL",
            "REF_D2O_31P_H3PO4_4H2O_MODEL",
        ):
            declaration = json.loads(
                (
                    self.project / "reference_jobs" / model_id / "model_declaration.json"
                ).read_text(encoding="utf-8")
            )
            self.assertIs(
                declaration["equivalent_to_85wt_percent_H3PO4_external_standard"], False
            )
            self.assertIn(
                "MODEL_ONLY_NOT_EQUIVALENT_TO_85WT_PERCENT_H3PO4",
                declaration["model_limitation"],
            )
            self.assertIn(
                "PRESERVE_EACH_REQUESTED_ATOM_VALUE",
                declaration["raw_shielding_aggregation_policy"],
            )

    def test_dry_run_validates_all_and_does_not_start_orca(self) -> None:
        with mock.patch("run_reference_serial.subprocess.Popen") as popen:
            evidence = dry_run_validate(self.project, self.manifest, self.gate)
        self.assertEqual(len(evidence), 6)
        popen.assert_not_called()

    def test_execute_rejects_draft_before_any_orca_process(self) -> None:
        with mock.patch("run_reference_serial.subprocess.Popen") as popen:
            with self.assertRaises(WorkflowValidationError):
                run_reference_serial(
                    project=self.project,
                    manifest_path=self.manifest,
                    gate_path=self.gate,
                    execute=True,
                )
        popen.assert_not_called()

    def test_input_or_method_tamper_fails_closed(self) -> None:
        first = self._manifest_rows()[0]
        opt = Path(first["job_directory"]) / first["opt_input_file"]
        opt.write_text(opt.read_text(encoding="utf-8") + "# intentional tamper\n", encoding="utf-8")
        with self.assertRaises(WorkflowValidationError):
            dry_run_validate(self.project, self.manifest, self.gate)


if __name__ == "__main__":
    unittest.main()
