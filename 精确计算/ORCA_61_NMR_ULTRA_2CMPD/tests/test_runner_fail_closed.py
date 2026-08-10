from __future__ import annotations

import csv
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch


PACKAGE = Path(__file__).resolve().parents[1]
SCRIPT = PACKAGE / "scripts/run_orca_manifest.py"


def load_module():
    spec = importlib.util.spec_from_file_location("run_orca_manifest_fail_closed_test", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class GenericRunnerFailClosedTests(unittest.TestCase):
    def test_default_is_hash_validating_dry_run(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            job = root / "jobs/job1"
            job.mkdir(parents=True)
            input_path = job / "job.inp"
            input_path.write_text("! HF\n* xyz 0 1\nH 0 0 0\n*\n", encoding="utf-8")
            manifest = root / "manifest.csv"
            fields = (
                "job_id",
                "stage",
                "job_directory",
                "input_file",
                "output_file",
                "expected_input_sha256",
                "enabled",
                "readiness",
                "approval_status",
                "depends_on_job_id",
            )
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerow(
                    {
                        "job_id": "JOB_1",
                        "stage": "goat",
                        "job_directory": "jobs/job1",
                        "input_file": "job.inp",
                        "output_file": "job.out",
                        "expected_input_sha256": module.sha256_file(input_path),
                        "enabled": "true",
                        "readiness": "ready",
                        "approval_status": "APPROVED_FOR_GOAT_CALCULATION",
                        "depends_on_job_id": "",
                    }
                )
            module.PROJECT = root
            stream = io.StringIO()
            with patch.object(sys, "argv", ["run_orca_manifest.py", str(manifest), "--stage", "goat"]):
                with redirect_stdout(stream):
                    self.assertIsNone(module.main())
            self.assertIn('"mode": "DRY_RUN"', stream.getvalue())
            self.assertFalse((job / "job.out").exists())

    def test_draft_or_missing_gate_is_rejected(self):
        module = load_module()
        row = {
            "job_id": "JOB_1",
            "stage": "goat",
            "enabled": "true",
            "readiness": "ready",
            "approval_status": "DRAFT_NOT_APPROVED",
        }
        self.assertTrue(module.approval_problem(row, "goat"))
        row.pop("approval_status")
        self.assertTrue(module.approval_problem(row, "goat"))
        row.pop("enabled")
        self.assertFalse(module.enabled_row(row))

    def test_dependency_order_and_cycle_detection(self):
        module = load_module()
        first = {"job_id": "ff", "depends_on_job_id": ""}
        second = {"job_id": "nmr", "depends_on_job_id": "ff"}
        self.assertEqual([row["job_id"] for row in module.topological_rows([second, first])], ["ff", "nmr"])
        with self.assertRaisesRegex(ValueError, "cycle"):
            module.topological_rows(
                [
                    {"job_id": "a", "depends_on_job_id": "b"},
                    {"job_id": "b", "depends_on_job_id": "a"},
                ]
            )

    def test_cross_stage_dependency_requires_hash_bound_completion(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            module.PROJECT = root
            job = root / "jobs/vpt2"
            job.mkdir(parents=True)
            output = job / "vpt2_ff.out"
            output.write_text("****ORCA TERMINATED NORMALLY****\n", encoding="utf-8")
            hessian = job / "vpt2_ff.hess"
            vpt2 = job / "vpt2_ff.vpt2"
            hessian.write_text("hessian\n", encoding="utf-8")
            vpt2.write_text("vpt2\n", encoding="utf-8")
            parent_id = "PARENT_FF"
            status_path = root / "logs/vpt2_force_field" / f"{parent_id}.json"
            status_path.parent.mkdir(parents=True)
            status_path.write_text(
                json.dumps(
                    {
                        "status": "complete",
                        "stage": "vpt2_force_field",
                        "job_id": parent_id,
                        "input_bundle_sha256": "A" * 64,
                        "orca_executable_sha256": "B" * 64,
                        "orca_version": "6.1.1",
                        "output": str(output),
                        "output_sha256": module.sha256_file(output),
                        "required_artifacts_complete": True,
                        "required_artifacts": [
                            {
                                "name": path.name,
                                "path": str(path),
                                "sha256": module.sha256_file(path),
                            }
                            for path in (hessian, vpt2)
                        ],
                    }
                ),
                encoding="utf-8",
            )
            child = {
                "job_id": "CHILD_NMR",
                "stage": "vpt2_nmr",
                "job_directory": "jobs/vpt2",
                "depends_on_job_id": parent_id,
                "depends_on_stage": "vpt2_force_field",
                "depends_on_artifacts": '["vpt2_ff.hess","vpt2_ff.vpt2"]',
            }
            self.assertEqual(module.topological_rows([child]), [child])
            hessian.write_text("tampered\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Unverified cross-stage dependency"):
                module.topological_rows([child])


if __name__ == "__main__":
    unittest.main(verbosity=2)
