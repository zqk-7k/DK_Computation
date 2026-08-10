from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


PACKAGE = Path(__file__).resolve().parents[1]
SCRIPT = PACKAGE / "scripts/activate_tms_reference_jobs.py"


def load_module():
    spec = importlib.util.spec_from_file_location("activate_tms_reference_jobs_test", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ReferenceActivationTests(unittest.TestCase):
    def test_activation_emits_only_requested_stage_rows_with_strict_tokens(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            tms = root / "reference_jobs/TMS_DMSO"
            production = tms / "production_revdsd_pcsseg3"
            crosscheck = tms / "crosscheck_pbe0_pcsseg3"
            production.mkdir(parents=True)
            crosscheck.mkdir(parents=True)
            geometry = "5\nTMS fixture\n" + "\n".join(
                [
                    "Si 0.000000 0.000000 0.000000",
                    "C 1.900000 0.000000 0.000000",
                    "C -1.900000 0.000000 0.000000",
                    "C 0.000000 1.900000 0.000000",
                    "C 0.000000 -1.900000 0.000000",
                ]
            ) + "\n"
            (tms / "tms_opt.xyz").write_text(geometry, encoding="utf-8")
            (tms / "tms_opt.out").write_text("ORCA TERMINATED NORMALLY\n", encoding="utf-8")
            manifest = root / "data/reference_job_manifest.csv"
            manifest.parent.mkdir(parents=True)
            fields = (
                "job_id",
                "stage",
                "method_id",
                "portable_job_directory",
                "coordinate_file",
                "enabled",
                "readiness",
                "approval_status",
            )
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(
                    [
                        {
                            "job_id": "TMS__opt",
                            "stage": "reference_optfreq",
                            "method_id": "PBE0_opt",
                            "portable_job_directory": "reference_jobs/TMS_DMSO",
                            "coordinate_file": "tms_start.xyz",
                            "enabled": "true",
                            "readiness": "ready",
                            "approval_status": "APPROVED_FOR_REFERENCE_GEOMETRY_CALCULATION",
                        },
                        {
                            "job_id": "TMS__prod",
                            "stage": "reference_shielding",
                            "method_id": "production_revdsd_pcsseg3",
                            "portable_job_directory": "reference_jobs/TMS_DMSO/production_revdsd_pcsseg3",
                            "coordinate_file": "tms_opt.xyz",
                            "enabled": "false",
                            "readiness": "blocked",
                            "approval_status": "DRAFT_NOT_APPROVED",
                        },
                        {
                            "job_id": "TMS__cross",
                            "stage": "reference_shielding",
                            "method_id": "crosscheck_pbe0_pcsseg3",
                            "portable_job_directory": "reference_jobs/TMS_DMSO/crosscheck_pbe0_pcsseg3",
                            "coordinate_file": "tms_opt.xyz",
                            "enabled": "false",
                            "readiness": "blocked",
                            "approval_status": "DRAFT_NOT_APPROVED",
                        },
                    ]
                )
            active = root / "data/reference_shielding_active_manifest.csv"
            module.PROJECT = root
            module.ROOT = tms
            module.MANIFEST = manifest
            module.ACTIVE_MANIFEST = active
            argv = [
                "activate_tms_reference_jobs.py",
                "--method-id",
                "production_revdsd_pcsseg3",
                "--reviewer",
                "TEST_REVIEWER",
                "--confirm-normal-termination",
                "--confirm-no-significant-imaginary-frequency",
            ]
            with patch.object(sys, "argv", argv):
                self.assertEqual(module.main(), 0)
            with active.open("r", encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["stage"], "reference_shielding")
            self.assertEqual(rows[0]["method_id"], "production_revdsd_pcsseg3")
            self.assertEqual(rows[0]["enabled"], "true")
            self.assertEqual(rows[0]["readiness"], "ready")
            self.assertEqual(
                rows[0]["approval_status"],
                "APPROVED_FOR_REFERENCE_SHIELDING_CALCULATION",
            )
            self.assertTrue((production / "tms_opt.xyz").is_file())
            self.assertFalse((crosscheck / "tms_opt.xyz").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
