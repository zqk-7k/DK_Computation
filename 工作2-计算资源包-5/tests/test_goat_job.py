"""No quantum programs are launched: only temporary fixtures and dry-run tests."""
from pathlib import Path
import json
import os
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import goat_job as g


class GoatJobTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.folder = self.root / "goat" / "j1"
        self.folder.mkdir(parents=True)
        (self.root / "registry").mkdir()
        self.orca = self.root / ("orca.exe" if os.name == "nt" else "orca")
        self.xtb = self.root / ("otool_xtb.exe" if os.name == "nt" else "otool_xtb")
        self.orca.write_bytes(b"FAKE TEST FIXTURE; NEVER EXECUTE; Program Version 6.1.0 - RELEASE")
        self.xtb.write_bytes(b"FAKE TEST FIXTURE; NEVER EXECUTE")
        self.orca.chmod(0o700)
        self.xtb.chmod(0o700)
        (self.root / "libiomp5md.dll").write_bytes(b"FAKE TEST FIXTURE")
        (self.folder / "goat.inp").write_text("! GOAT-ENTROPY XTB ALPB(DMSO)\n%pal nprocs 8 end\n%maxcore 2000\n* xyzfile 0 1 start.xyz\n", encoding="utf-8")
        (self.folder / "start.xyz").write_text("1\nfixture only\nHe 0 0 0\n", encoding="utf-8")
        (self.folder / "cluster.mol").write_text("fixture only", encoding="utf-8")
        (self.folder / "atom_map.json").write_text("[]", encoding="utf-8")
        (self.folder / "provenance.json").write_text("{}", encoding="utf-8")
        self.row = {"job_id": "j1", "job_directory": "goat/j1", "initial_release": True,
                    "status": "PREPARED_NOT_RUN", **{field: g.digest(self.folder / name) for name, field in g.INPUT_HASH_FIELDS.items()}}
        self.save()
        self.mpi = patch.object(g.shutil, "which", return_value="fixture-mpiexec")
        self.mpi.start()
        self.addCleanup(self.mpi.stop)
        self.child = patch.object(g.subprocess, "Popen", side_effect=AssertionError("No subprocess may launch in tests"))
        self.child.start()
        self.addCleanup(self.child.stop)
        self.environment = patch.dict(os.environ, {}, clear=True)
        self.environment.start()
        self.addCleanup(self.environment.stop)

    def save(self):
        (self.root / "registry" / "jobs.json").write_text(json.dumps([self.row]), encoding="utf-8")

    def dry(self):
        return g.run_one(self.root, "j1", self.orca, self.xtb, acknowledge_models=True, dry_run=True)

    def test_preflight_reads_marker_without_execution(self):
        result = g.preflight(self.orca, self.xtb)
        self.assertTrue(result["passed"])
        self.assertEqual(result["runtime_goat_smoke"], "UNVALIDATED_NOT_EXECUTED")
        self.assertEqual(result["quantum_jobs_launched"], 0)

    def test_dry_run_writes_nothing(self):
        before = {p.relative_to(self.root): p.read_bytes() for p in self.root.rglob("*") if p.is_file()}
        result = self.dry()
        self.assertEqual(result["status"], "DRY_RUN_NO_EXECUTION")
        after = {p.relative_to(self.root): p.read_bytes() for p in self.root.rglob("*") if p.is_file()}
        self.assertEqual(before, after)

    def test_wrong_and_ambiguous_version_blocked(self):
        for content in (b"Program Version 6.1.1", b"Program Version 6.1.0 Program Version 6.1.1", b"no version"):
            self.orca.write_bytes(content)
            self.assertFalse(g.preflight(self.orca, self.xtb)["passed"])

    def test_missing_xtb_or_mpi_blocked(self):
        self.assertFalse(g.preflight(self.orca, self.root / "missing")["passed"])
        with patch.object(g.shutil, "which", return_value=None):
            self.assertFalse(g.preflight(self.orca, self.xtb)["passed"])

    def test_mislocated_xtb_blocked(self):
        wrong = self.root / "wrong_xtb"
        wrong.write_bytes(b"fixture")
        wrong.chmod(0o700)
        self.assertFalse(g.preflight(self.orca, wrong)["passed"])

    def test_acknowledgement_required_even_dry_run(self):
        with self.assertRaisesRegex(ValueError, "acknowledge"):
            g.run_one(self.root, "j1", self.orca, self.xtb, dry_run=True)

    def test_held_and_nonprepared_blocked(self):
        self.row["initial_release"] = False
        self.save()
        with self.assertRaisesRegex(ValueError, "Held"):
            self.dry()
        self.row["initial_release"] = True
        self.row["status"] = "RUNNING"
        self.save()
        with self.assertRaisesRegex(ValueError, "PREPARED_NOT_RUN"):
            self.dry()

    def test_each_bound_input_tamper_blocked(self):
        for name in g.INPUT_HASH_FIELDS:
            with self.subTest(name=name):
                path = self.folder / name
                original = path.read_bytes()
                path.write_bytes(original + b"altered")
                with self.assertRaisesRegex(ValueError, "hash mismatch"):
                    self.dry()
                path.write_bytes(original)

    def test_existing_output_receipt_and_lock_blocked(self):
        for name in ("goat.out", "goat.finalensemble.xyz", "receipt.json", ".goat-run.lock"):
            path = self.folder / name
            path.write_text("fixture", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Existing output"):
                self.dry()
            path.unlink()

    def test_escape_and_duplicate_id_blocked(self):
        self.row["job_directory"] = "../outside"
        self.save()
        with self.assertRaisesRegex(ValueError, "confined"):
            self.dry()
        (self.root / "registry" / "jobs.json").write_text(json.dumps([self.row, self.row]), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "duplicated"):
            self.dry()

    def test_xyz_fifth_column_rejected(self):
        path = self.folder / "start.xyz"
        path.write_text("1\nfixture\nH 0 0 0 M=2.01410177812\n", encoding="utf-8")
        self.row["xyz_sha256"] = g.digest(path)
        self.save()
        with self.assertRaisesRegex(ValueError, "exactly"):
            self.dry()

    def test_scheduler_resource_mismatch_blocked(self):
        with patch.dict(os.environ, {"SLURM_CPUS_PER_TASK": "4"}):
            with self.assertRaisesRegex(ValueError, "CPUs"):
                self.dry()


if __name__ == "__main__":
    unittest.main()
