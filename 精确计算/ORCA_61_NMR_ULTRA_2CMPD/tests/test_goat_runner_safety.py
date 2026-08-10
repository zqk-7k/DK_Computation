from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
RUNNER_PATH = PROJECT / "scripts" / "run_goat_batch.py"
SPEC = importlib.util.spec_from_file_location("goat_runner_under_test", RUNNER_PATH)
assert SPEC is not None and SPEC.loader is not None
runner = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="")


class FakeKernel32:
    def __init__(self, initial_priority: int = 0x20, open_ok: bool = True):
        self.priority = initial_priority
        self.open_ok = open_ok
        self.closed = False

    def OpenProcess(self, access: int, inherit: bool, pid: int) -> int:
        return 1234 if self.open_ok else 0

    def GetPriorityClass(self, handle: int) -> int:
        return self.priority

    def SetPriorityClass(self, handle: int, value: int) -> int:
        self.priority = value
        return 1

    def CloseHandle(self, handle: int) -> int:
        self.closed = True
        return 1


class LegacyMigrationTests(unittest.TestCase):
    def make_case(self, root: Path, job_id: str = runner.LEGACY_TWO_FILE_PILOT_JOB_ID):
        source = root / "source"
        run_root = root / "runs"
        work = run_root / job_id / "attempt_001" / "work"
        write_text(source / "goat.inp", "! GOAT\n* xyzfile 0 1 start.xyz\n")
        write_text(source / "start.xyz", "1\nmodel\nH 0 0 0\n")
        write_text(source / "atom_map.csv", "atom\n1\n")
        write_text(source / "isotope_map.csv", "atom,mass\n")
        write_text(source / "structure_validation.json", "{}\n")
        write_text(work / "goat.inp", (source / "goat.inp").read_text(encoding="utf-8"))
        write_text(work / "start.xyz", (source / "start.xyz").read_text(encoding="utf-8"))
        job = runner.Job(
            row={}, job_id=job_id, figure_id="F3_001", source_directory=source
        )
        state = {
            "status": "completed",
            "input_bundle_sha256": "legacy-hash",
            "attempts": [
                {
                    "status": "completed",
                    "input_bundle_sha256": "legacy-attempt-hash",
                    "run_directory": str(work),
                }
            ],
        }
        files = [
            source / "goat.inp",
            source / "start.xyz",
            source / "atom_map.csv",
            source / "isotope_map.csv",
            source / "structure_validation.json",
        ]
        return job, state, files, run_root, work

    def test_named_pilot_migrates_only_after_per_file_sha_match(self):
        with tempfile.TemporaryDirectory() as temporary:
            job, state, files, run_root, _ = self.make_case(Path(temporary))
            record, reason = runner.verified_legacy_pilot_migration(
                job, state, files, "canonical-five-file-hash", run_root
            )
            self.assertIsNotNone(record, reason)
            assert record is not None
            self.assertEqual(record["verification_basis"], "PER_FILE_SHA256_MATCH")
            self.assertTrue(record["orca_physics_inputs_unchanged"])
            self.assertEqual(len(record["verified_physics_inputs"]), 2)
            for item in record["verified_physics_inputs"]:
                self.assertEqual(item["current_source_sha256"], item["legacy_staged_sha256"])

    def test_changed_physics_input_rejects_migration(self):
        with tempfile.TemporaryDirectory() as temporary:
            job, state, files, run_root, work = self.make_case(Path(temporary))
            write_text(work / "start.xyz", "1\nchanged\nH 1 0 0\n")
            record, reason = runner.verified_legacy_pilot_migration(
                job, state, files, "canonical-five-file-hash", run_root
            )
            self.assertIsNone(record)
            self.assertIn("differs", reason)

    def test_any_other_job_rejects_migration(self):
        with tempfile.TemporaryDirectory() as temporary:
            job, state, files, run_root, _ = self.make_case(
                Path(temporary), job_id="F3_002_DMSO_01_primary_s01"
            )
            record, reason = runner.verified_legacy_pilot_migration(
                job, state, files, "canonical-five-file-hash", run_root
            )
            self.assertIsNone(record)
            self.assertIn("restricted", reason)

    def test_work_directory_outside_job_root_rejects_migration(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            job, state, files, run_root, work = self.make_case(root)
            outside = root / "outside"
            write_text(outside / "goat.inp", (work / "goat.inp").read_text(encoding="utf-8"))
            write_text(outside / "start.xyz", (work / "start.xyz").read_text(encoding="utf-8"))
            state["attempts"][-1]["run_directory"] = str(outside)
            record, reason = runner.verified_legacy_pilot_migration(
                job, state, files, "canonical-five-file-hash", run_root
            )
            self.assertIsNone(record)
            self.assertIn("outside", reason)


class LockPriorityAndAggregateTests(unittest.TestCase):
    def test_process_global_lock_excludes_second_runner(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "runner_global.lock"
            first = runner.ExclusiveFileLock(path, stale_hours=192, wait_seconds=0)
            second = runner.ExclusiveFileLock(path, stale_hours=192, wait_seconds=0)
            self.assertTrue(first.acquire())
            self.assertFalse(second.acquire())
            first.release()
            self.assertTrue(second.acquire())
            second.release()

    def test_live_same_host_lock_never_becomes_age_stale(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "runner_global.lock"
            lock = runner.ExclusiveFileLock(path, stale_hours=0.000001, wait_seconds=0)
            self.assertTrue(lock.acquire())
            record = json.loads(path.read_text(encoding="utf-8"))
            record["created_utc"] = "2000-01-01T00:00:00Z"
            path.write_text(json.dumps(record), encoding="utf-8")
            self.assertFalse(lock._existing_is_stale())
            lock.release()

    def test_priority_is_verified_and_set_with_fake_windows_api(self):
        api = FakeKernel32(initial_priority=0x20)
        record = runner.ensure_windows_below_normal_priority(999, kernel32=api)
        self.assertEqual(record["before"], "NORMAL_PRIORITY_CLASS")
        self.assertEqual(record["after"], "BELOW_NORMAL_PRIORITY_CLASS")
        self.assertTrue(record["set_attempted"])
        self.assertTrue(record["set_succeeded"])
        self.assertTrue(record["verified"])
        self.assertEqual(record["warning"], "")
        self.assertTrue(api.closed)

    def test_priority_open_failure_is_a_nonfatal_recorded_warning(self):
        record = runner.ensure_windows_below_normal_priority(
            999, kernel32=FakeKernel32(open_ok=False)
        )
        self.assertFalse(record["verified"])
        self.assertIn("could not be verified", record["warning"])

    def test_ram_floor_cannot_be_weakened(self):
        self.assertEqual(runner.effective_minimum_free_ram_gib(1), 24.0)
        self.assertEqual(runner.effective_minimum_free_ram_gib(24), 24.0)
        self.assertEqual(runner.effective_minimum_free_ram_gib(32), 32.0)

    def test_default_mode_is_dry_run_and_does_not_take_global_lock(self):
        with tempfile.TemporaryDirectory() as temporary:
            state_root = Path(temporary) / "state"
            original = runner.execute_runner
            observed = {}

            def fake_execute(args):
                observed["execute"] = args.execute
                observed["lock_exists"] = (state_root / "runner_global.lock").exists()
                return 0

            runner.execute_runner = fake_execute
            try:
                code = runner.main(["--state-root", str(state_root)])
            finally:
                runner.execute_runner = original
            self.assertEqual(code, 0)
            self.assertFalse(observed["execute"])
            self.assertFalse(observed["lock_exists"])

    def test_main_holds_global_lock_for_the_execute_runner_lifetime(self):
        with tempfile.TemporaryDirectory() as temporary:
            state_root = Path(temporary) / "state"
            original = runner.execute_runner
            observed = {}

            def fake_execute(args):
                observed["lock_exists"] = (state_root / "runner_global.lock").is_file()
                observed["ram_floor"] = args.minimum_free_ram_gib
                return 0

            runner.execute_runner = fake_execute
            try:
                code = runner.main(
                    [
                        "--state-root",
                        str(state_root),
                        "--minimum-free-ram-gib",
                        "1",
                        "--execute",
                    ]
                )
            finally:
                runner.execute_runner = original
            self.assertEqual(code, 0)
            self.assertTrue(observed["lock_exists"])
            self.assertEqual(observed["ram_floor"], 24.0)
            self.assertFalse((state_root / "runner_global.lock").exists())

    def test_main_refuses_when_another_runner_owns_global_lock(self):
        with tempfile.TemporaryDirectory() as temporary:
            state_root = Path(temporary) / "state"
            owner = runner.ExclusiveFileLock(
                state_root / "runner_global.lock", stale_hours=192, wait_seconds=0
            )
            self.assertTrue(owner.acquire())
            original = runner.execute_runner
            called = False

            def should_not_execute(args):
                nonlocal called
                called = True
                return 0

            runner.execute_runner = should_not_execute
            try:
                code = runner.main(["--state-root", str(state_root), "--execute"])
            finally:
                runner.execute_runner = original
                owner.release()
            self.assertEqual(code, 3)
            self.assertFalse(called)

    def test_aggregate_separates_canonical_and_attempt_hashes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            job = runner.Job(
                row={}, job_id="F3_001_test", figure_id="F3_001", source_directory=root
            )
            store = runner.StateStore(root / "state", root / "project", [job])
            store.jobs_root.mkdir(parents=True, exist_ok=True)
            state = {
                "status": "completed",
                "attempt_count": 1,
                "input_bundle_sha256": "canonical-five-file-hash",
                "attempts": [
                    {
                        "attempt_no": 1,
                        "status": "completed",
                        "input_bundle_sha256": "legacy-two-file-hash",
                        "artifacts": {},
                    }
                ],
            }
            store.state_path(job).write_text(json.dumps(state), encoding="utf-8")
            row = store._row_for(job)
            self.assertEqual(
                row["canonical_state_input_bundle_sha256"], "canonical-five-file-hash"
            )
            self.assertEqual(
                row["latest_attempt_input_bundle_sha256"], "legacy-two-file-hash"
            )
            self.assertNotIn("input_bundle_sha256", runner.AGGREGATE_FIELDS)


class ProductionManifestStaticTests(unittest.TestCase):
    def test_all_production_rows_match_the_five_file_canonical_hash(self):
        manifest = PROJECT / "data" / "GOAT_explicit_job_manifest.csv"
        with manifest.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 200)
        self.assertEqual(len({row["job_id"] for row in rows}), 200)
        self.assertEqual({row["figure_id"] for row in rows}, {"F3_005", "F3_015"})
        self.assertEqual(len({(row["figure_id"], row["scenario_id"]) for row in rows}), 10)
        names = (
            "goat.inp",
            "start.xyz",
            "atom_map.csv",
            "isotope_map.csv",
            "structure_validation.json",
        )
        for row in rows:
            directory = Path(row["job_directory"])
            files = runner.source_input_files(directory, "goat.inp")
            self.assertEqual([path.name for path in files], list(names))
            actual = runner.input_bundle_sha256(files)
            self.assertEqual(actual.casefold(), row["input_bundle_sha256"].casefold())
            lock_hash = (directory / "INPUT_BUNDLE_SHA256.txt").read_text().strip()
            self.assertEqual(lock_hash.casefold(), actual.casefold())


if __name__ == "__main__":
    unittest.main()
