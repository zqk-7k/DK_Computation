from __future__ import annotations

import csv
import fnmatch
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
HPC = HERE.parent
BIN = HPC / "bin"
if str(BIN) not in sys.path:
    sys.path.insert(0, str(BIN))

from hpc_common import verify_materialized_dependency_receipt


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


class HpcLayerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        (self.root / "hpc/templates").mkdir(parents=True)
        for template in (HPC / "templates").glob("*.in"):
            shutil.copy2(template, self.root / "hpc/templates" / template.name)
        (self.root / "inputs/job1").mkdir(parents=True)
        self.input = self.root / "inputs/job1/goat.inp"
        (self.root / "inputs/job1/atom_map.csv").write_text(
            "atom_index_1based,element,component,source_atom_index_1based,environment_atom,exchangeable,include_in_analyte_nmr,isotope_mass_u\n"
            "1,H,analyte,1,false,false,true,\n",
            encoding="utf-8",
        )
        (self.root / "inputs/job1/isotope_map.csv").write_text(
            "atom_index_1based,element,mass_u,component,goat_xyz_mass_annotation,downstream_coordinate_source,downstream_inline_mass_syntax,orca_6_1_manual_verified\n"
            "1,H,2.01410177812,analyte,NONE_PLAIN_H,REVIEWED_GOAT_SELECTED_GEOMETRY_ONLY,H {x_A} {y_A} {z_A} M=2.01410177812,true\n",
            encoding="utf-8",
        )
        self.write_fake_input(nprocs=1)
        self.manifest = self.root / "manifest.csv"
        self.write_manifest(expected=sha256(self.input))
        self.site = self.root / "site.json"
        self.site.write_text(
            json.dumps(
                {
                    "orca_executable": sys.executable,
                    "orca_sha256": sha256(Path(sys.executable)),
                    "xtb_executable": sys.executable,
                    "python_executable": sys.executable,
                    "module_commands": [],
                    "results_root": "results/hpc",
                    "strict_hash": True,
                    "keep_failed_scratch": False,
                    "slurm": {"account": "", "partition": "", "qos": "", "extra_directives": []},
                    "pbs": {"project": "", "queue": "", "extra_directives": []},
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        self.resources = self.root / "resources.json"
        self.resources.write_text(
            json.dumps(
                {
                    "stages": {
                        "goat": {
                            "cpus": 1,
                            "memory_mb": 1024,
                            "walltime": "00:10:00",
                            "max_concurrent": 2,
                            "input_default": "fake.inp",
                            "required_artifact_globs": ["*.globalminimum.xyz", "*.finalensemble.xyz"],
                            "copy_back_globs": ["*.inp", "*.out", "*.stderr.log", "*.xyz"],
                        }
                    }
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        (self.root / "config").mkdir()
        shutil.copy2(HPC.parent / "config/resource_policy.json", self.root / "config/resource_policy.json")
        conditions = self.root / "config/nucleus_conditions.csv"
        shutil.copy2(HPC.parent / "config/nucleus_conditions.csv", conditions)
        with conditions.open("a", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                ["F3_TEST", "1H", "DMSO-d6", "298.15", "400.13", "VERIFIED_TEST_FIXTURE", "HPC unit test"]
            )
            writer.writerow(
                ["F3_TEST", "13C", "DMSO-d6", "298.15", "100.61", "VERIFIED_TEST_FIXTURE", "HPC unit test"]
            )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write_fake_input(self, nprocs: int) -> None:
        self.input.write_text(
            "# fake ORCA executable input interpreted by Python\n"
            f"# %pal nprocs {nprocs} end\n"
            "from pathlib import Path\n"
            "Path('fake.globalminimum.xyz').write_text('1\\n-100.000000000000 converged=true\\nH 0 0 0\\n')\n"
            "Path('fake.finalensemble.xyz').write_text('1\\n-100.000000000000 converged=true\\nH 0 0 0\\n')\n"
            "print('****ORCA TERMINATED NORMALLY****')\n",
            encoding="utf-8",
        )

    def write_manifest(self, expected: str) -> None:
        with self.manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=(
                    "job_id",
                    "figure_id",
                    "compound_number",
                    "display_name",
                    "stage",
                    "scenario_id",
                    "primary_solvent",
                    "continuum_keyword",
                    "explicit_solvent_count",
                    "charge",
                    "multiplicity",
                    "portable_job_directory",
                    "input_file",
                    "output_file",
                    "expected_input_sha256",
                    "enabled",
                    "readiness",
                    "approval_status",
                ),
            )
            writer.writeheader()
            writer.writerow(
                {
                    "job_id": "TEST_001",
                    "figure_id": "F3_TEST",
                    "compound_number": "1",
                    "display_name": "test hydrogen",
                    "stage": "goat",
                    "scenario_id": "DMSO_TEST",
                    "primary_solvent": "DMSO-d6",
                    "continuum_keyword": "DMSO",
                    "explicit_solvent_count": "0",
                    "charge": "0",
                    "multiplicity": "1",
                    "portable_job_directory": "inputs/job1",
                    "input_file": "goat.inp",
                    "output_file": "fake.out",
                    "expected_input_sha256": expected,
                    "enabled": "true",
                    "readiness": "READY",
                    "approval_status": "APPROVED_FOR_GOAT_CALCULATION",
                }
            )

    def run_task(self) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(BIN / "run_orca_task.py"),
                "--package-root",
                str(self.root),
                "--manifest",
                str(self.manifest),
                "--stage",
                "goat",
                "--index",
                "1",
                "--site-config",
                str(self.site),
                "--stage-config",
                str(self.resources),
                "--allocated-cpus",
                "1",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

    def test_run_resume_status_and_collect(self) -> None:
        first = self.run_task()
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        state_path = self.root / "state/tasks/goat/TEST_001.json"
        state = json.loads(state_path.read_text(encoding="utf-8"))
        self.assertEqual(state["status"], "COMPLETED")
        self.assertFalse(Path(state["attempt_path"]).is_absolute())
        attempts = list((self.root / "results/hpc/goat/TEST_001").glob("attempt_*"))
        self.assertEqual(len(attempts), 1)

        second = self.run_task()
        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
        self.assertIn("RESUME", second.stdout)
        attempts_after = list((self.root / "results/hpc/goat/TEST_001").glob("attempt_*"))
        self.assertEqual(len(attempts_after), 1)

        status = subprocess.run(
            [
                sys.executable,
                str(BIN / "status.py"),
                "--package-root",
                str(self.root),
                "--manifest",
                str(self.manifest),
                "--stage",
                "goat",
                "--require-complete",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        self.assertEqual(status.returncode, 0, status.stdout + status.stderr)
        collect = subprocess.run(
            [
                sys.executable,
                str(BIN / "collect_results.py"),
                "--package-root",
                str(self.root),
                "--manifest",
                str(self.manifest),
                "--stage",
                "goat",
                "--site-config",
                str(self.site),
                "--require-all",
                "--materialize-job-directories",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        self.assertEqual(collect.returncode, 0, collect.stdout + collect.stderr)
        self.assertTrue((self.root / "inputs/job1/fake.out").is_file())
        receipts = list((self.root / "inputs/job1").glob(".hpc_materialization__goat__TEST_001.json"))
        self.assertEqual(len(receipts), 1)
        receipt = json.loads(receipts[0].read_text(encoding="utf-8"))
        self.assertFalse(Path(receipt["attempt_path"]).is_absolute())
        dependency_row = {
            "depends_on_job_id": "TEST_001",
            "depends_on_stage": "goat",
            "depends_on_artifacts": '["fake.globalminimum.xyz","fake.finalensemble.xyz"]',
        }
        verified = verify_materialized_dependency_receipt(
            self.root,
            dependency_row,
            self.root / "inputs/job1",
            "synthetic_child",
        )
        self.assertEqual(verified["parent_job_id"], "TEST_001")
        (self.root / "inputs/job1/fake.globalminimum.xyz").write_text(
            "tampered\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(ValueError, "absent or altered"):
            verify_materialized_dependency_receipt(
                self.root,
                dependency_row,
                self.root / "inputs/job1",
                "synthetic_child",
            )

    def test_materialization_conflict_is_fail_closed(self) -> None:
        first = self.run_task()
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        conflict = self.root / "inputs/job1/fake.out"
        conflict.write_text("different bytes\n", encoding="utf-8")
        collect = subprocess.run(
            [
                sys.executable,
                str(BIN / "collect_results.py"),
                "--package-root",
                str(self.root),
                "--manifest",
                str(self.manifest),
                "--stage",
                "goat",
                "--site-config",
                str(self.site),
                "--materialize-job-directories",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        self.assertEqual(collect.returncode, 3, collect.stdout + collect.stderr)
        self.assertEqual(conflict.read_text(encoding="utf-8"), "different bytes\n")

    def test_missing_approval_is_fail_closed(self) -> None:
        with self.manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = list(reader.fieldnames or [])
            rows = list(reader)
        rows[0]["approval_status"] = ""
        with self.manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        result = self.run_task()
        self.assertEqual(result.returncode, 3, result.stdout + result.stderr)
        state = json.loads((self.root / "state/tasks/goat/TEST_001.json").read_text(encoding="utf-8"))
        self.assertEqual(state["status"], "BLOCKED_SCIENTIFIC_GATE")

    def test_submit_cannot_skip_preflight(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(BIN / "submit_stage.py"),
                "--package-root",
                str(self.root),
                "--manifest",
                str(self.manifest),
                "--stage",
                "goat",
                "--scheduler",
                "slurm",
                "--site-config",
                str(self.site),
                "--stage-config",
                str(self.resources),
                "--submit",
                "--skip-preflight",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("never with --submit", result.stdout + result.stderr)

    def test_hash_mismatch_is_blocked(self) -> None:
        self.write_manifest(expected="0" * 64)
        result = self.run_task()
        self.assertEqual(result.returncode, 4, result.stdout + result.stderr)
        state = json.loads((self.root / "state/tasks/goat/TEST_001.json").read_text(encoding="utf-8"))
        self.assertEqual(state["status"], "BLOCKED_HASH_VALIDATION")

    def test_strict_site_requires_exact_orca_binary_hash(self) -> None:
        site = json.loads(self.site.read_text(encoding="utf-8"))
        site["orca_sha256"] = ""
        self.site.write_text(json.dumps(site), encoding="utf-8")
        result = self.run_task()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("strict_hash=true requires orca_sha256", result.stdout + result.stderr)

    def test_resume_rejects_tampered_attempt_artifact(self) -> None:
        first = self.run_task()
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        state_path = self.root / "state/tasks/goat/TEST_001.json"
        state = json.loads(state_path.read_text(encoding="utf-8"))
        attempt = self.root / state["attempt_path"]
        (attempt / "fake.finalensemble.xyz").write_text("tampered\n", encoding="utf-8")
        second = self.run_task()
        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
        self.assertNotIn("RESUME", second.stdout)
        attempts = [p for p in (self.root / "results/hpc/goat/TEST_001").glob("attempt_*") if p.is_dir()]
        self.assertEqual(len(attempts), 2)

    def test_cpu_oversubscription_is_blocked(self) -> None:
        self.write_fake_input(nprocs=2)
        self.write_manifest(expected=sha256(self.input))
        result = self.run_task()
        self.assertNotEqual(result.returncode, 0)
        state = json.loads((self.root / "state/tasks/goat/TEST_001.json").read_text(encoding="utf-8"))
        self.assertEqual(state["status"], "FAILED_RUNNER")
        self.assertIn("requests more parallel workers", state["error"])

    def test_slurm_and_pbs_render_are_bounded_and_do_not_use_mpirun(self) -> None:
        with self.manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = list(reader.fieldnames or [])
            rows = list(reader)
        second = dict(rows[0])
        second["job_id"] = "TEST_002"
        with self.manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows([rows[0], second])
        for scheduler, suffix, array_line in (
            ("slurm", "sbatch", "#SBATCH --array=1-2%2"),
            ("pbs", "pbs", "#PBS -J 1-2%2"),
        ):
            result = subprocess.run(
                [
                    sys.executable,
                    str(BIN / "submit_stage.py"),
                    "--package-root",
                    str(self.root),
                    "--manifest",
                    str(self.manifest),
                    "--stage",
                    "goat",
                    "--scheduler",
                    scheduler,
                    "--site-config",
                    str(self.site),
                    "--stage-config",
                    str(self.resources),
                    "--array-range",
                    "1-2",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            scripts = sorted((self.root / "state/generated_submissions").glob(f"*.{suffix}"))
            self.assertTrue(scripts)
            text = scripts[-1].read_text(encoding="utf-8")
            self.assertIn(array_line, text)
            self.assertNotIn("mpirun", text.casefold())
            self.assertNotIn("mpiexec", text.casefold())
            self.assertIn("run_orca_task.py", text)

    def test_preflight_static_checks_pass_with_fake_scheduler(self) -> None:
        fake_bin = self.root / "fake_bin"
        fake_bin.mkdir()
        if os.name == "nt":
            scheduler = fake_bin / "sbatch.bat"
            scheduler.write_text("@exit /b 0\n", encoding="ascii")
        else:
            scheduler = fake_bin / "sbatch"
            scheduler.write_text("#!/bin/sh\nexit 0\n", encoding="ascii")
            scheduler.chmod(0o755)
        environment = os.environ.copy()
        environment["PATH"] = str(fake_bin) + os.pathsep + environment.get("PATH", "")
        result = subprocess.run(
            [
                sys.executable,
                str(BIN / "preflight.py"),
                "--package-root",
                str(self.root),
                "--manifest",
                str(self.manifest),
                "--stage",
                "goat",
                "--scheduler",
                "slurm",
                "--site-config",
                str(self.site),
                "--stage-config",
                str(self.resources),
                "--allow-non-linux",
                "--no-execute",
            ],
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Preflight PASS", result.stdout)

    def test_normalize_manifest_converts_internal_absolute_paths(self) -> None:
        generated = self.root / "generated.csv"
        with generated.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=("job_id", "job_directory", "input_file"))
            writer.writeheader()
            writer.writerow(
                {
                    "job_id": "TEST_001",
                    "job_directory": str((self.root / "inputs/job1").resolve()),
                    "input_file": "goat.inp",
                }
            )
        portable = self.root / "portable.csv"
        result = subprocess.run(
            [
                sys.executable,
                str(BIN / "normalize_manifest_paths.py"),
                "--package-root",
                str(self.root),
                "--input",
                str(generated),
                "--output",
                str(portable),
                "--stage",
                "goat",
                "--stage-config",
                str(self.resources),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with portable.open("r", encoding="utf-8-sig", newline="") as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(row["portable_job_directory"], "inputs/job1")
        self.assertEqual(row["input_file"], "goat.inp")
        self.assertEqual(row["expected_input_sha256"], sha256(self.input))
        self.assertEqual(row["enabled"], "false")
        self.assertFalse(Path(row["portable_job_directory"]).is_absolute())

    def test_goat_hpc_state_bridge_draft_gate_and_dft_builder(self) -> None:
        run = self.run_task()
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        (self.root / "data").mkdir(exist_ok=True)
        scenario_registry = self.root / "data/scenarios.csv"
        scenario_registry.write_text(
            "figure_id,compound_number,display_name,scenario_id,primary_solvent,continuum_keyword,explicit_solvent_count,minimum_goat_starts,species_smiles,assumed_pd,species_note,status\n"
            "F3_TEST,1,test hydrogen,DMSO_TEST,DMSO-d6,DMSO,0,1,[H],NOT_APPLICABLE,test fixture,AUTHOR_APPROVED_ASSUMPTION_FOR_CALCULATION\n",
            encoding="utf-8",
        )
        compat = self.root / "data/goat_gate_compatible.csv"
        bridge = subprocess.run(
            [
                sys.executable,
                str(BIN / "export_workflow_states.py"),
                "--package-root",
                str(self.root),
                "--stage",
                "goat",
                "--manifest",
                str(self.manifest),
                "--output-state-root",
                "logs/goat_explicit_batch/jobs",
                "--compat-manifest",
                str(compat),
            ],
            cwd=self.root,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        self.assertEqual(bridge.returncode, 0, bridge.stdout + bridge.stderr)
        legacy_state = json.loads(
            (self.root / "logs/goat_explicit_batch/jobs/TEST_001.json").read_text(encoding="utf-8")
        )
        self.assertEqual(legacy_state["status"], "completed")
        self.assertEqual(legacy_state["attempts"][-1]["status"], "completed")
        with compat.open("r", encoding="utf-8-sig", newline="") as handle:
            compat_row = next(csv.DictReader(handle))
        self.assertEqual(compat_row["job_directory"], "inputs/job1")

        draft = subprocess.run(
            [
                sys.executable,
                str(HPC.parent / "scripts/draft_goat_selection_gate.py"),
                "--project",
                str(self.root),
                "--manifest",
                str(compat),
                "--scenario-registry",
                str(scenario_registry),
                "--state-root",
                str(self.root / "logs/goat_explicit_batch/jobs"),
                "--output-root",
                str(self.root / "jobs"),
                "--figure-id",
                "F3_TEST",
            ],
            cwd=self.root,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        self.assertEqual(draft.returncode, 0, draft.stdout + draft.stderr)
        draft_path = self.root / "jobs/F3_TEST/goat_ensemble_selection/goat_ensemble_selection_DRAFT.json"
        gate = json.loads(draft_path.read_text(encoding="utf-8"))
        self.assertEqual(gate["approval_status"], "DRAFT_NOT_APPROVED")

        # Synthetic test fixture only: emulate a completed human review to exercise the builder contract.
        gate["approval_status"] = "APPROVED_FOR_OPT_FREQ"
        gate["selection_completeness_check"] = "PASS"
        gate["protonation_and_cluster_integrity_review"] = "PASS"
        gate["reviewer"] = "AUTOMATED_TEST_FIXTURE_NOT_A_SCIENTIFIC_APPROVAL"
        gate["reviewed_at"] = "2000-01-01T00:00:00Z"
        for selection in gate["selected_conformers"]:
            selection["approval_status"] = "APPROVED_FOR_OPT_FREQ"
            selection["connectivity_stereochemistry_check"] = "PASS"
            selection["cluster_composition_check"] = "PASS"
            selection["manual_extraction_check"] = "PASS"
        approved = draft_path.with_name("goat_ensemble_selection_approved.json")
        approved.write_text(json.dumps(gate, indent=2) + "\n", encoding="utf-8")
        dft_manifest = self.root / "data/DFT_fixture.csv"
        build = subprocess.run(
            [
                sys.executable,
                str(HPC.parent / "scripts/build_dft_jobs.py"),
                "--project",
                str(self.root),
                "--gate",
                str(approved),
                "--manifest",
                str(compat),
                "--scenario-registry",
                str(scenario_registry),
                "--resource-policy",
                str(self.root / "config/resource_policy.json"),
                "--output-root",
                str(self.root / "02_dft_opt_freq_explicit"),
                "--output-manifest",
                str(dft_manifest),
            ],
            cwd=self.root,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        self.assertEqual(build.returncode, 0, build.stdout + build.stderr)
        with dft_manifest.open("r", encoding="utf-8-sig", newline="") as handle:
            dft_row = next(csv.DictReader(handle))
        dft_input = Path(dft_row["job_directory"]) / dft_row["input_file"]
        self.assertIn("nprocs 16", dft_input.read_text(encoding="utf-8"))

    def test_packaged_goat_manifest_static_contract(self) -> None:
        package_root = HPC.parent
        manifest = package_root / "data/GOAT_explicit_job_manifest.csv"
        if not manifest.is_file():
            self.skipTest("Portable GOAT manifest is not present in this development checkout")
        with manifest.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 200)
        for row in rows:
            raw_directory = row["portable_job_directory"]
            self.assertNotRegex(raw_directory, r"^[A-Za-z]:[\\/]")
            directory = package_root / raw_directory
            input_path = directory / row["input_file"]
            coordinate_path = directory / row["coordinate_file"]
            self.assertTrue(input_path.is_file(), input_path)
            self.assertTrue(coordinate_path.is_file(), coordinate_path)
            self.assertEqual(sha256(input_path), row["expected_input_sha256"].upper())
            self.assertEqual(sha256(coordinate_path), row["expected_coordinate_sha256"].upper())
            digest = hashlib.sha256()
            for name in (
                "goat.inp",
                "start.xyz",
                "atom_map.csv",
                "isotope_map.csv",
                "structure_validation.json",
            ):
                digest.update(name.encode())
                digest.update((directory / name).read_bytes())
            self.assertEqual(digest.hexdigest().upper(), row["expected_input_bundle_sha256"].upper())
            input_text = input_path.read_text(encoding="utf-8")
            self.assertRegex(input_text, r"(?i)\bnprocs\s+8\b")
            self.assertRegex(input_text, r"(?i)\bNWorkers\s+8\b")
            self.assertEqual(row.get("approval_status"), "APPROVED_FOR_GOAT_CALCULATION")

    def test_high_level_and_vpt2_scheduler_stages_are_packaged(self) -> None:
        package_root = HPC.parent
        resources = json.loads(
            (package_root / "hpc/config/stage_resources.json").read_text(encoding="utf-8")
        )["stages"]
        for stage in (
            "j_method_benchmark",
            "high_level_benchmark",
            "vpt2_force_field",
            "vpt2_nmr",
        ):
            self.assertIn(stage, resources)
            self.assertGreater(int(resources[stage]["cpus"]), 0)
        for name in (
            "submit_high_level_benchmark.sh",
            "submit_j_method_benchmark.sh",
            "submit_vpt2_force_field.sh",
            "submit_vpt2_nmr.sh",
        ):
            self.assertTrue((package_root / "hpc/submit" / name).is_file())
        self.assertIn("*.vpt2", resources["vpt2_force_field"]["copy_back_globs"])
        self.assertIn("vpt2_nmr.vpt2", resources["vpt2_nmr"]["required_artifact_globs"])
        self.assertIn("*.autociresp_*", resources["high_level_benchmark"]["copy_back_globs"])

        expected = {
            "optfreq": {
                "dft_opt_freq.gbw",
                "dft_opt_freq.hess",
                "dft_opt_freq.property.txt",
                "dft_opt_freq.xyz",
            },
            "reference_optfreq": {
                "tms_opt.gbw",
                "tms_opt.hess",
                "tms_opt.property.txt",
                "tms_opt.xyz",
            },
        }
        for stage, names in expected.items():
            patterns = set(resources[stage]["required_artifact_globs"])
            self.assertEqual(patterns, names)
        for stage, resource in resources.items():
            input_name = resource["input_default"]
            self.assertFalse(
                any(fnmatch.fnmatch(input_name, pattern) for pattern in resource["required_artifact_globs"]),
                f"{stage} required artifacts can be falsely satisfied by its input file",
            )

    def test_row_level_autoci_outputs_and_helper_are_enforced_and_copied(self) -> None:
        job = self.root / "inputs/job1"
        benchmark = job / "benchmark.inp"
        benchmark.write_text(
            "# %pal nprocs 1 end\n"
            "from pathlib import Path\n"
            "Path('benchmark.gbw').write_text('gbw')\n"
            "Path('benchmark.property.txt').write_text('property')\n"
            "Path('benchmark.autociresp_0').write_text('response')\n"
            "print('****ORCA TERMINATED NORMALLY****')\n",
            encoding="utf-8",
        )
        self.resources.write_text(
            json.dumps(
                {
                    "stages": {
                        "high_level_benchmark": {
                            "cpus": 1,
                            "memory_mb": 1024,
                            "walltime": "00:10:00",
                            "max_concurrent": 1,
                            "input_default": "benchmark.inp",
                            "required_artifact_globs": [
                                "benchmark.gbw",
                                "benchmark.property.txt",
                            ],
                            "copy_back_globs": ["*.inp", "*.out", "*.stderr.log", "*.gbw", "*.property.txt"],
                        }
                    }
                }
            ),
            encoding="utf-8",
        )
        site = json.loads(self.site.read_text(encoding="utf-8"))
        site["helper_executables"] = {"orca_autoci": sys.executable}
        self.site.write_text(json.dumps(site), encoding="utf-8")
        fields = (
            "job_id",
            "stage",
            "portable_job_directory",
            "input_file",
            "output_file",
            "expected_input_sha256",
            "enabled",
            "readiness",
            "approval_status",
            "required_artifacts",
            "required_output_globs",
            "required_helpers",
        )
        with self.manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerow(
                {
                    "job_id": "AUTOCI_001",
                    "stage": "high_level_benchmark",
                    "portable_job_directory": "inputs/job1",
                    "input_file": "benchmark.inp",
                    "output_file": "benchmark.out",
                    "expected_input_sha256": sha256(benchmark),
                    "enabled": "true",
                    "readiness": "ready",
                    "approval_status": "APPROVED_FOR_HIGH_LEVEL_BENCHMARK",
                    "required_artifacts": '["benchmark.gbw","benchmark.property.txt"]',
                    "required_output_globs": '["benchmark.autociresp_*"]',
                    "required_helpers": '["orca_autoci"]',
                }
            )
        command = [
            sys.executable,
            str(BIN / "run_orca_task.py"),
            "--package-root",
            str(self.root),
            "--manifest",
            str(self.manifest),
            "--stage",
            "high_level_benchmark",
            "--index",
            "1",
            "--site-config",
            str(self.site),
            "--stage-config",
            str(self.resources),
            "--allocated-cpus",
            "1",
        ]
        result = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        state = json.loads(
            (self.root / "state/tasks/high_level_benchmark/AUTOCI_001.json").read_text(
                encoding="utf-8"
            )
        )
        attempt = self.root / state["attempt_path"]
        self.assertTrue((attempt / "benchmark.autociresp_0").is_file())
        ledger = (attempt / "FILE_MANIFEST_SHA256.csv").read_text(encoding="utf-8-sig")
        self.assertIn("benchmark.autociresp_0", ledger)
        (attempt / "benchmark.autociresp_0").unlink()
        rerun = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False
        )
        self.assertEqual(rerun.returncode, 0, rerun.stdout + rerun.stderr)
        self.assertNotIn("RESUME", rerun.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
