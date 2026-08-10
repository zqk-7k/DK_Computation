from __future__ import annotations

import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


PROJECT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import process_nmr_predictions as processor  # noqa: E402
from nmr_method_validation import (  # noqa: E402
    MethodValidationError,
    validate_nmr_method_evidence,
)
from runner_provenance import (  # noqa: E402
    RunnerProvenanceError,
    canonical_row_hash,
    sha256_file,
    validate_runner_output,
)


METHOD_CASES = (
    (
        "production_revdsd_pcsseg3",
        "revDSD-PBEP86-D4/2021",
        "revDSD-PBEP86-D4/2021|pcSseg-3|GIAO|NoFrozenCore|ORCA_6.1.1",
    ),
    (
        "benchmark_ri_mp2_pcsseg3",
        "RI-MP2",
        "RI-MP2|pcSseg-3|GIAO|NoFrozenCore|ORCA_6.1.1",
    ),
    (
        "benchmark_dlpno_mp2_pcsseg3",
        "DLPNO-MP2",
        "DLPNO-MP2|pcSseg-3|GIAO|NoFrozenCore|ORCA_6.1.1",
    ),
)


def method_input(method_token: str, *, relaxed: bool = True) -> str:
    density = "Density Relaxed" if relaxed else "Density Unrelaxed"
    return (
        f"! {method_token} pcSseg-3 NMR NoFrozenCore\n"
        "%mp2\n"
        f"  {density}\n"
        "end\n"
        "* xyzfile 0 1 molecule.xyz\n"
    )


def method_output(input_text: str) -> str:
    return (
        input_text
        + "\nCHEMICAL SHIELDING SUMMARY (ppm)\n"
        + "****ORCA TERMINATED NORMALLY****\n"
    )


class CorrelatedMethodValidationTests(unittest.TestCase):
    def test_relaxed_density_and_method_identity_are_checked_for_all_mp2_families(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for ordinal, (method_id, token, fingerprint) in enumerate(METHOD_CASES):
                input_path = root / f"method_{ordinal}.inp"
                output_path = root / f"method_{ordinal}.out"
                input_text = method_input(token)
                input_path.write_text(input_text, encoding="utf-8")
                output_path.write_text(method_output(input_text), encoding="utf-8")
                result = validate_nmr_method_evidence(
                    input_path, output_path, method_id, fingerprint
                )
                self.assertEqual(result["relaxed_mp2_density"], "verified_in_input_and_output")

    def test_unrelaxed_density_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_path = root / "method.inp"
            output_path = root / "method.out"
            text = method_input("RI-MP2", relaxed=False)
            input_path.write_text(text, encoding="utf-8")
            output_path.write_text(method_output(text), encoding="utf-8")
            with self.assertRaisesRegex(MethodValidationError, "relaxed MP2 density"):
                validate_nmr_method_evidence(
                    input_path,
                    output_path,
                    "benchmark_ri_mp2_pcsseg3",
                    "RI-MP2|pcSseg-3|GIAO|NoFrozenCore|ORCA_6.1.1",
                )

    def test_vpt2_is_diagnostic_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_path = root / "vpt2.inp"
            output_path = root / "vpt2.out"
            input_path.write_text("! VPT2\n", encoding="utf-8")
            output_path.write_text("VPT2\n", encoding="utf-8")
            with self.assertRaisesRegex(MethodValidationError, "diagnostic-only"):
                validate_nmr_method_evidence(
                    input_path,
                    output_path,
                    "diagnostic_vpt2",
                    "VPT2|ORCA_6.1.1",
                )

    def test_gas_phase_autoci_is_not_a_production_delta_branch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_path = root / "autoci.inp"
            output_path = root / "autoci.out"
            input_path.write_text("! AUTOCI pcSseg-4 NMR\n", encoding="utf-8")
            output_path.write_text("AUTOCI CCSD(T)\n", encoding="utf-8")
            with self.assertRaisesRegex(MethodValidationError, "diagnostic-only"):
                validate_nmr_method_evidence(
                    input_path,
                    output_path,
                    "diagnostic_autoci_ccsdt_pcsseg4_gas",
                    "AUTOCI-CCSD(T)|pcSseg-4|GIAO|GAS_PHASE_DELTA_DIAGNOSTIC|ORCA_6.1.1",
                )


class PropertyWorkflowProfileProvenanceTests(unittest.TestCase):
    @staticmethod
    def _fixture(root: Path) -> tuple[dict[str, str], dict[str, object]]:
        digest = "A" * 64
        job = root / "property_job"
        job.mkdir()
        row = {
            "job_id": "F3_005__DMSO_02_W00_primary__conf01__shielding",
            "figure_id": "F3_005",
            "scenario_id": "DMSO_02_W00_primary",
            "method_fingerprint": "METHOD_FP",
            "workflow_profile": "PRACTICAL_HIGH_ACCURACY",
            "workflow_property_scope": "profile_all_conformers_PBE0_plus_cumulative95_revDSD",
            "workflow_profile_config_sha256": digest,
            "workflow_profile_manifest_sha256": digest,
            "workflow_profile_scenario_registry_sha256": digest,
            "source_dft_manifest_sha256": digest,
            "nmr_atom_map_sha256": digest,
            "source_isotope_map_sha256": digest,
            "source_opt_freq_state_sha256": digest,
            "approved_minimum_gate_sha256": digest,
            "conformer_xyz_sha256": digest,
            "job_directory": str(job),
        }
        payload: dict[str, object] = dict(row)
        payload["optimized_coordinate_sha256"] = payload.pop("conformer_xyz_sha256")
        (job / "provenance.json").write_text(json.dumps(payload), encoding="utf-8")
        return row, payload

    def test_property_row_and_provenance_are_profile_hash_bound(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            row, _payload = self._fixture(Path(temporary))
            validated = processor.validate_property_provenance(row)
            binding = processor.resolve_property_profile_binding([validated])
            self.assertEqual(binding["name"], "PRACTICAL_HIGH_ACCURACY")
            self.assertEqual(
                binding["binding_status"],
                "PASS_ROW_PROVENANCE_AND_ALL_PROPERTY_JOBS_IDENTICAL",
            )

    def test_property_profile_hash_substitution_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            row, payload = self._fixture(Path(temporary))
            payload["workflow_profile_manifest_sha256"] = "B" * 64
            (Path(row["job_directory"]) / "provenance.json").write_text(
                json.dumps(payload), encoding="utf-8"
            )
            with self.assertRaisesRegex(
                processor.ProcessingError, "workflow_profile_manifest_sha256"
            ):
                processor.validate_property_provenance(row)

    def test_mixed_property_profile_bindings_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            row, _payload = self._fixture(root)
            first = processor.validate_property_provenance(row)
            second_payload = dict(first["payload"])
            second_payload["source_dft_manifest_sha256"] = "B" * 64
            with self.assertRaisesRegex(
                processor.ProcessingError, "source_dft_manifest_sha256"
            ):
                processor.resolve_property_profile_binding(
                    [first, {"payload": second_payload}]
                )

    @staticmethod
    def _live_binding_fixture(root: Path) -> tuple[dict[str, str], Path]:
        profile = "PRACTICAL_HIGH_ACCURACY"
        config = root / "config" / "workflow_profiles.json"
        manifest = (
            root
            / "data"
            / "workflow_profiles"
            / profile
            / "GOAT_explicit_job_manifest.csv"
        )
        registry = manifest.with_name("explicit_scenario_registry.csv")
        dft_manifest = root / "data" / "DFT_explicit_job_manifest.csv"
        for path, text in (
            (config, "profile config\n"),
            (manifest, "profile manifest\n"),
            (registry, "scenario registry\n"),
            (dft_manifest, "dft manifest\n"),
        ):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding="utf-8")
        return (
            {
                "name": profile,
                "config_sha256": processor.sha256_file(config).upper(),
                "materialized_manifest_sha256": processor.sha256_file(manifest).upper(),
                "scenario_registry_sha256": processor.sha256_file(registry).upper(),
                "source_dft_manifest_sha256": processor.sha256_file(dft_manifest).upper(),
                "binding_status": "PASS_ROW_PROVENANCE_AND_ALL_PROPERTY_JOBS_IDENTICAL",
            },
            dft_manifest,
        )

    def test_live_profile_and_dft_files_are_hash_bound(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            binding, dft_manifest = self._live_binding_fixture(root)
            checked = processor.validate_live_workflow_profile_binding(
                binding, dft_manifest, root
            )
            self.assertEqual(
                checked["live_file_binding_status"],
                "PASS_DFT_CONFIG_PROFILE_MANIFEST_AND_REGISTRY",
            )

    def test_live_dft_manifest_substitution_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            binding, dft_manifest = self._live_binding_fixture(root)
            dft_manifest.write_text("substituted dft manifest\n", encoding="utf-8")
            with self.assertRaisesRegex(
                processor.ProcessingError, "source_dft_manifest_sha256"
            ):
                processor.validate_live_workflow_profile_binding(
                    binding, dft_manifest, root
                )

    def test_live_profile_config_substitution_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            binding, dft_manifest = self._live_binding_fixture(root)
            (root / "config" / "workflow_profiles.json").write_text(
                "substituted config\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(processor.ProcessingError, "config_sha256"):
                processor.validate_live_workflow_profile_binding(
                    binding, dft_manifest, root
                )

    def test_legacy_property_binding_remains_explicitly_marked_compatibility(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            row, payload = self._fixture(root)
            row["workflow_profile"] = "ULTRA_LEGACY_DEFAULT"
            payload["workflow_profile"] = "ULTRA_LEGACY_DEFAULT"
            for field in (
                "workflow_profile_config_sha256",
                "workflow_profile_manifest_sha256",
                "workflow_profile_scenario_registry_sha256",
            ):
                row[field] = ""
                payload[field] = ""
            (Path(row["job_directory"]) / "provenance.json").write_text(
                json.dumps(payload), encoding="utf-8"
            )
            validated = processor.validate_property_provenance(row)
            binding = processor.resolve_property_profile_binding([validated])
            self.assertEqual(
                binding["binding_status"],
                "LEGACY_COMPATIBILITY_NO_PROFILE_CONFIG_HASHES",
            )

    def test_legacy_binding_rejects_pseudo_profile_file_hashes(self) -> None:
        payload = {
            "workflow_profile": "ULTRA_LEGACY_DEFAULT",
            "workflow_profile_config_sha256": "A" * 64,
            "workflow_profile_manifest_sha256": "",
            "workflow_profile_scenario_registry_sha256": "",
            "source_dft_manifest_sha256": "B" * 64,
        }
        with self.assertRaisesRegex(
            processor.ProcessingError, "must not carry pseudo profile-file hashes"
        ):
            processor.resolve_property_profile_binding([{"payload": payload}])


class RunnerAndReferenceProvenanceTests(unittest.TestCase):
    @staticmethod
    def _local_runner_fixture(root: Path) -> tuple[dict[str, str], Path, Path, Path]:
        job = root / "reference_job"
        job.mkdir(parents=True)
        coordinate = job / "molecule.xyz"
        coordinate.write_text("1\nmodel\nH 0 0 0\n", encoding="utf-8")
        input_path = job / "tms_shield.inp"
        input_text = method_input("revDSD-PBEP86-D4/2021")
        input_path.write_text(input_text, encoding="utf-8")
        output = job / "tms_shield.out"
        output.write_text(method_output(input_text), encoding="utf-8")
        row = {
            "job_id": "TMS_DMSO__production_revdsd_pcsseg3",
            "stage": "reference_shielding",
            "method_id": "production_revdsd_pcsseg3",
            "method_fingerprint": (
                "revDSD-PBEP86-D4/2021|pcSseg-3|GIAO|NoFrozenCore|ORCA_6.1.1"
            ),
            "job_directory": str(job),
            "input_file": input_path.name,
            "output_file": output.name,
            "expected_input_sha256": sha256_file(input_path),
        }
        state_root = root / "state"
        state_root.mkdir()
        state_path = state_root / f"{row['job_id']}.json"
        state = {
            "job_id": row["job_id"],
            "stage": "reference_shielding",
            "status": "complete",
            "manifest_hash_validation": {
                "checks": [
                    {
                        "field": "expected_input_sha256",
                        "expected": sha256_file(input_path).lower(),
                        "actual": sha256_file(input_path).lower(),
                        "matches": True,
                    }
                ],
                "matches": True,
            },
            "input_bundle": {
                "files": [
                    {
                        "role": "input",
                        "path": str(input_path),
                        "sha256": sha256_file(input_path),
                    }
                ]
            },
            "output": str(output),
            "output_sha256": sha256_file(output),
            "required_artifacts": [],
            "required_artifacts_complete": True,
            "orca_version": "6.1.1",
        }
        state_path.write_text(json.dumps(state), encoding="utf-8")
        return row, state_root, state_path, coordinate

    def test_local_runner_output_hash_is_revalidated(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            row, state_root, _state_path, _coordinate = self._local_runner_fixture(root)
            result = validate_runner_output(row, "reference_shielding", root, state_root)
            output = Path(result["output_path"])
            self.assertEqual(result["output_sha256"], sha256_file(output))
            output.write_text(output.read_text(encoding="utf-8") + "tampered\n", encoding="utf-8")
            with self.assertRaisesRegex(RunnerProvenanceError, "output hash mismatch"):
                validate_runner_output(row, "reference_shielding", root, state_root)

    def test_hpc_runner_attempt_and_ledger_are_hash_bound(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            job = root / "job"
            job.mkdir()
            input_path = job / "shielding.inp"
            input_path.write_text("! PBE0 NMR\n", encoding="utf-8")
            row = {
                "job_id": "HPC_JOB",
                "job_directory": str(job),
                "input_file": input_path.name,
                "output_file": "shielding.out",
                "input_sha256": sha256_file(input_path),
            }
            attempt = root / "results" / "attempt_001"
            attempt.mkdir(parents=True)
            output = attempt / "shielding.out"
            output.write_text("****ORCA TERMINATED NORMALLY****\n", encoding="utf-8")
            metadata = attempt / "run_metadata.json"
            metadata.write_text(
                json.dumps({"job_id": row["job_id"], "stage": "shielding"}),
                encoding="utf-8",
            )
            ledger = attempt / "FILE_MANIFEST_SHA256.csv"
            ledger.write_text(
                "file,size_bytes,sha256\n"
                + f"{output.name},{output.stat().st_size},{sha256_file(output)}\n"
                + f"{metadata.name},{metadata.stat().st_size},{sha256_file(metadata)}\n",
                encoding="utf-8",
            )
            state_root = root / "state"
            state_root.mkdir()
            state_path = state_root / "HPC_JOB.json"
            state_path.write_text(
                json.dumps(
                    {
                        "status": "COMPLETED",
                        "job_id": "HPC_JOB",
                        "stage": "shielding",
                        "manifest_row_sha256": canonical_row_hash(row),
                        "input_sha256": sha256_file(input_path),
                        "hash_validation": {
                            "strict_checks_present": True,
                            "matches": True,
                        },
                        "normal_termination_marker": True,
                        "required_artifacts_ok": True,
                        "attempt_path": str(attempt),
                        "output_file": output.name,
                        "file_manifest_sha256": sha256_file(ledger),
                    }
                ),
                encoding="utf-8",
            )
            result = validate_runner_output(row, "shielding", root, state_root)
            self.assertEqual(result["runner_kind"], "portable_hpc")
            ledger.write_text(ledger.read_text(encoding="utf-8") + "\n", encoding="utf-8")
            with self.assertRaisesRegex(RunnerProvenanceError, "file-manifest hash mismatch"):
                validate_runner_output(row, "shielding", root, state_root)

    def test_processor_requires_approved_hash_bound_reference_registry(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            row, state_root, state_path, coordinate = self._local_runner_fixture(root)
            input_path = Path(row["job_directory"]) / row["input_file"]
            output = Path(row["job_directory"]) / row["output_file"]
            entries = []
            for nucleus in ("1H", "13C"):
                entries.append(
                    {
                        "reference_id": f"TMS_{nucleus}",
                        "approval_status": "APPROVED_FOR_NMR_PROPERTIES",
                        "nucleus": nucleus,
                        "method_id": row["method_id"],
                        "method_fingerprint": row["method_fingerprint"],
                        "input_file": str(input_path),
                        "coordinate_file": str(coordinate),
                        "output_file": str(output),
                        "input_sha256": sha256_file(input_path),
                        "coordinate_sha256": sha256_file(coordinate),
                        "output_sha256": sha256_file(output),
                        "runner_state_file": str(state_path),
                        "runner_state_sha256": sha256_file(state_path),
                        "runner_kind": "local_manifest_runner",
                        "vibrational_correction_policy": (
                            "DIAGNOSTIC_ONLY_NOT_APPLIED_TO_MAIN_DELTA"
                        ),
                        "vpt2_correction_applied_to_main_delta": False,
                        "reviewer": "Reviewer A",
                    }
                )
            registry = root / "reference_registry.json"
            registry.write_text(
                json.dumps(
                    {
                        "status": "APPROVED_FOR_NMR_PROPERTIES",
                        "reviewer": "Reviewer A",
                        "approved_at_utc": "2026-08-10T00:00:00Z",
                        "reference_jobs": entries,
                    }
                ),
                encoding="utf-8",
            )
            with patch.object(processor, "PROJECT", root):
                selected, approval = processor.validate_approved_reference_registry(
                    registry,
                    row,
                    row["method_id"],
                    row["method_fingerprint"],
                    None,
                    state_root,
                )
            self.assertEqual(selected, output)
            self.assertEqual(approval["registry_status"], "APPROVED_FOR_NMR_PROPERTIES")

            entries[0]["vpt2_correction_applied_to_main_delta"] = True
            registry.write_text(
                json.dumps(
                    {
                        "status": "APPROVED_FOR_NMR_PROPERTIES",
                        "reviewer": "Reviewer A",
                        "approved_at_utc": "2026-08-10T00:00:00Z",
                        "reference_jobs": entries,
                    }
                ),
                encoding="utf-8",
            )
            with patch.object(processor, "PROJECT", root):
                with self.assertRaisesRegex(processor.ProcessingError, "may not enter"):
                    processor.validate_approved_reference_registry(
                        registry,
                        row,
                        row["method_id"],
                        row["method_fingerprint"],
                        None,
                        state_root,
                    )


class WeightGateBindingTests(unittest.TestCase):
    def test_weight_table_is_bound_to_named_reviewer_approved_gate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gate = root / "minimum_boltzmann_approved.json"
            weight_sets = {}
            gate_sets = []
            for nucleus, temperature, values in (
                ("1H", 312.609, (0.7, 0.3)),
                ("13C", 304.0651, (0.4, 0.6)),
            ):
                conformers = {}
                gate_rows = []
                for conformer, energy, weight in zip(
                    ("conf01", "conf02"), (-100.0, -99.99), values
                ):
                    conformers[conformer] = {
                        "weight": weight,
                        "qrrho_gibbs_free_energy_hartree": energy,
                    }
                    gate_rows.append(
                        {
                            "conformer_id": conformer,
                            "weight": weight,
                            "qrrho_gibbs_free_energy_hartree": energy,
                        }
                    )
                weight_sets[nucleus] = {
                    "temperature_K": temperature,
                    "conformers": conformers,
                }
                gate_sets.append(
                    {
                        "condition_scenario_id": "DMSO_01",
                        "nucleus": nucleus,
                        "temperature_K": temperature,
                        "approval_status": "APPROVED_FOR_NMR",
                        "weights": gate_rows,
                    }
                )
            gate.write_text(
                json.dumps(
                    {
                        "approval_status": "APPROVED_FOR_NMR",
                        "figure_id": "F3_005",
                        "reviewer": "Reviewer B",
                        "reviewed_at": "2026-08-10T00:00:00Z",
                        "upstream_state_hashes": [
                            {"job_id": "conf01", "sha256": "B" * 64},
                            {"job_id": "conf02", "sha256": "C" * 64},
                        ],
                        "boltzmann_weight_sets": gate_sets,
                    }
                ),
                encoding="utf-8",
            )
            gate_hash = sha256_file(gate)
            for weight_set in weight_sets.values():
                weight_set["approved_gate_sha256"] = gate_hash
            manifests = []
            for ordinal in range(2):
                directory = root / f"property_{ordinal}"
                directory.mkdir()
                (directory / "provenance.json").write_text(
                    json.dumps({"approved_minimum_gate_file": str(gate)}),
                    encoding="utf-8",
                )
                manifests.append(
                    {
                        f"conf{ordinal + 1:02d}": {
                            "job_id": f"property_{ordinal}",
                            "job_directory": str(directory),
                            "approved_minimum_gate_sha256": gate_hash,
                            "source_opt_freq_state_sha256": (
                                "B" * 64 if ordinal == 0 else "C" * 64
                            ),
                        }
                    }
                )
            result = processor.validate_weight_gate_binding(
                "F3_005", "DMSO_01", weight_sets, *manifests
            )
            self.assertEqual(
                result["status"],
                "PASS_PROPERTY_AND_WEIGHT_TABLE_SHARE_APPROVED_GATE",
            )
            dft_records = {
                "conf01": {"runner_provenance": {"runner_state_sha256": "B" * 64}},
                "conf02": {"runner_provenance": {"runner_state_sha256": "C" * 64}},
            }
            state_result = processor.validate_property_dft_state_binding(
                dft_records, *manifests
            )
            self.assertEqual(
                state_result["status"],
                "PASS_PROPERTY_JOBS_BIND_CURRENT_TRUSTED_DFT_RUNNER_STATES",
            )
            weight_sets["1H"]["conformers"]["conf01"]["weight"] = 0.6
            with self.assertRaisesRegex(processor.ProcessingError, "value mismatch"):
                processor.validate_weight_gate_binding(
                    "F3_005", "DMSO_01", weight_sets, *manifests
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
