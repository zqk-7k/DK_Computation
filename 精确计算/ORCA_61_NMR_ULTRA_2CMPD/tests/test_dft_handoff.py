from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "scripts"))

from build_dft_jobs import build_dft_jobs  # noqa: E402
from build_property_jobs import (  # noqa: E402
    build_property_jobs,
    shielding_branches,
    shielding_fingerprint,
)
from draft_goat_selection_gate import build_draft_gates  # noqa: E402
from draft_minimum_boltzmann_gate import build_draft_minimum_gates  # noqa: E402
from run_orca_manifest import build_input_bundle  # noqa: E402
from run_opt_freq_serial import ExclusiveLock, dry_run_validate  # noqa: E402
from workflow_schema import (  # noqa: E402
    DEUTERIUM_MASS_U,
    WorkflowValidationError,
    runner_input_bundle_sha256,
    sha256_file,
    source_input_files,
)


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def frame_text(symbols: list[str], xyz: list[tuple[float, float, float]], energy: float) -> str:
    lines = [str(len(symbols)), f"{energy:.12f} converged=true"]
    lines.extend(
        f"{symbol:<2s} {x: .12f} {y: .12f} {z: .12f}"
        for symbol, (x, y, z) in zip(symbols, xyz, strict=True)
    )
    return "\n".join(lines) + "\n"


class ExplicitDftHandoffTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.project = Path(self.temp.name) / "project"
        self.project.mkdir()
        self.manifest = self.project / "data" / "GOAT_explicit_job_manifest.csv"
        self.registry = self.project / "data" / "explicit_scenario_registry.csv"
        self.state_root = self.project / "logs" / "goat_explicit_batch" / "jobs"
        self.gate_root = self.project / "jobs"
        self.resource_policy = self.project / "config" / "resource_policy.json"
        write_csv(
            self.project / "config" / "nucleus_conditions.csv",
            [
                {
                    "figure_id": "F3_001",
                    "nucleus": nucleus,
                    "solvent": "D2O",
                    "temperature_K": temperature,
                    "spectrometer_frequency_MHz": frequency,
                    "condition_status": "LOCKED_TEST_FIXTURE",
                    "source_note": "synthetic unit-test fixture",
                }
                for nucleus, frequency, temperature in (
                    ("1H", "400.0", "310.0"),
                    ("13C", "100.0", "300.0"),
                )
            ],
            [
                "figure_id",
                "nucleus",
                "solvent",
                "temperature_K",
                "spectrometer_frequency_MHz",
                "condition_status",
                "source_note",
            ],
        )
        write(
            self.resource_policy,
            json.dumps(
                {
                    "execution": {
                        "one_compound_at_a_time": True,
                        "absolute_mpi_rank_cap": 8,
                    },
                    "stage_resources": {
                        "opt_freq": {"mpi_ranks": 6, "maxcore_mb_per_rank": 4000},
                        "shielding": {"mpi_ranks": 6, "maxcore_mb_per_rank": 4000},
                        "j_coupling": {"mpi_ranks": 6, "maxcore_mb_per_rank": 4500},
                    },
                },
                indent=2,
            )
            + "\n",
        )
        write(
            self.project / "data" / "real_nmr_calibration_audit" / "calibration.json",
            '{"role":"fixture only"}\n',
        )
        self.rows: list[dict[str, object]] = []
        scenarios = [("D2O_01_primary", 1), ("D2O_02_primary", 2)]
        registry_rows: list[dict[str, object]] = []
        for scenario_id, water_count in scenarios:
            registry_rows.append(
                {
                    "figure_id": "F3_001",
                    "compound_number": 1,
                    "display_name": "formamide fixture",
                    "scenario_id": scenario_id,
                    "primary_solvent": "D2O",
                    "continuum_keyword": "WATER",
                    "explicit_solvent_count": water_count,
                    "minimum_goat_starts": 1,
                    "species_smiles": "NC=O",
                    "assumed_pd": "7.0",
                    "species_note": "synthetic unit-test fixture",
                    "status": "TEST_ONLY",
                }
            )
            self.rows.append(self._make_completed_job(scenario_id, water_count))
        write_csv(self.manifest, self.rows, list(self.rows[0].keys()))
        write_csv(self.registry, registry_rows, list(registry_rows[0].keys()))

    def tearDown(self) -> None:
        self.temp.cleanup()

    def _make_completed_job(self, scenario_id: str, water_count: int) -> dict[str, object]:
        job_id = f"F3_001_{scenario_id}_s01"
        job_dir = self.project / "01_goat_explicit" / "F3_001" / scenario_id / "start_01"
        symbols = ["N", "C", "O", "H", "H", "H"]
        xyz = [
            (-0.70, 0.00, 0.00),
            (0.55, 0.00, 0.00),
            (1.20, 1.00, 0.00),
            (-1.10, 0.85, 0.00),
            (-1.10, -0.85, 0.00),
            (1.05, -0.95, 0.00),
        ]
        components = ["analyte"] * 6
        source_indices: list[object] = [1, 2, 3, 4, 5, 6]
        exchangeable = ["false", "false", "false", "true", "true", "false"]
        include = ["true"] * 6
        for number in range(1, water_count + 1):
            base = 3.0 * number
            symbols.extend(["O", "H", "H"])
            xyz.extend([(base, 0.0, 0.0), (base + 0.75, 0.5, 0.0), (base - 0.75, 0.5, 0.0)])
            components.extend([f"D2O_{number}"] * 3)
            source_indices.extend(["", "", ""])
            exchangeable.extend(["false", "true", "true"])
            include.extend(["false", "false", "false"])
        write(job_dir / "start.xyz", frame_text(symbols, xyz, -20.0))
        write(
            job_dir / "goat.inp",
            "! XTB2 GOAT\n%pal nprocs 1 end\n* xyzfile 0 1 start.xyz\n",
        )
        atom_rows = []
        isotope_rows = []
        deuterium_indices = {4, 5}
        deuterium_indices.update(index for index in range(7, len(symbols) + 1) if symbols[index - 1] == "H")
        for index, symbol in enumerate(symbols, start=1):
            atom_rows.append(
                {
                    "atom_index_1based": index,
                    "element": symbol,
                    "component": components[index - 1],
                    "source_atom_index_1based": source_indices[index - 1],
                    "environment_atom": "false" if include[index - 1] == "true" else "true",
                    "exchangeable": exchangeable[index - 1],
                    "include_in_analyte_nmr": include[index - 1],
                    "isotope_mass_u": (
                        f"{DEUTERIUM_MASS_U:.11f}" if index in deuterium_indices else ""
                    ),
                }
            )
            if index in deuterium_indices:
                isotope_rows.append(
                    {
                        "atom_index_1based": index,
                        "element": "H",
                        "mass_u": f"{DEUTERIUM_MASS_U:.11f}",
                        "component": components[index - 1],
                        "goat_xyz_mass_annotation": "NONE_PLAIN_H",
                        "downstream_coordinate_source": "REVIEWED_GOAT_SELECTED_GEOMETRY_ONLY",
                        "downstream_inline_mass_syntax": f"H {{x_A}} {{y_A}} {{z_A}} M={DEUTERIUM_MASS_U:.11f}",
                        "orca_6_1_manual_verified": "true",
                    }
                )
        write_csv(
            job_dir / "atom_map.csv",
            atom_rows,
            [
                "atom_index_1based",
                "element",
                "component",
                "source_atom_index_1based",
                "environment_atom",
                "exchangeable",
                "include_in_analyte_nmr",
                "isotope_mass_u",
            ],
        )
        write_csv(
            job_dir / "isotope_map.csv",
            isotope_rows,
            [
                "atom_index_1based",
                "element",
                "mass_u",
                "component",
                "goat_xyz_mass_annotation",
                "downstream_coordinate_source",
                "downstream_inline_mass_syntax",
                "orca_6_1_manual_verified",
            ],
        )
        bundle_hash = runner_input_bundle_sha256(source_input_files(job_dir))
        run_dir = self.project / "02_goat_explicit_runs" / job_id / "attempt_001" / "work"
        output = run_dir / "goat.out"
        ensemble = run_dir / "goat.finalensemble.xyz"
        minimum = run_dir / "goat.globalminimum.xyz"
        write(output, "synthetic fixture\nORCA TERMINATED NORMALLY\n")
        write(
            ensemble,
            frame_text(symbols, xyz, -25.0)
            + frame_text(
                symbols,
                [(x, y + (0.40 if index == 0 else 0.0), z) for index, (x, y, z) in enumerate(xyz)],
                -24.998,
            ),
        )
        write(minimum, frame_text(symbols, xyz, -25.0))
        state_path = self.state_root / f"{job_id}.json"
        state = {
            "schema_version": 1,
            "job_id": job_id,
            "figure_id": "F3_001",
            "source_job_directory": str(job_dir),
            "status": "completed",
            "attempt_count": 1,
            "attempts": [
                {
                    "attempt_no": 1,
                    "status": "completed",
                    "return_code": 0,
                    "input_bundle_sha256": bundle_hash,
                    "source_input_unchanged": True,
                    "source_input_bundle_sha256_after": bundle_hash,
                    "output_file": str(output),
                    "output_sha256": sha256_file(output),
                    "artifacts": {
                        "globalminimum": str(minimum),
                        "finalensemble": str(ensemble),
                    },
                }
            ],
            "reason": "",
            "input_bundle_sha256": bundle_hash,
        }
        write(state_path, json.dumps(state, indent=2) + "\n")
        return {
            "queue_order": 1,
            "job_id": job_id,
            "figure_id": "F3_001",
            "compound_number": 1,
            "display_name": "formamide fixture",
            "scenario_id": scenario_id,
            "primary_solvent": "D2O",
            "continuum_keyword": "WATER",
            "explicit_solvent_count": water_count,
            "start_index": 1,
            "charge": 0,
            "multiplicity": 1,
            "preparation_seed": 1,
            "job_directory": str(job_dir),
            "goat_input_sha256": sha256_file(job_dir / "goat.inp"),
            "start_xyz_sha256": sha256_file(job_dir / "start.xyz"),
            "input_bundle_sha256": bundle_hash,
            "notice": "CALCULATED_MODEL_NOT_EXPERIMENTAL",
        }

    def _draft_and_approve(self) -> tuple[Path, Path]:
        drafts = build_draft_gates(
            project=self.project,
            manifest_path=self.manifest,
            scenario_registry_path=self.registry,
            state_root=self.state_root,
            output_root=self.gate_root,
            figure_ids=["F3_001"],
            energy_window_kcal_mol=12.0,
            duplicate_rmsd_A=0.05,
        )
        self.assertEqual(len(drafts), 1)
        draft_path = drafts[0]
        draft = json.loads(draft_path.read_text(encoding="utf-8"))
        self.assertEqual(draft["approval_status"], "DRAFT_NOT_APPROVED")
        self.assertTrue(all(item["approval_status"] == "DRAFT_NOT_APPROVED" for item in draft["selected_conformers"]))
        draft["approval_status"] = "APPROVED_FOR_OPT_FREQ"
        draft["selection_completeness_check"] = "PASS"
        draft["protonation_and_cluster_integrity_review"] = "PASS"
        draft["reviewer"] = "UNIT_TEST_REVIEWER"
        draft["reviewed_at"] = "2026-07-20T00:00:00Z"
        for selection in draft["selected_conformers"]:
            selection["approval_status"] = "APPROVED_FOR_OPT_FREQ"
            selection["connectivity_stereochemistry_check"] = "PASS"
            selection["cluster_composition_check"] = "PASS"
            selection["manual_extraction_check"] = "PASS"
        approved = draft_path.with_name("goat_ensemble_selection_approved.json")
        write(approved, json.dumps(draft, indent=2) + "\n")
        return draft_path, approved

    def _build_dft_and_fake_completed_opt(self) -> tuple[list[dict[str, object]], Path]:
        _draft, approved = self._draft_and_approve()
        manifest_path = self.project / "data" / "DFT_explicit_job_manifest.csv"
        rows = build_dft_jobs(
            self.project,
            [approved],
            self.manifest,
            self.registry,
            self.resource_policy,
            self.project / "02_dft_opt_freq_explicit",
            manifest_path,
        )
        validated = dry_run_validate(self.project, manifest_path)
        self.assertEqual(len(validated), len(rows))
        state_root = self.project / "logs" / "opt_freq_serial" / "jobs"
        for row in rows:
            job_id = str(row["job_id"])
            job_dir = Path(str(row["job_directory"]))
            ordinal = int(row["conformer_ordinal"])
            scenario_offset = 0.0 if row["scenario_id"] == "D2O_01_primary" else -10.0
            electronic = -100.100000 - scenario_offset - 0.001 * ordinal
            gibbs = -100.000000 - scenario_offset - 0.0005 * ordinal
            run_dir = self.project / "02_dft_opt_freq_runs" / job_id / "attempt_fixture" / "work"
            output = run_dir / "dft_opt_freq.out"
            optimized = run_dir / "dft_opt_freq.xyz"
            hess = run_dir / "dft_opt_freq.hess"
            gbw = run_dir / "dft_opt_freq.gbw"
            prop = run_dir / "dft_opt_freq.property.txt"
            write(
                output,
                "\n".join(
                    [
                        f"FINAL SINGLE POINT ENERGY       {electronic:.12f}",
                        "Vibrational entropy computed according to the QRRHO of S. Grimme",
                        "Chem.Eur.J. 2012 18 9955 using a reference frequency of 100.0 cm-1",
                        "     0:       0.00 cm**-1",
                        "     1:       0.00 cm**-1",
                        "     2:       0.00 cm**-1",
                        "     3:       0.00 cm**-1",
                        "     4:       0.00 cm**-1",
                        "     5:       0.00 cm**-1",
                        "     6:     100.00 cm**-1",
                        "     7:    1200.00 cm**-1",
                        "Temperature         ...   300.00 K",
                        f"Final Gibbs free energy         ...    {gibbs:.12f} Eh",
                        "Temperature         ...   310.00 K",
                        f"Final Gibbs free energy         ...    {gibbs + 0.0001 * ordinal:.12f} Eh",
                        "****ORCA TERMINATED NORMALLY****",
                        "",
                    ]
                ),
            )
            optimized.parent.mkdir(parents=True, exist_ok=True)
            optimized.write_bytes((job_dir / "approved_selected_cluster.xyz").read_bytes())
            write(hess, "$hessian\nfixture\n$end\n")
            write(gbw, "fixture gbw\n")
            write(prop, "fixture property\n")
            artifacts = {
                role: {
                    "path": str(path),
                    "exists_nonempty": True,
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
                for role, path in {
                    "gbw": gbw,
                    "hess": hess,
                    "property": prop,
                    "optimized_xyz": optimized,
                }.items()
            }
            state = {
                "schema_version": 2,
                "notice": "CALCULATED MODEL - NOT EXPERIMENTAL DATA",
                "job_id": job_id,
                "figure_id": row["figure_id"],
                "scenario_id": row["scenario_id"],
                "source_job_directory": str(job_dir),
                "status": "completed",
                "input_bundle_sha256": row["expected_input_bundle_sha256"],
                "provenance_sha256": sha256_file(job_dir / "provenance.json"),
                "attempt_count": 1,
                "attempts": [
                    {
                        "attempt_no": 1,
                        "status": "completed",
                        "return_code": 0,
                        "normal_termination_marker": True,
                        "source_input_unchanged": True,
                        "input_sha256": row["input_sha256"],
                        "input_bundle_sha256": row["expected_input_bundle_sha256"],
                        "output_file": str(output),
                        "output_sha256": sha256_file(output),
                        "artifacts": artifacts,
                    }
                ],
            }
            write(state_root / f"{job_id}.json", json.dumps(state, indent=2) + "\n")
        return rows, manifest_path

    def _draft_and_approve_minimum(self) -> tuple[Path, Path, list[dict[str, object]], Path]:
        rows, manifest = self._build_dft_and_fake_completed_opt()
        drafts = build_draft_minimum_gates(
            project=self.project,
            manifest_path=manifest,
            state_root=self.project / "logs" / "opt_freq_serial" / "jobs",
            output_root=self.project / "jobs",
            figure_ids=["F3_001"],
        )
        draft_path = drafts[0]
        payload = json.loads(draft_path.read_text(encoding="utf-8"))
        self.assertEqual(payload["approval_status"], "DRAFT_NOT_APPROVED")
        payload["approval_status"] = "APPROVED_FOR_NMR"
        payload["same_composition_only_check"] = "PASS"
        payload["reviewer"] = "UNIT_TEST_REVIEWER"
        payload["reviewed_at"] = "2026-07-20T00:00:00Z"
        for decision in payload["conformer_decisions"]:
            decision["stationary_point_check"] = "PASS"
            decision["cluster_integrity_check"] = "PASS"
            decision["approval_status"] = "APPROVED_FOR_NMR"
        for weight_set in payload["boltzmann_weight_sets"]:
            weight_set["approval_status"] = "APPROVED_FOR_NMR"
        approved = draft_path.with_name("minimum_boltzmann_approved.json")
        write(approved, json.dumps(payload, indent=2) + "\n")
        return draft_path, approved, rows, manifest

    def _approved_reference_registry(self, mismatched: bool = False) -> Path:
        ref_input = self.project / "reference_fixture" / "reference.inp"
        ref_output = self.project / "reference_fixture" / "reference.out"
        write(ref_input, "reference fixture input\n")
        write(ref_output, "reference fixture output\n****ORCA TERMINATED NORMALLY****\n")
        entries = []
        for branch in shielding_branches("WATER"):
            fingerprint = branch["fingerprint"]
            for nucleus in ("1H", "13C"):
                entries.append(
                    {
                        "reference_id": f"D2O_{nucleus}_{branch['branch']}_FIXTURE",
                        "approval_status": "APPROVED_FOR_NMR_PROPERTIES",
                        "nucleus": nucleus,
                        "solvent_domain": "D2O",
                        "raw_calculation_reference": "fixture",
                        "experimental_reference_convention": "fixture",
                        "method_fingerprint": (
                            fingerprint + "|INTENTIONAL_MISMATCH" if mismatched else fingerprint
                        ),
                        "input_file": str(ref_input),
                        "output_file": str(ref_output),
                        "input_sha256": sha256_file(ref_input),
                        "output_sha256": sha256_file(ref_output),
                    }
                )
        path = self.project / "config" / (
            "reference_registry_mismatch.json" if mismatched else "reference_registry_approved.json"
        )
        write(
            path,
            json.dumps(
                {
                    "schema_version": "2.0",
                    "status": "APPROVED_FOR_NMR_PROPERTIES",
                    "reference_jobs": entries,
                },
                indent=2,
            )
            + "\n",
        )
        return path

    def _approved_j_gate(self) -> Path:
        files = {}
        for method in ("production", "benchmark"):
            for role in ("input", "output"):
                path = self.project / "j_benchmark_fixture" / f"{method}_{role}.txt"
                write(
                    path,
                    f"{method} {role} fixture\n"
                    + ("****ORCA TERMINATED NORMALLY****\n" if role == "output" else ""),
                )
                files[f"{method}_{role}_file"] = str(path)
                files[f"{method}_{role}_sha256"] = sha256_file(path)
        path = self.project / "jobs" / "j_method_benchmark_approved.json"
        write(
            path,
            json.dumps(
                {
                    "schema_version": "2.0",
                    "approval_status": "APPROVED_FOR_NMR_PROPERTIES",
                    "stage": "j_method_benchmark",
                    "production_method_fingerprint": "PBE|pcJ-3|AutoAux|CPCM(actual_solvent)|VeryTightSCF|DefGrid3|EPRNMR_Tol=1e-10|Z_Tol=1e-12|SpinSpinRThresh=6.05A|ORCA_6.1.1",
                    "benchmark_method_fingerprint": "PBE0|pcJ-2|AutoAux|CPCM(actual_solvent)|VeryTightSCF|DefGrid3|EPRNMR_Tol=1e-10|Z_Tol=1e-12|SpinSpinRThresh=6.05A|ORCA_6.1.1",
                    "representative_set_coverage_check": "PASS",
                    "numerical_agreement_review": "PASS",
                    "production_method_lock_check": "PASS",
                    "pair_screen_sensitivity_review": "PASS",
                    "reviewer": "UNIT_TEST_REVIEWER",
                    "reviewed_at": "2026-07-20T00:00:00Z",
                    "representative_case_results": [
                        {"approval_status": "APPROVED_FOR_NMR_PROPERTIES", **files}
                    ],
                },
                indent=2,
            )
            + "\n",
        )
        return path

    def _property_build(
        self,
        gate: Path,
        manifest: Path,
        reference_registry: Path,
        tag: str,
    ) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
        return build_property_jobs(
            project=self.project,
            gate_paths=[gate],
            dft_manifest_path=manifest,
            reference_registry_path=reference_registry,
            j_gate_path=self._approved_j_gate(),
            iodine_gate_path=self.project / "jobs" / "unused_iodine_gate.json",
            resource_policy_path=self.resource_policy,
            output_root=self.project / f"03_nmr_properties_{tag}",
            shielding_manifest_path=self.project / "data" / f"shielding_{tag}.csv",
            j_manifest_path=self.project / "data" / f"j_{tag}.csv",
            weights_path=self.project / "data" / f"weights_{tag}.csv",
        )

    def test_draft_is_never_approved_and_builder_rejects_it(self) -> None:
        draft, _approved = self._draft_and_approve()
        with self.assertRaises(WorkflowValidationError):
            build_dft_jobs(
                self.project,
                [draft],
                self.manifest,
                self.registry,
                self.resource_policy,
                self.project / "bad_output",
                self.project / "data" / "bad_manifest.csv",
            )

    def test_approved_gate_builds_scenario_isolated_explicit_jobs(self) -> None:
        _draft, approved = self._draft_and_approve()
        output_root = self.project / "02_dft_opt_freq_explicit"
        rows = build_dft_jobs(
            self.project,
            [approved],
            self.manifest,
            self.registry,
            self.resource_policy,
            output_root,
            self.project / "data" / "DFT_explicit_job_manifest.csv",
        )
        self.assertGreaterEqual(len(rows), 2)
        self.assertEqual({row["scenario_id"] for row in rows}, {"D2O_01_primary", "D2O_02_primary"})
        for row in rows:
            job_dir = Path(str(row["job_directory"]))
            inp = (job_dir / "dft_opt_freq.inp").read_text(encoding="utf-8")
            self.assertIn(
                "! PBE0 D4 def2-TZVPP def2/J RIJCOSX CPCM(WATER) VeryTightOpt Freq VeryTightSCF DefGrid3",
                inp,
            )
            self.assertIn("%pal nprocs 6 end", inp)
            self.assertIn("Temp 300, 310", inp)
            self.assertEqual(
                build_input_bundle(job_dir / "dft_opt_freq.inp")["sha256"],
                row["expected_input_bundle_sha256"],
            )
            with (job_dir / "source_isotope_map.csv").open("r", encoding="utf-8-sig") as handle:
                isotope_count = len(list(csv.DictReader(handle)))
            self.assertEqual(inp.count(f"M={DEUTERIUM_MASS_U:.11f}"), isotope_count)
            provenance = json.loads((job_dir / "provenance.json").read_text(encoding="utf-8"))
            self.assertEqual(provenance["condition_isolation_key"], ["F3_001", row["scenario_id"]])
            self.assertEqual(provenance["checks"]["explicit_cluster_atom_order_preserved"], "PASS")
            self.assertTrue(provenance["downstream_real_nmr_calibration_reference"])
            source_selected = Path(provenance["selected_coordinate_source_file"]).read_bytes()
            self.assertEqual((job_dir / "approved_selected_cluster.xyz").read_bytes(), source_selected)
            with (job_dir / "nmr_atom_map_v2.csv").open("r", encoding="utf-8-sig") as handle:
                map_rows = list(csv.DictReader(handle))
            self.assertTrue(any(item["nmr_scope"] == "ENVIRONMENT_EXCLUDED_FROM_ANALYTE_NMR" for item in map_rows))
            self.assertEqual(len(map_rows), provenance["cluster_atom_count"])

    def test_opt_runner_lock_and_provenance_tamper_gate(self) -> None:
        _draft, approved = self._draft_and_approve()
        manifest = self.project / "data" / "DFT_explicit_job_manifest.csv"
        build_dft_jobs(
            self.project,
            [approved],
            self.manifest,
            self.registry,
            self.resource_policy,
            self.project / "02_dft_opt_freq_explicit",
            manifest,
        )
        lock = self.project / "logs" / "opt_freq_serial" / "RUNNER.lock"
        with ExclusiveLock(lock):
            with self.assertRaises(WorkflowValidationError):
                with ExclusiveLock(lock):
                    pass
        payload = json.loads(approved.read_text(encoding="utf-8"))
        payload["decision_notes"].append("intentional unit-test tamper")
        write(approved, json.dumps(payload, indent=2) + "\n")
        with self.assertRaises(WorkflowValidationError):
            dry_run_validate(self.project, manifest)

    def test_minimum_weights_and_property_requests_remain_scenario_and_analyte_only(self) -> None:
        _draft, approved, _rows, manifest = self._draft_and_approve_minimum()
        references = self._approved_reference_registry()
        shielding, coupling, weights = self._property_build(
            approved, manifest, references, "valid"
        )
        self.assertTrue(shielding)
        self.assertTrue(coupling)
        self.assertEqual(
            {row["scenario_id"] for row in shielding},
            {"D2O_01_primary", "D2O_02_primary"},
        )
        grouped: dict[tuple[str, str], float] = {}
        for row in weights:
            key = (str(row["scenario_id"]), str(row["nucleus"]))
            grouped[key] = grouped.get(key, 0.0) + float(row["weight"])
        self.assertEqual(
            set(grouped),
            {
                ("D2O_01_primary", "1H"),
                ("D2O_01_primary", "13C"),
                ("D2O_02_primary", "1H"),
                ("D2O_02_primary", "13C"),
            },
        )
        self.assertTrue(all(abs(total - 1.0) < 1e-10 for total in grouped.values()))
        self.assertEqual(
            {row["nucleus"]: float(row["temperature_K"]) for row in weights},
            {"1H": 310.0, "13C": 300.0},
        )
        for row in shielding:
            directory = Path(str(row["job_directory"]))
            with (directory / "nmr_atom_map_v2.csv").open(
                "r", encoding="utf-8-sig"
            ) as handle:
                atom_map = list(csv.DictReader(handle))
            environment = {
                int(item["atom_index_1based"])
                for item in atom_map
                if item["nmr_scope"] != "ANALYTE"
            }
            with (directory / "shielding_request_map.csv").open(
                "r", encoding="utf-8-sig"
            ) as handle:
                requests = list(csv.DictReader(handle))
            requested = {
                int(value)
                for item in requests
                for value in item["atom_indices_1based_nuclei"].split(",")
            }
            self.assertTrue(requested)
            self.assertTrue(requested.isdisjoint(environment))
            self.assertIn("CPCM(WATER)", (directory / "shielding.inp").read_text())
            self.assertEqual(
                build_input_bundle(directory / "shielding.inp")["sha256"],
                row["expected_input_bundle_sha256"],
            )
        for row in coupling:
            directory = Path(str(row["job_directory"]))
            j_input = (directory / "j_coupling.inp").read_text(encoding="utf-8")
            self.assertIn("%method Z_Tol 1e-12 end", j_input)
            self.assertIn("Tol = 1e-10", j_input)
            self.assertIn("SpinSpinRThresh 6.05", j_input)
            with (directory / "j_pair_allowlist.csv").open(
                "r", encoding="utf-8-sig"
            ) as handle:
                pairs = list(csv.DictReader(handle))
            self.assertTrue(pairs)
            self.assertTrue(
                all(
                    int(item["left_atom_index_1based_nuclei"]) <= 6
                    and int(item["right_atom_index_1based_nuclei"]) <= 6
                    and int(item["left_atom_index_0based_spinspin"]) <= 5
                    and int(item["right_atom_index_0based_spinspin"]) <= 5
                    for item in pairs
                )
            )

    def test_property_builder_rejects_cross_scenario_decision(self) -> None:
        _draft, approved, _rows, manifest = self._draft_and_approve_minimum()
        payload = json.loads(approved.read_text(encoding="utf-8"))
        original = payload["conformer_decisions"][0]["condition_scenario_id"]
        payload["conformer_decisions"][0]["condition_scenario_id"] = (
            "D2O_02_primary" if original == "D2O_01_primary" else "D2O_01_primary"
        )
        crossed = approved.with_name("minimum_boltzmann_crossed.json")
        write(crossed, json.dumps(payload, indent=2) + "\n")
        with self.assertRaises(WorkflowValidationError):
            self._property_build(
                crossed,
                manifest,
                self._approved_reference_registry(),
                "crossed",
            )

    def test_property_builder_rejects_reference_method_mismatch(self) -> None:
        _draft, approved, _rows, manifest = self._draft_and_approve_minimum()
        with self.assertRaises(WorkflowValidationError):
            self._property_build(
                approved,
                manifest,
                self._approved_reference_registry(mismatched=True),
                "reference_mismatch",
            )

    def test_property_builder_rejects_draft_reference_registry(self) -> None:
        _draft, approved, _rows, manifest = self._draft_and_approve_minimum()
        registry = self._approved_reference_registry()
        payload = json.loads(registry.read_text(encoding="utf-8"))
        payload["status"] = "DRAFT_NOT_APPROVED"
        draft_registry = registry.with_name("reference_registry_draft.json")
        write(draft_registry, json.dumps(payload, indent=2) + "\n")
        with self.assertRaises(WorkflowValidationError):
            self._property_build(
                approved,
                manifest,
                draft_registry,
                "draft_reference",
            )

    def test_draft_j_gate_blocks_only_j_jobs(self) -> None:
        _draft, approved, _rows, manifest = self._draft_and_approve_minimum()
        j_gate = self._approved_j_gate()
        payload = json.loads(j_gate.read_text(encoding="utf-8"))
        payload["approval_status"] = "DRAFT_NOT_APPROVED"
        draft_j = j_gate.with_name("j_method_benchmark_DRAFT.json")
        write(draft_j, json.dumps(payload, indent=2) + "\n")
        shielding, coupling, _weights = build_property_jobs(
            project=self.project,
            gate_paths=[approved],
            dft_manifest_path=manifest,
            reference_registry_path=self._approved_reference_registry(),
            j_gate_path=draft_j,
            iodine_gate_path=self.project / "jobs" / "unused_iodine_gate.json",
            resource_policy_path=self.resource_policy,
            output_root=self.project / "03_nmr_properties_draft_j",
            shielding_manifest_path=self.project / "data" / "shielding_draft_j.csv",
            j_manifest_path=self.project / "data" / "j_draft_j.csv",
            weights_path=self.project / "data" / "weights_draft_j.csv",
        )
        self.assertTrue(shielding)
        self.assertEqual(coupling, [])

    def test_read_only_orca_deuterium_smoke_evidence(self) -> None:
        input_path = PROJECT / "scratch" / "syntax_test_deuterium.inp"
        output_path = PROJECT / "scratch" / "syntax_test_deuterium.out"
        if not input_path.is_file() or not output_path.is_file():
            self.skipTest("Parent-run ORCA 6.1.1 deuterium smoke evidence is absent")
        input_text = input_path.read_text(encoding="utf-8", errors="replace")
        output_bytes = output_path.read_bytes()
        output_text = output_bytes.decode(
            "utf-16" if output_bytes.startswith((b"\xff\xfe", b"\xfe\xff")) else "utf-8",
            errors="replace",
        )
        self.assertGreaterEqual(input_text.count(f"M={DEUTERIUM_MASS_U:.11f}"), 1)
        self.assertIn("Total Mass", output_text)
        self.assertIn("20.03 AMU", output_text)
        self.assertIn("ORCA TERMINATED NORMALLY", output_text)


if __name__ == "__main__":
    unittest.main()
