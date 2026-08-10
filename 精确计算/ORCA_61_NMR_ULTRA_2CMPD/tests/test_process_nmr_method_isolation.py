from __future__ import annotations

import csv
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from process_nmr_predictions import (  # noqa: E402
    DEFAULT_SHIELDING_METHOD,
    ProcessingError,
    group_manifest_by_state_and_conformer,
    load_locked_nucleus_conditions,
    load_nucleus_weight_sets,
    load_scenario_atom_map,
    parse_atom_map,
    parser,
    property_coordinate_path,
    resolve_same_method_tms_output,
    select_state_keys,
    select_shielding_method_rows,
    validate_scenario_isotope_maps,
)


PRODUCTION = "production_revdsd_pcsseg3"
CROSSCHECK = "crosscheck_pbe0_pcsseg3"
PRODUCTION_FP = (
    "revDSD-PBEP86-D4/2021|pcSseg-3|AutoAux|RIJCOSX|GIAO|NoFrozenCore|"
    "CPCM(DMSO)|ExtremeSCF|DefGrid3|ORCA_6.1.1"
)
CROSSCHECK_FP = (
    "PBE0|pcSseg-3|AutoAux|RIJCOSX|GIAO|CPCM(DMSO)|ExtremeSCF|"
    "DefGrid3|ORCA_6.1.1"
)


def shielding_row(conformer: str, branch: str, fingerprint: str) -> dict[str, str]:
    return {
        "job_id": f"{conformer}__{branch}",
        "figure_id": "F3_005",
        "conformer_id": conformer,
        "property_branch": branch,
        "method_fingerprint": fingerprint,
    }


class ShieldingMethodIsolationTests(unittest.TestCase):
    def test_cli_defaults_to_production_and_resolves_tms_later(self) -> None:
        args = parser().parse_args([])
        self.assertEqual(args.shielding_method, DEFAULT_SHIELDING_METHOD)
        self.assertIsNone(args.tms_output)
        self.assertEqual(args.dft_manifest.name, "DFT_explicit_job_manifest.csv")
        self.assertEqual(args.weights.name, "NMR_conformer_weights.csv")
        self.assertEqual(args.coupling_manifest.name, "NMR_j_coupling_job_manifest.csv")

    def test_two_branches_are_filtered_before_conformer_grouping(self) -> None:
        rows = [
            shielding_row("conf01", PRODUCTION, PRODUCTION_FP),
            shielding_row("conf01", CROSSCHECK, CROSSCHECK_FP),
            shielding_row("conf02", PRODUCTION, PRODUCTION_FP),
            shielding_row("conf02", CROSSCHECK, CROSSCHECK_FP),
        ]
        selected, fingerprint = select_shielding_method_rows(rows, PRODUCTION)
        self.assertEqual([row["conformer_id"] for row in selected], ["conf01", "conf02"])
        self.assertTrue(all(row["property_branch"] == PRODUCTION for row in selected))
        self.assertEqual(fingerprint, PRODUCTION_FP)

    def test_inconsistent_fingerprints_within_branch_fail_closed(self) -> None:
        rows = [
            shielding_row("conf01", PRODUCTION, PRODUCTION_FP),
            shielding_row("conf02", PRODUCTION, CROSSCHECK_FP),
        ]
        with self.assertRaisesRegex(ProcessingError, "multiple fingerprints"):
            select_shielding_method_rows(rows, PRODUCTION)

    def test_missing_branch_columns_fail_closed(self) -> None:
        with self.assertRaisesRegex(ProcessingError, "missing columns"):
            select_shielding_method_rows(
                [{"job_id": "legacy", "figure_id": "F3_005", "conformer_id": "c1"}],
                PRODUCTION,
            )

    def test_tms_default_is_selected_by_exact_method_fingerprint(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference_directory = root / "reference_jobs" / "TMS" / PRODUCTION
            manifest = root / "reference_job_manifest.csv"
            fields = [
                "job_id",
                "stage",
                "method_id",
                "method_fingerprint",
                "job_directory",
                "output_file",
            ]
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerow(
                    {
                        "job_id": "TMS_DMSO__production",
                        "stage": "reference_shielding",
                        "method_id": PRODUCTION,
                        "method_fingerprint": PRODUCTION_FP,
                        "job_directory": str(reference_directory),
                        "output_file": "tms_shield.out",
                    }
                )
            output, row = resolve_same_method_tms_output(
                manifest, PRODUCTION, PRODUCTION_FP
            )
            self.assertEqual(output, (reference_directory / "tms_shield.out").resolve())
            self.assertEqual(row["job_id"], "TMS_DMSO__production")

    def test_tms_fingerprint_mismatch_fails_even_with_explicit_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "reference_job_manifest.csv"
            fields = [
                "job_id",
                "stage",
                "method_id",
                "method_fingerprint",
                "job_directory",
                "output_file",
            ]
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerow(
                    {
                        "job_id": "TMS_DMSO__production",
                        "stage": "reference_shielding",
                        "method_id": PRODUCTION,
                        "method_fingerprint": CROSSCHECK_FP,
                        "job_directory": str(root),
                        "output_file": "tms_shield.out",
                    }
                )
            with self.assertRaisesRegex(ProcessingError, "fingerprint mismatch"):
                resolve_same_method_tms_output(
                    manifest,
                    PRODUCTION,
                    PRODUCTION_FP,
                    root / "relocated_tms.out",
                )


class ScenarioIsolationAndPortableMapTests(unittest.TestCase):
    @staticmethod
    def _write_map(path: Path, final_element: str = "O") -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        fields = [
            "atom_index_1based",
            "element",
            "component",
            "environment_atom",
            "exchangeable",
            "include_in_analyte_nmr",
        ]
        rows = [
            (1, "N", "analyte", "false", "false", "true"),
            (2, "H", "analyte", "false", "true", "true"),
            (3, "C", "analyte", "false", "false", "true"),
            (4, final_element, "dmso_1", "true", "false", "false"),
        ]
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for index, element, component, environment, exchangeable, include in rows:
                writer.writerow(
                    {
                        "atom_index_1based": index,
                        "element": element,
                        "component": component,
                        "environment_atom": environment,
                        "exchangeable": exchangeable,
                        "include_in_analyte_nmr": include,
                    }
                )

    @staticmethod
    def _orca_output(path: Path, gibbs: float) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            "\n".join(
                (
                    "Quasi RRHO ... True",
                    "THERMOCHEMISTRY AT 298.15 K",
                    f"Final Gibbs free energy ... {gibbs:.12f} Eh",
                    f"FINAL SINGLE POINT ENERGY {gibbs:.12f}",
                    "VIBRATIONAL FREQUENCIES",
                    "  0:       100.0000 cm**-1",
                    "NORMAL MODES",
                    "****ORCA TERMINATED NORMALLY****",
                    "",
                )
            ),
            encoding="utf-8",
        )

    def test_same_conformer_token_in_two_scenarios_is_never_grouped_together(self) -> None:
        rows = [
            {
                "figure_id": "F3_005",
                "scenario_id": "DMSO_00",
                "conformer_id": "conf01",
            },
            {
                "figure_id": "F3_005",
                "scenario_id": "DMSO_01",
                "conformer_id": "conf01",
            },
        ]
        grouped = group_manifest_by_state_and_conformer(rows)
        self.assertEqual(set(grouped), {("F3_005", "DMSO_00"), ("F3_005", "DMSO_01")})
        self.assertEqual(len(grouped[("F3_005", "DMSO_00")]), 1)
        self.assertEqual(len(grouped[("F3_005", "DMSO_01")]), 1)

    def test_scenario_selector_keeps_states_explicit(self) -> None:
        available = {
            ("F3_005", "DMSO_00"),
            ("F3_005", "DMSO_01"),
            ("F3_015", "DMSO_00"),
        }
        self.assertEqual(
            select_state_keys(available, ["F3_005"], ["DMSO_01"]),
            {("F3_005", "DMSO_01")},
        )

    def test_portable_map_excludes_environment_and_accepts_no_adjacency(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "nmr_atom_map_v2.csv"
            self._write_map(path)
            atom_map = parse_atom_map(path)
            self.assertEqual(len(atom_map), 4)
            self.assertEqual(atom_map[1]["exchangeable_H_flag"], 1)
            self.assertEqual(atom_map[1]["exchangeability_structural_confidence"], "medium")
            self.assertEqual(atom_map[3]["environment_atom"], 1)
            self.assertEqual(atom_map[3]["include_in_analyte_nmr"], 0)
            self.assertNotEqual(atom_map[3]["nmr_scope"], "ANALYTE")

    def test_enriched_legacy_orca_index_and_adjacency_schema_remains_supported(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "nmr_atom_map_v2.csv"
            fields = [
                "orca_index_0based",
                "orca_index_1based",
                "element",
                "formal_charge",
                "canonical_symmetry_rank",
                "attached_atom_indices_0based",
                "attached_elements",
                "exchangeable_H_flag",
            ]
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerow(
                    {
                        "orca_index_0based": 0,
                        "orca_index_1based": 1,
                        "element": "N",
                        "formal_charge": 0,
                        "canonical_symmetry_rank": 1,
                        "attached_atom_indices_0based": "1",
                        "attached_elements": "H",
                        "exchangeable_H_flag": 0,
                    }
                )
                writer.writerow(
                    {
                        "orca_index_0based": 1,
                        "orca_index_1based": 2,
                        "element": "H",
                        "formal_charge": 0,
                        "canonical_symmetry_rank": 2,
                        "attached_atom_indices_0based": "0",
                        "attached_elements": "N",
                        "exchangeable_H_flag": 1,
                    }
                )
            atom_map = parse_atom_map(path)
            self.assertEqual(atom_map[1]["exchangeable_H_parent_index_0based"], 0)
            self.assertEqual(atom_map[1]["exchangeable_H_class"], "N-H")
            self.assertEqual(atom_map[1]["exchangeability_structural_confidence"], "high")

    def test_all_job_local_maps_must_have_one_hash_per_scenario(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            directories = [root / name for name in ("s1", "s2", "j1", "j2")]
            for directory in directories:
                self._write_map(directory / "nmr_atom_map_v2.csv")
            shielding = {
                "c1": {"job_directory": str(directories[0])},
                "c2": {"job_directory": str(directories[1])},
            }
            coupling = {
                "c1": {"job_directory": str(directories[2])},
                "c2": {"job_directory": str(directories[3])},
            }
            _path, atom_map, digest, checked = load_scenario_atom_map(
                "F3_005", "DMSO_00", shielding, coupling
            )
            self.assertEqual(len(atom_map), 4)
            self.assertEqual(len(digest), 64)
            self.assertEqual(len(checked), 4)
            self._write_map(directories[3] / "nmr_atom_map_v2.csv", final_element="S")
            with self.assertRaisesRegex(ProcessingError, "hash mismatch"):
                load_scenario_atom_map("F3_005", "DMSO_00", shielding, coupling)

    def test_optimized_coordinate_is_preferred_with_legacy_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            optimized = directory / "optimized.xyz"
            legacy = directory / "conformer.xyz"
            legacy.write_text("0\nlegacy\n", encoding="utf-8")
            self.assertEqual(
                property_coordinate_path({"job_directory": str(directory)}), legacy.resolve()
            )
            optimized.write_text("0\noptimized\n", encoding="utf-8")
            self.assertEqual(
                property_coordinate_path({"job_directory": str(directory)}), optimized.resolve()
            )

    def test_explicit_dmso_d6_mass_map_survives_xyz_property_handoff(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            directories = [root / name for name in ("s1", "s2", "j1", "j2")]
            payload = "atom_index_1based,element,isotope,mass_u,component\n4,H,2,2.01410177812,explicit_DMSO_d6_01\n"
            rows = []
            for ordinal, directory in enumerate(directories):
                directory.mkdir()
                isotope_path = directory / "source_isotope_map.csv"
                isotope_path.write_text(payload, encoding="utf-8")
                digest = hashlib.sha256(isotope_path.read_bytes()).hexdigest()
                rows.append(
                    {
                        "job_id": f"property_{ordinal}",
                        "job_directory": str(directory),
                        "source_isotope_map_sha256": digest,
                    }
                )
            atom_map = {
                0: {"element": "C", "component": "analyte", "environment_atom": 0},
                1: {"element": "H", "component": "analyte", "environment_atom": 0},
                2: {"element": "S", "component": "explicit_DMSO_d6_01", "environment_atom": 1},
                3: {"element": "H", "component": "explicit_DMSO_d6_01", "environment_atom": 1},
            }
            result = validate_scenario_isotope_maps(
                "F3_005",
                "DMSO_01",
                atom_map,
                {"c1": rows[0], "c2": rows[1]},
                {"c1": rows[2], "c2": rows[3]},
            )
            self.assertEqual(result["deuterium_atom_indices_1based"], [4])
            self.assertEqual(
                result["status"], "PASS_EXPLICIT_DMSO_D6_MASS_MAPPING_PRESERVED"
            )

    def test_nucleus_specific_weights_and_locked_temperatures_are_authoritative(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            conditions = root / "nucleus_conditions.csv"
            with conditions.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=(
                        "figure_id",
                        "nucleus",
                        "temperature_K",
                        "condition_status",
                        "source_note",
                    ),
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "figure_id": "F3_005",
                        "nucleus": "1H",
                        "temperature_K": "312.609",
                        "condition_status": "VERIFIED_FROM_LOCAL_BRUKER_ACQUS",
                        "source_note": "Sample 1-1-H",
                    }
                )
                writer.writerow(
                    {
                        "figure_id": "F3_005",
                        "nucleus": "13C",
                        "temperature_K": "304.0651",
                        "condition_status": "VERIFIED_FROM_LOCAL_BRUKER_ACQUS",
                        "source_note": "Sample 1-1-C",
                    }
                )
            rows = []
            for nucleus, temperature, weights in (
                ("1H", "312.609", ("0.75", "0.25")),
                ("13C", "304.0651", ("0.20", "0.80")),
            ):
                for conformer, energy, weight in zip(
                    ("conf01", "conf02"), ("-100.0", "-99.999"), weights
                ):
                    rows.append(
                        {
                            "figure_id": "F3_005",
                            "scenario_id": "DMSO_01",
                            "nucleus": nucleus,
                            "temperature_K": temperature,
                            "conformer_id": conformer,
                            "qrrho_gibbs_free_energy_hartree": energy,
                            "weight": weight,
                            "approved_gate_sha256": "A" * 64,
                        }
                    )
            loaded = load_nucleus_weight_sets(
                rows,
                {("F3_005", "DMSO_01")},
                load_locked_nucleus_conditions(conditions),
            )[("F3_005", "DMSO_01")]
            self.assertEqual(loaded["1H"]["temperature_K"], 312.609)
            self.assertEqual(loaded["13C"]["temperature_K"], 304.0651)
            self.assertEqual(loaded["1H"]["conformers"]["conf01"]["weight"], 0.75)
            self.assertEqual(loaded["13C"]["conformers"]["conf01"]["weight"], 0.20)

            rows[0]["temperature_K"] = "298.15"
            with self.assertRaisesRegex(ProcessingError, "mixes temperatures|does not match"):
                load_nucleus_weight_sets(
                    rows,
                    {("F3_005", "DMSO_01")},
                    load_locked_nucleus_conditions(conditions),
                )



if __name__ == "__main__":
    unittest.main()
