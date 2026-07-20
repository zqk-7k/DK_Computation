#!/usr/bin/env python3
"""Offline QA: registry, input syntax, HRMS engine, post-processing and audit."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from common import (  # noqa: E402
    build_shield_input,
    read_csv,
    read_json,
    read_xyz,
    runtime_environment_identity,
    sha256,
    write_csv,
    write_json,
)
from postprocess_nmr import first_order_lines, lorentzian_area  # noqa: E402


def assert_registry() -> None:
    species = read_csv(ROOT / "config" / "primary_species.csv")
    assert len(species) == 52
    assert [row["figure_id"] for row in species] == [f"F3_{i:03d}" for i in range(1, 53)]
    by_id = {row["figure_id"]: row for row in species}
    assert by_id["F3_016"]["charge"] == "1"
    assert by_id["F3_017"]["charge"] == "1"
    assert by_id["F3_021"]["charge"] == "2"
    assert by_id["F3_035"]["observable_13C_count"] == "0"
    assert by_id["F3_043"]["observable_13C_count"] == "0"
    assert by_id["F3_044"]["observable_13C_count"] == "0"
    for row in species:
        figure_id = row["figure_id"]
        atom_map = read_csv(ROOT / "structures" / figure_id / "atom_map.csv")
        _, atoms = read_xyz(ROOT / "structures" / figure_id / "start.xyz")
        assert len(atom_map) == len(atoms)
        assert all(atom["environment_atom"] == "false" for atom in atom_map)


def assert_input_order() -> None:
    text = build_shield_input(
        "PBE0 pcSseg-1 AutoAux RIJCOSX TightSCF CPCM(DMSO)",
        0, 1, "geometry.xyz", 4, 4800, ["C", "O", "N", "H"],
    )
    assert text.index("* xyzfile") < text.index("%eprnmr")
    assert "Nuclei = ALL H { SHIFT }" in text
    assert "Nuclei = ALL C { SHIFT }" in text
    assert "mpirun" not in text.lower()


def assert_profile_numerics() -> None:
    couplings = [
        {
            "equivalent_partner_count": 3,
            "J_Hz": value,
        }
        for value in (12.0, 10.0, 8.0, 6.0, 4.0, 2.0)
    ]
    lines = first_order_lines(5.0, 2.0, couplings, 400.132470802)
    assert len(lines) <= 512
    assert abs(sum(area for _, area in lines) - 2.0) < 1.0e-10
    centroid = sum(ppm * area for ppm, area in lines) / sum(area for _, area in lines)
    assert abs(centroid - 5.0) < 1.0e-10

    carbon_axis = np.linspace(220.0, -10.0, 32768)
    profile = lorentzian_area(carbon_axis, 123.456, 1.0 / 100.622829803, 1.0)
    captured = abs(float(np.trapezoid(profile, carbon_axis)))
    assert 0.95 <= captured <= 1.05


def shielding_rows(
    figure_id: str,
    method_id: str,
    atom_map_path: Path,
    raw_output: Path,
    run_dir: Path,
    reference: bool = False,
) -> list[dict[str, object]]:
    retained_input = next(raw_output.parent.glob("*.inp"))
    retained_geometry = next(raw_output.parent.glob("*.xyz"))
    rows = []
    for atom in read_csv(atom_map_path):
        if atom["element"] not in {"H", "C"}:
            continue
        index = int(atom["atom_index_1based"])
        if atom["element"] == "H":
            shielding = 31.0 if reference else 24.0 + 0.12 * index
        else:
            shielding = 185.0 if reference else 150.0 - 1.5 * index
        rows.append(
            {
                "figure_id": figure_id,
                "atom_index_1based": index,
                "element": atom["element"],
                "isotropic_shielding_ppm": shielding,
                "anisotropy_ppm": 0.0,
                "include_in_analyte_nmr": atom.get("include_in_analyte_nmr", "true"),
                "exchangeable": atom.get("exchangeable", "false"),
                "reference_signal": atom.get("reference_signal", "false"),
                "method_id": method_id,
                "method_quality_rank": 20,
                "method_quality_label": "TIER1_DFT_SHIFT",
                "geometry_tier": "BASELINE_GEOMETRY",
                "implicit_solvent": "DMSO",
                "source_output": str(raw_output),
                "source_output_relative_to_run_dir": raw_output.relative_to(run_dir).as_posix(),
                "source_output_sha256": sha256(raw_output),
                "source_input_sha256": sha256(retained_input),
                "source_geometry_sha256": sha256(retained_geometry),
                "workflow_fingerprint": "OFFLINE_TEST_FINGERPRINT",
                "run_uuid": "OFFLINE_TEST_UUID",
            }
        )
    return rows


def assert_postprocess_and_audit() -> None:
    method = "TIER1_PBE0_PCSSEG1"
    with tempfile.TemporaryDirectory(prefix="orca52_offline_") as temporary:
        run_dir = Path(temporary)
        write_json(
            run_dir / "run_config.json",
            {
                "run_uuid": "OFFLINE_TEST_UUID",
                "workflow_fingerprint": "OFFLINE_TEST_FINGERPRINT",
                "runtime_environment_identity": runtime_environment_identity(),
            },
        )
        def retained_attempt(job: str, geometry_source: Path, text: str) -> Path:
            folder = run_dir / "raw_orca" / job / "attempt_001"
            folder.mkdir(parents=True, exist_ok=True)
            raw = folder / "shield.out"
            raw.write_text(text, encoding="utf-8")
            retained_input = folder / "shield.inp"
            retained_input.write_text("! OFFLINE PROVENANCE FIXTURE\n", encoding="utf-8")
            retained_geometry = folder / "geometry.xyz"
            shutil.copy2(geometry_source, retained_geometry)
            write_json(
                folder / "run_record.json",
                {
                    "status": "COMPLETE",
                    "input_sha256": sha256(retained_input),
                    "geometry_sha256": sha256(retained_geometry),
                    "output_file": str(raw),
                },
            )
            return raw

        sample_raw = retained_attempt(
            "F3_001_shield_TIER1_PBE0_PCSSEG1",
            ROOT / "structures" / "F3_001" / "start.xyz",
            "OFFLINE TEST SAMPLE RAW OUTPUT\n",
        )
        reference_raw = retained_attempt(
            "REF_TMS_DMSO_shield_TIER1_PBE0_PCSSEG1",
            ROOT / "references" / "TMS_DMSO" / "start.xyz",
            "OFFLINE TEST REFERENCE RAW OUTPUT\n",
        )
        sample_path = run_dir / "parsed_shieldings" / "compounds" / "F3_001" / f"{method}.csv"
        sample_rows = shielding_rows(
            "F3_001", method, ROOT / "structures" / "F3_001" / "atom_map.csv",
            sample_raw, run_dir,
        )
        write_csv(sample_path, sample_rows, list(sample_rows[0].keys()))
        reference_path = run_dir / "parsed_shieldings" / "references" / "TMS_DMSO" / f"{method}.csv"
        reference_rows = shielding_rows(
            "TMS_DMSO", method, ROOT / "references" / "TMS_DMSO" / "atom_map.csv",
            reference_raw, run_dir, reference=True,
        )
        write_csv(reference_path, reference_rows, list(reference_rows[0].keys()))
        geometry = run_dir / "geometries" / "baseline" / "F3_001.xyz"
        geometry.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / "structures" / "F3_001" / "start.xyz", geometry)

        result = subprocess.run(
            [
                sys.executable, str(ROOT / "scripts" / "postprocess_nmr.py"),
                "--run-dir", str(run_dir), "--only", "F3_001,F3_035",
            ],
            text=True, capture_output=True,
        )
        if result.returncode != 0:
            raise AssertionError(f"postprocess failed\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")

        hrms = subprocess.run(
            [
                sys.executable, str(ROOT / "scripts" / "generate_hrms_all.py"),
                "--output-dir", str(run_dir / "HRMS_CALCULATED_NOT_EXPERIMENTAL"),
                "--allowed-root", str(run_dir),
                "--only", "F3_001,F3_035",
                "--overwrite",
            ],
            text=True,
            capture_output=True,
        )
        if hrms.returncode != 0:
            raise AssertionError(f"HRMS fixture failed\nSTDOUT:\n{hrms.stdout}\nSTDERR:\n{hrms.stderr}")
        audit = subprocess.run(
            [
                sys.executable, str(ROOT / "scripts" / "audit_outputs.py"),
                "--run-dir", str(run_dir), "--only", "F3_001,F3_035",
            ],
            text=True, capture_output=True,
        )
        if audit.returncode != 0:
            raise AssertionError(f"audit failed\nSTDOUT:\n{audit.stdout}\nSTDERR:\n{audit.stderr}")
        final = read_json(run_dir / "FINAL_AUDIT.json")
        assert final["status"] == "PASS"


def main() -> int:
    assert_registry()
    print("PASS registry/species/implicit-only checks")
    assert_input_order()
    print("PASS ORCA NMR input ordering checks")
    assert_profile_numerics()
    print("PASS multiplet area/centroid and digital-profile checks")
    assert_postprocess_and_audit()
    print("PASS end-to-end offline postprocess/audit fixture")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
