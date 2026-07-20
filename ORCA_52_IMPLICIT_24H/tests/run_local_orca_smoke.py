#!/usr/bin/env python3
"""Optional real-ORCA smoke test (xTB opt + fast NMR shielding)."""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from common import (  # noqa: E402
    build_opt_input,
    build_shield_input,
    find_final_geometry,
    parse_shielding_summary,
    read_csv,
    read_xyz,
    run_orca_job,
    read_json,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orca", type=Path, required=True)
    parser.add_argument("--figure-id", default="F3_001")
    parser.add_argument("--timeout", type=float, default=900.0)
    args = parser.parse_args()
    species = {row["figure_id"]: row for row in read_csv(ROOT / "config" / "primary_species.csv")}
    if args.figure_id not in species:
        raise ValueError(f"Unknown Figure ID: {args.figure_id}")
    row = species[args.figure_id]
    start = ROOT / "structures" / args.figure_id / "start.xyz"
    methods = read_json(ROOT / "config" / "methods_24h.json")
    with tempfile.TemporaryDirectory(prefix="orca52_smoke_") as temp:
        temp_root = Path(temp)
        raw = temp_root / "raw"
        scratch = temp_root / "scratch"
        xtb_input = build_opt_input(
            f"XTB Opt ALPB({row['orca_continuum_keyword']})",
            int(row["charge"]), int(row["multiplicity"]), start.name, 2, 1800
        )
        xtb, xtb_dir = run_orca_job(
            job_id=f"{args.figure_id}_xtb", input_text=xtb_input, geometry_path=start,
            orca_path=args.orca, scratch_root=scratch, raw_root=raw,
            cpu_ids=[0, 1], timeout_seconds=args.timeout, basename="xtb_opt",
        )
        if xtb.status != "COMPLETE":
            raise RuntimeError(f"xTB smoke step failed: {xtb}")
        geometry = find_final_geometry(xtb_dir, "xtb_opt")
        if geometry is None:
            raise RuntimeError("xTB smoke step produced no geometry")
        elements = [atom[0] for atom in read_xyz(geometry)[1]]
        iodine = methods["iodine_override"]["basis_block"] if args.figure_id == "F3_006" else None
        shield_input = build_shield_input(
            f"PBE0 pcSseg-1 AutoAux RIJCOSX TightSCF DefGrid3 CPCM({row['orca_continuum_keyword']})",
            int(row["charge"]), int(row["multiplicity"]), geometry.name, 4, 3000, elements,
            iodine_override=iodine,
        )
        shield, shield_dir = run_orca_job(
            job_id=f"{args.figure_id}_shield", input_text=shield_input, geometry_path=geometry,
            orca_path=args.orca, scratch_root=scratch, raw_root=raw,
            cpu_ids=[0, 1, 2, 3], timeout_seconds=args.timeout, basename="shield",
        )
        if shield.status != "COMPLETE":
            raise RuntimeError(f"shielding smoke step failed: {shield}")
        values = parse_shielding_summary(shield_dir / "shield.out")
        expected = sum(element in {"H", "C"} for element in elements)
        if len(values) != expected:
            raise RuntimeError(f"shielding count {len(values)} != expected {expected}")
        print(f"PASS real ORCA smoke test {args.figure_id}: {len(values)} H/C shieldings parsed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
