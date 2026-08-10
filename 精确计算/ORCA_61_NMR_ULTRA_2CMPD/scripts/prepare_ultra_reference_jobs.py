"""Prepare same-level TMS reference jobs for every shielding branch."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "reference_jobs" / "TMS_DMSO"
MANIFEST = PROJECT / "data" / "reference_job_manifest.csv"


METHODS = (
    (
        "production_revdsd_pcsseg3",
        "revDSD-PBEP86-D4/2021 pcSseg-3 AutoAux NMR NoFrozenCore CPCM(DMSO) RIJCOSX DefGrid3 ExtremeSCF",
        "%mp2\n  Density Relaxed\nend",
        "shielding",
        "revDSD-PBEP86-D4/2021|pcSseg-3|AutoAux|RIJCOSX|GIAO|NoFrozenCore|CPCM(DMSO)|ExtremeSCF|DefGrid3|ORCA_6.1.1",
    ),
    (
        "crosscheck_pbe0_pcsseg3",
        "PBE0 pcSseg-3 AutoAux NMR CPCM(DMSO) RIJCOSX DefGrid3 ExtremeSCF",
        "",
        "shielding",
        "PBE0|pcSseg-3|AutoAux|RIJCOSX|GIAO|CPCM(DMSO)|ExtremeSCF|DefGrid3|ORCA_6.1.1",
    ),
    (
        "benchmark_ri_mp2_pcsseg3",
        "RIJK RI-MP2 pcSseg-3 AutoAux NMR NoFrozenCore CPCM(DMSO) DefGrid3 ExtremeSCF",
        "%mp2\n  Density Relaxed\nend",
        "ri_mp2",
        "RI-MP2|pcSseg-3|NoFrozenCore|GIAO|CPCM(DMSO)|ExtremeSCF|ORCA_6.1.1",
    ),
    (
        "benchmark_dlpno_mp2_pcsseg3",
        "DLPNO-MP2 TightPNO pcSseg-3 AutoAux NMR NoFrozenCore CPCM(DMSO) RIJCOSX DefGrid3 ExtremeSCF",
        "%mp2\n  Density Relaxed\n  PAOOverlapThresh 1e-5\n  TScalePNO_LShift 0.1\nend",
        "dlpno_mp2",
        "DLPNO-MP2|TightPNO|pcSseg-3|NoFrozenCore|GIAO|CPCM(DMSO)|ExtremeSCF|ORCA_6.1.1",
    ),
    (
        "diagnostic_autoci_ccsdt_pcsseg4_gas",
        "AUTOCI-CCSD(T) NMR pcSseg-4 cc-pwCV5Z/C ExtremeSCF NoFrozenCore",
        "%autoci\n  STol 1e-10\n  Density Relaxed\nend",
        "autoci_ccsdt",
        "AUTOCI-CCSD(T)|pcSseg-4|cc-pwCV5Z/C|NoFrozenCore|GIAO|GAS_PHASE_DELTA_DIAGNOSTIC|ORCA_6.1.1",
    ),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def tms_xyz() -> str:
    mol = Chem.AddHs(Chem.MolFromSmiles("C[Si](C)(C)C"))
    params = AllChem.ETKDGv3()
    params.randomSeed = 6102026
    params.enforceChirality = True
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError("RDKit failed to embed TMS")
    # MMFF does not cover every silicon environment in every RDKit build.
    # Prefer MMFF when parameters are available and otherwise use UFF, which
    # supports the TMS starting geometry used only for the subsequent DFT opt.
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
    elif AllChem.UFFHasAllMoleculeParams(mol):
        AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
    else:
        raise RuntimeError("RDKit has neither MMFF nor UFF parameters for TMS")
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), "TMS same-level NMR reference; calculated model"]
    for atom in mol.GetAtoms():
        point = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():<2s} {point.x: .10f} {point.y: .10f} {point.z: .10f}")
    return "\n".join(lines) + "\n"


def main() -> int:
    ROOT.mkdir(parents=True, exist_ok=True)
    start = ROOT / "tms_start.xyz"
    start.write_text(tms_xyz(), encoding="utf-8", newline="\n")
    geometry = ROOT / "tms_opt.inp"
    geometry.write_text(
        """# CALCULATED MODEL - NOT EXPERIMENTAL DATA
! PBE0 D4 def2-TZVPP def2/J RIJCOSX CPCM(DMSO) VeryTightOpt Freq VeryTightSCF DefGrid3
%pal nprocs 16 end
%maxcore 6000
%geom MaxIter 500 end
%scf MaxIter 500 end
* xyzfile 0 1 tms_start.xyz
""",
        encoding="utf-8",
        newline="\n",
    )
    rows: list[dict[str, object]] = [
        {
            "job_id": "TMS_DMSO__opt_freq",
            "stage": "reference_optfreq",
            "method_id": "PBE0_D4_def2_TZVPP",
            "method_fingerprint": "PBE0-D4|def2-TZVPP|CPCM(DMSO)|VeryTightOpt|Freq|ORCA_6.1.1",
            "job_directory": "reference_jobs/TMS_DMSO",
            "portable_job_directory": "reference_jobs/TMS_DMSO",
            "input_file": "tms_opt.inp",
            "coordinate_file": "tms_start.xyz",
            "output_file": "tms_opt.out",
            "expected_input_sha256": sha256(geometry),
            "expected_coordinate_sha256": sha256(start),
            "enabled": "true",
            "readiness": "ready",
            "approval_status": "APPROVED_FOR_REFERENCE_GEOMETRY_CALCULATION",
            "notice": "CALCULATED_MODEL_NOT_EXPERIMENTAL",
        }
    ]
    for method_id, keywords, extra, resource, fingerprint in METHODS:
        directory = ROOT / method_id
        directory.mkdir(parents=True, exist_ok=True)
        input_path = directory / "tms_shield.inp"
        text = [
            "# CALCULATED MODEL - NOT EXPERIMENTAL DATA",
            "# Copy the reviewed reference_jobs/TMS_DMSO/tms_opt.xyz into this directory before activation.",
            f"! {keywords}",
            "%pal nprocs 16 end",
            "%maxcore 12000",
        ]
        if extra:
            text.extend((extra, ""))
        # ORCA 6.1 requires the EPR/NMR block after the coordinate specification.
        text.extend(("* xyzfile 0 1 tms_opt.xyz", "%eprnmr", "  Ori = GIAO", "  Tol = 1e-10", "end", ""))
        input_path.write_text("\n".join(text), encoding="utf-8", newline="\n")
        rel = directory.relative_to(PROJECT).as_posix()
        rows.append(
            {
                "job_id": f"TMS_DMSO__{method_id}",
                "stage": "reference_shielding",
                "method_id": method_id,
                "method_fingerprint": fingerprint,
                "resource_profile": resource,
                "job_directory": rel,
                "portable_job_directory": rel,
                "input_file": "tms_shield.inp",
                "coordinate_file": "tms_opt.xyz",
                "output_file": "tms_shield.out",
                "expected_input_sha256": sha256(input_path),
                "expected_coordinate_sha256": "TO_LOCK_AFTER_REFERENCE_GEOMETRY",
                "required_output_globs": (
                    '["tms_shield.autociresp_*"]'
                    if method_id.startswith("diagnostic_autoci")
                    else "[]"
                ),
                "required_helpers": (
                    '["orca_autoci"]'
                    if method_id.startswith("diagnostic_autoci")
                    else "[]"
                ),
                "result_use": (
                    "GAS_PHASE_DIAGNOSTIC_ONLY"
                    if method_id.startswith("diagnostic_autoci")
                    else "METHOD_MATCHED_REFERENCE"
                ),
                "enabled": "false",
                "readiness": "blocked_until_reference_geometry_review",
                "approval_status": "DRAFT_NOT_APPROVED",
                "notice": "CALCULATED_MODEL_NOT_EXPERIMENTAL",
            }
        )
    fields = list(dict.fromkeys(key for row in rows for key in row))
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"reference_jobs": len(rows), "manifest": str(MANIFEST), "status": "PREPARED"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
