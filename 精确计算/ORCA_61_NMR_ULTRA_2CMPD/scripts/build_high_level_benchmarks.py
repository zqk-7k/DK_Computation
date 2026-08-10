"""Build fail-closed RI-MP2, DLPNO-MP2, AUTOCI-CCSD(T) and VPT2 jobs.

The jobs are generated from APPROVED_FOR_NMR minimum gates but are emitted with
``enabled=false`` and ``readiness=manual_opt_in``.  They are diagnostics, not a
mandatory production dependency.  The script never launches ORCA and never
promotes a method to any validation status.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Mapping

from build_property_jobs import _validate_approved_minimum_gate, analyte_index_map, continuum_token
from workflow_schema import (
    NOTICE,
    WorkflowValidationError,
    atomic_write_bytes,
    atomic_write_csv,
    ensure_within,
    parse_xyz_frames,
    read_csv,
    read_json,
    read_text_auto,
    sha256_file,
)


PROJECT = Path(__file__).resolve().parents[1]
ISOTOPE = {"H": 1, "C": 13}


BRANCHES = (
    {
        "id": "benchmark_ri_mp2_pcsseg3",
        "top_n": 3,
        "solute_only": False,
        "keywords": "RIJK RI-MP2 pcSseg-3 AutoAux NMR NoFrozenCore CPCM({SOLVENT}) DefGrid3 ExtremeSCF",
        "extra": "%mp2\n  Density Relaxed\nend",
        "resource": "ri_mp2",
        "fingerprint": (
            "RI-MP2|pcSseg-3|NoFrozenCore|GIAO|CPCM({SOLVENT})|ExtremeSCF|"
            "EPRNMR_Tol=1e-12|Z_Tol=1e-12|ORCA_6.1.1"
        ),
    },
    {
        "id": "benchmark_dlpno_mp2_pcsseg3",
        "top_n": 1,
        "solute_only": False,
        "keywords": "DLPNO-MP2 TightPNO pcSseg-3 AutoAux NMR NoFrozenCore CPCM({SOLVENT}) RIJCOSX DefGrid3 ExtremeSCF",
        "extra": "%mp2\n  Density Relaxed\n  PAOOverlapThresh 1e-5\n  TScalePNO_LShift 0.1\nend",
        "resource": "dlpno_mp2",
        "fingerprint": (
            "DLPNO-MP2|TightPNO|pcSseg-3|NoFrozenCore|GIAO|CPCM({SOLVENT})|"
            "ExtremeSCF|EPRNMR_Tol=1e-12|Z_Tol=1e-12|ORCA_6.1.1"
        ),
    },
    {
        "id": "diagnostic_autoci_ccsdt_pcsseg4_gas",
        "top_n": 1,
        "solute_only": True,
        "keywords": "AUTOCI-CCSD(T) NMR pcSseg-4 cc-pwCV5Z/C ExtremeSCF NoFrozenCore",
        "extra": "%autoci\n  STol 1e-10\n  Density Relaxed\nend",
        "resource": "autoci_ccsdt",
        "fingerprint": (
            "AUTOCI-CCSD(T)|pcSseg-4|cc-pwCV5Z/C|NoFrozenCore|GIAO|"
            "AUTOCI_STol=1e-10|EPRNMR_Tol=1e-12|Z_Tol=1e-12|"
            "GAS_PHASE_DELTA_DIAGNOSTIC|ORCA_6.1.1"
        ),
        "required_helpers": ["orca_autoci"],
        "required_output_globs": ["benchmark.autociresp_*"],
    },
)

DIAGNOSTIC_APPROVED = "APPROVED_FOR_DIAGNOSTIC_EXECUTION"
VPT2_METHOD_ID = "vpt2_composite_diagnostic"


def manifest_activation_fields(
    stage: str, enabled: bool, diagnostic_gate_status: str
) -> dict[str, str]:
    """Map the scientific diagnostic gate to the HPC stage-token contract."""
    if not enabled:
        return {
            "enabled": "false",
            "readiness": "manual_opt_in",
            "approval_status": "DRAFT_NOT_APPROVED",
            "diagnostic_gate_status": diagnostic_gate_status,
        }
    approval = (
        "APPROVED_FOR_VPT2"
        if stage in {"vpt2_force_field", "vpt2_nmr"}
        else "APPROVED_FOR_HIGH_LEVEL_BENCHMARK"
    )
    return {
        "enabled": "true",
        "readiness": "ready",
        "approval_status": approval,
        "diagnostic_gate_status": diagnostic_gate_status,
    }


def _branch_status(gate: Mapping[str, object], branch_id: str) -> str:
    branches = gate.get("branches")
    if not isinstance(branches, dict) or branch_id not in branches:
        raise WorkflowValidationError(f"High-level gate lacks branch {branch_id}")
    value = branches[branch_id]
    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        return str(value.get("approval_status", ""))
    raise WorkflowValidationError(f"Malformed high-level gate branch {branch_id}")


def _validate_branch_scientific_reviews(gate: Mapping[str, object], branch_id: str) -> None:
    """Require every declared numerical/scientific review before activation.

    A single approval token is not enough for these deliberately expensive
    diagnostics.  The branch must also contain a non-empty mapping of named
    reviews, and every review must have been explicitly approved.  This keeps
    VPT2 geometry/displacement checks and correlated-method pilot checks
    machine-readable and fail-closed.
    """
    branches = gate.get("branches")
    branch = branches.get(branch_id) if isinstance(branches, dict) else None
    if not isinstance(branch, dict):
        raise WorkflowValidationError(
            f"Approved high-level branch {branch_id} lacks structured scientific reviews"
        )
    reviews = branch.get("scientific_reviews")
    if not isinstance(reviews, dict) or not reviews:
        raise WorkflowValidationError(
            f"Approved high-level branch {branch_id} lacks scientific_reviews"
        )
    unapproved = sorted(
        str(name)
        for name, status in reviews.items()
        if str(status) != "APPROVED_AFTER_REVIEW"
    )
    if unapproved:
        raise WorkflowValidationError(
            f"Approved high-level branch {branch_id} has unapproved scientific reviews: "
            + ", ".join(unapproved)
        )


def activation_state(
    project: Path, gate_path: Path, branch_id: str, require_autoci_helper: bool = False
) -> tuple[bool, str]:
    """Fail closed unless a branch is approved for one hash-bound ORCA machine."""
    gate = read_json(gate_path)
    status = _branch_status(gate, branch_id)
    if status != DIAGNOSTIC_APPROVED:
        return False, status or "DRAFT_NOT_APPROVED"
    if gate.get("approval_status") != DIAGNOSTIC_APPROVED:
        raise WorkflowValidationError(f"Branch {branch_id} approved under an unapproved parent gate")
    _validate_branch_scientific_reviews(gate, branch_id)
    if not str(gate.get("reviewer", "")).strip() or not str(gate.get("reviewed_at", "")).strip():
        raise WorkflowValidationError(f"Approved high-level gate lacks reviewer/time: {gate_path}")
    machine = gate.get("machine_binding")
    if not isinstance(machine, dict):
        raise WorkflowValidationError(f"Approved high-level gate lacks machine_binding: {gate_path}")
    if machine.get("orca_version_exact") != "6.1.1":
        raise WorkflowValidationError("High-level diagnostics require an exact ORCA 6.1.1 machine lock")
    for field in ("platform", "orca_executable_sha256", "preflight_file", "preflight_sha256"):
        if not str(machine.get(field, "")).strip():
            raise WorkflowValidationError(f"Machine binding lacks {field}")
    if not re.fullmatch(r"[0-9A-Fa-f]{64}", str(machine["orca_executable_sha256"])):
        raise WorkflowValidationError("Machine ORCA executable SHA-256 is malformed")
    preflight = ensure_within(Path(str(machine["preflight_file"])), project, "machine preflight")
    if sha256_file(preflight) != str(machine["preflight_sha256"]).upper():
        raise WorkflowValidationError("Machine preflight hash mismatch")
    preflight_payload = read_json(preflight)
    if preflight_payload.get("orca_version_exact") != "6.1.1":
        raise WorkflowValidationError("Machine preflight did not verify ORCA 6.1.1")
    helpers = set(preflight_payload.get("available_helpers", []))
    if require_autoci_helper and "orca_autoci" not in helpers:
        raise WorkflowValidationError("AutoCI activation requires orca_autoci on the bound machine")
    return True, DIAGNOSTIC_APPROVED


def isotope_masses_from_map(path: Path) -> dict[int, float]:
    masses: dict[int, float] = {}
    for row in read_csv(path):
        token = row.get("effective_isotope_mass_u", "").strip()
        if not token or token == "NATURAL_ABUNDANCE":
            continue
        index = int(row["atom_index_1based"])
        mass = float(token)
        if row.get("element", "").upper() != "H" or abs(mass - 2.01410177812) > 1e-9:
            raise WorkflowValidationError(f"Unsupported isotope override at atom {index} in {path}")
        masses[index] = mass
    return masses


def validated_vpt2_reference_record(
    project: Path, record: object
) -> tuple[bool, dict[str, object]]:
    """Verify the complete paired TMS static/VPT2 correction provenance."""
    if not isinstance(record, dict) or record.get("approval_status") != DIAGNOSTIC_APPROVED:
        return False, {}
    expected_fingerprint = (
        "PBE0-D4/def2-TZVPP_FF_then_PBE0/pcSseg-3_AvgProp_NMR|"
        "AnharmDisp=0.05|PropDisp=0.05|EPRNMR_Tol=1e-10|Z_Tol=1e-12|"
        "ISOTOPE_MASS_LOCKED|ORCA_6.1.1"
    )
    if record.get("method_fingerprint") != expected_fingerprint:
        raise WorkflowValidationError("Paired TMS VPT2 method fingerprint mismatch")
    for field in ("reference_id", "reference_geometry_sha256"):
        if not str(record.get(field, "")).strip():
            raise WorkflowValidationError(f"Paired TMS VPT2 record lacks {field}")
    if not re.fullmatch(r"[0-9A-Fa-f]{64}", str(record["reference_geometry_sha256"])):
        raise WorkflowValidationError("Paired TMS VPT2 geometry SHA-256 is malformed")
    roles = (
        "static_input",
        "static_output",
        "vpt2_ff_input",
        "vpt2_ff_output",
        "vpt2_nmr_input",
        "vpt2_nmr_output",
        "runner_receipt",
    )
    resolved: dict[str, Path] = {}
    for role in roles:
        path_field = f"{role}_file"
        hash_field = f"{role}_sha256"
        expected_hash = str(record.get(hash_field, ""))
        if not re.fullmatch(r"[0-9A-Fa-f]{64}", expected_hash):
            raise WorkflowValidationError(f"Paired TMS VPT2 {hash_field} is malformed")
        resolved[role] = ensure_within(
            Path(str(record.get(path_field, ""))), project, f"paired TMS VPT2 {role}"
        )
        if sha256_file(resolved[role]) != expected_hash.upper():
            raise WorkflowValidationError(f"Paired TMS VPT2 {role} hash mismatch")
    for role in ("static_output", "vpt2_ff_output", "vpt2_nmr_output"):
        if "ORCA TERMINATED NORMALLY" not in read_text_auto(resolved[role]):
            raise WorkflowValidationError(f"Paired TMS VPT2 {role} lacks normal termination")
    receipt = read_json(resolved["runner_receipt"])
    if receipt.get("status") != "completed" or receipt.get("orca_version_exact") != "6.1.1":
        raise WorkflowValidationError("Paired TMS VPT2 runner receipt is not completed ORCA 6.1.1")
    return True, dict(record)


def coordinate_block(
    optimized_xyz: Path,
    charge: int,
    multiplicity: int,
    isotope_masses: Mapping[int, float],
) -> list[str]:
    frames = parse_xyz_frames(optimized_xyz)
    if len(frames) != 1:
        raise WorkflowValidationError(f"Expected one optimized XYZ frame: {optimized_xyz}")
    frame = frames[0]
    lines = [f"* xyz {charge} {multiplicity}"]
    for index, (symbol, xyz) in enumerate(
        zip(frame.symbols, frame.coordinates, strict=True), start=1
    ):
        suffix = f" M={isotope_masses[index]:.11f}" if index in isotope_masses else ""
        lines.append(
            f"  {symbol:<2s} {xyz[0]: .12f} {xyz[1]: .12f} {xyz[2]: .12f}{suffix}"
        )
    lines.extend(("*", ""))
    return lines


def csv_bytes(rows: list[dict[str, object]], fields: tuple[str, ...]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def render_nmr(
    keywords: str,
    extra: str,
    charge: int,
    multiplicity: int,
    index_map: dict[str, list[int]],
    ranks: int,
    maxcore: int,
    optimized_xyz: Path,
    isotope_masses: Mapping[int, float],
) -> str:
    lines = [
        f"# {NOTICE}",
        "# Manual high-level diagnostic; complete approved electronic model retained.",
        f"! {keywords}",
        f"%pal nprocs {ranks} end",
        f"%maxcore {maxcore}",
        "%method Z_Tol 1e-12 end",
    ]
    if extra:
        lines.extend((extra, ""))
    lines.extend(coordinate_block(optimized_xyz, charge, multiplicity, isotope_masses))
    lines.extend(("%eprnmr", "  Ori = GIAO", "  Tol = 1e-12"))
    for element in ("H", "C"):
        indices = index_map.get(element, [])
        if indices:
            lines.append(
                f"  Nuclei = {','.join(map(str, indices))} {{ shift, ist = {ISOTOPE[element]} }};"
            )
    lines.extend(("end", ""))
    return "\n".join(lines)


def render_vpt2_ff(
    charge: int,
    multiplicity: int,
    solvent: str,
    ranks: int,
    maxcore: int,
    optimized_xyz: Path,
    isotope_masses: Mapping[int, float],
) -> str:
    lines = f"""# {NOTICE}
# DIAGNOSTIC ONLY: cannot enter main delta without paired analyte/reference corrections.
! PBE0 D4 def2-TZVPP def2/J RIJCOSX CPCM({solvent}) DefGrid3 ExtremeSCF VPT2
%pal nprocs {ranks} end
%maxcore {maxcore}
%vpt2
  VPT2 On
  AnharmDisp 0.05
  HessianCutoff 1e-12
end
%method
  Z_Tol 1e-12
end
""".splitlines()
    lines.extend(coordinate_block(optimized_xyz, charge, multiplicity, isotope_masses))
    return "\n".join(lines)


def render_vpt2_nmr(
    charge: int,
    multiplicity: int,
    solvent: str,
    index_map: dict[str, list[int]],
    ranks: int,
    maxcore: int,
    optimized_xyz: Path,
    isotope_masses: Mapping[int, float],
) -> str:
    lines = [
        f"# {NOTICE}",
        f"! PBE0 pcSseg-3 AutoAux CPCM({solvent}) DefGrid3 ExtremeSCF NMR",
        f"%pal nprocs {ranks} end",
        f"%maxcore {maxcore}",
        "%geom InHess Read InHessName \"vpt2_ff.hess\" end",
        "%vpt2",
        "  VPT2 On",
        "  AvgProp NMR",
        "  PropDisp 0.05",
        "  HessianCutoff 1e-12",
        "end",
        "%method Z_Tol 1e-12 end",
        "# Requires both vpt2_ff.hess and vpt2_ff.vpt2 from the same isotope-aware force-field job.",
    ]
    lines.extend(coordinate_block(optimized_xyz, charge, multiplicity, isotope_masses))
    lines.extend([
        "%eprnmr",
        "  Ori = GIAO",
        "  Tol = 1e-10",
    ])
    for element in ("H", "C"):
        indices = index_map.get(element, [])
        if indices:
            lines.append(
                f"  Nuclei = {','.join(map(str, indices))} {{ shift, ist = {ISOTOPE[element]} }};"
            )
    lines.extend(("end", ""))
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=PROJECT)
    parser.add_argument("--gate", type=Path, action="append", required=True)
    parser.add_argument("--dft-manifest", type=Path, default=Path("data/DFT_explicit_job_manifest.csv"))
    parser.add_argument(
        "--high-level-gate",
        type=Path,
        default=Path("jobs/high_level_manual_opt_in_gate.json"),
    )
    parser.add_argument("--output-root", type=Path, default=Path("04_high_level_benchmarks"))
    parser.add_argument("--manifest", type=Path, default=Path("data/high_level_benchmark_manifest.csv"))
    args = parser.parse_args()

    project = args.project.resolve()
    dft_manifest = ensure_within(project / args.dft_manifest, project, "DFT manifest")
    output_root = ensure_within(project / args.output_root, project, "benchmark output root")
    manifest_path = ensure_within(project / args.manifest, project, "benchmark manifest")
    high_level_gate = ensure_within(
        project / args.high_level_gate, project, "high-level activation gate"
    )
    high_level_gate_payload = read_json(high_level_gate)
    high_level_gate_hash = sha256_file(high_level_gate)
    dft_rows = read_csv(dft_manifest)
    rows_by_id = {row["job_id"]: row for row in dft_rows}
    resources = json.loads((project / "config/resource_policy.json").read_text(encoding="utf-8"))["stage_resources"]
    manifest_rows: list[dict[str, object]] = []

    for raw_gate in args.gate:
        gate_path = ensure_within(project / raw_gate, project, "minimum gate")
        gate, states, approved, _weight_sets = _validate_approved_minimum_gate(project, gate_path, dft_rows)
        grouped: dict[str, list[str]] = defaultdict(list)
        for job_id in approved:
            grouped[rows_by_id[job_id]["scenario_id"]].append(job_id)
        rank: dict[str, int] = {}
        for scenario_id, job_ids in grouped.items():
            ordered = sorted(
                job_ids,
                key=lambda item: float(
                    approved[item]["thermochemistry_by_nucleus"][
                        "1H"
                        if "1H" in approved[item]["thermochemistry_by_nucleus"]
                        else sorted(approved[item]["thermochemistry_by_nucleus"])[0]
                    ]["qrrho_gibbs_free_energy_hartree"]
                ),
            )
            rank.update({job_id: position for position, job_id in enumerate(ordered, start=1)})

        for job_id, decision in approved.items():
            source = rows_by_id[job_id]
            evidence = states[job_id]
            scenario_id = source["scenario_id"]
            optimized = Path(evidence["artifacts"]["optimized_xyz"]["path"])
            atom_map = Path(evidence["job_directory"]) / "nmr_atom_map_v2.csv"
            index_map = analyte_index_map(atom_map)
            isotope_masses = isotope_masses_from_map(atom_map)
            charge = int(source["charge"])
            multiplicity = int(source["multiplicity"])
            solvent = continuum_token(source["continuum_keyword"])
            explicit_count = int(source.get("explicit_solvent_count", "0") or 0)

            for branch in BRANCHES:
                if rank[job_id] > int(branch["top_n"]):
                    continue
                if branch["solute_only"] and explicit_count != 0:
                    continue
                method_id = str(branch["id"])
                branch_enabled, branch_status = activation_state(
                    project,
                    high_level_gate,
                    method_id,
                    require_autoci_helper=method_id.startswith("diagnostic_autoci"),
                )
                resource = resources[str(branch["resource"])]
                keywords = str(branch["keywords"]).replace("{SOLVENT}", solvent)
                fingerprint = str(branch["fingerprint"]).replace("{SOLVENT}", solvent)
                text = render_nmr(
                    keywords,
                    str(branch["extra"]),
                    charge,
                    multiplicity,
                    index_map,
                    int(resource["mpi_ranks"]),
                    int(resource["maxcore_mb_per_rank"]),
                    optimized,
                    isotope_masses,
                )
                directory = output_root / gate["figure_id"] / scenario_id / job_id / method_id
                directory.mkdir(parents=True, exist_ok=True)
                (directory / "optimized.xyz").write_bytes(optimized.read_bytes())
                (directory / "nmr_atom_map_v2.csv").write_bytes(atom_map.read_bytes())
                (directory / "source_isotope_map.csv").write_bytes(
                    (Path(evidence["job_directory"]) / "source_isotope_map.csv").read_bytes()
                )
                (directory / "benchmark.inp").write_text(text, encoding="utf-8", newline="\n")
                rel = directory.relative_to(project).as_posix()
                manifest_rows.append(
                    {
                        "job_id": f"{job_id}__{method_id}",
                        "figure_id": gate["figure_id"],
                        "scenario_id": scenario_id,
                        "source_job_id": job_id,
                        "fixed_composition_rank": rank[job_id],
                        "stage": "high_level_benchmark",
                        "method_id": method_id,
                        "method_fingerprint": fingerprint,
                        "job_directory": rel,
                        "portable_job_directory": rel,
                        "input_file": "benchmark.inp",
                        "coordinate_file": "optimized.xyz",
                        "output_file": "benchmark.out",
                        "expected_input_sha256": sha256_file(directory / "benchmark.inp"),
                        "expected_coordinate_sha256": sha256_file(directory / "optimized.xyz"),
                        "source_isotope_map_sha256": sha256_file(
                            directory / "source_isotope_map.csv"
                        ),
                        "required_artifacts": json.dumps(
                            (
                                ["benchmark.gbw"]
                                if method_id.startswith("diagnostic_autoci")
                                else ["benchmark.gbw", "benchmark.property.txt"]
                            ),
                            separators=(",", ":"),
                        ),
                        "required_output_globs": json.dumps(
                            branch.get("required_output_globs", []), separators=(",", ":")
                        ),
                        "required_helpers": json.dumps(
                            branch.get("required_helpers", []), separators=(",", ":")
                        ),
                        "high_level_gate_sha256": high_level_gate_hash,
                        **manifest_activation_fields(
                            "high_level_benchmark", branch_enabled, branch_status
                        ),
                        "result_use": "DIAGNOSTIC_ONLY_NEVER_AUTOMATICALLY_SUPERSEDES_SOLVATED_PRODUCTION",
                        "notice": NOTICE,
                    }
                )

            if rank[job_id] == 1 and explicit_count <= 1:
                resource = resources["vpt2"]
                vpt2_enabled, vpt2_status = activation_state(
                    project, high_level_gate, VPT2_METHOD_ID
                )
                reference_ready, reference_correction = validated_vpt2_reference_record(
                    project,
                    high_level_gate_payload.get("paired_reference_vpt2_correction", {}),
                )
                directory = output_root / gate["figure_id"] / scenario_id / job_id / VPT2_METHOD_ID
                directory.mkdir(parents=True, exist_ok=True)
                (directory / "optimized.xyz").write_bytes(optimized.read_bytes())
                (directory / "nmr_atom_map_v2.csv").write_bytes(atom_map.read_bytes())
                (directory / "source_isotope_map.csv").write_bytes(
                    (Path(evidence["job_directory"]) / "source_isotope_map.csv").read_bytes()
                )
                (directory / "vpt2_ff.inp").write_text(
                    render_vpt2_ff(
                        charge,
                        multiplicity,
                        solvent,
                        int(resource["mpi_ranks"]),
                        int(resource["maxcore_mb_per_rank"]),
                        optimized,
                        isotope_masses,
                    ),
                    encoding="utf-8",
                    newline="\n",
                )
                (directory / "vpt2_nmr.inp").write_text(
                    render_vpt2_nmr(
                        charge,
                        multiplicity,
                        solvent,
                        index_map,
                        int(resource["mpi_ranks"]),
                        int(resource["maxcore_mb_per_rank"]),
                        optimized,
                        isotope_masses,
                    ),
                    encoding="utf-8",
                    newline="\n",
                )
                rel = directory.relative_to(project).as_posix()
                for stage, input_file, output_file, dependency in (
                    ("vpt2_force_field", "vpt2_ff.inp", "vpt2_ff.out", ""),
                    ("vpt2_nmr", "vpt2_nmr.inp", "vpt2_nmr.out", f"{job_id}__vpt2_force_field"),
                ):
                    manifest_rows.append(
                        {
                            "job_id": f"{job_id}__{stage}",
                            "figure_id": gate["figure_id"],
                            "scenario_id": scenario_id,
                            "source_job_id": job_id,
                            "fixed_composition_rank": 1,
                            "stage": stage,
                            "method_id": VPT2_METHOD_ID,
                            "method_fingerprint": "PBE0-D4/def2-TZVPP_FF_then_PBE0/pcSseg-3_AvgProp_NMR|AnharmDisp=0.05|PropDisp=0.05|EPRNMR_Tol=1e-10|Z_Tol=1e-12|ISOTOPE_MASS_LOCKED|ORCA_6.1.1",
                            "job_directory": rel,
                            "portable_job_directory": rel,
                            "input_file": input_file,
                            "coordinate_file": "optimized.xyz",
                            "output_file": output_file,
                            "depends_on_job_id": dependency,
                            "depends_on_stage": (
                                "vpt2_force_field" if stage == "vpt2_nmr" else ""
                            ),
                            "depends_on_artifacts": json.dumps(
                                ["vpt2_ff.hess", "vpt2_ff.vpt2"] if stage == "vpt2_nmr" else [],
                                separators=(",", ":"),
                            ),
                            "required_artifacts": json.dumps(
                                (
                                    ["vpt2_ff.hess", "vpt2_ff.vpt2"]
                                    if stage == "vpt2_force_field"
                                    else ["vpt2_nmr.vpt2", "vpt2_nmr.property.txt"]
                                ),
                                separators=(",", ":"),
                            ),
                            "expected_input_sha256": sha256_file(directory / input_file),
                            "expected_coordinate_sha256": sha256_file(directory / "optimized.xyz"),
                            "source_isotope_map_sha256": sha256_file(
                                directory / "source_isotope_map.csv"
                            ),
                            "high_level_gate_sha256": high_level_gate_hash,
                            **manifest_activation_fields(stage, vpt2_enabled, vpt2_status),
                            "main_delta_eligibility": (
                                "ELIGIBLE_ONLY_FOR_REVIEWED_COMPOSITE"
                                if reference_ready
                                else "BLOCKED_PAIRED_TMS_VPT2_CORRECTION_REQUIRED"
                            ),
                            "composite_equation": (
                                "delta_final = delta_static_high + "
                                "[(sigma_ref_vib_low-sigma_ref_static_low) - "
                                "(sigma_analyte_vib_low-sigma_analyte_static_low)]"
                            ),
                            "analyte_static_and_vibrational_low_level_values_required": "true",
                            "paired_reference_vpt2_record": json.dumps(
                                reference_correction if reference_ready else {},
                                separators=(",", ":"),
                                sort_keys=True,
                            ),
                            "notice": NOTICE,
                        }
                    )

    if not manifest_rows:
        raise WorkflowValidationError("No high-level benchmark jobs were generated")
    fields = tuple(dict.fromkeys(key for row in manifest_rows for key in row))
    atomic_write_csv(manifest_path, manifest_rows, fields)
    print(json.dumps({"jobs": len(manifest_rows), "manifest": str(manifest_path), "status": "MANUAL_OPT_IN_REQUIRED"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
