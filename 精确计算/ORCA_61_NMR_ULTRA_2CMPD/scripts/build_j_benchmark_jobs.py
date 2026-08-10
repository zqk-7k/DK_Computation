"""Build the representative PBE/pcJ-3 versus PBE0/pcJ-2 J benchmark.

This closes the J-method gate without pre-approving it: jobs are executable
only after the conformer/minimum gate, while numerical agreement and the final
production-method decision still require a named review of hash-bound outputs.
The 13C presentation remains proton-decoupled; H-C pairs here are method
diagnostics, not instructions to draw 13C-H splitting in the main spectrum.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

from build_property_jobs import (
    _validate_approved_minimum_gate,
    analyte_index_map,
    continuum_token,
    render_j_input,
)
from workflow_schema import (
    NOTICE,
    WorkflowValidationError,
    atomic_write_csv,
    atomic_write_json,
    ensure_within,
    parse_xyz_frames,
    read_csv,
    read_json,
    sha256_file,
)


PROJECT = Path(__file__).resolve().parents[1]
RADIUS_A = 6.05
BRANCHES = (
    {
        "method_id": "production_pbe_pcj3",
        "fingerprint": (
            "PBE|pcJ-3|AutoAux|CPCM(actual_solvent)|VeryTightSCF|DefGrid3|"
            "EPRNMR_Tol=1e-10|Z_Tol=1e-12|SpinSpinRThresh=6.05A|ORCA_6.1.1"
        ),
        "keywords": "PBE pcJ-3 AutoAux VeryTightSCF DefGrid3 CPCM({CONTINUUM})",
    },
    {
        "method_id": "benchmark_pbe0_pcj2",
        "fingerprint": (
            "PBE0|pcJ-2|AutoAux|CPCM(actual_solvent)|VeryTightSCF|DefGrid3|"
            "EPRNMR_Tol=1e-10|Z_Tol=1e-12|SpinSpinRThresh=6.05A|ORCA_6.1.1"
        ),
        "keywords": "PBE0 pcJ-2 AutoAux VeryTightSCF DefGrid3 CPCM({CONTINUUM})",
    },
)


def j_benchmark_execution_state(
    project: Path,
    gate_path: Path,
    minimum_gate_hashes: list[str],
    resource_policy_hash: str,
) -> tuple[bool, str]:
    """Require a named, machine-bound gate before expensive benchmark execution."""
    gate = read_json(gate_path)
    approved_token = "APPROVED_FOR_J_BENCHMARK_CALCULATION"
    if gate.get("approval_status") != approved_token:
        return False, "DRAFT_NOT_APPROVED"
    if not str(gate.get("reviewer", "")).strip() or not str(gate.get("reviewed_at", "")).strip():
        raise WorkflowValidationError("Approved J benchmark execution gate lacks reviewer/time")
    supplied_minimum_hashes = gate.get("approved_minimum_gate_sha256")
    if not isinstance(supplied_minimum_hashes, list) or {
        str(value).upper() for value in supplied_minimum_hashes
    } != {value.upper() for value in minimum_gate_hashes}:
        raise WorkflowValidationError("J benchmark execution gate minimum-gate hashes mismatch")
    if str(gate.get("resource_policy_sha256", "")).upper() != resource_policy_hash.upper():
        raise WorkflowValidationError("J benchmark execution gate resource-policy hash mismatch")
    machine = gate.get("machine_binding")
    if not isinstance(machine, dict) or machine.get("orca_version_exact") != "6.1.1":
        raise WorkflowValidationError("J benchmark execution gate lacks ORCA 6.1.1 machine binding")
    for field in ("preflight_file", "preflight_sha256", "orca_executable_sha256"):
        if not str(machine.get(field, "")).strip():
            raise WorkflowValidationError(f"J benchmark machine binding lacks {field}")
    if not re.fullmatch(r"[0-9A-Fa-f]{64}", str(machine["orca_executable_sha256"])):
        raise WorkflowValidationError("J benchmark ORCA executable hash is malformed")
    preflight = ensure_within(Path(str(machine["preflight_file"])), project, "J benchmark preflight")
    if sha256_file(preflight) != str(machine["preflight_sha256"]).upper():
        raise WorkflowValidationError("J benchmark machine preflight hash mismatch")
    if read_json(preflight).get("orca_version_exact") != "6.1.1":
        raise WorkflowValidationError("J benchmark preflight did not verify ORCA 6.1.1")
    return True, approved_token


def benchmark_pairs(index_map: dict[str, list[int]], xyz_path: Path) -> list[dict[str, object]]:
    """Include all analyte H-H and H-C pairs inside the locked distance screen."""
    frames = parse_xyz_frames(xyz_path)
    if len(frames) != 1:
        raise WorkflowValidationError(f"Expected one optimized frame: {xyz_path}")
    frame = frames[0]
    candidates: list[tuple[int, int, str, str]] = []
    hydrogens = sorted(index_map.get("H", []))
    carbons = sorted(index_map.get("C", []))
    candidates.extend((a, b, "H", "H") for pos, a in enumerate(hydrogens) for b in hydrogens[pos + 1 :])
    candidates.extend((h, c, "H", "C") for h in hydrogens for c in carbons)
    rows: list[dict[str, object]] = []
    for left, right, left_element, right_element in candidates:
        a = frame.coordinates[left - 1]
        b = frame.coordinates[right - 1]
        distance = math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b, strict=True)))
        if distance > RADIUS_A:
            continue
        rows.append(
            {
                "left_atom_index_1based_nuclei": left,
                "right_atom_index_1based_nuclei": right,
                "left_atom_index_0based_spinspin": left - 1,
                "right_atom_index_0based_spinspin": right - 1,
                "left_element": left_element,
                "right_element": right_element,
                "pair_class": f"{left_element}-{right_element}",
                "distance_A": f"{distance:.6f}",
                "scope": "ANALYTE_ONLY_METHOD_BENCHMARK",
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=PROJECT)
    parser.add_argument("--gate", type=Path, action="append", required=True)
    parser.add_argument("--dft-manifest", type=Path, default=Path("data/DFT_explicit_job_manifest.csv"))
    parser.add_argument("--output-root", type=Path, default=Path("04_j_method_benchmark"))
    parser.add_argument("--manifest", type=Path, default=Path("data/J_method_benchmark_job_manifest.csv"))
    parser.add_argument("--coverage", type=Path, default=Path("data/J_method_benchmark_coverage.json"))
    parser.add_argument(
        "--execution-gate",
        type=Path,
        default=Path("jobs/j_benchmark_execution_gate.json"),
    )
    args = parser.parse_args()

    project = args.project.resolve()
    dft_manifest = (project / args.dft_manifest).resolve()
    dft_rows = read_csv(dft_manifest)
    rows_by_id = {row["job_id"]: row for row in dft_rows}
    resource_policy_path = project / "config" / "resource_policy.json"
    resource_policy_hash = sha256_file(resource_policy_path)
    resources = read_json(resource_policy_path)["stage_resources"]["j_coupling"]
    ranks = int(resources["mpi_ranks"])
    maxcore = int(resources["maxcore_mb_per_rank"])
    output_root = (project / args.output_root).resolve()
    manifest_rows: list[dict[str, object]] = []
    coverage_rows: list[dict[str, object]] = []
    minimum_gate_paths = [(project / value).resolve() for value in args.gate]
    minimum_gate_hashes = [sha256_file(path) for path in minimum_gate_paths]
    execution_gate = (project / args.execution_gate).resolve()
    execution_enabled, execution_status = j_benchmark_execution_state(
        project, execution_gate, minimum_gate_hashes, resource_policy_hash
    )
    execution_gate_hash = sha256_file(execution_gate)

    for gate_path in minimum_gate_paths:
        gate, states, approved, _weights = _validate_approved_minimum_gate(
            project, gate_path, dft_rows
        )
        by_scenario: dict[str, list[str]] = defaultdict(list)
        for job_id in approved:
            by_scenario[rows_by_id[job_id]["scenario_id"]].append(job_id)
        selected: list[str] = []
        for scenario_id, job_ids in sorted(by_scenario.items()):
            selected.append(
                min(
                    job_ids,
                    key=lambda item: float(
                        approved[item]["thermochemistry_by_nucleus"][
                            "1H"
                        ]["qrrho_gibbs_free_energy_hartree"]
                    ),
                )
            )
        for job_id in selected:
            source = rows_by_id[job_id]
            evidence = states[job_id]
            xyz_path = Path(evidence["artifacts"]["optimized_xyz"]["path"])
            atom_map = Path(evidence["job_directory"]) / "nmr_atom_map_v2.csv"
            index_map = analyte_index_map(atom_map)
            pairs = benchmark_pairs(index_map, xyz_path)
            pair_counts = Counter(str(row["pair_class"]) for row in pairs)
            if not pairs or not pair_counts.get("H-H"):
                raise WorkflowValidationError(f"J benchmark lacks H-H coverage: {job_id}")
            coverage_rows.append(
                {
                    "figure_id": gate["figure_id"],
                    "scenario_id": source["scenario_id"],
                    "source_job_id": job_id,
                    "pair_counts": dict(pair_counts),
                    "distance_screen_A": RADIUS_A,
                }
            )
            for branch in BRANCHES:
                text = render_j_input(
                    int(source["charge"]),
                    int(source["multiplicity"]),
                    continuum_token(source["continuum_keyword"]),
                    index_map,
                    pairs,
                    ranks,
                    maxcore,
                    str(branch["keywords"]),
                )
                directory = (
                    output_root
                    / str(gate["figure_id"])
                    / source["scenario_id"]
                    / job_id
                    / str(branch["method_id"])
                )
                directory.mkdir(parents=True, exist_ok=True)
                (directory / "optimized.xyz").write_bytes(xyz_path.read_bytes())
                (directory / "nmr_atom_map_v2.csv").write_bytes(atom_map.read_bytes())
                (directory / "j_benchmark.inp").write_text(text, encoding="utf-8", newline="\n")
                fields = tuple(pairs[0])
                atomic_write_csv(directory / "j_benchmark_pair_allowlist.csv", pairs, fields)
                rel = directory.relative_to(project).as_posix()
                manifest_rows.append(
                    {
                        "job_id": f"{job_id}__{branch['method_id']}",
                        "figure_id": gate["figure_id"],
                        "scenario_id": source["scenario_id"],
                        "source_job_id": job_id,
                        "stage": "j_method_benchmark",
                        "method_id": branch["method_id"],
                        "method_fingerprint": branch["fingerprint"],
                        "job_directory": rel,
                        "portable_job_directory": rel,
                        "input_file": "j_benchmark.inp",
                        "coordinate_file": "optimized.xyz",
                        "output_file": "j_benchmark.out",
                        "required_artifacts": '["j_benchmark.gbw","j_benchmark.property.txt"]',
                        "expected_input_sha256": sha256_file(directory / "j_benchmark.inp"),
                        "expected_coordinate_sha256": sha256_file(directory / "optimized.xyz"),
                        "approved_minimum_gate_sha256": sha256_file(gate_path),
                        "resource_policy_sha256": resource_policy_hash,
                        "j_benchmark_execution_gate_sha256": execution_gate_hash,
                        "pair_count": len(pairs),
                        "pair_class_counts": json.dumps(dict(pair_counts), separators=(",", ":")),
                        "distance_screen_A": RADIUS_A,
                        "enabled": "true" if execution_enabled else "false",
                        "readiness": "ready" if execution_enabled else "manual_opt_in",
                        "approval_status": execution_status,
                        "result_use": "METHOD_SELECTION_ONLY_NOT_13C_PROTON_SPLITTING",
                        "notice": NOTICE,
                    }
                )

    figures = {str(row["figure_id"]) for row in coverage_rows}
    required_figures = {"F3_005", "F3_015"}
    if figures != required_figures:
        raise WorkflowValidationError(
            f"Representative J benchmark must cover both compounds; got {sorted(figures)}"
        )
    if not any(row["pair_counts"].get("H-C", 0) for row in coverage_rows):
        raise WorkflowValidationError("Representative J benchmark lacks H-C diagnostic coverage")
    manifest = (project / args.manifest).resolve()
    fields = tuple(dict.fromkeys(key for row in manifest_rows for key in row))
    atomic_write_csv(manifest, manifest_rows, fields)
    coverage = (project / args.coverage).resolve()
    atomic_write_json(
        coverage,
        {
            "schema_version": "1.0",
            "notice": NOTICE,
            "status": "INPUTS_PREPARED_RESULTS_NOT_REVIEWED",
            "locked_distance_screen_A": RADIUS_A,
            "coverage": coverage_rows,
            "required_review": [
                "Compare PBE/pcJ-3 and PBE0/pcJ-2 on identical pairs.",
                "Report sensitivity to the 6.05 A pair screen before approving truncation.",
                "Do not show H-C splitting in the routine proton-decoupled 13C spectrum.",
            ],
        },
    )
    print(json.dumps({"jobs": len(manifest_rows), "coverage": str(coverage)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
