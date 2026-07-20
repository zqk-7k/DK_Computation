#!/usr/bin/env python3
"""Deadline-aware whole-node ORCA workflow for the 52 Figure 3 compounds.

The scheduler uses a physical-core token pool.  It first produces a complete
low-cost shielding layer, then replaces successful compounds with a higher-level
result while time remains.  It never labels a fallback result as a refined one.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import time
import uuid
from collections import deque
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

from rdkit import Chem
from rdkit.Chem import Lipinski

from common import (
    ROOT,
    CorePool,
    artifact_manifest_issues,
    build_opt_input,
    build_shield_input,
    find_final_geometry,
    frozen_geometry_issues,
    geometry_same_composition,
    parse_shielding_summary,
    physical_cpu_ids,
    read_banner_csv,
    read_csv,
    read_json,
    read_xyz,
    runtime_environment_identity,
    run_orca_job,
    sha256,
    terminate_active_processes,
    utc_stamp,
    write_csv,
    write_json,
)


STOP_REQUESTED = False


def captured_text(value: str | bytes | None) -> str:
    """Normalize TimeoutExpired capture fields across CPython platforms."""

    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


def deadline_completed_process(
    command: Sequence[str], exc: subprocess.TimeoutExpired, marker: str
) -> subprocess.CompletedProcess[str]:
    return subprocess.CompletedProcess(
        list(command),
        124,
        captured_text(exc.stdout),
        captured_text(exc.stderr) + f"\n{marker}\n",
    )


def request_stop(signum: int, _frame: Any) -> None:
    global STOP_REQUESTED
    STOP_REQUESTED = True
    print(
        f"[{utc_stamp()}] received signal {signum}; stopping new work and terminating active ORCA job steps",
        flush=True,
    )
    terminated = terminate_active_processes(grace_seconds=12.0)
    print(f"[{utc_stamp()}] termination requested for process groups {terminated}", flush=True)


@dataclass
class Task:
    task_id: str
    cores: int
    description: str
    runner: Callable[[list[int]], dict[str, Any]]


def method_by_id(methods: dict[str, Any], method_id: str) -> dict[str, Any]:
    matches = [item for item in methods["shielding_methods"] if item["method_id"] == method_id]
    if len(matches) != 1:
        raise KeyError(f"Unknown shielding method {method_id}")
    return matches[0]


def workflow_fingerprint() -> tuple[str, dict[str, str]]:
    relative_paths = [
        "config/hardware_104physical_24h.json",
        "config/methods_24h.json",
        "config/appearance_400MHz.json",
        "config/primary_species.csv",
        "config/input_manifest_sha256.csv",
        "requirements.txt",
        "scripts/common.py",
        "scripts/run_24h_pipeline.py",
        "scripts/postprocess_nmr.py",
        "scripts/audit_outputs.py",
        "scripts/finalize_manifest.py",
        "scripts/generate_formula_hrms.py",
        "scripts/generate_hrms_all.py",
    ]
    hashes = {relative: sha256(ROOT / relative) for relative in relative_paths}
    runtime_identity = runtime_environment_identity()
    for key, value in sorted(runtime_identity.items()):
        hashes[f"runtime_environment/{key}"] = value
    digest = hashlib.sha256()
    for relative, value in sorted(hashes.items()):
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(value.encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest().upper(), hashes


def latest_matching_success(raw_root: Path, job_id: str, input_text: str, geometry: Path) -> tuple[dict[str, Any], Path] | None:
    job_root = raw_root / job_id
    input_hash = hashlib.sha256(input_text.encode("utf-8")).hexdigest().upper()
    geometry_hash = sha256(geometry)
    for attempt in sorted(job_root.glob("attempt_*"), reverse=True):
        record_path = attempt / "run_record.json"
        if not record_path.is_file():
            continue
        record = read_json(record_path)
        recorded_output = Path(str(record.get("output_file", "")))
        portable_output = attempt / recorded_output.name if recorded_output.name else Path()
        output_exists = portable_output.is_file() or recorded_output.is_file()
        if (
            record.get("status") == "COMPLETE"
            and record.get("input_sha256") == input_hash
            and record.get("geometry_sha256") == geometry_hash
            and output_exists
        ):
            if portable_output.is_file():
                record["output_file"] = str(portable_output)
            return record, attempt
    return None


def select_job_timeout(stage_limit: float, absolute_deadline: float) -> float:
    return max(0.0, min(stage_limit, absolute_deadline - time.monotonic() - 30.0))


def run_or_reuse(
    *,
    job_id: str,
    input_text: str,
    geometry: Path,
    orca: Path,
    scratch: Path,
    raw_root: Path,
    cpu_ids: Sequence[int],
    timeout: float,
    basename: str,
) -> tuple[dict[str, Any], Path]:
    reused = latest_matching_success(raw_root, job_id, input_text, geometry)
    if reused:
        record, path = reused
        record = dict(record)
        record["reused"] = True
        return record, path
    if timeout < 60:
        return {
            "job_id": job_id,
            "status": "SKIPPED_INSUFFICIENT_TIME",
            "reused": False,
            "elapsed_seconds": 0.0,
            "error": f"remaining timeout {timeout:.1f} s",
        }, raw_root / job_id
    record, path = run_orca_job(
        job_id=job_id,
        input_text=input_text,
        geometry_path=geometry,
        orca_path=orca,
        scratch_root=scratch,
        raw_root=raw_root,
        cpu_ids=cpu_ids,
        timeout_seconds=timeout,
        basename=basename,
    )
    value = asdict(record)
    value["reused"] = False
    return value, path


def run_task_batch(
    *,
    stage_name: str,
    tasks: list[Task],
    cpu_ids: Sequence[int],
    run_dir: Path,
    stop_launch_at: float,
    status_context: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    """Launch tasks whenever enough physical-core tokens are available."""

    pool = CorePool(cpu_ids)
    pending = deque(tasks)
    running: dict[Future[dict[str, Any]], tuple[Task, list[int]]] = {}
    outcomes: dict[str, dict[str, Any]] = {}
    stage_start = time.monotonic()
    last_report = 0.0
    with ThreadPoolExecutor(max_workers=max(1, len(tasks))) as executor:
        while pending or running:
            now = time.monotonic()
            launch_allowed = now < stop_launch_at and not STOP_REQUESTED
            if launch_allowed and pending:
                made_progress = True
                while made_progress and pending:
                    made_progress = False
                    for index, task in enumerate(pending):
                        allocation = pool.acquire(task.cores)
                        if allocation is None:
                            continue
                        del pending[index]
                        print(
                            f"[{utc_stamp()}] START {stage_name}/{task.task_id} "
                            f"cores={allocation} {task.description}",
                            flush=True,
                        )
                        future = executor.submit(task.runner, allocation)
                        running[future] = (task, allocation)
                        made_progress = True
                        break

            if not running:
                if pending:
                    reason = "STOP_REQUESTED" if STOP_REQUESTED else "STAGE_LAUNCH_CUTOFF"
                    while pending:
                        task = pending.popleft()
                        outcomes[task.task_id] = {
                            "task_id": task.task_id,
                            "status": f"NOT_LAUNCHED_{reason}",
                            "description": task.description,
                        }
                break

            done, _ = wait(running, timeout=5.0, return_when=FIRST_COMPLETED)
            for future in done:
                task, allocation = running.pop(future)
                pool.release(allocation)
                try:
                    outcome = future.result()
                except Exception as exc:
                    outcome = {
                        "task_id": task.task_id,
                        "status": "PYTHON_EXCEPTION",
                        "error": f"{type(exc).__name__}: {exc}",
                    }
                outcomes[task.task_id] = outcome
                print(
                    f"[{utc_stamp()}] END {stage_name}/{task.task_id} "
                    f"status={outcome.get('status')} free_cores={pool.available}",
                    flush=True,
                )

            if now - last_report >= 30:
                stage_report = {
                    **status_context,
                    "updated_utc": utc_stamp(),
                    "current_stage": stage_name,
                    "stage_elapsed_seconds": now - stage_start,
                    "pending_count": len(pending),
                    "running_count": len(running),
                    "completed_count": len(outcomes),
                    "free_physical_cores": pool.available,
                    "running": [task.task_id for task, _ in running.values()],
                    "pending_preview": [task.task_id for task in list(pending)[:20]],
                }
                write_json(run_dir / "pipeline_status.json", stage_report)
                last_report = now
    write_json(run_dir / "stage_status" / f"{stage_name}.json", outcomes)
    return outcomes


def parse_and_validate_shielding(
    *,
    output: Path,
    atom_map_path: Path,
    parsed_path: Path,
    method: dict[str, Any],
    geometry_tier: str,
    solvent: str,
    figure_id: str,
    run_dir: Path,
    run_record: dict[str, Any],
    workflow_fingerprint_value: str,
    run_uuid: str,
) -> list[dict[str, Any]]:
    atom_map = read_csv(atom_map_path)
    expected = [row for row in atom_map if row["element"] in {"H", "C"}]
    values = parse_shielding_summary(output)
    observed_lookup = {(int(row["orca_atom_index_0based"]), row["element"]): row for row in values}
    try:
        source_output_relative = output.resolve().relative_to(run_dir.resolve()).as_posix()
    except ValueError:
        source_output_relative = ""
    rows: list[dict[str, Any]] = []
    for atom in expected:
        index0 = int(atom["atom_index_1based"]) - 1
        key = (index0, atom["element"])
        if key not in observed_lookup:
            raise ValueError(f"Missing shielding for {figure_id} atom {index0 + 1} {atom['element']}")
        value = observed_lookup[key]
        rows.append(
            {
                "figure_id": figure_id,
                "atom_index_1based": index0 + 1,
                "element": atom["element"],
                "isotropic_shielding_ppm": value["isotropic_shielding_ppm"],
                "anisotropy_ppm": value["anisotropy_ppm"],
                "include_in_analyte_nmr": atom["include_in_analyte_nmr"],
                "exchangeable": atom["exchangeable"],
                "reference_signal": atom.get("reference_signal", "false"),
                "method_id": method["method_id"],
                "method_quality_rank": method["quality_rank"],
                "method_quality_label": method["quality_label"],
                "geometry_tier": geometry_tier,
                "implicit_solvent": solvent,
                "source_output": str(output),
                "source_output_relative_to_run_dir": source_output_relative,
                "source_output_sha256": sha256(output),
                "source_input_sha256": run_record.get("input_sha256", ""),
                "source_geometry_sha256": run_record.get("geometry_sha256", ""),
                "workflow_fingerprint": workflow_fingerprint_value,
                "run_uuid": run_uuid,
            }
        )
    if len(rows) != len(expected):
        raise ValueError(f"Shielding count mismatch for {figure_id}: {len(rows)} != {len(expected)}")
    write_csv(parsed_path, rows, list(rows[0].keys()) if rows else ["figure_id"])
    return rows


def composition_cost(row: dict[str, str]) -> tuple[int, int, int]:
    mol = Chem.AddHs(Chem.MolFromSmiles(row["model_smiles"]))
    heavy = mol.GetNumHeavyAtoms()
    rotatable = Lipinski.NumRotatableBonds(Chem.RemoveHs(mol))
    score = int(heavy * heavy * (1.0 + 0.25 * rotatable))
    return score, heavy, rotatable


def hrms_validation_issues(output: Path, expected_ids: set[str]) -> list[str]:
    validation_path = output / "VALIDATION_AND_PROVENANCE.json"
    summary_path = output / "Figure3_formula_HRMS_prediction_summary.csv"
    manifest_path = output / "FILE_MANIFEST_SHA256.csv"
    if not validation_path.is_file():
        return ["validation_missing"]
    audit = read_json(validation_path)
    issues: list[str] = []
    expected_count = len(expected_ids)
    if audit.get("prediction_status") != "CALCULATED HRMS PREDICTION - NOT EXPERIMENTAL":
        issues.append("prediction_label_invalid")
    if audit.get("expected_compound_count") != expected_count:
        issues.append("expected_compound_count_mismatch")
    if audit.get("observed_compound_count") != expected_count:
        issues.append("observed_compound_count_mismatch")
    if audit.get("missing_ids") or audit.get("extra_ids"):
        issues.append("missing_or_extra_ids")
    if audit.get("all_self_tests_pass") is not True:
        issues.append("element_self_tests_failed")
    crosscheck = audit.get("independent_reference_mass_crosscheck", {})
    if crosscheck.get("pass") is not True:
        issues.append("independent_mass_crosscheck_failed")
    if crosscheck.get("comparison_count") != expected_count:
        issues.append("independent_mass_crosscheck_count_mismatch")
    if audit.get("fragment_ions_generated") is not False:
        issues.append("fragment_scope_flag_invalid")
    if audit.get("experimental_peaks_generated") is not False:
        issues.append("experimental_scope_flag_invalid")
    if audit.get("input_sdf_sha256") != sha256(
        ROOT / "registry" / "Figure3_52_ORCA_ordered_reference_states.sdf"
    ):
        issues.append("input_sdf_hash_mismatch")
    if audit.get("input_registry_sha256") != sha256(
        ROOT / "registry" / "Figure3_52_ORCA_reference_state_registry.csv"
    ):
        issues.append("input_registry_hash_mismatch")
    if audit.get("pinned_engine_sha256") != sha256(
        ROOT / "scripts" / "generate_formula_hrms.py"
    ):
        issues.append("pinned_engine_hash_mismatch")
    if not summary_path.is_file():
        issues.append("summary_missing")
    else:
        observed_ids = {row["figure_id"] for row in read_banner_csv(summary_path)}
        if observed_ids != expected_ids:
            issues.append("summary_id_set_mismatch")
    issues.extend(
        f"internal_{item}"
        for item in artifact_manifest_issues(
            manifest_path, workflow_root=ROOT, generated_root=output
        )
    )
    return issues


def run_hrms(run_dir: Path, only: str) -> dict[str, Any]:
    output = run_dir / "HRMS_CALCULATED_NOT_EXPERIMENTAL"
    validation = output / "VALIDATION_AND_PROVENANCE.json"
    selected = {item.strip() for item in only.split(",") if item.strip()}
    expected_ids = selected or {f"F3_{index:03d}" for index in range(1, 53)}
    reuse_issues = hrms_validation_issues(output, expected_ids) if validation.is_file() else ["validation_missing"]
    if not reuse_issues:
        return {"status": "REUSED_VALIDATED", "output": str(output), "compound_ids": sorted(expected_ids)}
    command = [
        sys.executable,
        str(ROOT / "scripts" / "generate_hrms_all.py"),
        "--output-dir",
        str(output),
        "--allowed-root",
        str(run_dir),
        "--overwrite",
    ]
    if only:
        command.extend(["--only", only])
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    record = {
        "status": "COMPLETE" if result.returncode == 0 else "FAILED",
        "return_code": result.returncode,
        "command": command,
        "stdout": result.stdout[-10000:],
        "stderr": result.stderr[-10000:],
        "output": str(output),
    }
    write_json(run_dir / "hrms_status.json", record)
    if result.returncode != 0:
        raise RuntimeError(f"HRMS generation failed: {result.stderr[-2000:]}")
    validation_issues = hrms_validation_issues(output, expected_ids)
    if validation_issues:
        raise RuntimeError(f"HRMS validation failed after generation: {validation_issues}")
    return record


def render_planned_inputs(
    run_dir: Path,
    species: list[dict[str, str]],
    methods: dict[str, Any],
    hardware: dict[str, Any],
) -> None:
    planned = run_dir / "planned_inputs"
    planned.mkdir(parents=True, exist_ok=True)
    geom_method = methods["geometry"]["baseline"]
    primary = method_by_id(methods, "TIER1_PBE0_PCSSEG1")
    refined = method_by_id(methods, "TIER2_PBE0_PCSSEG2")
    for row in species:
        figure_id = row["figure_id"]
        structure = ROOT / "structures" / figure_id / "start.xyz"
        elements = [atom[0] for atom in read_xyz(structure)[1]]
        nprocs = 8 if len(elements) >= int(hardware["per_job"]["large_atom_threshold"]) else 4
        folder = planned / figure_id
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "01_xtb_opt.inp").write_text(
            build_opt_input(
                geom_method["simple_keywords"].format(solvent=row["orca_continuum_keyword"]),
                int(row["charge"]),
                int(row["multiplicity"]),
                "start.xyz",
                2,
                2000,
            ),
            encoding="utf-8",
        )
        iodine = methods["iodine_override"]["basis_block"] if figure_id == "F3_006" else None
        (folder / "02_tier1_shield.inp").write_text(
            build_shield_input(
                primary["simple_keywords"].format(solvent=row["orca_continuum_keyword"]),
                int(row["charge"]), int(row["multiplicity"]), "geometry.xyz", nprocs,
                int(hardware["per_job"]["maxcore_mb_per_rank_baseline"]), elements, iodine,
            ), encoding="utf-8"
        )
        (folder / "03_tier2_shield.inp").write_text(
            build_shield_input(
                refined["simple_keywords"].format(solvent=row["orca_continuum_keyword"]),
                int(row["charge"]), int(row["multiplicity"]), "geometry.xyz", nprocs,
                int(hardware["per_job"]["maxcore_mb_per_rank_refined"]), elements, iodine,
            ), encoding="utf-8"
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orca", type=Path, help="Absolute path to ORCA 6.1.1 driver")
    parser.add_argument("--scratch", type=Path, help="Node-local NVMe scratch")
    parser.add_argument("--run-dir", type=Path, default=ROOT / "results" / "production_24h")
    parser.add_argument("--core-cap", type=int, default=0)
    parser.add_argument("--deadline-hours", type=float, default=0.0)
    parser.add_argument(
        "--budget-start-epoch",
        type=float,
        default=0.0,
        help="Unix epoch when the node workflow started; includes preflight in the hard budget",
    )
    parser.add_argument("--only", default="", help="Comma-separated F3 IDs for pilot/debug")
    parser.add_argument("--disable-goat", action="store_true")
    parser.add_argument("--skip-hrms", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    for sig in (signal.SIGINT, signal.SIGTERM):
        signal.signal(sig, request_stop)
    if hasattr(signal, "SIGUSR1"):
        signal.signal(signal.SIGUSR1, request_stop)

    hardware = read_json(ROOT / "config" / "hardware_104physical_24h.json")
    methods = read_json(ROOT / "config" / "methods_24h.json")
    all_species = read_csv(ROOT / "config" / "primary_species.csv")
    selected_ids = {item.strip() for item in args.only.split(",") if item.strip()}
    species = [row for row in all_species if not selected_ids or row["figure_id"] in selected_ids]
    if selected_ids and selected_ids != {row["figure_id"] for row in species}:
        missing = sorted(selected_ids - {row["figure_id"] for row in species})
        raise ValueError(f"Unknown --only IDs: {missing}")
    args.run_dir.mkdir(parents=True, exist_ok=True)
    if args.dry_run:
        render_planned_inputs(args.run_dir, species, methods, hardware)
        print(f"Rendered planned inputs for {len(species)} compounds: {args.run_dir / 'planned_inputs'}")
        return 0
    if args.orca is None or not args.orca.is_absolute():
        raise ValueError("--orca must be the absolute path to the ORCA driver")
    if not args.orca.is_file():
        raise FileNotFoundError(args.orca)

    scratch = args.scratch
    if scratch is None:
        scratch_env = os.environ.get("SLURM_TMPDIR") or os.environ.get("TMPDIR") or "/tmp"
        scratch = Path(scratch_env) / f"orca52_{os.environ.get('USER', 'user')}"
    scratch.mkdir(parents=True, exist_ok=True)

    cap = args.core_cap or int(hardware["resource_limits"]["physical_core_cap"])
    cpu_ids = physical_cpu_ids(cap)
    if len(cpu_ids) < 4:
        raise RuntimeError(f"At least four physical cores are required; found {len(cpu_ids)}")
    configured_seconds = float(hardware["deadline_policy"]["total_budget_seconds"])
    total_seconds = args.deadline_hours * 3600 if args.deadline_hours else configured_seconds
    elapsed_before_pipeline = (
        max(0.0, time.time() - args.budget_start_epoch) if args.budget_start_epoch else 0.0
    )
    start_clock = time.monotonic() - elapsed_before_pipeline
    hard_deadline = start_clock + total_seconds
    scale = total_seconds / configured_seconds
    geometry_deadline = start_clock + float(hardware["deadline_policy"]["geometry_stage_cutoff_seconds"]) * scale
    reference_deadline = start_clock + float(hardware["deadline_policy"]["reference_stage_cutoff_seconds"]) * scale
    baseline_deadline = start_clock + float(hardware["deadline_policy"]["baseline_stage_cutoff_seconds"]) * scale
    baseline_recovery_deadline = start_clock + float(hardware["deadline_policy"]["baseline_recovery_cutoff_seconds"]) * scale
    refine_launch_cutoff = start_clock + float(hardware["deadline_policy"]["refinement_launch_cutoff_seconds"]) * scale
    optional_deadline = start_clock + float(hardware["deadline_policy"]["terminate_optional_jobs_seconds"]) * scale
    finalization_deadline = start_clock + float(hardware["deadline_policy"]["finalization_target_seconds"]) * scale
    if time.monotonic() >= hard_deadline:
        raise RuntimeError("The workflow budget was already exhausted during preflight")

    fingerprint, fingerprint_inputs = workflow_fingerprint()
    run_config_path = args.run_dir / "run_config.json"
    previous: dict[str, Any] = {}
    if run_config_path.is_file():
        previous = read_json(run_config_path)
        previous_fingerprint = previous.get("workflow_fingerprint")
        if previous_fingerprint and previous_fingerprint != fingerprint:
            raise RuntimeError(
                "This run directory was created with a different workflow fingerprint. "
                "Use a new --run-dir instead of mixing stale and new parsed results."
            )
    run_uuid = str(previous.get("run_uuid") or uuid.uuid4())
    status_context = {
        "schema_version": "1.0",
        "run_start_utc": utc_stamp(),
        "run_dir": str(args.run_dir.resolve()),
        "orca": str(args.orca.resolve()),
        "scratch": str(scratch.resolve()),
        "physical_cpu_ids": cpu_ids,
        "total_budget_seconds": total_seconds,
        "budget_start_epoch": args.budget_start_epoch or time.time(),
        "preflight_elapsed_seconds": elapsed_before_pipeline,
        "compound_count": len(species),
        "implicit_solvent_only": True,
        "prediction_only": True,
        "workflow_fingerprint": fingerprint,
        "workflow_fingerprint_inputs": fingerprint_inputs,
        "run_uuid": run_uuid,
        "runtime_environment_identity": runtime_environment_identity(),
    }
    write_json(run_config_path, status_context)
    raw_root = args.run_dir / "raw_orca"
    parsed_root = args.run_dir / "parsed_shieldings"
    geometry_root = args.run_dir / "geometries"
    raw_root.mkdir(parents=True, exist_ok=True)
    parsed_root.mkdir(parents=True, exist_ok=True)
    geometry_root.mkdir(parents=True, exist_ok=True)

    if not args.skip_hrms:
        print(f"[{utc_stamp()}] generating formula-level HRMS", flush=True)
        run_hrms(args.run_dir, args.only)

    observable_species = [
        row for row in species
        if int(row["observable_1H_count_after_exchange_policy"]) > 0 or int(row["observable_13C_count"]) > 0
    ]

    # Stage 1: baseline geometries for analytes and both references.
    geometry_tasks: list[Task] = []
    baseline_geom_method = methods["geometry"]["baseline"]
    xtb_timeout = float(hardware["stage_timeouts_seconds"]["xtb_opt"])

    def make_xtb_task(row: dict[str, str]) -> Task:
        figure_id = row["figure_id"]
        start_xyz = ROOT / "structures" / figure_id / "start.xyz"
        destination = geometry_root / "baseline" / f"{figure_id}.xyz"

        def runner(cores: list[int]) -> dict[str, Any]:
            input_text = build_opt_input(
                baseline_geom_method["simple_keywords"].format(solvent=row["orca_continuum_keyword"]),
                int(row["charge"]), int(row["multiplicity"]), start_xyz.name,
                len(cores), 1800,
            )
            timeout = select_job_timeout(xtb_timeout, geometry_deadline)
            record, attempt = run_or_reuse(
                job_id=f"{figure_id}_xtb_opt", input_text=input_text, geometry=start_xyz,
                orca=args.orca, scratch=scratch, raw_root=raw_root, cpu_ids=cores,
                timeout=timeout, basename="xtb_opt",
            )
            source = find_final_geometry(attempt, "xtb_opt") if record.get("status") == "COMPLETE" else None
            integrity_issues: list[str] = []
            if source is not None:
                integrity_issues = frozen_geometry_issues(
                    start_xyz,
                    source,
                    ROOT / "structures" / figure_id / "bonds.csv",
                    row["model_smiles"],
                )
            fallback = False
            if source is None or integrity_issues:
                source = start_xyz
                fallback = True
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            result = {
                "task_id": figure_id,
                "status": "BASELINE_GEOMETRY_START_FALLBACK" if fallback else "COMPLETE",
                "geometry": str(destination),
                "geometry_sha256": sha256(destination),
                "integrity_issues": integrity_issues,
                "orca_run": record,
            }
            write_json(destination.with_suffix(".json"), result)
            return result

        return Task(figure_id, min(2, len(cpu_ids)), "xTB/ALPB baseline optimization", runner)

    references = read_csv(ROOT / "references" / "reference_registry.csv")
    for reference in references:
        ref_id = reference["reference_id"]
        start_xyz = ROOT / "references" / ref_id / "start.xyz"
        destination = geometry_root / "references" / "baseline" / f"{ref_id}.xyz"

        def make_ref_xtb(ref: dict[str, str], src: Path, dst: Path) -> Task:
            def runner(cores: list[int]) -> dict[str, Any]:
                input_text = build_opt_input(
                    baseline_geom_method["simple_keywords"].format(solvent=ref["orca_continuum_keyword"]),
                    int(ref["charge"]), int(ref["multiplicity"]), src.name, len(cores), 1800,
                )
                timeout = select_job_timeout(xtb_timeout, geometry_deadline)
                record, attempt = run_or_reuse(
                    job_id=f"REF_{ref['reference_id']}_xtb_opt", input_text=input_text, geometry=src,
                    orca=args.orca, scratch=scratch, raw_root=raw_root, cpu_ids=cores,
                    timeout=timeout, basename="xtb_opt",
                )
                final = find_final_geometry(attempt, "xtb_opt") if record.get("status") == "COMPLETE" else None
                fallback = final is None or not geometry_same_composition(src, final)
                final = src if fallback else final
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(final, dst)
                return {
                    "task_id": ref["reference_id"],
                    "status": "REFERENCE_START_FALLBACK" if fallback else "COMPLETE",
                    "geometry": str(dst),
                    "geometry_sha256": sha256(dst),
                    "orca_run": record,
                }
            return Task(ref["reference_id"], min(2, len(cpu_ids)), "reference xTB/ALPB optimization", runner)
        geometry_tasks.append(make_ref_xtb(reference, start_xyz, destination))
    # Reference geometries are deliberately first in the queue because every
    # reported chemical shift depends on them.  Analytes fill the remaining
    # cores immediately and retain the same four-hour full-coverage fallback.
    geometry_tasks.extend(make_xtb_task(row) for row in observable_species)

    geometry_outcomes = run_task_batch(
        stage_name="01_baseline_geometry", tasks=geometry_tasks, cpu_ids=cpu_ids,
        run_dir=args.run_dir, stop_launch_at=geometry_deadline, status_context=status_context,
    )
    # Full-coverage invariant: every analyte/reference gets at least its frozen,
    # audited start geometry even if the xTB stage reached its four-hour gate.
    for row in observable_species:
        destination = geometry_root / "baseline" / f"{row['figure_id']}.xyz"
        if not destination.is_file():
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / "structures" / row["figure_id"] / "start.xyz", destination)
            write_json(
                destination.with_suffix(".json"),
                {
                    "task_id": row["figure_id"],
                    "status": "BASELINE_GEOMETRY_START_FALLBACK_AFTER_STAGE_GATE",
                    "geometry": str(destination),
                    "geometry_sha256": sha256(destination),
                    "stage_outcome": geometry_outcomes.get(row["figure_id"], {}),
                },
            )
    for ref in references:
        destination = geometry_root / "references" / "baseline" / f"{ref['reference_id']}.xyz"
        if not destination.is_file():
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / "references" / ref["reference_id"] / "start.xyz", destination)

    # Stage 2: baseline references are a hard prerequisite for every reported
    # chemical shift.  Run them before any optional Tier 2 reference work so a
    # slow r2SCAN optimization can never consume their six-hour gate.
    reference_atom_maps = {
        ref["reference_id"]: ROOT / "references" / ref["reference_id"] / "atom_map.csv"
        for ref in references
    }

    def make_reference_shield(
        ref: dict[str, str],
        meth: dict[str, Any],
        geom: Path,
        parsed_path: Path,
        geom_tier: str,
        stage_deadline: float,
    ) -> Task:
        def runner(cores: list[int]) -> dict[str, Any]:
            elements = [atom[0] for atom in read_xyz(geom)[1]]
            input_text = build_shield_input(
                meth["simple_keywords"].format(solvent=ref["orca_continuum_keyword"]),
                int(ref["charge"]), int(ref["multiplicity"]), geom.name, len(cores),
                int(
                    hardware["per_job"][
                        "maxcore_mb_per_rank_refined"
                        if geom_tier == "REFINED_GEOMETRY"
                        else "maxcore_mb_per_rank_baseline"
                    ]
                ),
                elements,
            )
            timeout = select_job_timeout(
                float(hardware["stage_timeouts_seconds"]["reference_job"]),
                stage_deadline,
            )
            job_id = f"REF_{ref['reference_id']}_shield_{meth['method_id']}"
            record, attempt = run_or_reuse(
                job_id=job_id, input_text=input_text, geometry=geom, orca=args.orca,
                scratch=scratch, raw_root=raw_root, cpu_ids=cores, timeout=timeout,
                basename="shield",
            )
            if record.get("status") != "COMPLETE":
                return {"task_id": job_id, "status": record.get("status"), "orca_run": record}
            output = attempt / "shield.out"
            rows = parse_and_validate_shielding(
                output=output, atom_map_path=reference_atom_maps[ref["reference_id"]],
                parsed_path=parsed_path, method=meth, geometry_tier=geom_tier,
                solvent=ref["orca_continuum_keyword"], figure_id=ref["reference_id"],
                run_dir=args.run_dir, run_record=record,
                workflow_fingerprint_value=fingerprint, run_uuid=run_uuid,
            )
            return {
                "task_id": job_id, "status": "COMPLETE", "parsed": str(parsed_path),
                "nuclei": len(rows), "orca_run": record,
            }

        return Task(
            f"{ref['reference_id']}_{meth['method_id']}",
            min(4, len(cpu_ids)),
            f"matching {meth['method_id']} computational reference",
            runner,
        )

    baseline_reference_methods = (
        "TIER1_PBE0_PCSSEG1", "FALLBACK_PBE_PCSSEG1", "EMERGENCY_PBE_DEF2SVP"
    )
    baseline_reference_tasks: list[Task] = []
    for reference in references:
        ref_id = reference["reference_id"]
        geometry = geometry_root / "references" / "baseline" / f"{ref_id}.xyz"
        for method_id in baseline_reference_methods:
            method = method_by_id(methods, method_id)
            parsed = parsed_root / "references" / ref_id / f"{method_id}.csv"
            baseline_reference_tasks.append(
                make_reference_shield(
                    reference, method, geometry, parsed, "BASELINE_GEOMETRY",
                    reference_deadline,
                )
            )
    run_task_batch(
        stage_name="02_baseline_reference_shieldings",
        tasks=baseline_reference_tasks,
        cpu_ids=cpu_ids,
        run_dir=args.run_dir,
        stop_launch_at=reference_deadline,
        status_context=status_context,
    )
    missing_reference_solvents = []
    for ref in references:
        available = [
            method_id for method_id in baseline_reference_methods
            if (parsed_root / "references" / ref["reference_id"] / f"{method_id}.csv").is_file()
        ]
        if not available:
            missing_reference_solvents.append(ref["reference_id"])
    if missing_reference_solvents:
        raise RuntimeError(
            "No matching baseline reference method finished for: "
            + ", ".join(missing_reference_solvents)
        )

    # Stage 3: optional r2SCAN reference geometries for Tier 2 upgrades.
    refined_geom_method = methods["geometry"]["refined"]
    ref_opt_tasks: list[Task] = []
    for reference in references:
        ref_id = reference["reference_id"]
        source = geometry_root / "references" / "baseline" / f"{ref_id}.xyz"
        destination = geometry_root / "references" / "refined" / f"{ref_id}.xyz"

        def make_ref_r2(ref: dict[str, str], src: Path, dst: Path) -> Task:
            def runner(cores: list[int]) -> dict[str, Any]:
                input_text = build_opt_input(
                    refined_geom_method["simple_keywords"].format(solvent=ref["orca_continuum_keyword"]),
                    int(ref["charge"]), int(ref["multiplicity"]), src.name, len(cores),
                    int(hardware["per_job"]["maxcore_mb_per_rank_refined"]),
                )
                timeout = select_job_timeout(float(hardware["stage_timeouts_seconds"]["r2scan_opt"]), reference_deadline)
                record, attempt = run_or_reuse(
                    job_id=f"REF_{ref['reference_id']}_r2scan_opt", input_text=input_text, geometry=src,
                    orca=args.orca, scratch=scratch, raw_root=raw_root, cpu_ids=cores,
                    timeout=timeout, basename="r2scan_opt",
                )
                final = find_final_geometry(attempt, "r2scan_opt") if record.get("status") == "COMPLETE" else None
                fallback = final is None or not geometry_same_composition(src, final)
                if fallback:
                    return {
                        "task_id": ref["reference_id"],
                        "status": "REFERENCE_REFINED_GEOMETRY_NOT_COMPLETED",
                        "geometry": "",
                        "orca_run": record,
                    }
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(final, dst)
                return {
                    "task_id": ref["reference_id"], "status": "COMPLETE",
                    "geometry": str(dst), "geometry_sha256": sha256(dst), "orca_run": record,
                }
            return Task(
                f"{ref['reference_id']}_r2scan", min(4, len(cpu_ids)),
                "reference r2SCAN-3c/CPCM optimization", runner,
            )
        ref_opt_tasks.append(make_ref_r2(reference, source, destination))
    run_task_batch(
        stage_name="03_reference_refined_geometry", tasks=ref_opt_tasks, cpu_ids=cpu_ids,
        run_dir=args.run_dir, stop_launch_at=reference_deadline, status_context=status_context,
    )

    # Stage 4: optional Tier 2 references. Baseline production does not wait for
    # these; if the reference gate has closed, sample refinements are skipped.
    reference_tasks: list[Task] = []
    for reference in references:
        ref_id = reference["reference_id"]
        method = method_by_id(methods, "TIER2_PBE0_PCSSEG2")
        geometry = geometry_root / "references" / "refined" / f"{ref_id}.xyz"
        if not geometry.is_file():
            continue
        parsed = parsed_root / "references" / ref_id / f"{method['method_id']}.csv"
        reference_tasks.append(
            make_reference_shield(
                reference, method, geometry, parsed, "REFINED_GEOMETRY",
                reference_deadline,
            )
        )
    run_task_batch(
        stage_name="04_optional_tier2_reference_shieldings", tasks=reference_tasks, cpu_ids=cpu_ids,
        run_dir=args.run_dir, stop_launch_at=reference_deadline, status_context=status_context,
    )

    # Stage 4: complete shielding coverage, with explicit method fallbacks.
    baseline_tasks: list[Task] = []
    method_sequence = [
        method_by_id(methods, "TIER1_PBE0_PCSSEG1"),
        method_by_id(methods, "FALLBACK_PBE_PCSSEG1"),
        method_by_id(methods, "EMERGENCY_PBE_DEF2SVP"),
    ]
    large_threshold = int(hardware["per_job"]["large_atom_threshold"])

    def make_baseline_shield(row: dict[str, str]) -> Task:
        figure_id = row["figure_id"]
        geometry = geometry_root / "baseline" / f"{figure_id}.xyz"
        atom_map = ROOT / "structures" / figure_id / "atom_map.csv"
        atom_count = len(read_xyz(geometry)[1])
        ncores = min(8 if atom_count >= large_threshold else 4, len(cpu_ids))

        def runner(cores: list[int]) -> dict[str, Any]:
            elements = [atom[0] for atom in read_xyz(geometry)[1]]
            attempts: list[dict[str, Any]] = []
            for index, method in enumerate(method_sequence):
                ref_id = "TMS_DMSO" if row["orca_continuum_keyword"] == "DMSO" else "TSP_MINUS_WATER"
                matching_ref = parsed_root / "references" / ref_id / f"{method['method_id']}.csv"
                if not matching_ref.is_file():
                    attempts.append({"method_id": method["method_id"], "status": "SKIPPED_REFERENCE_MISSING"})
                    continue
                iodine = methods["iodine_override"]["basis_block"] if figure_id == "F3_006" and "PCSSEG" in method["method_id"] else None
                input_text = build_shield_input(
                    method["simple_keywords"].format(solvent=row["orca_continuum_keyword"]),
                    int(row["charge"]), int(row["multiplicity"]), geometry.name, len(cores),
                    int(hardware["per_job"]["maxcore_mb_per_rank_baseline"]), elements,
                    iodine_override=iodine, slowconv=index > 0,
                )
                stage_limit = float(
                    hardware["stage_timeouts_seconds"]["baseline_shield" if index == 0 else "baseline_retry"]
                )
                timeout = select_job_timeout(stage_limit, baseline_deadline)
                job_id = f"{figure_id}_shield_{method['method_id']}"
                record, attempt = run_or_reuse(
                    job_id=job_id, input_text=input_text, geometry=geometry, orca=args.orca,
                    scratch=scratch, raw_root=raw_root, cpu_ids=cores, timeout=timeout,
                    basename="shield",
                )
                attempts.append({"method_id": method["method_id"], "orca_run": record})
                if record.get("status") != "COMPLETE":
                    continue
                parsed = parsed_root / "compounds" / figure_id / f"{method['method_id']}.csv"
                try:
                    rows = parse_and_validate_shielding(
                        output=attempt / "shield.out", atom_map_path=atom_map, parsed_path=parsed,
                        method=method, geometry_tier="BASELINE_GEOMETRY",
                        solvent=row["orca_continuum_keyword"], figure_id=figure_id,
                        run_dir=args.run_dir, run_record=record,
                        workflow_fingerprint_value=fingerprint, run_uuid=run_uuid,
                    )
                except Exception as exc:
                    attempts[-1]["parse_error"] = f"{type(exc).__name__}: {exc}"
                    continue
                return {
                    "task_id": figure_id,
                    "status": "COMPLETE" if method["method_id"] == "TIER1_PBE0_PCSSEG1" else "COMPLETE_WITH_FALLBACK",
                    "selected_method_id": method["method_id"],
                    "selected_method_quality_rank": method["quality_rank"],
                    "parsed": str(parsed),
                    "nuclei": len(rows),
                    "attempts": attempts,
                }
            return {"task_id": figure_id, "status": "ALL_BASELINE_METHODS_FAILED", "attempts": attempts}

        return Task(figure_id, ncores, "PBE0/pcSseg-1 shielding with audited fallbacks", runner)

    scored_species = sorted(observable_species, key=lambda row: composition_cost(row)[0], reverse=True)
    baseline_tasks.extend(make_baseline_shield(row) for row in scored_species)
    baseline_outcomes = run_task_batch(
        stage_name="04_baseline_shieldings", tasks=baseline_tasks, cpu_ids=cpu_ids,
        run_dir=args.run_dir, stop_launch_at=baseline_deadline, status_context=status_context,
    )
    baseline_missing = sorted(
        row["figure_id"] for row in observable_species
        if baseline_outcomes.get(row["figure_id"], {}).get("status") not in {"COMPLETE", "COMPLETE_WITH_FALLBACK"}
    )
    if baseline_missing:
        emergency_method = method_by_id(methods, "EMERGENCY_PBE_DEF2SVP")
        recovery_tasks: list[Task] = []

        def make_recovery(row: dict[str, str]) -> Task:
            figure_id = row["figure_id"]
            geometry = geometry_root / "baseline" / f"{figure_id}.xyz"
            atom_map = ROOT / "structures" / figure_id / "atom_map.csv"
            desired_cores = 8 if len(read_xyz(geometry)[1]) >= large_threshold else 4

            def runner(cores: list[int]) -> dict[str, Any]:
                ref_id = "TMS_DMSO" if row["orca_continuum_keyword"] == "DMSO" else "TSP_MINUS_WATER"
                matching_ref = parsed_root / "references" / ref_id / f"{emergency_method['method_id']}.csv"
                if not matching_ref.is_file():
                    return {
                        "task_id": figure_id,
                        "status": "RECOVERY_REFERENCE_MISSING",
                        "reference": str(matching_ref),
                    }
                elements = [atom[0] for atom in read_xyz(geometry)[1]]
                input_text = build_shield_input(
                    emergency_method["simple_keywords"].format(solvent=row["orca_continuum_keyword"]),
                    int(row["charge"]), int(row["multiplicity"]), geometry.name, len(cores),
                    int(hardware["per_job"]["maxcore_mb_per_rank_baseline"]), elements,
                    slowconv=True,
                )
                timeout = select_job_timeout(
                    float(hardware["stage_timeouts_seconds"]["baseline_retry"]),
                    baseline_recovery_deadline,
                )
                record, attempt = run_or_reuse(
                    job_id=f"{figure_id}_shield_EMERGENCY_RECOVERY",
                    input_text=input_text,
                    geometry=geometry,
                    orca=args.orca,
                    scratch=scratch,
                    raw_root=raw_root,
                    cpu_ids=cores,
                    timeout=timeout,
                    basename="shield",
                )
                if record.get("status") != "COMPLETE":
                    return {"task_id": figure_id, "status": "RECOVERY_FAILED", "orca_run": record}
                parsed = parsed_root / "compounds" / figure_id / f"{emergency_method['method_id']}.csv"
                rows = parse_and_validate_shielding(
                    output=attempt / "shield.out",
                    atom_map_path=atom_map,
                    parsed_path=parsed,
                    method=emergency_method,
                    geometry_tier="BASELINE_GEOMETRY",
                    solvent=row["orca_continuum_keyword"],
                    figure_id=figure_id,
                    run_dir=args.run_dir,
                    run_record=record,
                    workflow_fingerprint_value=fingerprint,
                    run_uuid=run_uuid,
                )
                return {
                    "task_id": figure_id,
                    "status": "COMPLETE_WITH_FALLBACK",
                    "selected_method_id": emergency_method["method_id"],
                    "parsed": str(parsed),
                    "nuclei": len(rows),
                    "orca_run": record,
                }

            return Task(
                figure_id,
                min(desired_cores, len(cpu_ids)),
                "deadline recovery with PBE/def2-SVP shielding",
                runner,
            )

        by_id = {row["figure_id"]: row for row in observable_species}
        recovery_tasks.extend(make_recovery(by_id[figure_id]) for figure_id in baseline_missing)
        recovery_outcomes = run_task_batch(
            stage_name="04b_baseline_recovery",
            tasks=recovery_tasks,
            cpu_ids=cpu_ids,
            run_dir=args.run_dir,
            stop_launch_at=baseline_recovery_deadline,
            status_context=status_context,
        )
        baseline_outcomes.update(recovery_outcomes)
        baseline_missing = sorted(
            figure_id for figure_id in baseline_missing
            if recovery_outcomes.get(figure_id, {}).get("status") != "COMPLETE_WITH_FALLBACK"
        )
    write_json(args.run_dir / "baseline_coverage.json", {
        "observable_compound_count": len(observable_species),
        "complete_count": len(observable_species) - len(baseline_missing),
        "missing": baseline_missing,
        "pass": not baseline_missing,
    })

    # Stage 5: opportunistic GOAT + r2SCAN geometry + Tier 2 shielding.
    refined_tasks: list[Task] = []
    goat_method = methods["geometry"]["goat"]
    refined_method = method_by_id(methods, "TIER2_PBE0_PCSSEG2")

    def make_refinement(row: dict[str, str]) -> Task:
        figure_id = row["figure_id"]
        baseline_geometry = geometry_root / "baseline" / f"{figure_id}.xyz"
        atom_map = ROOT / "structures" / figure_id / "atom_map.csv"
        score, heavy, rotatable = composition_cost(row)
        atom_count = len(read_xyz(baseline_geometry)[1])
        ncores = min(8 if atom_count >= large_threshold else 4, len(cpu_ids))

        def runner(cores: list[int]) -> dict[str, Any]:
            history: list[dict[str, Any]] = []
            current = baseline_geometry
            current_tier = "BASELINE_GEOMETRY"
            frozen_start = ROOT / "structures" / figure_id / "start.xyz"
            frozen_bonds = ROOT / "structures" / figure_id / "bonds.csv"
            if not args.disable_goat and rotatable >= 2:
                goat_options = dict(goat_method)
                goat_options["nworkers"] = 8 if len(cores) >= 8 else 4
                goat_input = build_opt_input(
                    goat_method["simple_keywords"].format(solvent=row["orca_continuum_keyword"]),
                    int(row["charge"]), int(row["multiplicity"]), current.name, len(cores), 1800,
                    goat=goat_options,
                )
                timeout = select_job_timeout(float(hardware["stage_timeouts_seconds"]["goat"]), optional_deadline)
                goat_record, goat_attempt = run_or_reuse(
                    job_id=f"{figure_id}_goat", input_text=goat_input, geometry=current,
                    orca=args.orca, scratch=scratch, raw_root=raw_root, cpu_ids=cores,
                    timeout=timeout, basename="goat",
                )
                history.append({"stage": "GOAT", "run": goat_record})
                candidate = goat_attempt / "goat.globalminimum.xyz"
                goat_issues = (
                    frozen_geometry_issues(frozen_start, candidate, frozen_bonds, row["model_smiles"])
                    if candidate.is_file() else ["goat_global_minimum_missing"]
                )
                history[-1]["geometry_integrity_issues"] = goat_issues
                if goat_record.get("status") == "COMPLETE" and candidate.is_file() and not goat_issues:
                    goat_geometry = geometry_root / "goat" / f"{figure_id}.xyz"
                    goat_geometry.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(candidate, goat_geometry)
                    current = goat_geometry
                    current_tier = "GOAT_GLOBAL_MINIMUM"

            r2_input = build_opt_input(
                refined_geom_method["simple_keywords"].format(solvent=row["orca_continuum_keyword"]),
                int(row["charge"]), int(row["multiplicity"]), current.name, len(cores),
                int(hardware["per_job"]["maxcore_mb_per_rank_refined"]),
            )
            timeout = select_job_timeout(float(hardware["stage_timeouts_seconds"]["r2scan_opt"]), optional_deadline)
            r2_record, r2_attempt = run_or_reuse(
                job_id=f"{figure_id}_r2scan_opt", input_text=r2_input, geometry=current,
                orca=args.orca, scratch=scratch, raw_root=raw_root, cpu_ids=cores,
                timeout=timeout, basename="r2scan_opt",
            )
            history.append({"stage": "R2SCAN_OPT", "run": r2_record})
            candidate = find_final_geometry(r2_attempt, "r2scan_opt") if r2_record.get("status") == "COMPLETE" else None
            r2_issues = (
                frozen_geometry_issues(frozen_start, candidate, frozen_bonds, row["model_smiles"])
                if candidate is not None else ["r2scan_final_geometry_missing"]
            )
            history[-1]["geometry_integrity_issues"] = r2_issues
            if candidate is not None and not r2_issues:
                refined_geometry = geometry_root / "refined" / f"{figure_id}.xyz"
                refined_geometry.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(candidate, refined_geometry)
                current = refined_geometry
                current_tier = "REFINED_GEOMETRY"
            else:
                return {
                    "task_id": figure_id,
                    "status": "REFINEMENT_GEOMETRY_NOT_COMPLETED",
                    "geometry_tier": current_tier,
                    "history": history,
                }

            elements = [atom[0] for atom in read_xyz(current)[1]]
            iodine = methods["iodine_override"]["basis_block"] if figure_id == "F3_006" else None
            shield_input = build_shield_input(
                refined_method["simple_keywords"].format(solvent=row["orca_continuum_keyword"]),
                int(row["charge"]), int(row["multiplicity"]), current.name, len(cores),
                int(hardware["per_job"]["maxcore_mb_per_rank_refined"]), elements, iodine,
            )
            timeout = select_job_timeout(float(hardware["stage_timeouts_seconds"]["refined_shield"]), optional_deadline)
            shield_record, shield_attempt = run_or_reuse(
                job_id=f"{figure_id}_shield_{refined_method['method_id']}", input_text=shield_input,
                geometry=current, orca=args.orca, scratch=scratch, raw_root=raw_root,
                cpu_ids=cores, timeout=timeout, basename="shield",
            )
            history.append({"stage": "TIER2_SHIELD", "run": shield_record})
            if shield_record.get("status") != "COMPLETE":
                return {
                    "task_id": figure_id, "status": "REFINEMENT_NOT_COMPLETED",
                    "geometry_tier": current_tier, "history": history,
                }
            parsed = parsed_root / "compounds" / figure_id / f"{refined_method['method_id']}.csv"
            try:
                rows = parse_and_validate_shielding(
                    output=shield_attempt / "shield.out", atom_map_path=atom_map, parsed_path=parsed,
                    method=refined_method, geometry_tier=current_tier,
                    solvent=row["orca_continuum_keyword"], figure_id=figure_id,
                    run_dir=args.run_dir, run_record=shield_record,
                    workflow_fingerprint_value=fingerprint, run_uuid=run_uuid,
                )
            except Exception as exc:
                return {
                    "task_id": figure_id, "status": "REFINEMENT_PARSE_FAILED",
                    "error": f"{type(exc).__name__}: {exc}", "history": history,
                }
            return {
                "task_id": figure_id, "status": "COMPLETE", "parsed": str(parsed),
                "nuclei": len(rows), "geometry_tier": current_tier,
                "rotatable_bonds": rotatable, "priority_score": score, "history": history,
            }

        return Task(figure_id, ncores, f"GOAT={rotatable >= 2 and not args.disable_goat}; r2SCAN-3c; PBE0/pcSseg-2", runner)

    refinable = [
        row for row in scored_species
        if row["figure_id"] not in baseline_missing
        and (
            parsed_root
            / "references"
            / ("TMS_DMSO" if row["orca_continuum_keyword"] == "DMSO" else "TSP_MINUS_WATER")
            / "TIER2_PBE0_PCSSEG2.csv"
        ).is_file()
    ]
    refined_tasks.extend(make_refinement(row) for row in refinable)
    refinement_outcomes = run_task_batch(
        stage_name="05_refinements", tasks=refined_tasks, cpu_ids=cpu_ids,
        run_dir=args.run_dir, stop_launch_at=refine_launch_cutoff,
        status_context=status_context,
    )

    # Final post-processing and audit are separate scripts so they can be rerun.
    # They are bounded by the same node-level clock; optional ORCA work can
    # therefore never consume the time reserved for packaging and audit.
    postprocess_command = [
        sys.executable, str(ROOT / "scripts" / "postprocess_nmr.py"),
        "--run-dir", str(args.run_dir),
    ]
    if args.only:
        postprocess_command.extend(["--only", args.only])
    print(f"[{utc_stamp()}] post-processing NMR spectra", flush=True)
    post_timeout = min(900.0, max(1.0, finalization_deadline - time.monotonic() - 240.0))
    try:
        post = subprocess.run(
            postprocess_command, cwd=ROOT, text=True, capture_output=True,
            timeout=post_timeout,
        )
    except subprocess.TimeoutExpired as exc:
        post = deadline_completed_process(postprocess_command, exc, "POSTPROCESS_DEADLINE_EXCEEDED")
    (args.run_dir / "postprocess.stdout.txt").write_text(post.stdout, encoding="utf-8")
    (args.run_dir / "postprocess.stderr.txt").write_text(post.stderr, encoding="utf-8")
    if post.returncode != 0:
        print(post.stdout)
        print(post.stderr, file=sys.stderr)

    audit_command = [
        sys.executable, str(ROOT / "scripts" / "audit_outputs.py"),
        "--run-dir", str(args.run_dir),
    ]
    if args.only:
        audit_command.extend(["--only", args.only])
    audit_timeout = min(240.0, max(1.0, finalization_deadline - time.monotonic() - 45.0))
    try:
        audit = subprocess.run(
            audit_command, cwd=ROOT, text=True, capture_output=True,
            timeout=audit_timeout,
        )
    except subprocess.TimeoutExpired as exc:
        audit = deadline_completed_process(audit_command, exc, "AUDIT_DEADLINE_EXCEEDED")
    (args.run_dir / "audit.stdout.txt").write_text(audit.stdout, encoding="utf-8")
    (args.run_dir / "audit.stderr.txt").write_text(audit.stderr, encoding="utf-8")

    final = {
        **status_context,
        "end_utc": utc_stamp(),
        "elapsed_seconds": time.monotonic() - start_clock,
        "hard_deadline_seconds": total_seconds,
        "baseline_missing": baseline_missing,
        "refined_complete_count": sum(item.get("status") == "COMPLETE" for item in refinement_outcomes.values()),
        "postprocess_return_code": post.returncode,
        "audit_return_code": audit.returncode,
        "manifest_generation_status": "PENDING_FINALIZER",
        "status": "COMPLETE" if post.returncode == 0 and audit.returncode == 0 else "INCOMPLETE",
    }
    write_json(args.run_dir / "pipeline_final.json", final)

    finalizer_command = [
        sys.executable, str(ROOT / "scripts" / "finalize_manifest.py"),
        "--run-dir", str(args.run_dir),
    ]
    finalizer_timeout = max(1.0, hard_deadline - time.monotonic() - 10.0)
    try:
        finalizer = subprocess.run(
            finalizer_command, cwd=ROOT, text=True, capture_output=True,
            timeout=finalizer_timeout,
        )
    except subprocess.TimeoutExpired as exc:
        finalizer = deadline_completed_process(finalizer_command, exc, "FINAL_MANIFEST_DEADLINE_EXCEEDED")
    # Do not write into run_dir after a successful finalizer: that would make
    # the just-created SHA256 manifest stale.  Slurm captures these messages.
    print(finalizer.stdout, end="", flush=True)
    if finalizer.returncode != 0:
        print(finalizer.stderr, file=sys.stderr, flush=True)
    print(json.dumps(final, ensure_ascii=False, indent=2), flush=True)
    return 0 if final["status"] == "COMPLETE" and finalizer.returncode == 0 else 3


if __name__ == "__main__":
    raise SystemExit(main())
