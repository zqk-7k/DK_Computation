#!/usr/bin/env python3
"""Fail-closed adaptive GOAT wave controller.

The controller releases HELD candidates only.  It never invokes scheduler cancellation and never
changes a RUNNING/COMPLETING job.  A later wave is eligible only after every
previously released candidate has reached a terminal state and an assessment
covering all those released jobs has been written.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from package_io import ACTIVE_STATES, KNOWN_STATES, TERMINAL_STATES, atomic_write_json, atomic_write_text, load_json


TRUE_VALUES = {"1", "true", "yes", "y", "accepted", "pass", "passed"}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def first(row: dict[str, str], *names: str, default: str = "") -> str:
    for name in names:
        value = row.get(name, "")
        if str(value).strip():
            return str(value).strip()
    return default


def normalise_row(row: dict[str, str], root: Path) -> dict[str, object]:
    compound = first(row, "compound_id", "figure_id")
    scenario = first(row, "scenario_id", "scenario", "scenario_name")
    job_id = first(row, "job_id", "candidate_id", "start_id")
    if not (compound and scenario and job_id):
        raise ValueError("manifest rows require compound_id/figure_id, scenario_id/scenario, and job_id/candidate_id")
    wave = int(first(row, "wave", "release_wave", default="1"))
    target = int(first(row, "target_success", "success_target", default="0"))
    if target <= 0:
        raise ValueError(f"{job_id}: target_success must be positive")
    raw_dir = first(row, "job_dir", "job_directory", "portable_job_directory", "relative_job_directory", "candidate_dir")
    if not raw_dir:
        raise ValueError(f"{job_id}: missing job_dir")
    job_dir = Path(raw_dir)
    if not job_dir.is_absolute():
        job_dir = (root / job_dir).resolve()
    state = first(row, "state", default="HELD").upper()
    if state in {"RELEASE_READY", "READY", "HOLD"}:
        state = "HELD"
    state_file = job_dir / "job_state.json"
    if state_file.is_file():
        disk = load_json(state_file)
        state = str(disk.get("state", state)).upper()
    if state not in KNOWN_STATES:
        raise ValueError(f"{job_id}: unknown state {state}")
    accepted = first(row, "accepted", default="").lower() in TRUE_VALUES
    if state_file.is_file():
        accepted = bool(load_json(state_file).get("accepted", accepted))
    return {
        "raw": row,
        "compound_id": compound,
        "scenario_id": scenario,
        "job_id": job_id,
        "wave": wave,
        "target_success": target,
        "job_dir": job_dir,
        "state": state,
        "accepted": accepted,
        "submit_script": first(row, "submit_script", "scheduler_script"),
    }


def load_manifest(path: Path, root: Path) -> tuple[list[dict[str, str]], list[dict[str, object]], list[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("manifest has no header")
        raw = list(reader)
        fields = list(reader.fieldnames)
    jobs = [normalise_row(row, root) for row in raw]
    if len({str(j["job_id"]) for j in jobs}) != len(jobs):
        raise ValueError("job_id values must be globally unique")
    return raw, jobs, fields


def manifest_digest(jobs: Iterable[dict[str, object]]) -> str:
    stable = [{k: str(j[k]) for k in ("compound_id", "scenario_id", "job_id", "wave", "state", "accepted")} for j in jobs]
    return hashlib.sha256(json.dumps(stable, sort_keys=True).encode("utf-8")).hexdigest()


def gate_path(gate_root: Path, compound: str, scenario: str) -> Path:
    safe = lambda x: re.sub(r"[^A-Za-z0-9_.-]+", "_", x)
    return gate_root / safe(compound) / f"{safe(scenario)}.json"


def evaluate_gate(gate: dict[str, object], jobs: list[dict[str, object]], target: int,
                  limits: dict[str, float] | None = None) -> tuple[bool, list[str]]:
    reasons: list[str] = []
    released = {str(j["job_id"]) for j in jobs if j["state"] != "HELD"}
    covered = {str(x) for x in gate.get("evaluated_job_ids", [])}
    if not released.issubset(covered):
        reasons.append("assessment does not cover every released job")
    accepted_completed = sum(j["state"] == "COMPLETED" and bool(j["accepted"]) for j in jobs)
    gate_count = int(gate.get("accepted_success_count", -1))
    if gate_count != accepted_completed:
        reasons.append(f"accepted_success_count mismatch: gate={gate_count}, states={accepted_completed}")
    if accepted_completed < target:
        reasons.append(f"accepted successes {accepted_completed} < target {target}")
    for name in ("identity_gate_pass", "composition_gate_pass", "stereochemistry_gate_pass"):
        if gate.get(name) is not True:
            reasons.append(f"{name} is not true")
    if gate.get("low_energy_basin_rediscovered") is not True:
        reasons.append("low-energy basin rediscovery gate not passed")
    if gate.get("motif_coverage_pass") is not True:
        reasons.append("motif coverage gate not passed")
    stability = gate.get("shift_stability", {})
    if not isinstance(stability, dict) or stability.get("evaluated") is not True:
        reasons.append("shift stability not evaluated")
    else:
        gate_limits = limits or {"h_nonexchangeable_max_delta_ppm": 0.05, "c_max_delta_ppm": 0.5, "nh_max_delta_ppm": 0.2}
        for name, limit in gate_limits.items():
            value = stability.get(name)
            if value is None or float(value) > limit:
                reasons.append(f"{name} is absent or exceeds {limit} ppm")
        if stability.get("peak_order_stable") is not True:
            reasons.append("peak ordering is not stable")
    return not reasons, reasons


def gate_is_current(gate: dict[str, object], jobs: list[dict[str, object]]) -> tuple[bool, list[str]]:
    """Separate record freshness from a scientifically negative assessment."""
    reasons = []
    released = {str(j["job_id"]) for j in jobs if j["state"] != "HELD"}
    covered = {str(x) for x in gate.get("evaluated_job_ids", [])}
    if not released.issubset(covered):
        reasons.append("assessment does not cover every released job")
    accepted_completed = sum(j["state"] == "COMPLETED" and bool(j["accepted"]) for j in jobs)
    try:
        gate_count = int(gate.get("accepted_success_count", -1))
    except (TypeError, ValueError):
        gate_count = -1
    if gate_count != accepted_completed:
        reasons.append(f"accepted_success_count mismatch: gate={gate_count}, states={accepted_completed}")
    return not reasons, reasons


def plan_group(jobs: list[dict[str, object]], gate_root: Path, limits: dict[str, float] | None = None) -> dict[str, object]:
    jobs = sorted(jobs, key=lambda j: (int(j["wave"]), str(j["job_id"])))
    compound, scenario = str(jobs[0]["compound_id"]), str(jobs[0]["scenario_id"])
    targets = {int(j["target_success"]) for j in jobs}
    if len(targets) != 1:
        raise ValueError(f"{compound}/{scenario}: inconsistent target_success")
    target = targets.pop()
    released = [j for j in jobs if j["state"] != "HELD"]
    active = [j for j in released if j["state"] in ACTIVE_STATES]
    terminal = [j for j in released if j["state"] in TERMINAL_STATES]
    held = [j for j in jobs if j["state"] == "HELD"]
    result: dict[str, object] = {
        "compound_id": compound, "scenario_id": scenario, "target_success": target,
        "counts": {"total": len(jobs), "held": len(held), "released": len(released), "active": len(active), "terminal": len(terminal)},
        "release_job_ids": [], "decision": "WAIT", "reasons": [],
    }
    if active:
        result["decision"] = "WAIT_RELEASED_TASKS_TERMINAL"
        result["reasons"] = ["at least one previously released job is PENDING/RUNNING/COMPLETING"]
        return result
    if not released:
        first_wave = min(int(j["wave"]) for j in held)
        to_release = [str(j["job_id"]) for j in held if int(j["wave"]) == first_wave]
        result.update(decision="RELEASE_WAVE", release_wave=first_wave, release_job_ids=to_release,
                      reasons=["initial wave; no prior released jobs"])
        return result
    assessment = gate_path(gate_root, compound, scenario)
    result["assessment_path"] = str(assessment)
    if not assessment.is_file():
        result["decision"] = "WAIT_GATE_ASSESSMENT"
        result["reasons"] = ["all released jobs are terminal, but no current scientific gate assessment exists"]
        return result
    gate = load_json(assessment)
    current, freshness_reasons = gate_is_current(gate, jobs)
    if not current:
        result["decision"] = "WAIT_GATE_ASSESSMENT"
        result["reasons"] = freshness_reasons
        return result
    passed, reasons = evaluate_gate(gate, jobs, target, limits)
    if passed:
        result["decision"] = "STOP_TARGET_AND_COVERAGE_PASSED"
        result["reasons"] = ["target, identity, basin, motif and shift-stability gates passed"]
        return result
    if held:
        max_released_wave = max(int(j["wave"]) for j in released)
        candidate_waves = sorted({int(j["wave"]) for j in held if int(j["wave"]) > max_released_wave})
        if not candidate_waves:
            candidate_waves = sorted({int(j["wave"]) for j in held})
        next_wave = candidate_waves[0]
        to_release = [str(j["job_id"]) for j in held if int(j["wave"]) == next_wave]
        result.update(decision="RELEASE_WAVE", release_wave=next_wave, release_job_ids=to_release,
                      reasons=reasons)
    else:
        result["decision"] = "EXHAUSTED_FAIL_CLOSED_MANUAL_ESCALATION"
        result["reasons"] = reasons + ["candidate ceiling exhausted"]
    return result


def extract_scheduler_id(stdout: str) -> str:
    match = re.search(r"(?:Submitted batch job\s+)?([0-9]+)", stdout)
    if not match:
        raise ValueError(f"could not parse scheduler job id from: {stdout!r}")
    return match.group(1)


def write_manifest(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    for needed in ("state", "scheduler_job_id"):
        if needed not in fields:
            fields.append(needed)
    from io import StringIO
    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
    atomic_write_text(path, stream.getvalue())


def apply_plan(plan: list[dict[str, object]], jobs: list[dict[str, object]], raw: list[dict[str, str]], fields: list[str],
               manifest: Path, root: Path, site: dict[str, object]) -> list[dict[str, object]]:
    by_id = {str(j["job_id"]): j for j in jobs}
    raw_by_id = {first(r, "job_id", "candidate_id", "start_id"): r for r in raw}
    submission_results = []
    submit_prefix = [str(x) for x in site.get("submit_command", [])]
    if not submit_prefix:
        raise ValueError("site config has no submit_command")
    for group in plan:
        decision = str(group["decision"])
        if decision == "STOP_TARGET_AND_COVERAGE_PASSED":
            c, s = str(group["compound_id"]), str(group["scenario_id"])
            base = root / "gates" / "goat" / re.sub(r"[^A-Za-z0-9_.-]+", "_", c) / re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
            released = [str(j["job_id"]) for j in jobs if j["compound_id"] == c and j["scenario_id"] == s and j["state"] != "HELD"]
            held = [str(j["job_id"]) for j in jobs if j["compound_id"] == c and j["scenario_id"] == s and j["state"] == "HELD"]
            request = {"schema_version": "1.0", "created_utc": utc_now(), "action": "DO_NOT_RELEASE_REMAINING_HELD", "released_terminal_job_ids": released, "held_job_ids": held, "explicitly_cancelled_job_ids": []}
            atomic_write_json(base / "stop_request.json", request)
            atomic_write_json(base / "stop_final_receipt.json", {**request, "verified_no_active_released_jobs": True, "note": "No RUNNING/COMPLETING task was cancelled."})
        if decision != "RELEASE_WAVE":
            continue
        for job_id in group["release_job_ids"]:
            job = by_id[str(job_id)]
            if job["state"] != "HELD":
                raise RuntimeError(f"refusing to submit non-HELD job {job_id}")
            script = Path(str(job["submit_script"]))
            if not script.is_absolute():
                script = (root / script).resolve()
            if not script.is_file():
                raise FileNotFoundError(f"submit script missing for {job_id}: {script}")
            proc = subprocess.run(submit_prefix + [str(script)], text=True, capture_output=True, check=False)
            record = {"job_id": job_id, "returncode": proc.returncode, "stdout": proc.stdout, "stderr": proc.stderr, "submitted_utc": utc_now()}
            if proc.returncode == 0:
                scheduler_id = extract_scheduler_id(proc.stdout)
                record["scheduler_job_id"] = scheduler_id
                state = {"schema_version": "1.0", "job_id": job_id, "state": "PENDING", "scheduler_job_id": scheduler_id, "updated_utc": utc_now(), "accepted": False}
                atomic_write_json(Path(job["job_dir"]) / "job_state.json", state)
                raw_by_id[str(job_id)]["state"] = "PENDING"
                raw_by_id[str(job_id)]["scheduler_job_id"] = scheduler_id
            submission_results.append(record)
    write_manifest(manifest, raw, fields)
    return submission_results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--gate-root", type=Path)
    parser.add_argument("--site-config", type=Path)
    parser.add_argument("--policy", type=Path)
    parser.add_argument("--plan-output", type=Path)
    parser.add_argument("--apply", action="store_true", help="Submit eligible HELD jobs. Without this flag the command is read-only.")
    args = parser.parse_args()
    root = args.package_root.resolve()
    manifest = args.manifest.resolve()
    gate_root = (args.gate_root or (root / "gates" / "goat_assessments")).resolve()
    raw, jobs, fields = load_manifest(manifest, root)
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for job in jobs:
        grouped[(str(job["compound_id"]), str(job["scenario_id"]))].append(job)
    policy_path = (args.policy or (root / "config" / "adaptive_goat_policy.json")).resolve()
    policy = load_json(policy_path)
    raw_limits = policy.get("suggested_ensemble_delta_gates", {})
    limits = {
        "h_nonexchangeable_max_delta_ppm": float(raw_limits["nonexchangeable_1H_ppm"]),
        "c_max_delta_ppm": float(raw_limits["13C_ppm"]),
        "nh_max_delta_ppm": float(raw_limits["exchangeable_NH_ppm"]),
    }
    decisions = [plan_group(grouped[key], gate_root, limits) for key in sorted(grouped)]
    result: dict[str, object] = {
        "schema_version": "1.0", "created_utc": utc_now(), "mode": "APPLY" if args.apply else "DRY_RUN",
        "manifest": str(manifest), "manifest_state_sha256": manifest_digest(jobs), "policy": str(policy_path), "limits": limits, "decisions": decisions,
        "safety": {"cancels_running_jobs": False, "releases_only_held": True, "requires_released_terminal_before_next_wave": True},
    }
    if args.apply:
        if not args.site_config:
            raise SystemExit("--apply requires --site-config")
        site = load_json(args.site_config.resolve())
        result["submissions"] = apply_plan(decisions, jobs, raw, fields, manifest, root, site)
    output = args.plan_output or (root / "package_metadata" / "last_adaptive_goat_plan.json")
    atomic_write_json(output, result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
