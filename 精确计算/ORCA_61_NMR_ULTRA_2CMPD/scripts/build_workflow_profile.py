#!/usr/bin/env python3
"""Materialize a runnable GOAT manifest for a declared workflow profile.

The source GOAT inputs remain immutable and hash locked.  A profile selects a
scientifically declared subset of fixed-composition scenarios and deterministic
starts.  It never combines energies or populations across solvent
stoichiometries and never launches ORCA.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT / "config" / "workflow_profiles.json"
DEFAULT_MANIFEST = PROJECT / "data" / "GOAT_explicit_job_manifest.csv"
DEFAULT_SCENARIOS = PROJECT / "data" / "explicit_scenario_registry.csv"
DEFAULT_OUTPUT_ROOT = PROJECT / "data" / "workflow_profiles"


class ProfileError(RuntimeError):
    """Raised when a profile cannot be materialized without ambiguity."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def safe_profile(value: str) -> str:
    token = value.strip().upper()
    if not re.fullmatch(r"[A-Z0-9_]+", token):
        raise ProfileError(f"Unsafe workflow profile identifier: {value!r}")
    return token


def display_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def load_config(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ProfileError(f"Cannot read workflow profile configuration {path}: {exc}") from exc
    if not isinstance(payload.get("profiles"), dict) or not payload["profiles"]:
        raise ProfileError("workflow_profiles.json must contain profiles")
    return payload


def materialize(
    config_path: Path,
    source_manifest: Path,
    source_scenarios: Path,
    output_root: Path,
    profile_name: str | None = None,
) -> dict[str, Any]:
    payload = load_config(config_path)
    selected_name = safe_profile(profile_name or str(payload.get("default_profile", "")))
    if selected_name not in payload["profiles"]:
        raise ProfileError(f"Unknown workflow profile: {selected_name}")
    profile = payload["profiles"][selected_name]
    starts_required = int(profile.get("goat", {}).get("starts_per_scenario", 0))
    if starts_required < 1:
        raise ProfileError(f"Profile {selected_name} has no valid starts_per_scenario")

    manifest_rows = read_csv(source_manifest)
    scenario_rows = read_csv(source_scenarios)
    if not manifest_rows or not scenario_rows:
        raise ProfileError("Source GOAT manifest/scenario registry is empty")
    jobs_by_key: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in manifest_rows:
        jobs_by_key.setdefault((row["figure_id"], row["scenario_id"]), []).append(row)
    registry_by_key = {(row["figure_id"], row["scenario_id"]): row for row in scenario_rows}

    selected_jobs: list[dict[str, Any]] = []
    selected_scenarios: list[dict[str, Any]] = []
    seen_job_ids: set[str] = set()
    profile_hash = sha256(config_path)
    source_manifest_hash = sha256(source_manifest)
    source_scenario_hash = sha256(source_scenarios)
    compounds = profile.get("compounds")
    if not isinstance(compounds, dict) or not compounds:
        raise ProfileError(f"Profile {selected_name} has no compounds")

    for figure_id, compound in compounds.items():
        declared = compound.get("scenarios") if isinstance(compound, dict) else None
        if not isinstance(declared, list) or not declared:
            raise ProfileError(f"Profile {selected_name}/{figure_id} has no scenarios")
        for scenario in declared:
            scenario_id = str(scenario.get("id", ""))
            key = (figure_id, scenario_id)
            if key not in registry_by_key or key not in jobs_by_key:
                raise ProfileError(f"Missing source scenario/input jobs for {figure_id}/{scenario_id}")
            candidates = sorted(jobs_by_key[key], key=lambda row: int(row["start_index"]))
            if len(candidates) < starts_required:
                raise ProfileError(
                    f"{figure_id}/{scenario_id} has {len(candidates)} starts; {starts_required} required"
                )
            role = str(scenario.get("role", "UNDECLARED"))
            eligible = bool(scenario.get("main_prediction_eligibility", False))
            registry_row: dict[str, Any] = dict(registry_by_key[key])
            registry_row.update(
                {
                    "minimum_goat_starts": starts_required,
                    "workflow_profile": selected_name,
                    "scenario_role": role,
                    "main_prediction_eligibility": str(eligible).lower(),
                    "cross_stoichiometry_weight": "NONE",
                    "profile_config_sha256": profile_hash,
                    "source_scenario_registry_sha256": source_scenario_hash,
                }
            )
            selected_scenarios.append(registry_row)
            for source_row in candidates[:starts_required]:
                row: dict[str, Any] = dict(source_row)
                if row["job_id"] in seen_job_ids:
                    raise ProfileError(f"Duplicate selected GOAT job: {row['job_id']}")
                seen_job_ids.add(row["job_id"])
                row.update(
                    {
                        "workflow_profile": selected_name,
                        "scenario_role": role,
                        "main_prediction_eligibility": str(eligible).lower(),
                        "cross_stoichiometry_weight": "NONE",
                        "profile_config_sha256": profile_hash,
                        "source_master_manifest_sha256": source_manifest_hash,
                    }
                )
                selected_jobs.append(row)

    selected_jobs.sort(
        key=lambda row: (int(row["compound_number"]), str(row["scenario_id"]), int(row["start_index"]))
    )
    selected_scenarios.sort(key=lambda row: (int(row["compound_number"]), str(row["scenario_id"])))
    destination = output_root / selected_name
    manifest_path = destination / "GOAT_explicit_job_manifest.csv"
    scenario_path = destination / "explicit_scenario_registry.csv"
    summary_path = destination / "profile_summary.json"
    write_csv(manifest_path, selected_jobs, list(selected_jobs[0]))
    write_csv(scenario_path, selected_scenarios, list(selected_scenarios[0]))
    summary = {
        "schema_version": "1.0",
        "notice": "CALCULATED WORKFLOW SELECTION - NOT EXPERIMENTAL DATA",
        "workflow_profile": selected_name,
        "profile_status": profile.get("status"),
        "default_profile": selected_name == safe_profile(str(payload.get("default_profile", ""))),
        "compounds": len({row["figure_id"] for row in selected_jobs}),
        "fixed_composition_scenarios": len(selected_scenarios),
        "goat_jobs": len(selected_jobs),
        "starts_per_scenario": starts_required,
        "mandatory_trace_water_scenarios": sum(
            row["scenario_role"] == "MANDATORY_TRACE_WATER_SENSITIVITY"
            for row in selected_scenarios
        ),
        "cross_stoichiometry_boltzmann_mixing": "FORBIDDEN",
        "source_manifest": display_path(source_manifest),
        "source_manifest_sha256": source_manifest_hash,
        "source_scenario_registry": display_path(source_scenarios),
        "source_scenario_registry_sha256": source_scenario_hash,
        "profile_config": display_path(config_path),
        "profile_config_sha256": profile_hash,
        "manifest": display_path(manifest_path),
        "manifest_sha256": sha256(manifest_path),
        "scenario_registry": display_path(scenario_path),
        "scenario_registry_sha256": sha256(scenario_path),
    }
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--profile", help="Profile ID; defaults to config default_profile")
    value.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    value.add_argument("--source-manifest", type=Path, default=DEFAULT_MANIFEST)
    value.add_argument("--source-scenarios", type=Path, default=DEFAULT_SCENARIOS)
    value.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    return value


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        summary = materialize(
            args.config.resolve(),
            args.source_manifest.resolve(),
            args.source_scenarios.resolve(),
            args.output_root.resolve(),
            args.profile,
        )
    except (ProfileError, OSError, ValueError, KeyError) as exc:
        print(f"Workflow profile build failed: {exc}")
        return 2
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
