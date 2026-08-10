"""Static, fail-closed audit for the two-compound ORCA 6.1 NMR package."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PROJECT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    checks: list[dict[str, Any]] = []
    errors: list[str] = []

    def check(name: str, passed: bool, detail: Any) -> None:
        checks.append({"name": name, "status": "PASS" if passed else "ERROR", "detail": detail})
        if not passed:
            errors.append(f"{name}: {detail}")

    registry = rows(PROJECT / "registry/Figure3_52_ORCA_reference_state_registry.csv")
    check("registry_scope", {row["figure_id"] for row in registry} == {"F3_005", "F3_015"}, registry)
    check(
        "identity_lock",
        {(row["figure_id"], row["formula"], row["net_charge_as_drawn"]) for row in registry}
        == {("F3_005", "C7H7NO", "0"), ("F3_015", "C5H10N2O", "0")},
        "benzamide and neutral (S)-prolinamide",
    )

    manifest_path = PROJECT / "data/GOAT_explicit_job_manifest.csv"
    manifest = rows(manifest_path)
    counts = Counter((row["figure_id"], row["scenario_id"]) for row in manifest)
    check("goat_job_count", len(manifest) == 200, len(manifest))
    check("goat_unique_job_ids", len({row["job_id"] for row in manifest}) == 200, len({row["job_id"] for row in manifest}))
    check(
        "goat_scenario_count",
        len(counts) == 10 and set(counts.values()) == {20},
        [
            {"figure_id": key[0], "scenario_id": key[1], "job_count": value}
            for key, value in sorted(counts.items())
        ],
    )
    check("goat_enabled", all(row.get("enabled") == "true" and row.get("readiness") == "ready" for row in manifest), "all 200 rows")

    bad_paths: list[str] = []
    bad_hashes: list[str] = []
    bad_resources: list[str] = []
    for row in manifest:
        raw = row["portable_job_directory"].replace("\\", "/")
        if Path(raw).is_absolute() or re.match(r"^[A-Za-z]:/", raw) or ".." in Path(raw).parts:
            bad_paths.append(row["job_id"])
            continue
        directory = PROJECT / raw
        input_path = directory / row["input_file"]
        coordinate = directory / row["coordinate_file"]
        if not input_path.is_file() or not coordinate.is_file():
            bad_paths.append(row["job_id"])
            continue
        if sha256(input_path) != row["expected_input_sha256"] or sha256(coordinate) != row["expected_coordinate_sha256"]:
            bad_hashes.append(row["job_id"])
        text = input_path.read_text(encoding="utf-8")
        if not (re.search(r"%pal\s+nprocs\s+8\s+end", text, re.I) and re.search(r"NWorkers\s+8", text, re.I) and re.search(r"%maxcore\s+2000", text, re.I)):
            bad_resources.append(row["job_id"])
    check("portable_paths_and_files", not bad_paths, bad_paths[:10])
    check("goat_hashes", not bad_hashes, bad_hashes[:10])
    check("goat_resources", not bad_resources, bad_resources[:10])

    reference = rows(PROJECT / "data/reference_job_manifest.csv")
    enabled_reference = [row["job_id"] for row in reference if row.get("enabled") == "true"]
    check("reference_jobs", len(reference) == 6, len(reference))
    check("reference_fail_closed", enabled_reference == ["TMS_DMSO__opt_freq"], enabled_reference)
    check(
        "same_level_tms_production",
        any(row.get("method_id") == "production_revdsd_pcsseg3" for row in reference),
        [row.get("method_id") for row in reference],
    )

    gate = json.loads((PROJECT / "jobs/high_level_manual_opt_in_gate.json").read_text(encoding="utf-8"))
    check(
        "high_level_gate",
        gate.get("approval_status") == "DRAFT_NOT_APPROVED",
        gate.get("approval_status"),
    )

    syntax_failures: list[str] = []
    for path in PROJECT.rglob("*.inp"):
        text = path.read_text(encoding="utf-8", errors="replace")
        if "%eprnmr" in text.lower():
            if text.lower().find("* xyzfile") > text.lower().find("%eprnmr"):
                syntax_failures.append(path.relative_to(PROJECT).as_posix() + ": EPRNMR before coordinates")
            for match in re.finditer(r"Nuclei\s*=\s*([0-9,]+)", text, re.I):
                if min(int(token) for token in match.group(1).split(",")) < 1:
                    syntax_failures.append(path.relative_to(PROJECT).as_posix() + ": Nuclei is not 1-based")
    check("orca_nmr_index_and_block_order", not syntax_failures, syntax_failures[:20])

    status = json.loads((PROJECT / "package_metadata/STATUS.json").read_text(encoding="utf-8"))
    check("truth_status", status.get("hpc_execution_state") == "NOT_RUN" and status.get("calculated_spectra_included") is False, status)

    secret_hits: list[str] = []
    # Never copy a user's literal credential into the repository, including into
    # the audit that is meant to detect it.  These patterns catch common secret
    # assignments, GitHub tokens, private keys, and credential-bearing URLs.
    secret_patterns = (
        re.compile(r"(?i)\b(?:password|passwd|pwd)\s*[:=]\s*[^\s,;]+"),
        re.compile(r"\bgh[pousr]_[A-Za-z0-9]{20,}\b"),
        re.compile(r"-----BEGIN (?:OPENSSH|RSA|EC|DSA) PRIVATE KEY-----"),
        re.compile(r"https?://[^\s/@:]+:[^\s/@]+@github\.com", re.I),
    )
    for path in PROJECT.rglob("*"):
        if path.is_file() and path.suffix.lower() in {".py", ".json", ".csv", ".md", ".txt", ".sh", ".inp"}:
            text = path.read_text(encoding="utf-8", errors="ignore")
            if any(pattern.search(text) for pattern in secret_patterns):
                secret_hits.append(path.relative_to(PROJECT).as_posix())
    check("no_user_credentials", not secret_hits, secret_hits)

    report = {
        "schema_version": "1.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if not errors else "FAIL",
        "project": ".",
        "checks": checks,
        "error_count": len(errors),
        "errors": errors,
    }
    output = PROJECT / "package_metadata/static_audit.json"
    output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps({"status": report["status"], "checks": len(checks), "errors": len(errors)}, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
