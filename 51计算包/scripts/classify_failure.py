#!/usr/bin/env python3
"""Classify ORCA/HPC failures without applying an unsafe generic retry."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


RULES = (
    ("FILESYSTEM_READ_ONLY", (r"read-only file system", r"input/output error", r"stale file handle", r"no space left on device"),
     "Do not retry in the same directory. Verify node-local scratch, quota and stage-out; move computation off the shared/spool filesystem."),
    ("MPI_FABRIC", (r"openfabrics", r"openib", r"ofi_init", r"btl.*failed", r"orte|orted", r"pmix.*error"),
     "Verify that ORCA 6.1.1 and its matching OpenMPI runtime/module are used. Do not wrap the ORCA driver in mpirun. Change fabric variables only under the site administrator's documented policy."),
    ("MEMORY_EXHAUSTION", (r"out of memory", r"oom-kill", r"cannot allocate memory", r"mp2 response", r"not enough memory", r"std::bad_alloc"),
     "Move to the declared memory retry ladder: fewer ORCA ranks and/or a higher-memory profile. Keep nprocs*maxcore within the package memory-fraction limit."),
    ("SCF_NONCONVERGENCE", (r"scf not converged", r"scf failed", r"no convergence in scf", r"error termination.*scf"),
     "Do not treat as a memory retry. Inspect charge/spin/geometry, restart with a documented SCF repair branch, and preserve the failed attempt."),
    ("VIRTUAL_FREQUENCY", (r"imaginary frequenc", r"negative frequenc", r"saddle point", r"not a minimum"),
     "Re-optimise from the displaced mode or a different basin. Do not accept the structure for thermochemical weighting until the frequency gate passes."),
    ("STALE_LOCK", (r"stale lock", r"\blocked\b", r"lock already exists", r"heartbeat.*stale"),
     "Query the scheduler and inspect the heartbeat. Clear the lock only with a signed/manual clearance receipt proving the former job is terminal."),
    ("WALLTIME", (r"time limit", r"due to time limit", r"walltime", r"cancelled.*time"),
     "Restart from a scientifically valid checkpoint when the stage supports it, or resubmit from the approved geometry with a longer declared walltime. GOAT is not assumed to have an exact restart."),
    ("ORCA_ABNORMAL_TERMINATION", (r"orca finished by error termination", r"orca terminated abnormally", r"error termination"),
     "Inspect the first upstream error and classify it before retrying; the final ORCA error line alone is not a root cause."),
)


def classify_text(text: str) -> dict[str, object]:
    lower = text.lower()
    matches = []
    for category, patterns, action in RULES:
        hit = [pattern for pattern in patterns if re.search(pattern, lower, re.IGNORECASE | re.MULTILINE)]
        if hit:
            matches.append({"category": category, "matched_patterns": hit, "recommended_action": action})
    primary = matches[0]["category"] if matches else "UNKNOWN"
    return {
        "schema_version": "1.0",
        "primary_category": primary,
        "matches": matches,
        "automatic_retry_permitted": False,
        "decision": "MANUAL_OR_RULE_SPECIFIC_REVIEW_REQUIRED",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("logs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    chunks = []
    missing = []
    for path in args.logs:
        if not path.is_file():
            missing.append(str(path))
        else:
            chunks.append(f"\n===== {path} =====\n" + path.read_text(encoding="utf-8", errors="replace"))
    result = classify_text("".join(chunks))
    result["source_logs"] = [str(x.resolve()) for x in args.logs if x.is_file()]
    result["missing_logs"] = missing
    payload = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        from package_io import atomic_write_text
        atomic_write_text(args.output, payload)
    else:
        print(payload, end="")
    return 2 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
