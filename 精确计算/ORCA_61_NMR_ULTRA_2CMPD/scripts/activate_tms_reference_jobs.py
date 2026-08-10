"""Activate same-method TMS shielding jobs after a reviewed geometry calculation.

This is deliberately fail-closed.  It does not infer that an optimization is
scientifically acceptable: a named reviewer must explicitly confirm that the
TMS optimization terminated normally and has no significant imaginary mode.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT / "data" / "reference_job_manifest.csv"
ACTIVE_MANIFEST = PROJECT / "data" / "reference_shielding_active_manifest.csv"
ROOT = PROJECT / "reference_jobs" / "TMS_DMSO"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def atomic_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument(
        "--active-manifest",
        type=Path,
        default=ACTIVE_MANIFEST,
        help="Stage-specific manifest containing only the explicitly activated shielding rows",
    )
    parser.add_argument("--method-id", action="append", required=True)
    parser.add_argument("--reviewer", required=True)
    parser.add_argument("--confirm-normal-termination", action="store_true")
    parser.add_argument("--confirm-no-significant-imaginary-frequency", action="store_true")
    args = parser.parse_args()

    if not args.reviewer.strip():
        raise SystemExit("A non-empty named reviewer is required")
    if not args.confirm_normal_termination or not args.confirm_no_significant_imaginary_frequency:
        raise SystemExit("Both explicit geometry review confirmations are required")

    output = ROOT / "tms_opt.out"
    geometry = ROOT / "tms_opt.xyz"
    if not output.is_file() or not geometry.is_file():
        raise SystemExit("Missing reference_jobs/TMS_DMSO/tms_opt.out or tms_opt.xyz")
    text = output.read_text(encoding="utf-8", errors="replace")
    if "ORCA TERMINATED NORMALLY" not in text:
        raise SystemExit("TMS optimization output did not terminate normally")
    if geometry.stat().st_size < 100:
        raise SystemExit("TMS optimized geometry is unexpectedly small")

    manifest = args.manifest.resolve()
    with manifest.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    requested = set(args.method_id)
    found: set[str] = set()
    digest = sha256(geometry)
    reviewed_at = datetime.now(timezone.utc).isoformat()
    active_rows: list[dict[str, str]] = []
    for row in rows:
        method = row.get("method_id", "")
        if method not in requested:
            continue
        if row.get("stage") != "reference_shielding":
            raise SystemExit(f"Requested method is not a reference shielding row: {method}")
        destination = PROJECT / row["portable_job_directory"] / row["coordinate_file"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(geometry, destination)
        row["expected_coordinate_sha256"] = digest
        row["enabled"] = "true"
        # Keep readiness machine-readable.  Scientific approval is carried by
        # the exact stage-specific approval token below, not by prose in this
        # field.
        row["readiness"] = "ready"
        row["approval_status"] = "APPROVED_FOR_REFERENCE_SHIELDING_CALCULATION"
        row["reviewer"] = args.reviewer.strip()
        row["reviewed_at"] = reviewed_at
        active_rows.append(dict(row))
        found.add(method)
    missing = sorted(requested - found)
    if missing:
        raise SystemExit(f"Unknown reference method id(s): {', '.join(missing)}")
    atomic_csv(manifest, rows)
    active_manifest = args.active_manifest.resolve()
    try:
        active_manifest.relative_to(PROJECT.resolve())
    except ValueError as exc:
        raise SystemExit("--active-manifest must remain inside the package root") from exc
    if not active_rows:
        raise SystemExit("No reference shielding rows were activated")
    atomic_csv(active_manifest, active_rows)
    print(
        json.dumps(
            {
                "status": "ACTIVATED_AFTER_MANUAL_REVIEW",
                "methods": sorted(found),
                "reviewer": args.reviewer.strip(),
                "reviewed_at": reviewed_at,
                "optimized_geometry_sha256": digest,
                "active_manifest": str(active_manifest.relative_to(PROJECT)).replace("\\", "/"),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
