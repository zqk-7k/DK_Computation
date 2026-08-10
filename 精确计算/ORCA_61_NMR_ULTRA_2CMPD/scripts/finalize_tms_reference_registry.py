"""Create the approved same-method TMS registry from completed ORCA outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from nmr_method_validation import MethodValidationError, validate_nmr_method_evidence
from runner_provenance import RunnerProvenanceError, validate_runner_output


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = PROJECT / "data/reference_job_manifest.csv"
DEFAULT_OUTPUT = PROJECT / "config/reference_job_registry.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--method-id", action="append", required=True)
    parser.add_argument("--reviewer", required=True)
    parser.add_argument("--runner-state-root", type=Path)
    parser.add_argument("--confirm-shielding-summary-reviewed", action="store_true")
    args = parser.parse_args()
    if not args.reviewer.strip() or not args.confirm_shielding_summary_reviewed:
        raise SystemExit("Named reviewer and --confirm-shielding-summary-reviewed are required")

    manifest = args.manifest.resolve()
    with manifest.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    index = {row["method_id"]: row for row in rows if row.get("stage") == "reference_shielding"}
    entries: list[dict[str, object]] = []
    for method_id in args.method_id:
        if method_id not in index:
            raise SystemExit(f"Unknown reference method: {method_id}")
        row = index[method_id]
        if row.get("enabled") != "true" or not row.get("readiness", "").startswith("ready"):
            raise SystemExit(f"Reference method is not activated: {method_id}")
        directory = PROJECT / row["portable_job_directory"]
        input_path = directory / row["input_file"]
        coordinate = directory / row["coordinate_file"]
        for path in (input_path, coordinate):
            if not path.is_file():
                raise SystemExit(f"Missing TMS reference artifact: {path}")
        if sha256(input_path) != row["expected_input_sha256"]:
            raise SystemExit(f"TMS input hash mismatch: {method_id}")
        if sha256(coordinate) != row["expected_coordinate_sha256"]:
            raise SystemExit(f"TMS coordinate hash mismatch: {method_id}")
        try:
            runner = validate_runner_output(
                row,
                "reference_shielding",
                PROJECT,
                args.runner_state_root,
            )
            output_path = Path(runner["output_path"])
            method_validation = validate_nmr_method_evidence(
                input_path,
                output_path,
                method_id,
                row["method_fingerprint"],
            )
        except (RunnerProvenanceError, MethodValidationError) as exc:
            raise SystemExit(f"TMS runner/method validation failed for {method_id}: {exc}") from exc
        for nucleus in ("1H", "13C"):
            entries.append(
                {
                    "reference_id": f"TMS_DMSO_{nucleus}_{method_id}",
                    "approval_status": "APPROVED_FOR_NMR_PROPERTIES",
                    "nucleus": nucleus,
                    "solvent_domain": "DMSO-d6",
                    "raw_calculation_reference": "TMS in CPCM(DMSO), calculated model",
                    "experimental_reference_convention": f"TMS delta({nucleus}) = 0 ppm",
                    "method_id": method_id,
                    "method_fingerprint": row["method_fingerprint"],
                    "input_file": str(input_path.resolve()),
                    "coordinate_file": str(coordinate.resolve()),
                    "output_file": str(output_path.resolve()),
                    "input_sha256": sha256(input_path),
                    "coordinate_sha256": sha256(coordinate),
                    "output_sha256": sha256(output_path),
                    "runner_state_file": runner["runner_state"],
                    "runner_state_sha256": runner["runner_state_sha256"],
                    "runner_kind": runner["runner_kind"],
                    "method_validation": method_validation,
                    "vibrational_correction_policy": (
                        "DIAGNOSTIC_ONLY_NOT_APPLIED_TO_MAIN_DELTA"
                    ),
                    "vpt2_correction_applied_to_main_delta": False,
                    "reviewer": args.reviewer.strip(),
                }
            )
    payload = {
        "schema_version": "2.0-ultra",
        "status": "APPROVED_FOR_NMR_PROPERTIES",
        "approved_at_utc": datetime.now(timezone.utc).isoformat(),
        "reviewer": args.reviewer.strip(),
        "reference_scale": "TMS delta(1H)=0 and delta(13C)=0",
        "rule": "Target and reference use the identical method fingerprint; no cross-method reference is permitted.",
        "vibrational_correction_policy": (
            "VPT2 or other vibrational shielding corrections remain diagnostic-only unless "
            "an independently approved, same-level target/reference correction protocol exists."
        ),
        "reference_jobs": entries,
        "notice": "CALCULATED MODEL - NOT EXPERIMENTAL DATA",
    }
    output = args.output.resolve()
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    temporary.replace(output)
    print(json.dumps({"status": payload["status"], "entries": len(entries), "output": str(output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
