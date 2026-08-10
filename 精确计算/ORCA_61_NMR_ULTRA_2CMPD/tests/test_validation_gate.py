from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT / "scripts" / "validate_nmr_predictions.py"
TARGETS = PROJECT / "validation" / "experimental_peak_targets.csv"


class ValidationGateTests(unittest.TestCase):
    def predicted_from_targets(self, destination: Path) -> None:
        with TARGETS.open("r", encoding="utf-8-sig", newline="") as handle:
            targets = list(csv.DictReader(handle))
        fields = (
            "figure_id",
            "nucleus",
            "signal_id",
            "ppm",
            "exchange_class",
            "method_id",
            "assignment_source",
            "assignment_group_id",
            "assignment_atom_indices_1based",
            "symmetry_source_sha256",
            "source_bundle_sha256",
            "review_status",
        )
        symmetry_hash = __import__("hashlib").sha256(
            (PROJECT / "config" / "solution_equivalence_groups.json").read_bytes()
        ).hexdigest().upper()
        with destination.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            for row in targets:
                writer.writerow(
                    {
                        "figure_id": row["figure_id"],
                        "nucleus": row["nucleus"],
                        "signal_id": row["peak_id"],
                        "ppm": row["ppm"],
                        "exchange_class": row["exchange_class"],
                        "method_id": "production_revdsd_pcsseg3",
                        "assignment_source": (
                            "EXCHANGEABLE_ATOM_RESOLVED"
                            if row["exchange_class"] == "exchangeable"
                            else "ATOM_MAP_UNIQUE_SITE"
                        ),
                        "assignment_group_id": row["peak_id"],
                        "assignment_atom_indices_1based": "1",
                        "symmetry_source_sha256": symmetry_hash,
                        "source_bundle_sha256": "A" * 64,
                        "review_status": "LOCKED_BEFORE_SCORING",
                    }
                )

    def test_perfect_values_remain_provisional_when_holdout_metadata_are_missing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            predicted = root / "predicted.csv"
            output = root / "report.json"
            self.predicted_from_targets(predicted)
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--predicted-signals",
                    str(predicted),
                    "--output",
                    str(output),
                    "--confirm-peak-order-and-multiplicity",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            report = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(report["status"], "PROVISIONAL_PASS_METADATA_PENDING")
            self.assertFalse(report["metrics"]["holdout_vendor_metadata_complete"])

    def test_missing_signal_fails_count_gate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            predicted = root / "predicted.csv"
            output = root / "report.json"
            self.predicted_from_targets(predicted)
            rows = predicted.read_text(encoding="utf-8").splitlines()
            predicted.write_text("\n".join(rows[:-1]) + "\n", encoding="utf-8")
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--predicted-signals",
                    str(predicted),
                    "--output",
                    str(output),
                    "--confirm-peak-order-and-multiplicity",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
            report = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(report["status"], "VALIDATION_NOT_PASSED")


if __name__ == "__main__":
    unittest.main()
