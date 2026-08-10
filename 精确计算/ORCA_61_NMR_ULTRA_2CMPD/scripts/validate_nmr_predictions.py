"""Score pre-declared calculated NMR signal centres against locked targets.

The input must contain one row per predicted observable signal, not one row per
atom or transition.  This script performs assignment and scoring only; it never
fits peak-specific offsets, linewidths, intensities, or a regression model.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linear_sum_assignment


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_TARGETS = PROJECT / "validation" / "experimental_peak_targets.csv"
DEFAULT_POLICY = PROJECT / "validation" / "validation_policy.json"
DEFAULT_SOLUTION_EQUIVALENCE = PROJECT / "config" / "solution_equivalence_groups.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def normalize_exchange(row: dict[str, str], nucleus: str) -> str:
    if nucleus != "1H":
        return "not_applicable"
    token = (row.get("exchange_class") or row.get("exchangeable") or "").strip().lower()
    if token in {"exchangeable", "true", "1", "yes"}:
        return "exchangeable"
    if token in {"nonexchangeable", "false", "0", "no"}:
        return "nonexchangeable"
    raise ValueError("Every 1H predicted signal must declare exchange_class")


def assign(targets: list[float], predicted: list[float]) -> tuple[list[dict[str, float]], float | None]:
    """Return an optimal assignment for an explicitly auxiliary diagnostic only."""
    if not targets or not predicted:
        return [], None
    cost = np.abs(np.asarray(targets)[:, None] - np.asarray(predicted)[None, :])
    target_idx, predicted_idx = linear_sum_assignment(cost)
    pairs = [
        {
            "target_ppm": float(targets[i]),
            "predicted_ppm": float(predicted[j]),
            "error_ppm": float(predicted[j] - targets[i]),
            "absolute_error_ppm": float(cost[i, j]),
        }
        for i, j in zip(target_idx.tolist(), predicted_idx.tolist())
    ]
    return pairs, float(np.mean([row["absolute_error_ppm"] for row in pairs]))


def fixed_assign(
    targets: list[dict[str, str]], predicted: list[dict[str, str]]
) -> tuple[list[dict[str, Any]], float | None, bool]:
    target_by_id: dict[str, dict[str, str]] = {}
    predicted_by_id: dict[str, dict[str, str]] = {}
    for row in targets:
        signal_id = row.get("peak_id", "").strip()
        if not signal_id or signal_id in target_by_id:
            raise ValueError(f"Blank or duplicate target peak_id: {signal_id!r}")
        target_by_id[signal_id] = row
    for row in predicted:
        signal_id = row.get("signal_id", "").strip()
        if not signal_id or signal_id in predicted_by_id:
            raise ValueError(f"Blank or duplicate predicted signal_id: {signal_id!r}")
        predicted_by_id[signal_id] = row
    pairs: list[dict[str, Any]] = []
    for signal_id in sorted(set(target_by_id) & set(predicted_by_id)):
        target_ppm = float(target_by_id[signal_id]["ppm"])
        predicted_ppm = float(predicted_by_id[signal_id]["ppm"])
        pairs.append(
            {
                "signal_id": signal_id,
                "target_ppm": target_ppm,
                "predicted_ppm": predicted_ppm,
                "error_ppm": predicted_ppm - target_ppm,
                "absolute_error_ppm": abs(predicted_ppm - target_ppm),
                "assignment_source": predicted_by_id[signal_id]["assignment_source"],
                "assignment_group_id": predicted_by_id[signal_id]["assignment_group_id"],
                "assignment_atom_indices_1based": predicted_by_id[signal_id][
                    "assignment_atom_indices_1based"
                ],
                "symmetry_source_sha256": predicted_by_id[signal_id]["symmetry_source_sha256"],
                "review_status": predicted_by_id[signal_id]["review_status"],
            }
        )
    mae = float(np.mean([row["absolute_error_ppm"] for row in pairs])) if pairs else None
    return pairs, mae, set(target_by_id) == set(predicted_by_id)


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predicted-signals", type=Path, required=True)
    parser.add_argument("--targets", type=Path, default=DEFAULT_TARGETS)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument(
        "--solution-equivalence-config",
        type=Path,
        default=DEFAULT_SOLUTION_EQUIVALENCE,
    )
    parser.add_argument("--method-id", default="production_revdsd_pcsseg3")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--request-promotion", action="store_true")
    parser.add_argument("--reviewer")
    parser.add_argument("--confirm-peak-order-and-multiplicity", action="store_true")
    args = parser.parse_args()

    predicted_path = args.predicted_signals.resolve()
    target_path = args.targets.resolve()
    policy_path = args.policy.resolve()
    symmetry_path = args.solution_equivalence_config.resolve()
    predicted_all = read_csv(predicted_path)
    targets = read_csv(target_path)
    policy = json.loads(policy_path.read_text(encoding="utf-8"))

    required = {
        "figure_id",
        "nucleus",
        "signal_id",
        "ppm",
        "method_id",
        "assignment_source",
        "assignment_group_id",
        "assignment_atom_indices_1based",
        "symmetry_source_sha256",
        "review_status",
        "source_bundle_sha256",
    }
    if not predicted_all or not required.issubset(predicted_all[0]):
        raise SystemExit(
            "Predicted signal CSV requires fixed signal IDs, assignment provenance, source hashes, "
            "figure_id,nucleus,ppm,method_id and exchange_class for 1H"
        )
    predicted = [row for row in predicted_all if row["method_id"] == args.method_id]
    if not predicted:
        raise SystemExit(f"No predicted rows for method_id={args.method_id}")
    methods = {row["method_id"] for row in predicted}
    if methods != {args.method_id}:
        raise SystemExit("Method mixing is forbidden")
    allowed_assignment_sources = {
        "DECLARED_SOLUTION_EQUIVALENCE",
        "ATOM_MAP_UNIQUE_SITE",
        "EXCHANGEABLE_ATOM_RESOLVED",
    }
    expected_symmetry_hash = sha256(symmetry_path)
    assignment_provenance_complete = all(
        row["assignment_source"] in allowed_assignment_sources
        and row["review_status"] == "LOCKED_BEFORE_SCORING"
        and bool(re.fullmatch(r"\d+(?:;\d+)*", row["assignment_atom_indices_1based"]))
        and row["symmetry_source_sha256"].upper() == expected_symmetry_hash
        and bool(re.fullmatch(r"[0-9A-Fa-f]{64}", row["source_bundle_sha256"]))
        for row in predicted
    )

    target_groups: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    source_statuses: set[str] = set()
    for row in targets:
        nucleus = row["nucleus"]
        exchange = normalize_exchange(row, nucleus)
        target_groups[(row["figure_id"], nucleus, exchange)].append(row)
        source_statuses.add(row.get("source_status", ""))

    predicted_groups: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in predicted:
        nucleus = row["nucleus"]
        if nucleus not in {"1H", "13C"}:
            raise SystemExit(f"Unsupported nucleus: {nucleus}")
        exchange = normalize_exchange(row, nucleus)
        predicted_groups[(row["figure_id"], nucleus, exchange)].append(row)

    group_results: list[dict[str, Any]] = []
    all_keys = sorted(set(target_groups) | set(predicted_groups))
    for key in all_keys:
        target_values = target_groups.get(key, [])
        predicted_values = predicted_groups.get(key, [])
        pairs, mae, fixed_id_match = fixed_assign(target_values, predicted_values)
        auxiliary_pairs, auxiliary_mae = assign(
            [float(row["ppm"]) for row in target_values],
            [float(row["ppm"]) for row in predicted_values],
        )
        figure_id, nucleus, exchange = key
        group_results.append(
            {
                "figure_id": figure_id,
                "nucleus": nucleus,
                "exchange_class": exchange,
                "target_count": len(target_values),
                "predicted_count": len(predicted_values),
                "count_match": len(target_values) == len(predicted_values),
                "fixed_signal_id_set_match": fixed_id_match,
                "mae_ppm": mae,
                "assignments": pairs,
                "hungarian_auxiliary_diagnostic": {
                    "mae_ppm": auxiliary_mae,
                    "assignments": auxiliary_pairs,
                    "controls_validation": False,
                },
            }
        )

    def aggregate(nucleus: str, exchange: str | None = None) -> float | None:
        errors = [
            pair["absolute_error_ppm"]
            for group in group_results
            if group["nucleus"] == nucleus
            and (exchange is None or group["exchange_class"] == exchange)
            for pair in group["assignments"]
        ]
        return float(np.mean(errors)) if errors else None

    h_nonexchange_mae = aggregate("1H", "nonexchangeable")
    h_exchange_mae = aggregate("1H", "exchangeable")
    c_mae = aggregate("13C")
    thresholds = policy["primary_metrics"]
    counts_pass = all(
        group["count_match"] and group["fixed_signal_id_set_match"]
        for group in group_results
    )
    numerical_pass = (
        h_nonexchange_mae is not None
        and h_nonexchange_mae <= float(thresholds["1H_nonexchangeable_MAE_ppm_target"])
        and c_mae is not None
        and c_mae <= float(thresholds["13C_MAE_ppm_target"])
    )
    review_pass = bool(args.confirm_peak_order_and_multiplicity)
    pass_all = counts_pass and numerical_pass and review_pass and assignment_provenance_complete
    metadata_complete = not any("processed_image_only" in token for token in source_statuses)

    if args.request_promotion and not (args.reviewer or "").strip():
        raise SystemExit("--request-promotion requires --reviewer")
    if args.request_promotion and pass_all and metadata_complete:
        status = "VALIDATED_FOR_THIS_TWO_COMPOUND_PROTOCOL_AND_LOCKED_CONDITIONS"
    elif pass_all:
        status = "PROVISIONAL_PASS_METADATA_PENDING"
    else:
        status = "VALIDATION_NOT_PASSED"

    payload = {
        "schema_version": "1.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "method_id": args.method_id,
        "reviewer": (args.reviewer or "").strip() or None,
        "promotion_requested": bool(args.request_promotion),
        "no_fit_statement": (
            "Primary scoring used pre-locked signal IDs and atom/symmetry provenance; no "
            "peak-specific fitting or regression was performed. Hungarian matching is auxiliary only."
        ),
        "input_provenance": {
            "predicted_signals": str(predicted_path),
            "predicted_signals_sha256": sha256(predicted_path),
            "targets": str(target_path),
            "targets_sha256": sha256(target_path),
            "policy": str(policy_path),
            "policy_sha256": sha256(policy_path),
            "solution_equivalence_config": str(symmetry_path),
            "solution_equivalence_config_sha256": expected_symmetry_hash,
        },
        "metrics": {
            "1H_nonexchangeable_MAE_ppm": h_nonexchange_mae,
            "1H_exchangeable_MAE_ppm_report_only": h_exchange_mae,
            "13C_MAE_ppm": c_mae,
            "all_peak_counts_match": counts_pass,
            "peak_order_and_multiplicity_review_confirmed": review_pass,
            "numerical_thresholds_pass": numerical_pass,
            "holdout_vendor_metadata_complete": metadata_complete,
            "fixed_assignment_provenance_complete": assignment_provenance_complete,
        },
        "groups": group_results,
        "limitations": [
            "Predicted input must contain reviewed observable signal centres; atom-level rows are not silently collapsed.",
            "The fixed signal-ID/atom/symmetry assignment controls validation; minimum-error Hungarian matching is diagnostic only.",
            "Exchangeable 1H errors are reported but do not control promotion.",
            "F3_015 cannot receive protocol-scoped validation status while only processed-image metadata are available.",
        ],
    }
    atomic_json(args.output.resolve(), payload)
    print(json.dumps({"status": status, "metrics": payload["metrics"]}, indent=2))
    return 0 if status != "VALIDATION_NOT_PASSED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
