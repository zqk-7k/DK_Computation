"""Source-bound multiplet comparisons, not automatic scientific approval.

Experimental inputs must be supplied by the user. Synthetic engineering tests
exercise the numeric helpers but are never installed as experimental panels.
Every production invocation re-collects real completed ORCA jobs read-only.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from .collect import finite, inside, read_json, sha256_file, unique_output
from .spectra import make_spectrum


DEFAULT_TOLERANCES = {
    "absolute_position_rmse_ppm": 0.15,
    "centered_position_rmse_Hz": 1.0,
    "adjacent_spacing_rmse_Hz": 1.0,
    "relative_line_area_mae": 0.10,
    "relative_region_integral_error": 0.15,
    "line_count_difference": 0,
    "strong_coupling_ratio_warning": 10.0,
}


def parameter_digest(atoms, pairs, nucleus, field_MHz):
    data = {"atom_rows": atoms, "j_pairs": pairs, "nucleus": nucleus,
            "field_MHz": float(field_MHz)}
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def _human_review(value):
    if not isinstance(value, dict):
        raise ValueError("A named experimental assignment review is required")
    for key in ("reviewer", "date", "notes"):
        if not isinstance(value.get(key), str) or not value[key].strip():
            raise ValueError("Review requires reviewer, date and notes")
    if value["reviewer"].strip().lower() in {"auto", "automatic", "codex", "claude", "test", "ai"}:
        raise ValueError("Software cannot certify experimental assignments")


def _artifact(root, row):
    if not isinstance(row, dict):
        raise ValueError("Artifact path and sha256 are required")
    path = inside(root, row.get("path"))
    if not path.is_file() or sha256_file(path) != row.get("sha256"):
        raise ValueError("Missing or changed source artifact")
    return path


def _positive(value, label, zero=False):
    result = finite(value, label)
    if result < 0 or (not zero and result == 0):
        raise ValueError(f"{label} must be {'nonnegative' if zero else 'positive'}")
    return result


def _range(value, label):
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{label} must be a two-value list")
    low, high = sorted(finite(x, label) for x in value)
    if low == high:
        raise ValueError(f"{label} must have nonzero width")
    return low, high


def resolve_transitions(lines, field_MHz, resolution_Hz, relative_area_floor=0.01):
    """Merge unresolved ideal transitions; NOT a fitted experimental peak picker.

    Each group spans at most resolution_Hz (not transitive nearest-neighbor
    chaining). Positions are integrated-area centroids. Comparison therefore
    requires experimental fitted/assigned component AREAS, not peak heights.
    """
    field = _positive(field_MHz, "field")
    resolution = _positive(resolution_Hz, "resolution")
    floor = _positive(relative_area_floor, "relative area floor", zero=True)
    if floor >= 1:
        raise ValueError("Relative area floor must be below one")
    data = sorted((finite(r["ppm"], "line ppm"), _positive(r["intensity"], "line area", zero=True))
                  for r in lines)
    groups = []
    for ppm, area in data:
        if area == 0:
            continue
        if not groups or (ppm - groups[-1][0][0]) * field > resolution:
            groups.append([])
        groups[-1].append((ppm, area))
    resolved = [{"ppm": sum(x*a for x, a in g)/sum(a for _, a in g),
                 "area": sum(a for _, a in g), "transition_count": len(g)} for g in groups]
    if not resolved:
        return []
    maximum = max(r["area"] for r in resolved)
    return [r for r in resolved if r["area"] >= maximum * floor]


def compare_lines(observed, predicted, field_MHz):
    """Numerical engineering helper. Inputs are resolved lines with ppm/area."""
    field = _positive(field_MHz, "field")
    def normalized(rows):
        rows = sorted(({"ppm": finite(r["ppm"], "ppm"), "area": _positive(r["area"], "area")}
                       for r in rows), key=lambda r: r["ppm"])
        total = sum(r["area"] for r in rows)
        return rows, np.array([r["ppm"] for r in rows]), np.array([r["area"]/total for r in rows])
    obs, ox, oa = normalized(observed)
    pred, px, pa = normalized(predicted)
    result = {"experimental_line_count": len(obs), "model_resolved_line_count": len(pred),
              "line_count_difference": abs(len(obs)-len(pred)),
              "absolute_position_rmse_ppm": None, "centered_position_rmse_Hz": None,
              "adjacent_spacing_rmse_Hz": None, "relative_line_area_mae": None,
              "experimental_adjacent_spacings_Hz": (np.diff(ox)*field).tolist(),
              "model_adjacent_spacings_Hz": (np.diff(px)*field).tolist(), "matched_lines": []}
    if not obs or not pred:
        return result
    ii, jj = linear_sum_assignment(abs(ox[:, None]-px[None, :]))
    result["absolute_position_rmse_ppm"] = float(np.sqrt(np.mean((ox[ii]-px[jj])**2)))
    result["matched_lines"] = [{"experimental_ppm": float(ox[i]), "model_ppm": float(px[j]),
                                 "position_error_Hz": float((px[j]-ox[i])*field)} for i, j in zip(ii, jj)]
    result["centroid_error_ppm"] = float(np.dot(px, pa)-np.dot(ox, oa))
    if len(obs) == len(pred):
        # Sorted one-to-one pattern comparison; no J refit or per-line shifting.
        result["centered_position_rmse_Hz"] = float(np.sqrt(np.mean(((px-px.mean())-(ox-ox.mean()))**2))*field)
        result["relative_line_area_mae"] = float(np.mean(abs(pa-oa)))
        if len(obs) > 1:
            result["adjacent_spacing_rmse_Hz"] = float(np.sqrt(np.mean((np.diff(px)-np.diff(ox))**2))*field)
    return result


def _verified_prediction(root, panel):
    prediction_path = _artifact(root, panel["prediction"])
    prediction = read_json(prediction_path)
    metadata = prediction.get("metadata", {})
    from .collect import collect_species
    actual = collect_species(root, metadata["species_id"], metadata["scenario_id"],
                             metadata["j_method"], exploratory=True, write=False)
    for key in ("atom_rows", "spectrum_atom_rows", "j_pairs"):
        if prediction.get(key) != actual.get(key):
            raise ValueError("Prediction numbers differ from independently re-collected ORCA results")
    for key in ("species_id", "compound_id", "scenario_id", "temperature_K", "solvent", "protocol_hash",
                "j_method", "shielding_method_fingerprint", "j_method_fingerprint", "weights",
                "reference_weights", "provenance", "ensemble_snapshot"):
        if metadata.get(key) != actual["metadata"].get(key):
            raise ValueError(f"Prediction provenance mismatch: {key}")
    frozen = metadata.get("frozen_ensembles", {})
    for name in ("target", "reference"):
        frozen_path = _artifact(root, frozen.get(name))
        if read_json(frozen_path) != metadata["ensemble_snapshot"][name]:
            raise ValueError("Frozen ensemble differs from prediction snapshot")
    return prediction, prediction_path


def _external_transitions(root, panel, prediction, nucleus, field):
    """Accept independently reviewed solver output, without claiming solver QA."""
    entry = panel["external_transitions"]
    path = _artifact(root, entry)
    data = read_json(path)
    _human_review(data.get("review"))
    expected = parameter_digest(prediction["spectrum_atom_rows"], prediction["j_pairs"], nucleus, field)
    if data.get("parameter_digest") != expected or data.get("prediction_sha256") != panel["prediction"]["sha256"]:
        raise ValueError("External transitions not bound to the actual prediction/field")
    for key in ("solver", "solver_version", "hamiltonian_description", "limitations"):
        if not isinstance(data.get(key), str) or not data[key].strip():
            raise ValueError("External solver metadata/limitations are required")
    if data.get("spin_model_scope") not in ("REVIEWED_FULL_SPIN_MODEL", "DECLARED_APPROXIMATION"):
        raise ValueError("External solver must declare full reviewed spin scope or an approximation")
    if (data["spin_model_scope"] == "REVIEWED_FULL_SPIN_MODEL"
            and any(r["element"] == "F" for r in prediction["spectrum_atom_rows"])):
        additional = data.get("additional_spin_parameters", {})
        # The current collector has H-H/H-F/C-F, NOT 19F shifts or F-F J.
        if not additional.get("description") or not additional.get("artifacts"):
            raise ValueError("A full F-containing spin model needs source-bound additional 19F shifts/F-F parameters or an explicit justified scope")
        for source in additional["artifacts"]:
            _artifact(root, source)
    sources = data.get("solver_artifacts", [])
    if not sources:
        raise ValueError("External solver input/output artifacts are required")
    for source in sources:
        _artifact(root, source)
    required = sorted(r["atom_index_zero_based"] for r in prediction["spectrum_atom_rows"]
                      if r["element"] == ("H" if nucleus == "1H" else "C"))
    if data.get("observed_atom_indices_zero_based") != required:
        raise ValueError("External transition observation scope is incomplete")
    lines = data.get("sticks", [])
    if not lines:
        raise ValueError("External solver contains no transitions")
    for line in lines:
        finite(line.get("ppm"), "external line ppm")
        _positive(line.get("intensity"), "external transition area", zero=True)
        ids = line.get("source_atom_indices")
        if not isinstance(ids, list) or not ids or any(type(i) is not int or i not in required for i in ids):
            raise ValueError("External transitions need valid observed-nucleus scope")
    if not math.isclose(sum(r["intensity"] for r in lines), len(required), abs_tol=1e-6):
        raise ValueError("External transition total area must equal the observed nucleus count")
    return lines, data


def validate_multiplets(root, panel_path, *, write=True, output_path=None):
    """Read experimental panel, validate provenance, report metrics for review.

    No scientific approval or replacement for real experimental validation is
    generated. Only corrupt/incomplete inputs are hard errors; discrepancies
    and model limitations generate recommendations, never alter original data.
    """
    root = Path(root).resolve(strict=True)
    path = Path(panel_path)
    if not path.is_absolute():
        path = inside(root, str(path))
    elif not path.resolve().is_relative_to(root):
        raise ValueError("Panel must be inside the resource package")
    panel = read_json(path)
    if panel.get("schema_version") != "1.0" or panel.get("data_kind") != "EXPERIMENTAL_NMR":
        raise ValueError("Versioned real experimental NMR panel required; no synthetic panel")
    _human_review(panel.get("assignment_review"))
    sources = panel.get("experimental_artifacts", [])
    if not sources:
        raise ValueError("Experimental source artifacts are required")
    for source in sources:
        _artifact(root, source)
    experiment = panel["experiment"]
    for key in ("sample_identity", "solvent", "reference", "processing", "decoupling"):
        if not isinstance(experiment.get(key), str) or not experiment[key].strip():
            raise ValueError(f"Experimental {key} must be declared")
    field = _positive(experiment.get("field_MHz"), "experimental field")
    temperature = _positive(experiment.get("temperature_K"), "experimental temperature")
    nucleus = experiment.get("nucleus")
    if nucleus not in ("1H", "13C"):
        raise ValueError("Nucleus must be 1H or 13C")
    if nucleus == "1H" and experiment["decoupling"] != "none":
        raise ValueError("Internal proton model expects no heteronuclear decoupling")
    if nucleus == "13C" and experiment["decoupling"] != "1H":
        raise ValueError("Internal carbon model expects proton decoupling only")
    prediction, prediction_path = _verified_prediction(root, panel)
    md = prediction["metadata"]
    if experiment["sample_identity"] not in (md["species_id"], md["compound_id"]):
        raise ValueError("Experimental sample and calculation identity disagree")
    settings = dict(panel.get("spectrum_settings", {}))
    nested_settings = settings.pop("spectrum", {})
    if not isinstance(nested_settings, dict) or "spectrum" in nested_settings:
        raise ValueError("Spectrum settings must contain at most one dictionary nesting layer")
    settings.update(nested_settings)
    # Field is always the experimental instrument field, not a fitted setting.
    settings["field_MHz_H" if nucleus == "1H" else "field_MHz_C"] = field
    settings.setdefault("linewidth_Hz_H" if nucleus == "1H" else "linewidth_Hz_C", 1.0 if nucleus == "1H" else 2.0)
    model = make_spectrum(prediction["spectrum_atom_rows"], prediction["j_pairs"], nucleus, settings)
    if model["field_MHz"] != field:
        raise ValueError("Spectrum field must match the experimental instrument field")
    lines = model["sticks"]
    external = None
    if panel.get("external_transitions"):
        lines, external = _external_transitions(root, panel, prediction, nucleus, field)
    tolerances = {**DEFAULT_TOLERANCES, **panel.get("tolerances", {})}
    if set(tolerances) != set(DEFAULT_TOLERANCES):
        raise ValueError("Unknown multiplet tolerance key")
    for key, value in tolerances.items():
        _positive(value, key, zero=True)
    observed_atoms = {r["atom_index_zero_based"]: r for r in prediction["spectrum_atom_rows"]
                      if r["element"] == ("H" if nucleus == "1H" else "C")}
    regions = panel.get("regions", [])
    if not isinstance(regions, list) or not regions:
        raise ValueError("At least one assigned region is required")
    reports, covered, ids = [], set(), set()
    all_warnings = list(model["warnings"])
    solvent = experiment["solvent"].lower().replace("_", "-")
    if solvent not in ("dmso", "dmso-d6") or abs(temperature-md["temperature_K"]) > 0.1:
        all_warnings.append("Experimental solvent/temperature differs from the calculated model; review transferability.")
    for region in regions:
        rid = region.get("id")
        if not isinstance(rid, str) or not rid.strip() or rid in ids:
            raise ValueError("Region IDs must be nonempty and unique")
        ids.add(rid)
        assigned = region.get("assigned_atom_indices_zero_based")
        if (not isinstance(assigned, list) or not assigned or any(type(i) is not int or i not in observed_atoms for i in assigned)
                or len(set(assigned)) != len(assigned) or covered.intersection(assigned)):
            raise ValueError("Regions need unique non-overlapping observed-atom assignments")
        covered.update(assigned)
        low, high = _range(region["experimental_window_ppm"], "experimental window")
        plow, phigh = _range(region.get("model_window_ppm", region["experimental_window_ppm"]), "model window")
        if not isinstance(region.get("overlap"), bool):
            raise ValueError("Region overlap must be explicitly declared")
        if region.get("second_order_observed") not in ("yes", "no", "unknown"):
            raise ValueError("Second-order observation must be yes/no/unknown")
        if not isinstance(region.get("multiplicity_description"), str) or not region["multiplicity_description"].strip():
            raise ValueError("Provide the experimental multiplicity description, including unresolved/overlap when appropriate")
        warnings = []
        observed = region.get("lines", [])
        if not observed and not region.get("unresolved_reason"):
            raise ValueError("Missing experimental lines require an explicit unresolved reason")
        for r in observed:
            if not low <= finite(r["ppm"], "observed ppm") <= high:
                raise ValueError("Experimental line lies outside its assigned window")
            _positive(r["area"], "observed component area")
        if observed and region.get("line_area_kind") != "fitted_or_assigned_component_area":
            raise ValueError("Line areas must be fitted/assigned integrated areas, not peak heights")
        if not isinstance(region.get("line_area_quantitative"), bool):
            raise ValueError("Declare whether relative line areas are quantitatively comparable")
        resolution = _positive(region.get("resolution_Hz"), "experimental resolution")
        floor = region.get("relative_area_floor", 0.01)
        selected = [r for r in lines if plow <= r["ppm"] <= phigh]
        resolved = resolve_transitions(selected, field, resolution, floor)
        metrics = compare_lines(observed, resolved, field)
        metrics["raw_model_window_area"] = sum(r["intensity"] for r in selected)
        metrics["resolved_model_area_after_detection_floor"] = sum(r["area"] for r in resolved)
        if not region["line_area_quantitative"]:
            metrics["relative_line_area_mae"] = None
            warnings.append("Relative line areas recorded but not scored: quantitative comparability was not established.")
        integral = region.get("integral", {})
        measured = _positive(integral.get("value"), "region integral", zero=True)
        per_nucleus = _positive(integral.get("normalization_per_nucleus"), "integral normalization")
        if not isinstance(integral.get("quantitative"), bool) or not isinstance(integral.get("normalization_description"), str):
            raise ValueError("Declare integral normalization and quantitative validity")
        expected = len(assigned)
        if not math.isclose(metrics["raw_model_window_area"], expected, rel_tol=0, abs_tol=0.05):
            warnings.append("Model-window transition area differs from assigned nucleus count: inspect truncation, shared blocks and contamination.")
        metrics["region_integral_in_nucleus_units"] = measured/per_nucleus
        metrics["expected_assigned_nucleus_count"] = expected
        metrics["relative_region_integral_error"] = abs(measured/per_nucleus-expected)/expected if integral["quantitative"] else None
        if not integral["quantitative"]:
            warnings.append("Integral recorded but not scored: acquisition/normalization is nonquantitative.")
        if nucleus == "13C":
            warnings.append("Calculated carbon site areas do not predict NOE/relaxation intensities; compare line areas only if experimentally justified.")
        if region["overlap"]:
            warnings.append("Overlapping region: line counts/areas/spacing require reviewer confirmation of deconvolution and assignments.")
        if not observed:
            warnings.append("Unresolved experimental region: splitting comparison remains incomplete, not passed.")
        if any(observed_atoms[i].get("exchangeable_H") for i in assigned):
            warnings.append("Exchangeable H present; exchange-dependent splitting is excluded by this model.")
        relevant = [b for b in model["block_report"] if set(b["atom_indices_zero_based"]).intersection(assigned)]
        approximate = any("first-order" in b["method"] for b in relevant)
        if approximate and (not external or external["spin_model_scope"] == "DECLARED_APPROXIMATION"):
            warnings.append("Internal model is first order for this block: full strong-coupling pattern validation is not established.")
        if external:
            warnings.append("External transitions used; solver spin-scope adequacy is independently reviewed, not certified by this comparison module.")
        ratios = []
        atom_map = {r["atom_index_zero_based"]: r for r in prediction["spectrum_atom_rows"]}
        for pair in prediction["j_pairs"]:
            i, j = pair["atom_i_zero_based"], pair["atom_j_zero_based"]
            if i not in assigned and j not in assigned:
                continue
            if atom_map[i]["element"] == atom_map[j]["element"] == "H" and pair["J_Hz"] != 0:
                ratio = abs((atom_map[i]["shift_ppm"]-atom_map[j]["shift_ppm"])*field/pair["J_Hz"])
                ratios.append({"pair_zero_based": [i, j], "delta_nu_over_abs_J": ratio})
        second_order_risk = any(r["delta_nu_over_abs_J"] < tolerances["strong_coupling_ratio_warning"] for r in ratios)
        if second_order_risk:
            warnings.append("Small delta-nu/|J| found: inspect second-order behavior; equivalent spins alone do not imply visible strong-coupling distortion.")
        outside_sources = sorted({i for r in selected for i in r.get("source_atom_indices", []) if i not in assigned})
        if outside_sources:
            warnings.append("Model window contains transitions sourced from other assigned nuclei or a shared exact block; review overlap/window completeness.")
        flags = [key for key in DEFAULT_TOLERANCES if key != "strong_coupling_ratio_warning"
                 and metrics.get(key) is not None and metrics[key] > tolerances[key]]
        reports.append({"id": rid, "assignment": region, "metrics": metrics, "resolved_model_lines": resolved,
                        "raw_model_transition_count": len(selected), "metric_flags": flags, "warnings": warnings,
                        "second_order_diagnostics": ratios, "internal_block_models": relevant,
                        "recommendation": "REVIEW_REQUIRED" if flags or warnings else "NUMERIC_AGREEMENT_REQUIRES_HUMAN_REVIEW"})
    excluded = panel.get("excluded_atoms", [])
    excluded_ids = set()
    for row in excluded:
        i = row.get("atom_index_zero_based")
        if type(i) is not int or i not in observed_atoms or i in excluded_ids or i in covered or not row.get("reason"):
            raise ValueError("Excluded atom needs a unique valid index and scientific reason")
        excluded_ids.add(i)
    if covered | excluded_ids != set(observed_atoms):
        raise ValueError("Every observed nucleus must be assigned or explicitly excluded")
    report = {"schema_version": "1.0", "status": "METRICS_COMPLETE_AWAITING_HUMAN_REVIEW",
              "notice": "A validation workflow/comparison report, NOT experimental agreement certification or automatic approval.",
              "panel": {"path": path.relative_to(root).as_posix(), "sha256": sha256_file(path)},
              "prediction": {"path": prediction_path.relative_to(root).as_posix(), "sha256": sha256_file(prediction_path)},
              "experimental_artifacts": sources, "experiment": experiment, "assignment_review": panel["assignment_review"],
              "tolerances_are_recommendations_only": tolerances, "regions": reports, "excluded_atoms": excluded,
              "coverage": {"observed_nuclei": len(observed_atoms), "assigned_nuclei": len(covered), "excluded_nuclei": len(excluded_ids)},
              "source_recollection": "ACTUAL_ORCA_OUTPUTS_RECOLLECTED_AND_MATCHED_READ_ONLY",
              "model": {k: model[k] for k in ("method", "field_MHz", "linewidth_FWHM_Hz", "threshold_Hz", "block_report")},
              "active_transition_source": "EXTERNALLY_REVIEWED_SOLVER_OUTPUT" if external else "INTERNAL_DECLARED_SPIN_MODEL",
              "external_solver": external, "warnings": all_warnings,
              "limitations": ["Line matching never fits or adjusts calculated shifts/J values.",
                              "Centered metrics assess pattern independent of a uniform shift error; absolute errors are retained.",
                              "Finite-resolution grouping uses area centroids, not an experimental lineshape fit.",
                              "Exchange kinetics, relaxation, quadrupolar couplings and 13C satellite patterns are not covered.",
                              "First-order/oversized/F-containing blocks require a reviewed adequate external spin model when strong coupling matters."]}
    if write:
        destination = inside(root, str(output_path)) if output_path else unique_output(root, "results/multiplet_validation") / "report.json"
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("x", encoding="utf-8") as handle:
            json.dump(report, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write("\n")
        report["report_path"] = destination.relative_to(root).as_posix()
    return report
