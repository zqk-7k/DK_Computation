"""Compare real computed J branches with an independently reviewed experiment.

Input panel JSON contains data_origin='experimental', protocol_hash,
temperature_K, solvent, scenario_id, chemical_domain, independence_statement,
review={reviewer,reviewed_at,decision:'APPROVED_INDEPENDENT_EXPERIMENTAL_PANEL'},
predeclared_selection={criterion:'lowest_MAE_Hz',max_MAE_Hz,max_abs_error_Hz,
min_sign_accuracy,min_compounds,min_pairs}, and cases. Each case contains
case_id/compound_id/species_id/scenario_id/temperature_K, experimental_evidence
[{path,sha256,source_type:'acquired_experimental',source_description}],
assignment_review={reviewer,reviewed_at,decision:'APPROVED_ATOM_PAIR_ASSIGNMENT'},
and measured_pairs=[{atoms_zero_based:[i,j],J_Hz,sign_known:true|false,
evidence_index:0,source_locator:'page/table/row/assignment'}].

No predicted values may be pasted into the panel. Predictions are parsed from
completed, source-validated ORCA cluster jobs. The generated report NEVER grants
approval. A separate reviewed method-selection gate is needed for production.
"""
from __future__ import annotations

from datetime import datetime
import hashlib
import json
import math
from pathlib import Path

from .collect import (J_METHODS, _context, _write_new, collect_candidate_couplings,
                      finite, inside, read_json, sha256_file, unique_output)

REPORT_BODY_KEYS = ("scope", "panel_path", "panel_sha256", "panel_review", "independence_statement",
                    "predeclared_selection", "method_comparison", "recommended_methods", "measurements", "source_records")


def _review(record, decision):
    if not isinstance(record, dict) or record.get("decision") != decision:
        raise ValueError(f"Explicit human review {decision} is required")
    if not isinstance(record.get("reviewer"), str) or not record["reviewer"].strip():
        raise ValueError("Named reviewer is required; software cannot approve itself")
    if record["reviewer"].strip().lower() in {"auto", "automatic", "codex", "agent", "system", "unknown", "todo"}:
        raise ValueError("A software-generated reviewer is not scientific approval")
    try:
        datetime.fromisoformat(record["reviewed_at"].replace("Z", "+00:00"))
    except (KeyError, AttributeError, TypeError, ValueError) as exc:
        raise ValueError("A real ISO-format review date/time is required") from exc


def _pair_class(atoms, pair):
    return "".join(sorted(atoms[i]["element"] for i in pair))


def evaluate_comparison(measurements, policy):
    """Pure comparison; signed pairs compare signed J, unsigned pairs compare |J|."""
    required = ("max_MAE_Hz", "max_abs_error_Hz", "min_sign_accuracy", "min_compounds", "min_pairs")
    if not isinstance(policy, dict) or policy.get("criterion") != "lowest_MAE_Hz" or any(k not in policy for k in required):
        raise ValueError("Predeclared explicit thresholds and lowest_MAE_Hz criterion are required")
    mae_limit = finite(policy["max_MAE_Hz"], "MAE threshold")
    maximum_limit = finite(policy["max_abs_error_Hz"], "maximum-error threshold")
    sign_limit = finite(policy["min_sign_accuracy"], "sign-accuracy threshold")
    for key, minimum in (("min_compounds", 2), ("min_pairs", 2)):
        if type(policy[key]) is not int or policy[key] < minimum:
            raise ValueError("Cross-compound benchmark requires at least two compounds and matched pairs")
    if mae_limit < 0 or maximum_limit < 0 or not 0 <= sign_limit <= 1:
        raise ValueError("Invalid predeclared validation thresholds")
    if len(measurements) < policy["min_pairs"] or len({r["compound_id"] for r in measurements}) < policy["min_compounds"]:
        raise ValueError("Independent benchmark panel has insufficient compound/pair coverage")
    results = {}
    for method in J_METHODS:
        errors, signed_errors, signs = [], [], []
        for row in measurements:
            exp = finite(row["experimental_J_Hz"], "experimental J")
            predicted = finite(row["predicted_J_Hz"][method], "predicted J")
            if type(row.get("sign_known")) is not bool:
                raise ValueError("Experimental sign-known status must be explicit")
            if row["sign_known"]:
                delta = predicted - exp
                signed_errors.append(delta)
                if exp != 0:
                    signs.append((predicted > 0) == (exp > 0) and predicted != 0)
            else:
                delta = abs(predicted) - abs(exp)
            errors.append(delta)
        mae = math.fsum(abs(x) for x in errors) / len(errors)
        maximum = max(abs(x) for x in errors)
        sign_accuracy = math.fsum(signs) / len(signs) if signs else None
        sign_pass = sign_limit == 0 or sign_accuracy is not None and sign_accuracy >= sign_limit
        results[method] = {"MAE_Hz": mae, "max_abs_error_Hz": maximum,
            "mean_comparison_error_Hz": math.fsum(errors) / len(errors),
            "mean_signed_error_Hz": math.fsum(signed_errors) / len(signed_errors) if signed_errors else None,
            "sign_accuracy": sign_accuracy, "sign_evaluated_pairs": len(signs),
            "pairs": len(errors), "per_pair_errors_Hz": errors,
            "passes_predeclared_numeric_thresholds": mae <= mae_limit and maximum <= maximum_limit and sign_pass}
    passing = [m for m in J_METHODS if results[m]["passes_predeclared_numeric_thresholds"]]
    best = min((results[m]["MAE_Hz"] for m in passing), default=None)
    recommended = [m for m in passing if math.isclose(results[m]["MAE_Hz"], best, abs_tol=1e-12)]
    return results, recommended


def _evidence(root, records):
    if not isinstance(records, list) or not records:
        raise ValueError("Acquired experimental source artifacts are required")
    evidence = []
    for entry in records:
        if entry.get("source_type") != "acquired_experimental" or not str(entry.get("source_description", "")).strip():
            raise ValueError("Simulated/unknown data cannot be the experimental benchmark")
        path = inside(root, entry.get("path"))
        # Never relabel this workflow's predictions or QM outputs as experiment.
        if path.relative_to(root).parts[0] in {"results", "jobs", "ensembles", "01_goat_explicit"}:
            raise ValueError("Calculated workflow artifacts cannot supply experimental J observations")
        if not path.is_file() or sha256_file(path) != entry.get("sha256"):
            raise ValueError("Experimental evidence artifact is missing or changed")
        if path.suffix.lower() in {".txt", ".csv", ".tsv", ".json"}:
            header = path.read_text(encoding="utf-8-sig", errors="replace")[:8192].upper()
            if any(token in header for token in ("CALCULATED_NOT_EXPERIMENTAL", "NOT EXPERIMENTAL DATA",
                                                 "SYNTHETIC UNIT TEST", "SYNTHETIC ENGINEERING", "SIMULATED DATA")):
                raise ValueError("Experimental artifact explicitly declares simulated/synthetic data")
        if path.suffix.lower() == ".json":
            payload = read_json(path)
            if isinstance(payload, dict):
                markers = [payload.get("data_origin", ""), payload.get("notice", ""), payload.get("status", "")]
                markers += [payload.get("metadata", {}).get(k, "") for k in ("notice", "status") if isinstance(payload.get("metadata", {}), dict)]
                if any(any(token in str(v).upper() for token in ("SIMULAT", "SYNTHETIC", "CALCULATED", "EXPLORATORY")) for v in markers):
                    raise ValueError("Experimental artifact declares calculated/synthetic data")
        evidence.append({"path": path.relative_to(root).as_posix(), "sha256": sha256_file(path),
                         "source_type": entry["source_type"], "source_description": entry["source_description"]})
    return evidence


def run_benchmark(root, panel_path, target_compounds, *, output_path=None, write=True):
    root = Path(root).resolve(strict=True)
    path = Path(panel_path)
    if path.is_absolute():
        path = path.resolve()
        if not path.is_relative_to(root):
            raise ValueError("Copy the experimental panel into the package before binding provenance")
    else:
        path = inside(root, str(path))
    panel = read_json(path)
    if panel.get("data_origin") != "experimental" or panel.get("simulated", False):
        raise ValueError("Only independently acquired experimental J panels are accepted")
    if any(k in panel for k in ("predictions", "predicted_J_Hz", "fitted_target_peaks")):
        raise ValueError("Predictions must come from real completed ORCA jobs, not panel entries")
    _review(panel.get("review"), "APPROVED_INDEPENDENT_EXPERIMENTAL_PANEL")
    if not str(panel.get("independence_statement", "")).strip() or not str(panel.get("chemical_domain", "")).strip():
        raise ValueError("Independent origin and predeclared chemical-domain scope must be stated")
    targets = sorted(set(target_compounds or []))
    if not targets or any(not isinstance(t, str) or not t.strip() for t in targets):
        raise ValueError("Explicit prediction targets are needed to exclude target-self validation")
    settings, ph = _context(root)
    if panel.get("protocol_hash") != ph or panel.get("solvent") != "DMSO":
        raise ValueError("Experimental panel must target the exact locked protocol/solvent")
    temperature = finite(panel.get("temperature_K"), "panel temperature")
    if not math.isclose(temperature, settings["temperature_K"], abs_tol=1e-8, rel_tol=0):
        raise ValueError("Experimental and calculated temperatures differ")
    scenario = panel.get("scenario_id")
    if not isinstance(scenario, str) or not scenario:
        raise ValueError("One fixed composition scenario is required per benchmark")
    cases = panel.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("No independent experimental cases")
    measurements, source_records, case_ids, case_compounds = [], [], set(), set()
    fingerprints = {m: set() for m in J_METHODS}
    pair_classes = set()
    for case in cases:
        cid, compound = case.get("case_id"), case.get("compound_id")
        if not isinstance(cid, str) or not cid or cid in case_ids:
            raise ValueError("Missing/duplicate benchmark case ID")
        case_ids.add(cid)
        if not isinstance(compound, str) or not compound or compound in targets:
            raise ValueError("Benchmark must be independent of every prediction target")
        if compound in case_compounds:
            raise ValueError("Use one case per independent compound; repeated spectra are not independent compounds")
        case_compounds.add(compound)
        if case.get("data_origin", "experimental") != "experimental" or case.get("simulated", False):
            raise ValueError("Simulated case rejected")
        if any(k in case for k in ("predictions", "predicted_J_Hz")):
            raise ValueError("Pasted predicted J data rejected")
        if case.get("scenario_id") != scenario or not math.isclose(finite(case.get("temperature_K"), "case temperature"), temperature, abs_tol=1e-8, rel_tol=0):
            raise ValueError("Case composition/temperature mismatch")
        _review(case.get("assignment_review"), "APPROVED_ATOM_PAIR_ASSIGNMENT")
        evidence = _evidence(root, case.get("experimental_evidence"))
        candidates = {m: collect_candidate_couplings(root, case.get("species_id"), scenario, m) for m in J_METHODS}
        for method, candidate in candidates.items():
            e = candidate["ensemble"]
            if e.get("compound_id") != compound or e.get("protocol_hash") != ph or e.get("temperature_K") != temperature:
                raise ValueError("Calculated case identity/protocol/temperature does not match experimental assignment")
            if not candidate["method_fingerprint"]:
                raise ValueError("J benchmark requires actual computed J pairs")
            fingerprints[method].add(candidate["method_fingerprint"])
        if candidates[J_METHODS[0]]["source"] != candidates[J_METHODS[1]]["source"]:
            raise ValueError("Both J methods must use the same reviewed geometry/population ensemble")
        pairs = case.get("measured_pairs")
        if not isinstance(pairs, list) or not pairs:
            raise ValueError("No matched experimental pairs in case")
        seen = set()
        for entry in pairs:
            pair = entry.get("atoms_zero_based")
            if not isinstance(pair, list) or len(pair) != 2 or any(type(i) is not int for i in pair):
                raise ValueError("Experimental assignments require two explicit zero-based atom indices")
            pair = tuple(sorted(pair))
            if pair[0] == pair[1] or pair in seen:
                raise ValueError("Duplicate/invalid experimental atom pair")
            seen.add(pair)
            if any(pair not in c["couplings"] for c in candidates.values()):
                raise ValueError("Experimental pair must match computed analyte pairs in BOTH methods")
            idx = entry.get("evidence_index")
            if type(idx) is not int or not 0 <= idx < len(evidence) or not str(entry.get("source_locator", "")).strip():
                raise ValueError("Every experimental assignment needs an artifact and exact source locator")
            if type(entry.get("sign_known")) is not bool:
                raise ValueError("Do not infer a coupling sign that was not measured")
            atoms = candidates[J_METHODS[0]]["atoms"]
            if any(atoms[i].get("exchangeable_H", False) for i in pair):
                raise ValueError("Exchangeable-H J needs a separately commissioned validation domain")
            pair_class = _pair_class(atoms, pair)
            pair_classes.add(pair_class)
            measurements.append({"case_id": cid, "compound_id": compound, "species_id": case["species_id"],
                "atoms_zero_based": list(pair), "pair_class": pair_class,
                "experimental_J_Hz": finite(entry.get("J_Hz"), "experimental J"), "sign_known": entry["sign_known"],
                "predicted_J_Hz": {m: c["couplings"][pair] for m, c in candidates.items()},
                "experimental_artifact": evidence[idx], "source_locator": entry["source_locator"]})
        source_records.append({"case_id": cid, "assignment_review": case["assignment_review"],
            "experimental_evidence": evidence, "methods": {m: {"source": c["source"], "provenance": c["provenance"],
            "conformer_ids": [r["id"] for r in c["ensemble"]["conformers"]], "weights": c["weights"]} for m, c in candidates.items()}})
    if any(len(fs) != 1 for fs in fingerprints.values()):
        raise ValueError("Benchmark contains mixed method/basis fingerprints; split the validation domains")
    comparison, recommended = evaluate_comparison(measurements, panel.get("predeclared_selection"))
    body = {"scope": {"protocol_hash": ph, "solvent": "DMSO", "scenario_id": scenario, "temperature_K": temperature,
            "chemical_domain": panel["chemical_domain"], "target_compounds": targets, "benchmark_compounds": sorted(case_compounds),
            "pair_classes": sorted(pair_classes), "method_fingerprints": {m: next(iter(fs)) for m, fs in fingerprints.items()}},
        "panel_path": path.relative_to(root).as_posix(), "panel_sha256": sha256_file(path),
        "panel_review": panel["review"], "independence_statement": panel["independence_statement"],
        "predeclared_selection": panel["predeclared_selection"], "method_comparison": comparison,
        "recommended_methods": recommended, "measurements": measurements, "source_records": source_records}
    digest = hashlib.sha256(json.dumps(body, sort_keys=True, ensure_ascii=False).encode()).hexdigest()
    result = {**body, "scientific_digest": digest,
        "status": "EVALUATED_REQUIRES_SELECTION_REVIEW" if recommended else "NO_CANDIDATE_PASSES_DECLARED_THRESHOLDS",
        "experimental_origin_status": "USER_ATTESTED_ACQUIRED_EXPERIMENTAL_WITH_REVIEW_AND_ARTIFACT_HASHES",
        "approval_granted": False,
        "limitations": ["Artifact hashes establish provenance, not experimental authenticity or independent scientific correctness.",
            "This comparison does not validate chemical shifts, all compound classes, all couplings, or sampling convergence.",
            "Numerical J agreement is not a multiplet-topology benchmark; both require suitable experimental assignments.",
            "A separate human method-selection gate is mandatory; do not auto-approve this report."]}
    if write:
        if output_path is None:
            destination = unique_output(root, f"benchmarks/{ph}/{scenario}") / "comparison.json"
        else:
            destination = inside(root, str(output_path))
            destination.parent.mkdir(parents=True, exist_ok=True)
        result["result_path"] = destination.relative_to(root).as_posix()
        _write_new(destination, json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    return result


def validate_benchmark_gate(root, gate_path, target_compound, j_method, scenario_id, temperature_K,
                            protocol_hash, *, target_pair_classes=None, method_fingerprint=None):
    root = Path(root).resolve(strict=True)
    path = inside(root, str(gate_path))
    gate = read_json(path)
    _review(gate, "APPROVED_FOR_STATED_SCOPE")
    if not str(gate.get("selection_rationale", "")).strip():
        raise ValueError("Human method selection requires a stated scientific rationale")
    report_path = inside(root, gate.get("report_path"))
    if sha256_file(report_path) != gate.get("report_sha256"):
        raise ValueError("Benchmark report missing/changed")
    report = read_json(report_path)
    try:
        bound_body = {k: report[k] for k in REPORT_BODY_KEYS}
    except KeyError as exc:
        raise ValueError("Incomplete experimental benchmark report") from exc
    bound_digest = hashlib.sha256(json.dumps(bound_body, sort_keys=True, ensure_ascii=False).encode()).hexdigest()
    if bound_digest != report.get("scientific_digest"):
        raise ValueError("Reviewed benchmark report does not match its scientific content digest")
    scope = report["scope"]
    if gate.get("scope") != scope:
        raise ValueError("Reviewed gate scope differs from benchmark scope")
    for key, expected in (("protocol_hash", protocol_hash), ("scenario_id", scenario_id), ("temperature_K", temperature_K)):
        if scope.get(key) != expected:
            raise ValueError(f"J gate {key} mismatch")
    if target_compound not in scope["target_compounds"] or target_compound in scope["benchmark_compounds"]:
        raise ValueError("Target is outside independent validation scope or was used for self-validation")
    if j_method not in J_METHODS or gate.get("selected_j_method") != j_method:
        raise ValueError("J branch does not match reviewed selection")
    if target_pair_classes is not None and not set(target_pair_classes).issubset(scope["pair_classes"]):
        raise ValueError("Target contains coupling classes not covered by experimental validation")
    if method_fingerprint is not None and scope["method_fingerprints"].get(j_method) != method_fingerprint:
        raise ValueError("Target method/basis differs from independently validated fingerprint")
    fresh = run_benchmark(root, report["panel_path"], scope["target_compounds"], write=False)
    if fresh["scientific_digest"] != report.get("scientific_digest"):
        raise ValueError("Benchmark experimental/computational evidence changed; repeat review")
    if j_method not in fresh["recommended_methods"]:
        raise ValueError("Selected J method does not pass the predeclared method-selection rule")
    return {"gate_path": path.relative_to(root).as_posix(), "gate_sha256": sha256_file(path),
            "report_path": report_path.relative_to(root).as_posix(), "report_sha256": sha256_file(report_path),
            "reviewer": gate["reviewer"], "reviewed_at": gate["reviewed_at"], "selected_j_method": j_method,
            "scope": scope, "status": "REVIEWED_J_METHOD_FOR_STATED_SCOPE_ONLY"}
