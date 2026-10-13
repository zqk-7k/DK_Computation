"""Collect one fixed-composition cluster ensemble; never pool solvent scenarios.

Both analyte and TMS use their own complete, source-validated Gibbs populations.
Exploratory output is explicitly unvalidated. Production output additionally
requires an independently reviewed experimental J-method gate. No QM is run.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import re
import uuid

from .orcaio import parse_output
from .common import _io_path

J_METHODS = ("PBE_pcJ3", "PBE0_pcJ2")


def read_json(path):
    return json.loads(_io_path(path).read_text(encoding="utf-8"))


def sha256_file(path):
    digest = hashlib.sha256()
    with _io_path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            digest.update(block)
    return digest.hexdigest()


def inside(root, relative):
    root = Path(root).resolve()
    if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
        raise ValueError("Expected a nonempty package-relative path")
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError("Path escapes package")
    return path


def finite(value, label):
    if isinstance(value, bool):
        raise ValueError(f"{label} must be finite numeric data")
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be finite numeric data") from exc
    if not math.isfinite(value):
        raise ValueError(f"{label} must be finite numeric data")
    return value


def identifier(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9_-]+", value):
        raise ValueError("Unsafe or missing identifier")
    return value


def _context(root):
    from .common import config, hash_protocol
    settings = config(root)
    return settings, hash_protocol(settings)


def validate_atom_map(rows):
    """Validate the FULL cluster graph; return unchanged global atom indices."""
    if not isinstance(rows, list) or not rows:
        raise ValueError("A full-cluster atom map is required")
    atoms = {}
    for expected, row in enumerate(rows):
        if not isinstance(row, dict) or type(row.get("index_zero_based")) is not int or row["index_zero_based"] != expected:
            raise ValueError("Atom map indices must preserve contiguous original order")
        if not re.fullmatch(r"[A-Z][a-z]?", row.get("element", "")):
            raise ValueError("Invalid atom element")
        if type(row.get("environment_atom")) is not bool or type(row.get("include_in_analyte_nmr")) is not bool:
            raise ValueError("Explicit environment and NMR flags are required")
        if not isinstance(row.get("component"), str) or not row["component"]:
            raise ValueError("Component assignment is required")
        if row["environment_atom"] == (row["component"] == "analyte"):
            raise ValueError("Environment/component assignment disagrees")
        iso = row.get("isotope", 0)
        if type(iso) is not int or iso < 0:
            raise ValueError("Invalid isotope annotation")
        observed = not row["environment_atom"] and (
            row["element"] == "H" and iso in (0, 1) or row["element"] == "C" and iso in (0, 13))
        if row["include_in_analyte_nmr"] != observed:
            raise ValueError("Analyte H/C observation flag disagrees with isotope/component")
        if row["element"] == "H" and type(row.get("exchangeable_H")) is not bool:
            raise ValueError("Hydrogen exchangeability must be explicit")
        atoms[expected] = row
    solute = [i for i, a in atoms.items() if not a["environment_atom"]]
    if not solute or solute != list(range(len(solute))):
        raise ValueError("This cluster implementation requires preserved analyte-first ordering")
    for i, row in atoms.items():
        ns = row.get("neighbours_zero_based")
        if not isinstance(ns, list) or any(type(n) is not int for n in ns) or len(ns) != len(set(ns)):
            raise ValueError("Invalid atom-map neighbors")
        for n in ns:
            if n not in atoms or n == i or i not in atoms[n].get("neighbours_zero_based", []):
                raise ValueError("Atom-map bonds must be symmetric and valid")
            if row["component"] != atoms[n]["component"]:
                raise ValueError("Frozen cluster map contains a cross-component bond")
    return atoms


def analyte_atoms(atoms):
    return {i: a for i, a in atoms.items() if not a["environment_atom"]}


def required_pairs(atoms):
    solute = analyte_atoms(atoms)
    observed = {i: a for i, a in solute.items() if a["element"] in ("C", "F") or
                a["element"] == "H" and a.get("isotope", 0) in (0, 1)}
    return {(i, j) for i in observed for j in observed if i < j and
            sorted((observed[i]["element"], observed[j]["element"])) in
            (["H", "H"], ["F", "H"], ["C", "F"])}


def weighted_average(values, weights):
    if not isinstance(values, (list, tuple)) or not values or len(values) != len(weights):
        raise ValueError("Complete matching conformer values and weights are required")
    ws = [finite(w, "weight") for w in weights]
    if any(w < 0 for w in ws) or not math.isclose(math.fsum(ws), 1, abs_tol=1e-8, rel_tol=0):
        raise ValueError("Weights must be nonnegative and already normalized")
    keys = set(values[0])
    if any(set(v) != keys for v in values):
        raise ValueError("Missing values or inconsistent conformer coverage; never fill with zero")
    return {k: math.fsum(w * finite(v[k], "prediction") for w, v in zip(ws, values)) for k in keys}


def extract_shieldings(parsed, atoms):
    needed = {i for i, a in atoms.items() if a["include_in_analyte_nmr"]}
    values = {}
    for row in parsed.get("shielding", []):
        i = row.get("atom_index_zero_based")
        if type(i) is not int or i not in needed or i in values or row.get("element") != atoms[i]["element"]:
            raise ValueError("Unexpected/duplicate shielding nucleus or solvent leakage")
        values[i] = finite(row["isotropic_ppm"], "shielding")
    if not needed or set(values) != needed:
        raise ValueError("Missing required analyte H/C shieldings")
    return values


def extract_couplings(parsed, metadata, atoms):
    declared = metadata.get("expected_j_pairs_zero_based")
    if not isinstance(declared, list):
        raise ValueError("Explicit exact J-pair whitelist is required")
    pairs = []
    for p in declared:
        if not isinstance(p, (list, tuple)) or len(p) != 2 or any(type(i) is not int for i in p):
            raise ValueError("Malformed J whitelist")
        pairs.append(tuple(sorted(p)))
    needed = required_pairs(atoms)
    if len(set(pairs)) != len(pairs) or set(pairs) != needed:
        raise ValueError("J whitelist must cover exactly the declared analyte HH/HF/CF pairs")
    values = {}
    for row in parsed.get("j_pairs", []):
        i, j = row.get("atom_i_zero_based"), row.get("atom_j_zero_based")
        if type(i) is not int or type(j) is not int:
            raise ValueError("Invalid J indices")
        pair = tuple(sorted((i, j)))
        if pair not in needed or pair in values:
            raise ValueError("Unexpected/duplicate J pair or environment leakage")
        values[pair] = finite(row["J_Hz"], "J")
    if set(values) != needed:
        raise ValueError("Missing J pairs; unavailable couplings cannot be set to zero")
    return values


def methyl_average(atoms, sigma, coupling):
    """Explicit fast-CH3-rotation assumption; never average CH2 or solvent nuclei."""
    solute = analyte_atoms(atoms)
    groups = []
    for a in solute.values():
        hs = sorted(i for i in a["neighbours_zero_based"] if solute[i]["element"] == "H"
                    and solute[i].get("isotope", 0) in (0, 1))
        if a["element"] == "C" and len(a["neighbours_zero_based"]) == 4 and len(hs) == 3:
            groups.append(hs)
    membership = {i: n for n, group in enumerate(groups) for i in group}
    shifts, js = dict(sigma), dict(coupling)
    for group in groups:
        if shifts:
            average = math.fsum(shifts[i] for i in group) / 3
            shifts.update({i: average for i in group})
    buckets = {}
    for pair, value in js.items():
        key = tuple(sorted((('m', membership[i]) if i in membership else ('a', i) for i in pair)))
        buckets.setdefault(key, []).append((pair, value))
    for (left, right), rows in buckets.items():
        expected = 3 if left == right and left[0] == 'm' else (
            (3 if left[0] == 'm' else 1) * (3 if right[0] == 'm' else 1))
        if len(rows) != expected:
            raise ValueError("Incomplete CH3 coupling symmetry orbit")
        average = math.fsum(v for _, v in rows) / len(rows)
        js.update({pair: average for pair, _ in rows})
    return shifts, js, groups


def _receipt(directory):
    receipt = read_json(directory / "receipt.json")
    if receipt.get("status") != "completed":
        raise ValueError("Completed runner receipt is required")
    for key in ("input_sha256", "output_sha256"):
        if not isinstance(receipt.get(key), dict) or not receipt[key]:
            raise ValueError("Runner receipt lacks artifact hashes")
        for name, expected in receipt[key].items():
            path = inside(directory, name)
            if not path.is_file() or sha256_file(path) != expected:
                raise ValueError("Runner receipt artifact changed")
    if not {"job.inp", "job.json"}.issubset(receipt["input_sha256"]) or "job.out" not in receipt["output_sha256"]:
        raise ValueError("Runner receipt does not bind input, metadata and output")
    for name in ("input.xyz", "input_structure.mol", "cluster.mol", "atom_map.json"):
        if (directory / name).exists() and name not in receipt["input_sha256"]:
            raise ValueError("Unbound structure or map in receipt")
    for name in ("job.property.txt", "propertyjob.property.txt"):
        if (directory / name).exists() and name not in receipt["output_sha256"]:
            raise ValueError("Unbound property output in receipt")
    return receipt


def _ensemble(root, species, scenario, settings, ph):
    from .pipeline import validate_ensemble_sources
    checked = validate_ensemble_sources(root, species, scenario)
    if not isinstance(checked, dict) or checked.get("accepted") is not True:
        raise ValueError("Source ensemble validation did not accept this ensemble")
    path = inside(root, f"ensembles/{ph}/{identifier(species)}/{identifier(scenario)}/ensemble.json")
    e = read_json(path)
    for key, expected in (("species_id", species), ("scenario_id", scenario), ("protocol_hash", ph)):
        if e.get(key) != expected:
            raise ValueError(f"Ensemble {key} mismatch")
    if e.get("solvent") != "DMSO":
        raise ValueError("This protocol only supports the locked DMSO scenarios")
    temperature = finite(e.get("temperature_K"), "temperature")
    if temperature <= 0 or not math.isclose(temperature, float(settings["temperature_K"]), abs_tol=1e-8):
        raise ValueError("Ensemble temperature differs from locked protocol")
    rows = e.get("conformers")
    if not isinstance(rows, list) or not rows:
        raise ValueError("Missing conformer ensemble")
    ids = [r.get("id") for r in rows]
    if len(set(ids)) != len(ids) or any(not isinstance(i, str) or not i for i in ids):
        raise ValueError("Missing/duplicate conformer identities")
    weights = [finite(r.get("weight"), "conformer weight") for r in rows]
    weighted_average([{} for _ in weights], weights)
    amap = inside(root, e["atom_map"])
    atoms = validate_atom_map(read_json(amap))
    return e, atoms, weights, {"ensemble_path": path.relative_to(root).as_posix(),
        "ensemble_sha256": sha256_file(path), "atom_map_sha256": sha256_file(amap)}


def _load_job(root, path, e, conformer, stage, settings, map_hash):
    from .gates import inspect_job
    directory = inside(root, path)
    receipt = _receipt(directory)
    meta = read_json(directory / "job.json")
    for key, expected in (("species_id", e["species_id"]), ("scenario_id", e["scenario_id"]),
                          ("protocol_hash", e["protocol_hash"]), ("solvent", e["solvent"]),
                          ("stage", stage), ("conformer_id", conformer["id"]), ("atom_map_sha256", map_hash)):
        if meta.get(key) != expected:
            raise ValueError(f"Property job {key} mismatch")
    fingerprint = meta.get("method_fingerprint")
    if not isinstance(fingerprint, str) or not fingerprint:
        raise ValueError("Property method fingerprint is required")
    check = inspect_job(directory, settings)
    if check.get("accepted") is not True:
        raise ValueError(f"Property job gate failed: {check.get('errors')}")
    text = (directory / "job.out").read_text(encoding="utf-8", errors="replace")
    prop = next((directory / n for n in ("job.property.txt", "propertyjob.property.txt") if (directory / n).exists()), None)
    parsed = parse_output(text, prop.read_text(encoding="utf-8", errors="replace") if prop else None)
    versions = re.findall(r"Program\s+Version\s+(\d+\.\d+(?:\.\d+)?)", text)
    if not parsed["normal_termination"] or not versions or versions[-1] != settings["orca_version"]:
        raise ValueError("Actual ORCA termination/version mismatch")
    record = {"directory": path, "stage": stage, "method_fingerprint": fingerprint,
              "job_json_sha256": sha256_file(directory / "job.json"), "output_sha256": sha256_file(directory / "job.out"),
              "receipt_sha256": sha256_file(directory / "receipt.json"), "actual_orca_version": versions[-1],
              "input_sha256": receipt["input_sha256"], "output_sha256_by_file": receipt["output_sha256"]}
    if prop:
        record["property_sha256"] = sha256_file(prop)
    return meta, parsed, record


def collect_candidate_couplings(root, species_id, scenario_id, j_method):
    """Read real completed J jobs only; shared by benchmark and final collection."""
    if j_method not in J_METHODS:
        raise ValueError("Unsupported J candidate")
    root = Path(root).resolve(strict=True)
    settings, ph = _context(root)
    e, atoms, weights, evidence = _ensemble(root, species_id, scenario_id, settings, ph)
    values, provenance, methods = [], [], set()
    for row in e["conformers"]:
        relative = row.get("j_jobs", {}).get(j_method)
        if relative is None:
            if required_pairs(atoms):
                raise ValueError("Missing J job for an analyte with required pairs")
            values.append({})
            continue
        meta, parsed, record = _load_job(root, relative, e, row, "j_" + j_method, settings, evidence["atom_map_sha256"])
        values.append(extract_couplings(parsed, meta, atoms))
        provenance.append(record)
        methods.add(meta["method_fingerprint"])
    if len(methods) > 1:
        raise ValueError("Mixed J methods within one conformer ensemble")
    raw = weighted_average(values, weights)
    _, js, groups = methyl_average(atoms, {}, raw)
    return {"ensemble": e, "atoms": atoms, "weights": weights, "couplings": js,
            "raw_couplings": raw, "method_fingerprint": next(iter(methods), None),
            "methyl_groups": groups, "provenance": provenance, "source": evidence}


def _shielding_ensemble(root, e, atoms, weights, settings, evidence):
    values, provenance, fingerprints = [], [], set()
    for row in e["conformers"]:
        meta, parsed, record = _load_job(root, row["shielding_job"], e, row, "shielding", settings, evidence["atom_map_sha256"])
        values.append(extract_shieldings(parsed, atoms))
        provenance.append(record)
        fingerprints.add(meta["method_fingerprint"])
    if len(fingerprints) != 1:
        raise ValueError("Mixed shielding methods within ensemble")
    return weighted_average(values, weights), values, provenance, next(iter(fingerprints))


def _validate_tms(atoms):
    solute = analyte_atoms(atoms)
    si = [i for i, a in solute.items() if a["element"] == "Si"]
    cs = {i for i, a in solute.items() if a["element"] == "C"}
    hs = {i for i, a in solute.items() if a["element"] == "H" and a.get("isotope", 0) in (0, 1)}
    if len(solute) != 17 or len(si) != 1 or len(cs) != 4 or len(hs) != 12:
        raise ValueError("Reference analyte is not explicit-H TMS")
    if set(solute[si[0]]["neighbours_zero_based"]) != cs:
        raise ValueError("TMS Si-C connectivity mismatch")
    for c in cs:
        ns = set(solute[c]["neighbours_zero_based"])
        if len(ns) != 4 or si[0] not in ns or not ns - {si[0]} <= hs:
            raise ValueError("TMS methyl connectivity mismatch")
    if any(len(solute[h]["neighbours_zero_based"]) != 1 or solute[h]["neighbours_zero_based"][0] not in cs for h in hs):
        raise ValueError("TMS hydrogen connectivity mismatch")


def _write_new(path, data):
    with _io_path(path).open("x", encoding="utf-8", newline="") as handle:
        handle.write(data)


def unique_output(root, relative):
    folder = inside(root, relative) / (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ") + "_" + uuid.uuid4().hex[:8])
    _io_path(folder).mkdir(parents=True, exist_ok=False)
    return folder


def collect_species(root, species_id, scenario_id, j_method=None, *, exploratory=False,
                    benchmark_gate=None, convergence_gate=None, write=True):
    root = Path(root).resolve(strict=True)
    from .policy import resolve_j_method
    j_method = resolve_j_method(root, species_id, j_method)
    candidate = collect_candidate_couplings(root, species_id, scenario_id, j_method)
    e, atoms, weights = candidate["ensemble"], candidate["atoms"], candidate["weights"]
    settings, ph = _context(root)
    ref_relative = f"ensembles/{ph}/TMS_DMSO/{identifier(scenario_id)}/ensemble.json"
    if e.get("reference_ensemble", ref_relative) != ref_relative:
        raise ValueError("Reference must use matching TMS scenario/protocol; no cross-composition reference")
    ref, ref_atoms, ref_weights, ref_evidence = _ensemble(root, "TMS_DMSO", scenario_id, settings, ph)
    _validate_tms(ref_atoms)
    raw_sigma, values, provenance, method = _shielding_ensemble(root, e, atoms, weights, settings, candidate["source"])
    ref_sigma, _, ref_provenance, ref_method = _shielding_ensemble(root, ref, ref_atoms, ref_weights, settings, ref_evidence)
    if ref_method != method or ref["temperature_K"] != e["temperature_K"]:
        raise ValueError("TMS method or temperature differs from analyte")
    reference = {el: math.fsum(v for i, v in ref_sigma.items() if ref_atoms[i]["element"] == el) /
                 sum(ref_atoms[i]["element"] == el for i in ref_sigma) for el in ("H", "C")}
    sigma, js, groups = methyl_average(atoms, raw_sigma, candidate["raw_couplings"])
    approval = convergence_approval = None
    if not exploratory:
        if not benchmark_gate:
            raise ValueError("Production collection requires independent reviewed J benchmark gate; use exploratory=True for unvalidated candidates")
        if not convergence_gate:
            raise ValueError("Production collection requires a source-bound target/reference sampling convergence gate")
        from .benchmark import validate_benchmark_gate
        pair_classes = {"".join(sorted(atoms[i]["element"] for i in pair)) for pair in candidate["couplings"]}
        approval = validate_benchmark_gate(root, benchmark_gate, e["compound_id"], j_method, scenario_id, e["temperature_K"], ph,
                                          target_pair_classes=pair_classes, method_fingerprint=candidate["method_fingerprint"])
        from .convergence import validate_gate
        convergence_approval = validate_gate(root, convergence_gate, species_id, scenario_id, j_method)
    atom_rows = []
    for i, value in sorted(sigma.items()):
        a = atoms[i]
        atom_rows.append({"atom_index_zero_based": i, "element": a["element"], "nucleus": "1H" if a["element"] == "H" else "13C",
            "shielding_ppm": value, "raw_weighted_shielding_ppm": raw_sigma[i], "reference_shielding_ppm": reference[a["element"]],
            "shift_ppm": reference[a["element"]] - value, "exchangeable_H": a.get("exchangeable_H", False),
            "sigma_weighted_sd_ppm": math.sqrt(math.fsum(w * (s[i] - raw_sigma[i]) ** 2 for w, s in zip(weights, values))),
            "neighbours_zero_based": a["neighbours_zero_based"]})
    j_rows = [{"atom_i_zero_based": i, "atom_j_zero_based": j, "J_Hz": v,
               "raw_weighted_J_Hz": candidate["raw_couplings"][(i, j)],
               "touches_exchangeable_H": bool(atoms[i].get("exchangeable_H") or atoms[j].get("exchangeable_H"))}
              for (i, j), v in sorted(js.items())]
    spectrum_atoms = atom_rows + [{"atom_index_zero_based": i, "element": "F", "nucleus": "19F", "shift_ppm": None,
                                  "exchangeable_H": False} for i, a in atoms.items() if not a["environment_atom"] and a["element"] == "F"]
    warnings = ["CALCULATED MODEL, NOT EXPERIMENTAL DATA; no accuracy guarantee.",
        "Each composition, species and J method remains a separate result; no pooled or spliced spectra.",
        "18O mass alone does not predict an experimental isotope chemical shift.",
        "Conformer SD is dispersion, not a prediction error bar; linewidths are display assumptions.",
        "Fast CH3 rotation assumed; CH2 protons never merged; exchangeable-H J excluded from coherent splitting."]
    if exploratory:
        warnings.insert(0, "UNVALIDATED_EXPLORATORY_CANDIDATE: no approved J-method selection or sampling-accuracy claim.")
    reviewed_approximations=[{'role':role,'conformer_id':row['id'],**row['scientific_review']}
                            for role,ensemble in [('target',e),('reference',ref)]
                            for row in ensemble['conformers'] if row.get('scientific_review')]
    if reviewed_approximations:
        warnings.append('Human-reviewed negative-frequency approximations are included. See reviewer, source hash and limitations; a true physical minimum is not certified.')
    metadata = {"status": "UNVALIDATED_EXPLORATORY_CANDIDATE" if exploratory else "CALCULATED_WITH_REVIEWED_J_SCOPE",
        "notice": "CALCULATED_NOT_EXPERIMENTAL", "species_id": species_id, "compound_id": e["compound_id"],
        "scenario_id": scenario_id, "solvent": e["solvent"], "temperature_K": e["temperature_K"], "protocol_hash": ph,
        "j_method": j_method, "shielding_method_fingerprint": method, "j_method_fingerprint": candidate["method_fingerprint"],
        "conformer_ids": [c["id"] for c in e["conformers"]], "weights": weights,
        "reference_conformer_ids": [c["id"] for c in ref["conformers"]], "reference_weights": ref_weights,
        "reference_shieldings_ppm": reference, "reference_model": "own complete TMS ensemble, same scenario/method/temperature",
        "source": candidate["source"], "reference_source": ref_evidence, "methyl_groups_zero_based": groups,
        "sampling_convergence_established": bool(e.get("sampling_convergence_established")),
        "reference_sampling_convergence_established": bool(ref.get("sampling_convergence_established")),
        "benchmark_approval": approval, "convergence_approval": convergence_approval,
        "human_reviewed_approximations":reviewed_approximations,
        "ensemble_snapshot": {"target": e, "reference": ref},
        "provenance": provenance + ref_provenance + candidate["provenance"]}
    result = {"atom_rows": atom_rows, "spectrum_atom_rows": spectrum_atoms, "j_pairs": j_rows, "warnings": warnings, "metadata": metadata}
    if write:
        mode = "exploratory" if exploratory else "reviewed_scope"
        destination = unique_output(root, f"results/{ph}/{identifier(species_id)}/{identifier(scenario_id)}/{j_method}/{mode}")
        result["output_directory"] = destination.relative_to(root).as_posix()
        result["result_path"] = (destination / "prediction.json").relative_to(root).as_posix()
        frozen = {}
        for label, evidence, snapshot in (("target", candidate["source"], e), ("reference", ref_evidence, ref)):
            source = inside(root, evidence["ensemble_path"])
            raw = _io_path(source).read_bytes()
            if hashlib.sha256(raw).hexdigest() != evidence["ensemble_sha256"] or json.loads(raw.decode("utf-8")) != snapshot:
                raise ValueError("Ensemble changed during collection; do not create an inconsistent snapshot")
            target = destination / (label + "_ensemble_snapshot.json")
            with _io_path(target).open("xb") as handle:
                handle.write(raw)
            frozen[label] = {"path": target.relative_to(root).as_posix(), "sha256": sha256_file(target)}
        metadata["frozen_ensembles"] = frozen
        _write_new(destination / "prediction.json", json.dumps(result, ensure_ascii=False, indent=2) + "\n")
        for name, rows in (("chemical_shifts.csv", atom_rows), ("j_couplings.csv", j_rows)):
            stream = io.StringIO(newline="")
            fields = list(rows[0]) if rows else ["atom_i_zero_based", "atom_j_zero_based", "J_Hz"]
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for row in rows:
                writer.writerow({k: json.dumps(v) if isinstance(v, list) else v for k, v in row.items()})
            _write_new(destination / name, stream.getvalue())
    return result
