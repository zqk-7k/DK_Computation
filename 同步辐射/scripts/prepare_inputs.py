#!/usr/bin/env python3
"""Prepare independent, geometry-based FEFF9/10 jobs from the supplied CIFs.

This module does not run DFT/FEFF, invent electrode potentials, remove H, reduce
cells, or infer oxidation states. Site indices are zero-based pymatgen indices;
the exact site order, CIF labels, coordinates and periodic images are recorded.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import inspect
import json
import math
from pathlib import Path
import re
import sys
import warnings


CORE_MODELS = [
    "CoFe-LDH.cif",
    "Top-O1.cif",
    "Top-O2.cif",
    "W1-path/CH2OH-CH2OH.cif",
    "W1-path/CH2OH-CHO.cif",
    "CoFe-LDH-path/CH2OH-CH2OH.cif",
    "CoFe-LDH-path/CH2OH-CHO.cif",
]
EDGES = {"W": "L3", "Co": "K", "Fe": "K"}
DEFAULT_FEFF = {
    "cluster_radius_angstrom": 8.0,
    "scf_radius_angstrom": 5.0,
    "fms_radius_angstrom": 6.0,
    "rpath_angstrom": 5.0,
    "exafs_kmax": 14.0,
    "xanes_card": [4.0, 0.05, 0.2],
    "scf_card_suffix": [0, 30, 0.2, 1],
    "fms_card_suffix": [0],
    "corehole": "FSR",
    "s02": 1.0,
    "nleg": 4,
    "print": [1, 0, 0, 0, 0, 3],
    "modes": ["xanes", "exafs"],
}


def fail(message: str) -> None:
    raise ValueError(message)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative_child(root: Path, raw: str, field: str) -> Path:
    """Keep portable package paths and prevent writing through an outside path."""
    p = (root / raw).resolve()
    if not p.is_relative_to(root.resolve()):
        fail(f"{field} must stay inside the workflow package: {raw}")
    return p


def number(value) -> str:
    x = float(value)
    if not math.isfinite(x):
        fail("Non-finite FEFF parameter")
    return f"{x:.10g}"


def tokens(value, name: str) -> str:
    values = value.split() if isinstance(value, str) else value
    if not isinstance(values, (list, tuple)) or not values:
        fail(f"{name} must be a nonempty numeric list or numeric string")
    return " ".join(number(v) for v in values)


def settings(config: dict, profile: str) -> tuple[dict, dict]:
    profiles = config.get("profiles", {})
    if not re.fullmatch(r"[A-Za-z0-9_-]+", profile):
        fail("Profile names may contain letters, digits, underscores and hyphens only")
    if profile not in {"pilot", "core", "full"} and profile not in profiles:
        fail(f"Unknown profile {profile}; define it in profiles first")
    selected = profiles.get(profile, {})
    if isinstance(selected, list):
        selected = {"models": selected}
    if not isinstance(selected, dict):
        fail(f"profiles.{profile} must be an object")
    feff = dict(DEFAULT_FEFF)
    feff.update(config.get("feff", {}))
    feff.update(selected.get("feff", {}))
    for key in ("cluster_radius_angstrom", "scf_radius_angstrom",
                "fms_radius_angstrom", "rpath_angstrom", "exafs_kmax"):
        value = float(feff[key])
        if not math.isfinite(value) or value <= 0:
            fail(f"feff.{key} must be finite and positive")
        feff[key] = value
    if feff["cluster_radius_angstrom"] < max(
        feff["scf_radius_angstrom"], feff["fms_radius_angstrom"], feff["rpath_angstrom"]
    ):
        fail("The atom-cluster radius must cover SCF, FMS and EXAFS path radii")
    if feff["scf_radius_angstrom"] > feff["fms_radius_angstrom"]:
        fail("Use SCF radius <= FMS radius for this baseline workflow")
    if not set(feff["modes"]).issubset({"xanes", "exafs"}) or not feff["modes"]:
        fail("feff.modes must contain xanes and/or exafs")
    if len(set(feff["modes"])) != len(feff["modes"]):
        fail("Duplicate FEFF modes")
    if str(feff["corehole"]).upper() not in {"FSR", "REGULAR", "RPA", "NONE"}:
        fail("corehole must be FSR, regular, RPA or none")
    if not 0 < float(feff["s02"]) <= 1:
        fail("s02 must be in (0, 1]; default 1 is an uncalibrated theory baseline")
    if int(feff["nleg"]) != float(feff["nleg"]) or int(feff["nleg"]) < 2:
        fail("nleg must be an integer >= 2")
    return selected, feff


def select_models(structures_dir: Path, selected: dict, profile: str) -> list[Path]:
    names = selected.get("models")
    if names == "all":
        names = [p.relative_to(structures_dir).as_posix()
                 for p in sorted(structures_dir.rglob("*.cif"))]
    if names is None:
        if profile == "pilot":
            names = ["Top-O1.cif"]
        elif profile == "core":
            names = CORE_MODELS
        elif profile == "full":
            names = [p.relative_to(structures_dir).as_posix()
                     for p in sorted(structures_dir.rglob("*.cif"))]
        else:
            fail(f"Custom profile {profile} needs an explicit models list")
    if not isinstance(names, list) or not names:
        fail(f"No model list for {profile}")
    paths = [relative_child(structures_dir, str(n), "model") for n in names]
    if len(set(paths)) != len(paths):
        fail("Duplicate CIF in selected model list")
    for p in paths:
        if p.suffix.lower() != ".cif" or not p.is_file():
            fail(f"CIF not found: {p}")
    return paths


def load_structure(path: Path):
    try:
        from pymatgen.io.cif import CifParser
    except ImportError as exc:
        raise RuntimeError(
            "pymatgen is required; activate the package environment before preparing inputs."
        ) from exc
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        # No rational-coordinate snapping, occupancy normalization or primitive cell.
        parser = CifParser(str(path), occupancy_tolerance=1.0,
                           site_tolerance=1e-5, frac_tolerance=0.0)
        kwargs = {"primitive": False, "symmetrized": False, "check_occu": True}
        if "on_error" in inspect.signature(parser.parse_structures).parameters:
            kwargs["on_error"] = "raise"
        structures = parser.parse_structures(**kwargs)
    if len(structures) != 1:
        fail(f"Expected exactly one structure in {path}, found {len(structures)}")
    structure = structures[0]
    if not structure.is_ordered:
        fail(f"Partial occupancy/disorder in {path}; explicit configurations are required. "
             "Occupancies will not be rounded or silently converted to full occupancy.")
    for i, site in enumerate(structure):
        total = sum(float(v) for v in site.species.values())
        if abs(total - 1.0) > 1e-8:
            fail(f"Site {i} has occupancy {total} in {path}; resolve explicitly")
    parsed_warnings = [str(w.message) for w in captured]
    return structure, list(dict.fromkeys(parsed_warnings))


def build_cluster(structure, absorber: int, radius: float) -> list[dict]:
    origin = structure[absorber]
    atoms = [{"site_index": absorber, "cif_label": str(origin.label),
              "element": origin.specie.symbol, "Z": int(origin.specie.Z),
              "image": [0, 0, 0], "xyz_angstrom": [0.0, 0.0, 0.0],
              "distance_angstrom": 0.0, "absorber": True}]
    seen = {(absorber, (0, 0, 0))}
    neighbors = structure.get_sites_in_sphere(origin.coords, radius,
                                             include_index=True, include_image=True)
    for neighbor in neighbors:
        idx = int(neighbor.index)
        image = tuple(int(round(float(v))) for v in neighbor.image)
        key = (idx, image)
        if key in seen:
            continue
        seen.add(key)
        xyz = [float(x) for x in neighbor.coords - origin.coords]
        distance = math.sqrt(sum(x*x for x in xyz))
        if distance < 1e-6:
            fail(f"Overlapping atoms at absorber {absorber}, site {idx}, image {image}")
        if distance > radius + 1e-6:
            continue
        site = structure[idx]
        atoms.append({"site_index": idx, "cif_label": str(site.label),
                      "element": site.specie.symbol, "Z": int(site.specie.Z),
                      "image": list(image), "xyz_angstrom": xyz,
                      "distance_angstrom": distance, "absorber": False})
    atoms[1:] = sorted(atoms[1:], key=lambda a: (
        round(a["distance_angstrom"], 10), a["site_index"], tuple(a["image"])))
    return atoms


def make_input(structure, atoms: list[dict], model: str, edge: str,
               absorber: int, mode: str, feff: dict) -> str:
    counts = Counter(site.specie.symbol for site in structure)
    central_element = atoms[0]["element"]
    # Non-central potentials are grouped by element. This approximation must be
    # checked for chemically distinct O and metal environments during convergence.
    scatterer_elements = sorted({a["element"] for a in atoms[1:]},
                                 key=lambda e: next(a["Z"] for a in atoms if a["element"] == e))
    ipots = {element: i+1 for i, element in enumerate(scatterer_elements)}
    absorber_stoich = 0.01 if central_element in ipots else float(counts[central_element])
    lines = [
        f"TITLE Geometry-based theoretical XAS: {model}",
        f"TITLE absorber_index_zero_based={absorber}; edge={central_element}_{edge}; mode={mode}",
        "* No electrode potential, DFT charge density or magnetic order inferred from CIF.",
        "* All hydrogen atoms retained. Atomic coordinates are in angstrom.",
        f"EDGE {edge}",
        f"S02 {number(feff['s02'])}",
        f"COREHOLE {str(feff['corehole'])}",
        "CONTROL 1 1 1 1 1 1",
        "PRINT " + tokens(feff["print"], "print"),
        f"SCF {number(feff['scf_radius_angstrom'])} " + tokens(feff["scf_card_suffix"], "scf_card_suffix"),
    ]
    if mode == "xanes":
        lines += ["XANES " + tokens(feff["xanes_card"], "xanes_card"),
                  f"FMS {number(feff['fms_radius_angstrom'])} " + tokens(feff["fms_card_suffix"], "fms_card_suffix"),
                  "RPATH -1"]
    else:
        lines += [f"EXAFS {number(feff['exafs_kmax'])}",
                  f"RPATH {number(feff['rpath_angstrom'])}",
                  f"NLEG {int(feff['nleg'])}"]
    lines += ["", "POTENTIALS", "* ipot Z tag lmax1 lmax2 xnatph",
              f"0 {atoms[0]['Z']} {central_element} -1 -1 {number(absorber_stoich)}"]
    for symbol, ipot in ipots.items():
        z = next(a["Z"] for a in atoms if a["element"] == symbol)
        stoich = float(counts[symbol]) - (absorber_stoich if symbol == central_element else 0.0)
        lines.append(f"{ipot} {z} {symbol} -1 -1 {number(stoich)}")
    lines += ["", "ATOMS", "* x y z ipot tag distance"]
    for index, atom in enumerate(atoms):
        ipot = 0 if index == 0 else ipots[atom["element"]]
        atom["ipot"] = ipot
        tag = f"{atom['element']}{atom['site_index']:04d}"[:6]
        x, y, z = atom["xyz_angstrom"]
        lines.append(f"{x: .9f} {y: .9f} {z: .9f} {ipot:3d} {tag:6s} {atom['distance_angstrom']:.9f}")
    lines += ["END", ""]
    return "\n".join(lines)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8", newline="\n")
    temporary.replace(path)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--profile", required=True,
                    help="pilot, core, full, or a named profile defined in config.json")
    args = ap.parse_args(argv)
    config_path = args.config.resolve()
    root = config_path.parent
    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    selected, feff = settings(config, args.profile)
    selected_edges = selected.get("edges", ["W_L3"] if args.profile == "pilot" else
                                  [f"{e}_{edge}" for e, edge in EDGES.items()])
    allowed_edges = {f"{e}_{edge}" for e, edge in EDGES.items()}
    if (not isinstance(selected_edges, list) or not selected_edges or
            not set(selected_edges).issubset(allowed_edges)):
        fail(f"Profile edges must be selected from {sorted(allowed_edges)}")
    site_selection = selected.get("site_selection", "all")
    if site_selection not in {"all", "representative_nearest_w"}:
        fail("site_selection must be all or representative_nearest_w")
    paths = config.get("paths", {})
    structures_dir = relative_child(root, paths.get("structures_dir", "structures"), "structures_dir")
    jobs_dir = relative_child(root, paths.get("jobs_dir", "jobs"), "jobs_dir")
    manifests_dir = relative_child(root, paths.get("manifests_dir", "manifests"), "manifests_dir")
    sources = select_models(structures_dir, selected, args.profile)
    jobs, prepared, model_ids = [], [], set()
    for source in sources:
        rel_model = source.relative_to(structures_dir).as_posix()
        model = re.sub(r"[^A-Za-z0-9_.-]+", "_", rel_model[:-4].replace("/", "__"))
        if model in model_ids:
            fail(f"Model ID collision for {source}; rename the copied CIFs explicitly")
        model_ids.add(model)
        structure, parser_warnings = load_structure(source)
        counts = Counter(site.specie.symbol for site in structure)
        selected_sites = {e: [i for i, s in enumerate(structure) if s.specie.symbol == e]
                          for e, edge in EDGES.items() if f"{e}_{edge}" in selected_edges}
        if site_selection == "representative_nearest_w":
            w_sites = [i for i, s in enumerate(structure) if s.specie.symbol == "W"]
            for e, indices in selected_sites.items():
                if not indices:
                    continue
                if e == "W" or not w_sites:
                    chosen = indices[0]
                else:
                    chosen = min(indices, key=lambda i: (
                        min(float(structure.get_distance(i, w)) for w in w_sites), i))
                selected_sites[e] = [chosen]
        selected_indices = {i for indices in selected_sites.values() for i in indices}
        absorber_count = len(selected_indices)
        if not absorber_count:
            fail(f"No Co/Fe/W absorber in {source}")
        site_table = [{"site_index": i, "cif_label": str(s.label),
                       "element": s.specie.symbol, "Z": int(s.specie.Z),
                       "occupancy": 1.0,
                       "frac_coords": [float(v) for v in s.frac_coords]}
                      for i, s in enumerate(structure)]
        for absorber, site in enumerate(structure):
            element = site.specie.symbol
            if absorber not in selected_indices:
                continue
            edge = f"{element}_{EDGES[element]}"
            if edge not in selected_edges:
                continue
            atoms = build_cluster(structure, absorber, feff["cluster_radius_angstrom"])
            if "max_atoms" in feff and len(atoms) > int(feff["max_atoms"]):
                fail(f"{model} site {absorber}: {len(atoms)} cluster atoms exceed "
                     f"feff.max_atoms={feff['max_atoms']}; check the FEFF build limits")
            for mode in feff["modes"]:
                job_path = jobs_dir / args.profile / model / edge / f"site_{absorber:04d}" / mode
                info = {
                    "id": f"{args.profile}/{model}/{edge}/site_{absorber:04d}/{mode}",
                    "model_id": model,
                    "source_cif": source.relative_to(root).as_posix(),
                    "source_cif_sha256": digest(source),
                    "edge": edge, "element": element, "absorber_index": absorber,
                    "site_index": absorber, "mode": mode,
                    "weight": 1.0/len(selected_sites[element]),
                    "path": job_path.relative_to(root).as_posix(),
                    "profile": args.profile,
                    "site_selection": site_selection,
                    "element_site_count": len(selected_sites[element]),
                    "element_total_site_count": counts[element],
                    "is_full_element_average": site_selection == "all",
                    "cluster_atom_count": len(atoms),
                    "interpretation": "geometry-based theoretical XAS; no assigned electrode potential",
                }
                text = make_input(structure, atoms, model, EDGES[element], absorber, mode, feff)
                info["feff_input_sha256"] = hashlib.sha256(text.encode("utf-8")).hexdigest()
                metadata = dict(info)
                metadata.update({
                    "site_index_basis": "zero-based order returned by pymatgen CifParser, primitive=False",
                    "pymatgen_version": importlib.metadata.version("pymatgen"),
                    "lattice_vectors_angstrom": structure.lattice.matrix.tolist(),
                    "composition": dict(counts), "source_site_table": site_table,
                    "cluster_atoms": [dict(a) for a in atoms],
                    "all_hydrogen_retained": True,
                    "cluster_contains_c_axis_periodic_images": any(a["image"][2] for a in atoms),
                    "potential_grouping": "element, with separate central core-hole potential 0",
                    "feff_settings": feff, "parser_warnings": parser_warnings,
                    "convergence_status": "not established; pilot radius/path/corehole convergence required",
                    "missing_from_cif": ["DFT density", "electrode potential", "DFT+U settings", "magnetic moments"],
                })
                existing = job_path / "feff.inp"
                if existing.exists() and existing.read_text(encoding="utf-8") != text:
                    fail(f"Existing input differs: {existing}. Use a new jobs_dir for changed settings; "
                         "stale FEFF outputs must not be mixed with new inputs.")
                prepared.append((job_path, text, metadata))
                jobs.append(info)
        print(f"{rel_model}: {len(structure)} atoms, {absorber_count} absorbers, "
              f"{absorber_count * len(feff['modes'])} independent jobs")
        for message in parser_warnings:
            print(f"  CIF parser warning: {message}", file=sys.stderr)
    # Complete validation before creating any inputs.
    for job_path, text, metadata in prepared:
        write_text(job_path / "feff.inp", text)
        write_text(job_path / "job.json", json.dumps(metadata, ensure_ascii=False, indent=2) + "\n")
    manifest = manifests_dir / f"{args.profile}_jobs.json"
    write_text(manifest, json.dumps(jobs, ensure_ascii=False, indent=2) + "\n")
    print(f"Prepared {len(jobs)} independent jobs. Manifest: {manifest}")
    if site_selection != "all":
        print("Representative-site pilot: these spectra are NOT full element-averaged spectra.")
    print("No FEFF calculation was run. No convergence or electrode-potential assignment is implied.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
