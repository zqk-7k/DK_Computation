#!/usr/bin/env python3
"""Prepare four finite-POM FEFF jobs; never replicate the storage CIF cell."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
Z = {"O": 8, "P": 15, "W": 74}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False)+"\n", encoding="utf-8", newline="\n")


def read_config(path=None):
    if path is None:
        path = ROOT/"config_reference.json"
    path = Path(path).resolve()
    if path.is_file():
        return path, json.loads(path.read_text(encoding="utf-8-sig"))
    if path != ROOT/"config_reference.json":
        raise FileNotFoundError(path)
    return path, {}


def read_xyz(path):
    lines = Path(path).read_text(encoding="utf-8-sig").splitlines()
    count = int(lines[0])
    atoms = []
    for line in lines[2:]:
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 4 or fields[0] not in Z:
            raise ValueError("Expected one P/W/O Cartesian XYZ atom per line")
        xyz = [float(v) for v in fields[1:]]
        if not all(math.isfinite(v) for v in xyz):
            raise ValueError("Nonfinite coordinate")
        atoms.append({"element": fields[0], "xyz_A": xyz, "source_index": len(atoms)})
    if len(atoms) != count or count != 53 or Counter(a["element"] for a in atoms) != {"P": 1, "W": 12, "O": 40}:
        raise ValueError("Full molecular reference must contain exactly P1 W12 O40 (53 atoms)")
    if min(math.dist(a["xyz_A"], b["xyz_A"]) for i, a in enumerate(atoms) for b in atoms[i+1:]) < 1.0:
        raise ValueError("Unexpected severe atomic overlap")
    return atoms


def assign_potentials(atoms, absorber):
    tungsten = [a for a in atoms if a["element"] == "W"]
    oxygen_degree = Counter()
    for atom in atoms:
        if atom["element"] == "O":
            degree = sum(math.dist(atom["xyz_A"], w["xyz_A"]) <= 2.6 for w in tungsten)
            if degree not in {1, 2, 3}:
                raise ValueError("Unrecognized oxygen W coordination in POM")
            oxygen_degree[degree] += 1
            atom["oxygen_W_degree"] = degree
            atom["ipot"] = {1: 3, 2: 4, 3: 5}[degree]
        elif atom["element"] == "P":
            atom["ipot"] = 2
        else:
            atom["ipot"] = 0 if atom["source_index"] == absorber else 1
    if oxygen_degree != {1: 12, 2: 24, 3: 4}:
        raise ValueError("POM oxygen-class counts differ from the crystallographic model")
    return oxygen_degree


def finite_grid(raw):
    if not isinstance(raw, (list, tuple)) or len(raw) != 3:
        raise ValueError("xanes_relative_egrid_eV must have three values")
    low, high, step = map(float, raw)
    if not all(math.isfinite(v) for v in [low, high, step]) or not low < 0 < high or step <= 0:
        raise ValueError("Invalid XANES relative EGRID")
    count = (high-low)/step
    if abs(count-round(count)) > 1e-7 or count+1 > 1800:
        raise ValueError("Use an integral EGRID span with no more than 1800 points")
    if high > 137.1594:
        raise ValueError("Upper XANES EGRID exceeds approximate k=6 A^-1 FMS range")
    return [low, high, step]


def input_text(atoms, absorber, profile, mode, radius, settings):
    centre = atoms[absorber]["xyz_A"]
    grid = finite_grid(settings.get("xanes_relative_egrid_eV", [-160, 130, .5]))
    lines = [
        "TITLE Calculated finite POM reference: crystallographic PW12O40 geometry",
        f"TITLE representative W source index {absorber}; profile {profile}; mode {mode}",
        "* All 53 atoms retained exactly; no periodic replication or surface assignment.",
        "* Formal anion charge is -3; counterions, solvent and excess charge are NOT explicitly modeled.",
        "* Neutral-atom FEFF electronic baseline; no inferred W(VI) occupation or applied voltage.",
        "EDGE L3", "S02 1.0", "COREHOLE FSR", "CONTROL 1 1 1 1 1 1",
        "PRINT 1 0 0 1 0 3", f"SCF {radius:.8f} 1 50 0.2 1",
    ]
    if mode == "xanes":
        lines += ["XANES 6 0.05 0.2", f"FMS {radius:.8f} 1", "RPATH -1",
                  "EGRID", "e_grid " + " ".join(f"{v:g}" for v in grid)]
    else:
        lines += ["EXAFS 14", f"RPATH {float(settings.get('rpath_A', 5.0)):g}",
                  f"NLEG {int(settings.get('nleg', 4))}"]
    counts = Counter(a["ipot"] for a in atoms)
    labels = {0: (74, "Wabs"), 1: (74, "Wrest"), 2: (15, "P"),
              3: (8, "Oterm"), 4: (8, "Obridg"), 5: (8, "Ocentr")}
    lines += ["", "POTENTIALS", "* ipot Z tag lmax1 lmax2 xnatph"]
    for ipot, (z, label) in labels.items():
        lines.append(f"{ipot} {z} {label} -1 -1 {counts[ipot]}")
    lines += ["", "ATOMS", "* x y z ipot tag distance"]
    ordered = sorted(atoms, key=lambda a: (a["source_index"] != absorber,
                                          math.dist(a["xyz_A"], centre), a["source_index"]))
    for atom in ordered:
        xyz = [x-y for x, y in zip(atom["xyz_A"], centre)]
        distance = math.sqrt(sum(x*x for x in xyz))
        lines.append(" ".join(f"{x:.9f}" for x in xyz)+f" {atom['ipot']} {atom['element']}{atom['source_index']:03d} {distance:.9f}")
    lines += ["END", ""]
    return "\n".join(lines)


def prepare(config_path=None):
    config_path, config = read_config(config_path)
    root = config_path.parent
    source = root/"references/alpha_PW12O40_crystallographic.xyz"
    meta_path = root/"references/reference_model_manifest.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    full = meta["full_model"]
    if not meta.get("nonperiodic") or not meta.get("pbc_forbidden_for_feff"):
        raise ValueError("Reference metadata must explicitly require a finite model")
    if not full.get("all_W_in_single_crystallographic_orbit"):
        raise ValueError("One-site reference calculation needs established W equivalence")
    absorber = int(full["representative_W_index_zero_based"])
    atoms = read_xyz(source)
    if atoms[absorber]["element"] != "W":
        raise ValueError("Representative absorber is not W")
    def radial_fingerprint(index):
        return {element: sorted(math.dist(a["xyz_A"], atoms[index]["xyz_A"])
                                for a in atoms if a["element"] == element)
                for element in ["P", "W", "O"]}
    reference_radial = radial_fingerprint(absorber)
    deviation = 0.0
    for i, atom in enumerate(atoms):
        if atom["element"] != "W":
            continue
        actual_radial = radial_fingerprint(i)
        deviation = max(deviation, max(abs(a-b) for element in reference_radial
                                      for a, b in zip(reference_radial[element], actual_radial[element])))
    if deviation > 2e-5:
        raise ValueError("W radial environments differ; cannot use one representative W for this reference")
    assign_potentials(atoms, absorber)
    centre = atoms[absorber]["xyz_A"]
    maximum = max(math.dist(a["xyz_A"], centre) for a in atoms)
    full_radius = math.ceil((maximum+.1)*100)/100
    settings = config.get("feff", {})
    small_radius = float(settings.get("small_radius_A", 6.0))
    if not 3 < small_radius < full_radius:
        raise ValueError("Small radius must include the first shell and be smaller than full-molecule radius")
    if int(settings.get("nleg", 4)) < 2 or float(settings.get("rpath_A", 5.0)) <= 0:
        raise ValueError("Invalid EXAFS path settings")
    grid = finite_grid(settings.get("xanes_relative_egrid_eV", [-160, 130, .5]))
    window = list(map(float, settings.get("absolute_window_eV", [10100, 10300])))
    if len(window) != 2 or not all(math.isfinite(x) for x in window) or window[0] >= window[1]:
        raise ValueError("Invalid absolute window")
    records = []
    for profile, radius in [("small", small_radius), ("full", full_radius)]:
        for mode in ["xanes", "exafs"]:
            directory = root/"jobs_reference"/profile/mode
            text = input_text(atoms, absorber, profile, mode, radius, settings)
            job = {"id": f"{profile}/{mode}", "profile": profile, "mode": mode,
                   "edge": "W_L3", "path": directory.relative_to(root).as_posix(),
                   "source_xyz": source.relative_to(root).as_posix(), "source_sha256": sha(source),
                   "source_metadata_sha256": sha(meta_path), "absorber_source_index": absorber,
                   "equivalent_W_sites": 12, "atom_count": len(atoms), "finite_molecule": True,
                   "radial_equivalence_max_difference_A": deviation,
                   "periodic_images_added": 0, "scf_radius_A": radius,
                   "fms_radius_A": radius if mode == "xanes" else None,
                   "full_geometry_maximum_absorber_distance_A": maximum,
                   "scf_cluster_atoms": sum(math.dist(a["xyz_A"], centre) <= radius for a in atoms),
                   "relative_egrid_eV": grid if mode == "xanes" else None,
                   "absolute_window_eV": window if mode == "xanes" else None,
                   "feff_input_sha256": hashlib.sha256(text.encode()).hexdigest(),
                   "scope": "Unrelaxed crystallographic molecular reference, not the supported catalyst; no imposed ionic charge, potential or active-site assignment"}
            directory.mkdir(parents=True, exist_ok=True)
            inp = directory/"feff.inp"
            if inp.exists() and inp.read_bytes() != text.encode():
                raise ValueError(f"Existing input differs: {inp}; preserve old calculations and use a new project for changed settings")
            inp.write_text(text, encoding="utf-8", newline="\n")
            write_json(directory/"job.json", job)
            records.append(job)
    write_json(root/"jobs_reference/manifest.json", records)
    print(f"Prepared four finite-reference jobs: small={small_radius:g} A, full={full_radius:g} A; all 53 atoms retained in every input")
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path)
    args = parser.parse_args()
    prepare(args.config)


if __name__ == "__main__":
    main()
