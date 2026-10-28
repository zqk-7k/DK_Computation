#!/usr/bin/env python3
"""Read completed FEFF10 calculations, average signed spectra, and transform with Larch.

No synthetic spectra or experimental fits are produced by the calculation path.
The isolated --self-test uses an OS temporary directory and is labelled TEST ONLY.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import re
import sys
import tempfile
import textwrap
import warnings

import numpy as np
import scipy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import Normalize

DEFAULTS = {
    "output_dir": "results", "energy_step_eV": 0.25,
    "kmin": 3.0, "kmax": 12.0, "kstep": 0.05, "dk": 1.0,
    "nfft": 2048, "rmax": 6.0, "native_k_gap_max": 0.25,
    "endpoint_tolerance": 0.002, "xmu_relation": "diagnostic",
    "exafs_site_average": "native_k_structural_approximation",
    "energy_shift_eV_by_edge": {}, "png_dpi": 450,
}
REFERENCES = {
    "FEFF10": "https://times-software.github.io/feff10/assets/doc/feff10_users_guide.pdf",
    "xftf": "https://xraypy.github.io/xraylarch/xafs_fourier.html",
    "cauchy_wavelet": "https://xraypy.github.io/xraylarch/xafs_wavelets.html",
    "wavelet_method": "Munoz, Argoul and Farges, American Mineralogist 88, 694-700 (2003)",
}
PALETTE = ["#2676A5", "#CC6539", "#438D75", "#8666A7", "#AF8A36", "#555555", "#B85376", "#5C9B9B"]


class DataError(RuntimeError):
    pass


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def dump_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def slug(value):
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("._")
    if not text:
        raise DataError(f"Empty safe name from {value!r}")
    return text


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def larch_api():
    try:
        from larch import Group
        from larch.xafs import xftf, cauchy_wavelet
    except ImportError as exc:
        raise DataError("xraylarch is required; activate the supplied HPC Python environment. No substitute transform is used.") from exc
    return Group, xftf, cauchy_wavelet


def parse_feff_table(path, mode, relation="diagnostic"):
    """Use explicit native FEFF columns; retain original comments and extra columns."""
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    headers, rows = [], []
    for lineno, line in enumerate(text.splitlines(), 1):
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "*", "!")):
            if stripped:
                headers.append(line)
            continue
        tokens = stripped.split()
        try:
            row = [float(v.replace("D", "E").replace("d", "e")) for v in tokens]
        except ValueError:
            if rows:
                raise DataError(f"Non-numeric row after data began: {path}:{lineno}")
            headers.append(line)
            continue
        rows.append(row)
    if len(rows) < 8:
        raise DataError(f"Too few numeric rows in {path}")
    widths = {len(row) for row in rows}
    if len(widths) != 1:
        raise DataError(f"Inconsistent numeric column counts in {path}: {widths}")
    data = np.asarray(rows, dtype=float)
    required = 6 if mode == "xanes" else 2
    if data.shape[1] < required or not np.isfinite(data).all():
        raise DataError(f"Invalid {mode} table in {path}; need >= {required} finite columns")
    if np.any(np.diff(data[:, 0]) <= 0):
        raise DataError(f"First column must increase strictly; duplicate/unsorted points in {path}")
    meta = {"source": str(path), "sha256": sha256(path), "rows": len(data),
            "columns": data.shape[1], "native_header": headers,
            "column_interpretation": (["photon_energy_eV", "relative_energy_eV", "k_A_inv", "mu", "mu0", "native_chi_or_scaled_chi"]
                                      if mode == "xanes" else ["k_A_inv", "signed_chi"]),
            "axis_min": float(data[0, 0]), "axis_max": float(data[-1, 0])}
    if mode == "xanes":
        # FEFF10 defines column 6 as a possibly rescaled chi. Neither identity is
        # assumed universally. The directly calculated mu column is authoritative.
        mu, mu0, chi = data[:, 3], data[:, 4], data[:, 5]
        scale = max(float(np.max(np.abs(mu))), 1e-12)
        additive = float(np.max(np.abs(mu - (mu0 + chi))) / scale)
        multiplicative = float(np.max(np.abs(mu - mu0 * (1 + chi))) / scale)
        meta["mu_identity_diagnostics"] = {"relative_max_error_mu0_plus_chi": additive,
            "relative_max_error_mu0_times_1_plus_chi": multiplicative, "enforced_relation": relation,
            "note": "FEFF8+ scaled chi does not universally obey mu=mu0*(1+chi). mu is read directly, never reconstructed."}
        if relation not in {"diagnostic", "additive", "multiplicative"}:
            raise DataError("xmu_relation must be diagnostic, additive or multiplicative")
        selected = additive if relation == "additive" else multiplicative
        if relation != "diagnostic" and selected > 5e-4:
            raise DataError(f"Explicitly requested xmu {relation} identity fails in {path}: {selected:g}")
        meta["energy_reference"] = "Native photon energy, column 1. No per-site Fermi subtraction or white-line alignment."
    return data, meta


def grid_inside(lo, hi, step):
    start = math.ceil((lo - 1e-9) / step) * step
    stop = math.floor((hi + 1e-9) / step) * step
    n = int(round((stop - start) / step)) + 1
    if n < 2:
        raise DataError("No usable common interpolation interval")
    return start + step * np.arange(n)


def checked_interp(x, y, new_x, tolerance=0.0):
    if new_x[0] < x[0] - tolerance or new_x[-1] > x[-1] + tolerance:
        raise DataError(f"Interpolation would extrapolate [{new_x[0]}, {new_x[-1]}] beyond [{x[0]}, {x[-1]}]")
    # Only roundoff/printed endpoint tolerance may use the nearest boundary value.
    return np.interp(np.clip(new_x, x[0], x[-1]), x, y)


def write_csv(path, headers, array):
    np.savetxt(path, np.asarray(array), delimiter=",", header=",".join(headers), comments="", fmt="%.12g")


def write_matrix(path, k, r, values):
    write_csv(path, ["R_A_phase_uncorrected_vs_k_A_inv"] + [f"{x:.8g}" for x in k],
              np.column_stack([r, values]))


def configure_font():
    names = {f.name for f in font_manager.fontManager.ttflist}
    font = "Arial" if "Arial" in names else "DejaVu Sans"
    note = "Arial located" if font == "Arial" else "Arial unavailable: explicitly using DejaVu Sans; install licensed Arial and rerun for Arial exports."
    if font != "Arial":
        warnings.warn(note, stacklevel=2)
    plt.rcParams.update({"font.family": font, "font.size": 7, "axes.labelsize": 7,
        "axes.titlesize": 8, "legend.fontsize": 6.5, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "axes.linewidth": 0.5, "lines.linewidth": 0.5, "xtick.major.width": 0.5,
        "ytick.major.width": 0.5, "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "legend.frameon": False,
        "mathtext.fontset": "custom", "mathtext.rm": font, "mathtext.it": f"{font}:italic",
        "mathtext.bf": f"{font}:bold"})
    return {"requested_font": "Arial", "actual_font": font, "message": note}


def export_figure(fig, stem, partial, p, footer="Model-derived; not experimental measurements"):
    status = "PARTIAL - DIAGNOSTIC ONLY | " if partial else ""
    fig.suptitle(status + "Calculated FEFF10", fontsize=9, y=0.995)
    if p.get("_representative_sites_only"):
        footer += " | Selected sites, not a full-element average"
    fig.text(0.5, 0.012, footer, ha="center", va="bottom", fontsize=6, color="#555555")
    fig.savefig(str(stem) + ".pdf")
    fig.savefig(str(stem) + ".svg")
    fig.savefig(str(stem) + ".png", dpi=int(p["png_dpi"]))
    fig.savefig(str(stem) + "_preview.png", dpi=160)
    plt.close(fig)


def transform_chi(k, chi, p):
    """Larch only; signed input first, complex arrays retained, no magnitude averaging."""
    Group, xftf, cauchy_wavelet = larch_api()
    if abs(k[0]) > 1e-9 or not np.allclose(np.diff(k), p["kstep"], atol=1e-9):
        raise DataError("Larch transform arrays must be uniform and start at zero")
    if len(k) >= p["nfft"] // 2:
        raise DataError("WT input must be shorter than nfft/2 to avoid Larch truncation")
    out = {"ft": {}}
    for weight in (1, 2, 3):
        group = Group()
        xftf(k, chi, group=group, kmin=p["kmin"], kmax=p["kmax"], kweight=weight,
             dk=p["dk"], dk2=p["dk"], window="hanning", nfft=p["nfft"],
             kstep=p["kstep"], rmax_out=p["rmax"])
        out["ft"][weight] = {"r": np.asarray(group.r), "complex": np.asarray(group.chir),
                              "window": np.asarray(group.kwin)}
    group = Group()
    cauchy_wavelet(k, chi, group=group, kweight=2, rmax_out=p["rmax"], nfft=p["nfft"])
    out["wt"] = {"r": np.asarray(group.wcauchy_r), "k": k.copy(), "complex": np.asarray(group.wcauchy)}
    if out["wt"]["complex"].shape != (len(out["wt"]["r"]), len(k)):
        raise DataError("Unexpected Larch Cauchy array shape; inspect installed Larch version")
    for transform in list(out["ft"].values()) + [out["wt"]]:
        if not np.isfinite(transform["complex"]).all():
            raise DataError("Non-finite Larch transform output")
    return out


def write_transform(directory, k, chi, transforms, p, is_mean=False):
    directory.mkdir(parents=True, exist_ok=True)
    display = (k >= p["kmin"] - 1e-8) & (k <= p["kmax"] + 1e-8)
    all_chi = np.column_stack([k, chi, k * chi, k**2 * chi, k**3 * chi])
    headers = ["k_A_inv", "chi_signed", "k_chi_A_inv", "k2_chi_A_inv2", "k3_chi_A_inv3"]
    write_csv(directory / "chi_full_grid.csv", headers, all_chi)
    write_csv(directory / "chi_common_3_to_12.csv", headers, all_chi[display])
    archive = {"k_A_inv": k, "chi_signed": chi, "display_k_mask": display}
    for weight, tr in transforms["ft"].items():
        value, r = tr["complex"], tr["r"]
        write_csv(directory / f"FT_kweight{weight}.csv", ["R_A_phase_uncorrected", "real", "imag", "magnitude"],
                  np.column_stack([r, value.real, value.imag, np.abs(value)]))
        archive[f"FT_kw{weight}_complex"] = value
        archive[f"FT_kw{weight}_r_A"] = r
        archive[f"FT_kw{weight}_window"] = tr["window"]
    wt = transforms["wt"]
    for part, value in (("REAL", wt["complex"].real), ("IMAG", wt["complex"].imag), ("MAG", np.abs(wt["complex"]))):
        write_matrix(directory / f"WT_Cauchy_{part}.csv", k, wt["r"], value)
    archive.update({"WT_Cauchy_complex": wt["complex"], "WT_k_A_inv": k, "WT_R_A": wt["r"]})
    np.savez_compressed(directory / "complex_transforms.npz", **archive)
    dump_json(directory / "processing.json", {
        "calculation": "FEFF10 model-derived", "aggregate": "site-weighted signed chi before transforms" if is_mean else "single absorber site",
        "exafs_site_average": p["exafs_site_average"],
        "energy_caveat": "EXAFS uses each site's native threshold k. This structural average does not mix site chemical shifts on absolute photon energy.",
        "ft": {"library": "larch.xafs.xftf", "kmin": p["kmin"], "kmax": p["kmax"],
               "dk": p["dk"], "dk2": p["dk"], "window": "hanning", "nfft": p["nfft"],
               "kweights": [1, 2, 3], "phase_corrected": False},
        "wt": {"library": "larch.xafs.cauchy_wavelet", "kweight": 2, "nfft": p["nfft"],
               "input_window": "none; common full available k grid", "input_k_range": [float(k[0]), float(k[-1])],
               "display_k_range": [p["kmin"], p["kmax"]], "color_scale": "shared across all plotted mean maps in this run",
               "edge_caution": "Finite-k boundary effects remain; no Morlet COI mask is applied to the Cauchy transform."},
        "fit": "none; calculated spectra are not fitted to themselves or to unavailable experiment"})


def load_completed_jobs(config_path, profile, allow_partial):
    config_path = Path(config_path).resolve()
    config = read_json(config_path)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", profile):
        raise DataError("Profile must be a safe alphanumeric name with optional underscores/hyphens")
    if profile not in config.get("profiles", {}):
        raise DataError(f"Profile {profile!r} is not defined in config.profiles")
    root = Path(config.get("project_root", config_path.parent))
    if not root.is_absolute():
        root = config_path.parent / root
    root = root.resolve()
    manifest_directory = Path(config.get("paths", {}).get("manifests_dir", "manifests"))
    if not manifest_directory.is_absolute():
        manifest_directory = root / manifest_directory
    manifest = manifest_directory / f"{profile}_jobs.json"
    jobs = read_json(manifest)
    if not isinstance(jobs, list) or not jobs:
        raise DataError(f"Expected a nonempty job list: {manifest}")
    required = {"id", "model_id", "source_cif", "edge", "element", "absorber_index", "mode", "weight", "path"}
    seen, present, problems, totals = set(), [], [], defaultdict(float)
    for listed in jobs:
        missing = required - set(listed)
        if missing:
            raise DataError(f"Manifest job is missing fields: {sorted(missing)}")
        if listed["id"] in seen:
            raise DataError(f"Duplicate job id {listed['id']}")
        seen.add(listed["id"])
        if listed["mode"] not in {"xanes", "exafs"}:
            raise DataError(f"Unknown mode {listed['mode']}")
        weight = float(listed["weight"])
        if not np.isfinite(weight) or weight <= 0:
            raise DataError("All site weights must be positive finite numbers")
        key = (str(listed["model_id"]), str(listed["element"]), str(listed["edge"]), listed["mode"])
        totals[key] += weight
        directory = Path(listed["path"])
        if not directory.is_absolute():
            directory = root / directory
        directory = directory.resolve()
        try:
            local = read_json(directory / "job.json")
            for field in required - {"path"}:
                if local.get(field) != listed[field]:
                    raise DataError(f"job.json mismatch on {field}")
            result = read_json(directory / "run_result.json")
            if result.get("status") != "complete":
                raise DataError(f"run status is {result.get('status')!r}")
            source = directory / ("xmu.dat" if listed["mode"] == "xanes" else "chi.dat")
            if not source.is_file() or source.stat().st_size == 0:
                raise DataError(f"Missing nonempty {source.name}")
            input_path = directory / "feff.inp"
            if result.get("input_sha256") and input_path.is_file() and result["input_sha256"] != sha256(input_path):
                raise DataError("feff.inp changed after successful run; rerun the calculation")
            present.append(dict(listed, directory=directory, source=source, run_result=result,
                                site_selection=local.get("site_selection", listed.get("site_selection", "as specified in manifest")),
                                is_full_element_average=local.get("is_full_element_average", listed.get("is_full_element_average", False))))
        except (OSError, ValueError, DataError) as exc:
            problems.append({"job_id": listed["id"], "path": str(directory), "reason": str(exc)})
    for key, total in totals.items():
        if not math.isclose(total, 1.0, abs_tol=1e-6):
            raise DataError(f"Manifest site weights for {key} sum to {total}, not 1")
    if problems and not allow_partial:
        details = "\n".join(f"  {v['job_id']}: {v['reason']}" for v in problems[:20])
        raise DataError(f"{len(problems)} jobs are incomplete. Refusing a partial site average.\n{details}\nUse --allow-partial only for labelled diagnostic output.")
    if not present:
        raise DataError("No completed FEFF jobs available")
    p = dict(DEFAULTS)
    p.update(config.get("postprocess", {}))
    # Fixed scientific contract; config cannot silently alter these transforms.
    for key, expected in {"kmin": 3.0, "kmax": 12.0, "kstep": 0.05, "dk": 1.0, "nfft": 2048}.items():
        if not np.isclose(float(p[key]), expected):
            raise DataError(f"This workflow contract requires {key}={expected}; revise the documented method before changing it")
    p["nfft"] = int(p["nfft"])
    if p["exafs_site_average"] != "native_k_structural_approximation":
        raise DataError("Absolute-energy EXAFS mixing is not implemented: provide and validate per-site E(k) mappings in a separate workflow")
    if p["energy_step_eV"] <= 0 or p["rmax"] <= 0 or p["endpoint_tolerance"] > 0.005:
        raise DataError("Invalid postprocessing step/range/tolerance")
    return root, config, p, jobs, present, problems, manifest


def process_xanes(edge_jobs, output, p, audit):
    parsed = []
    for job in edge_jobs:
        data, meta = parse_feff_table(job["source"], "xanes", p["xmu_relation"])
        parsed.append((job, data, meta))
    edge_id = edge_key(edge_jobs[0])
    shift = float(p["energy_shift_eV_by_edge"].get(edge_id, 0.0))
    lo, hi = max(d[0, 0] for _, d, _ in parsed), min(d[-1, 0] for _, d, _ in parsed)
    grid = grid_inside(lo, hi, float(p["energy_step_eV"]))
    if len(grid) < 20 or grid[-1] - grid[0] < 10:
        raise DataError(f"Insufficient shared absolute-energy coverage for {edge_id}")
    by_model = defaultdict(list)
    for job, data, meta in parsed:
        columns = np.column_stack([checked_interp(data[:, 0], data[:, i], grid) for i in (3, 4, 5)])
        directory = output / "sites" / slug(job["model_id"]) / edge_id / f"site_{job['absorber_index']}" / "xanes"
        directory.mkdir(parents=True, exist_ok=True)
        write_csv(directory / "native_xmu.csv", ["photon_energy_eV", "relative_energy_eV", "k_A_inv", "mu", "mu0", "native_chi_or_scaled_chi"] + [f"extra_native_col_{i+1}" for i in range(6, data.shape[1])], data)
        write_csv(directory / "common_energy_mu.csv", ["photon_energy_eV", "mu", "mu0", "native_chi_or_scaled_chi"], np.column_stack([grid + shift, columns]))
        meta.update({"job_id": job["id"], "common_energy_step_eV": p["energy_step_eV"],
                     "shared_edge_shift_eV": shift, "interpolation": "linear; no extrapolation",
                     "common_grid_range_eV": [float(grid[0] + shift), float(grid[-1] + shift)]})
        dump_json(directory / "source_and_processing.json", meta)
        by_model[str(job["model_id"])].append((job, columns))
        audit.append(meta)
    results = []
    for model, entries in sorted(by_model.items()):
        original = np.asarray([float(j["weight"]) for j, _ in entries])
        weights = original / original.sum()
        mean = np.tensordot(weights, np.asarray([cols for _, cols in entries]), axes=(0, 0))
        directory = output / "means" / slug(model) / edge_id / "xanes"
        directory.mkdir(parents=True, exist_ok=True)
        write_csv(directory / "mean_mu.csv", ["photon_energy_eV", "mean_mu", "mean_mu0", "mean_native_chi_or_scaled_chi"], np.column_stack([grid + shift, mean]))
        dump_json(directory / "averaging.json", {"job_ids": [j["id"] for j, _ in entries], "original_weight_sum": float(original.sum()),
                    "normalized_weights_used": weights.tolist(), "operation": "Average mu on shared absolute photon energy before any transform; no peak alignment or separate normalization.",
                    "is_full_element_average": all(j.get("is_full_element_average", False) for j, _ in entries) and math.isclose(float(original.sum()), 1.0, abs_tol=1e-6),
                    "number_of_available_sites": len(entries), "original_weights": original.tolist(),
                    "site_selection": sorted(set(j.get("site_selection", "unspecified") for j, _ in entries)),
                    "shared_edge_shift_eV": shift})
        results.append({"model": model, "edge": edge_id, "energy": grid + shift, "mu": mean[:, 0]})
    return results


def edge_key(job):
    edge = str(job["edge"]).replace(" ", "_")
    element = str(job["element"])
    return slug(edge if edge.lower().startswith(element.lower() + "_") else element + "_" + edge)


def process_exafs(edge_jobs, output, p, audit):
    parsed = []
    tolerance = float(p["endpoint_tolerance"])
    support_min, support_max = p["kmin"] - p["dk"] / 2, p["kmax"] + p["dk"] / 2
    for job in edge_jobs:
        data, meta = parse_feff_table(job["source"], "exafs")
        native_data = data.copy()
        if np.any(data[:, 0] < -tolerance):
            raise DataError(f"Negative k in chi.dat: {job['source']}")
        if np.max(np.diff(data[:, 0])) > p["native_k_gap_max"] + tolerance:
            raise DataError(f"Large native k gap in {job['source']}; refusing interpolation")
        if data[0, 0] > p["kstep"] + tolerance:
            raise DataError(f"Native chi must begin at k=0 or near one grid step for full-grid Cauchy WT: {job['source']}")
        if data[-1, 0] < support_max - tolerance or data[0, 0] > support_min + tolerance:
            raise DataError(f"Hanning kmin=3,kmax=12,dk=1 needs native k support 2.5-12.5; insufficient in {job['source']}")
        meta["zero_origin_added"] = bool(data[0, 0] > 1e-9)
        if data[0, 0] > 1e-9:
            data = np.vstack([np.zeros(data.shape[1]), data])
        else:
            data[0, 0] = 0.0
        parsed.append((job, data, meta, native_data))
    end = math.floor((min(d[-1, 0] for _, d, _, _ in parsed) + tolerance) / p["kstep"]) * p["kstep"]
    k = p["kstep"] * np.arange(int(round(end / p["kstep"])) + 1)
    by_model = defaultdict(list)
    edge_id = edge_key(edge_jobs[0])
    for job, data, meta, native_data in parsed:
        chi = checked_interp(data[:, 0], data[:, 1], k, tolerance)
        directory = output / "sites" / slug(job["model_id"]) / edge_id / f"site_{job['absorber_index']}" / "exafs"
        transforms = transform_chi(k, chi, p)
        write_transform(directory, k, chi, transforms, p)
        write_csv(directory / "native_chi.csv", ["k_A_inv", "chi_signed"] + [f"extra_native_col_{i+1}" for i in range(2, native_data.shape[1])], native_data)
        meta.update({"job_id": job["id"], "common_full_k_range": [0.0, float(k[-1])],
                     "kstep": p["kstep"], "interpolation": "linear; endpoint clipping only within printed-coordinate tolerance",
                     "endpoint_tolerance_A_inv": tolerance,
                     "zero_origin_note": "If native k starts near 0.05, (0,0) is prepended; k^w chi at zero is zero for w=1,2,3.",
                     "phase_information": "Signed native chi preserved. FT/WT real and imaginary arrays retained."})
        dump_json(directory / "source_and_processing.json", meta)
        by_model[str(job["model_id"])].append((job, chi, transforms))
        audit.append(meta)
    results = []
    for model, entries in sorted(by_model.items()):
        original = np.asarray([float(j["weight"]) for j, _, _ in entries])
        weights = original / original.sum()
        mean_chi = np.tensordot(weights, np.asarray([chi for _, chi, _ in entries]), axes=(0, 0))
        transforms = transform_chi(k, mean_chi, p)
        # Linearity QA detects accidental averaging of |FT|/|WT| instead of chi.
        residuals = {}
        for weight in (1, 2, 3):
            average_complex = np.tensordot(weights, np.asarray([tr["ft"][weight]["complex"] for _, _, tr in entries]), axes=(0, 0))
            delta = float(np.max(np.abs(average_complex - transforms["ft"][weight]["complex"])))
            if not np.allclose(average_complex, transforms["ft"][weight]["complex"], rtol=1e-9, atol=1e-10):
                raise DataError("FT linearity check failed")
            residuals[f"FT_kw{weight}"] = delta
        average_complex = np.tensordot(weights, np.asarray([tr["wt"]["complex"] for _, _, tr in entries]), axes=(0, 0))
        if not np.allclose(average_complex, transforms["wt"]["complex"], rtol=1e-9, atol=1e-10):
            raise DataError("WT linearity check failed")
        residuals["WT"] = float(np.max(np.abs(average_complex - transforms["wt"]["complex"])))
        directory = output / "means" / slug(model) / edge_id / "exafs"
        write_transform(directory, k, mean_chi, transforms, p, is_mean=True)
        dump_json(directory / "averaging.json", {"job_ids": [j["id"] for j, _, _ in entries], "original_weight_sum": float(original.sum()),
            "normalized_weights_used": weights.tolist(), "operation": "Average signed chi on the common native-threshold k grid, then apply xftf/Cauchy WT.",
            "is_full_element_average": all(j.get("is_full_element_average", False) for j, _, _ in entries) and math.isclose(float(original.sum()), 1.0, abs_tol=1e-6),
            "number_of_available_sites": len(entries), "original_weights": original.tolist(),
            "site_selection": sorted(set(j.get("site_selection", "unspecified") for j, _, _ in entries)),
            "approximation": p["exafs_site_average"], "complex_linearity_max_errors": residuals,
            "not_done": "No magnitude averaging; no phase correction; no experimental fitting; no invented energy-to-k mappings."})
        results.append({"model": model, "edge": edge_id, "k": k, "chi": mean_chi, "transforms": transforms})
    return results


def make_plots(xanes, exafs, output, p, partial):
    figures = output / "figures"
    figures.mkdir()
    by_edge = defaultdict(list)
    for item in xanes:
        by_edge[item["edge"]].append(item)
    for edge, entries in sorted(by_edge.items()):
        for page, start in enumerate(range(0, len(entries), 6), 1):
            fig, ax = plt.subplots(figsize=(183 / 25.4, 95 / 25.4))
            for i, item in enumerate(entries[start:start+6]):
                ax.plot(item["energy"], item["mu"], color=PALETTE[i % len(PALETTE)], label=item["model"])
            ax.set(xlabel="Photon energy (eV)", ylabel=r"Calculated $\mu$ (FEFF units)", title=edge.replace("_", " "))
            ax.ticklabel_format(axis="x", useOffset=False, style="plain")
            handles, labels = ax.get_legend_handles_labels()
            fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.95), ncol=2, fontsize=6.5)
            fig.subplots_adjust(left=0.10, right=0.98, bottom=0.17, top=0.72)
            export_figure(fig, figures / f"{edge}_XANES_page{page}", partial, p)
    by_edge = defaultdict(list)
    for item in exafs:
        by_edge[item["edge"]].append(item)
    # One global raw WT scale for every mean map, without per-map normalization.
    visible_values = []
    for item in exafs:
        wt = item["transforms"]["wt"]
        keep = (wt["k"] >= p["kmin"] - 1e-8) & (wt["k"] <= p["kmax"] + 1e-8)
        visible_values.append(float(np.max(np.abs(wt["complex"][:, keep]))))
    vmax = max(visible_values, default=1.0) or 1.0
    pages = [(edge, entries[start:start+6], start//6+1) for edge, entries in sorted(by_edge.items()) for start in range(0, len(entries), 6)]
    for edge, entries, page in pages:
        fig, axes = plt.subplots(3, 2, figsize=(183 / 25.4, 165 / 25.4))
        for i, item in enumerate(entries):
            color = PALETTE[i % len(PALETTE)]
            keep = (item["k"] >= p["kmin"] - 1e-8) & (item["k"] <= p["kmax"] + 1e-8)
            for row, weight in enumerate((1, 2, 3)):
                tr = item["transforms"]["ft"][weight]
                axes[row, 0].plot(item["k"][keep], (item["k"]**weight * item["chi"])[keep], color=color, label=item["model"])
                axes[row, 1].plot(tr["r"], np.abs(tr["complex"]), color=color)
        for row, weight in enumerate((1, 2, 3)):
            axes[row, 0].set(xlabel=r"$k$ ($\mathrm{\AA}^{-1}$)", ylabel=rf"$k^{weight}\chi(k)$ ($\mathrm{{\AA}}^{{-{weight}}}$)")
            axes[row, 1].set(xlabel=r"$R$ ($\mathrm{\AA}$), phase uncorrected", ylabel=rf"$|\chi(R)|$ ($\mathrm{{\AA}}^{{-{weight+1}}}$)")
            for col in range(2):
                axes[row, col].text(-0.13, 1.03, "abcdef"[2 * row + col], transform=axes[row, col].transAxes, fontweight="bold", fontsize=9)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.963), ncol=2)
        fig.subplots_adjust(left=0.10, right=0.98, bottom=0.095, top=0.84, wspace=0.32, hspace=0.40)
        export_figure(fig, figures / f"{edge}_EXAFS_FT_page{page}", partial, p, f"{edge.replace('_', ' ')} | Model-derived; native-k structural average; no experimental fit")
        for start in (0,):
            chunk = entries
            ncols = min(3, len(chunk)); nrows = math.ceil(len(chunk) / ncols)
            fig, axes = plt.subplots(nrows, ncols, figsize=(183 / 25.4, (25 + 64 * nrows) / 25.4), squeeze=False)
            image = None
            for i, item in enumerate(chunk):
                ax = axes.flat[i]; wt = item["transforms"]["wt"]
                keep = (wt["k"] >= p["kmin"] - 1e-8) & (wt["k"] <= p["kmax"] + 1e-8)
                image = ax.pcolormesh(wt["k"][keep], wt["r"], np.abs(wt["complex"][:, keep]), shading="auto", cmap="viridis", norm=Normalize(0, vmax), rasterized=True)
                title_lines = [line for segment in item["model"].split("__")
                               for line in textwrap.wrap(segment, width=27, break_long_words=True, break_on_hyphens=False)]
                ax.set(xlabel=r"$k$ ($\mathrm{\AA}^{-1}$)", ylabel=r"$R$ ($\mathrm{\AA}$)", title="\n".join(title_lines))
            for ax in list(axes.flat)[len(chunk):]:
                ax.set_visible(False)
            fig.subplots_adjust(left=0.08, right=0.86, bottom=0.17 if nrows == 1 else 0.10, top=0.80 if nrows == 1 else 0.88, wspace=0.40, hspace=0.45)
            cax = fig.add_axes([0.89, 0.25, 0.016, 0.48])
            fig.colorbar(image, cax=cax, label=r"Cauchy WT magnitude (a.u.)")
            export_figure(fig, figures / f"{edge}_WT_page{page}", partial, p, f"{edge.replace('_', ' ')} | Model-derived; common raw color scale; finite-k boundary effects")
    return {"WT_common_vmin": 0.0, "WT_common_vmax": vmax, "WT_normalization": "none",
            "WT_scale_scope": "All plotted mean spectra, all available edges, visible k=3-12 range"}


def run(config_path, profile, allow_partial=False):
    root, config, p, jobs, completed, problems, manifest = load_completed_jobs(config_path, profile, allow_partial)
    larch_api()  # Fail before making results if the actual transformation dependency is missing.
    partial = bool(problems)
    p["_representative_sites_only"] = not all(j.get("is_full_element_average", False) for j in completed)
    base = Path(p["output_dir"])
    if not base.is_absolute():
        base = root / base
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output = base / profile / (("PARTIAL_DIAGNOSTIC_" if partial else "COMPLETE_") + stamp)
    output.mkdir(parents=True, exist_ok=False)
    report = {"status": "processing", "data_type": "Calculated FEFF10; model-derived; not experimental measurements",
        "profile": profile, "partial": partial, "expected_jobs": len(jobs), "completed_jobs": len(completed),
        "unavailable_jobs": problems, "manifest": str(manifest), "manifest_sha256": sha256(manifest),
        "site_selection": sorted(set(j.get("site_selection", "unspecified") for j in completed)),
        "is_full_element_average": all(j.get("is_full_element_average", False) for j in completed) and not partial,
        "config_sha256": sha256(config_path), "parameters": p, "references": REFERENCES,
        "figure_contract": {"claim": "Compare model-derived spectra from supplied CIF structures without assigning electrochemical states or proving exclusive active sites.",
            "evidence": "FEFF native outputs -> common-grid signed site averages -> Larch complex FT/WT -> source-linked figures",
            "archetype": "quantitative grid", "backend": "Python/matplotlib", "exports": "editable-text PDF/SVG and PNG; all source arrays and complex transforms retained",
            "risks": "Structural labels are not OCP/potentials. Native-k EXAFS averaging is approximate; absolute-energy chemical shifts are retained only for XANES."},
        "versions": {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
                     "xraylarch": importlib.metadata.version("xraylarch")}, "font": configure_font()}
    audit, xanes, exafs = [], [], []
    try:
        by_type = defaultdict(list)
        for job in completed:
            by_type[(edge_key(job), job["mode"])].append(job)
        for (_, mode), selected in sorted(by_type.items()):
            if mode == "xanes":
                xanes.extend(process_xanes(selected, output, p, audit))
            else:
                exafs.extend(process_exafs(selected, output, p, audit))
        report["plot_settings"] = make_plots(xanes, exafs, output, p, partial)
        report["source_audit"] = audit
        report["status"] = "PARTIAL_DIAGNOSTIC" if partial else "complete"
        report["models_and_edges"] = {"xanes": [{"model": v["model"], "edge": v["edge"]} for v in xanes],
                                      "exafs": [{"model": v["model"], "edge": v["edge"]} for v in exafs]}
        report["interpretation_limits"] = ["Calculated spectra do not establish an exclusive catalytic active site.",
            "No calculated model is labelled OCP or assigned electrode potential.",
            "No experimental fit, coordination-number fit or oxidation-state calibration has been performed.",
            "FT R is phase uncorrected; peak locations are not direct bond lengths.",
            "EXAFS native-k averaging preserves signed oscillations but omits absolute-energy mixing of site chemical shifts."]
        dump_json(output / "postprocess_report.json", report)
        files = [{"path": str(f.relative_to(output)), "bytes": f.stat().st_size, "sha256": sha256(f)}
                 for f in sorted(output.rglob("*")) if f.is_file()]
        dump_json(output / "file_inventory.json", files)
        dump_json(base / profile / "latest_run.json", {"output": str(output), "status": report["status"], "profile": profile})
    except Exception as exc:
        report["status"] = "failed"
        report["error"] = str(exc)
        dump_json(output / "postprocess_report.json", report)
        raise
    print(json.dumps({"status": report["status"], "output": str(output), "completed_jobs": len(completed)}, ensure_ascii=False))
    return output


def self_test():
    """Isolated numerical fixture, never included in FEFF calculation results."""
    p = dict(DEFAULTS)
    larch_api()
    with tempfile.TemporaryDirectory(prefix="FEFF_POSTPROCESS_TEST_ONLY_") as dirname:
        temp = Path(dirname)
        energy = np.arange(10000, 10020.01, 0.25)
        mu0 = 0.8 + 0.002 * (energy - energy[0])
        chi = 0.1 * np.sin((energy - energy[0]) / 3)
        table = np.column_stack([energy, energy - 10005, np.sqrt(np.maximum(energy - 10005, 0) / 3.81), mu0 + chi, mu0, chi])
        np.savetxt(temp / "xmu.dat", table, header="TEST ONLY; omega e k mu mu0 scaled_chi", fmt="%.12g")
        loaded, meta = parse_feff_table(temp / "xmu.dat", "xanes", "additive")
        assert np.allclose(loaded, table)
        rejected = False
        try:
            parse_feff_table(temp / "xmu.dat", "xanes", "multiplicative")
        except DataError:
            rejected = True
        assert rejected, "Inapplicable mu relation must not silently pass"
        rejected = False
        try:
            checked_interp(energy, table[:, 3], np.array([9999.0, 10010.0]))
        except DataError:
            rejected = True
        assert rejected, "Extrapolation guard must work"
        k = np.arange(281) * 0.05
        first = 0.1 * np.sin(4 * k) * np.exp(-k / 9)
        second = 0.075 * np.sin(5.7 * k + 0.2) * np.exp(-k / 8)
        np.savetxt(temp / "chi.dat", np.column_stack([k, first]), header="TEST ONLY; k chi")
        loaded, _ = parse_feff_table(temp / "chi.dat", "exafs")
        assert np.allclose(loaded[:, 1], first)
        t1, t2 = transform_chi(k, first, p), transform_chi(k, second, p)
        mixed = transform_chi(k, 0.4 * first + 0.6 * second, p)
        for weight in (1, 2, 3):
            assert np.allclose(mixed["ft"][weight]["complex"], 0.4 * t1["ft"][weight]["complex"] + 0.6 * t2["ft"][weight]["complex"], atol=1e-10)
        assert np.allclose(mixed["wt"]["complex"], 0.4 * t1["wt"]["complex"] + 0.6 * t2["wt"]["complex"], atol=1e-10)
        assert not np.allclose(np.abs(mixed["ft"][2]["complex"]), 0.4 * np.abs(t1["ft"][2]["complex"]) + 0.6 * np.abs(t2["ft"][2]["complex"]))
    print(json.dumps({"self_test": "PASS", "data": "TEST ONLY analytical fixture, not FEFF results", "fixture_location": "Temporary directory removed", "checks": ["native columns", "mu identity applicability", "no extrapolation", "signed average", "complex FT/WT linearity", "no magnitude averaging"]}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--profile", help="Safe profile name defined in config.profiles (pilot/core/full or a convergence profile)")
    parser.add_argument("--allow-partial", action="store_true", help="Explicitly labelled diagnostics only; missing site weights are renormalized and reported")
    parser.add_argument("--self-test", action="store_true", help="Run an isolated TEST ONLY numerical fixture; does not read or create calculation results")
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
        elif args.config and args.profile:
            run(args.config, args.profile, args.allow_partial)
        else:
            parser.error("Provide --config and --profile, or --self-test")
    except (DataError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
