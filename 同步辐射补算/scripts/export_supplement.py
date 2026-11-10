#!/usr/bin/env python3
"""Source-linked TXT tables and figures from completed, calculated FEFF spectra.

This module never creates spectra, adjusts individual peak positions, changes FT
distances, or infers electrochemical potentials from a structural model.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys
import textwrap

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

from postprocess import DataError, configure_font, slug

EDGES = ("W_L3", "Co_K", "Fe_K")
EDGE_TITLES = {"W_L3": r"W $L_3$ edge", "Co_K": "Co K edge", "Fe_K": "Fe K edge"}
COLORS = ("#555555", "#3175A3", "#C56A3D")
STYLES = ("-", "--", "-.")
FORMATS = ("png", "pdf", "svg", "tiff")


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2,
                                   allow_nan=False) + "\n", encoding="utf-8")


def load_csv(path, min_cols=2):
    path = Path(path)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream)
        columns = next(reader)
        rows = [[float(x) for x in row] for row in reader if row]
    array = np.asarray(rows, dtype=float)
    if array.ndim != 2 or array.shape[1] < min_cols or not np.isfinite(array).all():
        raise DataError(f"Invalid finite numeric table: {path.name}")
    if len(columns) != array.shape[1]:
        raise DataError(f"Header/data column mismatch: {path.name}")
    return columns, array


def save_txt(path, columns, array, notes):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    header = "Calculated FEFF10 data; structure/adsorbate models, not measured spectra.\n"
    header += "\n".join(str(v) for v in notes) + "\n" + "\t".join(columns)
    np.savetxt(path, np.asarray(array), delimiter="\t", header=header,
               comments="# ", fmt="%.12g", encoding="utf-8")


def validate_source(result_dir):
    report = read_json(result_dir / "postprocess_report.json")
    if report.get("status") != "complete" or report.get("partial"):
        raise DataError("A complete postprocess_report is required; partial diagnostics are not final data.")
    if not report.get("is_full_element_average", False):
        raise DataError("Supplement figures require full-element site averages, not representative pilot spectra.")
    inventory = read_json(result_dir / "file_inventory.json")
    lookup = {}
    for entry in inventory:
        relative = Path(entry["path"])
        path = (result_dir / relative).resolve()
        if relative.is_absolute() or not path.is_relative_to(result_dir.resolve()):
            raise DataError("Source inventory path escapes the result directory")
        if relative.suffix in {".csv", ".npz", ".json"}:
            if not path.is_file() or sha(path) != entry["sha256"]:
                raise DataError(f"Source result changed or missing: {relative.as_posix()}")
            lookup[relative.as_posix()] = entry
    for tree in ("means", "sites"):
        for path in (result_dir / tree).rglob("*"):
            if path.is_file() and path.suffix in {".csv", ".npz"}:
                if path.relative_to(result_dir).as_posix() not in lookup:
                    raise DataError("Scientific table is absent from the verified source inventory")
    return report, lookup


def validate_comparisons(config, source_report, config_dir):
    export = config.get("export", {})
    comparisons = export.get("comparisons", [])
    if not comparisons:
        raise DataError("Define export.comparisons with exactly three model/label series per comparison.")
    available = {}
    for mode in ("xanes", "exafs"):
        available[mode] = {(x["model"], x["edge"]) for x in source_report["models_and_edges"][mode]}
    used = set()
    for comparison in comparisons:
        key = slug(comparison["id"])
        if key in used:
            raise DataError("Duplicate comparison ID")
        used.add(key)
        series = comparison.get("series", [])
        if len(series) != 3 or len({s["model"] for s in series}) != 3:
            raise DataError("Each comparison needs three distinct models; missing states must not be duplicated.")
        for item in series:
            label = str(item["label"])
            # Use structural labels until independently documented potential mapping exists.
            if re.search(r"\bocp\b|\b\d+(?:\.\d+)?\s*v\b|versus\s+rhe|vs\.?\s*rhe", label, re.I):
                raise DataError("Use structural state labels. OCP/voltage assignment is not established by this workflow.")
            for mode in available:
                for edge in EDGES:
                    if (item["model"], edge) not in available[mode]:
                        raise DataError(f"Missing {mode}: model={item['model']}, edge={edge}")
    shifts = export.get("energy_reference_shift_eV", {})
    prior_shifts = source_report.get("parameters", {}).get("energy_shift_eV_by_edge", {})
    if any(abs(float(v)) > 1e-12 for v in prior_shifts.values()):
        raise DataError("Export expects native-energy postprocessing; apply at most one documented global edge calibration here.")
    if any(abs(float(v)) > 1e-12 for v in shifts.values()):
        evidence = export.get("energy_reference_evidence")
        if not isinstance(evidence, dict) or not evidence.get("description") or not evidence.get("file"):
            raise DataError("Nonzero common edge shifts require energy_reference_evidence with description and file.")
        evidence_path = (config_dir / evidence["file"]).resolve()
        if not evidence_path.is_file():
            raise DataError("Energy reference evidence file is missing")
    if any(k not in EDGES for k in shifts) or not all(np.isfinite(float(v)) for v in shifts.values()):
        raise DataError("Energy calibration must be a finite single value per edge")
    return export, comparisons


def mean_path(result_dir, model, edge, mode, name):
    path = result_dir / "means" / slug(model) / edge / mode / name
    if not path.is_file():
        raise DataError(f"Missing processed table for {model}/{edge}/{mode}/{name}")
    return path


def common_grid(arrays, window=None):
    if any(np.any(np.diff(a[:, 0]) <= 0) for a in arrays):
        raise DataError("Spectral axes must be strictly increasing")
    low = max(a[0, 0] for a in arrays)
    high = min(a[-1, 0] for a in arrays)
    if window is not None:
        lo, hi = map(float, window)
        if not np.isfinite([lo, hi]).all() or hi <= lo:
            raise DataError("Invalid requested plotting window")
        if low > lo + 1e-7 or high < hi - 1e-7:
            raise DataError(f"Actual spectra cover only {low:.4f}–{high:.4f}; requested {lo:g}–{hi:g}. Recalculate; no padding or extrapolation.")
        low, high = lo, hi
    step = max(float(np.median(np.diff(a[:, 0]))) for a in arrays)
    if high <= low or step <= 0:
        raise DataError("No common spectral interval")
    # Include requested endpoints using interpolation strictly inside actual data.
    count = int(np.ceil((high - low) / step)) + 1
    return np.linspace(low, high, count)


def export_all_tables(result_dir, output, inventory):
    records = []
    # Means and individual absorber data are both retained in tab-delimited TXT.
    for relative in sorted(inventory):
        rel = Path(relative)
        if rel.suffix != ".csv" or rel.parts[0] not in {"means", "sites"}:
            continue
        source = result_dir / rel
        columns, array = load_csv(source)
        target = output / "txt" / rel.with_suffix(".txt")
        save_txt(target, columns, array, [f"source={rel.as_posix()}", f"source_sha256={sha(source)}",
                 "Numeric values retained; no peak shifts, smoothing, or noise added."])
        records.append({"source": rel.as_posix(), "output": target.relative_to(output).as_posix(),
                        "rows": len(array), "columns": len(columns)})
    # Long-format wavelets for ordinary plotting software; matrix TXT files above
    # remain available as well. Restrict duplicated long form to model averages.
    for source in sorted((result_dir / "means").rglob("complex_transforms.npz")):
        if source.relative_to(result_dir).as_posix() not in inventory:
            raise DataError("Wavelet NPZ is absent from verified inventory")
        with np.load(source, allow_pickle=False) as z:
            k, radius = z["WT_k_A_inv"], z["WT_R_A"]
            value = z["WT_Cauchy_complex"]
            if value.shape != (len(radius), len(k)) or not np.isfinite(value).all():
                raise DataError("Invalid WT array shape or values")
            kk, rr = np.meshgrid(k, radius)
            array = np.column_stack([kk.ravel(), rr.ravel(), value.real.ravel(), value.imag.ravel(), np.abs(value).ravel()])
        rel = source.relative_to(result_dir).with_name("WT_Cauchy_long.txt")
        target = output / "txt" / rel
        save_txt(target, ["k_A_inv", "R_A_phase_uncorrected", "WT_real", "WT_imag", "WT_magnitude"], array,
                 [f"source={source.relative_to(result_dir).as_posix()}", f"source_sha256={sha(source)}",
                  "Cauchy WT; no per-map normalization; boundary effects remain."])
        records.append({"source": source.relative_to(result_dir).as_posix(), "output": target.relative_to(output).as_posix(),
                        "rows": len(array), "columns": 5})
    return records


def save_figure(fig, stem, settings, footer=None):
    footer = footer or "Calculated structural models; no electrochemical-potential assignment"
    fig.text(0.5, 0.012, footer, ha="center", va="bottom", fontsize=6, color="#555555")
    stem.parent.mkdir(parents=True, exist_ok=True)
    for fmt in FORMATS:
        kwargs = {}
        if fmt in {"png", "tiff"}:
            kwargs["dpi"] = int(settings.get("dpi", 600))
        if fmt == "tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(stem.with_suffix("." + fmt), **kwargs)
    plt.close(fig)


def stacked_figure(entries, comparison, mode, output, settings):
    fig, axes = plt.subplots(3, 1, figsize=(110 / 25.4, 177 / 25.4))
    for row, edge in enumerate(EDGES):
        entry = entries[edge]
        for i, item in enumerate(comparison["series"]):
            axes[row].plot(entry["x"], entry["y"][:, i], color=COLORS[i], linestyle=STYLES[i],
                           linewidth=0.5, label=item["label"])
        axes[row].text(0.02, 0.92, EDGE_TITLES[edge], transform=axes[row].transAxes, va="top", fontsize=8)
        axes[row].text(-0.15, 1.02, "abc"[row], transform=axes[row].transAxes, fontweight="bold", fontsize=8)
        if mode == "XANES":
            axes[row].set(xlabel="Photon energy (eV)", ylabel=r"Calculated $\mu$ (FEFF units)")
            axes[row].ticklabel_format(axis="x", style="plain", useOffset=False)
        elif mode == "FT":
            axes[row].set(xlabel=r"$R$ ($\mathrm{\AA}$), phase uncorrected", ylabel=r"$|\chi(R)|$ ($\mathrm{\AA}^{-3}$)", xlim=(0, 6))
        else:
            axes[row].set(xlabel=r"$k$ ($\mathrm{\AA}^{-1}$)", ylabel=r"$k^3\chi(k)$ ($\mathrm{\AA}^{-3}$)", xlim=(3, 12))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.52, 0.95), ncol=1, fontsize=7)
    title = comparison.get("title", comparison["id"])
    fig.suptitle("\n".join(textwrap.wrap(title, 65)), fontsize=8, y=0.99)
    fig.subplots_adjust(left=0.18, right=0.97, bottom=0.09, top=0.84, hspace=0.52)
    footer = None if mode == "XANES" else "Calculated structural models; FT distances are not direct bond lengths"
    save_figure(fig, output / "figures" / f"{slug(comparison['id'])}_{mode}_W_Co_Fe", settings, footer)


def load_comparison(result_dir, comparison, settings):
    xanes, ft, oscillations, peak_rows, metrics = {}, {}, {}, [], []
    window_w = settings.get("w_energy_window_eV", [10100.0, 10300.0])
    for edge in EDGES:
        shift = float(settings.get("energy_reference_shift_eV", {}).get(edge, 0.0))
        arrays = []
        for item in comparison["series"]:
            _, table = load_csv(mean_path(result_dir, item["model"], edge, "xanes", "mean_mu.csv"), 4)
            table = table.copy(); table[:, 0] += shift
            arrays.append(table)
        window = window_w if edge == "W_L3" else settings.get("energy_window_eV_by_edge", {}).get(edge)
        x = common_grid(arrays, window)
        y = np.column_stack([np.interp(x, a[:, 0], a[:, 1]) for a in arrays])
        xanes[edge] = {"x": x, "y": y, "shift": shift}
        reference_scale = float(np.sqrt(np.mean(y[:, 0]**2)))
        for i in (1, 2):
            rms = float(np.sqrt(np.mean((y[:, i] - y[:, 0])**2)))
            metrics.append({"edge": edge, "model": comparison["series"][i]["model"],
                            "reference": comparison["series"][0]["model"], "window_eV": [float(x[0]), float(x[-1])],
                            "raw_mu_rms_difference": rms,
                            "relative_to_reference_rms": rms/reference_scale if reference_scale > 1e-15 else None,
                            "normalization": "None; same-edge native FEFF mu, common energy grid"})
        raw_ft = [load_csv(mean_path(result_dir, item["model"], edge, "exafs", "FT_kweight2.csv"), 4)[1]
                  for item in comparison["series"]]
        r = common_grid(raw_ft, [0.0, min(6.0, min(a[-1, 0] for a in raw_ft))])
        ft[edge] = {"x": r, "y": np.column_stack([np.interp(r, a[:, 0], a[:, 3]) for a in raw_ft])}
        raw_chi = [load_csv(mean_path(result_dir, item["model"], edge, "exafs", "chi_common_3_to_12.csv"), 5)[1]
                   for item in comparison["series"]]
        k = common_grid(raw_chi, [3.0, 12.0])
        oscillations[edge] = {"x": k, "y": np.column_stack([np.interp(k, a[:, 0], a[:, 4]) for a in raw_chi])}
        if edge == "W_L3":
            for item, a in zip(comparison["series"], raw_ft):
                keep = (a[:, 0] >= 1.0) & (a[:, 0] <= 2.5)
                if not keep.any():
                    raise DataError("Missing W FT first-shell display region")
                subset = a[keep]; idx = int(np.argmax(subset[:, 3]))
                peak_rows.append({"model": item["model"], "R_peak_A_phase_uncorrected": float(subset[idx, 0]),
                                  "FT_magnitude": float(subset[idx, 3]), "search_range_A": [1.0, 2.5],
                                  "at_search_boundary": idx in {0, len(subset)-1},
                                  "interpretation": "Grid maximum only; no target enforced and no bond-length/coordination fit"})
    return xanes, ft, oscillations, peak_rows, metrics


def wavelet_figure(result_dir, comparison, output, settings):
    values = {}
    for edge in EDGES:
        values[edge] = []
        for item in comparison["series"]:
            path = mean_path(result_dir, item["model"], edge, "exafs", "complex_transforms.npz")
            with np.load(path, allow_pickle=False) as z:
                k, r, wavelet = z["WT_k_A_inv"], z["WT_R_A"], z["WT_Cauchy_complex"]
                keep = (k >= 3.0 - 1e-8) & (k <= 12.0 + 1e-8)
                values[edge].append((k[keep], r, np.abs(wavelet[:, keep])))
    fig, axes = plt.subplots(3, 3, figsize=(183/25.4, 175/25.4))
    scales = {}
    for row, edge in enumerate(EDGES):
        maximum = max(float(v[2].max()) for v in values[edge]) or 1.0
        scales[edge] = [0.0, maximum]
        for col, (k, r, value) in enumerate(values[edge]):
            ax = axes[row, col]
            image = ax.pcolormesh(k, r, value, shading="auto", cmap="viridis",
                                 norm=Normalize(0, maximum), rasterized=True)
            ax.set(xlabel=r"$k$ ($\mathrm{\AA}^{-1}$)", ylabel=r"$R$ ($\mathrm{\AA}$)", ylim=(0, 6))
            ax.text(-0.19, 1.05, "abcdefghi"[row*3+col], transform=ax.transAxes, fontweight="bold", fontsize=8)
            if row == 0:
                ax.set_title("\n".join(textwrap.wrap(comparison["series"][col]["label"], 24)), fontsize=7)
            if col == 0:
                ax.text(0.02, 0.94, EDGE_TITLES[edge], transform=ax.transAxes, color="white", va="top", fontsize=7)
        cax = fig.add_axes([0.91, [0.65, 0.365, 0.08][row], 0.012, 0.20])
        fig.colorbar(image, cax=cax, label="WT magnitude (a.u.)")
    fig.subplots_adjust(left=0.08, right=0.87, bottom=0.08, top=0.91, wspace=0.43, hspace=0.48)
    save_figure(fig, output/"figures"/f"{slug(comparison['id'])}_WT_W_Co_Fe", settings,
                "Calculated Cauchy WT; common scale within each edge; phase uncorrected R")
    return scales


def run(config_path, result_dir, output_dir=None):
    config_path = Path(config_path).resolve(); result_dir = Path(result_dir).resolve()
    config = read_json(config_path)
    report, inventory = validate_source(result_dir)
    settings, comparisons = validate_comparisons(config, report, config_path.parent)
    # Validate all requested means/windows before writing any final output.
    prepared = [(c, load_comparison(result_dir, c, settings)) for c in comparisons]
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output = Path(output_dir).resolve() if output_dir else result_dir / ("FINAL_EXPORT_" + timestamp)
    output.mkdir(parents=True, exist_ok=False)
    fonts = configure_font()
    records = export_all_tables(result_dir, output, inventory)
    comparison_reports = []
    for comparison, (xanes, ft, oscillations, peaks, metrics) in prepared:
        comp_dir = output / "txt" / "comparisons" / slug(comparison["id"])
        for name, collection, axis in [("XANES", xanes, "photon_energy_eV"),
                                       ("FT_kweight2", ft, "R_A_phase_uncorrected"),
                                       ("k3_chi", oscillations, "k_A_inv")]:
            for edge, entry in collection.items():
                save_txt(comp_dir / f"{edge}_{name}.txt", [axis] + [s["model"] for s in comparison["series"]],
                         np.column_stack([entry["x"], entry["y"]]),
                         [f"comparison={comparison['id']}", "Common plotting grid; linear interpolation inside actual coverage only.",
                          f"same_edge_energy_shift_eV={entry.get('shift',0.0)}" if name == "XANES" else "No phase correction or R-axis shift."])
        write_json(comp_dir / "W_first_shell_peak_report.json", peaks)
        save_txt(comp_dir / "W_first_shell_peak_report.txt", ["series_index", "R_peak_A_phase_uncorrected", "FT_magnitude"],
                 [[i, p["R_peak_A_phase_uncorrected"], p["FT_magnitude"]] for i, p in enumerate(peaks)],
                 ["model_order=" + ", ".join(p["model"] for p in peaks),
                  "Grid maximum in 1.0-2.5 A only; not a measured/fitted W-O bond length; no target peak imposed."])
        save_txt(comp_dir / "XANES_difference_metrics.txt",
                 ["edge_index", "series_index", "raw_mu_RMS_difference", "reference_mu_RMS", "energy_min_eV", "energy_max_eV"],
                 [[EDGES.index(v["edge"]), next(i for i,s in enumerate(comparison["series"]) if s["model"] == v["model"]),
                   v["raw_mu_rms_difference"], float(np.sqrt(np.mean(xanes[v["edge"]]["y"][:, 0]**2))),
                   *v["window_eV"]] for v in metrics],
                 ["edge_order=" + ", ".join(EDGES),
                  "model_order=" + ", ".join(s["model"] for s in comparison["series"]),
                  "Reference is series 0. No separate curve normalization or oxidation-state/active-site assignment."])
        for mode, data in [("XANES", xanes), ("FT", ft), ("k3_EXAFS", oscillations)]:
            stacked_figure(data, comparison, mode, output, settings)
        scales = wavelet_figure(result_dir, comparison, output, settings)
        comparison_reports.append({"id": comparison["id"], "series": comparison["series"],
                                   "W_first_shell_region_peaks": peaks, "XANES_differences": metrics,
                                   "WT_scales_by_edge": scales})
    summary = {"status": "complete", "data_kind": "Calculated FEFF10 structural-model spectra",
               "source_result_directory_name": result_dir.name,
               "source_postprocess_report_sha256": sha(result_dir/"postprocess_report.json"),
               "config_sha256": sha(config_path), "font": fonts, "formats": list(FORMATS),
               "txt_tables": records, "comparisons": comparison_reports,
               "normalization": "None; all raw means retained. No individual peak alignment, smoothing or noise added.",
               "energy_reference_shift_eV": settings.get("energy_reference_shift_eV", {}),
               "energy_reference_evidence": settings.get("energy_reference_evidence"),
               "limits": ["No direct electrode-potential mapping from static structures.",
                          "Changes in spectra do not prove an exclusive active site.",
                          "FT R is phase uncorrected and is not a direct W-O bond distance.",
                          "Cauchy WT is a transform of the same EXAFS signal, not independent experimental evidence."]}
    write_json(output/"export_report.json", summary)
    files = [{"path": p.relative_to(output).as_posix(), "bytes": p.stat().st_size, "sha256": sha(p)}
             for p in sorted(output.rglob("*")) if p.is_file()]
    write_json(output/"export_inventory.json", files)
    print(json.dumps({"status":"complete", "output": str(output), "txt_tables":len(records),
                      "comparisons":len(comparison_reports)}, ensure_ascii=False))
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--result-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    try:
        run(args.config, args.result_dir, args.output_dir)
    except (DataError, OSError, ValueError, KeyError) as exc:
        print("EXPORT STOPPED: " + str(exc), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
