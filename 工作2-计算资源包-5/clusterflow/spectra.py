"""Auditable, noise-free CALCULATED NMR browsing spectra.

Small isolated 1H blocks use a spin-1/2 isotropic Hamiltonian, diagonalized
within conserved-magnetization sectors. F-containing or oversized H blocks
are explicitly first order. 13C{1H} sites are independent isotopomers with
first-order C-F splitting, not a molecule containing 100% 13C.

Figure contract: display positions and the declared model's fine structure;
one quantitative panel, Python/matplotlib only, editable SVG plus PNG and
source CSVs. Lorentzian width is a user-chosen display parameter. No noise,
relaxation prediction, exchange kinetics, isotope-shift claim, or experimental
13C intensity model is introduced here. Fixtures used in tests are synthetic
engineering cases, not molecular predictions.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import textwrap

import numpy as np


def _finite(value, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite number") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    return result


def _setting(settings: dict, names: tuple[str, ...], default):
    for name in names:
        if name in settings:
            return settings[name]
    return default


def _merge(lines: list[dict], decimals: int = 10) -> list[dict]:
    merged: dict[float, dict] = {}
    for row in lines:
        key = round(row["ppm"], decimals)
        if key not in merged:
            merged[key] = {"ppm": key, "intensity": 0.0, "source_atom_indices": set(),
                           "models": set()}
        merged[key]["intensity"] += float(row["intensity"])
        merged[key]["source_atom_indices"].update(row["source_atom_indices"])
        merged[key]["models"].add(row["model"])
    return [{"ppm": row["ppm"], "intensity": row["intensity"],
             "source_atom_indices": sorted(row["source_atom_indices"]),
             "model": "; ".join(sorted(row["models"]))}
            for _, row in sorted(merged.items())]


def _components(ids: list[int], couplings: dict[tuple[int, int], float]) -> list[list[int]]:
    neighbors = {i: set() for i in ids}
    for i, j in couplings:
        if i in neighbors and j in neighbors:
            neighbors[i].add(j)
            neighbors[j].add(i)
    unseen, components = set(ids), []
    while unseen:
        stack, component = [min(unseen)], []
        while stack:
            i = stack.pop()
            if i not in unseen:
                continue
            unseen.remove(i)
            component.append(i)
            stack.extend(sorted(neighbors[i] & unseen, reverse=True))
        components.append(sorted(component))
    return components


def _exact_h(component: list[int], atoms: dict[int, dict], pairs: dict, field: float) -> list[dict]:
    """High-temperature total transverse absorption; exact within this H block.

    I_z=+1/2 for bit=1. Isotropic J contributes J*Iz*Iz and J/2
    flip-flops. Only adjacent magnetization sectors are connected by sum I-.
    Relative energies retain their signs, so negative chemical shifts survive.
    """
    n = len(component)
    if n == 1:
        return [{"ppm": atoms[component[0]]["shift_ppm"], "intensity": 1.0,
                 "source_atom_indices": component, "model": "isolated singlet"}]
    ppm = np.array([atoms[i]["shift_ppm"] for i in component], dtype=float)
    carrier = float(ppm.mean() * field)
    offsets = ppm * field - carrier
    local = {atom: k for k, atom in enumerate(component)}
    local_pairs = [(local[i], local[j], value) for (i, j), value in pairs.items()
                   if i in local and j in local]
    sectors = []
    for k in range(n + 1):
        states = [state for state in range(1 << n) if state.bit_count() == k]
        index = {state: position for position, state in enumerate(states)}
        hamiltonian = np.zeros((len(states), len(states)), dtype=float)
        for row, state in enumerate(states):
            z = np.array([0.5 if state & (1 << i) else -0.5 for i in range(n)])
            diagonal = float(np.dot(offsets, z))
            for i, j, coupling in local_pairs:
                diagonal += coupling * z[i] * z[j]
                if bool(state & (1 << i)) != bool(state & (1 << j)):
                    hamiltonian[row, index[state ^ (1 << i) ^ (1 << j)]] += coupling / 2
            hamiltonian[row, row] += diagonal
        energies, vectors = np.linalg.eigh(hamiltonian)
        sectors.append((states, index, energies, vectors))
    lines = []
    for k in range(1, n + 1):
        states, _, high_e, high_v = sectors[k]
        low_states, low_index, low_e, low_v = sectors[k - 1]
        lowering = np.zeros((len(low_states), len(states)), dtype=float)
        for column, state in enumerate(states):
            for i in range(n):
                if state & (1 << i):
                    lowering[low_index[state ^ (1 << i)], column] += 1.0
        intensity = np.square(low_v.T @ lowering @ high_v)
        frequencies = high_e[np.newaxis, :] - low_e[:, np.newaxis] + carrier
        rows, columns = np.nonzero(intensity > 1e-12)
        lines.extend({"ppm": float(frequencies[r, c] / field),
                      "intensity": float(intensity[r, c]),
                      "source_atom_indices": component, "model": "exact H-H block"}
                     for r, c in zip(rows, columns))
    area = sum(line["intensity"] for line in lines)
    if area <= 0:
        raise ValueError("Exact spin block produced no positive transitions")
    for line in lines:
        line["intensity"] *= n / area
    return lines


def _magnetically_equivalent(i: int, j: int, atoms: dict, pairs: dict) -> bool:
    # Equality here is algebraic equivalence within the supplied parameter
    # model, not an assertion of experimental or graph-based equivalence.
    if abs(atoms[i]["shift_ppm"] - atoms[j]["shift_ppm"]) > 1e-10:
        return False
    for other in atoms:
        if other in (i, j) or atoms[other]["element"] not in ("H", "F"):
            continue
        if abs(pairs.get(tuple(sorted((i, other))), 0.0)
               - pairs.get(tuple(sorted((j, other))), 0.0)) > 1e-10:
            return False
    return True


def _first_order(atom: int, couplings: list[float], atoms: dict, field: float,
                 bin_hz: float, model: str) -> list[dict]:
    # Each spin-1/2 partner contributes ±J/2 at equal population. Fine bins
    # bound combinatorial line growth; their size is recorded in metadata.
    distribution = {0: 1.0}
    for coupling in sorted(couplings, key=lambda v: (-abs(v), v)):
        new = {}
        for key, intensity in distribution.items():
            for sign in (-1, 1):
                position = round((key * bin_hz + sign * coupling / 2) / bin_hz)
                new[position] = new.get(position, 0.0) + intensity / 2
        if len(new) > 200000:
            raise ValueError("First-order line count exceeds the explicit safety limit")
        distribution = new
    return [{"ppm": atoms[atom]["shift_ppm"] + offset * bin_hz / field,
             "intensity": intensity, "source_atom_indices": [atom], "model": model}
            for offset, intensity in sorted(distribution.items())]


def make_spectrum(atom_rows, j_pairs, nucleus, settings) -> dict:
    """Return a JSON-serializable spectrum with per-observed-nucleus unit area.

    atom_rows must include metadata for every J endpoint (including F with
    shift_ppm=None). Observed H/C shifts must be finite. J keys follow orcaio:
    atom_i_zero_based, atom_j_zero_based, J_Hz. Unknown/duplicate endpoints
    fail instead of silently losing a coupling. Both protocol.json spectrum
    keys and field_MHz_H/C, max_exact_H_spins aliases are accepted.
    """
    symbol = {"H": "H", "1H": "H", "C": "C", "13C": "C"}.get(str(nucleus))
    if symbol is None:
        raise ValueError("Only H/1H and C/13C spectra are supported")
    settings = dict(settings or {})
    if isinstance(settings.get("spectrum"), dict):
        settings = {**settings, **settings["spectrum"]}
    if settings.get("noise", False):
        raise ValueError("Synthetic noise is not supported by the calculated-spectrum layer")
    field = _finite(_setting(settings, (f"field_MHz_{symbol}",
                    "proton_frequency_MHz" if symbol == "H" else "carbon_frequency_MHz"),
                    400.0 if symbol == "H" else 100.61), "field MHz")
    linewidth = _finite(_setting(settings, (f"linewidth_Hz_{symbol}", "linewidth_Hz",
                        "proton_linewidth_Hz" if symbol == "H" else "carbon_linewidth_Hz"),
                        1.0 if symbol == "H" else 1.5), "linewidth Hz")
    threshold = _finite(_setting(settings, ("threshold_Hz", "J_graph_threshold_Hz"), 0.1), "J threshold")
    maximum = _setting(settings, ("max_exact_H_spins", "max_exact_proton_spins"), 10)
    if not isinstance(maximum, int) or isinstance(maximum, bool) or not 0 <= maximum <= 10:
        raise ValueError("max_exact_H_spins must be an integer between 0 and 10")
    bin_hz = _finite(settings.get("first_order_line_bin_Hz", min(0.01, linewidth / 100)), "line bin Hz")
    if field <= 0 or linewidth <= 0 or threshold < 0 or bin_hz <= 0 or bin_hz > linewidth / 20:
        raise ValueError("Invalid field, linewidth, J threshold or first-order line-bin size")
    atoms, original_atoms, warnings = {}, [], []
    for raw in atom_rows:
        row = dict(raw)
        atom = row.get("atom_index_zero_based")
        if not isinstance(atom, int) or isinstance(atom, bool) or atom < 0 or atom in atoms:
            raise ValueError("Atom indices must be unique nonnegative integers")
        row["element"] = str(row.get("element", ""))
        if not row["element"]:
            raise ValueError("An atom element is required")
        if row["element"] == symbol:
            row["shift_ppm"] = _finite(row.get("shift_ppm"), "observed chemical shift")
        # Validate other H shifts too: H-model equivalence uses them only for H.
        if symbol == "H" and row["element"] == "H":
            row["shift_ppm"] = _finite(row.get("shift_ppm"), "proton chemical shift")
        if not isinstance(row.get("exchangeable_H", False), bool):
            raise ValueError("exchangeable_H must be a boolean")
        row["exchangeable_H"] = row.get("exchangeable_H", False)
        atoms[atom] = row
        original_atoms.append(row)
    observed = sorted(i for i, row in atoms.items() if row["element"] == symbol)
    if not observed:
        raise ValueError(f"No observed {symbol} atoms were supplied")
    pairs, original_pairs, seen_pairs, omitted_exchange = {}, [], set(), 0
    for raw in j_pairs:
        row = dict(raw)
        i, j = row.get("atom_i_zero_based"), row.get("atom_j_zero_based")
        if (any(not isinstance(index, int) or isinstance(index, bool) for index in (i, j))
                or i not in atoms or j not in atoms or i == j):
            raise ValueError("Every J endpoint must identify a distinct atom_rows entry")
        pair = tuple(sorted((i, j)))
        if pair in seen_pairs:
            raise ValueError("Duplicate J pair; do not count a symmetric pair twice")
        seen_pairs.add(pair)
        coupling = _finite(row.get("J_Hz"), "J Hz")
        row["J_Hz"] = coupling
        original_pairs.append(row)
        if any(atoms[k]["element"] == "H" and atoms[k]["exchangeable_H"] for k in pair):
            omitted_exchange += 1
            continue
        if abs(coupling) >= threshold and coupling != 0:
            pairs[pair] = coupling
    if omitted_exchange:
        warnings.append(f"Excluded {omitted_exchange} J pairs involving exchangeable H; its shift peak is retained, without an exchange/relaxation model.")
    lines, methods, block_report = [], set(), []
    if symbol == "H":
        hh = {pair: value for pair, value in pairs.items()
              if all(atoms[i]["element"] == "H" for i in pair)}
        hf_atoms = {i for pair in pairs for i in pair if atoms[i]["element"] == "H"
                    and any(atoms[j]["element"] == "F" for j in pair)}
        for component in _components(observed, hh):
            exact = len(component) <= maximum and not (set(component) & hf_atoms)
            if exact:
                lines.extend(_exact_h(component, atoms, hh, field))
                method = "exact H-H spin-1/2 block"
            else:
                reason = "H-F present" if set(component) & hf_atoms else "H block exceeds exact-size limit"
                warnings.append(f"H block {component}: {reason}; H-H/H-F rendered first order, not a full strong-coupling spectrum.")
                method = "first-order H-H/H-F approximation"
                for atom in component:
                    coupling_list = []
                    for pair, value in pairs.items():
                        if atom not in pair:
                            continue
                        other = pair[1] if pair[0] == atom else pair[0]
                        element = atoms[other]["element"]
                        if element == "F" or (element == "H" and not _magnetically_equivalent(atom, other, atoms, pairs)):
                            coupling_list.append(value)
                    lines.extend(_first_order(atom, coupling_list, atoms, field, bin_hz, method))
            methods.add(method)
            block_report.append({"atom_indices_zero_based": component, "method": method})
        warnings.append("13C satellites and couplings to quadrupolar nuclei are not modelled in the proton spectrum.")
    else:
        method = "13C{1H}: independent sites; first-order C-F only"
        methods.add(method)
        for atom in observed:
            values = []
            for pair, value in pairs.items():
                if atom in pair:
                    other = pair[1] if pair[0] == atom else pair[0]
                    if atoms[other]["element"] == "F":
                        values.append(value)
            lines.extend(_first_order(atom, values, atoms, field, bin_hz, method))
        warnings.append("13C sites have idealized unit areas; experimental relaxation/NOE intensities and 13C-13C satellites are not predicted.")
        if any(atoms[i]["element"] == "F" for pair in pairs for i in pair):
            warnings.append("C-F splitting is first order; F-F strong coupling and 19F decoupling are not modelled.")
    sticks = _merge(lines)
    total = sum(row["intensity"] for row in sticks)
    if not math.isclose(total, len(observed), rel_tol=1e-10, abs_tol=1e-10):
        raise ValueError("Spectrum area conservation failed")
    fwhm_ppm = linewidth / field
    padding = max(0.5, 40 * fwhm_ppm)
    limits = settings.get("x_range_ppm", (sticks[0]["ppm"] - padding, sticks[-1]["ppm"] + padding))
    if len(limits) != 2:
        raise ValueError("x_range_ppm must contain exactly two limits")
    low, high = sorted(_finite(v, "ppm limit") for v in limits)
    if low >= high or low > sticks[0]["ppm"] or high < sticks[-1]["ppm"]:
        raise ValueError("ppm range must include all stick positions and have positive width")
    requested = settings.get("n_points", 16384)
    if not isinstance(requested, int) or isinstance(requested, bool) or requested < 128:
        raise ValueError("n_points must be an integer >= 128")
    # At least eight samples/FWHM. FFT convolution uses linear stick deposition;
    # it does not snap sticks in the source-data output to a plotting grid.
    points = max(requested, int(math.ceil((high - low) / (fwhm_ppm / 8))) + 1)
    if points > 500000:
        raise ValueError("Requested linewidth/range needs >500000 samples; narrow the range or increase linewidth")
    x = np.linspace(low, high, points)
    step = float(x[1] - x[0])
    impulse = np.zeros(points)
    for row in sticks:
        location = (row["ppm"] - low) / step
        left = min(points - 2, max(0, int(math.floor(location))))
        fraction = min(1.0, max(0.0, location - left))
        impulse[left] += row["intensity"] * (1 - fraction)
        impulse[left + 1] += row["intensity"] * fraction
    distance = np.arange(-(points - 1), points) * step
    halfwidth = fwhm_ppm / 2
    kernel = halfwidth / (math.pi * (distance * distance + halfwidth * halfwidth))
    convolution_length = len(impulse) + len(kernel) - 1
    fft_length = 1 << (convolution_length - 1).bit_length()
    full = np.fft.irfft(np.fft.rfft(impulse, fft_length) * np.fft.rfft(kernel, fft_length), fft_length)
    y = np.maximum(full[points - 1:2 * points - 1], 0.0)
    # Explicit trapezoidal integration also supports the declared numpy>=1.26.
    sampled_area = float(np.sum((y[:-1] + y[1:]) * np.diff(x) / 2))
    if sampled_area <= 0:
        raise ValueError("Lorentzian trace has zero area")
    y *= len(observed) / sampled_area
    return {"sticks": sticks, "x_ppm": x.tolist(), "y": y.tolist(),
            "multiplet_readiness": ("APPROXIMATION_REVIEW_REQUIRED" if any('first-order' in m for m in methods) and symbol == 'H'
                                    else "DECLARED_MODEL_REQUIRES_THRESHOLD_AND_EXPERIMENTAL_REVIEW"),
            "full_strong_coupling_for_all_64_verified": False,
            "warnings": warnings, "method": "; ".join(sorted(methods)),
            "nucleus": "1H" if symbol == "H" else "13C", "n_observed_nuclei": len(observed),
            "data_status": "CALCULATED_MODEL_NOT_EXPERIMENTAL",
            "field_MHz": field, "linewidth_FWHM_Hz": linewidth,
            "threshold_Hz": threshold, "first_order_line_bin_Hz": bin_hz,
            "exact_transition_intensity_floor": 1e-12,
            "normalization": "stick area and finite-window trace area equal number of observed nuclei",
            "broadening": "Lorentzian; FFT convolution with linear grid deposition; >=8 grid points/FWHM; no noise",
            "block_report": block_report, "source_atom_rows": original_atoms,
            "source_j_pairs": original_pairs}


def _write_rows(path: Path, rows: list[dict], columns: list[str]) -> None:
    from .common import _io_path
    with _io_path(path).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(row.get(key), ensure_ascii=False)
                             if isinstance(row.get(key), (list, dict)) else row.get(key)
                             for key in columns})


def export_spectrum(result, out_prefix, title) -> dict[str, str]:
    """Save transparent-provenance SVG/PNG, trace/sticks/raw-source CSV and JSON.

    No outputs are overwritten: choose a fresh prefix for a rerun. Metadata
    accompanies the figure; the title/footer always identify a calculated model.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .common import _io_path

    prefix = Path(out_prefix)
    outputs = {key: Path(str(prefix) + suffix) for key, suffix in {
        "svg": ".svg", "png": ".png", "trace_csv": "_trace.csv",
        "sticks_csv": "_sticks.csv", "source_atoms_csv": "_source_atoms.csv",
        "source_j_csv": "_source_J.csv", "metadata_json": "_metadata.json"}.items()}
    if any(_io_path(path).exists() for path in outputs.values()):
        raise FileExistsError("Spectrum export refuses to overwrite an existing artifact prefix")
    _io_path(prefix.parent).mkdir(parents=True, exist_ok=True)
    x, y = np.array(result["x_ppm"], dtype=float), np.array(result["y"], dtype=float)
    if x.ndim != 1 or x.shape != y.shape or len(x) < 2 or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Spectrum trace must contain matched finite one-dimensional arrays")
    provenance = {"data_status": "CALCULATED_MODEL_NOT_EXPERIMENTAL", "nucleus": result["nucleus"],
                  "field_MHz": result["field_MHz"], "linewidth_FWHM_Hz": result["linewidth_FWHM_Hz"],
                  "method": result["method"]}
    _write_rows(outputs["trace_csv"], [{"ppm": float(a), "intensity_per_ppm": float(b), **provenance}
                                     for a, b in zip(x, y)], ["ppm", "intensity_per_ppm", *provenance])
    _write_rows(outputs["sticks_csv"], [{**row, **provenance, "stick_model": row["model"]}
                                      for row in result["sticks"]],
                ["ppm", "intensity", "source_atom_indices", "stick_model", *provenance])
    for key, rows in (("source_atoms_csv", result["source_atom_rows"]), ("source_j_csv", result["source_j_pairs"])):
        columns = sorted(set().union(*(row.keys() for row in rows))) if rows else ["no_source_rows"]
        _write_rows(outputs[key], rows, columns)
    metadata = {key: value for key, value in result.items() if key not in ("x_ppm", "y", "source_atom_rows", "source_j_pairs")}
    metadata["title"] = str(title)
    metadata["figure_contract"] = "One quantitative calculated-spectrum browsing panel; no experiment/accuracy assertion; SVG text editable, CSV source values preserved."
    _io_path(outputs["metadata_json"]).write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    style = {"font.family": "sans-serif", "font.sans-serif": ["Arial", "Microsoft YaHei", "DejaVu Sans"],
             "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 9,
             "axes.spines.right": False, "axes.spines.top": False, "axes.linewidth": 0.8,
             "text.usetex": False}
    with plt.rc_context(style):
        fig, ax = plt.subplots(figsize=(7.2, 3.2))
        fig.subplots_adjust(left=0.105, right=0.98, top=0.76, bottom=0.27)
        ax.plot(x, y, color="#204F70", linewidth=0.85)
        ax.set_xlim(float(x[-1]), float(x[0]))
        ax.set_ylim(0, float(y.max()) * 1.08)
        ax.set_xlabel("Chemical shift / ppm")
        ax.set_ylabel("Calculated intensity / a.u.")
        ax.set_yticks([])
        ax.tick_params(direction="out", width=0.7, length=3)
        fig.text(0.105, 0.93, "CALCULATED | " + "\n".join(textwrap.wrap(str(title), 82)),
                 ha="left", va="top", fontsize=10, weight="bold")
        fig.text(0.105, 0.16, f"{result['nucleus']} | {result['field_MHz']:g} MHz | Lorentzian FWHM {result['linewidth_FWHM_Hz']:g} Hz | no noise",
                 ha="left", va="top", fontsize=8)
        fig.text(0.105, 0.095, "\n".join(textwrap.wrap(result["method"], 107)),
                 ha="left", va="top", fontsize=7.5, color="#4C5961")
        fig.savefig(_io_path(outputs["svg"]), facecolor="white")
        fig.savefig(_io_path(outputs["png"]), facecolor="white", dpi=300)
        plt.close(fig)
    return {key: str(value) for key, value in outputs.items()}
