#!/usr/bin/env python3
"""Render a transparent, condition-dependent NMR *appearance* layer.

The input is a clean transition/line table.  ORCA is not called here.  This
module adds line broadening, named solvent/water components, a deterministic
baseline and seeded noise.  It therefore models how a calculated line list
might look under a declared acquisition profile; it is not a first-principles
prediction of an acquired spectrum.

Required CSV columns: figure_id,nucleus,ppm,intensity
Optional columns: linewidth_hz,exchangeable_H_flag,multiplicity,J_hz
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import wofz


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT / "config" / "appearance_profiles.json"
DISCLAIMER = (
    "SYNTHETIC, CONDITIONED APPEARANCE - NOT AN ACQUIRED SPECTRUM. "
    "Noise, baseline, solvent intensity, peak height and linewidth are not "
    "static-ORCA first-principles observables."
)


class AppearanceError(RuntimeError):
    """Raised for invalid input or an internally inconsistent profile."""


@dataclass(frozen=True)
class CleanLine:
    figure_id: str
    nucleus: str
    ppm: float
    intensity: float
    linewidth_hz: float | None
    exchangeable: bool
    multiplicity: str
    j_hz: tuple[float, ...]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_nucleus(value: str) -> str:
    token = value.strip().upper().replace("^", "").replace(" ", "")
    aliases = {"H": "1H", "H1": "1H", "1H": "1H", "C": "13C", "C13": "13C", "13C": "13C"}
    try:
        return aliases[token]
    except KeyError as exc:
        raise AppearanceError(f"Unsupported nucleus {value!r}; use 1H or 13C") from exc


def parse_bool(value: str | bool | None) -> bool:
    if isinstance(value, bool):
        return value
    token = str(value or "").strip().lower()
    if token in {"", "0", "false", "no", "n"}:
        return False
    if token in {"1", "true", "yes", "y"}:
        return True
    raise AppearanceError(f"Invalid Boolean value: {value!r}")


def parse_j_values(value: str | None) -> tuple[float, ...]:
    if value is None or not value.strip():
        return ()
    try:
        values = tuple(float(item) for item in re.split(r"[;,/|\s]+", value.strip()) if item)
    except ValueError as exc:
        raise AppearanceError(f"Invalid J_hz value: {value!r}") from exc
    if any(item <= 0 for item in values):
        raise AppearanceError("Every J_hz value must be positive")
    return values


def read_clean_lines(path: Path) -> list[CleanLine]:
    required = {"figure_id", "nucleus", "ppm", "intensity"}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise AppearanceError(f"Clean-line CSV is missing columns: {sorted(missing)}")
        lines: list[CleanLine] = []
        for row_number, row in enumerate(reader, start=2):
            try:
                ppm = float(row["ppm"])
                intensity = float(row["intensity"])
                width_text = (row.get("linewidth_hz") or "").strip()
                linewidth = float(width_text) if width_text else None
            except (TypeError, ValueError) as exc:
                raise AppearanceError(f"Invalid numeric value at CSV row {row_number}") from exc
            if not math.isfinite(ppm) or not math.isfinite(intensity) or intensity < 0:
                raise AppearanceError(f"ppm/intensity is invalid at CSV row {row_number}")
            if linewidth is not None and (not math.isfinite(linewidth) or linewidth <= 0):
                raise AppearanceError(f"linewidth_hz must be positive at CSV row {row_number}")
            figure_id = (row["figure_id"] or "").strip()
            if not figure_id:
                raise AppearanceError(f"figure_id is empty at CSV row {row_number}")
            lines.append(
                CleanLine(
                    figure_id=figure_id,
                    nucleus=normalize_nucleus(row["nucleus"]),
                    ppm=ppm,
                    intensity=intensity,
                    linewidth_hz=linewidth,
                    exchangeable=parse_bool(row.get("exchangeable_H_flag")),
                    multiplicity=(row.get("multiplicity") or "s").strip().lower(),
                    j_hz=parse_j_values(row.get("J_hz") or row.get("j_hz")),
                )
            )
    if not lines:
        raise AppearanceError("Clean-line CSV contains no transitions")
    identities = {(line.figure_id, line.nucleus) for line in lines}
    if len(identities) != 1:
        raise AppearanceError(
            "One invocation must contain exactly one figure_id/nucleus pair; split mixed CSV files first"
        )
    return lines


def load_profiles(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AppearanceError(f"Cannot read appearance configuration {path}: {exc}") from exc
    if not isinstance(payload.get("profiles"), dict) or not payload["profiles"]:
        raise AppearanceError("appearance_profiles.json must contain a non-empty 'profiles' object")
    return payload


def select_profile(payload: dict[str, Any], profile_id: str | None, nucleus: str) -> tuple[str, dict[str, Any]]:
    if profile_id is None:
        profile_id = str(payload.get("default_profile_by_nucleus", {}).get(nucleus, ""))
    if not profile_id or profile_id not in payload["profiles"]:
        raise AppearanceError(f"Unknown or missing appearance profile: {profile_id!r}")
    profile = payload["profiles"][profile_id]
    if normalize_nucleus(str(profile.get("nucleus", ""))) != nucleus:
        raise AppearanceError(f"Profile {profile_id!r} does not match input nucleus {nucleus}")
    status = str(profile.get("status", "")).upper()
    if not status:
        raise AppearanceError(f"Profile {profile_id!r} must declare a provenance status")
    return profile_id, profile


def _binomial_pattern(order: int) -> tuple[np.ndarray, np.ndarray]:
    offsets = np.arange(order, dtype=float) - (order - 1) / 2.0
    weights = np.asarray([math.comb(order - 1, index) for index in range(order)], dtype=float)
    return offsets, weights / weights.sum()


def _multiplicity_orders(label: str, j_count: int) -> list[int]:
    if j_count == 0:
        return []
    token = label.strip().lower().replace(" ", "")
    aliases = {
        "s": 1,
        "singlet": 1,
        "d": 2,
        "doublet": 2,
        "t": 3,
        "triplet": 3,
        "q": 4,
        "quartet": 4,
        "quint": 5,
        "quintet": 5,
        "sext": 6,
        "sextet": 6,
        "sept": 7,
        "septet": 7,
        "m": 1,
        "multiplet": 1,
    }
    if token in aliases:
        order = aliases[token]
        return [] if order == 1 else [order] * j_count
    character_orders = {"d": 2, "t": 3, "q": 4}
    if token and all(character in character_orders for character in token):
        values = [character_orders[character] for character in token]
        if len(values) != j_count:
            raise AppearanceError(
                f"Multiplicity {label!r} describes {len(values)} couplings but J_hz has {j_count}"
            )
        return values
    raise AppearanceError(f"Unsupported multiplicity {label!r}")


def multiplet_components(multiplicity: str, j_values: Iterable[float]) -> tuple[np.ndarray, np.ndarray]:
    j_tuple = tuple(j_values)
    orders = _multiplicity_orders(multiplicity, len(j_tuple))
    offsets = np.asarray([0.0])
    weights = np.asarray([1.0])
    for order, j_hz in zip(orders, j_tuple):
        local_offsets, local_weights = _binomial_pattern(order)
        offsets = (offsets[:, None] + local_offsets[None, :] * j_hz).ravel()
        weights = (weights[:, None] * local_weights[None, :]).ravel()
    return offsets, weights / weights.sum()


def line_profile(
    axis_ppm: np.ndarray,
    center_ppm: float,
    fwhm_ppm: float,
    shape: str,
    voigt_lorentz_fraction: float,
) -> np.ndarray:
    """Return an area-normalized line profile on a ppm axis."""
    if fwhm_ppm <= 0:
        raise AppearanceError("Linewidth must be positive")
    delta = axis_ppm - center_ppm
    shape_token = shape.strip().lower()
    if shape_token in {"lorentz", "lorentzian"}:
        gamma = fwhm_ppm / 2.0
        return gamma / (math.pi * (delta * delta + gamma * gamma))
    if shape_token in {"gauss", "gaussian"}:
        sigma = fwhm_ppm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        return np.exp(-0.5 * (delta / sigma) ** 2) / (sigma * math.sqrt(2.0 * math.pi))
    if shape_token == "voigt":
        fraction = float(voigt_lorentz_fraction)
        if not 0 <= fraction <= 1:
            raise AppearanceError("voigt_lorentz_fraction must be between 0 and 1")
        if fraction <= 1e-12:
            return line_profile(axis_ppm, center_ppm, fwhm_ppm, "gaussian", fraction)
        if fraction >= 1 - 1e-12:
            return line_profile(axis_ppm, center_ppm, fwhm_ppm, "lorentzian", fraction)
        lorentz_fwhm = fraction * fwhm_ppm
        gaussian_fwhm = (1.0 - fraction) * fwhm_ppm
        gamma = lorentz_fwhm / 2.0
        sigma = gaussian_fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        z_value = (delta + 1j * gamma) / (sigma * math.sqrt(2.0))
        return np.real(wofz(z_value)) / (sigma * math.sqrt(2.0 * math.pi))
    raise AppearanceError(f"Unsupported lineshape {shape!r}; use lorentzian, gaussian or voigt")


def _add_component(
    axis: np.ndarray,
    signal: np.ndarray,
    center_ppm: float,
    area: float,
    linewidth_hz: float,
    spectrometer_mhz: float,
    shape: str,
    voigt_fraction: float,
    multiplicity: str = "s",
    j_values: Iterable[float] = (),
    explicit_pattern: list[dict[str, float]] | None = None,
) -> None:
    if explicit_pattern:
        offsets = np.asarray([float(item["offset_hz"]) for item in explicit_pattern])
        weights = np.asarray([float(item["weight"]) for item in explicit_pattern])
        if np.any(weights < 0) or weights.sum() <= 0:
            raise AppearanceError("A named-component pattern has invalid weights")
        weights = weights / weights.sum()
    else:
        offsets, weights = multiplet_components(multiplicity, j_values)
    width_ppm = linewidth_hz / spectrometer_mhz
    for offset_hz, weight in zip(offsets, weights):
        signal += area * float(weight) * line_profile(
            axis,
            center_ppm + float(offset_hz) / spectrometer_mhz,
            width_ppm,
            shape,
            voigt_fraction,
        )


def synthesize(lines: list[CleanLine], profile: dict[str, Any], seed: int) -> dict[str, np.ndarray]:
    nucleus = lines[0].nucleus
    axis_cfg = profile.get("axis", {})
    ppm_min = float(axis_cfg["ppm_min"])
    ppm_max = float(axis_cfg["ppm_max"])
    points = int(axis_cfg["points"])
    if not ppm_min < ppm_max or points < 512:
        raise AppearanceError("Profile axis requires ppm_min < ppm_max and at least 512 points")
    axis = np.linspace(ppm_min, ppm_max, points)

    spectrometer_mhz = float(profile["spectrometer_mhz"])
    if spectrometer_mhz <= 0:
        raise AppearanceError("spectrometer_mhz must be positive")
    shape = str(profile.get("lineshape", "voigt"))
    voigt_fraction = float(profile.get("voigt_lorentz_fraction", 0.65))
    default_width = float(profile["default_linewidth_hz"])
    exchange_extra = float(profile.get("exchangeable_H_extra_linewidth_hz", 0.0))

    analyte = np.zeros(points)
    for line in lines:
        linewidth = line.linewidth_hz if line.linewidth_hz is not None else default_width
        if nucleus == "1H" and line.exchangeable:
            linewidth += exchange_extra
        _add_component(
            axis,
            analyte,
            line.ppm,
            line.intensity,
            linewidth,
            spectrometer_mhz,
            shape,
            voigt_fraction,
            line.multiplicity,
            line.j_hz,
        )

    named = np.zeros(points)
    enabled_components: list[dict[str, Any]] = []
    for component in profile.get("named_components", []):
        if not parse_bool(component.get("enabled", True)):
            continue
        component_nucleus = normalize_nucleus(str(component.get("nucleus", nucleus)))
        if component_nucleus != nucleus:
            continue
        center = float(component["ppm"])
        if not ppm_min <= center <= ppm_max:
            raise AppearanceError(f"Named component {component.get('id')} lies outside the ppm axis")
        _add_component(
            axis,
            named,
            center,
            float(component["relative_integral"]),
            float(component["linewidth_hz"]),
            spectrometer_mhz,
            str(component.get("lineshape", shape)),
            float(component.get("voigt_lorentz_fraction", voigt_fraction)),
            str(component.get("multiplicity", "s")),
            parse_j_values(str(component.get("J_hz", ""))),
            component.get("pattern"),
        )
        enabled_components.append(component)

    signal_reference = max(float(np.max(np.abs(analyte + named))), np.finfo(float).eps)
    normalized_axis = 2.0 * (axis - ppm_min) / (ppm_max - ppm_min) - 1.0
    baseline_cfg = profile.get("baseline", {})
    coefficients = [float(value) for value in baseline_cfg.get("polynomial_coefficients", [0.0])]
    baseline_relative = np.polynomial.polynomial.polyval(normalized_axis, coefficients)
    roll = baseline_cfg.get("roll", {})
    baseline_relative += float(roll.get("amplitude", 0.0)) * np.sin(
        2.0 * math.pi * float(roll.get("cycles_across_axis", 1.0)) * (normalized_axis + 1.0) / 2.0
        + float(roll.get("phase_radians", 0.0))
    )
    baseline = signal_reference * baseline_relative
    deterministic = analyte + named + baseline

    noise_cfg = profile.get("noise", {})
    sigma = float(noise_cfg.get("sigma_relative_to_signal_max", 0.0)) * signal_reference
    rng = np.random.default_rng(int(seed))
    noise = rng.normal(0.0, sigma, points)
    correlation_points = float(noise_cfg.get("correlation_points", 0.0))
    if correlation_points > 0.5 and sigma > 0:
        radius = max(2, int(math.ceil(4 * correlation_points)))
        grid = np.arange(-radius, radius + 1, dtype=float)
        kernel = np.exp(-0.5 * (grid / correlation_points) ** 2)
        kernel /= kernel.sum()
        noise = np.convolve(noise, kernel, mode="same")
        rms = float(np.std(noise))
        if rms > 0:
            noise *= sigma / rms
    stochastic = deterministic + noise
    return {
        "ppm": axis,
        "analyte": analyte,
        "named_components": named,
        "baseline": baseline,
        "deterministic": deterministic,
        "noise": noise,
        "stochastic": stochastic,
        "enabled_components": np.asarray(enabled_components, dtype=object),
    }


def safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return token.strip("_.") or "spectrum"


def write_csv(path: Path, columns: list[tuple[str, np.ndarray]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow([name for name, _ in columns])
        for values in zip(*(array for _, array in columns)):
            writer.writerow([f"{float(value):.12g}" for value in values])


def plot_spectrum(
    path: Path,
    ppm: np.ndarray,
    intensity: np.ndarray,
    figure_id: str,
    nucleus: str,
    profile_id: str,
    layer: str,
    dpi: int,
) -> None:
    figure, axis = plt.subplots(figsize=(11.0, 5.2), constrained_layout=True)
    axis.plot(ppm, intensity, color="#222222", linewidth=0.85)
    axis.set_xlim(float(ppm.max()), float(ppm.min()))
    axis.set_xlabel(f"{nucleus} chemical shift (ppm)")
    axis.set_ylabel("arbitrary intensity")
    axis.set_title(f"{figure_id} | {layer} | {profile_id}", fontsize=11)
    axis.text(
        0.01,
        0.98,
        "SYNTHETIC APPEARANCE - NOT EXPERIMENTAL",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        color="#9c2f2f",
        fontweight="bold",
    )
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    figure.savefig(path, dpi=dpi, facecolor="white")
    plt.close(figure)


def render(
    input_path: Path,
    config_path: Path,
    output_dir: Path,
    profile_id: str | None = None,
    seed_override: int | None = None,
) -> dict[str, Path]:
    lines = read_clean_lines(input_path)
    payload = load_profiles(config_path)
    figure_id, nucleus = lines[0].figure_id, lines[0].nucleus
    selected_id, profile = select_profile(payload, profile_id, nucleus)
    seed = int(profile.get("noise", {}).get("seed", 20260810) if seed_override is None else seed_override)
    result = synthesize(lines, profile, seed)
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{safe_token(figure_id)}_{nucleus}"
    deterministic_csv = output_dir / f"{stem}_deterministic.csv"
    stochastic_csv = output_dir / f"{stem}_stochastic.csv"
    deterministic_png = output_dir / f"{stem}_deterministic.png"
    stochastic_png = output_dir / f"{stem}_stochastic.png"
    metadata_path = output_dir / f"{stem}_appearance_metadata.json"

    write_csv(
        deterministic_csv,
        [
            ("ppm", result["ppm"]),
            ("analyte_signal", result["analyte"]),
            ("named_solvent_water_signal", result["named_components"]),
            ("baseline_signal", result["baseline"]),
            ("deterministic_total", result["deterministic"]),
        ],
    )
    write_csv(
        stochastic_csv,
        [
            ("ppm", result["ppm"]),
            ("deterministic_total", result["deterministic"]),
            ("seeded_noise", result["noise"]),
            ("stochastic_total", result["stochastic"]),
        ],
    )
    dpi = int(profile.get("plot_dpi", 180))
    plot_spectrum(
        deterministic_png,
        result["ppm"],
        result["deterministic"],
        figure_id,
        nucleus,
        selected_id,
        "deterministic conditioned layer",
        dpi,
    )
    plot_spectrum(
        stochastic_png,
        result["ppm"],
        result["stochastic"],
        figure_id,
        nucleus,
        selected_id,
        "seeded stochastic conditioned layer",
        dpi,
    )
    metadata = {
        "schema_version": "1.0",
        "figure_id": figure_id,
        "nucleus": nucleus,
        "profile_id": selected_id,
        "profile_status": profile["status"],
        "seed": seed,
        "disclaimer": DISCLAIMER,
        "scientific_scope": {
            "calls_orca": False,
            "first_principles_prediction": False,
            "experimental_fit_performed": False,
            "role": "conditional instrument/experimental appearance model",
        },
        "input": {"file": input_path.name, "sha256": sha256(input_path), "line_count": len(lines)},
        "configuration": {"file": config_path.name, "sha256": sha256(config_path)},
        "lineshape": profile.get("lineshape", "voigt"),
        "spectrometer_mhz": profile["spectrometer_mhz"],
        "named_components": [
            {
                "id": item.get("id"),
                "identity": item.get("identity"),
                "ppm": item.get("ppm"),
                "status": item.get("status", profile["status"]),
            }
            for item in result["enabled_components"].tolist()
        ],
        "outputs": {
            "deterministic_csv": deterministic_csv.name,
            "stochastic_csv": stochastic_csv.name,
            "deterministic_png": deterministic_png.name,
            "stochastic_png": stochastic_png.name,
        },
    }
    metadata_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return {
        "deterministic_csv": deterministic_csv,
        "stochastic_csv": stochastic_csv,
        "deterministic_png": deterministic_png,
        "stochastic_png": stochastic_png,
        "metadata": metadata_path,
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Render declared deterministic and seeded-stochastic NMR appearance layers."
    )
    result.add_argument("--input", type=Path, help="Clean-line CSV")
    result.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    result.add_argument("--profile", help="Profile ID; otherwise use the nucleus default")
    result.add_argument("--output-dir", type=Path, help="Destination directory")
    result.add_argument("--seed", type=int, help="Override the profile noise seed")
    result.add_argument("--list-profiles", action="store_true")
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        payload = load_profiles(args.config)
        if args.list_profiles:
            for profile_name, profile in payload["profiles"].items():
                print(f"{profile_name}\t{profile.get('nucleus')}\t{profile.get('status')}")
            return 0
        if args.input is None or args.output_dir is None:
            raise AppearanceError("--input and --output-dir are required unless --list-profiles is used")
        outputs = render(args.input, args.config, args.output_dir, args.profile, args.seed)
    except (AppearanceError, OSError, KeyError, ValueError) as exc:
        print(f"Appearance rendering failed: {exc}", file=sys.stderr)
        return 2
    print(DISCLAIMER)
    for name, path in outputs.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
