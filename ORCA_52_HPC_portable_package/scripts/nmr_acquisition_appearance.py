from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import uuid
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_CALIBRATION_DIR = PROJECT / "data" / "real_nmr_calibration_audit"
DEFAULT_OUTPUT_ROOT = PROJECT / "results" / "NMR_acquisition_appearance_SYNTHETIC_NOT_ACQUIRED"

SCHEMA_VERSION = "1.0.0"
LAYER_A_NOTICE = "CALCULATED MODEL - NOT EXPERIMENTAL"
LAYER_B_NOTICE = "SYNTHETIC APPEARANCE - NOT ACQUIRED"
GLOBAL_WARNING = (
    "These spectra are calculations/synthetic visualizations. They are not acquired data "
    "and cannot establish identity, purity, yield, or experimental reproducibility."
)
SUPPORTED_NUCLEI = {"1H", "13C", "19F", "31P"}
REQUIRED_STICK_FIELDS = {
    "transition_id",
    "component_id",
    "nucleus",
    "shift_ppm",
    "intensity",
    "linewidth_hz",
    "assignment",
    "data_status",
}

PALETTE = {
    "calculated": "#0F4D92",
    "synthetic": "#272727",
    "solvent": "#42949E",
    "water": "#9A4D8E",
    "neutral": "#767676",
    "light": "#D8D8D8",
    "warning": "#B64342",
}


class AppearanceError(RuntimeError):
    """Raised when a synthetic spectrum would be scientifically ambiguous."""


@dataclass(frozen=True)
class Transition:
    transition_id: str
    component_id: str
    nucleus: str
    shift_ppm: float
    intensity: float
    linewidth_hz: float
    assignment: str
    data_status: str


@dataclass
class AcquisitionProfile:
    solvent: str
    nucleus: str
    calibration_status: str
    calibration_group: str
    calibration_dataset_id: str
    frequency_mhz: float
    spectral_width_hz: float
    acquired_complex_points: int
    processed_points: int
    acquisition_time_s: float
    line_broadening_hz: float
    high_ppm: float
    low_ppm: float
    temperature_k: float | None
    scans: int | None
    delay_s: float | None
    pulse_us: float | None
    noise_stationary_rel_max_percent: float
    noise_target_lag1: float
    baseline_span_rel_max_percent: float
    baseline_slow_rms_rel_max_percent: float
    parameter_sources: dict[str, str]


@dataclass(frozen=True)
class NamedComponent:
    component_id: str
    component_name: str
    component_class: str
    center_ppm: float | None
    fwhm_hz: float | None
    target_height_rel_analyte_percent: float | None
    provenance_status: str
    calibration_record: str
    generation_status: str
    multiplet_offsets_hz: tuple[float, ...] = (0.0,)
    multiplet_intensities: tuple[float, ...] = (1.0,)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def derive_seed(*parts: Any) -> int:
    payload = "\x1f".join(str(part) for part in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big", signed=False)


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{uuid.uuid4().hex}")
    try:
        temporary.write_text(text, encoding="utf-8", newline="\n")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_write_json(path: Path, payload: Any) -> None:
    atomic_write_text(
        path,
        json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
    )


def atomic_write_csv(path: Path, rows: Iterable[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{uuid.uuid4().hex}")
    try:
        with temporary.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="raise")
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{uuid.uuid4().hex}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(fields), delimiter="\t", extrasaction="raise"
            )
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _float(value: Any, field: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise AppearanceError(f"Invalid numeric value for {field}: {value!r}") from exc
    if not math.isfinite(result):
        raise AppearanceError(f"Non-finite numeric value for {field}: {value!r}")
    return result


def _optional_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def _optional_int(value: Any) -> int | None:
    parsed = _optional_float(value)
    return None if parsed is None else int(round(parsed))


def normalize_solvent(value: str) -> str:
    token = value.strip().upper().replace("₆", "6").replace("₂", "2").replace("-", "")
    if token in {"DMSO", "DMSOD6", "D6DMSO"}:
        return "DMSO-d6"
    if token in {"D2O", "2HO", "DEUTERIUMOXIDE"}:
        return "D2O"
    raise AppearanceError("Only DMSO-d6 and D2O appearance presets are supported")


def load_transitions(path: Path, nucleus: str) -> list[Transition]:
    if nucleus not in SUPPORTED_NUCLEI:
        raise AppearanceError(f"Unsupported nucleus {nucleus!r}; choose {sorted(SUPPORTED_NUCLEI)}")
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_STICK_FIELDS - set(reader.fieldnames or [])
        if missing:
            raise AppearanceError(f"Stick CSV lacks required columns: {sorted(missing)}")
        rows = list(reader)
    selected: list[Transition] = []
    seen: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        row_nucleus = str(row["nucleus"]).strip()
        if row_nucleus != nucleus:
            continue
        transition_id = str(row["transition_id"]).strip()
        component_id = str(row["component_id"]).strip()
        assignment = str(row["assignment"]).strip()
        status = str(row["data_status"]).strip().upper()
        if not transition_id or not component_id or not assignment:
            raise AppearanceError(
                f"Row {row_number}: transition_id, component_id and assignment must be named"
            )
        if transition_id in seen:
            raise AppearanceError(f"Duplicate transition_id for {nucleus}: {transition_id}")
        if status not in {"CALCULATED", "CALCULATED_PREDICTION"}:
            raise AppearanceError(
                f"Row {row_number}: data_status must be CALCULATED or CALCULATED_PREDICTION; "
                "experimental and anonymous impurity rows are not accepted by this module"
            )
        transition = Transition(
            transition_id=transition_id,
            component_id=component_id,
            nucleus=row_nucleus,
            shift_ppm=_float(row["shift_ppm"], "shift_ppm"),
            intensity=_float(row["intensity"], "intensity"),
            linewidth_hz=_float(row["linewidth_hz"], "linewidth_hz"),
            assignment=assignment,
            data_status="CALCULATED",
        )
        if transition.intensity <= 0.0:
            raise AppearanceError(f"Row {row_number}: intensity must be positive")
        if transition.linewidth_hz <= 0.0:
            raise AppearanceError(f"Row {row_number}: linewidth_hz must be positive")
        selected.append(transition)
        seen.add(transition_id)
    if not selected:
        raise AppearanceError(f"No {nucleus} transitions found in {path}")
    return sorted(selected, key=lambda row: (row.shift_ppm, row.transition_id))


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def calibration_hashes(calibration_dir: Path) -> list[dict[str, Any]]:
    names = [
        "nmr_acquisition_processing_audit.csv",
        "nmr_spectral_calibration_per_dataset.csv",
        "nmr_representative_calibration_by_group.csv",
    ]
    return [
        {"file": name, "sha256": sha256_file(calibration_dir / name)} for name in names
    ]


def _select_empirical_row(rows: Sequence[dict[str, str]], seed: int) -> dict[str, str]:
    if not rows:
        raise AppearanceError("No empirical calibration rows satisfy the requested group")
    ordered = sorted(rows, key=lambda row: str(row.get("dataset_id", "")))
    index = int(np.random.default_rng(seed).integers(0, len(ordered)))
    return ordered[index]


def _preset_profile(solvent: str, nucleus: str) -> AcquisitionProfile:
    presets: dict[str, dict[str, float | int]] = {
        "1H": {
            "frequency_mhz": 400.13,
            "spectral_width_hz": 8012.82,
            "acquired_complex_points": 32768,
            "processed_points": 65536,
            "acquisition_time_s": 4.0894,
            "line_broadening_hz": 0.30,
            "high_ppm": 16.18,
            "low_ppm": -3.845,
            "noise_percent": 0.00269,
            "noise_lag1": 0.766,
            "baseline_span_percent": 0.217,
            "baseline_rms_percent": 0.117,
            "scans": 8,
            "delay_s": 1.0,
            "pulse_us": 10.0,
        },
        "13C": {
            "frequency_mhz": 100.613,
            "spectral_width_hz": 23809.52,
            "acquired_complex_points": 32768,
            "processed_points": 32768,
            "acquisition_time_s": 1.3763,
            "line_broadening_hz": 1.00,
            "high_ppm": 218.0,
            "low_ppm": -18.7,
            "noise_percent": 0.109,
            "noise_lag1": 0.619,
            "baseline_span_percent": 0.0364,
            "baseline_rms_percent": 0.0122,
            "scans": 512,
            "delay_s": 2.0,
            "pulse_us": 8.0,
        },
        "19F": {
            "frequency_mhz": 376.50,
            "spectral_width_hz": 75300.0,
            "acquired_complex_points": 32768,
            "processed_points": 65536,
            "acquisition_time_s": 0.435,
            "line_broadening_hz": 1.00,
            "high_ppm": 50.0,
            "low_ppm": -150.0,
            "noise_percent": 0.015,
            "noise_lag1": 0.65,
            "baseline_span_percent": 0.08,
            "baseline_rms_percent": 0.03,
            "scans": 32,
            "delay_s": 1.0,
            "pulse_us": 12.0,
        },
        "31P": {
            "frequency_mhz": 161.98,
            "spectral_width_hz": 32396.0,
            "acquired_complex_points": 32768,
            "processed_points": 65536,
            "acquisition_time_s": 1.0115,
            "line_broadening_hz": 1.00,
            "high_ppm": 100.0,
            "low_ppm": -100.0,
            "noise_percent": 0.020,
            "noise_lag1": 0.65,
            "baseline_span_percent": 0.08,
            "baseline_rms_percent": 0.03,
            "scans": 64,
            "delay_s": 1.0,
            "pulse_us": 10.0,
        },
    }
    p = presets[nucleus]
    status = (
        "ASSUMED_NOT_LOCALLY_CALIBRATED"
        if solvent == "D2O"
        else "ASSUMED_NUCLEUS_PRESET_NOT_LOCALLY_CALIBRATED"
    )
    source = "D2O appearance preset; no local D2O calibration data"
    if solvent != "D2O":
        source = f"{nucleus} appearance preset; local DMSO audit has no {nucleus} group"
    return AcquisitionProfile(
        solvent=solvent,
        nucleus=nucleus,
        calibration_status=status,
        calibration_group=f"{solvent}_{nucleus}",
        calibration_dataset_id="NONE",
        frequency_mhz=float(p["frequency_mhz"]),
        spectral_width_hz=float(p["spectral_width_hz"]),
        acquired_complex_points=int(p["acquired_complex_points"]),
        processed_points=int(p["processed_points"]),
        acquisition_time_s=float(p["acquisition_time_s"]),
        line_broadening_hz=float(p["line_broadening_hz"]),
        high_ppm=float(p["high_ppm"]),
        low_ppm=float(p["low_ppm"]),
        temperature_k=298.15,
        scans=int(p["scans"]),
        delay_s=float(p["delay_s"]),
        pulse_us=float(p["pulse_us"]),
        noise_stationary_rel_max_percent=float(p["noise_percent"]),
        noise_target_lag1=float(p["noise_lag1"]),
        baseline_span_rel_max_percent=float(p["baseline_span_percent"]),
        baseline_slow_rms_rel_max_percent=float(p["baseline_rms_percent"]),
        parameter_sources={"all": source},
    )


def load_acquisition_profile(
    calibration_dir: Path,
    solvent: str,
    nucleus: str,
    seed: int,
    acquired_points_override: int | None = None,
    processed_points_override: int | None = None,
) -> AcquisitionProfile:
    if solvent != "DMSO-d6" or nucleus not in {"1H", "13C"}:
        profile = _preset_profile(solvent, nucleus)
    else:
        acquisition_rows = read_csv(calibration_dir / "nmr_acquisition_processing_audit.csv")
        spectral_rows = read_csv(calibration_dir / "nmr_spectral_calibration_per_dataset.csv")
        candidates = [
            row
            for row in acquisition_rows
            if str(row.get("primary_dataset_for_filename", "")).lower() == "true"
            and str(row.get("solvent_parameter", "")).upper() == "DMSO"
            and str(row.get("nucleus", "")) == nucleus
        ]
        selected = _select_empirical_row(candidates, seed)
        spectral_by_id = {row["dataset_id"]: row for row in spectral_rows}
        spectral = spectral_by_id.get(selected["dataset_id"])
        if spectral is None:
            raise AppearanceError(
                f"Calibration audit lacks spectral row for {selected['dataset_id']}"
            )
        frequency = _float(selected["frequency_MHz"], "frequency_MHz")
        sw = _float(selected["spectral_width_Hz"], "spectral_width_Hz")
        acquired = int(_float(selected["acquired_complex_points"], "acquired_complex_points"))
        processed = int(_float(selected["SI_processed_points"], "SI_processed_points"))
        high = _float(selected["high_ppm_document_axis"], "high_ppm_document_axis")
        low = _float(selected["low_ppm_nominal_document_axis"], "low_ppm_nominal_document_axis")
        profile = AcquisitionProfile(
            solvent=solvent,
            nucleus=nucleus,
            calibration_status="LOCAL_EMPIRICAL_DMSO_DISTRIBUTION",
            calibration_group=f"DMSO_{nucleus}",
            calibration_dataset_id=selected["dataset_id"],
            frequency_mhz=frequency,
            spectral_width_hz=sw,
            acquired_complex_points=acquired,
            processed_points=processed,
            acquisition_time_s=_float(selected["AQ_s"], "AQ_s"),
            line_broadening_hz=_float(selected["LB_Hz"], "LB_Hz"),
            high_ppm=high,
            low_ppm=low,
            temperature_k=_optional_float(selected.get("temperature_K")),
            scans=_optional_int(selected.get("NS")),
            delay_s=_optional_float(selected.get("D1_s")),
            pulse_us=_optional_float(selected.get("P1_us")),
            noise_stationary_rel_max_percent=_float(
                spectral["noise_stationary_rel_max_percent"],
                "noise_stationary_rel_max_percent",
            ),
            noise_target_lag1=_float(
                spectral["noise_lag1_correlation_clipped"],
                "noise_lag1_correlation_clipped",
            ),
            baseline_span_rel_max_percent=_float(
                spectral["baseline_span_rel_max_percent"],
                "baseline_span_rel_max_percent",
            ),
            baseline_slow_rms_rel_max_percent=_float(
                spectral["baseline_slow_rms_rel_max_percent"],
                "baseline_slow_rms_rel_max_percent",
            ),
            parameter_sources={
                "acquisition": "local empirical DMSO audit row",
                "noise_and_baseline": "matched local processed-spectrum audit row",
                "selection": "deterministic empirical draw from eligible primary datasets",
            },
        )
    if acquired_points_override is not None:
        if acquired_points_override < 256:
            raise AppearanceError("acquired-points override must be at least 256")
        profile.acquired_complex_points = int(acquired_points_override)
        profile.acquisition_time_s = profile.acquired_complex_points / profile.spectral_width_hz
        profile.parameter_sources["acquired_complex_points"] = "explicit command-line override"
    if processed_points_override is not None:
        if processed_points_override < profile.acquired_complex_points:
            raise AppearanceError("processed-points cannot be smaller than acquired-points")
        profile.processed_points = int(processed_points_override)
        profile.parameter_sources["processed_points"] = "explicit command-line override"
    if profile.processed_points < profile.acquired_complex_points:
        raise AppearanceError("Calibration profile has processed_points < acquired_complex_points")
    return profile


def _sample_named_peak(
    calibration_dir: Path,
    kind: str,
    seed: int,
) -> tuple[float, float, float, str, str]:
    rows = read_csv(calibration_dir / "nmr_spectral_calibration_per_dataset.csv")
    prefix = "residual_dmso" if kind == "residual_dmso" else "water_in_dmso"
    candidates = [
        row
        for row in rows
        if str(row.get("primary_dataset_for_filename", "")).lower() == "true"
        and str(row.get("solvent_parameter", "")).upper() == "DMSO"
        and str(row.get("nucleus", "")) == "1H"
        and str(row.get(f"{prefix}_confidence", "")).lower() == "high"
        and row.get(f"{prefix}_center_ppm") not in (None, "")
        and row.get(f"{prefix}_fwhm_hz") not in (None, "")
        and row.get(f"{prefix}_height_rel_tallest_peak_percent") not in (None, "")
    ]
    selected = _select_empirical_row(candidates, seed)
    return (
        _float(selected[f"{prefix}_center_ppm"], f"{prefix}_center_ppm"),
        _float(selected[f"{prefix}_fwhm_hz"], f"{prefix}_fwhm_hz"),
        _float(
            selected[f"{prefix}_height_rel_tallest_peak_percent"],
            f"{prefix}_height_rel_tallest_peak_percent",
        ),
        selected["dataset_id"],
        "LOCAL_EMPIRICAL_DMSO_DISTRIBUTION",
    )


def named_components(
    calibration_dir: Path,
    solvent: str,
    nucleus: str,
    spectrum_seed: int,
) -> list[NamedComponent]:
    if solvent == "DMSO-d6" and nucleus == "1H":
        dmso = _sample_named_peak(
            calibration_dir,
            "residual_dmso",
            derive_seed(spectrum_seed, "residual_dmso_d6_1h"),
        )
        water = _sample_named_peak(
            calibration_dir,
            "water_in_dmso",
            derive_seed(spectrum_seed, "water_in_dmso_d6_1h"),
        )
        return [
            NamedComponent(
                "residual_dmso_d6_1h",
                "residual DMSO-d6 proton",
                "named_residual_solvent",
                dmso[0],
                dmso[1],
                dmso[2],
                dmso[4],
                dmso[3],
                "generated",
            ),
            NamedComponent(
                "water_in_dmso_d6_1h",
                "water in DMSO-d6",
                "named_water",
                water[0],
                water[1],
                water[2],
                water[4],
                water[3],
                "generated",
            ),
        ]
    if solvent == "DMSO-d6" and nucleus == "13C":
        intensities = (1.0, 3.0, 6.0, 7.0, 6.0, 3.0, 1.0)
        offsets = tuple((index - 3) * 19.8 for index in range(7))
        return [
            NamedComponent(
                "dmso_d6_13c_septet",
                "DMSO-d6 carbon (1:3:6:7:6:3:1 septet)",
                "named_residual_solvent",
                39.52,
                1.5,
                35.0,
                "ASSUMED_STANDARD_SOLVENT_COMPONENT_NOT_LOCALLY_PEAK_CALIBRATED",
                "no local 13C solvent-peak fit; conventional center and 1J(C,D) used",
                "generated",
                offsets,
                intensities,
            ),
            NamedComponent(
                "water_13c_not_applicable",
                "water has no ordinary 13C resonance",
                "named_water",
                None,
                None,
                None,
                "NOT_APPLICABLE",
                "nucleus-specific exclusion",
                "not_generated",
            ),
        ]
    if solvent == "D2O" and nucleus == "1H":
        return [
            NamedComponent(
                "hdo_in_d2o_1h",
                "HDO in D2O",
                "named_water",
                4.79,
                8.0,
                80.0,
                "ASSUMED_NOT_LOCALLY_CALIBRATED",
                "D2O/HDO preset; position and width depend on temperature, pD and exchange",
                "generated",
            )
        ]
    reason = (
        "D2O has no ordinary solvent resonance for this observed nucleus"
        if solvent == "D2O"
        else "DMSO-d6 solvent/water appearance component not defined for this nucleus"
    )
    return [
        NamedComponent(
            f"solvent_component_{nucleus}_not_generated",
            reason,
            "named_solvent_or_water",
            None,
            None,
            None,
            "NOT_APPLICABLE",
            "nucleus-specific exclusion",
            "not_generated",
        )
    ]


def fid_from_lines(
    shifts_ppm: np.ndarray,
    amplitudes: np.ndarray,
    linewidths_hz: np.ndarray,
    profile: AcquisitionProfile,
) -> np.ndarray:
    n = profile.acquired_complex_points
    time = np.arange(n, dtype=float) / profile.spectral_width_hz
    carrier_ppm = (profile.high_ppm + profile.low_ppm) / 2.0
    result = np.zeros(n, dtype=np.complex128)
    for shift, amplitude, linewidth in zip(shifts_ppm, amplitudes, linewidths_hz):
        offset_hz = (float(shift) - carrier_ppm) * profile.frequency_mhz
        result += float(amplitude) * np.exp(
            2j * np.pi * offset_hz * time - np.pi * float(linewidth) * time
        )
    return result


def fft_spectrum(fid: np.ndarray, profile: AcquisitionProfile, window: np.ndarray | None) -> np.ndarray:
    processed = fid if window is None else fid * window
    return np.fft.fftshift(np.fft.fft(processed, n=profile.processed_points))


def ppm_axis(profile: AcquisitionProfile) -> np.ndarray:
    dwell = 1.0 / profile.spectral_width_hz
    frequency_hz = np.fft.fftshift(np.fft.fftfreq(profile.processed_points, d=dwell))
    carrier_ppm = (profile.high_ppm + profile.low_ppm) / 2.0
    return carrier_ppm + frequency_hz / profile.frequency_mhz


def exponential_window(profile: AcquisitionProfile) -> np.ndarray:
    time = np.arange(profile.acquired_complex_points, dtype=float) / profile.spectral_width_hz
    return np.exp(-np.pi * profile.line_broadening_hz * time)


def apply_phase(
    spectrum: np.ndarray, phase0_deg: float, phase1_deg: float
) -> np.ndarray:
    normalized = np.linspace(-1.0, 1.0, spectrum.size, endpoint=True)
    radians = np.deg2rad(phase0_deg + phase1_deg * normalized)
    return spectrum * np.exp(1j * radians)


def _complex_ar1_from_innovations(
    innovations: np.ndarray, coefficient: float
) -> np.ndarray:
    """Return a circular complex Gaussian AR(1) realization.

    Reusing one innovation vector while solving for the coefficient makes the
    lag-matching step deterministic and monotonic enough for a bounded search.
    """
    coefficient = float(np.clip(coefficient, -0.98, 0.98))
    result = np.empty(innovations.size, dtype=np.complex128)
    result[0] = innovations[0]
    innovation_scale = math.sqrt(1.0 - coefficient**2)
    for index in range(1, innovations.size):
        result[index] = (
            coefficient * result[index - 1]
            + innovation_scale * innovations[index]
        )
    return result


def calibrated_correlated_noise_fid(
    profile: AcquisitionProfile,
    window: np.ndarray,
    target_frequency_lag1: float,
    seed: int,
    phase0_deg: float,
    phase1_deg: float,
) -> tuple[np.ndarray, float, float]:
    """Construct a complex Gaussian FID whose processed noise matches lag-1.

    The local audit reports correlation after truncation, apodization, zero
    filling and Fourier transformation.  We therefore create a deterministic
    circular complex Gaussian AR(1) realization on the processed frequency grid,
    project it to a complex time-domain FID, truncate it to the acquired TD, and
    run the same processing chain.  A bounded coefficient search matches the
    *processed* lag-1 rather than pretending it is a raw-FID parameter.
    """
    rng = np.random.default_rng(seed)
    innovations = (
        rng.normal(size=profile.processed_points)
        + 1j * rng.normal(size=profile.processed_points)
    ) / math.sqrt(2.0)

    def project(coefficient: float) -> tuple[np.ndarray, float]:
        frequency_noise = _complex_ar1_from_innovations(innovations, coefficient)
        full_time = np.fft.ifft(np.fft.ifftshift(frequency_noise))
        raw_fid = full_time[: profile.acquired_complex_points].copy()
        processed = apply_phase(
            fft_spectrum(raw_fid, profile, window), phase0_deg, phase1_deg
        )
        return raw_fid, lag1_correlation(np.real(processed))

    target = float(np.clip(target_frequency_lag1, -0.95, 0.95))
    low, high = -0.98, 0.98
    low_fid, low_value = project(low)
    high_fid, high_value = project(high)
    if target <= low_value:
        return low_fid, low, low_value
    if target >= high_value:
        return high_fid, high, high_value
    selected_fid = low_fid
    selected_coefficient = low
    selected_value = low_value
    for _ in range(18):
        middle = (low + high) / 2.0
        middle_fid, middle_value = project(middle)
        if abs(middle_value - target) < abs(selected_value - target):
            selected_fid = middle_fid
            selected_coefficient = middle
            selected_value = middle_value
        if middle_value < target:
            low = middle
        else:
            high = middle
    return selected_fid, selected_coefficient, selected_value


def lag1_correlation(values: np.ndarray) -> float:
    centered = np.asarray(values, dtype=float) - float(np.mean(values))
    if centered.size < 3 or np.std(centered) == 0.0:
        return 0.0
    return float(np.corrcoef(centered[:-1], centered[1:])[0, 1])


def slow_baseline(size: int, target_span: float, seed: int) -> np.ndarray:
    if target_span <= 0.0:
        return np.zeros(size, dtype=float)
    rng = np.random.default_rng(seed)
    x = np.linspace(-1.0, 1.0, size)
    coefficients = rng.normal(0.0, 1.0, size=4)
    baseline = (
        coefficients[0] * x
        + coefficients[1] * (x**2 - np.mean(x**2))
        + coefficients[2] * (x**3 - np.mean(x**3))
        + 0.5 * coefficients[3] * np.sin(2.0 * np.pi * (0.65 * x + 0.17))
    )
    baseline -= np.mean(baseline)
    span = float(np.percentile(baseline, 95) - np.percentile(baseline, 5))
    if span == 0.0:
        return np.zeros(size, dtype=float)
    return baseline * (target_span / span)


def component_fid(component: NamedComponent, profile: AcquisitionProfile) -> np.ndarray:
    if component.generation_status != "generated":
        return np.zeros(profile.acquired_complex_points, dtype=np.complex128)
    assert component.center_ppm is not None
    assert component.fwhm_hz is not None
    offsets = np.asarray(component.multiplet_offsets_hz, dtype=float)
    intensities = np.asarray(component.multiplet_intensities, dtype=float)
    intensities /= np.sum(intensities)
    shifts = component.center_ppm + offsets / profile.frequency_mhz
    linewidths = np.full_like(shifts, component.fwhm_hz, dtype=float)
    return fid_from_lines(shifts, intensities, linewidths, profile)


def _peak_height(spectrum: np.ndarray) -> float:
    value = float(np.max(np.abs(np.real(spectrum))))
    if not math.isfinite(value) or value <= 0.0:
        raise AppearanceError("Spectrum has no finite non-zero absorption signal")
    return value


def build_appearance(
    transitions: Sequence[Transition],
    profile: AcquisitionProfile,
    components: Sequence[NamedComponent],
    spectrum_seed: int,
    phase0_deg: float | None,
    phase1_deg: float | None,
) -> dict[str, Any]:
    shifts = np.asarray([row.shift_ppm for row in transitions], dtype=float)
    intensities = np.asarray([row.intensity for row in transitions], dtype=float)
    linewidths = np.asarray([row.linewidth_hz for row in transitions], dtype=float)
    if np.min(shifts) < profile.low_ppm or np.max(shifts) > profile.high_ppm:
        raise AppearanceError(
            f"Calculated transition range {np.min(shifts):.3f} to {np.max(shifts):.3f} ppm "
            f"falls outside acquisition range {profile.low_ppm:.3f} to {profile.high_ppm:.3f} ppm"
        )

    analyte_fid = fid_from_lines(shifts, intensities, linewidths, profile)
    layer_a_complex = fft_spectrum(analyte_fid, profile, window=None)
    layer_a_scale = _peak_height(layer_a_complex)
    layer_a = np.real(layer_a_complex) / layer_a_scale

    window = exponential_window(profile)
    analyte_processed = fft_spectrum(analyte_fid, profile, window)
    analyte_peak = _peak_height(analyte_processed)
    named_fid = np.zeros_like(analyte_fid)
    component_scales: dict[str, float] = {}
    for component in components:
        if component.generation_status != "generated":
            component_scales[component.component_id] = 0.0
            continue
        unit_fid = component_fid(component, profile)
        unit_spectrum = fft_spectrum(unit_fid, profile, window)
        unit_peak = _peak_height(unit_spectrum)
        target_fraction = float(component.target_height_rel_analyte_percent or 0.0) / 100.0
        scale = target_fraction * analyte_peak / unit_peak
        named_fid += unit_fid * scale
        component_scales[component.component_id] = scale

    phase_rng = np.random.default_rng(derive_seed(spectrum_seed, "phase"))
    actual_phase0 = float(phase0_deg) if phase0_deg is not None else float(phase_rng.normal(0, 0.30))
    actual_phase1 = float(phase1_deg) if phase1_deg is not None else float(phase_rng.normal(0, 0.12))

    signal_fid = analyte_fid + named_fid
    signal_processed_unphased = fft_spectrum(signal_fid, profile, window)
    signal_processed = apply_phase(signal_processed_unphased, actual_phase0, actual_phase1)
    analyte_processed_phased = apply_phase(analyte_processed, actual_phase0, actual_phase1)
    named_processed_phased = apply_phase(
        fft_spectrum(named_fid, profile, window), actual_phase0, actual_phase1
    )
    signal_peak = _peak_height(signal_processed)

    target_lag1 = float(np.clip(profile.noise_target_lag1, -0.95, 0.95))
    unit_noise_fid, spectral_ar1, unscaled_achieved_lag1 = calibrated_correlated_noise_fid(
        profile,
        window,
        target_lag1,
        derive_seed(spectrum_seed, "complex_fid_noise"),
        actual_phase0,
        actual_phase1,
    )
    unit_noise_processed = apply_phase(
        fft_spectrum(unit_noise_fid, profile, window), actual_phase0, actual_phase1
    )
    unit_sigma = float(np.std(np.real(unit_noise_processed)))
    target_sigma = profile.noise_stationary_rel_max_percent / 100.0 * signal_peak
    noise_scale = 0.0 if unit_sigma == 0.0 else target_sigma / unit_sigma
    noise_fid = unit_noise_fid * noise_scale
    noise_processed = unit_noise_processed * noise_scale

    baseline_target_span = profile.baseline_span_rel_max_percent / 100.0 * signal_peak
    baseline = slow_baseline(
        profile.processed_points,
        baseline_target_span,
        derive_seed(spectrum_seed, "slow_baseline"),
    )
    final_complex = signal_processed + noise_processed
    final_real = np.real(final_complex) + baseline
    normalization = max(float(np.max(np.abs(final_real))), np.finfo(float).tiny)
    axis = ppm_axis(profile)
    total_fid = signal_fid + noise_fid

    return {
        "axis_ppm": axis,
        "layer_a": layer_a,
        "layer_a_complex": layer_a_complex / layer_a_scale,
        "analyte_fid": analyte_fid,
        "named_fid": named_fid,
        "noise_fid": noise_fid,
        "total_fid": total_fid,
        "window": window,
        "analyte_processed": analyte_processed_phased,
        "named_processed": named_processed_phased,
        "noise_processed": noise_processed,
        "baseline": baseline,
        "final_real": final_real,
        "final_imag": np.imag(final_complex),
        "final_normalized": final_real / normalization,
        "normalization": normalization,
        "component_scales": component_scales,
        "phase0_deg": actual_phase0,
        "phase1_deg": actual_phase1,
        "noise_scale": noise_scale,
        "noise_target_sigma": target_sigma,
        "noise_achieved_sigma": float(np.std(np.real(noise_processed))),
        "noise_achieved_lag1": lag1_correlation(np.real(noise_processed)),
        "baseline_achieved_span": float(np.percentile(baseline, 95) - np.percentile(baseline, 5)),
        "baseline_achieved_rms": float(np.std(baseline)),
        "noise_spectral_ar1_coefficient_before_time_projection": spectral_ar1,
        "noise_unscaled_achieved_lag1": unscaled_achieved_lag1,
    }


def configure_matplotlib() -> str:
    try:
        arial = font_manager.findfont("Arial", fallback_to_default=False)
    except ValueError as exc:
        raise AppearanceError("Arial is required for NMR figure export") from exc
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )
    return arial


def nucleus_math(nucleus: str) -> str:
    return {
        "1H": r"$^{1}\mathrm{H}$",
        "13C": r"$^{13}\mathrm{C}$",
        "19F": r"$^{19}\mathrm{F}$",
        "31P": r"$^{31}\mathrm{P}$",
    }[nucleus]


def default_display_limits(
    nucleus: str,
    profile: AcquisitionProfile,
    transitions: Sequence[Transition],
    high_override: float | None,
    low_override: float | None,
) -> tuple[float, float]:
    defaults = {
        "1H": (12.0, -0.5),
        "13C": (218.0, -10.0),
        "19F": (30.0, -140.0),
        "31P": (80.0, -80.0),
    }
    high, low = defaults[nucleus]
    max_shift = max(row.shift_ppm for row in transitions)
    min_shift = min(row.shift_ppm for row in transitions)
    high = max(high, max_shift + 0.5)
    low = min(low, min_shift - 0.5)
    high = min(profile.high_ppm, high_override if high_override is not None else high)
    low = max(profile.low_ppm, low_override if low_override is not None else low)
    if high <= low:
        raise AppearanceError("Display ppm high limit must exceed low limit")
    return high, low


def render_figure(
    output_stem: Path,
    spectrum_id: str,
    display_name: str,
    nucleus: str,
    profile: AcquisitionProfile,
    transitions: Sequence[Transition],
    components: Sequence[NamedComponent],
    result: dict[str, Any],
    display_high: float,
    display_low: float,
    dpi: int,
) -> list[Path]:
    configure_matplotlib()
    fig, axes = plt.subplots(2, 1, figsize=(7.09, 3.55), sharex=True)
    axis = result["axis_ppm"]
    axes[0].plot(axis, result["layer_a"], color=PALETTE["calculated"], lw=0.75)
    axes[1].plot(axis, result["final_normalized"], color=PALETTE["synthetic"], lw=0.62)
    for index, ax in enumerate(axes):
        ax.set_xlim(display_high, display_low)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="x", direction="out", length=3.0, width=0.7, labelsize=7)
        ax.text(
            -0.035,
            1.01,
            "a" if index == 0 else "b",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            fontweight="bold",
        )
    axes[0].set_title(
        f"{nucleus_math(nucleus)} NMR | calculated analyte model",
        loc="left",
        fontsize=8,
        fontweight="bold",
        pad=5,
    )
    axes[0].text(
        1.0,
        1.02,
        LAYER_A_NOTICE,
        transform=axes[0].transAxes,
        ha="right",
        va="bottom",
        fontsize=6.4,
        color=PALETTE["calculated"],
        fontweight="bold",
    )
    axes[1].set_title(
        f"{nucleus_math(nucleus)} NMR | synthetic acquisition appearance",
        loc="left",
        fontsize=8,
        fontweight="bold",
        pad=5,
    )
    axes[1].text(
        1.0,
        1.02,
        LAYER_B_NOTICE,
        transform=axes[1].transAxes,
        ha="right",
        va="bottom",
        fontsize=6.4,
        color=PALETTE["warning"],
        fontweight="bold",
    )
    generated = [row for row in components if row.generation_status == "generated"]
    y_top = max(1.0, float(np.max(result["final_normalized"])))
    for component_index, component in enumerate(generated):
        assert component.center_ppm is not None
        if not (display_low <= component.center_ppm <= display_high):
            continue
        color = PALETTE["water"] if "water" in component.component_class else PALETTE["solvent"]
        axes[1].axvline(component.center_ppm, color=color, lw=0.55, ls=(0, (2, 2)), alpha=0.75)
        axes[1].text(
            component.center_ppm,
            y_top * (0.91 - 0.13 * (component_index % 2)),
            component.component_name,
            color=color,
            rotation=90,
            ha="right",
            va="top",
            fontsize=5.7,
        )
    axes[1].set_xlabel("Chemical shift (ppm)", fontsize=7.5)
    axes[0].set_ylabel("Calculated\nintensity", fontsize=6.5, labelpad=8)
    axes[1].set_ylabel("Relative\nintensity", fontsize=6.5, labelpad=8)
    axes[0].set_ylim(min(-0.03, float(np.min(result["layer_a"])) * 1.03), 1.08)
    lower = min(-0.05, float(np.percentile(result["final_normalized"], 0.1)) * 1.2)
    upper = max(1.05, float(np.max(result["final_normalized"])) * 1.05)
    axes[1].set_ylim(lower, upper)
    fig.suptitle(f"{spectrum_id} | {display_name} | {profile.solvent}", x=0.09, ha="left", fontsize=8.2)
    fig.text(
        0.5,
        0.012,
        "Not acquired. Identity and purity require experimental characterization.",
        ha="center",
        va="bottom",
        fontsize=6.2,
        color=PALETTE["warning"],
    )
    fig.subplots_adjust(left=0.09, right=0.985, top=0.88, bottom=0.17, hspace=0.36)
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for suffix, fmt in ((".svg", "svg"), (".png", "png")):
        path = output_stem.with_suffix(suffix)
        temporary = path.with_name(f".{path.stem}.tmp.{uuid.uuid4().hex}{suffix}")
        try:
            fig.savefig(
                temporary,
                format=fmt,
                dpi=dpi if fmt == "png" else None,
                facecolor="white",
                metadata={
                    "Title": f"{spectrum_id} {nucleus} NMR calculated and synthetic layers",
                    "Description": f"{LAYER_A_NOTICE}; {LAYER_B_NOTICE}",
                },
            )
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)
        outputs.append(path)
    plt.close(fig)
    return outputs


def profile_rows(axis: np.ndarray, result: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    layer_a: list[dict[str, Any]] = []
    layer_b: list[dict[str, Any]] = []
    normalization = float(result["normalization"])
    for index, shift in enumerate(axis):
        layer_a.append(
            {
                "notice": LAYER_A_NOTICE,
                "point_index": index,
                "chemical_shift_ppm": f"{shift:.9f}",
                "calculated_analyte_intensity_normalized": f"{result['layer_a'][index]:.12g}",
                "data_status": "CALCULATED",
            }
        )
        layer_b.append(
            {
                "notice": LAYER_B_NOTICE,
                "point_index": index,
                "chemical_shift_ppm": f"{shift:.9f}",
                "analyte_absorption_arbitrary_units": f"{np.real(result['analyte_processed'][index]):.12g}",
                "named_components_absorption_arbitrary_units": f"{np.real(result['named_processed'][index]):.12g}",
                "complex_noise_real_arbitrary_units": f"{np.real(result['noise_processed'][index]):.12g}",
                "complex_noise_imaginary_arbitrary_units": f"{np.imag(result['noise_processed'][index]):.12g}",
                "slow_baseline_arbitrary_units": f"{result['baseline'][index]:.12g}",
                "final_real_arbitrary_units": f"{result['final_real'][index]:.12g}",
                "final_imaginary_arbitrary_units": f"{result['final_imag'][index]:.12g}",
                "final_real_normalized": f"{result['final_real'][index] / normalization:.12g}",
                "data_status": "SYNTHETIC_NOT_ACQUIRED",
            }
        )
    return layer_a, layer_b


def fid_rows(profile: AcquisitionProfile, result: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index in range(profile.acquired_complex_points):
        rows.append(
            {
                "notice": LAYER_B_NOTICE,
                "point_index": index,
                "time_s": f"{index / profile.spectral_width_hz:.12g}",
                "analyte_real": f"{np.real(result['analyte_fid'][index]):.12g}",
                "analyte_imaginary": f"{np.imag(result['analyte_fid'][index]):.12g}",
                "named_components_real": f"{np.real(result['named_fid'][index]):.12g}",
                "named_components_imaginary": f"{np.imag(result['named_fid'][index]):.12g}",
                "complex_noise_real": f"{np.real(result['noise_fid'][index]):.12g}",
                "complex_noise_imaginary": f"{np.imag(result['noise_fid'][index]):.12g}",
                "total_real_before_window": f"{np.real(result['total_fid'][index]):.12g}",
                "total_imaginary_before_window": f"{np.imag(result['total_fid'][index]):.12g}",
                "window_value": f"{result['window'][index]:.12g}",
                "data_status": "SYNTHETIC_NOT_ACQUIRED",
            }
        )
    return rows


def ledger_rows(
    spectrum_id: str,
    nucleus: str,
    transitions: Sequence[Transition],
    components: Sequence[NamedComponent],
    result: dict[str, Any],
    seeds: dict[str, int],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for transition in transitions:
        rows.append(
            {
                "notice": LAYER_A_NOTICE,
                "spectrum_id": spectrum_id,
                "nucleus": nucleus,
                "component_id": transition.component_id,
                "component_name": transition.assignment,
                "component_class": "analyte_calculated_transition",
                "transition_id": transition.transition_id,
                "center_ppm": transition.shift_ppm,
                "fwhm_hz": transition.linewidth_hz,
                "relative_input_intensity": transition.intensity,
                "target_height_rel_analyte_percent": "",
                "fid_amplitude_scale": transition.intensity,
                "provenance_status": "CALCULATED",
                "calibration_record": "input stick/transition CSV",
                "generation_status": "generated",
                "deterministic_seed": "",
            }
        )
    for component in components:
        rows.append(
            {
                "notice": LAYER_B_NOTICE,
                "spectrum_id": spectrum_id,
                "nucleus": nucleus,
                "component_id": component.component_id,
                "component_name": component.component_name,
                "component_class": component.component_class,
                "transition_id": "",
                "center_ppm": "" if component.center_ppm is None else component.center_ppm,
                "fwhm_hz": "" if component.fwhm_hz is None else component.fwhm_hz,
                "relative_input_intensity": "",
                "target_height_rel_analyte_percent": ""
                if component.target_height_rel_analyte_percent is None
                else component.target_height_rel_analyte_percent,
                "fid_amplitude_scale": result["component_scales"].get(component.component_id, 0.0),
                "provenance_status": component.provenance_status,
                "calibration_record": component.calibration_record,
                "generation_status": component.generation_status,
                "deterministic_seed": seeds.get(component.component_id, ""),
            }
        )
    for component_id, component_name, component_class, seed_key, provenance in (
        (
            "complex_fid_noise",
            "correlated circular complex Gaussian acquisition noise",
            "acquisition_noise",
            "complex_fid_noise",
            "local DMSO amplitude/processed-lag target; complex Gaussian frequency AR(1) projected to an acquisition FID",
        ),
        (
            "exponential_window",
            "exponential apodization",
            "processing_window",
            "",
            "acquisition profile line-broadening parameter",
        ),
        (
            "zero_filling",
            "zero filling before FFT",
            "digital_processing",
            "",
            "acquisition profile TD/SI parameters",
        ),
        (
            "residual_phase",
            "residual zero/first-order phase",
            "processing_phase",
            "phase",
            "deterministic assumed post-phasing residual",
        ),
        (
            "slow_baseline",
            "named slow processing baseline",
            "processing_baseline",
            "slow_baseline",
            "local DMSO baseline span target or explicit assumed preset",
        ),
    ):
        rows.append(
            {
                "notice": LAYER_B_NOTICE,
                "spectrum_id": spectrum_id,
                "nucleus": nucleus,
                "component_id": component_id,
                "component_name": component_name,
                "component_class": component_class,
                "transition_id": "",
                "center_ppm": "",
                "fwhm_hz": "",
                "relative_input_intensity": "",
                "target_height_rel_analyte_percent": "",
                "fid_amplitude_scale": result["noise_scale"] if component_id == "complex_fid_noise" else "",
                "provenance_status": provenance,
                "calibration_record": "named processing component; never an impurity peak",
                "generation_status": "generated",
                "deterministic_seed": seeds.get(seed_key, "") if seed_key else "",
            }
        )
    return rows


def run_generation(args: argparse.Namespace) -> Path:
    input_path = args.input.expanduser().resolve()
    calibration_dir = args.calibration_dir.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()
    solvent = normalize_solvent(args.solvent)
    nucleus = args.nucleus
    transitions = load_transitions(input_path, nucleus)
    input_hash = sha256_file(input_path)
    spectrum_seed = derive_seed(
        SCHEMA_VERSION,
        args.seed,
        args.spectrum_id,
        nucleus,
        solvent,
        input_hash,
    )
    seeds = {
        key: derive_seed(spectrum_seed, key)
        for key in (
            "acquisition_profile",
            "residual_dmso_d6_1h",
            "water_in_dmso_d6_1h",
            "hdo_in_d2o_1h",
            "complex_fid_noise",
            "phase",
            "slow_baseline",
        )
    }
    profile = load_acquisition_profile(
        calibration_dir,
        solvent,
        nucleus,
        seeds["acquisition_profile"],
        args.acquired_points,
        args.processed_points,
    )
    components = named_components(calibration_dir, solvent, nucleus, spectrum_seed)
    result = build_appearance(
        transitions,
        profile,
        components,
        spectrum_seed,
        args.phase0_deg,
        args.phase1_deg,
    )

    safe_nucleus = nucleus.replace("/", "_")
    output_dir = output_root / args.spectrum_id / safe_nucleus
    if output_dir.exists() and any(output_dir.iterdir()) and not args.overwrite:
        raise AppearanceError(f"Output directory is not empty; use --overwrite: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    layer_a_rows, layer_b_rows = profile_rows(result["axis_ppm"], result)
    layer_a_fields = list(layer_a_rows[0])
    layer_b_fields = list(layer_b_rows[0])
    raw_fid_rows = fid_rows(profile, result)
    raw_fid_fields = list(raw_fid_rows[0])
    ledger = ledger_rows(args.spectrum_id, nucleus, transitions, components, result, seeds)
    ledger_fields = list(ledger[0])
    transition_rows = [
        {
            "notice": LAYER_A_NOTICE,
            **asdict(row),
            "source_csv_sha256": input_hash,
        }
        for row in transitions
    ]
    transition_fields = list(transition_rows[0])

    written: list[Path] = []
    for path, rows, fields, writer in (
        (output_dir / "layer_A_calculated_model.csv", layer_a_rows, layer_a_fields, atomic_write_csv),
        (output_dir / "layer_A_calculated_model.txt", layer_a_rows, layer_a_fields, atomic_write_tsv),
        (output_dir / "layer_B_synthetic_acquisition.csv", layer_b_rows, layer_b_fields, atomic_write_csv),
        (output_dir / "layer_B_synthetic_acquisition.txt", layer_b_rows, layer_b_fields, atomic_write_tsv),
        (output_dir / "layer_B_complex_FID.csv", raw_fid_rows, raw_fid_fields, atomic_write_csv),
        (output_dir / "component_ledger.csv", ledger, ledger_fields, atomic_write_csv),
        (output_dir / "input_transitions_normalized.csv", transition_rows, transition_fields, atomic_write_csv),
    ):
        writer(path, rows, fields)
        written.append(path)

    display_high, display_low = default_display_limits(
        nucleus,
        profile,
        transitions,
        args.display_high_ppm,
        args.display_low_ppm,
    )
    figure_paths = render_figure(
        output_dir / f"{args.spectrum_id}_{nucleus}_NMR_layers",
        args.spectrum_id,
        args.display_name,
        nucleus,
        profile,
        transitions,
        components,
        result,
        display_high,
        display_low,
        args.dpi,
    )
    written.extend(figure_paths)

    metadata = {
        "schema_version": SCHEMA_VERSION,
        "generated_utc": utc_now(),
        "generated_by": "scripts/nmr_acquisition_appearance.py",
        "spectrum_id": args.spectrum_id,
        "display_name": args.display_name,
        "nucleus": nucleus,
        "solvent": solvent,
        "layer_A_notice": LAYER_A_NOTICE,
        "layer_B_notice": LAYER_B_NOTICE,
        "global_warning": GLOBAL_WARNING,
        "scientific_scope": {
            "layer_A": "calculated analyte transitions convolved only with named input linewidths",
            "layer_B": (
                "synthetic acquisition appearance from complex analyte/named-component FID, "
                "complex Gaussian noise, exponential apodization, zero filling, FFT, residual "
                "phase and named slow baseline"
            ),
            "forbidden_interpretation": [
                "not experimental NMR",
                "not identity or purity evidence",
                "not a measured impurity profile",
                "not a replacement for raw vendor FID",
            ],
            "anonymous_impurity_peaks_added": False,
        },
        "input": {
            "file": input_path.name,
            "sha256": input_hash,
            "transition_count": len(transitions),
        },
        "deterministic_seeding": {
            "user_seed": args.seed,
            "spectrum_seed": spectrum_seed,
            "independent_component_seeds": seeds,
        },
        "acquisition_profile": asdict(profile),
        "named_components": [asdict(row) for row in components],
        "processing": {
            "window": "exponential exp(-pi*LB*t)",
            "zero_fill_from_points": profile.acquired_complex_points,
            "zero_fill_to_points": profile.processed_points,
            "fft": "numpy complex FFT with fftshift",
            "phase0_deg": result["phase0_deg"],
            "phase1_deg": result["phase1_deg"],
            "complex_fid_noise_model": (
                "circular complex Gaussian frequency-domain AR(1), inverse-FFT projected "
                "to a truncated acquisition FID and passed through the full processing chain"
            ),
            "noise_spectral_ar1_coefficient_before_time_projection": result[
                "noise_spectral_ar1_coefficient_before_time_projection"
            ],
            "noise_target_stationary_sigma_rel_max_percent": profile.noise_stationary_rel_max_percent,
            "noise_achieved_sigma_arbitrary_units": result["noise_achieved_sigma"],
            "noise_achieved_frequency_domain_lag1": result["noise_achieved_lag1"],
            "noise_target_frequency_domain_lag1_from_calibration": profile.noise_target_lag1,
            "noise_lag_note": (
                "The local lag-1 value is a calibration target, not forced exactly; the achieved "
                "value emerges from acquisition noise, truncation, apodization and zero filling."
            ),
            "baseline_target_span_rel_max_percent": profile.baseline_span_rel_max_percent,
            "baseline_achieved_span_arbitrary_units": result["baseline_achieved_span"],
            "baseline_achieved_rms_arbitrary_units": result["baseline_achieved_rms"],
        },
        "calibration_sources": calibration_hashes(calibration_dir),
        "display": {
            "archetype": "quantitative two-layer comparison",
            "core_conclusion": "computed transitions and synthetic acquisition appearance are visibly and provenance-wise distinct",
            "display_high_ppm": display_high,
            "display_low_ppm": display_low,
            "font": "Arial",
            "png_dpi": args.dpi,
        },
    }
    metadata_path = output_dir / "spectrum_metadata.json"
    atomic_write_json(metadata_path, metadata)
    written.append(metadata_path)

    readme = "\n".join(
        [
            LAYER_A_NOTICE,
            LAYER_B_NOTICE,
            "",
            f"Spectrum: {args.spectrum_id} | {args.display_name}",
            f"Nucleus/solvent: {nucleus} / {solvent}",
            f"Calibration status: {profile.calibration_status}",
            f"Deterministic spectrum seed: {spectrum_seed}",
            "",
            "Layer A contains only the calculated analyte transitions supplied in the input CSV.",
            "Layer B is synthetic and adds only named solvent/water components plus named acquisition/processing effects.",
            "No anonymous impurity peaks are generated.",
            GLOBAL_WARNING,
            "",
            "Numeric outputs:",
            "- layer_A_calculated_model.csv/.txt",
            "- layer_B_synthetic_acquisition.csv/.txt",
            "- layer_B_complex_FID.csv",
            "- component_ledger.csv",
            "- input_transitions_normalized.csv",
            "",
        ]
    )
    readme_path = output_dir / "README_CALCULATED_AND_SYNTHETIC_ONLY.txt"
    atomic_write_text(readme_path, readme)
    written.append(readme_path)

    hash_rows = [
        {
            "file": path.name,
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "notice": LAYER_A_NOTICE if "layer_A" in path.name else LAYER_B_NOTICE,
        }
        for path in sorted(written, key=lambda item: item.name)
    ]
    hash_path = output_dir / "output_file_hashes.csv"
    atomic_write_csv(hash_path, hash_rows, list(hash_rows[0]))
    completion = {
        "status": "COMPLETE_SYNTHETIC_NOT_ACQUIRED",
        "generated_utc": utc_now(),
        "output_hash_manifest": hash_path.name,
        "output_hash_manifest_sha256": sha256_file(hash_path),
        "file_count_excluding_completion_marker": len(hash_rows) + 1,
        "layer_A_notice": LAYER_A_NOTICE,
        "layer_B_notice": LAYER_B_NOTICE,
    }
    atomic_write_json(output_dir / "run_complete.json", completion)
    return output_dir


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description=(
            "Generate auditable Layer-A calculated NMR and Layer-B synthetic acquisition "
            "appearance from a standard transition CSV. This script never calls ORCA."
        )
    )
    result.add_argument("--input", type=Path, required=True, help="Standard stick/transition CSV")
    result.add_argument("--spectrum-id", required=True, help="Stable ID such as F3_019")
    result.add_argument("--display-name", required=True)
    result.add_argument("--nucleus", choices=sorted(SUPPORTED_NUCLEI), required=True)
    result.add_argument("--solvent", choices=["DMSO-d6", "D2O"], required=True)
    result.add_argument("--calibration-dir", type=Path, default=DEFAULT_CALIBRATION_DIR)
    result.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    result.add_argument("--seed", type=int, default=20260720)
    result.add_argument("--phase0-deg", type=float, default=None)
    result.add_argument("--phase1-deg", type=float, default=None)
    result.add_argument("--display-high-ppm", type=float, default=None)
    result.add_argument("--display-low-ppm", type=float, default=None)
    result.add_argument("--acquired-points", type=int, default=None, help="Explicit test/QA override")
    result.add_argument("--processed-points", type=int, default=None, help="Explicit test/QA override")
    result.add_argument("--dpi", type=int, default=600)
    result.add_argument("--overwrite", action="store_true")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        output = run_generation(args)
    except (AppearanceError, FileNotFoundError, OSError, ValueError) as exc:
        print(f"NMR appearance generation failed: {exc}", file=sys.stderr)
        return 2
    print(LAYER_A_NOTICE)
    print(LAYER_B_NOTICE)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
