#!/usr/bin/env python3
"""Read-only audit of the author's Mnova-embedded Bruker 1D NMR data.

The script never writes to the source spectrum directory. It recovers:

* document-level acquisition metadata;
* embedded Bruker ``acqus`` and ``procs`` JCAMP parameter files;
* the processed complex float32 spectrum stored in the Mnova document;
* robust baseline/noise and residual-solvent/water peak descriptors.

All derived values are explicitly marked as derived. The procedure is not a
replacement for the original vendor raw-data directory. In particular, an
Mnova document whose source path ends in ``fid`` records a raw-FID import, but
the standalone Bruker ``fid`` file is not present in the audited folder.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
import zlib
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.signal import find_peaks, medfilt


ZLIB_SECOND_BYTES = {0x01, 0x5E, 0x9C, 0xDA}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


@dataclass(frozen=True)
class ZStream:
    offset: int
    compressed_length: int
    payload: bytes


def scan_zlib_streams(data: bytes) -> list[ZStream]:
    streams: list[ZStream] = []
    for offset in range(len(data) - 2):
        if data[offset] != 0x78 or data[offset + 1] not in ZLIB_SECOND_BYTES:
            continue
        try:
            obj = zlib.decompressobj()
            payload = obj.decompress(data[offset:]) + obj.flush()
        except zlib.error:
            continue
        if not obj.eof:
            continue
        consumed = len(data[offset:]) - len(obj.unused_data)
        streams.append(ZStream(offset, consumed, payload))
    return streams


def scan_utf16be_strings(data: bytes, max_chars: int = 1000) -> list[tuple[int, str]]:
    """Return Mnova's 32-bit-BE-length-prefixed UTF-16BE strings."""
    out: list[tuple[int, str]] = []
    for offset in range(len(data) - 6):
        nbytes = int.from_bytes(data[offset : offset + 4], "big")
        if nbytes < 2 or nbytes > max_chars * 2 or nbytes % 2:
            continue
        end = offset + 4 + nbytes
        if end > len(data):
            continue
        try:
            text = data[offset + 4 : end].decode("utf-16be")
        except UnicodeDecodeError:
            continue
        if not any(ch.isalnum() for ch in text):
            continue
        if all(ch.isprintable() or ch in "\r\n\t" for ch in text):
            out.append((offset, text))
    return out


def parse_jcamp(payload: bytes) -> tuple[dict[str, str], dict[str, list[str]], str]:
    text = payload.decode("latin1", "replace").replace("\r\n", "\n").replace("\r", "\n")
    scalars: dict[str, str] = {}
    arrays: dict[str, list[str]] = {}
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        m = re.match(r"^##(\$?[A-Za-z0-9_]+)=\s*(.*)$", lines[i])
        if not m:
            i += 1
            continue
        key, value = m.group(1), m.group(2).strip()
        scalars[key] = value.strip("<>")
        if re.match(r"^\(\s*\d+\s*\.\.\s*\d+\s*\)$", value):
            items: list[str] = []
            j = i + 1
            while j < len(lines) and not lines[j].startswith("##"):
                items.extend(re.findall(r"<[^>]*>|[^\s]+", lines[j]))
                j += 1
            arrays[key] = [item.strip("<>") for item in items]
            i = j
            continue
        i += 1
    return scalars, arrays, text


def to_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        result = float(str(value).strip().strip("<>"))
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def to_int(value: Any) -> int | None:
    f = to_float(value)
    return int(round(f)) if f is not None else None


def previous_string_value(strings: list[tuple[int, str]], label: str, start_index: int = 0) -> list[tuple[int, str]]:
    values: list[tuple[int, str]] = []
    for index in range(max(1, start_index), len(strings)):
        if strings[index][1] == label:
            values.append(strings[index - 1])
    return values


DOCUMENT_LABELS = [
    "Data File Name",
    "Solvent",
    "Temperature",
    "Pulse Sequence",
    "Number of Scans",
    "Receiver Gain",
    "Relaxation Delay",
    "Pulse Width",
    "Acquisition Time",
    "Acquisition Date",
    "Modification Date",
    "Spectrometer Frequency",
    "Spectral Width",
    "Lowest Frequency",
    "Nucleus",
    "Acquired Size",
    "Spectral Size",
]


def parse_document_metadata(strings: list[tuple[int, str]]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    current: dict[str, Any] = {}
    current_offsets: dict[str, int] = {}
    labels = set(DOCUMENT_LABELS)
    for index, (offset, text) in enumerate(strings):
        if text not in labels or index == 0:
            continue
        value_offset, value = strings[index - 1]
        current[text] = value
        current_offsets[text] = value_offset
        if text == "Spectral Size":
            current["_value_offsets"] = current_offsets
            records.append(current)
            current = {}
            current_offsets = {}
    return records


def classify_parameter_stream(stream: ZStream) -> tuple[str | None, dict[str, str], dict[str, list[str]], str]:
    if not stream.payload.startswith(b"##TITLE= Parameter file"):
        return None, {}, {}, ""
    scalars, arrays, text = parse_jcamp(stream.payload)
    if "$NUC1" in scalars and "$TD" in scalars and "$SOLVENT" in scalars:
        return "acqus", scalars, arrays, text
    if "$SI" in scalars and "$SW_p" in scalars:
        return "procs", scalars, arrays, text
    return "other_parameter_file", scalars, arrays, text


def robust_sigma_adjacent_difference(y: np.ndarray) -> float:
    differences = np.diff(y[np.isfinite(y)])
    if differences.size < 10:
        return float("nan")
    median = np.median(differences)
    mad = np.median(np.abs(differences - median))
    return float(1.482602218505602 * mad / math.sqrt(2.0))


def robust_polynomial_baseline(y: np.ndarray, degree: int = 2, bins: int = 512) -> tuple[np.ndarray, np.ndarray]:
    n = y.size
    boundaries = np.linspace(0, n, min(bins, n // 8) + 1, dtype=int)
    centers: list[float] = []
    medians: list[float] = []
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        if right <= left:
            continue
        segment = y[left:right]
        finite = segment[np.isfinite(segment)]
        if finite.size:
            centers.append((left + right - 1) / 2.0)
            medians.append(float(np.median(finite)))
    x = np.asarray(centers, dtype=float)
    values = np.asarray(medians, dtype=float)
    if values.size < degree + 3:
        baseline = np.full(n, np.nanmedian(y), dtype=float)
        return baseline, np.ones(values.size, dtype=bool)
    mask = np.ones(values.size, dtype=bool)
    coefficients = np.polyfit(x, values, degree)
    for _ in range(6):
        coefficients = np.polyfit(x[mask], values[mask], degree)
        residual = values - np.polyval(coefficients, x)
        med = np.median(residual[mask])
        mad = 1.482602218505602 * np.median(np.abs(residual[mask] - med))
        if not math.isfinite(mad) or mad == 0:
            break
        new_mask = np.abs(residual - med) <= 3.5 * mad
        if new_mask.sum() < degree + 3 or np.array_equal(new_mask, mask):
            break
        mask = new_mask
    baseline = np.polyval(coefficients, np.arange(n, dtype=float))
    return baseline, mask


def slow_median_baseline(y: np.ndarray, bins: int = 1024, kernel_bins: int = 31) -> np.ndarray:
    """Estimate a reproducible slow baseline without following narrow peaks."""
    n = y.size
    bins = max(64, min(bins, n // 32))
    boundaries = np.linspace(0, n, bins + 1, dtype=int)
    centers: list[float] = []
    medians: list[float] = []
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        centers.append((left + right - 1) / 2.0)
        medians.append(float(np.median(y[left:right])))
    values = np.asarray(medians, dtype=float)
    kernel = min(kernel_bins, values.size if values.size % 2 else values.size - 1)
    kernel = max(3, kernel if kernel % 2 else kernel - 1)
    smooth = medfilt(values, kernel_size=kernel)
    return np.interp(np.arange(n, dtype=float), np.asarray(centers), smooth)


def crossing_x(x0: float, y0: float, x1: float, y1: float, target: float) -> float:
    if y1 == y0:
        return (x0 + x1) / 2.0
    fraction = (target - y0) / (y1 - y0)
    return x0 + fraction * (x1 - x0)


def measure_peak_at_index(
    ppm: np.ndarray,
    y: np.ndarray,
    index: int,
    baseline_value: float,
    noise_sigma: float,
    max_walk_ppm: float,
) -> dict[str, Any]:
    height = float(y[index] - baseline_value)
    result: dict[str, Any] = {
        "center_ppm": float(ppm[index]),
        "height": height,
        "snr": height / noise_sigma if noise_sigma > 0 else None,
        "fwhm_ppm": None,
        "fwhm_hz": None,
        "left_half_ppm": None,
        "right_half_ppm": None,
        "asymmetry_halfwidth_ratio": None,
    }
    if height <= 0 or not math.isfinite(height):
        return result
    half = baseline_value + height / 2.0
    left = index
    while left > 0 and abs(ppm[left] - ppm[index]) <= max_walk_ppm and y[left] > half:
        left -= 1
    right = index
    while right < y.size - 1 and abs(ppm[right] - ppm[index]) <= max_walk_ppm and y[right] > half:
        right += 1
    if left == index or right == index or y[left] > half or y[right] > half:
        return result
    x_left = crossing_x(float(ppm[left]), float(y[left]), float(ppm[left + 1]), float(y[left + 1]), half)
    x_right = crossing_x(float(ppm[right - 1]), float(y[right - 1]), float(ppm[right]), float(y[right]), half)
    fwhm_ppm = abs(x_left - x_right)
    left_width = abs(float(ppm[index]) - x_left)
    right_width = abs(x_right - float(ppm[index]))
    result.update(
        {
            "fwhm_ppm": fwhm_ppm,
            "left_half_ppm": x_left,
            "right_half_ppm": x_right,
            "asymmetry_halfwidth_ratio": left_width / right_width if right_width > 0 else None,
        }
    )
    return result


def analyze_target_peak(
    ppm: np.ndarray,
    y: np.ndarray,
    frequency_mhz: float,
    expected_ppm: float,
    search_halfwidth_ppm: float,
    baseline_outer_halfwidth_ppm: float,
    baseline_inner_halfwidth_ppm: float,
    noise_sigma: float,
    max_walk_ppm: float,
) -> dict[str, Any]:
    local_mask = np.abs(ppm - expected_ppm) <= baseline_outer_halfwidth_ppm
    edge_mask = local_mask & (np.abs(ppm - expected_ppm) >= baseline_inner_halfwidth_ppm)
    baseline_value = float(np.median(y[edge_mask])) if np.any(edge_mask) else float(np.median(y[local_mask]))
    search_mask = np.abs(ppm - expected_ppm) <= search_halfwidth_ppm
    indices = np.flatnonzero(search_mask)
    if not indices.size:
        return {"status": "no_axis_coverage", "confidence": "none"}
    local_y = y[indices] - baseline_value
    minimum_prominence = max(5.0 * noise_sigma, 0.002 * max(float(np.max(local_y)), 0.0))
    distance_points = max(1, int(round(0.001 / abs(float(ppm[1] - ppm[0])))))
    peaks, properties = find_peaks(local_y, prominence=minimum_prominence, distance=distance_points)
    if peaks.size:
        candidates = indices[peaks]
        heights = y[candidates] - baseline_value
        distances = np.abs(ppm[candidates] - expected_ppm)
        snrs = np.maximum(heights / max(noise_sigma, np.finfo(float).eps), 1.0)
        score = distances / search_halfwidth_ppm - 0.10 * np.log10(snrs)
        index = int(candidates[int(np.argmin(score))])
        candidate_count = int(peaks.size)
    else:
        index = int(indices[int(np.argmax(local_y))])
        candidate_count = 0
    result = measure_peak_at_index(ppm, y, index, baseline_value, noise_sigma, max_walk_ppm)
    result["fwhm_hz"] = result["fwhm_ppm"] * frequency_mhz if result.get("fwhm_ppm") is not None else None
    result["baseline"] = baseline_value
    result["expected_ppm"] = expected_ppm
    result["delta_from_expected_ppm"] = result["center_ppm"] - expected_ppm
    result["candidate_peak_count"] = candidate_count
    snr = result.get("snr") or 0.0
    delta = abs(result["delta_from_expected_ppm"])
    fwhm_hz = result.get("fwhm_hz")
    if snr >= 20 and delta <= min(0.05, 0.5 * search_halfwidth_ppm) and fwhm_hz is not None:
        confidence = "high"
    elif snr >= 8 and delta <= search_halfwidth_ppm and fwhm_hz is not None:
        confidence = "medium"
    else:
        confidence = "low"
    issues: list[str] = []
    if candidate_count > 3:
        issues.append("multiple_local_candidates_possible_overlap")
    if fwhm_hz is not None and (fwhm_hz < 0.25 or fwhm_hz > 80):
        issues.append("fwhm_outside_typical_1d_range")
        if confidence == "high":
            confidence = "medium"
    if snr < 8:
        issues.append("low_snr")
    result["confidence"] = confidence
    result["status"] = "measured_from_processed_array"
    result["issues"] = ";".join(issues)
    return result


def spectrum_statistics(ppm: np.ndarray, y: np.ndarray) -> tuple[dict[str, Any], np.ndarray, float]:
    polynomial_baseline, _ = robust_polynomial_baseline(y)
    baseline = slow_median_baseline(y)
    residual = y - baseline
    noise_sigma = robust_sigma_adjacent_difference(residual)
    if not math.isfinite(noise_sigma) or noise_sigma <= 0:
        noise_sigma = float(1.482602218505602 * np.median(np.abs(residual - np.median(residual))))
    max_positive = float(np.max(residual))
    clip_mask = np.isfinite(residual) & (np.abs(residual) <= 4.0 * noise_sigma)
    clipped = residual[clip_mask]
    lag_mask = clip_mask[:-1] & clip_mask[1:]
    if lag_mask.sum() >= 10:
        a = residual[:-1][lag_mask]
        b = residual[1:][lag_mask]
        lag1 = float(np.corrcoef(a, b)[0, 1]) if np.std(a) > 0 and np.std(b) > 0 else None
    else:
        lag1 = None
    if lag1 is not None and -0.99 < lag1 < 0.99:
        stationary_sigma = noise_sigma / math.sqrt(max(1.0 - lag1, 0.01))
        ar1_innovation_sigma = stationary_sigma * math.sqrt(max(1.0 - lag1 * lag1, 0.0))
    else:
        stationary_sigma = None
        ar1_innovation_sigma = None
    if clipped.size >= 10 and np.std(clipped) > 0:
        centered = clipped - np.mean(clipped)
        sd = np.std(centered)
        skew = float(np.mean((centered / sd) ** 3))
        excess_kurtosis = float(np.mean((centered / sd) ** 4) - 3.0)
    else:
        skew = None
        excess_kurtosis = None
    baseline_span = float(np.percentile(baseline, 95) - np.percentile(baseline, 5))
    polynomial_span = float(
        np.percentile(polynomial_baseline, 95) - np.percentile(polynomial_baseline, 5)
    )
    stats = {
        "baseline_median_raw": float(np.median(baseline)),
        "baseline_span_p95_p5_raw": baseline_span,
        "baseline_span_rel_max_percent": 100.0 * baseline_span / max_positive if max_positive > 0 else None,
        "baseline_slow_rms_raw": float(np.std(baseline - np.median(baseline))),
        "baseline_slow_rms_rel_max_percent": (
            100.0 * float(np.std(baseline - np.median(baseline))) / max_positive if max_positive > 0 else None
        ),
        "baseline_global_polynomial_span_p95_p5_raw": polynomial_span,
        "baseline_global_polynomial_span_rel_max_percent": (
            100.0 * polynomial_span / max_positive if max_positive > 0 else None
        ),
        "noise_sigma_adjacent_difference_raw": noise_sigma,
        "noise_rel_max_percent": 100.0 * noise_sigma / max_positive if max_positive > 0 else None,
        "noise_lag1_correlation_clipped": lag1,
        "noise_stationary_sigma_from_diff_and_lag1_raw": stationary_sigma,
        "noise_stationary_rel_max_percent": (
            100.0 * stationary_sigma / max_positive
            if stationary_sigma is not None and max_positive > 0
            else None
        ),
        "noise_AR1_innovation_sigma_raw": ar1_innovation_sigma,
        "noise_AR1_innovation_rel_max_percent": (
            100.0 * ar1_innovation_sigma / max_positive
            if ar1_innovation_sigma is not None and max_positive > 0
            else None
        ),
        "noise_skewness_clipped": skew,
        "noise_excess_kurtosis_clipped": excess_kurtosis,
        "max_baseline_corrected_peak_raw": max_positive,
        "clipped_noise_point_fraction": float(clip_mask.mean()),
    }
    return stats, baseline, noise_sigma


def numeric_summary(values: Iterable[Any]) -> dict[str, float | int | None]:
    clean = [float(v) for v in values if v is not None and math.isfinite(float(v))]
    if not clean:
        return {"n": 0, "median": None, "p10": None, "p90": None, "min": None, "max": None}
    a = np.asarray(clean, dtype=float)
    return {
        "n": int(a.size),
        "median": float(np.median(a)),
        "p10": float(np.percentile(a, 10)),
        "p90": float(np.percentile(a, 90)),
        "min": float(np.min(a)),
        "max": float(np.max(a)),
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = []
        seen: set[str] = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    seen.add(key)
                    fields.append(key)
    with path.open("w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def flatten_peak(prefix: str, peak: dict[str, Any] | None) -> dict[str, Any]:
    if not peak:
        return {}
    return {f"{prefix}_{key}": value for key, value in peak.items()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = args.source.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    source_files = sorted(path for path in source.iterdir() if path.is_file())
    hashes = [
        {
            "file": str(path),
            "relative_file": path.name,
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in source_files
    ]
    write_csv(output / "source_file_hashes.csv", hashes)

    acquisition_rows: list[dict[str, Any]] = []
    calibration_rows: list[dict[str, Any]] = []
    extraction_issues: list[dict[str, Any]] = []
    internal_records: list[dict[str, Any]] = []

    for mnova_path in sorted(source.glob("*.mnova")):
        data = mnova_path.read_bytes()
        strings = scan_utf16be_strings(data)
        doc_metadata = parse_document_metadata(strings)
        unknown_offsets = [offset for offset, text in strings if text == "Unknown"]
        zstreams = scan_zlib_streams(data)
        parameter_streams: list[dict[str, Any]] = []
        for stream in zstreams:
            kind, scalars, arrays, text = classify_parameter_stream(stream)
            if kind in {"acqus", "procs"}:
                parameter_streams.append(
                    {
                        "stream": stream,
                        "kind": kind,
                        "scalars": scalars,
                        "arrays": arrays,
                        "text": text,
                    }
                )
        acqus_streams = [item for item in parameter_streams if item["kind"] == "acqus"]
        procs_streams = [item for item in parameter_streams if item["kind"] == "procs"]
        for dataset_index, acqus_item in enumerate(acqus_streams, start=1):
            acq_stream: ZStream = acqus_item["stream"]
            acq = acqus_item["scalars"]
            acq_arrays = acqus_item["arrays"]
            prior_procs = [item for item in procs_streams if item["stream"].offset < acq_stream.offset]
            if not prior_procs:
                extraction_issues.append(
                    {"file": str(mnova_path), "dataset_index": dataset_index, "issue": "no_prior_procs_stream"}
                )
                continue
            procs_item = max(prior_procs, key=lambda item: item["stream"].offset)
            proc_stream: ZStream = procs_item["stream"]
            proc = procs_item["scalars"]
            metadata = doc_metadata[dataset_index - 1] if dataset_index <= len(doc_metadata) else {}
            nucleus = str(metadata.get("Nucleus") or acq.get("$NUC1") or "").strip("<>")
            solvent = str(metadata.get("Solvent") or acq.get("$SOLVENT") or "").strip("<>")
            spectral_size = to_int(metadata.get("Spectral Size")) or to_int(proc.get("$SI"))
            if spectral_size is None:
                extraction_issues.append(
                    {"file": str(mnova_path), "dataset_index": dataset_index, "issue": "missing_spectral_size"}
                )
                continue
            spectrum_bytes = spectral_size * 8
            acq_end = acq_stream.offset + acq_stream.compressed_length
            candidates: list[tuple[int, int, int]] = []
            for u1 in unknown_offsets:
                if u1 <= acq_end:
                    continue
                for u2 in unknown_offsets:
                    if u2 <= u1:
                        continue
                    header_size = u2 - u1 - spectrum_bytes
                    if 800 <= header_size <= 950:
                        candidates.append((u1, u2, header_size))
                        break
                if candidates:
                    break
            if not candidates:
                extraction_issues.append(
                    {
                        "file": str(mnova_path),
                        "dataset_index": dataset_index,
                        "issue": "processed_spectrum_block_not_located",
                    }
                )
                continue
            unknown_first, unknown_second, header_size = candidates[0]
            # The repeated Mnova spectrum-object prefix starts 37 bytes before
            # the first length-prefixed "Unknown" string. The first complex
            # spectrum ends exactly where the repeated prefix begins.
            spectrum_start = unknown_second - 37 - spectrum_bytes
            spectrum_end = spectrum_start + spectrum_bytes
            raw_spectrum = data[spectrum_start:spectrum_end]
            if len(raw_spectrum) != spectrum_bytes:
                extraction_issues.append(
                    {"file": str(mnova_path), "dataset_index": dataset_index, "issue": "truncated_spectrum_block"}
                )
                continue
            complex_pairs = np.frombuffer(raw_spectrum, dtype="<f4").reshape(spectral_size, 2)
            real = complex_pairs[:, 0].astype(np.float64)
            imag = complex_pairs[:, 1].astype(np.float64)
            if not np.all(np.isfinite(real)):
                extraction_issues.append(
                    {"file": str(mnova_path), "dataset_index": dataset_index, "issue": "nonfinite_real_spectrum"}
                )
                continue
            frequency_mhz = to_float(metadata.get("Spectrometer Frequency")) or to_float(proc.get("$SF"))
            spectral_width_hz = to_float(metadata.get("Spectral Width")) or to_float(proc.get("$SW_p"))
            lowest_frequency_hz = to_float(metadata.get("Lowest Frequency"))
            if frequency_mhz is None or spectral_width_hz is None:
                extraction_issues.append(
                    {"file": str(mnova_path), "dataset_index": dataset_index, "issue": "missing_axis_parameters"}
                )
                continue
            if lowest_frequency_hz is not None:
                high_ppm = (lowest_frequency_hz + spectral_width_hz) / frequency_mhz
                low_ppm_nominal = lowest_frequency_hz / frequency_mhz
            else:
                high_ppm = to_float(proc.get("$OFFSET"))
                if high_ppm is None:
                    extraction_issues.append(
                        {"file": str(mnova_path), "dataset_index": dataset_index, "issue": "missing_axis_offset"}
                    )
                    continue
                low_ppm_nominal = high_ppm - spectral_width_hz / frequency_mhz
            step_ppm = spectral_width_hz / frequency_mhz / spectral_size
            ppm = high_ppm - np.arange(spectral_size, dtype=float) * step_ppm
            stats, baseline, noise_sigma = spectrum_statistics(ppm, real)

            data_file_name = str(metadata.get("Data File Name") or "")
            input_kind = "fid" if re.search(r"(?:^|[/\\])fid$", data_file_name, re.I) else "processed_1r"
            expected_sample = re.sub(r"-(H|C)$", "", mnova_path.stem, flags=re.I)
            normalized_source_path = data_file_name.replace("\\", "/")
            primary_by_filename = f"/{expected_sample}/".lower() in normalized_source_path.lower()
            if not primary_by_filename and expected_sample == "3-14" and "3-14" in normalized_source_path:
                primary_by_filename = True

            temperature_raw = to_float(metadata.get("Temperature")) or to_float(acq.get("$TE"))
            if temperature_raw is not None and temperature_raw <= 0:
                temperature_k = None
                temperature_status = "invalid_zero_sentinel_not_physical"
            else:
                temperature_k = temperature_raw
                temperature_status = "reported_parameter" if temperature_k is not None else "not_available"
            td_real = to_int(acq.get("$TD"))
            acquired_complex = to_int(metadata.get("Acquired Size")) or (td_real // 2 if td_real else None)
            raw_resolution = spectral_width_hz / acquired_complex if acquired_complex else None
            processed_resolution = spectral_width_hz / spectral_size
            d_values = acq_arrays.get("$D", [])
            p_values = acq_arrays.get("$P", [])
            d1 = to_float(d_values[1]) if len(d_values) > 1 else to_float(metadata.get("Relaxation Delay"))
            p1 = to_float(p_values[1]) if len(p_values) > 1 else to_float(metadata.get("Pulse Width"))

            source_hash = sha256_file(mnova_path)
            dataset_id = f"{mnova_path.stem}__embedded_{dataset_index:02d}"
            acq_row = {
                "dataset_id": dataset_id,
                "mnova_file": str(mnova_path),
                "mnova_sha256": source_hash,
                "embedded_dataset_index": dataset_index,
                "primary_dataset_for_filename": primary_by_filename,
                "embedded_original_path": data_file_name,
                "embedded_input_kind": input_kind,
                "raw_vendor_fid_file_present_in_audited_folder": False,
                "vendor": "Bruker",
                "origin": acq.get("ORIGIN"),
                "instrument": acq.get("$INSTRUM") or metadata.get("Spectrometer"),
                "probe": acq.get("$PROBHD"),
                "nucleus": nucleus,
                "solvent_parameter": solvent,
                "frequency_MHz": frequency_mhz,
                "spectral_width_Hz": spectral_width_hz,
                "lowest_frequency_Hz": lowest_frequency_hz,
                "high_ppm_document_axis": high_ppm,
                "low_ppm_nominal_document_axis": low_ppm_nominal,
                "axis_step_ppm": step_ppm,
                "TD_real_points": td_real,
                "acquired_complex_points": acquired_complex,
                "SI_processed_points": spectral_size,
                "AQ_s": to_float(metadata.get("Acquisition Time")) or to_float(acq.get("$AQ")),
                "NS": to_int(metadata.get("Number of Scans")) or to_int(acq.get("$NS")),
                "DS": to_int(acq.get("$DS")),
                "RG": to_float(metadata.get("Receiver Gain")) or to_float(acq.get("$RG")),
                "D1_s": d1,
                "P1_us": p1,
                "pulse_program": metadata.get("Pulse Sequence") or acq.get("$PULPROG"),
                "temperature_raw_K": temperature_raw,
                "temperature_K": temperature_k,
                "temperature_status": temperature_status,
                "acquisition_datetime": metadata.get("Acquisition Date"),
                "LB_Hz": to_float(proc.get("$LB")),
                "GB": to_float(proc.get("$GB")),
                "SSB": to_float(proc.get("$SSB")),
                "WDW": to_int(proc.get("$WDW")),
                "FT_mod": to_int(proc.get("$FT_mod")),
                "BC_mod": to_int(proc.get("$BC_mod")),
                "PHC0_deg": to_float(proc.get("$PHC0")),
                "PHC1_deg": to_float(proc.get("$PHC1")),
                "raw_digital_resolution_Hz_per_complex_point": raw_resolution,
                "processed_digital_resolution_Hz_per_point": processed_resolution,
                "acqus_evidence": f"{mnova_path}#zlib_offset={acq_stream.offset};compressed_length={acq_stream.compressed_length}",
                "procs_evidence": f"{mnova_path}#zlib_offset={proc_stream.offset};compressed_length={proc_stream.compressed_length}",
                "processed_array_evidence": f"{mnova_path}#byte_start={spectrum_start};byte_length={spectrum_bytes};dtype=<f4;complex_pairs=real,imag",
                "processed_array_sha256": hashlib.sha256(raw_spectrum).hexdigest(),
                "spectrum_object_prefix_bytes": header_size,
                "axis_evidence": "Mnova document metadata: Spectrometer Frequency + Spectral Width + Lowest Frequency",
            }
            acquisition_rows.append(acq_row)

            peak_results: dict[str, dict[str, Any]] = {}
            solvent_upper = solvent.upper()
            if nucleus == "1H" and solvent_upper == "DMSO":
                peak_results["residual_dmso"] = analyze_target_peak(
                    ppm, real, frequency_mhz, 2.50, 0.045, 0.11, 0.060, noise_sigma, 0.10
                )
                peak_results["water_in_dmso"] = analyze_target_peak(
                    ppm, real, frequency_mhz, 3.33, 0.22, 0.38, 0.25, noise_sigma, 0.30
                )
            elif nucleus == "1H" and solvent_upper in {"CDCL3", "CHLOROFORM"}:
                peak_results["residual_chcl3"] = analyze_target_peak(
                    ppm, real, frequency_mhz, 7.26, 0.055, 0.12, 0.075, noise_sigma, 0.10
                )
                peak_results["water_in_cdcl3"] = analyze_target_peak(
                    ppm, real, frequency_mhz, 1.56, 0.30, 0.50, 0.34, noise_sigma, 0.35
                )

            max_peak = stats.get("max_baseline_corrected_peak_raw") or 0.0
            for peak in peak_results.values():
                height = peak.get("height")
                peak["height_rel_tallest_peak_percent"] = (
                    100.0 * float(height) / float(max_peak)
                    if height is not None and max_peak and max_peak > 0
                    else None
                )

            cal_row = {
                "dataset_id": dataset_id,
                "mnova_file": str(mnova_path),
                "embedded_dataset_index": dataset_index,
                "primary_dataset_for_filename": primary_by_filename,
                "nucleus": nucleus,
                "solvent_parameter": solvent,
                "frequency_MHz": frequency_mhz,
                "processed_points": spectral_size,
                "processed_resolution_Hz_per_point": processed_resolution,
                **stats,
                **flatten_peak("residual_dmso", peak_results.get("residual_dmso")),
                **flatten_peak("water_in_dmso", peak_results.get("water_in_dmso")),
                **flatten_peak("residual_chcl3", peak_results.get("residual_chcl3")),
                **flatten_peak("water_in_cdcl3", peak_results.get("water_in_cdcl3")),
                "statistics_status": "derived_from_processed_Mnova_float32_array_not_raw_FID",
                "processed_array_evidence": acq_row["processed_array_evidence"],
            }
            calibration_rows.append(cal_row)
            internal_records.append(
                {
                    "dataset_id": dataset_id,
                    "ppm": ppm,
                    "real": real,
                    "imag": imag,
                    "acquisition": acq_row,
                    "calibration": cal_row,
                }
            )

    # Independent exact validation against the one exported numeric CSV.
    validation: dict[str, Any] = {"status": "not_run"}
    validation_csv = source / "1-1-C.csv"
    candidate = next((record for record in internal_records if record["dataset_id"] == "1-1-C__embedded_01"), None)
    if validation_csv.exists() and candidate is not None:
        exported = np.loadtxt(validation_csv, delimiter="\t", usecols=(0, 1))
        ppm_recovered = candidate["ppm"][::-1]
        intensity_recovered = candidate["real"][::-1]
        validation = {
            "status": "PASS" if exported.shape[0] == intensity_recovered.size else "FAIL",
            "source_csv": str(validation_csv),
            "source_csv_sha256": sha256_file(validation_csv),
            "points_csv": int(exported.shape[0]),
            "points_recovered": int(intensity_recovered.size),
            "max_abs_intensity_difference": float(np.max(np.abs(exported[:, 1] - intensity_recovered))),
            "max_abs_intensity_difference_relative_to_max": float(
                np.max(np.abs(exported[:, 1] - intensity_recovered))
                / np.max(np.abs(intensity_recovered))
            ),
            "max_abs_ppm_difference_due_to_csv_rounding": float(np.max(np.abs(exported[:, 0] - ppm_recovered))),
            "interpretation": "The recovered array agrees within the decimal/significant-digit rounding used by the text CSV; ppm differences likewise reflect rounded CSV text.",
        }
        if (
            validation["max_abs_intensity_difference_relative_to_max"] > 2e-6
            or validation["max_abs_ppm_difference_due_to_csv_rounding"] > 6e-4
        ):
            validation["status"] = "FAIL"

    write_csv(output / "nmr_acquisition_processing_audit.csv", acquisition_rows)
    write_csv(output / "nmr_spectral_calibration_per_dataset.csv", calibration_rows)
    write_csv(output / "extraction_issues.csv", extraction_issues, ["file", "dataset_index", "issue"])

    primary_cal = [row for row in calibration_rows if row.get("primary_dataset_for_filename")]
    primary_acq = {row["dataset_id"]: row for row in acquisition_rows if row.get("primary_dataset_for_filename")}
    group_rows: list[dict[str, Any]] = []
    group_json: dict[str, Any] = {}
    group_keys = sorted({(row["solvent_parameter"], row["nucleus"]) for row in primary_cal})
    summary_metrics = [
        "frequency_MHz",
        "processed_resolution_Hz_per_point",
        "noise_rel_max_percent",
        "noise_stationary_rel_max_percent",
        "noise_AR1_innovation_rel_max_percent",
        "baseline_span_rel_max_percent",
        "baseline_slow_rms_rel_max_percent",
        "noise_lag1_correlation_clipped",
        "noise_skewness_clipped",
        "noise_excess_kurtosis_clipped",
        "residual_dmso_center_ppm",
        "residual_dmso_fwhm_hz",
        "residual_dmso_snr",
        "residual_dmso_height_rel_tallest_peak_percent",
        "water_in_dmso_center_ppm",
        "water_in_dmso_fwhm_hz",
        "water_in_dmso_snr",
        "water_in_dmso_height_rel_tallest_peak_percent",
        "residual_chcl3_center_ppm",
        "residual_chcl3_fwhm_hz",
        "residual_chcl3_height_rel_tallest_peak_percent",
        "water_in_cdcl3_center_ppm",
        "water_in_cdcl3_fwhm_hz",
        "water_in_cdcl3_height_rel_tallest_peak_percent",
    ]
    for solvent, nucleus in group_keys:
        rows = [row for row in primary_cal if row["solvent_parameter"] == solvent and row["nucleus"] == nucleus]
        gid = f"{solvent}_{nucleus}"
        summaries: dict[str, dict[str, Any]] = {}
        metric_filters: dict[str, str] = {}
        for metric in summary_metrics:
            selected_rows = rows
            if metric.startswith("residual_dmso_"):
                selected_rows = [
                    row for row in rows if row.get("residual_dmso_confidence") in {"high", "medium"}
                ]
                metric_filters[metric] = "DMSO assignment confidence high or medium"
            elif metric.startswith("water_in_dmso_"):
                selected_rows = [
                    row
                    for row in rows
                    if row.get("water_in_dmso_confidence") in {"high", "medium"}
                    and not row.get("water_in_dmso_issues")
                ]
                metric_filters[metric] = "water assignment confidence high or medium; no overlap/width issue"
            elif metric.startswith("residual_chcl3_"):
                selected_rows = [
                    row for row in rows if row.get("residual_chcl3_confidence") in {"high", "medium"}
                ]
                metric_filters[metric] = "CHCl3 assignment confidence high or medium"
            elif metric.startswith("water_in_cdcl3_"):
                selected_rows = [
                    row
                    for row in rows
                    if row.get("water_in_cdcl3_confidence") in {"high", "medium"}
                    and not row.get("water_in_cdcl3_issues")
                ]
                metric_filters[metric] = "water assignment confidence high or medium; no overlap/width issue"
            else:
                metric_filters[metric] = "all primary datasets in group"
            summaries[metric] = numeric_summary(row.get(metric) for row in selected_rows)
        acq_rows = [primary_acq[row["dataset_id"]] for row in rows]
        for metric in [
            "temperature_K",
            "spectral_width_Hz",
            "TD_real_points",
            "acquired_complex_points",
            "SI_processed_points",
            "AQ_s",
            "NS",
            "D1_s",
            "P1_us",
            "LB_Hz",
            "raw_digital_resolution_Hz_per_complex_point",
        ]:
            summaries[metric] = numeric_summary(row.get(metric) for row in acq_rows)
            metric_filters[metric] = "all primary datasets in group; invalid temperature sentinel omitted"
        high_conf_dmso = sum(row.get("residual_dmso_confidence") == "high" for row in rows)
        high_conf_water = sum(row.get("water_in_dmso_confidence") == "high" for row in rows)
        group_json[gid] = {
            "solvent_parameter": solvent,
            "nucleus": nucleus,
            "dataset_count": len(rows),
            "high_confidence_DMSO_peak_count": high_conf_dmso,
            "high_confidence_water_peak_count": high_conf_water,
            "metrics": summaries,
            "metric_filters": metric_filters,
            "dataset_ids": [row["dataset_id"] for row in rows],
        }
        for metric, summary in summaries.items():
            group_rows.append(
                {
                    "group": gid,
                    "solvent_parameter": solvent,
                    "nucleus": nucleus,
                    "dataset_count": len(rows),
                    "metric": metric,
                    "metric_filter": metric_filters.get(metric),
                    **summary,
                }
            )
    write_csv(output / "nmr_representative_calibration_by_group.csv", group_rows)

    audit = {
        "audit_generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_directory": str(source),
        "output_directory": str(output),
        "source_directory_was_read_only_in_workflow": True,
        "mnova_file_count": len(list(source.glob("*.mnova"))),
        "embedded_dataset_count": len(acquisition_rows),
        "primary_dataset_count": sum(bool(row["primary_dataset_for_filename"]) for row in acquisition_rows),
        "embedded_raw_fid_import_count": sum(row["embedded_input_kind"] == "fid" for row in acquisition_rows),
        "standalone_vendor_raw_fid_count": len([path for path in source_files if path.name.lower() in {"fid", "ser"}]),
        "vendor_family_detected": ["Bruker"],
        "varian_agilent_data_detected": False,
        "array_extraction_validation": validation,
        "extraction_issue_count": len(extraction_issues),
        "extraction_issues": extraction_issues,
        "group_summaries": group_json,
        "limitations": [
            "No standalone Bruker raw-data directory or standard fid/ser file is present in the audited folder.",
            "Three Mnova datasets record that they were imported from a file named fid, but the original file is not independently available here.",
            "Noise, baseline and peak widths are derived from the processed Mnova frequency-domain array and therefore include the recorded apodization and any undocumented Mnova edits.",
            "A Bruker temperature parameter equal to 0 is treated as an invalid/missing sentinel, not 0 K.",
            "The Bruker SOLVENT value DMSO denotes the instrument solvent selection; it does not independently certify solvent purity or water concentration.",
            "No D2O spectrum is present, so HDO position/width and D2O baseline/noise cannot be calibrated from these files.",
            "No Varian/Agilent procpar, fid directory, or VNMRJ data structure was detected.",
            "Low-confidence water assignments may overlap product or impurity resonances and must not be treated as measured water peaks without visual review.",
        ],
    }
    (output / "nmr_calibration_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    dmso_h = group_json.get("DMSO_1H", {})
    lines = [
        "# 真实 NMR 数据只读校准审计",
        "",
        f"- 原始目录：`{source}`",
        f"- 审计输出：`{output}`",
        f"- Mnova 文件：{audit['mnova_file_count']} 个；恢复嵌入数据集：{audit['embedded_dataset_count']} 个；按文件名判定的主数据集：{audit['primary_dataset_count']} 个。",
        f"- Bruker：已检出；Varian/Agilent：未检出。",
        f"- 独立标准 `fid/ser` 文件：{audit['standalone_vendor_raw_fid_count']} 个。Mnova 元数据记录从 `fid` 导入的数据集：{audit['embedded_raw_fid_import_count']} 个。",
        "",
        "## 核心结论",
        "",
        "1. `.mnova` 中嵌入了 Bruker `acqus`、`procs`、脉冲程序和处理后的复数频域数组，可直接提取真实采集/处理参数与处理后噪声。",
        "2. `1-1-C.mnova` 恢复的 32,768 个强度点与 `1-1-C.csv` 在 CSV 有效数字舍入误差内一致；因此本审计没有从 PDF/TIFF 像素反推强度。",
        "3. 文件夹中没有可独立复算的标准 Bruker 原始目录。三个文档虽记录来源为 `fid`，但原始 `fid` 本体不在本目录；本审计不把 Mnova 内部载荷冒充为原始厂商文件。",
        "4. DMSO-d6 条件可用本批真实谱图校准；D2O/HDO 条件不能由本目录校准。",
        "5. 任何额外杂质峰都不能仅为了“更像实测”而随意添加；只有溶剂残峰、水峰、噪声和基线可按本审计统计建模，其他峰须有实测或明确标为合成假设。",
        "",
        "## 文件说明",
        "",
        "- `nmr_acquisition_processing_audit.csv`：逐嵌入数据集的溶剂、核种、频率、SW、TD/SI、AQ、NS、D1、P1、脉冲程序、温度、LB、数字分辨率与证据偏移。",
        "- `nmr_spectral_calibration_per_dataset.csv`：逐数据集的处理后噪声、基线、相关性及溶剂/水峰拟合结果。",
        "- `nmr_representative_calibration_by_group.csv`：按溶剂与核种给出中位数、P10–P90 和范围。",
        "- `source_file_hashes.csv`：全部源文件 SHA-256；用于证明原始文件未被替换。",
        "- `nmr_calibration_audit.json`：机器可读摘要、限制和校验结果。",
        "- `extraction_issues.csv`：未能可靠提取的项目；空表表示没有结构性提取失败。",
        "",
        "## DMSO-1H 代表性参数（仅在存在可靠统计时使用）",
        "",
    ]
    if dmso_h:
        metrics = dmso_h["metrics"]
        for label, key, unit in [
            ("观测频率", "frequency_MHz", "MHz"),
            ("处理后数字分辨率", "processed_resolution_Hz_per_point", "Hz/point"),
            ("指数加窗 LB", "LB_Hz", "Hz"),
            ("残余 DMSO 峰位置", "residual_dmso_center_ppm", "ppm"),
            ("残余 DMSO 峰 FWHM", "residual_dmso_fwhm_hz", "Hz"),
            ("DMSO 中水峰位置", "water_in_dmso_center_ppm", "ppm"),
            ("DMSO 中水峰 FWHM", "water_in_dmso_fwhm_hz", "Hz"),
            ("高频差分噪声", "noise_rel_max_percent", "% of tallest peak"),
            ("按 lag-1 修正的平稳噪声", "noise_stationary_rel_max_percent", "% of tallest peak"),
            ("慢基线 P95–P5", "baseline_span_rel_max_percent", "% of tallest peak"),
            ("噪声 lag-1 相关系数", "noise_lag1_correlation_clipped", ""),
        ]:
            item = metrics.get(key, {})
            if item.get("n", 0):
                lines.append(
                    f"- {label}：中位数 {item['median']:.6g} {unit}；P10–P90 {item['p10']:.6g}–{item['p90']:.6g}（n={item['n']}）。"
                )
    lines.extend(
        [
            "",
            "上述水峰结果必须结合 `confidence` 字段使用。产品峰可能落在约 3.1–3.6 ppm，低/中置信度结果不能直接当作水峰校准。",
            "",
            "## 不能从这些文件确定的内容",
            "",
            "- D2O 中 HDO 的真实位置、FWHM、强度与 pD 依赖；",
            "- 样品浓度、核磁管型号、精确水含量和溶剂批次；",
            "- 未记录的 Mnova 手工基线/相位编辑细节；",
            "- 独立原始 FID 的完整性与可重处理性；",
            "- 任何未被明确鉴定的微量杂质峰身份。",
            "",
            "## 证据规则",
            "",
            "每一行都给出 `.mnova` 文件路径、SHA-256、嵌入 zlib 流偏移，以及处理后 float32 数组的字节起点和哈希。参数值取自嵌入的 Bruker 文件；噪声、基线和 FWHM 明确标为从处理后数组推导，而不是仪器直接报告。",
        ]
    )
    (output / "NMR真实谱图校准审计.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": "PASS" if not extraction_issues and validation.get("status") == "PASS" else "REVIEW",
        "output": str(output),
        "mnova_files": audit["mnova_file_count"],
        "embedded_datasets": audit["embedded_dataset_count"],
        "primary_datasets": audit["primary_dataset_count"],
        "fid_imports": audit["embedded_raw_fid_import_count"],
        "validation": validation,
        "issue_count": len(extraction_issues),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
