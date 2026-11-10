#!/usr/bin/env python3
"""Validated FEFF EGRID settings and native-output energy-window checks.

EGRID energies are relative to the FEFF Fermi level. They are not absolute
photon energies. This module never shifts spectra, extrapolates or pads tails.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def validate_egrid(grid, *, allow_high_energy_fms=False):
    if not isinstance(grid, (list, tuple)) or len(grid) != 3:
        raise ValueError("An EGRID must be [minimum_eV, maximum_eV, step_eV]")
    low, high, step = (float(v) for v in grid)
    if not all(math.isfinite(v) for v in (low, high, step)):
        raise ValueError("EGRID values must be finite")
    if not low < 0 < high or step <= 0:
        raise ValueError("XANES EGRID must span the Fermi level, with positive step")
    count = (high-low)/step
    if abs(count-round(count)) > 1.e-7:
        raise ValueError("EGRID span must be an integer multiple of its step")
    if int(round(count))+1 > 1800:
        raise ValueError("EGRID exceeds 1800 points; FEFF nex=2000 also needs a vertical contour")
    if high > 3.80998212*6.0**2 and not allow_high_energy_fms:
        raise ValueError("EGRID upper energy exceeds the manual's approximate k=6 A^-1 FMS range; review the physics before enabling allow_high_energy_fms")
    return [low, high, step]


def egrid_for_edge(feff, edge):
    grids = feff.get("egrid_by_edge", {})
    if edge not in grids:
        return None
    return validate_egrid(grids[edge], allow_high_energy_fms=bool(feff.get("allow_high_energy_fms", False)))


def absolute_window_for_edge(feff, edge):
    windows = feff.get("absolute_windows_eV", {"W_L3": [10100.0, 10300.0]})
    if edge not in windows:
        return None
    window = windows[edge]
    if not isinstance(window, (list, tuple)) or len(window) != 2:
        raise ValueError("An absolute energy window must be [minimum_eV, maximum_eV]")
    low, high = (float(v) for v in window)
    if not math.isfinite(low) or not math.isfinite(high) or not low < high:
        raise ValueError("Invalid absolute photon-energy window")
    return [low, high]


def check_absolute_window(rows, window, *, expected_step=None):
    """Require actual xmu.dat column-1 support; return machine-readable evidence."""
    energy = [float(row[0]) for row in rows]
    if len(energy) < 2 or not all(math.isfinite(v) for v in energy):
        raise ValueError("Insufficient/nonfinite native photon-energy samples")
    gaps = [b-a for a, b in zip(energy, energy[1:])]
    if min(gaps) <= 0:
        raise ValueError("Native photon energy must increase strictly")
    low, high = (float(v) for v in window)
    tolerance = 1.e-5
    if energy[0] > low+tolerance or energy[-1] < high-tolerance:
        raise ValueError(
            f"Native xmu.dat covers {energy[0]:.6f}..{energy[-1]:.6f} eV, "
            f"not the requested {low:g}..{high:g} eV. Archive this attempt and "
            "rerun a justified wider relative EGRID in a new profile. Do not pad, "
            "extrapolate, or shift individual spectra to manufacture coverage."
        )
    if expected_step is not None and max(gaps) > float(expected_step)*1.05 + tolerance:
        raise ValueError("Native output has gaps exceeding the requested EGRID step; check grid truncation and FEFF logs")
    return {"status": "PASS", "requested_absolute_window_eV": [low, high],
            "native_absolute_window_eV": [energy[0], energy[-1]],
            "native_maximum_step_eV": max(gaps), "native_rows": len(rows),
            "energy_source": "xmu.dat column 1 (omega), without per-spectrum shift",
            "extrapolation_or_tail_padding": False}


def check_job_energy_coverage(job, rows):
    metadata_path = Path(job)/"job.json"
    if not metadata_path.is_file():
        return None
    metadata = json.loads(metadata_path.read_text(encoding="utf-8-sig"))
    window = metadata.get("requested_absolute_window_eV")
    if window is None:
        return None
    grid = metadata.get("relative_egrid_eV")
    return check_absolute_window(rows, window, expected_step=grid[2] if grid else None)


def self_test():
    assert validate_egrid([-160, 130, .5]) == [-160.0, 130.0, .5]
    for bad in [[-500, 130, .1], [-160, 160, .5], [-160, 130, 0], [1, 100, .5]]:
        try:
            validate_egrid(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"Accepted invalid grid {bad}")
    rows = [[10050+i*.5, 0, 0, 1, 1, 0] for i in range(581)]
    assert check_absolute_window(rows, [10100, 10300], expected_step=.5)["status"] == "PASS"
    try:
        check_absolute_window(rows[150:], [10100, 10300])
    except ValueError:
        pass
    else:
        raise AssertionError("Accepted missing pre-edge coverage")
    print("EGRID and native-window logic self-test PASS; no physical spectrum was calculated")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.job is None:
        parser.error("Specify --job or --self-test")
    from run_feff_job import numeric_table
    rows, _ = numeric_table(args.job/"xmu.dat", 6)
    result = check_job_energy_coverage(args.job, rows)
    if result is None:
        raise ValueError("job.json has no requested_absolute_window_eV")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
