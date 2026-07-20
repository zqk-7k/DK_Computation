#!/usr/bin/env python3
"""Generate auditable formula-level HRMS candidate-ion predictions for Figure 3.

This wrapper deliberately keeps HRMS formula work separate from ORCA.  It uses
the previously audited exact-mass/isotope convolution engine, verifies that
engine by SHA-256, includes all 52 registered structures, and writes additional
CSV/TXT provenance and validation products.

The results are CALCULATED predictions.  They are not experimental HRMS scans,
do not predict ion abundance, and contain no fragmentation or contaminants.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import platform
import shutil
import sys
import textwrap
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


PREDICTION_LABEL = "CALCULATED PREDICTION — NOT EXPERIMENTAL"
SCOPE_NOTICE = (
    "Formula-level intact precursor/adduct exact masses and representative "
    "natural-isotope envelopes only; no fragment ions, experimental noise, "
    "contaminants, ion-abundance prediction or measured peaks."
)
NIST_ISOTOPE_URL = (
    "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?all=all"
    "&ascii=ascii2&isotype=all"
)

SERIAL_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_ROOT = SERIAL_ROOT.parent
ENGINE_PATH = (
    REFERENCE_ROOT
    / "ORCA_51_NMR_HRMS_predictions"
    / "scripts"
    / "generate_formula_hrms.py"
)
EXPECTED_ENGINE_SHA256 = (
    "40dc9cf406aa39adc9723c6cbab51904835a3d1f1a9bd0cac25f9e339905cf38"
)
DEFAULT_SDF = REFERENCE_ROOT / "Figure3_52_ORCA_ordered_reference_states.sdf"
DEFAULT_REGISTRY = REFERENCE_ROOT / "Figure3_52_ORCA_reference_state_registry.csv"
DEFAULT_OUTPUT = (
    SERIAL_ROOT / "results" / "HRMS_formula_predictions_CALCULATED_NOT_EXPERIMENTAL"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def ensure_safe_output(path: Path) -> None:
    resolved = path.resolve()
    allowed = (SERIAL_ROOT / "results").resolve()
    try:
        resolved.relative_to(allowed)
    except ValueError as exc:
        raise ValueError(f"Output must remain below {allowed}: {resolved}") from exc
    if "CALCULATED_NOT_EXPERIMENTAL" not in resolved.name:
        raise ValueError(
            "Output directory name must contain CALCULATED_NOT_EXPERIMENTAL"
        )


def load_locked_engine():
    if not ENGINE_PATH.is_file():
        raise FileNotFoundError(f"Pinned HRMS engine is missing: {ENGINE_PATH}")
    observed = sha256_file(ENGINE_PATH).lower()
    if observed != EXPECTED_ENGINE_SHA256.lower():
        raise RuntimeError(
            "Pinned HRMS engine hash changed; audit the engine before use. "
            f"Expected {EXPECTED_ENGINE_SHA256}, observed {observed}."
        )
    specification = importlib.util.spec_from_file_location(
        "locked_formula_hrms_engine", ENGINE_PATH
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Cannot import pinned HRMS engine: {ENGINE_PATH}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)

    # Correct the old display label and include benzamide in this 52-compound run.
    module.PREDICTION_LABEL = PREDICTION_LABEL
    module.NO_FRAGMENT_NOTICE = SCOPE_NOTICE
    module.SKIP_IDS = frozenset()

    # Br and K are not present in the current 52 formulas, but are included so
    # the workflow fails safely only outside this declared elemental coverage.
    module.ISOTOPES["Br"] = (
        module.Isotope(79, 78.91833760000, 0.5069),
        module.Isotope(81, 80.91628970000, 0.4931),
    )
    module.ISOTOPES["K"] = (
        module.Isotope(39, 38.96370648640, 0.932581),
        module.Isotope(40, 39.96399816600, 0.000117),
        module.Isotope(41, 40.96182525790, 0.067302),
    )
    install_unclipped_plotter(module)
    return module


def install_unclipped_plotter(engine) -> None:
    """Replace the legacy plot layout with a long-title-safe implementation."""

    def plot_profile(
        path,
        figure_id,
        display_name,
        ion,
        envelope,
        x,
        y,
        resolving_power,
    ) -> None:
        engine.plt.rcParams.update(
            {
                "font.family": "Arial",
                "font.size": 8.5,
                "axes.linewidth": 0.8,
                "xtick.major.width": 0.8,
                "ytick.major.width": 0.8,
            }
        )
        figure, axis = engine.plt.subplots(figsize=(7.2, 4.6))
        axis.plot(
            x,
            y,
            color="#202020",
            linewidth=1.0,
            label=f"Gaussian profile, R={resolving_power:g}",
        )
        axis.vlines(
            [stick.mz for stick in envelope.sticks],
            0,
            engine.np.array([stick.relative_intensity for stick in envelope.sticks]),
            color="#2574A9",
            linewidth=0.75,
            alpha=0.80,
            label="Natural-abundance isotope sticks",
        )
        axis.set_ylim(0, 108)
        axis.set_xlabel("Calculated m/z")
        axis.set_ylabel("Calculated relative intensity (%)")
        wrapped_name = "\n".join(
            textwrap.wrap(f"{figure_id} | {display_name}", width=66)
        )
        axis.set_title(
            f"{wrapped_name}\n{ion.label} | "
            f"{engine.charged_formula(ion.formula_counts, ion.charge)}",
            fontsize=9.0,
            fontweight="bold",
            pad=7.0,
        )
        figure.suptitle(
            PREDICTION_LABEL,
            fontsize=10.2,
            fontweight="bold",
            y=0.985,
        )
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.legend(frameon=False, fontsize=7.5, loc="upper right")
        figure.text(
            0.105,
            0.035,
            "Formula isotope envelope only. No fragment ions or experimental peaks.",
            fontsize=7.2,
            color="#555555",
        )
        figure.subplots_adjust(left=0.105, right=0.985, bottom=0.16, top=0.76)
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(
            path,
            dpi=300,
            metadata={
                "Title": PREDICTION_LABEL,
                "Description": SCOPE_NOTICE,
                "Author": "Formula-level Python calculation",
            },
        )
        engine.plt.close(figure)

    engine.plot_profile = plot_profile


def read_commented_csv(path: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    metadata: dict[str, str] = {}
    data_lines: list[str] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for line in handle:
            if line.startswith("#"):
                text = line[1:].strip()
                if ":" in text:
                    key, value = text.split(":", 1)
                    metadata[key.strip()] = value.strip()
                continue
            data_lines.append(line)
    return metadata, list(csv.DictReader(data_lines))


def write_banner_csv(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, Any]],
    metadata: Mapping[str, Any],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        handle.write(f"# {PREDICTION_LABEL}\n")
        handle.write(f"# scope: {SCOPE_NOTICE}\n")
        for key, value in metadata.items():
            handle.write(f"# {key}: {value}\n")
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_tab_text(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        handle.write(f"{PREDICTION_LABEL}\n")
        handle.write(f"{SCOPE_NOTICE}\n")
        handle.write("\t".join(fieldnames) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(key, "")) for key in fieldnames) + "\n")


def load_plain_registry(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def make_element_self_tests(engine) -> dict[str, Any]:
    tests: dict[str, Any] = {}
    expected_m2 = {
        "Cl": 100.0 * 0.2424 / 0.7576,
        "Br": 100.0 * 0.4931 / 0.5069,
        "S": 100.0 * 0.0425 / 0.9499,
    }
    for element, expected in expected_m2.items():
        result = engine.calculate_envelope(
            {element: 1}, charge=+1, max_states=100, min_relative_percent=0.0
        )
        mono = min(result.sticks, key=lambda item: abs(item.mass_shift_u))
        m2_candidates = [
            item for item in result.sticks if 1.8 <= item.mass_shift_u <= 2.2
        ]
        observed = max(m2_candidates, key=lambda item: item.relative_intensity)
        tests[element] = {
            "expected_M_plus_2_to_M_percent": expected,
            "observed_M_plus_2_to_M_percent": observed.relative_intensity,
            "absolute_difference_percentage_points": abs(
                expected - observed.relative_intensity
            ),
            "pass_tolerance_1e-8": abs(expected - observed.relative_intensity) < 1e-8,
            "monoisotopic_relative_intensity_percent": mono.relative_intensity,
        }
    tests["P"] = {
        "stable_isotopes_in_model": len(engine.ISOTOPES["P"]),
        "pass_single_stable_isotope": len(engine.ISOTOPES["P"]) == 1,
    }
    tests["I"] = {
        "stable_isotopes_in_model": len(engine.ISOTOPES["I"]),
        "pass_single_stable_isotope": len(engine.ISOTOPES["I"]) == 1,
    }
    return tests


def build_derivative_outputs(
    engine,
    output: Path,
    registry_path: Path,
    sdf_path: Path,
    selected_only: set[str],
) -> dict[str, Any]:
    summary_path = output / "Figure3_formula_HRMS_prediction_summary.csv"
    _, summary = read_commented_csv(summary_path)
    registry = load_plain_registry(registry_path)
    if selected_only:
        registry = [row for row in registry if row["figure_id"] in selected_only]

    expected_ids = [row["figure_id"] for row in registry]
    observed_ids = sorted({row["figure_id"] for row in summary})
    missing_ids = sorted(set(expected_ids) - set(observed_ids))
    extra_ids = sorted(set(observed_ids) - set(expected_ids))
    if missing_ids or extra_ids:
        raise RuntimeError(
            f"HRMS summary ID mismatch: missing={missing_ids}, extra={extra_ids}"
        )

    summary_fields = list(summary[0].keys()) if summary else []
    write_tab_text(output / "ALL_candidate_ions.txt", summary_fields, summary)

    per_compound = defaultdict(list)
    for row in summary:
        per_compound[row["figure_id"]].append(row)
    registry_by_id = {row["figure_id"]: row for row in registry}

    # Independent numerical cross-check: compare the embedded-isotope reference
    # formula mass with RDKit's separately maintained exact-mass implementation.
    # Small sub-microdalton differences are expected from isotope-table rounding.
    reference_mass_by_id = {
        figure_id: float(rows[0]["reference_exact_mass_u"])
        for figure_id, rows in per_compound.items()
    }
    rdkit_mass_checks: list[dict[str, Any]] = []
    for figure_id, molecule in engine.load_ordered_sdf(sdf_path):
        if figure_id not in reference_mass_by_id:
            continue
        rdkit_mass = float(engine.rdMolDescriptors.CalcExactMolWt(molecule))
        embedded_mass = reference_mass_by_id[figure_id]
        rdkit_mass_checks.append(
            {
                "figure_id": figure_id,
                "embedded_isotope_mass_u": embedded_mass,
                "rdkit_exact_mass_u": rdkit_mass,
                "absolute_difference_u": abs(embedded_mass - rdkit_mass),
            }
        )
    max_rdkit_difference = max(
        (row["absolute_difference_u"] for row in rdkit_mass_checks), default=0.0
    )
    if max_rdkit_difference > 1.0e-5:
        raise RuntimeError(
            "Independent RDKit exact-mass cross-check exceeded 1e-5 u: "
            f"{max_rdkit_difference:.12g}"
        )
    for figure_id in expected_ids:
        row = registry_by_id[figure_id]
        lines = [
            PREDICTION_LABEL,
            SCOPE_NOTICE,
            "",
            f"Figure ID: {figure_id}",
            f"Compound: {row['display_name']}",
            f"Registered formula/reference state: {row['formula']}",
            f"InChIKey: {row['inchikey']}",
            f"Registered isomeric SMILES: {row['isomeric_smiles']}",
            "",
            "Candidate intact ions (candidate does not mean experimentally observed):",
        ]
        for ion in per_compound[figure_id]:
            lines.append(
                f"- {ion['ion_mode']} {ion['ion_label']}: "
                f"{ion['ion_formula']}, calculated m/z {ion['monoisotopic_mz']}"
            )
        lines.extend(
            [
                "",
                "The corresponding isotope-stick CSV, Gaussian profile CSV and PNG "
                "are stored in this compound directory.",
                "Experimental reporting must use a measured 'found' m/z; these files "
                "supply only the calculated value.",
            ]
        )
        (output / figure_id / f"{figure_id}_HRMS_summary.txt").write_text(
            "\n".join(lines) + "\n", encoding="utf-8-sig"
        )

    combined_sticks: list[dict[str, Any]] = []
    for figure_id in expected_ids:
        for path in sorted((output / figure_id).glob("*_isotope_sticks.csv")):
            metadata, rows = read_commented_csv(path)
            for row in rows:
                combined_sticks.append(
                    {
                        "figure_id": figure_id,
                        "display_name": registry_by_id[figure_id]["display_name"],
                        "ion_label": metadata.get("ion_label", ""),
                        "ion_mode": metadata.get("ion_mode", ""),
                        "ion_formula": metadata.get("ion_formula", ""),
                        "source_file": str(path.relative_to(output)),
                        **row,
                    }
                )
    stick_fields = [
        "figure_id",
        "display_name",
        "ion_label",
        "ion_mode",
        "ion_formula",
        "isotopologue_index",
        "calculated_mz",
        "natural_abundance_probability",
        "relative_intensity_percent",
        "mass_shift_from_monoisotopic_mz_u",
        "source_file",
    ]
    write_banner_csv(
        output / "ALL_isotope_sticks.csv",
        stick_fields,
        combined_sticks,
        {"row_count": len(combined_sticks)},
    )
    write_tab_text(output / "ALL_isotope_sticks.txt", stick_fields, combined_sticks)

    isotope_rows: list[dict[str, Any]] = []
    for element, isotopes in sorted(engine.ISOTOPES.items()):
        for isotope in isotopes:
            isotope_rows.append(
                {
                    "element": element,
                    "mass_number": isotope.mass_number,
                    "exact_isotope_mass_u": isotope.exact_mass,
                    "representative_natural_abundance_fraction": isotope.abundance,
                    "source_url": NIST_ISOTOPE_URL,
                }
            )
    write_banner_csv(
        output / "isotope_mass_and_abundance_table.csv",
        list(isotope_rows[0].keys()),
        isotope_rows,
        {
            "source_note": (
                "Embedded representative terrestrial isotope compositions; natural "
                "materials can vary within documented ranges."
            )
        },
    )

    self_tests = make_element_self_tests(engine)
    audit = {
        "prediction_status": PREDICTION_LABEL,
        "scope": SCOPE_NOTICE,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_sdf": str(sdf_path.resolve()),
        "input_sdf_sha256": sha256_file(sdf_path),
        "input_registry": str(registry_path.resolve()),
        "input_registry_sha256": sha256_file(registry_path),
        "pinned_engine": str(ENGINE_PATH.resolve()),
        "pinned_engine_sha256": sha256_file(ENGINE_PATH),
        "expected_compound_count": len(expected_ids),
        "observed_compound_count": len(observed_ids),
        "candidate_ion_count": len(summary),
        "combined_isotope_stick_count": len(combined_sticks),
        "missing_ids": missing_ids,
        "extra_ids": extra_ids,
        "fragment_ions_generated": False,
        "experimental_peaks_generated": False,
        "independent_reference_mass_crosscheck": {
            "method": "RDKit CalcExactMolWt on each registered SDF molecule",
            "comparison_count": len(rdkit_mass_checks),
            "maximum_absolute_difference_u": max_rdkit_difference,
            "tolerance_u": 1.0e-5,
            "pass": max_rdkit_difference <= 1.0e-5,
            "records": rdkit_mass_checks,
        },
        "element_self_tests": self_tests,
        "all_self_tests_pass": all(
            all(value for key, value in test.items() if key.startswith("pass_"))
            for test in self_tests.values()
        ),
        "software": {
            "python": sys.version,
            "python_executable": sys.executable,
            "platform": platform.platform(),
            "rdkit": getattr(engine.Chem.rdBase, "rdkitVersion", "unknown"),
            "numpy": engine.np.__version__,
            "matplotlib": engine.matplotlib.__version__,
        },
    }
    (output / "VALIDATION_AND_PROVENANCE.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    if not audit["all_self_tests_pass"]:
        raise RuntimeError("Element isotope self-test failed; see validation JSON")
    return audit


def write_hash_manifest(output: Path, extra_inputs: Sequence[Path]) -> None:
    rows: list[dict[str, Any]] = []
    seen: set[Path] = set()
    for path in [*extra_inputs, *sorted(output.rglob("*"))]:
        if not path.is_file() or path.name == "FILE_MANIFEST_SHA256.csv":
            continue
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        try:
            relative = str(path.relative_to(SERIAL_ROOT))
            role = "generated_output" if output in path.parents else "workflow_file"
        except ValueError:
            relative = str(resolved)
            role = "locked_input"
        rows.append(
            {
                "role": role,
                "path": relative,
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    write_banner_csv(
        output / "FILE_MANIFEST_SHA256.csv",
        ("role", "path", "size_bytes", "sha256"),
        rows,
        {"manifest_excludes_itself": True, "file_count": len(rows)},
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate all-52 formula-level HRMS candidate-ion predictions."
    )
    parser.add_argument("--sdf", type=Path, default=DEFAULT_SDF)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--only", default="", help="Optional comma-separated F3 IDs for a test run."
    )
    parser.add_argument("--resolution", type=float, default=60000.0)
    parser.add_argument("--min-relative-percent", type=float, default=0.01)
    parser.add_argument("--max-isotope-states", type=int, default=20000)
    parser.add_argument("--points-per-fwhm", type=int, default=12)
    parser.add_argument("--max-profile-points", type=int, default=250000)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    ensure_safe_output(args.output_dir)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        if not args.overwrite:
            raise RuntimeError(
                f"Output is not empty: {args.output_dir}. Use --overwrite to rebuild."
            )
        shutil.rmtree(args.output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for key in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[key] = "1"

    engine = load_locked_engine()
    engine_args = argparse.Namespace(
        sdf=args.sdf,
        registry=args.registry,
        output_dir=args.output_dir,
        only=args.only,
        resolution=args.resolution,
        min_relative_percent=args.min_relative_percent,
        max_isotope_states=args.max_isotope_states,
        points_per_fwhm=args.points_per_fwhm,
        max_profile_points=args.max_profile_points,
    )
    result = engine.run(engine_args)
    if result != 0:
        return int(result)

    selected = {item.strip() for item in args.only.split(",") if item.strip()}
    audit = build_derivative_outputs(
        engine, args.output_dir, args.registry, args.sdf, selected
    )
    write_hash_manifest(
        args.output_dir,
        (
            args.sdf,
            args.registry,
            ENGINE_PATH,
            Path(__file__).resolve(),
            SERIAL_ROOT / "run_hrms_formula.ps1",
            SERIAL_ROOT / "README_HRMS_如何调用.md",
        ),
    )
    print(
        f"{PREDICTION_LABEL}: validated {audit['observed_compound_count']} "
        f"compounds and {audit['candidate_ion_count']} candidate ions."
    )
    print(f"Output: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
