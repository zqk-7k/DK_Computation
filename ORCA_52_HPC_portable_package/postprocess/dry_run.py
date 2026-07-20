"""Dependency and path smoke test for the portable HPC post-processing layer.

This script never executes ORCA and never reads a completed property file.  It
only validates imports, adjacent package resources, the CLI parser and relative
paths in the shipped example configuration.
"""

from __future__ import annotations

import ast
import importlib.metadata
import json
import sys
from pathlib import Path


POSTPROCESS = Path(__file__).resolve().parent
PACKAGE_ROOT = POSTPROCESS.parent
sys.dont_write_bytecode = True
if str(POSTPROCESS) not in sys.path:
    sys.path.insert(0, str(POSTPROCESS))

import matplotlib
import numpy
import rdkit

import production_spectra


def inside_package(path: Path) -> bool:
    try:
        path.resolve().relative_to(PACKAGE_ROOT.resolve())
        return True
    except ValueError:
        return False


def main() -> int:
    example_path = POSTPROCESS / "examples" / "benzamide_demo_config.json"
    example = json.loads(example_path.read_text(encoding="utf-8"))
    path_fields = ["atom_map", "shielding_output", "j_output", "output_root"]
    resolved = {
        field: (example_path.parent / example[field]).resolve() for field in path_fields
    }
    problems: list[str] = []
    for field in path_fields:
        raw = str(example[field])
        if Path(raw).is_absolute() or "\\" in raw:
            problems.append(f"{field} is not a portable POSIX-style relative path: {raw}")
        if not inside_package(resolved[field]):
            problems.append(f"{field} resolves outside the unpacked HPC package: {resolved[field]}")

    required_resources = [
        PACKAGE_ROOT / "scripts" / "nmr_acquisition_appearance.py",
        PACKAGE_ROOT / "scripts" / "generate_formula_hrms.py",
        production_spectra.CALIBRATION,
    ]
    missing = [str(path) for path in required_resources if not path.exists()]
    if missing:
        problems.extend(f"missing resource: {path}" for path in missing)

    for source in (POSTPROCESS / "orca_io.py", POSTPROCESS / "production_spectra.py"):
        ast.parse(source.read_text(encoding="utf-8"), filename=str(source))
    parser = production_spectra.build_parser()
    if not parser:
        problems.append("CLI parser did not initialize")

    report = {
        "status": "PASS" if not problems else "FAIL",
        "scope": "POSTPROCESS_DEPENDENCY_AND_RELATIVE_PATH_DRY_RUN_ONLY",
        "orca_started": False,
        "package_root": str(PACKAGE_ROOT),
        "python": sys.version.split()[0],
        "dependencies": {
            "numpy": numpy.__version__,
            "matplotlib": matplotlib.__version__,
            "rdkit": getattr(rdkit, "__version__", None)
            or importlib.metadata.version("rdkit"),
        },
        "calibration_directory": str(production_spectra.CALIBRATION),
        "example_paths": {
            field: {
                "raw": example[field],
                "resolved": str(resolved[field]),
                "exists_now": resolved[field].exists(),
                "note": (
                    "The completed-result input is expected to be absent before the HPC jobs finish."
                    if field != "output_root"
                    else "The output directory is created only when post-processing runs."
                ),
            }
            for field in path_fields
        },
        "problems": problems,
    }
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())
