"""Write a reproducible SHA-256 manifest for package source/input files."""

from __future__ import annotations

import csv
import hashlib
import re
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
OUTPUT = PROJECT / "package_metadata/file_manifest_sha256.csv"
EXCLUDED_PARTS = {"__pycache__", ".pytest_cache", "state", "results", "logs"}
RUNTIME_SUFFIXES = {
    ".out",
    ".gbw",
    ".hess",
    ".tmp",
    ".pyc",
    ".bibtex",
    ".cpcm",
    ".cpcm_corr",
    ".densities",
    ".densitiesinfo",
    ".loc",
    ".ges",
    ".hostnames",
    ".property.txt",
    ".vpt2",
    ".allxyz",
}


def is_runtime_product(path: Path) -> bool:
    """Return True for ORCA/Python runtime products, not source fixtures."""

    name = path.name
    return (
        any(name.endswith(suffix) for suffix in RUNTIME_SUFFIXES)
        or bool(re.search(r"\.bas[0-9]+$", name))
        or ".mp2resp_" in name
        or ".autociresp_" in name
        or name.endswith(".cpscfdata.tmp.0")
        or name.endswith(".propint.tmp.0")
        or name.endswith(".autociinp.tmp")
        or name.endswith(".SHARKINP.tmp")
        or ".SHARK." in name
        or name.endswith(".shark_grid.tmp")
    )


def main() -> int:
    paths = [
        path
        for path in PROJECT.rglob("*")
        if path.is_file()
        and path != OUTPUT
        and not EXCLUDED_PARTS.intersection(path.relative_to(PROJECT).parts)
        and not is_runtime_product(path)
    ]
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUT.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("relative_path", "size_bytes", "sha256"))
        for path in sorted(paths, key=lambda value: value.relative_to(PROJECT).as_posix()):
            digest = hashlib.sha256(path.read_bytes()).hexdigest().upper()
            writer.writerow((path.relative_to(PROJECT).as_posix(), path.stat().st_size, digest))
    temporary.replace(OUTPUT)
    print(f"Wrote {len(paths)} checksums to {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
