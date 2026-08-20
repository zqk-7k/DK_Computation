#!/usr/bin/env python3
"""Write or verify deterministic package SHA-256 checksums."""

from __future__ import annotations

import argparse
from pathlib import Path

from package_io import atomic_write_text, iter_package_files, sha256_file


DEFAULT_OUTPUT = "package_metadata/SHA256SUMS"
DEFAULT_EXCLUDES = ("package_metadata/SHA256SUMS", "package_metadata/audit_report.json", "package_metadata/last_adaptive_goat_plan.json")


def build(root: Path, output_rel: str = DEFAULT_OUTPUT) -> str:
    excludes = set(DEFAULT_EXCLUDES) | {output_rel.replace("\\", "/")}
    lines = [f"{sha256_file(path)}  {rel}" for path, rel in iter_package_files(root, excludes)]
    return "\n".join(lines) + "\n"


def verify(root: Path, checksum_file: Path) -> list[str]:
    errors = []
    for number, line in enumerate(checksum_file.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip(): continue
        try: digest, rel = line.split("  ", 1)
        except ValueError:
            errors.append(f"line {number}: malformed"); continue
        path = (root / rel).resolve()
        try: path.relative_to(root.resolve())
        except ValueError:
            errors.append(f"line {number}: path escapes package"); continue
        if not path.is_file(): errors.append(f"missing: {rel}")
        elif sha256_file(path) != digest: errors.append(f"hash mismatch: {rel}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    root = args.package_root.resolve(); output = (root / args.output).resolve()
    if args.verify:
        errors = verify(root, output)
        for item in errors: print(item)
        return 1 if errors else 0
    atomic_write_text(output, build(root, args.output))
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
