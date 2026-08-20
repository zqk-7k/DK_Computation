#!/usr/bin/env python3
"""Verify all files named in an ORCA runner receipt."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from package_io import sha256_file


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()
    receipt = json.loads(args.receipt.read_text(encoding="utf-8"))
    root = args.receipt.parent.resolve(); errors = []
    for rel, expected in receipt.get("files_sha256", {}).items():
        path = (root / rel).resolve()
        try: path.relative_to(root)
        except ValueError: errors.append(f"unsafe path: {rel}"); continue
        if not path.is_file(): errors.append(f"missing: {rel}")
        elif sha256_file(path) != expected: errors.append(f"hash mismatch: {rel}")
    print(json.dumps({"receipt": str(args.receipt), "valid": not errors, "errors": errors}, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
