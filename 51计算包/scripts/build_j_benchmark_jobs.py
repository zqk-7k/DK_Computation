#!/usr/bin/env python3
"""Fail-closed placeholder for a 51-compound cross-validated J-method gate.

The two-compound source implementation hard-coded its pilot compounds and is therefore
not scientifically transferable to this package. This command deliberately refuses
job generation until a new representative benchmark panel is approved.
"""

from __future__ import annotations

import argparse


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--approved-panel-gate", help="reserved for the commissioned 51-compound implementation")
    parser.parse_args()
    parser.error(
        "J benchmark generation is FAIL_CLOSED: approve a representative cross-compound "
        "panel and commission a hash-bound implementation; do not reuse the old two-compound gate"
    )
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
