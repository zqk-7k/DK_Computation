#!/usr/bin/env bash
set -euo pipefail
PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
exec python3 -m unittest -v "$PACKAGE_ROOT/hpc/tests/test_hpc_layer.py"
