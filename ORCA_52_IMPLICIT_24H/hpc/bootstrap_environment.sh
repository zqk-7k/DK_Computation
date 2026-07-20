#!/usr/bin/env bash
set -euo pipefail

PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PACKAGE_ROOT}"
PYTHON_EXE="${PYTHON_EXE:-python3}"

"${PYTHON_EXE}" -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip wheel
python -m pip install -r requirements.txt
python -m pip check
python -m py_compile scripts/*.py
echo "Environment ready: ${PACKAGE_ROOT}/.venv"
