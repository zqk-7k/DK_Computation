#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
bash scripts/bootstrap_python.sh
export PATH="$ROOT/.venv/bin:$PATH"
python scripts/supplement.py doctor --config config.json
bash scripts/install_qe.sh
python scripts/install_pseudos.py
if command -v ifort >/dev/null 2>&1; then export FEFF_FC="${FEFF_FC:-ifort}"; fi
bash scripts/install_feff10.sh
python scripts/supplement.py execute --config config.json --scheduler "${1:-auto}"
