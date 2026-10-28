#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
PROFILE="${1:-core}"
SCHEDULER="${2:-auto}"
bash scripts/bootstrap_python.sh
PY="$ROOT/.venv/bin/python"
export PATH="$ROOT/.venv/bin:$PATH"
if "$PY" -c 'import json; c=json.load(open("config.json")); raise SystemExit(0 if c.get("software",{}).get("feff_command") else 1)'; then
  echo 'Using configured existing FEFF command.'
else
  bash scripts/install_feff10.sh
fi
"$PY" scripts/workflow.py doctor --config config.json
"$PY" scripts/postprocess.py --self-test
exec "$PY" scripts/workflow.py execute --config config.json --profile "$PROFILE" --scheduler "$SCHEDULER"
