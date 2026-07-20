#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 STAGE MANIFEST [submit_stage.py options]" >&2
  exit 64
fi

STAGE="$1"
MANIFEST="$2"
shift 2
PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
SCHEDULER="${HPC_SCHEDULER:-slurm}"

exec python3 "$PACKAGE_ROOT/hpc/bin/submit_stage.py" \
  --package-root "$PACKAGE_ROOT" \
  --scheduler "$SCHEDULER" \
  --stage "$STAGE" \
  --manifest "$MANIFEST" \
  "$@"
