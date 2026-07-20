#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: bash hpc/run_direct_24h.sh /absolute/path/to/orca [run_directory]" >&2
  exit 2
fi
ORCA_EXE="$1"
RUN_DIR="${2:-}"
PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PACKAGE_ROOT}"
source .venv/bin/activate

if [[ -z "${RUN_DIR}" ]]; then
  RUN_DIR="${PACKAGE_ROOT}/results/direct_$(date -u +%Y%m%dT%H%M%SZ)"
fi
: "${NVME_SCRATCH_ROOT:?Set NVME_SCRATCH_ROOT to a node-local NVMe directory with at least 300 GiB free}"
SCRATCH_ROOT="${NVME_SCRATCH_ROOT%/}/orca52_${USER}_$(date +%s)"
WORKFLOW_START_EPOCH="$(date +%s)"
mkdir -p "${RUN_DIR}" "${SCRATCH_ROOT}"

python scripts/preflight.py --orca "${ORCA_EXE}" --scratch "${SCRATCH_ROOT}" --production --report "${RUN_DIR}/preflight.json"
python scripts/run_24h_pipeline.py \
  --orca "${ORCA_EXE}" \
  --scratch "${SCRATCH_ROOT}" \
  --run-dir "${RUN_DIR}" \
  --core-cap 92 \
  --budget-start-epoch "${WORKFLOW_START_EPOCH}"
