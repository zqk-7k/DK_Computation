#!/usr/bin/env bash
set -euo pipefail
cd "${WORK2_PACKAGE_ROOT:?}"
source scripts/scnet_environment.sh
export WORK2_STAGED_RELEASE=1
# Resolve each allocation's actual deadline before the concurrent node launch.
"$WORK2_PYTHON" -c 'import os; from pathlib import Path; from scripts.goat_bulk_worker import allocation_deadline; print(allocation_deadline(Path(os.environ["WORK2_GOAT_BULK_BASE"])))'
trap 'touch "${WORK2_GOAT_BULK_BASE:?}/DRAIN_${SLURM_JOB_ID}"' USR1
srun --kill-on-bad-exit=0 --nodes="$SLURM_NNODES" --ntasks="$SLURM_NNODES" --ntasks-per-node=1 \
  "$WORK2_PYTHON" scripts/goat_bulk_worker.py --scope "$WORK2_GOAT_BULK_SCOPE" \
  --scope-sha256 "$WORK2_GOAT_BULK_SCOPE_SHA256" &
driver=$!
set +e
while true; do
  wait "$driver"
  result=$?
  if kill -0 "$driver" 2>/dev/null; then continue; fi
  exit "$result"
done
