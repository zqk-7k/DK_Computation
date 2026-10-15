#!/usr/bin/env bash
# Finite start-level GOAT pool. Existing healthy parent jobs finish separately.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=114000M
#SBATCH --signal=B:USR1@3600
set -euo pipefail
cd "${WORK2_PACKAGE_ROOT:?}"
source scripts/scnet_environment.sh
export WORK2_STAGED_RELEASE=1
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
