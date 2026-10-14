#!/usr/bin/env bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=114000M
#SBATCH --signal=B:USR1@3600
#SBATCH --job-name=work2-64-pilot
set -euo pipefail
cd "${WORK2_PACKAGE_ROOT:?}"
source scripts/scnet_environment.sh
export WORK2_PYTHON="${WORK2_PYTHON:-/public/home/qikai/.conda/envs/orca-nmr-py310/bin/python}"
export WORK2_RUNTIME_EVIDENCE="$WORK2_PACKAGE_ROOT/operations/runtime_${SLURM_JOB_ID}.json"
# One canonical record for all nodes, so identical inputs retain one cache identity.
"$WORK2_PYTHON" scripts/site_check.py --execute-version-check --record "$WORK2_RUNTIME_EVIDENCE"
srun --nodes="$SLURM_NNODES" --ntasks="$SLURM_NNODES" --ntasks-per-node=1 \
  "$WORK2_PYTHON" scripts/check_packed_mpi.py
trap 'mkdir -p operations/packed; touch "operations/packed/DRAIN_${SLURM_JOB_ID}"' USR1
srun --nodes="$SLURM_NNODES" --ntasks="$SLURM_NNODES" --ntasks-per-node=1 \
  "$WORK2_PYTHON" scripts/packed_worker.py --execute --scope config/pilot_scope_20261014.json &
driver=$!
set +e
while true; do
  wait "$driver"
  result=$?
  if kill -0 "$driver" 2>/dev/null; then continue; fi
  exit "$result"
done
