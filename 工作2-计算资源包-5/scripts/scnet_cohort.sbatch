#!/usr/bin/env bash
# One finite, disjoint addition to the existing project allocation.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=114000M
#SBATCH --signal=B:USR1@3600
set -euo pipefail
cd "${WORK2_PACKAGE_ROOT:?}"
source scripts/scnet_environment.sh
export WORK2_STAGED_RELEASE=1
export WORK2_PRIORITY_POLICY=remaining_work
"$WORK2_PYTHON" -c 'import os; from clusterflow.common import sha; assert sha(os.environ["WORK2_COHORT_SCOPE"]) == os.environ["WORK2_COHORT_SCOPE_SHA256"], "Scope changed after submission"'
srun --nodes="$SLURM_NNODES" --ntasks="$SLURM_NNODES" --ntasks-per-node=1 \
  "$WORK2_PYTHON" scripts/check_packed_mpi.py
trap 'mkdir -p operations/packed; touch "operations/packed/DRAIN_${SLURM_JOB_ID}"' USR1
srun --nodes="$SLURM_NNODES" --ntasks="$SLURM_NNODES" --ntasks-per-node=1 \
  "$WORK2_PYTHON" scripts/packed_worker.py --execute --scope "$WORK2_COHORT_SCOPE" &
driver=$!
set +e
while true; do
  wait "$driver"
  result=$?
  if kill -0 "$driver" 2>/dev/null; then continue; fi
  exit "$result"
done
