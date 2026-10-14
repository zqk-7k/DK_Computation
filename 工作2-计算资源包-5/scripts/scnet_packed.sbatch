#!/usr/bin/env bash
# Supply account/partition/QOS/nodes/time at submission after live verification.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=114000M
#SBATCH --signal=B:USR1@1800
#SBATCH --job-name=work2-64-v5
set -euo pipefail
cd "${WORK2_PACKAGE_ROOT:?Set an independent production copy of 工作2-计算资源包-5}"
source scripts/scnet_environment.sh
test -r "${WORK2_RUNTIME_EVIDENCE:?Run scripts/site_check.py first}"
# Advance notice drains new claims on all nodes. Current ORCA calls continue.
trap 'mkdir -p operations/packed; touch "operations/packed/DRAIN_${SLURM_JOB_ID}"' USR1
srun --nodes="${SLURM_NNODES:?}" --ntasks="${SLURM_NNODES}" --ntasks-per-node=1 \
  /public/home/qikai/.conda/envs/orca-nmr-py310/bin/python scripts/packed_worker.py --execute &
work2_srun_pid=$!
wait "$work2_srun_pid" || work2_status=$?
# wait may be interrupted by USR1 while srun is still healthy; wait again.
if kill -0 "$work2_srun_pid" 2>/dev/null; then work2_status=0; wait "$work2_srun_pid" || work2_status=$?; fi
exit "${work2_status:-0}"
