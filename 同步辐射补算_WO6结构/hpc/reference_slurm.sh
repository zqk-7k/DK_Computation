#!/usr/bin/env bash
#SBATCH --job-name=pom_xas
#SBATCH --array=0-3%4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=72:00:00
#SBATCH --output=pom_xas_%A_%a.out
#SBATCH --error=pom_xas_%A_%a.err
set -euo pipefail
# Submit while standing in the package root. Add actual site partition/account
# and module commands if required; none are guessed here.
cd "${SLURM_SUBMIT_DIR:?Submit with sbatch from the package root}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"${XAS_PYTHON:-python3}" scripts/run_reference_job.py --config config_reference.json --task-index "${SLURM_ARRAY_TASK_ID}"
