#!/usr/bin/env bash
#PBS -N pom_xas
#PBS -J 0-3%4
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -V
#PBS -S /bin/bash
set -euo pipefail
# OpenPBS/PBS Pro example. For Torque, use -t 0-3%4 and site-appropriate
# nodes/ppn resource syntax. Add real queue/account and module commands.
cd "${PBS_O_WORKDIR:?Submit with qsub from the package root}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
TASK_INDEX="${PBS_ARRAY_INDEX:-${PBS_ARRAYID:?No PBS array index}}"
"${XAS_PYTHON:-python3}" scripts/run_reference_job.py --config config_reference.json --task-index "$TASK_INDEX"
