#!/usr/bin/env bash
#SBATCH --job-name=orca52_nmr_hrms_24h
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=92
#SBATCH --hint=nomultithread
#SBATCH --mem=680G
#SBATCH --time=23:55:00
#SBATCH --exclusive
#SBATCH --signal=B:USR1@2100
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err

set -euo pipefail

PACKAGE_ROOT="${PACKAGE_ROOT:-${SLURM_SUBMIT_DIR:-}}"
if [[ -z "${PACKAGE_ROOT}" ]] || [[ ! -f "${PACKAGE_ROOT}/config/hardware_104physical_24h.json" ]]; then
  echo "PACKAGE_ROOT is invalid. Submit through hpc/submit_24h.sh or export its absolute path." >&2
  exit 2
fi
PACKAGE_ROOT="$(cd "${PACKAGE_ROOT}" && pwd)"
cd "${PACKAGE_ROOT}"
mkdir -p logs results

: "${ORCA_EXE:?Set ORCA_EXE to the absolute path of the ORCA 6.1.1 driver}"
if [[ "${ORCA_EXE}" != /* ]] || [[ ! -x "${ORCA_EXE}" ]]; then
  echo "ORCA_EXE must be an absolute executable path: ${ORCA_EXE}" >&2
  exit 2
fi

if [[ -n "${VENV_DIR:-}" ]]; then
  source "${VENV_DIR}/bin/activate"
elif [[ -f "${PACKAGE_ROOT}/.venv/bin/activate" ]]; then
  source "${PACKAGE_ROOT}/.venv/bin/activate"
fi

PYTHON_EXE="${PYTHON_EXE:-python3}"
WORKFLOW_START_EPOCH="$(date +%s)"
if [[ -n "${SLURM_TMPDIR:-}" ]]; then
  SCRATCH_ROOT="${SLURM_TMPDIR}/orca52_${SLURM_JOB_ID:-direct}"
elif [[ -n "${NVME_SCRATCH_ROOT:-}" ]]; then
  SCRATCH_ROOT="${NVME_SCRATCH_ROOT%/}/orca52_${SLURM_JOB_ID:-direct}"
else
  echo "Production run requires SLURM_TMPDIR or NVME_SCRATCH_ROOT on node-local NVMe; refusing /tmp fallback." >&2
  exit 2
fi
RUN_DIR="${RUN_DIR:-${PACKAGE_ROOT}/results/production_${SLURM_JOB_ID:-direct}}"
mkdir -p "${SCRATCH_ROOT}" "${RUN_DIR}"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export ORCA_USE_SLURM_STEPS=1
ulimit -s unlimited || true

"${PYTHON_EXE}" scripts/preflight.py \
  --orca "${ORCA_EXE}" \
  --scratch "${SCRATCH_ROOT}" \
  --production \
  --report "${RUN_DIR}/preflight.json"

set +e
"${PYTHON_EXE}" scripts/run_24h_pipeline.py \
  --orca "${ORCA_EXE}" \
  --scratch "${SCRATCH_ROOT}" \
  --run-dir "${RUN_DIR}" \
  --core-cap 92 \
  --budget-start-epoch "${WORKFLOW_START_EPOCH}" &
PIPELINE_PID=$!

forward_usr1() {
  echo "Forwarding USR1 to pipeline ${PIPELINE_PID}" >&2
  kill -USR1 "${PIPELINE_PID}" 2>/dev/null || true
}
forward_term() {
  echo "Forwarding TERM to pipeline ${PIPELINE_PID}" >&2
  kill -TERM "${PIPELINE_PID}" 2>/dev/null || true
}
trap forward_usr1 USR1
trap forward_term TERM INT
wait "${PIPELINE_PID}"
STATUS=$?
set -e

echo "Pipeline exit code: ${STATUS}"
echo "Result directory: ${RUN_DIR}"
exit "${STATUS}"
