#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: bash hpc/submit_24h.sh /absolute/path/to/orca [extra sbatch options]" >&2
  echo "Example: bash hpc/submit_24h.sh /opt/orca_6_1_1/orca --partition=cpu" >&2
  exit 2
fi

ORCA_EXE="$1"
shift
if [[ "${ORCA_EXE}" != /* ]] || [[ ! -x "${ORCA_EXE}" ]]; then
  echo "ORCA path must be absolute and executable: ${ORCA_EXE}" >&2
  exit 2
fi

PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PACKAGE_ROOT}"
mkdir -p logs results

sbatch "$@" --export="ALL,ORCA_EXE=${ORCA_EXE},PACKAGE_ROOT=${PACKAGE_ROOT}" hpc/run_24h.slurm
