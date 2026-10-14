#!/usr/bin/env bash
# Source inside a Slurm compute allocation, not on a login node.
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo 'SCNET compute allocation required' >&2
  return 2 2>/dev/null || exit 2
fi
module load compiler/intel/oneapi/2021.3.0 || return 2
source /public/home/qikai/apprepo/orca/6.1.0-glibc229/scripts/env.sh || return 2
export XTB_ROOT=/public/home/qikai/apprepo/xtb/v6.7.0-intelmpi2021/app
export LD_LIBRARY_PATH="$XTB_ROOT/lib:$XTB_ROOT/lib64:${LD_LIBRARY_PATH:-}"
export PATH="$XTB_ROOT/bin:$PATH"
export PROGLIST="$XTB_ROOT/bin/xtb"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
export APPTAINERENV_OMPI_MCA_pml=ob1 APPTAINERENV_OMPI_MCA_btl=self,vader,tcp
export SINGULARITYENV_OMPI_MCA_pml=ob1 SINGULARITYENV_OMPI_MCA_btl=self,vader,tcp
