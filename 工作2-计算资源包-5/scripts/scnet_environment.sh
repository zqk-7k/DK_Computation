#!/usr/bin/env bash
# Source inside a Slurm compute allocation, not on a login node.
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo 'SCNET compute allocation required' >&2
  return 2 2>/dev/null || exit 2
fi
# sbatch inherits the login environment, including a conflicting default HPC-X.
# Purge only this job's module environment; no other job/session is changed.
module purge || return 2
module load compiler/intel/oneapi/2021.3.0 || return 2
# This site's module shell function can print Tcl errors yet return zero.
case ":${LOADEDMODULES:-}:" in
  *:compiler/intel/oneapi/2021.3.0:*) ;;
  *) echo 'Required oneAPI module was not loaded' >&2; return 2 ;;
esac
source /public/home/qikai/apprepo/orca/6.1.0-glibc229/scripts/env.sh || return 2
export XTB_ROOT=/public/home/qikai/apprepo/xtb/v6.7.0-intelmpi2021/app
export LD_LIBRARY_PATH="$XTB_ROOT/lib:$XTB_ROOT/lib64:${LD_LIBRARY_PATH:-}"
export PATH="$XTB_ROOT/bin:$PATH"
export PROGLIST="$XTB_ROOT/bin/xtb"
work2_xtb_libraries=$(ldd "$XTB_ROOT/bin/xtb" 2>&1) || {
  echo "$work2_xtb_libraries" >&2; return 2;
}
if [[ "$work2_xtb_libraries" == *'not found'* ]]; then
  echo "xTB dynamic libraries unresolved: $work2_xtb_libraries" >&2
  return 2
fi
unset work2_xtb_libraries
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
export APPTAINERENV_OMPI_MCA_pml=ob1 APPTAINERENV_OMPI_MCA_btl=self,vader,tcp
export SINGULARITYENV_OMPI_MCA_pml=ob1 SINGULARITYENV_OMPI_MCA_btl=self,vader,tcp
