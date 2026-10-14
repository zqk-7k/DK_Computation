"""Exercise inherited module conflict and misleading module success, without QM."""
import os
from pathlib import Path
import subprocess
import unittest

SCRIPT=Path(__file__).resolve().parents[1]/'scripts/scnet_environment.sh'
HARNESS=r'''
set -euo pipefail
export SLURM_JOB_ID=123 LOADEDMODULES=mpi/hpcx/2.11.0/gcc-7.3.1
module() {
  if [[ "$1" == purge ]]; then LOADEDMODULES=''; return 0; fi
  if [[ "$MODE" != masked_failure && -z "$LOADEDMODULES" ]]; then
    LOADEDMODULES=compiler/intel/oneapi/2021.3.0
  fi
  return 0
}
source() { echo ORCA_ENV_SOURCED; }
ldd() {
  if [[ "$MODE" == missing_library ]]; then echo 'libmkl_intel_lp64.so.1 => not found';
  else echo 'libmkl_intel_lp64.so.1 => /mock/mkl/libmkl_intel_lp64.so.1'; fi
}
. "$SCRIPT"
echo "READY:$LOADEDMODULES:$OMP_NUM_THREADS"
'''

class ScnetEnvironmentTests(unittest.TestCase):
    def run_mode(self,mode):
        return subprocess.run(['bash','-c',HARNESS],env={**os.environ,'SCRIPT':str(SCRIPT),'MODE':mode},capture_output=True,text=True)

    def test_inherited_hpcx_removed_before_oneapi(self):
        result=self.run_mode('healthy')
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('READY:compiler/intel/oneapi/2021.3.0:1',result.stdout)

    def test_zero_exit_module_failure_cannot_proceed(self):
        result=self.run_mode('masked_failure')
        self.assertNotEqual(result.returncode,0)
        self.assertNotIn('ORCA_ENV_SOURCED',result.stdout)

    def test_missing_mkl_prevents_quantum_start(self):
        result=self.run_mode('missing_library')
        self.assertNotEqual(result.returncode,0)
        self.assertIn('xTB dynamic libraries unresolved',result.stderr)
        self.assertNotIn('READY:',result.stdout)
