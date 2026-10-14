import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from scripts import check_packed_mpi

class PackedMPIProbeTests(unittest.TestCase):
    def probe(self,stdout,code=0):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            with patch.object(check_packed_mpi,'ROOT',root),patch.object(check_packed_mpi.socket,'gethostname',return_value='node-a'), \
                 patch.dict(os.environ,{'SLURM_JOB_ID':'123'}), \
                 patch.object(check_packed_mpi.subprocess,'run',return_value=subprocess.CompletedProcess([],code,stdout,'')) as run:
                if code==0 and stdout=='node-a\n'*8:check_packed_mpi.main()
                else:
                    with self.assertRaises(RuntimeError):check_packed_mpi.main()
                env=run.call_args.kwargs['env']
                self.assertEqual(env['APPTAINERENV_OMPI_MCA_plm'],'isolated')
                self.assertEqual(env['APPTAINERENV_OMPI_MCA_ras'],'^slurm')
                return json.loads(next((root/'operations/mpi_probes').glob('*.json')).read_text())

    def test_eight_local_ranks_pass(self):self.assertEqual(self.probe('node-a\n'*8)['status'],'PASS')
    def test_other_node_is_rejected(self):self.assertEqual(self.probe('node-b\n'*8)['status'],'FAIL')
    def test_empty_output_with_zero_exit_is_rejected(self):self.assertEqual(self.probe('')['status'],'FAIL')
