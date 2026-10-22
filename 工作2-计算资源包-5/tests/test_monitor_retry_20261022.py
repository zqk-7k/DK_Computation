import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from clusterflow.common import writej, sha
from clusterflow.observation import observe
from clusterflow.goat_retry import retry_token

class ObservationTests(unittest.TestCase):
    def test_zero_filled_observation_is_unknown_and_reported(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'heartbeat.json';p.write_bytes(b'\0'*80);errors=[]
            with patch('clusterflow.observation.time.sleep'):
                self.assertIsNone(observe(p, errors))
            self.assertEqual(errors[0]['status'],'UNREADABLE_OBSERVATION')
            self.assertEqual(p.read_bytes(),b'\0'*80)
    def test_atomic_replacement_transient_is_retried(self):
        with patch('clusterflow.observation.readj',side_effect=[ValueError('partial'),{'updated_epoch':1}]),patch('clusterflow.observation.time.sleep'):
            errors=[];self.assertEqual(observe('heartbeat',errors),{'updated_epoch':1});self.assertEqual(errors,[])
    def test_missing_optional_is_not_success(self):
        with tempfile.TemporaryDirectory() as d:
            errors=[];self.assertIsNone(observe(Path(d)/'missing',errors,missing_ok=True));self.assertEqual(errors,[])
    def test_missing_required_reported(self):
        with tempfile.TemporaryDirectory() as d,patch('clusterflow.observation.time.sleep'):
            errors=[];self.assertIsNone(observe(Path(d)/'missing',errors));self.assertEqual(len(errors),1)

class RetryTests(unittest.TestCase):
    def setUp(self):
        self.t=tempfile.TemporaryDirectory();self.addCleanup(self.t.cleanup);self.r=Path(self.t.name);self.d=self.r/'jobs/j';self.d.mkdir(parents=True)
        self.old=self.r/'old';writej(self.old/'receipt.json',{'status':'RUNNING'})
        self.proof=self.d/'recovery_history/attempt/proof.json';writej(self.proof,{'job_id':'j','scheduler_state':'NODE_FAIL','stage':str(self.old),'receipt_sha256':sha(self.old/'receipt.json')})
        self.auth={'status':'PREPARED_VERIFIED_NODE_FAIL','job_id':'j','manifest_sha256':'manifest','inputs':{'goat.inp':'input'},'attempt_token':'a'*32,'proof_path':str(self.proof.relative_to(self.r)),'proof_sha256':sha(self.proof)}
        writej(self.d/'goat_retry_authorization.json',self.auth)
    def get(self):return retry_token(self.r,self.d,'j','manifest',{'goat.inp':'input'})
    def test_authorized_retry_has_distinct_token(self):self.assertEqual(self.get(),'a'*32)
    def test_changed_inputs_rejected(self):
        with self.assertRaises(ValueError):retry_token(self.r,self.d,'j','manifest',{'goat.inp':'changed'})
    def test_changed_previous_result_rejected(self):
        writej(self.old/'receipt.json',{'status':'COMPLETED'})
        with self.assertRaises(ValueError):self.get()
    def test_changed_proof_rejected(self):
        writej(self.proof,{'scheduler_state':'RUNNING'})
        with self.assertRaises(ValueError):self.get()
    def test_mpi_retry_requires_absent_processes_and_failed_receipt(self):
        writej(self.old/'receipt.json',{'status':'FAILED'})
        (self.old/'goat.stderr').write_text('mpirun exited improperly')
        proof={'job_id':'j','stage':str(self.old),'receipt_sha256':sha(self.old/'receipt.json'),
               'failure_class':'MPI_ABORT_BEFORE_FIRST_GLOBAL_ITERATION',
               'producer_exit_proof':{'matches':[],'failures':[]},
               'mpi_stderr_sha256':sha(self.old/'goat.stderr')}
        writej(self.proof,proof)
        self.auth.update(status='PREPARED_VERIFIED_MPI_ABORT',proof_sha256=sha(self.proof))
        writej(self.d/'goat_retry_authorization.json',self.auth)
        self.assertEqual(self.get(),'a'*32)
        proof['producer_exit_proof']['matches']=[{'pid':123}];writej(self.proof,proof)
        self.auth['proof_sha256']=sha(self.proof);writej(self.d/'goat_retry_authorization.json',self.auth)
        with self.assertRaises(ValueError):self.get()

    def test_fresh_job_has_no_retry_token(self):
        (self.d/'goat_retry_authorization.json').unlink();self.assertIsNone(self.get())

class SchedulerEvidenceTests(unittest.TestCase):
    def test_failed_sibling_does_not_authorize_live_parent(self):
        from scripts.recover_goat_nodes_20261022 import exact_failed
        with self.assertRaises(ValueError):exact_failed([['10','a_0','RUNNING'],['11','a_1','NODE_FAIL']],'10')
    def test_active_child_blocks_retry(self):
        from scripts.recover_goat_nodes_20261022 import exact_failed
        with self.assertRaises(ValueError):exact_failed([['10','a_0','NODE_FAIL'],['10.0','a_0.0','RUNNING']],'10')
    def test_exact_terminal_allocation_accepted(self):
        from scripts.recover_goat_nodes_20261022 import exact_failed
        self.assertEqual(exact_failed([['10','a_0','NODE_FAIL'],['11','a_1','RUNNING']],'10')[0],'10')
