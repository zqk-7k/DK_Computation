import tempfile,time,unittest,json
from pathlib import Path
from unittest.mock import patch
from clusterflow.common import writej,sha
from scripts import prepare_goat_recovery_20261020 as m

class BoundedRecoveryTests(unittest.TestCase):
 def setUp(self):
  self.t=tempfile.TemporaryDirectory();self.addCleanup(self.t.cleanup);self.r=Path(self.t.name);self.d=self.r/'job';self.d.mkdir();self.s=self.r/'scratch';self.s.mkdir();self.row={'job_directory':'job'}
 def call(self,j,q=None):
  with patch.object(m,'validate_job'),patch.object(m,'verify_original_inputs'):return m.inspect(self.r,j,self.row,q or {'active':{}})
 def stage(self,status):
  writej(self.d/'stage_intents/i.json',{'stage':str(self.s)});writej(self.s/'receipt.json',{'status':status})
 def test_only_named_fresh_seeds_allowed(self):
  self.assertEqual(self.call(m.JOBS[-1])['retry_number'],0)
  with self.assertRaisesRegex(ValueError,'Unexpected fresh'):self.call(m.JOBS[0])
 def test_old_outputs_allowed_but_changed_input_rejected(self):
  row=dict(self.row,initial_release=True)
  for name,field in dict(m.INPUT_HASH_FIELDS,**{'provenance.json':'provenance_sha256'}).items():
   (self.d/name).write_text('synthetic '+name);row[field]=sha(self.d/name)
  (self.d/'goat.out').write_text('retained failed attempt')
  m.verify_original_inputs(self.r,row)
  (self.d/'goat.inp').write_text('changed')
  with self.assertRaisesRegex(ValueError,'identity changed'):m.verify_original_inputs(self.r,row)
 def test_existing_delegation_prevents_second_attempt(self):
  writej(self.d/'delegated_goat.json',{})
  with self.assertRaisesRegex(ValueError,'Already delegated'):self.call(m.JOBS[-1])
 def test_completed_output_never_retried(self):
  self.stage('COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW')
  with self.assertRaisesRegex(ValueError,'Unknown producer'):self.call(m.JOBS[0])
 def test_wrapper_only_failure_never_retried(self):
  self.stage('FAILED');writej(self.s/'receipt.json',{'status':'FAILED','errors':['ORCA nonzero or absent exit code']})
  with patch.object(m,'validate_terminal'):
   with self.assertRaisesRegex(ValueError,'Recover valid output'):self.call(m.JOBS[0])
 def test_fresh_heartbeat_blocks_orphan_recovery(self):
  j=m.JOBS[6];self.stage('RUNNING');writej(self.s/'runtime_heartbeat.json',{'updated_epoch':time.time()})
  with self.assertRaisesRegex(ValueError,'Fresh producer'):self.call(j,{'active':{j:{}}})
 def test_live_orphan_child_blocks_retry(self):
  j=m.JOBS[6];self.stage('RUNNING');writej(self.s/'runtime_heartbeat.json',{'updated_epoch':0,'pid':123})
  q={'active':{j:{'hostname':'i15r4n18','slurm_job_id':'124489624'}}}
  with patch.object(m.subprocess,'run') as run:
   run.return_value.stdout=json.dumps({'matches':[{'pid':456}],'failures':[]})
   with self.assertRaisesRegex(ValueError,'not conclusively absent'):self.call(j,q)
