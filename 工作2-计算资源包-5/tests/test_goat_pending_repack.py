import tempfile,unittest
from pathlib import Path
from clusterflow.common import writej
from scripts.repack_pending_goat_20261020 import no_started_work

class RepackTests(unittest.TestCase):
 def setUp(self):
  self.t=tempfile.TemporaryDirectory();self.addCleanup(self.t.cleanup);self.root=Path(self.t.name);self.base=self.root/'operations/pool';self.job=self.root/'job';self.job.mkdir();writej(self.root/'registry/jobs.json',[{'job_id':'j','job_directory':'job'}]);writej(self.base/'scope.json',{'job_ids':['j']})
 def test_pending_without_any_execution_evidence_can_be_repacked(self):
  self.assertEqual(len(no_started_work(self.root,self.base,'10_0|PENDING\n10_1|PENDING\n')),2)
 def test_any_running_element_blocks_cancellation(self):
  with self.assertRaisesRegex(ValueError,'PENDING'):no_started_work(self.root,self.base,'10_0|RUNNING\n10_1|PENDING\n')
 def test_missing_scheduler_records_are_not_permission_to_replace(self):
  with self.assertRaises(ValueError):no_started_work(self.root,self.base,'')
 def test_consumed_queue_blocks_even_if_scheduler_pending(self):
  writej(self.base/'queue.json',{'cursor':1,'active':{}})
  with self.assertRaisesRegex(ValueError,'consumed'):no_started_work(self.root,self.base,'10|PENDING')
 def test_durable_run_evidence_blocks_repack(self):
  writej(self.job/'.goat-run.lock',{'token':'producer'})
  with self.assertRaisesRegex(ValueError,'execution evidence'):no_started_work(self.root,self.base,'10|PENDING')
