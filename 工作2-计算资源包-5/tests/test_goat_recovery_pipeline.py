"""Synthetic fixtures: technical recovery must retain all chemistry gates."""
import shutil,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
import test_cluster_pipeline as fixture
from clusterflow import pipeline,autoflow
from clusterflow.common import readj,writej,sha
from clusterflow.goat_cleanup_review import evidence
from clusterflow.goat_delegation import FILES
from scripts.expand_goat_20261020 import resource_input

class RecoveryPipelineTests(unittest.TestCase):
 setUp=fixture.ClusterPipelineTests.setUp
 completed_goat=fixture.ClusterPipelineTests.completed_goat
 def snapshot(self):
  for name in ('config/automation.json','registry/batch_plan_v5.json'):
   shutil.copy2(fixture.ROOT/name,self.root/name)
  wave={'workspace':'automation/synthetic/wave_001/workspace'}
  workspace=autoflow._workspace(self.root,wave,{'fixture':{'job_ids':[self.row['job_id']]}},{self.row['job_id']:self.row})
  self.assertEqual(autoflow._verify_workspace(self.root,wave),workspace)
  pipeline._goat_source(workspace,self.row['job_id'])
  return wave
 def failed_cleanup(self):
  self.completed_goat();(self.goat/'goat.stderr').write_text('rm: cannot remove /var/tmp/work2-runtime-test/orca610_sandbox_test/a: Read-only file system\n')
  r=readj(self.goat/'receipt.json');r.update(status='FAILED',returncode=1,errors=['ORCA nonzero or absent exit code']);r['output_sha256']['goat.stderr']=sha(self.goat/'goat.stderr');writej(self.goat/'receipt.json',r)
 def review(self):
  writej(self.goat/'goat_cleanup_review.json',dict(status='NORMAL_ORCA_OUTPUT_WITH_FAILED_SANDBOX_CLEANUP',evidence=evidence(self.goat,self.goat/'receipt.json')))
 def test_failed_without_technical_review_stays_failed(self):
  self.failed_cleanup();self.assertEqual(autoflow._job_status(self.root,self.row,goat=True)['status'],'FAILED_REVIEW')
 def test_technical_review_preserves_failed_receipt_and_enters_chemistry_gates(self):
  self.failed_cleanup();before=sha(self.goat/'receipt.json');self.review();source=pipeline._goat_source(self.root,self.row['job_id'])
  self.assertIn('cleanup_review_sha256',source['bindings']);self.assertEqual(sha(self.goat/'receipt.json'),before)
  status=autoflow._job_status(self.root,self.row,goat=True);self.assertEqual(status['status'],'COMPLETE');self.assertEqual(status['original_receipt_status'],'FAILED')
  self.snapshot()
 def test_hash_bound_wrong_atom_order_still_rejected(self):
  self.failed_cleanup();p=self.goat/'goat.finalensemble.xyz';lines=p.read_text().splitlines();lines[2]='Xe '+ ' '.join(lines[2].split()[1:]);p.write_text('\n'.join(lines)+'\n')
  r=readj(self.goat/'receipt.json');r['output_sha256'][p.name]=sha(p);writej(self.goat/'receipt.json',r);self.review()
  self.assertEqual(autoflow._job_status(self.root,self.row,goat=True)['status'],'INVALID_REVIEW')
 def test_delegated_completed_output_is_recognized_despite_canonical_reservation(self):
  self.completed_goat();(self.goat/'isotope_map.json').write_text('{}')
  c=self.root/'operations/cohort';c.mkdir(parents=True)
  for name in ('config','registry','01_goat_explicit'):shutil.copytree(self.root/name,c/name)
  d=c/self.row['job_directory'];(d/'goat.inp').write_text(resource_input((d/'goat.inp').read_text()))
  row=dict(self.row,input_sha256=sha(d/'goat.inp'));writej(c/'registry/jobs.json',[row]);r=readj(d/'receipt.json');r['input_sha256']['goat.inp']=sha(d/'goat.inp');writej(d/'receipt.json',r)
  record=dict(kind='RESOURCE_COHORT_DELEGATION',token='synthetic',job_id=self.row['job_id'],cpus=16,cohort=str(c),canonical_manifest_sha256=sha(self.root/'registry/jobs.json'),cohort_manifest_sha256=sha(c/'registry/jobs.json'),canonical_inputs={n:sha(self.goat/n) for n in FILES})
  writej(self.goat/'delegated_goat.json',record);writej(self.goat/'.goat-run.lock',record)
  s=autoflow._job_status(self.root,self.row,goat=True);self.assertEqual(s['status'],'COMPLETE');self.assertEqual(s['canonical_directory'],str(self.goat))
  self.assertIn('delegation_sha256',pipeline._goat_source(self.root,self.row['job_id'])['bindings'])
  wave=self.snapshot()
  record['token']='changed';writej(self.goat/'delegated_goat.json',record)
  with self.assertRaisesRegex(ValueError,'delegation binding changed'):autoflow._verify_workspace(self.root,wave)
