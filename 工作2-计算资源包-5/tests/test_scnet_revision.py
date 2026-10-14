"""Regression evidence only: no ORCA, SSH, sbatch or production mutation."""
import copy
import json
import os
from pathlib import Path
import socket
import shutil
import tempfile
import unittest
from unittest.mock import patch
from rdkit import Chem
from clusterflow import runtime,operations,pipeline,autoflow
from clusterflow.common import readj,writej,sha
from scripts import goat_job,goat_portable
import test_goat_portable as portable_fixture

ROOT=Path(__file__).resolve().parents[1]


class ReleaseRegression(portable_fixture.GoatPortableTests):
    def test_held_release_survives_real_staged_validator(self):
        self.row.update(initial_release=False,species_id='example')
        writej(self.root/'registry/jobs.json',[self.row])
        writej(self.root/'config/automation.json',{'automatic_handoff':{'enabled':True}})
        writej(self.root/'config/protocol.json',{})
        writej(self.root/'registry/species.json',[])
        writej(self.root/'registry/batch_plan_v5.json',{'species':[{'species_id':'example','goat_start_ids':['j1']}]})
        names=('registry/batch_plan_v5.json','registry/jobs.json','registry/species.json','config/protocol.json','config/automation.json')
        writej(self.root/'automation/example/state.json',{'species_id':'example',
               'identity_sha256':{n:sha(self.root/n) for n in names},'released':{'example':['j1']}})
        def real_validation_mock_quantum(stage,*args,**kwargs):
            goat_job.validate_job(stage,'j1')  # regression: formerly rejected here
            self.assertFalse(readj(stage/'registry/jobs.json')[0]['initial_release'])
            return self.mock_run(stage,*args,**kwargs)
        with patch.object(goat_job,'run_one',side_effect=real_validation_mock_quantum):
            result=self.call()
        self.assertIn('COMPLETED',result['status'])

    def test_unreleased_held_task_still_rejected(self):
        self.row['initial_release']=False
        writej(self.root/'registry/jobs.json',[self.row])
        with self.assertRaisesRegex(ValueError,'Held job'):self.call()


class RuntimeRegression(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)

    def test_live_owner_not_mistaken_for_dead(self):
        path=self.root/'lock';writej(path,runtime.owner_record('test'))
        self.assertEqual(runtime.lock_state(path),'RUNNING')

    def test_dead_pid_is_interrupted_not_running(self):
        path=self.root/'lock';record=runtime.owner_record('test');record['pid']=999999999
        writej(path,record)
        self.assertEqual(runtime.lock_state(path),'INTERRUPTED_REVIEW')
        self.assertTrue(path.exists())

    def test_unknown_owner_not_deleted_or_claimed_running(self):
        path=self.root/'lock';writej(path,{})
        self.assertEqual(runtime.lock_state(path),'UNKNOWN_REVIEW')
        self.assertFalse(operations.retire_terminal_lock(path))
        self.assertTrue(path.exists())

    def test_remote_running_parent_does_not_prove_child_alive(self):
        path=self.root/'lock';writej(path,{'hostname':'other','slurm_job_id':'123'})
        with patch('subprocess.run',return_value=type('R',(),{'returncode':0,'stdout':'123|RUNNING\n'})()):
            self.assertEqual(runtime.lock_state(path),'REMOTE_OWNER_CHECK_REQUIRED')

    def test_terminal_scheduler_lock_archived_not_deleted(self):
        path=self.root/'lock';writej(path,{'hostname':'other','slurm_job_id':'123'})
        with patch('subprocess.run',return_value=type('R',(),{'returncode':0,'stdout':'123|TIMEOUT\n'})()):
            self.assertTrue(operations.retire_terminal_lock(path))
        self.assertEqual(len(list((self.root/'lock_history').glob('*.json'))),1)

    def test_readonly_mpi_not_mislabeled_oom(self):
        self.assertEqual(runtime.classify_failure('MPI_Win_allocate_shared Read-only file system errno 30'),'NODE_READ_ONLY')
        self.assertEqual(runtime.classify_failure('PROPINT'),'MODULE_FAILURE_REVIEW')
        self.assertEqual(runtime.classify_failure('TMatrixContainers::AddMatrix'),'MATRIX_IO_REVIEW')
        self.assertEqual(runtime.classify_failure('No valid structures found'),'GOAT_NO_VALID_ENSEMBLE')
        self.assertEqual(runtime.classify_failure('out of memory'),'OOM')

    def test_slot_memory_instead_of_whole_parent_memory(self):
        from clusterflow.runner import allocation
        with patch.dict(os.environ,{'SLURM_MEM_PER_NODE':'114000','WORK2_SLOT_MEMORY_MB':'20000','SLURM_NTASKS':'1','SLURM_CPUS_PER_TASK':'16'}):
            with self.assertRaisesRegex(ValueError,'memory'):allocation(8,3500)
            allocation(4,3500)

    def test_version_evidence_bound_to_wrapper_container_xtb_and_site(self):
        site=readj(ROOT/'config/site_scnet.json')
        for field in ('orca_wrapper','container','xtb'):
            path=self.root/field;path.write_text('synthetic '+field);path.chmod(0o700);site[field]=str(path)
        profile=self.root/'site.json';writej(profile,site)
        evidence=self.root/'evidence.json'
        writej(evidence,{'version':'6.1.0','wrapper_sha256':sha(site['orca_wrapper']),
                        'container_sha256':sha(site['container']),'xtb_sha256':sha(site['xtb']),'site_sha256':sha(profile)})
        with patch.dict(os.environ,{'WORK2_SITE_CONFIG':str(profile),'WORK2_RUNTIME_EVIDENCE':str(evidence)}):
            self.assertTrue(runtime.wrapper_preflight(site['orca_wrapper'],site['xtb'])['passed'])
            Path(site['container']).write_text('changed')
            self.assertFalse(runtime.wrapper_preflight(site['orca_wrapper'],site['xtb'])['passed'])

    def test_missing_runtime_evidence_blocks_wrapper(self):
        site=readj(ROOT/'config/site_scnet.json')
        with patch.dict(os.environ,{},clear=True):
            self.assertFalse(runtime.wrapper_preflight(site['orca_wrapper'])['passed'])

    def test_64_products_66_models_remain_and_previous_execution_not_live(self):
        report=operations.coverage(ROOT)
        self.assertEqual((report['product_count'],report['model_count']),(64,66))
        self.assertTrue(all(not b['final_delivery_approved'] for rows in report['products'].values() for b in rows))

    def test_water_equivalent_hydrogen_permutation_deduplicates(self):
        row=readj(ROOT/'registry/jobs.json')[0];directory=ROOT/row['job_directory']
        mol=Chem.MolFromMolBlock((directory/'cluster.mol').read_text(),removeHs=False)
        amap=readj(directory/'atom_map.json');other=Chem.Mol(mol)
        a,b=[r['index_zero_based'] for r in amap if r['component'].startswith('H2O_') and r['element']=='H']
        other.GetConformer().SetAtomPosition(a,mol.GetConformer().GetAtomPosition(b))
        other.GetConformer().SetAtomPosition(b,mol.GetConformer().GetAtomPosition(a))
        self.assertLess(pipeline.cluster_rms(mol,other,amap),1e-5)

    def test_stageout_recovery_does_not_rerun_or_overwrite(self):
        original=self.root/'jobs/test';stage=self.root/'scratch/test'
        original.mkdir(parents=True);stage.mkdir(parents=True)
        for directory in (original,stage):(directory/'job.inp').write_text('synthetic input')
        (stage/'job.out').write_text('synthetic output')
        inputs={'job.inp':sha(original/'job.inp')}
        writej(stage/'receipt.json',{'id':'test','status':'completed','input_sha256':inputs,
               'output_sha256':{'job.out':sha(stage/'job.out')}})
        intent=original/'stage_intents/one.json'
        writej(intent,{'original':str(original),'stage':str(stage),'kind':'DFT','inputs':inputs,'job_identity':'test'})
        with patch.object(operations,'scheduler_terminal',return_value=True),patch.object(runtime,'site_config',return_value={'shared_scratch':str(self.root/'scratch')}),patch('subprocess.Popen',side_effect=AssertionError('No QM')):
            result=operations.recover_stageout(self.root,intent,execute=True)
            self.assertEqual(result['quantum_jobs_launched'],0)
            self.assertEqual((original/'receipt.json').read_bytes(),(stage/'receipt.json').read_bytes())
            self.assertTrue(stage.exists())
            (original/'job.out').write_text('conflict')
            with self.assertRaisesRegex(ValueError,'Conflicting'):operations.recover_stageout(self.root,intent,execute=True)

    def test_no_stageout_recovery_with_live_or_unknown_scheduler(self):
        original=self.root/'jobs/test';stage=self.root/'scratch/test'
        intent=original/'stage_intents/one.json'
        writej(intent,{'original':str(original),'stage':str(stage),'kind':'DFT','inputs':{},'job_identity':'test'})
        with patch.object(operations,'scheduler_terminal',return_value=False),patch.object(runtime,'site_config',return_value={'shared_scratch':str(self.root/'scratch')}):
            with self.assertRaisesRegex(ValueError,'terminal'):operations.recover_stageout(self.root,intent,execute=True)


class ExactReuseRegression(unittest.TestCase):
    def setUp(self):
        import test_cluster_pipeline as fixture
        self.fixture=fixture.ClusterPipelineTests();self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root=self.fixture.root
        writej(self.root/'registry/batch_plan_v5.json',{})
        self.source=self.root/self.fixture.prepared()['jobs'][0]
        self.fixture.completed_opt(self.source)
        receipt=readj(self.source/'receipt.json');receipt.update(orca_executable_sha256='synthetic-executable',runtime_evidence_sha256=None)
        writej(self.source/'receipt.json',receipt)
        self.target=self.source.parent/'new_wave'
        self.target.mkdir()
        for name in ('job.inp','job.json','input.xyz','input_structure.mol','atom_map.json'):
            shutil.copyfile(self.source/name,self.target/name)
        metadata=readj(self.target/'job.json');metadata['id']='synthetic-new-wave'
        writej(self.target/'job.json',metadata)

    def test_reuse_real_gates_accepts_bound_exact_input_and_detects_tampering(self):
        from clusterflow import cache,gates
        with patch.dict(os.environ,{},clear=True):
            cache.publish(self.root,self.source)
            receipt=cache.try_reuse(self.root,self.target,'synthetic-executable')
        self.assertEqual(receipt['quantum_jobs_launched'],0)
        self.assertEqual(receipt['receipt_kind'],'VERIFIED_REUSE_NOT_NEW_ORCA_EXECUTION')
        self.assertTrue(gates.inspect_job(self.target,self.fixture.cfg)['accepted'])
        (self.source/'job.out').write_text('tampered')
        self.assertFalse(gates.inspect_job(self.target,self.fixture.cfg)['accepted'])

    def test_different_input_not_reused(self):
        from clusterflow import cache
        cache.publish(self.root,self.source)
        (self.target/'job.inp').write_text('different input')
        self.assertIsNone(cache.try_reuse(self.root,self.target,'synthetic-executable'))


class PackedRegression(unittest.TestCase):
    def setUp(self):
        from scripts import packed_worker
        self.module=packed_worker
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.base=self.root/'operations/packed'
        self.site=readj(ROOT/'config/site_scnet.json')
        self.root_patch=patch.object(packed_worker,'ROOT',self.root);self.root_patch.start();self.addCleanup(self.root_patch.stop)
        self.storage=patch.object(packed_worker,'check_storage',return_value={'free_bytes':10000*1024**3});self.storage.start();self.addCleanup(self.storage.stop)
        self.terminal=patch.object(packed_worker,'scheduler_terminal',return_value=False);self.terminal.start();self.addCleanup(self.terminal.stop)

    def test_cross_worker_exclusive_claims(self):
        a=self.module.claim(self.base,self.site,['a','b'])
        b=self.module.claim(self.base,self.site,['a','b'])
        self.assertEqual({a[1]['species_id'],b[1]['species_id']},{'a','b'})
        self.assertIsNone(self.module.claim(self.base,self.site,['a','b']))

    def test_shared_disk_growth_reservation_blocks_overcommit(self):
        with patch.object(self.module,'check_storage',return_value={'free_bytes':225*1024**3}):
            self.assertIsNotNone(self.module.claim(self.base,self.site,['a','b']))
            self.assertIsNone(self.module.claim(self.base,self.site,['a','b']))

    def test_durable_scientific_block_not_repeatedly_claimed(self):
        writej(self.base/'blocked/a.json',{'retry_on_different_node':False})
        self.assertIsNone(self.module.claim(self.base,self.site,['a']))

    def test_unresolved_legacy_lock_is_not_cleared(self):
        path=self.root/'automation/a/.autoflow.lock';writej(path,{})
        self.assertIsNone(self.module.claim(self.base,self.site,['a']))
        self.assertTrue(path.exists())

    def test_operator_pause_is_not_treated_as_completion(self):
        writej(self.root/'automation/a/state.json',{'phase':'REVIEW'})
        self.assertIsNone(self.module.claim(self.base,self.site,['a']))

    def test_nonterminal_receipt_without_lock_is_not_running(self):
        job=self.root/'goat/j1';writej(job/'receipt.json',{'status':'RUNNING'})
        report=autoflow._job_status(self.root,{'job_directory':'goat/j1'},goat=True)
        self.assertEqual(report['status'],'INCOMPLETE_REVIEW')


if __name__=='__main__':unittest.main()
