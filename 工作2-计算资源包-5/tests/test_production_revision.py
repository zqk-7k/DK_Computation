"""Regression and integration checks; temporary fixtures, no quantum execution."""
import contextlib
import copy
import io
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest
from unittest.mock import patch

from clusterflow import scheduling, runtime, gates, pipeline, autoflow, review_control
from clusterflow.common import readj,writej,sha
from scripts import packed_worker,submit_packed
import test_cluster_pipeline as fixtures

ROOT=Path(__file__).resolve().parents[1]
SID=fixtures.SPECIES
SCENARIO=fixtures.SCENARIO


class AdmissionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        writej(self.root/'registry/batch_plan_v5.json',{'species':[{'species_id':SID,'complexity_class':'X'}]})
        writej(self.root/'config/protocol.json',{})

    def test_long_task_deferred_before_short_parent_deadline(self):
        with patch.dict(os.environ,{'WORK2_DEADLINE_EPOCH':str(time.time()+7200)}):
            with self.assertRaises(scheduling.DeferredWork) as caught:scheduling.admit_quantum(self.root,'optfreq',SID)
        self.assertEqual(caught.exception.report['status'],'WAITING_WALLTIME')
        self.assertGreater(caught.exception.report['required_seconds'],7200)

    def test_30_day_parent_accepts_unmeasured_complex_task_without_timeout(self):
        with patch.dict(os.environ,{'WORK2_DEADLINE_EPOCH':str(time.time()+720*3600)}):
            self.assertIsNone(scheduling.admit_quantum(self.root,'optfreq',SID))

    def test_observed_slow_job_increases_budget(self):
        before,_=scheduling.duration_budget(self.root,'optfreq',SID)
        scheduling.record_timing(self.root,'optfreq',SID,before*2)
        after,_=scheduling.duration_budget(self.root,'optfreq',SID)
        self.assertGreaterEqual(after,3*before)

    def test_cache_key_execution_lock_defers_second_owner_and_releases(self):
        with scheduling.single_execution(self.root,'same_tms_input'):
            with self.assertRaises(scheduling.DeferredWork) as caught:
                with scheduling.single_execution(self.root,'same_tms_input'):self.fail('Second owner entered')
            self.assertEqual(caught.exception.report['status'],'WAITING_SHARED_RESULT')
        with scheduling.single_execution(self.root,'same_tms_input'):pass

    def test_crashed_driver_marker_prevents_duplicate_from_other_workspace(self):
        marker=self.root/'operations/cache_locks/orphan.owner.json'
        writej(marker,{'source_directory':str(self.root/'first'),'slurm_job_id':'123'})
        for terminal in (False,True):
            with patch('clusterflow.operations.scheduler_terminal',return_value=terminal), \
                 patch('clusterflow.runtime.lock_state',return_value='INTERRUPTED_REVIEW'):
                with self.assertRaises(scheduling.DeferredWork) as caught:
                    with scheduling.single_execution(self.root,'orphan',self.root/'second'):self.fail('Orphan duplicated')
            self.assertEqual(caught.exception.report['status'],'WAITING_SHARED_RESULT_RECOVERY')
            self.assertTrue(marker.exists())

    def test_reconciled_terminal_original_archives_owner_before_retry(self):
        marker=self.root/'operations/cache_locks/orphan.owner.json'
        writej(marker,{'source_directory':str(self.root/'first'),'slurm_job_id':'123'})
        with patch('clusterflow.operations.scheduler_terminal',return_value=True), \
             patch('clusterflow.runtime.lock_state',return_value='INTERRUPTED_REVIEW'):
            with scheduling.single_execution(self.root,'orphan',self.root/'first'):
                self.assertIn('token',readj(marker))
        self.assertFalse(marker.exists())
        self.assertEqual(len(list((marker.parent/'owner_history').glob('*.json'))),1)

    def test_shared_quota_or_disk_wait_does_not_quarantine_node(self):
        site=readj(ROOT/'config/site_scnet.json')
        writej(self.root/'registry/batch_plan_v5.json',{'species':[{'species_id':SID,'compound_id':'F2_E01'}]})
        with patch.object(packed_worker,'ROOT',self.root),patch.object(packed_worker,'site_config',return_value=site), \
             patch.object(packed_worker,'check_storage',side_effect=runtime.StorageAdmissionError('QUOTA_ADMISSION_BLOCKED','/public')), \
             patch('os.sched_getaffinity',return_value=set(range(32))),patch('os.sched_setaffinity'),patch('time.sleep'), \
             patch.dict(os.environ,{'SLURM_JOB_ID':'123456789'},clear=False):
            packed_worker.worker(0,4,time.time()+720*3600,1)
        self.assertFalse((self.root/'operations/packed/quarantine').exists())
        self.assertTrue(list((self.root/'operations/packed').glob('admission_*.json')))

    def test_parallel_branch_reserves_space_for_both_simultaneous_children(self):
        site=readj(ROOT/'config/site_scnet.json')
        site['packed'].update(site['packed_profiles']['branch_parallel'])
        base=self.root/'operations/packed'
        with patch.object(packed_worker,'ROOT',self.root), \
             patch.object(packed_worker,'check_storage',side_effect=[{'free_bytes':400*1024**3},{'free_bytes':350*1024**3}]):
            _,owner=packed_worker.claim(base,site,[SID,'second_branch'])
            self.assertEqual(owner['growth_gib'],128)
            # 128 GiB floor + 128 existing + 128 new > 350 GiB free.
            self.assertIsNone(packed_worker.claim(base,site,[SID,'second_branch']))

    def test_site_profiles_fit_cpu_and_memory_for_prescreen(self):
        from clusterflow.runner import allocation
        for profile in ('balanced','large_memory','branch_parallel'):
            with patch.dict(os.environ,{'WORK2_PACKED_PROFILE':profile,'SLURM_NTASKS':'1'}):
                site=runtime.site_config();slots=site['packed'];parallel=slots['branch_parallel_jobs']
                self.assertLessEqual(slots['slots_per_node']*slots['slot_memory_mb'],114000)
                with patch.dict(os.environ,{'SLURM_CPUS_PER_TASK':str(slots['slot_allocated_cpus']//parallel),
                                           'WORK2_SLOT_MEMORY_MB':str(slots['slot_memory_mb']//parallel)}):allocation(8,2500)

    def test_submit_default_is_30_days_and_dry_run_has_no_subprocess(self):
        output=io.StringIO()
        with patch('sys.argv',['submit_packed.py','--nodes','1']),patch.dict(os.environ,{},clear=False), \
             patch.object(submit_packed.subprocess,'run',side_effect=AssertionError('No scheduler mutation')),contextlib.redirect_stdout(output):
            submit_packed.main()
        self.assertIn('"walltime_hours": 720',output.getvalue())

    def test_execute_refuses_failed_scheduler_test_without_submitting(self):
        with patch('sys.argv',['submit_packed.py','--nodes','1','--execute']), \
             patch.dict(os.environ,{'WORK2_RUNTIME_EVIDENCE':'synthetic'},clear=False), \
             patch.object(submit_packed.subprocess,'run',side_effect=subprocess.CalledProcessError(1,'sbatch')) as run:
            with self.assertRaises(subprocess.CalledProcessError):submit_packed.main()
        self.assertEqual(run.call_count,1)
        self.assertIn('--test-only',run.call_args.args[0])

    def test_active_or_unknown_previous_chain_cannot_be_resubmitted(self):
        writej(self.root/'operations/submission_old.json',{'status':'SUBMITTED_NOT_RUNNING','jobs':['123']})
        with patch.object(submit_packed,'ROOT',self.root),patch.object(submit_packed,'scheduler_terminal',return_value=False), \
             patch.object(submit_packed.subprocess,'run') as run:
            with self.assertRaisesRegex(ValueError,'active or not confirmed'):
                submit_packed.submit_chain(['sbatch'],{},2,{})
        run.assert_not_called()

    def test_successful_chain_dependencies_and_terminal_predecessor(self):
        writej(self.root/'operations/submission_old.json',{'status':'SUBMITTED_NOT_RUNNING','jobs':['123']})
        replies=[subprocess.CompletedProcess([],0,'456\n',''),subprocess.CompletedProcess([],0,'457\n','')]
        with patch.object(submit_packed,'ROOT',self.root),patch.object(submit_packed,'scheduler_terminal',return_value=True), \
             patch.object(submit_packed.subprocess,'run',side_effect=replies) as run:
            report=submit_packed.submit_chain(['sbatch'],{},2,{})
        self.assertEqual(report['jobs'],['456','457'])
        self.assertIn('--dependency=afterany:456',run.call_args_list[1].args[0])
        self.assertEqual(readj(report['journal'])['status'],'SUBMITTED_NOT_RUNNING')

    def test_ambiguous_submission_is_preserved_and_blocks_blind_retry(self):
        reply=subprocess.CompletedProcess([],0,'not a job id\n','')
        with patch.object(submit_packed,'ROOT',self.root),patch.object(submit_packed.subprocess,'run',return_value=reply):
            with self.assertRaisesRegex(ValueError,'Unrecognized'):
                submit_packed.submit_chain(['sbatch'],{},2,{})
        journal=next((self.root/'operations').glob('submission_*.json'))
        self.assertIn('UNCERTAIN',readj(journal)['status'])
        with self.assertRaisesRegex(ValueError,'Unresolved'):submit_packed.verify_previous_submissions(self.root)


class IncrementalTests(unittest.TestCase):
    def setUp(self):
        self.f=fixtures.ClusterPipelineTests();self.f.setUp();self.addCleanup(self.f.doCleanups)
        self.f.cfg['prescreen']['frequency_required']=False
        writej(self.f.root/'config/protocol.json',self.f.cfg)
        self.job=self.f.root/self.f.prepared()['jobs'][0]

    def test_prescreen_uses_real_electronic_energy_without_freq_or_fake_G(self):
        self.assertNotIn(' Freq',(self.job/'job.inp').read_text())
        self.f.completed_opt(self.job)
        out=self.job/'job.out'
        out.write_text('Program Version 6.1.0\nTHE OPTIMIZATION HAS CONVERGED\nFINAL SINGLE POINT ENERGY -123.4\nORCA TERMINATED NORMALLY\n')
        receipt=readj(self.job/'receipt.json');receipt['output_sha256']['job.out']=sha(out);writej(self.job/'receipt.json',receipt)
        result=gates.inspect_job(self.job,self.f.cfg)
        self.assertTrue(result['accepted'],result['errors'])
        self.assertIsNone(result['parsed']['gibbs_free_energy_Eh'])
        pipeline.prepare_healthy_downstream(self.f.root,SID,SCENARIO)
        self.assertEqual(len(pipeline.manifest(self.f.root,'optfreq',SID,SCENARIO)),1)
        self.assertFalse((self.job.parent.parent/'prescreen_gate.json').exists())

    def test_final_negative_frequency_still_blocks_properties(self):
        self.f.completed_opt(self.job,negative=True)
        self.assertTrue(gates.inspect_job(self.job,self.f.cfg)['accepted'])
        pipeline.prepare_healthy_downstream(self.f.root,SID,SCENARIO)
        opt=self.f.root/pipeline.manifest(self.f.root,'optfreq',SID,SCENARIO)[0]['job_dir']
        self.f.completed_opt(opt,negative=True)
        self.assertFalse(gates.inspect_job(opt,self.f.cfg)['accepted'])
        pipeline.prepare_healthy_downstream(self.f.root,SID,SCENARIO)
        self.assertFalse(pipeline.manifest(self.f.root,'shielding',SID,SCENARIO))

    def test_healthy_properties_prepared_without_claiming_final_ensemble(self):
        self.f.completed_opt(self.job)
        pipeline.prepare_healthy_downstream(self.f.root,SID,SCENARIO)
        opt=self.f.root/pipeline.manifest(self.f.root,'optfreq',SID,SCENARIO)[0]['job_dir']
        self.f.completed_opt(opt)
        pipeline.prepare_healthy_downstream(self.f.root,SID,SCENARIO)
        self.assertTrue(pipeline.manifest(self.f.root,'shielding',SID,SCENARIO))
        self.assertFalse(list((self.f.root/'ensembles').glob('**/ensemble.json')))
        pipeline.advance(self.f.root,SID,SCENARIO,'prescreen')
        self.assertEqual(pipeline.advance(self.f.root,SID,SCENARIO,'optfreq')['status'],'advanced')

    def test_failed_member_does_not_prevent_healthy_property_work(self):
        from test_v5_integration import IndependentFinalStableMinimumGate
        fixture=IndependentFinalStableMinimumGate();fixture.setUp()
        try:
            directories=fixture.prepared_final_stage(3)
            receipt=readj(directories[1]/'receipt.json');receipt['status']='failed'
            writej(directories[1]/'receipt.json',receipt)
            with patch.object(autoflow,'_run_pending',return_value={'status':'SYNTHETIC_JOB_DISPATCH','quantum_jobs_launched':0}) as dispatch:
                result=autoflow._stage_step(fixture.root,{SID:{'stable_conformer_target':3}},execute=True,orca='unused',scratch_root='unused')
            self.assertEqual(result['status'],'SYNTHETIC_JOB_DISPATCH')
            rows=dispatch.call_args.args[1]
            self.assertTrue(any(row['stage']=='shielding' for row in rows))
            self.assertNotIn(directories[1].name,{row['conformer_id'] for row in rows})
            self.assertFalse(list((fixture.root/'ensembles').glob('**/ensemble.json')))
        finally:fixture.doCleanups()

    def test_later_prescreen_representative_supersedes_early_opt_without_deleting_it(self):
        from rdkit import Chem
        from clusterflow import autoselect
        from test_v5_integration import IndependentFinalStableMinimumGate
        fixture=IndependentFinalStableMinimumGate();fixture.setUp()
        try:
            cfg=readj(fixture.root/'config/protocol.json');cfg['prescreen']['frequency_required']=False
            writej(fixture.root/'config/protocol.json',cfg)
            models=[fixture.real_initial_model(i) for i in (1,2)]
            fixture.completed(fixture.rows[0],models,[-100,-99.999])
            selected=autoselect.screen_batch(fixture.root,[fixture.rows[0]['job_id']],2)
            autoselect.import_selection(fixture.root,selected)
            batch=pipeline.prepare_prescreen(fixture.root,SID,SCENARIO)
            first,second=[fixture.root/path for path in batch['jobs']]
            fixture.fixture.completed_opt(first)
            pipeline.prepare_healthy_downstream(fixture.root,SID,SCENARIO)
            early=fixture.root/pipeline.manifest(fixture.root,'optfreq',SID,SCENARIO)[0]['job_dir']
            # The later pre-screen converges into the same basin at lower E.
            model=Chem.MolFromMolBlock((first/'input_structure.mol').read_text(),removeHs=False)
            fixture.edit_completed_output(second,model=model)
            output=second/'job.out';output.write_text(output.read_text().replace('-124.000000','-124.000100'))
            receipt=readj(second/'receipt.json');receipt['output_sha256']['job.out']=sha(output)
            writej(second/'receipt.json',receipt)
            result=pipeline.advance(fixture.root,SID,SCENARIO,'prescreen')
            self.assertEqual(result['status'],'advanced')
            rows=pipeline.manifest(fixture.root,'optfreq',SID,SCENARIO)
            self.assertEqual(len(rows),1)
            self.assertNotEqual(rows[0]['conformer_id'],early.name)
            fixture.fixture.completed_opt(fixture.root/rows[0]['job_dir'])
            audit=pipeline.final_stable_candidates(fixture.root,SID,SCENARIO)
            self.assertTrue(audit['all_prepared_optfreq_pass'],audit['issues'])
            self.assertEqual(audit['superseded_prescreen_representatives_retained_on_disk'],[early.name])
            self.assertTrue((early/'job.json').exists())
        finally:fixture.doCleanups()


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        from test_autoflow_v5 import AutoflowTests
        self.f=AutoflowTests();self.f.setUp();self.addCleanup(self.f.doCleanups)
        self.root=self.f.root
        autoflow.step(self.root,SID)

    def test_review_resume_requeues_without_changing_targets_or_approving(self):
        path=self.root/'automation'/SID/'state.json'
        state=readj(path);state['phase']='REVIEW';writej(path,state)
        writej(self.root/'operations/packed/blocked'/(SID+'.json'),{'cause':'fixture'})
        report=review_control.resume_branch(self.root,SID,'Synthetic repair evidence',execute=True)
        self.assertEqual(report['phase'],'SAMPLING')
        self.assertFalse(report['scientific_convergence_approved'])
        self.assertEqual(readj(path)['targets'],state['targets'])
        with patch.object(packed_worker,'ROOT',self.root),patch.object(packed_worker,'check_storage',return_value={'free_bytes':10000*1024**3}):
            self.assertIsNotNone(packed_worker.claim(self.root/'operations/packed',runtime.site_config(),[SID]))

    def test_resume_rejects_active_claim(self):
        writej(self.root/'operations/packed/claims'/(SID+'.json'),{'token':'active'})
        with self.assertRaisesRegex(ValueError,'claimed'):review_control.resume_branch(self.root,SID,'repair',True)

    def test_resume_rejects_changed_frozen_protocol(self):
        writej(self.root/'config/protocol.json',{})
        with self.assertRaisesRegex(ValueError,'Frozen'):review_control.resume_branch(self.root,SID,'repair',True)

    def test_resume_dry_run_does_not_mutate_state(self):
        path=self.root/'automation'/SID/'state.json';before=sha(path)
        review_control.resume_branch(self.root,SID,'repair')
        self.assertEqual(sha(path),before)

    def test_shortfall_requires_verified_pool_and_bound_review(self):
        evidence=self.root/'review_evidence.txt';evidence.write_text('SYNTHETIC TEST ONLY')
        current={'species_id':SID,'scenario_id':SCENARIO,'workspace':'.','protocol_hash':'test',
                 'original_target':3,'actual_stable_count':2,'stable_audit_sha256':'audit',
                 'sampling_sha256':'sampling','all_pool_goat_verified':True,'all_final_jobs_pass':True,
                 'checks':{'rediscovery_and_motifs_reviewed':False}}
        review={**current,'status':'APPROVED_LIMITED_SHORTFALL_SCOPE','reviewer':'Synthetic fixture reviewer',
                'reason':'Synthetic unit fixture','limitations':'Not real scientific evidence',
                'checks':{'rediscovery_and_motifs_reviewed':True},
                'evidence':[{'path':evidence.name,'sha256':sha(evidence)}]}
        path=self.root/'review.json';writej(path,review)
        with patch.object(review_control,'draft_shortfall',return_value=current):
            report=review_control.install_shortfall(self.root,SID,SCENARIO,path,execute=True)
            self.assertEqual(report['original_target'],3)
            self.assertEqual(review_control.target_for(self.root,SID,SCENARIO,3)[0],2)
            evidence.write_text('tampered')
            with self.assertRaisesRegex(ValueError,'evidence changed'):review_control.target_for(self.root,SID,SCENARIO,3)

    def test_shortfall_cannot_approve_unfinished_sampling(self):
        current={'species_id':SID,'scenario_id':SCENARIO,'workspace':'.','protocol_hash':'test',
                 'original_target':3,'actual_stable_count':2,'stable_audit_sha256':'audit',
                 'sampling_sha256':'sampling','all_pool_goat_verified':False,'all_final_jobs_pass':True}
        with patch.object(review_control,'draft_shortfall',return_value=current):
            with self.assertRaisesRegex(ValueError,'verified sampling'):
                review_control._validate_shortfall(self.root,SID,SCENARIO,{**current,'status':'APPROVED_LIMITED_SHORTFALL_SCOPE'})


if __name__=='__main__':unittest.main()
