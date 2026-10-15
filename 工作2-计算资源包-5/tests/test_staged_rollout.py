import copy
import os
from pathlib import Path
import tempfile
import subprocess
import unittest
from unittest.mock import patch

from clusterflow.common import readj, writej, sha
from clusterflow import priority, stage_release, scheduling, autoflow
from scripts import submit_cohort, packed_worker

ROOT=Path(__file__).resolve().parents[1]


class StagedRolloutTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)

    def test_scopes_partition_all_models_and_fit_node_budget(self):
        pilot=set(packed_worker.selected_species(ROOT,ROOT/'config/pilot_scope_20261014.json'))
        first=set(packed_worker.selected_species(ROOT,ROOT/'config/complex_first_scope_20261015.json'))
        remaining=set(packed_worker.selected_species(ROOT,ROOT/'config/remaining_scope_20261015.json'))
        all_ids={r['species_id'] for r in readj(ROOT/'registry/batch_plan_v5.json')['species']}
        self.assertEqual((len(pilot),len(first),len(remaining)),(5,8,54))
        self.assertFalse(pilot&first or first&remaining or pilot&remaining)
        self.assertEqual(pilot|first|remaining,all_ids)
        policy=readj(ROOT/'config/rollout_20261015.json')
        self.assertEqual(2+sum(c['nodes'] for c in policy['cohorts']),12)
        classes={r['species_id']:r['complexity_class'] for r in readj(ROOT/'registry/batch_plan_v5.json')['species']}
        self.assertEqual({classes[s] for s in first},{'X'})

    def test_overlap_and_overbudget_refused_even_if_idle(self):
        active=[{'job_id':'123','nodes':2,'species_ids':['a'],'state':'RUNNING'}]
        submit_cohort.validate_expansion(active,['b'],10,12)
        for species,nodes in [(['a'],1),(['b'],11)]:
            with self.assertRaises(ValueError):submit_cohort.validate_expansion(active,species,nodes,12)

    def test_live_queue_covers_accounting_delay(self):
        with patch.object(submit_cohort.subprocess,'run',return_value=subprocess.CompletedProcess([],0,'123|RUNNING\n','')) as run:
            self.assertEqual(submit_cohort.scheduler_state('123'),'RUNNING')
        self.assertEqual(run.call_count,1)
        self.assertEqual(run.call_args.args[0][0],'squeue')

    def test_missing_queue_requires_accounting_evidence_of_terminal(self):
        replies=[subprocess.CompletedProcess([],1,'','Invalid job id'),
                 subprocess.CompletedProcess([],0,'123|COMPLETED|\n','')]
        with patch.object(submit_cohort.subprocess,'run',side_effect=replies):
            self.assertEqual(submit_cohort.scheduler_state('123'),'COMPLETED')

    def test_absent_from_both_scheduler_sources_stays_unknown(self):
        with patch.object(submit_cohort.subprocess,'run',return_value=subprocess.CompletedProcess([],0,'','')):
            with self.assertRaises(ValueError):submit_cohort.scheduler_state('123')

    def test_unknown_submission_journal_blocks_expansion(self):
        writej(self.root/'operations/submission_uncertain.json',{'status':'SUBMISSION_INTENT','jobs':[]})
        with self.assertRaises(ValueError):submit_cohort.active_inventory(self.root)

    def test_changed_active_scope_blocks_expansion(self):
        scope=self.root/'scope.json';writej(scope,{'species_ids':['a']})
        writej(self.root/'operations/submission_old.json',{'status':'SUBMITTED_NOT_RUNNING','jobs':['123'],
               'scope_file':'scope.json','scope_sha256':'wrong','nodes_per_parent':2})
        with patch.object(submit_cohort,'scheduler_state',return_value='RUNNING'):
            with self.assertRaises(ValueError):submit_cohort.active_inventory(self.root)

    def test_stage_gate_is_opt_in_and_does_not_block_goat(self):
        writej(self.root/'registry/batch_plan_v5.json',{'species':[]});writej(self.root/'config/protocol.json',{})
        with patch.dict(os.environ,{'WORK2_STAGED_RELEASE':'1'},clear=True):
            scheduling.admit_quantum(self.root,'goat','a')
            with self.assertRaises(scheduling.DeferredWork) as caught:scheduling.admit_quantum(self.root,'optfreq','a')
            self.assertEqual(caught.exception.report['status'],'WAITING_STAGE_VALIDATION')
        with patch.dict(os.environ,{},clear=True):scheduling.admit_quantum(self.root,'optfreq','a')

    def test_release_requires_unchanged_receipt_not_only_allowlist(self):
        writej(self.root/'config/protocol.json',{});writej(self.root/'registry/batch_plan_v5.json',{'species':[]})
        receipt=self.root/'proof/receipt.json';writej(receipt,{'status':'completed'})
        writej(self.root/'operations/stage_release.json',{'allowed_stages':['optfreq'],
               'evidence':{'optfreq':{'receipt':'proof/receipt.json','sha256':sha(receipt)}}})
        stage_release.check(self.root,'optfreq')
        writej(receipt,{'status':'failed'})
        with self.assertRaises(scheduling.DeferredWork):stage_release.check(self.root,'optfreq')

    def test_no_goat_acceptance_from_zero_exit_failed_receipt(self):
        writej(self.root/'registry/jobs.json',[{'job_id':'a','job_directory':'goat/a'}])
        writej(self.root/'goat/a/receipt.json',{'status':'FAILED','returncode':0})
        writej(self.root/'config/pilot_scope_20261014.json',{'species_ids':[]})
        with patch('clusterflow.pipeline._goat_source',side_effect=AssertionError('Failed source must not validate')):
            report=stage_release.audit(self.root)
        self.assertFalse(report['validated_goat']);self.assertEqual(report['allowed_stages'],['goat'])

    def test_complete_goat_still_requires_source_validation(self):
        writej(self.root/'registry/jobs.json',[{'job_id':'a','job_directory':'goat/a'}])
        writej(self.root/'goat/a/receipt.json',{'status':'COMPLETED','returncode':0})
        writej(self.root/'config/pilot_scope_20261014.json',{'species_ids':[]})
        with patch('clusterflow.pipeline._goat_source',side_effect=ValueError('Missing ensemble')):
            self.assertFalse(stage_release.audit(self.root)['validated_goat'])
        with patch('clusterflow.pipeline._goat_source',return_value={}):
            self.assertTrue(stage_release.audit(self.root)['validated_goat'])

    def test_failed_fast_attempt_never_shortens_priority_estimate(self):
        row=copy.deepcopy(readj(ROOT/'registry/batch_plan_v5.json')['species'][0]);sid=row['species_id']
        row['goat_start_ids']=['a'];row['initial_goat_starts']=1
        jobs={'a':{'job_directory':'goat/a'}}
        before=priority.estimate(self.root,row,jobs)
        writej(self.root/'goat/a/receipt.json',{'status':'FAILED','elapsed_seconds':1})
        self.assertEqual(priority.estimate(self.root,row,jobs)['remaining_serial_hours_estimate'],before['remaining_serial_hours_estimate'])
        writej(self.root/'goat/a/receipt.json',{'status':'COMPLETED','elapsed_seconds':3600})
        after=priority.estimate(self.root,row,jobs)
        self.assertEqual(after['pending_counts_estimate']['goat'],0)
        self.assertLess(after['remaining_serial_hours_estimate'],before['remaining_serial_hours_estimate'])

    def test_successful_optfreq_reduces_remaining_dft_work(self):
        row=copy.deepcopy(readj(ROOT/'registry/batch_plan_v5.json')['species'][0]);sid=row['species_id']
        workspace='automation/'+sid+'/wave_001/workspace'
        writej(self.root/'automation'/sid/'state.json',{'phase':'DFT','waves':[{'workspace':workspace}]})
        before=priority.estimate(self.root,row,{})
        folder=self.root/workspace/'jobs/hash'/sid/'scenario/optfreq/c001'
        writej(folder/'job.json',{'species_id':sid,'stage':'optfreq'})
        writej(folder/'receipt.json',{'status':'completed','elapsed_seconds':100*3600})
        after=priority.estimate(self.root,row,{})
        self.assertEqual(after['pending_counts_estimate']['optfreq'],before['pending_counts_estimate']['optfreq']-1)
        self.assertEqual(after['per_job_hours_estimate']['optfreq'],100)

    def test_unvalidated_stage_does_not_hide_ready_work(self):
        rows=[]
        for stage in ['shielding','prescreen']:
            writej(self.root/stage/'job.json',{'stage':stage});rows.append({'job_dir':stage})
        def check(root,stage):
            if stage=='shielding':raise scheduling.DeferredWork('WAITING_STAGE_VALIDATION')
        with patch.dict(os.environ,{'WORK2_STAGED_RELEASE':'1'}),patch.object(stage_release,'check',side_effect=check), \
             patch.object(autoflow,'_execute_dft',return_value={'status':'completed'}) as run:
            autoflow._run_pending(self.root,rows,'orca','scratch')
        self.assertEqual(run.call_args.args[1]['job_dir'],'prescreen')

    def test_priority_claim_retains_exclusion_of_active_and_blocked_branches(self):
        site=readj(ROOT/'config/site_scnet.json');base=self.root/'operations/packed'
        writej(base/'blocked/b.json',{'retry_on_different_node':False})
        writej(base/'claims/a.json',{'species_id':'a','growth_gib':64,'slurm_job_id':'1'})
        with patch.object(packed_worker,'ROOT',self.root),patch.object(packed_worker,'scheduler_terminal',return_value=False), \
             patch.object(packed_worker,'check_storage',return_value={'free_bytes':10000*1024**3}), \
             patch.object(priority,'ranked_species',return_value=['a','b','c']), \
             patch.dict(os.environ,{'WORK2_PRIORITY_POLICY':'remaining_work'}):
            _,record=packed_worker.claim(base,site,['a','b','c'])
        self.assertEqual(record['species_id'],'c')

    def cohort_fixture(self):
        site=readj(ROOT/'config/site_scnet.json')
        writej(self.root/'config/site_scnet.json',site)
        ids=['a','b','c','d']
        writej(self.root/'registry/batch_plan_v5.json',{'species':[{'species_id':s} for s in ids]})
        scope=self.root/'scope.json';writej(scope,{'species_ids':ids})
        runtime=self.root/'runtime.json';writej(runtime,{'version':'6.1.0','site_sha256':sha(self.root/'config/site_scnet.json')})
        return scope,runtime

    def test_no_bulk_submission_without_real_goat(self):
        scope,runtime=self.cohort_fixture()
        with patch.dict(os.environ,{'WORK2_RUNTIME_EVIDENCE':str(runtime)}), \
             patch.object(stage_release,'audit',return_value={'validated_goat':False}), \
             patch.object(submit_cohort.subprocess,'run',side_effect=AssertionError('No scheduler call permitted')):
            result=submit_cohort.submit(self.root,'test',scope,1,execute=True,require_goat=True)
        self.assertEqual(result['status'],'WAITING_VALID_GOAT')
        self.assertFalse((self.root/'operations/submission_cohort_test.json').exists())

    def test_test_only_failure_never_submits_or_records_intent(self):
        scope,runtime=self.cohort_fixture()
        with patch.dict(os.environ,{'WORK2_RUNTIME_EVIDENCE':str(runtime)}), \
             patch.object(submit_cohort,'check_storage',return_value={'free_bytes':10000*1024**3}), \
             patch.object(submit_cohort.subprocess,'run',side_effect=subprocess.CalledProcessError(1,'sbatch')) as run:
            with self.assertRaises(subprocess.CalledProcessError):submit_cohort.submit(self.root,'test',scope,1,execute=True)
        self.assertEqual(run.call_count,1)
        self.assertIn('--test-only',run.call_args.args[0])
        self.assertFalse((self.root/'operations/submission_cohort_test.json').exists())

    def test_same_cohort_reuses_submission_not_quantum_work(self):
        scope,runtime=self.cohort_fixture()
        returns=[subprocess.CompletedProcess([],0,'','test passed'),subprocess.CompletedProcess([],0,'',''),
                 subprocess.CompletedProcess([],0,'12345\n','')]
        with patch.dict(os.environ,{'WORK2_RUNTIME_EVIDENCE':str(runtime)}), \
             patch.object(submit_cohort,'check_storage',return_value={'free_bytes':10000*1024**3}), \
             patch.object(submit_cohort.subprocess,'run',side_effect=returns) as run:
            result=submit_cohort.submit(self.root,'test',scope,1,execute=True)
            self.assertEqual(result['jobs'],['12345']);self.assertEqual(run.call_count,3)
        with patch.dict(os.environ,{'WORK2_RUNTIME_EVIDENCE':str(runtime)}), \
             patch.object(submit_cohort,'scheduler_state',return_value='RUNNING'), \
             patch.object(submit_cohort.subprocess,'run',side_effect=AssertionError('Duplicate submission')):
            self.assertEqual(submit_cohort.submit(self.root,'test',scope,1,execute=True)['status'],'ALREADY_SUBMITTED')


if __name__=='__main__':unittest.main()
