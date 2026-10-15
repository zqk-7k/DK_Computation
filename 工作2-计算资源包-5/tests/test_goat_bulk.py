import multiprocessing
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch, Mock

from clusterflow.common import readj, writej, sha
from clusterflow.runtime import StorageAdmissionError
from scripts.goat_bulk_worker import claim, finish, verify_scope, allocation_deadline
from scripts.submit_goat_bulk import build_scope, IDENTITIES
from scripts.recover_goat_bulk_startup import recover
from scripts.goat_bulk_parent_guard import supervisor


def concurrent_claim(root, base, policy, indexed, results):
    result = claim(Path(root), Path(base), policy, indexed, storage_check=lambda *a: None)
    results.put(result['job_id'] if result else None)


class BulkTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.base = self.root/'operations/goat_bulk_20261015'
        self.rows = [{'job_id': 'j'+str(i), 'species_id': 's'+str(i//7),
                      'job_directory': 'jobs/j'+str(i), 'initial_release': i % 7 < 6}
                     for i in range(14)]
        self.indexed = {r['job_id']: r for r in self.rows}
        for row in self.rows:
            (self.root/row['job_directory']).mkdir(parents=True)
        plans = [{'species_id': 's'+str(i), 'compound_id': 'c'+str(i),
                  'complexity_class': 'X', 'initial_goat_starts': 6,
                  'initial_candidate_batch_size': 30,
                  'goat_start_ids': ['j'+str(j) for j in range(i*7, (i+1)*7)]}
                 for i in range(2)]
        for name in IDENTITIES:
            writej(self.root/name, {})
        writej(self.root/'registry/jobs.json', self.rows)
        writej(self.root/'registry/batch_plan_v5.json', {'species': plans})
        self.policy = build_scope(self.root)

    def get(self):
        return claim(self.root, self.base, self.policy, self.indexed, storage_check=lambda *a: None)

    def test_initial_starts_only(self):
        self.assertEqual(len(self.policy['job_ids']), 12)
        self.assertNotIn('j6', self.policy['job_ids'])
        self.assertNotIn('j13', self.policy['job_ids'])

    def test_occupied_and_blocked_are_excluded(self):
        writej(self.root/'jobs/j0/.goat-run.lock', {'pid': 1})
        writej(self.root/'jobs/j1/receipt.json', {'status': 'FAILED'})
        writej(self.root/'operations/packed/blocked/s1.json', {})
        policy = build_scope(self.root)
        self.assertEqual(policy['legacy_active_jobs'], 1)
        self.assertEqual(set(policy['job_ids']), {'j2', 'j3', 'j4', 'j5'})

    def test_concurrent_claims_unique(self):
        ctx = multiprocessing.get_context('fork')
        results = ctx.Queue()
        children = [ctx.Process(target=concurrent_claim, args=(str(self.root), str(self.base), self.policy, self.indexed, results)) for _ in range(16)]
        for child in children: child.start()
        for child in children:
            child.join(10)
            self.assertEqual(child.exitcode, 0)
        claimed = [results.get(timeout=2) for _ in children]
        self.assertEqual(len(set(c for c in claimed if c)), 12)
        self.assertEqual(sum(c is None for c in claimed), 4)

    def test_finish_failure_blocks_unstarted_siblings(self):
        record = self.get()
        finish(self.root, self.base, record, {'status': 'FAILED'})
        next_record = self.get()
        self.assertNotEqual(record['species_id'], next_record['species_id'])
        self.assertNotIn(record['job_id'], readj(self.base/'queue.json')['active'])

    def test_success_is_not_scientific_acceptance(self):
        record = self.get()
        finish(self.root, self.base, record, {'status': 'COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW'})
        self.assertFalse((self.base/'blocked'/('s0.json')).exists())
        self.assertNotIn('scientific_acceptance', readj(self.base/'results'/(record['job_id']+'.json')))

    def test_pause_does_not_claim(self):
        writej(self.base/'PAUSE_NEW_CLAIMS', {})
        self.assertIsNone(self.get())
        self.assertFalse((self.base/'queue.json').exists())

    def test_source_lock_created_after_snapshot_skipped(self):
        writej(self.root/'jobs/j0/.goat-run.lock', {'pid': 1})
        self.assertNotEqual(self.get()['job_id'], 'j0')

    def test_unknown_stageout_attempt_never_restarted(self):
        writej(self.root/'jobs/j0/stage_intents/intent.json', {})
        self.assertNotEqual(self.get()['job_id'], 'j0')

    def test_storage_growth_reserves_legacy_and_active(self):
        reserve = []
        self.policy['legacy_active_jobs'] = 12
        claim(self.root, self.base, self.policy, self.indexed, storage_check=lambda p, n: reserve.append(n))
        claim(self.root, self.base, self.policy, self.indexed, storage_check=lambda p, n: reserve.append(n))
        self.assertEqual(reserve, [128+13*8, 128+14*8])

    def test_low_space_cannot_advance_queue(self):
        def blocked(*args): raise StorageAdmissionError('DISK_ADMISSION_BLOCKED', '/scratch')
        with self.assertRaises(StorageAdmissionError):
            claim(self.root, self.base, self.policy, self.indexed, storage_check=blocked)
        self.assertFalse((self.base/'queue.json').exists())

    def test_frozen_inputs_verified(self):
        path = self.root/'scope.json'; writej(path, self.policy)
        verify_scope(self.root, path, sha(path))
        writej(self.root/'config/protocol.json', {'changed': True})
        with self.assertRaisesRegex(ValueError, 'identity changed'):
            verify_scope(self.root, path, sha(path))

    def test_held_start_cannot_be_injected(self):
        self.policy['job_ids'].append('j6')
        path = self.root/'scope.json'; writej(path, self.policy)
        with self.assertRaisesRegex(ValueError, 'unreleased start'):
            verify_scope(self.root, path, sha(path))

    def test_deadline_queried_once_for_shared_parent(self):
        result = Mock(stdout='JobId=123 JobState=RUNNING EndTime=2099-01-01T00:00:00 ArrayJobId=120 ArrayTaskId=3')
        with patch.dict(os.environ, {'SLURM_JOB_ID': '123'}), patch('scripts.goat_bulk_worker.subprocess.run', return_value=result) as run:
            first = allocation_deadline(self.base)
            self.assertEqual(first, allocation_deadline(self.base))
            self.assertEqual(run.call_count, 1)
            self.assertEqual(readj(self.base/'allocations/123/deadline.json')['array_job_id'], '120')

    def test_deadline_rejects_wrong_allocation(self):
        result = Mock(stdout='JobId=999 JobState=RUNNING EndTime=2099-01-01T00:00:00')
        with patch.dict(os.environ, {'SLURM_JOB_ID': '123'}), patch('scripts.goat_bulk_worker.subprocess.run', return_value=result):
            with self.assertRaisesRegex(ValueError, 'requested running parent'):
                allocation_deadline(self.base)

    def test_recovery_selects_only_missing_nodes_without_submitting(self):
        scope = self.base/'scope.json'; writej(scope, self.policy)
        writej(self.base/'submission.json', {'job_id': '120', 'array_tasks': 1, 'scope_sha256': sha(scope)})
        (self.base/'slurm-120_0.err').write_text('end_time TimeoutExpired')
        for host in ('good', 'missing'):
            writej(self.root/('operations/mpi_probes/123_'+host+'.json'), {'hostname': host, 'status': 'PASS'})
        writej(self.base/'nodes/123_good.json', {'status': 'RUNNING'})
        responses = [Mock(stdout='120_0|RUNNING\n'), Mock(stdout='JobId=123 ArrayJobId=120 ArrayTaskId=0 JobState=RUNNING')]
        with patch('scripts.recover_goat_bulk_startup.subprocess.run', side_effect=responses), patch('scripts.recover_goat_bulk_startup.subprocess.Popen') as spawn:
            report = recover(self.root)
            self.assertEqual(report[0]['nodes'], ['missing'])
            spawn.assert_not_called()

    def test_parent_guard_selects_supervisor_not_quantum_worker(self):
        proc = self.root/'proc'
        (proc/'20').mkdir(parents=True); (proc/'10').mkdir()
        (proc/'20/stat').write_text('20 (python) S 10 0 0')
        (proc/'10/cmdline').write_bytes(b'python\0scripts/goat_bulk_worker.py\0--scope\0x\0')
        self.assertEqual(supervisor({'workers': [{'pid': 20, 'alive': True}]}, proc), 10)
        (proc/'10/cmdline').write_bytes(b'orca\0goat.inp\0')
        with self.assertRaisesRegex(ValueError, 'Refusing to pause'):
            supervisor({'workers': [{'pid': 20, 'alive': True}]}, proc)


if __name__ == '__main__':
    unittest.main()
