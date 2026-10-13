"""Synthetic review fixtures only: no ORCA, no real approval or production data."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_cluster_pipeline as fixtures
from clusterflow.common import writej, readj, sha
from clusterflow import gates, issue_review as ir


class IssueReviewTests(unittest.TestCase):
    def setUp(self):
        self.f = fixtures.ClusterPipelineTests()
        self.f.setUp()
        self.addCleanup(self.f.doCleanups)
        self.f.cfg['scientific_issue_policy'] = 'user_review'
        writej(self.f.root / 'config/protocol.json', self.f.cfg)
        prepared = self.f.prepared()
        self.job = self.f.root / prepared['jobs'][0]

    def complete(self, negative=True):
        self.f.completed_opt(self.job, negative=negative)

    def approve_fixture(self):
        evidence = self.job / 'SYNTHETIC_UNIT_TEST_MODE_REVIEW.txt'
        evidence.write_text('SYNTHETIC UNIT TEST ONLY. Not a scientific review.', encoding='utf-8')
        review = ir.draft_scientific_review(self.job, self.f.cfg)
        review.update(status='ACCEPT_AS_REVIEWED_APPROXIMATION',
                      reviewer='SYNTHETIC_UNIT_TEST_NOT_REAL_APPROVAL',
                      reviewed_utc='2026-10-14T00:00:00Z',
                      rationale='Synthetic unit fixture only; not scientific evidence.',
                      limitations='Synthetic unit fixture only; does not prove a physical minimum.',
                      checks={key: True for key in ir.REVIEW_CHECKS},
                      mode_review_evidence=[{'path': evidence.name, 'sha256': sha(evidence),
                                             'description': 'Synthetic test evidence only'}])
        path = self.job / 'scientific_review.json'
        writej(path, review)
        return path

    def rewrite_receipt(self):
        receipt = readj(self.job / 'receipt.json')
        receipt['output_sha256']['job.out'] = sha(self.job / 'job.out')
        writej(self.job / 'receipt.json', receipt)

    def test_missing_output_has_diagnosis_but_not_approval(self):
        diagnosis = ir.diagnose_job(self.job, self.f.cfg)
        self.assertFalse(diagnosis['normal_termination_observed'])
        self.assertEqual(diagnosis['findings'][0]['code'], 'not_run_or_output_missing')
        self.assertFalse(gates.inspect_job(self.job, self.f.cfg)['accepted'])

    def test_negative_modes_without_review_are_pending_not_autoaccepted(self):
        self.complete()
        inspected = gates.inspect_job(self.job, self.f.cfg)
        self.assertFalse(inspected['accepted'])
        self.assertTrue(inspected['needs_review'])
        self.assertIn('scientific_review_missing', inspected['errors'])

    def test_draft_never_approves(self):
        self.complete()
        writej(self.job / 'scientific_review.json', ir.draft_scientific_review(self.job, self.f.cfg))
        inspected = gates.inspect_job(self.job, self.f.cfg)
        self.assertFalse(inspected['accepted'])
        self.assertIn('scientific_review_not_approved', inspected['errors'])

    def test_source_bound_review_accepts_approximation_not_minimum_claim(self):
        self.complete()
        self.approve_fixture()
        inspected = gates.inspect_job(self.job, self.f.cfg)
        self.assertTrue(inspected['accepted'], inspected['errors'])
        self.assertEqual(inspected['scientific_status'], 'HUMAN_REVIEWED_APPROXIMATION')
        self.assertFalse(inspected['scientific_review']['physical_minimum_certified'])

    def test_no_arbitrary_frequency_threshold_approves_small_negative(self):
        self.complete()
        output = self.job / 'job.out'
        output.write_text(output.read_text().replace('-20.00', '-0.01'), encoding='utf-8')
        self.rewrite_receipt()
        self.assertFalse(gates.inspect_job(self.job, self.f.cfg)['accepted'])
        self.assertEqual(ir.diagnose_job(self.job, self.f.cfg)['negative_modes'][0]['frequency_cm1'], -.01)

    def test_unconverged_output_cannot_be_overridden(self):
        self.complete()
        output = self.job / 'job.out'
        output.write_text(output.read_text().replace('THE OPTIMIZATION HAS CONVERGED', 'OPTIMIZATION NOT CONVERGED'), encoding='utf-8')
        self.rewrite_receipt()
        self.approve_fixture()
        inspected = gates.inspect_job(self.job, self.f.cfg)
        self.assertFalse(inspected['accepted'])
        self.assertIn('optimization_not_converged', inspected['errors'])

    def test_timeout_is_advisory_and_not_overridable(self):
        self.complete()
        receipt = readj(self.job / 'receipt.json')
        receipt.update(status='failed', timed_out=True, failure_reason='timeout')
        writej(self.job / 'receipt.json', receipt)
        self.approve_fixture()
        self.assertFalse(gates.inspect_job(self.job, self.f.cfg)['accepted'])
        codes = [row['code'] for row in ir.diagnose_job(self.job, self.f.cfg)['findings']]
        self.assertIn('timeout', codes)

    def test_review_invalid_after_source_change_even_rebound_receipt(self):
        self.complete()
        self.approve_fixture()
        output = self.job / 'job.out'
        output.write_text(output.read_text() + '\nSOURCE MODIFIED\n', encoding='utf-8')
        self.rewrite_receipt()
        self.assertIn('scientific_review_source_hash_mismatch', gates.inspect_job(self.job, self.f.cfg)['errors'])

    def test_mode_evidence_tampering_rejected(self):
        self.complete()
        self.approve_fixture()
        (self.job / 'SYNTHETIC_UNIT_TEST_MODE_REVIEW.txt').write_text('changed', encoding='utf-8')
        self.assertIn('scientific_review_mode_evidence_invalid', gates.inspect_job(self.job, self.f.cfg)['errors'])

    def test_software_name_not_human_approval(self):
        self.complete()
        path = self.approve_fixture()
        review = readj(path)
        review['reviewer'] = 'codex'
        writej(path, review)
        self.assertIn('scientific_review_requires_named_human', gates.inspect_job(self.job, self.f.cfg)['errors'])

    def test_incomplete_frequencies_not_overridable(self):
        self.complete()
        output = self.job / 'job.out'
        lines = output.read_text().splitlines()
        freq_indices = [i for i, line in enumerate(lines) if 'cm**-1' in line]
        del lines[freq_indices[-1]]
        output.write_text('\n'.join(lines), encoding='utf-8')
        self.rewrite_receipt()
        self.approve_fixture()
        self.assertIn('frequency_count_not_3N', gates.inspect_job(self.job, self.f.cfg)['errors'])

    def test_clean_job_needs_no_scientific_exception(self):
        self.complete(negative=False)
        self.assertTrue(gates.inspect_job(self.job, self.f.cfg)['accepted'])

    def test_legacy_strict_mode_never_uses_exception(self):
        self.f.cfg['scientific_issue_policy'] = 'strict'
        # Reprepare in a new fixture to bind the strict protocol hash correctly.
        self.f.doCleanups()
        self.f = fixtures.ClusterPipelineTests()
        self.f.setUp()
        self.addCleanup(self.f.doCleanups)
        self.f.cfg['scientific_issue_policy'] = 'strict'
        writej(self.f.root / 'config/protocol.json', self.f.cfg)
        self.job = self.f.root / self.f.prepared()['jobs'][0]
        self.complete()
        self.approve_fixture()
        inspected = gates.inspect_job(self.job, self.f.cfg)
        self.assertFalse(inspected['accepted'])
        self.assertIn('imaginary_frequencies_require_mode_review', inspected['errors'])


if __name__ == '__main__':
    unittest.main()
