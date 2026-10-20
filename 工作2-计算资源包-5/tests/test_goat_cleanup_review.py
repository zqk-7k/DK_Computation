import tempfile
import unittest
from pathlib import Path
from clusterflow.common import readj, writej, sha
from clusterflow.goat_cleanup_review import evidence, validate


class CleanupReviewTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root/'goat.out').write_text('SYNTHETIC TEST ONLY\nProgram Version 6.1.0\nORCA TERMINATED NORMALLY\n')
        (self.root/'goat.finalensemble.xyz').write_text('1\nsynthetic\nH 0 0 0\n')
        (self.root/'goat.stderr').write_text(
            'INFO:    Starting build...\n'
            'rm: cannot remove /var/tmp/work2-runtime-fixture/orca610_sandbox_test/a: Read-only file system\n')
        self.receipt = self.root/'receipt.json'
        self.rebind()

    def rebind(self):
        writej(self.receipt, {'status': 'FAILED', 'returncode': 1,
                             'errors': ['ORCA nonzero or absent exit code'],
                             'output_sha256': {p.name: sha(p) for p in self.root.glob('goat.*')}})

    def test_review_keeps_original_failed_receipt(self):
        before = sha(self.receipt)
        writej(self.root/'goat_cleanup_review.json', {
            'status': 'NORMAL_ORCA_OUTPUT_WITH_FAILED_SANDBOX_CLEANUP',
            'evidence': evidence(self.root, self.receipt)})
        self.assertEqual(validate(self.root, self.receipt), sha(self.root/'goat_cleanup_review.json'))
        self.assertEqual(before, sha(self.receipt))
        self.assertEqual(readj(self.receipt)['status'], 'FAILED')

    def test_changed_output_rejected(self):
        (self.root/'goat.out').write_text('changed')
        with self.assertRaisesRegex(ValueError, 'binding'):
            evidence(self.root, self.receipt)

    def test_quantum_error_rejected_even_with_normal_end(self):
        with (self.root/'goat.out').open('a') as f:
            f.write('GOAT ERROR: No structure was left')
        self.rebind()
        with self.assertRaisesRegex(ValueError, 'normal-output'):
            evidence(self.root, self.receipt)

    def test_other_stderr_errors_rejected(self):
        with (self.root/'goat.stderr').open('a') as f:
            f.write('FATAL I/O ERROR\n')
        self.rebind()
        with self.assertRaisesRegex(ValueError, 'Unreviewed'):
            evidence(self.root, self.receipt)

    def test_other_receipt_errors_rejected(self):
        r = readj(self.receipt)
        r['errors'].append('Normal termination missing')
        writej(self.receipt, r)
        with self.assertRaisesRegex(ValueError, 'isolated'):
            evidence(self.root, self.receipt)
