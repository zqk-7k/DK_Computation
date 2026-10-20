import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch
from clusterflow.common import readj
from scripts.goat_bulk_worker import allocation_deadline


class ArrayDeadlineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)

    def query(self, output):
        with patch.dict(os.environ, {'SLURM_JOB_ID': '120'}), patch(
                'scripts.goat_bulk_worker.subprocess.run', return_value=Mock(stdout=output)):
            return allocation_deadline(self.base)

    def test_array_master_line_before_siblings(self):
        self.query('JobId=120 JobState=RUNNING EndTime=2099-01-01T00:00:00 ArrayTaskId=4\n'
                   'JobId=121 JobState=RUNNING EndTime=2099-02-01T00:00:00 ArrayTaskId=0\n')
        fields = readj(self.base/'allocations/120/deadline.json')['scheduler_fields']
        self.assertEqual(fields['ArrayTaskId'], '4')
        self.assertEqual(fields['EndTime'], '2099-01-01T00:00:00')

    def test_array_master_line_after_siblings(self):
        self.query('JobId=121 JobState=FAILED EndTime=2020-01-01T00:00:00\n'
                   'JobId=120 JobState=RUNNING EndTime=2099-01-01T00:00:00\n')

    def test_running_sibling_does_not_authorize_failed_master(self):
        with self.assertRaisesRegex(ValueError, 'requested running parent'):
            self.query('JobId=120 JobState=FAILED EndTime=2099-01-01T00:00:00\n'
                       'JobId=121 JobState=RUNNING EndTime=2099-01-01T00:00:00\n')

    def test_future_sibling_cannot_hide_expired_actual_deadline(self):
        with self.assertRaisesRegex(ValueError, 'deadline already passed'):
            self.query('JobId=120 JobState=RUNNING EndTime=2020-01-01T00:00:00\n'
                       'JobId=121 JobState=RUNNING EndTime=2099-01-01T00:00:00\n')

    def test_ambiguous_duplicate_identity_rejected(self):
        line = 'JobId=120 JobState=RUNNING EndTime=2099-01-01T00:00:00\n'
        with self.assertRaisesRegex(ValueError, 'requested running parent'):
            self.query(line+line)
