"""Mock-only portable staging tests; no external executable is ever launched."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import goat_portable as portable
from scripts import goat_job as runner


class GoatPortableTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='goat_portable_')
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.root = self.base / '中文原始包'
        self.job = self.root / 'goat' / 'j1'
        self.job.mkdir(parents=True)
        self.scratch = self.base / 'scratch'
        for name in portable.IMMUTABLE_FILES:
            (self.job / name).write_text('{}', encoding='utf-8')
        (self.job / 'goat.inp').write_text('! GOAT-ENTROPY XTB ALPB(DMSO)\n%pal nprocs 8 end\n%maxcore 2000\n* xyzfile 0 1 start.xyz\n', encoding='utf-8')
        (self.job / 'start.xyz').write_text('1\nSYNTHETIC UNIT TEST ONLY\nHe 0 0 0\n', encoding='utf-8')
        self.row = {'job_id': 'j1', 'job_directory': 'goat/j1', 'initial_release': True, 'status': 'PREPARED_NOT_RUN',
                    **{field: runner.digest(self.job / name) for name, field in runner.INPUT_HASH_FIELDS.items()}}
        (self.root / 'registry').mkdir()
        runner.atomic_json(self.root / 'registry/jobs.json', [self.row])
        self.pre = patch.object(runner, 'preflight', return_value={'passed': True, 'reasons': [], 'runtime_goat_smoke': 'UNVALIDATED_NOT_EXECUTED'})
        self.pre.start()
        self.addCleanup(self.pre.stop)
        self.proc = patch('subprocess.Popen', side_effect=AssertionError('No actual process may launch'))
        self.proc.start()
        self.addCleanup(self.proc.stop)
        self.env = patch.dict(os.environ, {}, clear=True)
        self.env.start()
        self.addCleanup(self.env.stop)

    def call(self, dry=False, scratch=None):
        return portable.run_one(self.root, 'j1', Path('unused_orca'), Path('unused_xtb'), scratch or self.scratch,
                                acknowledge_models=True, dry_run=dry)

    def mock_run(self, stage_root, job_id, *args, **kwargs):
        stage = stage_root / 'job'
        self.assertTrue(str(stage).isascii())
        self.assertFalse(any(c.isspace() for c in str(stage)))
        self.assertEqual(runner.read_json(stage_root / 'registry/jobs.json')[0]['job_directory'], 'job')
        (stage / 'goat.out').write_text('SYNTHETIC UNIT TEST ONLY', encoding='utf-8')
        (stage / 'goat.finalensemble.xyz').write_text('1\nSYNTHETIC UNIT TEST ONLY\nHe 0 0 0\n', encoding='utf-8')
        (stage / 'worker').mkdir()
        (stage / 'worker' / 'evidence.txt').write_text('synthetic worker evidence', encoding='utf-8')
        value = {'job_id': job_id, 'status': 'COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW',
                 'job_directory': str(stage), 'manifest_sha256': runner.digest(stage_root / 'registry/jobs.json'),
                 'quantum_jobs_launched': 0,
                 'input_sha256': {name: runner.digest(stage / name) for name in portable.IMMUTABLE_FILES},
                 'output_sha256': {name: runner.digest(stage / name) for name in ('goat.out', 'goat.finalensemble.xyz')}}
        runner.atomic_json(stage / 'receipt.json', value)
        return value

    def test_dry_run_creates_nothing_and_never_calls_runner(self):
        before = {str(p): p.read_bytes() for p in self.root.rglob('*') if p.is_file()}
        with patch.object(runner, 'run_one', side_effect=AssertionError('dry-run must not launch')):
            result = self.call(dry=True)
        self.assertEqual(result['quantum_jobs_launched'], 0)
        self.assertFalse(self.scratch.exists())
        self.assertEqual(before, {str(p): p.read_bytes() for p in self.root.rglob('*') if p.is_file()})

    def test_reviewed_retry_retains_old_scratch_and_uses_new_identity(self):
        from clusterflow.common import writej, sha
        first = self.call(dry=True)
        old = Path(first['stage_job_directory']); old.mkdir(parents=True)
        writej(old/'receipt.json', {'status':'RUNNING'})
        proof = self.job/'recovery_history/a/proof.json'
        writej(proof, {'job_id':'j1','scheduler_state':'NODE_FAIL','stage':str(old),
                      'receipt_sha256':sha(old/'receipt.json')})
        writej(self.job/'goat_retry_authorization.json',
               {'status':'PREPARED_VERIFIED_NODE_FAIL','job_id':'j1','attempt_token':'b'*32,
                'manifest_sha256':runner.digest(self.root/'registry/jobs.json'),
                'inputs':portable._hashes(self.job),'proof_path':str(proof.relative_to(self.root)),
                'proof_sha256':sha(proof)})
        second = self.call(dry=True)
        self.assertNotEqual(first['stage_root'],second['stage_root'])
        self.assertEqual(runner.read_json(old/'receipt.json')['status'],'RUNNING')
        self.assertFalse(Path(second['stage_root']).exists())

    def test_invalid_scratch_names_rejected_without_creation(self):
        for value in (self.base / 'with space', self.base / '中文'):
            with self.assertRaisesRegex(ValueError, 'ASCII'):
                self.call(dry=True, scratch=value)
            self.assertFalse(value.exists())

    def test_mocked_stageout_retains_six_inputs_and_receipt_unchanged(self):
        original = portable._hashes(self.job)
        with patch.object(runner, 'run_one', side_effect=self.mock_run) as launch:
            result = self.call()
        self.assertEqual(launch.call_count, 1)
        stage = Path(result['stage_job_directory'])
        self.assertEqual((stage / 'receipt.json').read_bytes(), (self.job / 'receipt.json').read_bytes())
        self.assertTrue((self.job / 'worker/evidence.txt').is_file())
        self.assertEqual(original, portable._hashes(self.job))
        evidence = runner.read_json(self.job / 'portable_stageout.json')
        self.assertEqual(evidence['original_manifest_sha256'], runner.digest(self.root / 'registry/jobs.json'))
        self.assertTrue(stage.is_dir())
        self.assertFalse((self.job / '.goat-run.lock').exists())

    def test_changed_original_inputs_prevent_stageout(self):
        def changed(*args, **kwargs):
            value = self.mock_run(*args, **kwargs)
            (self.job / 'isotope_map.json').write_text('tampered', encoding='utf-8')
            return value
        with patch.object(runner, 'run_one', side_effect=changed):
            with self.assertRaisesRegex(ValueError, 'immutable inputs changed'):
                self.call()
        self.assertFalse((self.job / 'receipt.json').exists())
        self.assertFalse((self.job / 'goat.out').exists())
        self.assertTrue(list(self.scratch.glob('work2goat_*')))

    def test_missing_receipt_is_not_invented_and_no_retry(self):
        with patch.object(runner, 'run_one', return_value={'status': 'FAILED'}):
            with self.assertRaisesRegex(ValueError, 'did not produce receipt'):
                self.call()
        self.assertFalse((self.job / 'receipt.json').exists())
        with self.assertRaisesRegex(ValueError, 'Scratch leaf already exists'):
            self.call(dry=True)

    def test_existing_original_output_never_overwritten(self):
        (self.job / 'goat.out').write_text('keep me', encoding='utf-8')
        with patch.object(runner, 'run_one', side_effect=AssertionError('must block before runner')):
            with self.assertRaisesRegex(ValueError, 'Existing output'):
                self.call()
        self.assertEqual((self.job / 'goat.out').read_text(), 'keep me')


if __name__ == '__main__':
    unittest.main()
