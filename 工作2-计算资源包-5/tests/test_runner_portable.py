"""Mock-only downstream portable runner tests. No QM executable may launch.

All data below are explicitly synthetic temporary transport fixtures, not
calculation results, scientific acceptance, or benchmark evidence.
"""
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from clusterflow import runner
from clusterflow.common import hash_protocol,readj,sha,writej


class RunnerPortableTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='runner_portable_')
        self.addCleanup(self.tmp.cleanup)
        self.base=Path(self.tmp.name)
        self.root=self.base/'中文原始包'
        self.relative='jobs/synthetic/prescreen/c001'
        self.job=self.root/self.relative
        self.job.mkdir(parents=True)
        self.scratch=self.base/'scratch'
        self.exe=self.base/'orca.exe'
        self.exe.write_bytes(b'SYNTHETIC TEST ONLY Program Version 6.1.1')
        self.cfg={'orca_version':'6.1.1','synthetic_test_fixture':True}
        writej(self.root/'config/protocol.json',self.cfg)
        for name in runner.IMMUTABLE_FILES:
            (self.job/name).write_text('SYNTHETIC UNIT TEST ONLY\n',encoding='utf-8')
        self.metadata={'id':'synthetic-c001','stage':'prescreen','nprocs':1,'maxcore_mb':1,
                       'protocol_hash':hash_protocol(self.cfg),'input_sha256':sha(self.job/'job.inp'),
                       'atom_map_sha256':sha(self.job/'atom_map.json'),'required_outputs':['job.xyz'],
                       'notice':'SYNTHETIC TEST FIXTURE, NOT SCIENTIFIC EVIDENCE'}
        writej(self.job/'job.json',self.metadata)
        self.env=patch.dict(os.environ,{},clear=True)
        self.env.start();self.addCleanup(self.env.stop)
        self.preflight=patch.object(runner,'static_preflight',return_value={
            'passed':True,'issues':[],'orca':str(self.exe),'orca_sha256':sha(self.exe)})
        self.preflight.start();self.addCleanup(self.preflight.stop)
        self.process=patch('subprocess.Popen',side_effect=AssertionError('No QM processes in transport tests'))
        self.process.start();self.addCleanup(self.process.stop)
        self.stage=None

    def call(self,**kwargs):
        args={'acknowledge':True,'scratch_root':self.scratch}
        args.update(kwargs)
        return runner.run_one(self.root,self.relative,self.exe,**args)

    def mock_run(self,directory,exe,**kwargs):
        directory=Path(directory);self.stage=directory
        self.assertTrue(str(directory).isascii())
        self.assertFalse(any(c.isspace() for c in str(directory)))
        self.assertTrue((self.job/'.clusterflow-run.lock').is_file())
        if directory!=self.job:self.assertTrue((self.job/'.orca-run.lock').is_file())
        (directory/'job.out').write_bytes(b'SYNTHETIC TEST ONLY\r\nNOT A QM RESULT\r\n')
        (directory/'job.stderr').write_bytes(b'')
        (directory/'job.xyz').write_bytes(b'1\nSYNTHETIC TEST ONLY\nHe 0 0 0\n')
        (directory/'worker').mkdir(exist_ok=True)
        (directory/'worker'/'evidence.txt').write_bytes(b'SYNTHETIC WORKER EVIDENCE')
        receipt={'id':self.metadata['id'],'stage':self.metadata['stage'],'status':'completed',
                 'input_sha256':runner._inputs(directory),'orca_executable_sha256':sha(exe),
                 'output_sha256':{n:sha(directory/n) for n in ('job.out','job.stderr','job.xyz')},
                 'notice':'SYNTHETIC UNIT TEST ONLY; NOT SCIENTIFIC EVIDENCE'}
        writej(directory/'receipt.json',receipt)
        return receipt

    def successful(self,**kwargs):
        with patch.object(runner,'run_job',side_effect=self.mock_run):return self.call(**kwargs)

    def no_receipt(self):
        self.assertFalse((self.job/'receipt.json').exists())
        self.assertFalse((self.job/'.orca-run.lock').exists())
        self.assertFalse((self.job/'.clusterflow-run.lock').exists())

    def test_dry_run_writes_nothing(self):
        before={p.relative_to(self.base).as_posix():p.read_bytes() for p in self.base.rglob('*') if p.is_file()}
        with patch.object(runner,'run_job',side_effect=AssertionError('dry-run executed')):
            result=self.call(dry_run=True)
        self.assertEqual(result['quantum_jobs_launched'],0)
        self.assertEqual(result['status'],'DRY_RUN_NO_EXECUTION')
        self.assertFalse(self.scratch.exists())
        self.assertEqual(before,{p.relative_to(self.base).as_posix():p.read_bytes() for p in self.base.rglob('*') if p.is_file()})

    def test_chinese_source_requires_ascii_scratch(self):
        with self.assertRaisesRegex(ValueError,'scratch-root'):self.call(scratch_root=None,dry_run=True)
        for path in (self.base/'with space',self.base/'中文'):
            with self.assertRaisesRegex(ValueError,'ASCII'):self.call(scratch_root=path,dry_run=True)
            self.assertFalse(path.exists())

    def test_scratch_drive_root_and_regular_file_rejected(self):
        for path in (Path(self.base.anchor),self.exe):
            with self.assertRaisesRegex(ValueError,'dedicated scratch'):self.call(scratch_root=path,dry_run=True)

    def test_different_packages_do_not_share_scratch_identity(self):
        a=self.call(dry_run=True)['ascii_scratch_directory']
        other=self.base/'another_package';shutil.copytree(self.root,other)
        b=runner.run_one(other,self.relative,self.exe,acknowledge=True,dry_run=True,
                        scratch_root=self.scratch)['ascii_scratch_directory']
        self.assertNotEqual(a,b)

    def test_receipt_exact_bytes_outputs_and_inputs_preserved(self):
        inputs=runner._inputs(self.job)
        value=self.successful()
        self.assertEqual(value,readj(self.job/'receipt.json'))
        self.assertEqual((self.stage/'receipt.json').read_bytes(),(self.job/'receipt.json').read_bytes())
        self.assertEqual(inputs,runner._inputs(self.job))
        evidence=readj(self.job/'portable_stageout.json')
        self.assertEqual(evidence['staged_receipt_sha256'],sha(self.job/'receipt.json'))
        self.assertEqual(evidence['original_input_sha256'],inputs)
        self.assertIn('worker/evidence.txt',evidence['copied_output_sha256'])
        for name,digest in evidence['copied_output_sha256'].items():
            self.assertEqual(sha(self.job/name),digest)
        self.assertTrue(self.stage.is_dir())
        self.assertFalse((self.job/'.orca-run.lock').exists())

    def test_changed_original_each_input_rejects_stageout(self):
        # Each immutable input gets its own independent job/scratch attempt.
        for name in runner.IMMUTABLE_FILES:
            with self.subTest(name=name):
                saved=(self.job/name).read_bytes()
                def changed(*args,**kwargs):
                    value=self.mock_run(*args,**kwargs)
                    (self.job/name).write_bytes(saved+b'changed')
                    return value
                with patch.object(runner,'run_job',side_effect=changed):
                    with self.assertRaises((ValueError,json.JSONDecodeError)):
                        self.call(scratch_root=self.base/('scratch_'+name.replace('.','_')))
                self.no_receipt()
                self.assertFalse((self.job/'job.out').exists())
                (self.job/name).write_bytes(saved)

    def test_staged_input_mutation_rejected(self):
        def changed(*args,**kwargs):
            value=self.mock_run(*args,**kwargs)
            (self.stage/'input.xyz').write_bytes(b'changed after execution')
            return value
        with patch.object(runner,'run_job',side_effect=changed):
            with self.assertRaisesRegex(ValueError,'Staged immutable'):self.call()
        self.no_receipt()

    def test_changed_protocol_and_executable_rejected(self):
        for target in ('protocol','executable'):
            with self.subTest(target=target):
                def changed(*args,**kwargs):
                    value=self.mock_run(*args,**kwargs)
                    if target=='protocol':writej(self.root/'config/protocol.json',{'changed':True})
                    else:self.exe.write_bytes(b'changed executable')
                    return value
                with patch.object(runner,'run_job',side_effect=changed):
                    with self.assertRaises(ValueError):self.call(scratch_root=self.base/target)
                self.no_receipt()
                writej(self.root/'config/protocol.json',self.cfg)
                self.exe.write_bytes(b'SYNTHETIC TEST ONLY Program Version 6.1.1')

    def test_missing_on_disk_receipt_not_invented(self):
        with patch.object(runner,'run_job',return_value={'status':'completed'}):
            with self.assertRaisesRegex(ValueError,'Missing regular'):self.call()
        self.no_receipt()

    def test_returned_receipt_must_match_disk(self):
        def changed(*args,**kwargs):
            value=self.mock_run(*args,**kwargs)
            return {**value,'different':True}
        with patch.object(runner,'run_job',side_effect=changed):
            with self.assertRaisesRegex(ValueError,'Returned receipt differs'):self.call()
        self.no_receipt()

    def test_receipt_identity_binding_rejects_forgery(self):
        for field,bad in [('id','another-job'),('stage','shielding'),('status','running'),
                          ('input_sha256',{}),('orca_executable_sha256','0'*64),('output_sha256',{})]:
            with self.subTest(field=field):
                def changed(*args,**kwargs):
                    value=self.mock_run(*args,**kwargs);value[field]=bad
                    writej(self.stage/'receipt.json',value);return value
                with patch.object(runner,'run_job',side_effect=changed):
                    with self.assertRaises(ValueError):self.call(scratch_root=self.base/field)
                self.no_receipt()

    def test_receipt_output_modified_after_hashing_rejected(self):
        def changed(*args,**kwargs):
            value=self.mock_run(*args,**kwargs)
            (self.stage/'job.out').write_bytes(b'tampered output')
            return value
        with patch.object(runner,'run_job',side_effect=changed):
            with self.assertRaisesRegex(ValueError,'output hash mismatch'):self.call()
        self.no_receipt()

    def test_receipt_output_traversal_rejected(self):
        def changed(*args,**kwargs):
            value=self.mock_run(*args,**kwargs)
            outside=self.stage.parent/'outside.txt';outside.write_bytes(b'external')
            value['output_sha256']['../outside.txt']=sha(outside)
            writej(self.stage/'receipt.json',value);return value
        with patch.object(runner,'run_job',side_effect=changed):
            with self.assertRaisesRegex(ValueError,'Unsafe|reserved'):self.call()
        self.no_receipt()

    def test_required_output_must_be_in_receipt(self):
        def changed(*args,**kwargs):
            value=self.mock_run(*args,**kwargs);del value['output_sha256']['job.xyz']
            writej(self.stage/'receipt.json',value);return value
        with patch.object(runner,'run_job',side_effect=changed):
            with self.assertRaisesRegex(ValueError,'required output'):self.call()
        self.no_receipt()

    def test_original_locks_block_before_stage_creation(self):
        for name in ('.clusterflow-run.lock','.orca-run.lock'):
            with self.subTest(name=name):
                (self.job/name).write_text('{"token":"someone-else"}',encoding='utf-8')
                with patch.object(runner,'run_job',side_effect=AssertionError('already locked')):
                    with self.assertRaises(FileExistsError):self.call()
                self.assertEqual(readj(self.job/name)['token'],'someone-else')
                self.assertFalse(self.scratch.exists())
                (self.job/name).unlink()

    def test_second_scratch_cannot_race_same_original_job(self):
        def concurrent(*args,**kwargs):
            with self.assertRaises(FileExistsError):self.call(scratch_root=self.base/'other_scratch')
            self.assertFalse((self.base/'other_scratch').exists())
            return self.mock_run(*args,**kwargs)
        with patch.object(runner,'run_job',side_effect=concurrent):self.call()

    def test_existing_scratch_is_not_silently_reused(self):
        stage=Path(self.call(dry_run=True)['ascii_scratch_directory']);stage.mkdir(parents=True)
        (stage/'retained.txt').write_text('preserve',encoding='utf-8')
        with patch.object(runner,'run_job',side_effect=AssertionError('existing scratch')):
            with self.assertRaises(FileExistsError):self.call()
        self.assertEqual((stage/'retained.txt').read_text(),'preserve')
        self.no_receipt()

    def test_unreceipted_products_blocked_without_overwrite(self):
        for name in ('job.out','job.stderr','job.gbw','job.hess','job.property.txt'):
            with self.subTest(name=name):
                path=self.job/name;path.write_bytes(b'preserve')
                with patch.object(runner,'run_job',side_effect=AssertionError('existing output')):
                    with self.assertRaisesRegex(ValueError,'Unreceipted'):self.call()
                self.assertEqual(path.read_bytes(),b'preserve');path.unlink()

    def test_matching_completed_receipt_resumes_without_execution(self):
        self.successful();before=(self.job/'receipt.json').read_bytes()
        with patch.object(runner,'run_job',side_effect=AssertionError('completed job relaunched')):
            value=self.call(scratch_root=self.base/'different_scratch')
        self.assertTrue(value['skipped'])
        self.assertEqual((self.job/'receipt.json').read_bytes(),before)
        self.assertFalse((self.base/'different_scratch').exists())

    def test_tampered_completed_output_requires_review_not_execution(self):
        self.successful();(self.job/'job.out').write_bytes(b'keep tampered for review')
        with patch.object(runner,'run_job',side_effect=AssertionError('must not overwrite')):
            with self.assertRaisesRegex(ValueError,'output hash mismatch'):self.call()
        self.assertEqual((self.job/'job.out').read_bytes(),b'keep tampered for review')

    def test_failed_receipt_requires_explicit_retry(self):
        def failure(*args,**kwargs):
            value=self.mock_run(*args,**kwargs);value['status']='failed'
            writej(self.stage/'receipt.json',value);return value
        with patch.object(runner,'run_job',side_effect=failure):self.call()
        with patch.object(runner,'run_job',side_effect=AssertionError('implicit retry')):
            with self.assertRaisesRegex(ValueError,'explicit --retry'):self.call()

    def test_explicit_retry_archives_entire_old_output_inventory(self):
        self.successful();old_stage=self.stage
        old_receipt=(self.job/'receipt.json').read_bytes()
        # Receipt-bound old property file absent from next attempt must not linger.
        (self.job/'job.property.txt').write_bytes(b'OLD SYNTHETIC PROPERTY')
        old=readj(self.job/'receipt.json');old['output_sha256']['job.property.txt']=sha(self.job/'job.property.txt')
        writej(self.job/'receipt.json',old);old_receipt=(self.job/'receipt.json').read_bytes()
        self.successful(retry=True)
        self.assertNotEqual(self.stage,old_stage)
        archives=list((self.job/'attempts').glob('before_explicit_retry_*'))
        self.assertEqual(len(archives),1)
        self.assertEqual((archives[0]/'receipt.json').read_bytes(),old_receipt)
        self.assertTrue((archives[0]/'worker/evidence.txt').is_file())
        self.assertTrue((archives[0]/'job.property.txt').is_file())
        self.assertFalse((self.job/'job.property.txt').exists())

    def test_output_appearing_during_run_is_never_overwritten(self):
        def conflict(*args,**kwargs):
            value=self.mock_run(*args,**kwargs)
            (self.job/'job.out').write_bytes(b'competing output')
            return value
        with patch.object(runner,'run_job',side_effect=conflict):
            with self.assertRaisesRegex(ValueError,'overwrite existing'):self.call()
        self.assertEqual((self.job/'job.out').read_bytes(),b'competing output')
        self.no_receipt()

    def test_copy_creation_race_is_exclusive(self):
        original=runner._copy_new
        def competing(source,destination,expected):
            if destination==self.job/'job.out':destination.write_bytes(b'competing during copy')
            return original(source,destination,expected)
        with patch.object(runner,'run_job',side_effect=self.mock_run),patch.object(runner,'_copy_new',side_effect=competing):
            with self.assertRaises(FileExistsError):self.call()
        self.assertEqual((self.job/'job.out').read_bytes(),b'competing during copy')
        self.no_receipt()

    def test_corrupted_copy_rejected_before_receipt(self):
        original=runner._copy_new
        def corrupt(source,destination,expected):
            original(source,destination,expected)
            if destination==self.job/'job.out':destination.write_bytes(b'corrupt copy')
        with patch.object(runner,'run_job',side_effect=self.mock_run),patch.object(runner,'_copy_new',side_effect=corrupt):
            with self.assertRaisesRegex(ValueError,'output hash mismatch'):self.call()
        self.no_receipt()

    def test_inputs_changing_during_copy_reject_receipt(self):
        original=runner._copy_new
        def change(source,destination,expected):
            original(source,destination,expected)
            if destination==self.job/'job.out':(self.job/'input_structure.mol').write_bytes(b'changed')
        with patch.object(runner,'run_job',side_effect=self.mock_run),patch.object(runner,'_copy_new',side_effect=change):
            with self.assertRaisesRegex(ValueError,'immutable inputs changed'):self.call()
        self.no_receipt()

    def test_receipt_commit_is_last_and_atomic(self):
        original=runner._commit_receipt
        def check(source,destination,expected,validate):
            self.assertFalse(destination.exists())
            self.assertTrue((self.job/'portable_stageout.json').is_file())
            for name,digest in readj(source)['output_sha256'].items():self.assertEqual(sha(self.job/name),digest)
            return original(source,destination,expected,validate)
        with patch.object(runner,'run_job',side_effect=self.mock_run),patch.object(runner,'_commit_receipt',side_effect=check):
            self.call()
        self.assertFalse(list(self.job.glob('.receipt.*.tmp')))

    def test_input_change_at_final_receipt_commit_is_rejected(self):
        original=runner._copy_new
        def change(source,destination,expected):
            original(source,destination,expected)
            if destination.parent==self.job and destination.name.startswith('.receipt.'):
                (self.job/'input.xyz').write_bytes(b'changed at commit')
        with patch.object(runner,'run_job',side_effect=self.mock_run),patch.object(runner,'_copy_new',side_effect=change):
            with self.assertRaisesRegex(ValueError,'immutable inputs changed'):self.call()
        self.no_receipt()
        self.assertFalse(list(self.job.glob('.receipt.*.tmp')))

    def test_competing_receipt_is_never_replaced_at_commit(self):
        original=runner._commit_receipt
        def competing(source,destination,expected,validate):
            destination.write_bytes(b'competing receipt must survive')
            return original(source,destination,expected,validate)
        with patch.object(runner,'run_job',side_effect=self.mock_run),patch.object(runner,'_commit_receipt',side_effect=competing):
            with self.assertRaises(FileExistsError):self.call()
        self.assertEqual((self.job/'receipt.json').read_bytes(),b'competing receipt must survive')
        self.assertFalse(list(self.job.glob('.receipt.*.tmp')))

    def test_failure_and_success_restore_process_environment(self):
        os.environ['PATH']='before';os.environ['OMP_NUM_THREADS']='8'
        with patch.object(runner,'run_job',side_effect=RuntimeError('synthetic failure')):
            with self.assertRaises(RuntimeError):self.call()
        self.assertEqual(os.environ['PATH'],'before');self.assertEqual(os.environ['OMP_NUM_THREADS'],'8')
        self.no_receipt()
        self.successful(scratch_root=self.base/'new_scratch')
        self.assertEqual(os.environ['PATH'],'before');self.assertEqual(os.environ['OMP_NUM_THREADS'],'8')


if __name__=='__main__':unittest.main()
