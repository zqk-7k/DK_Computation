import hashlib
from pathlib import Path
import tempfile
import unittest
from clusterflow.common import writej,readj,sha
from clusterflow.scheduling import single_execution,DeferredWork
from scripts.goat_portable import IMMUTABLE_FILES
from scripts.recover_goat_stageout_20261020 import recover


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.base=Path(self.tmp.name);self.root=self.base/'root'
        self.original=self.root/'jobs/j';self.stage=self.base/'scratch/attempt/job'
        self.original.mkdir(parents=True);self.stage.mkdir(parents=True)
        writej(self.root/'registry/jobs.json',[dict(job_id='j',species_id='s',job_directory='jobs/j')])
        writej(self.root/'config/site_scnet.json',dict(shared_scratch=str(self.base/'scratch')))
        for n in IMMUTABLE_FILES:
            (self.original/n).write_text('input '+n)
            (self.stage/n).write_text('input '+n)
        inputs={n:sha(self.original/n) for n in IMMUTABLE_FILES}
        self.intent=dict(original=str(self.original),stage=str(self.stage),job_identity='j',inputs=inputs)
        writej(self.original/'stage_intents/one.json',self.intent)
        (self.stage/'goat.out').write_text('real completed output')
        (self.stage/'goat.finalensemble.xyz').write_text('1\nenergy\nH 0 0 0\n')
        self.receipt=dict(job_id='j',status='COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW',
                          input_sha256={n:v for n,v in inputs.items() if n!='isotope_map.json'},
                          output_sha256={n:sha(self.stage/n) for n in ['goat.out','goat.finalensemble.xyz']})
        writej(self.stage/'receipt.json',self.receipt)
        writej(self.root/'operations/goat_bulk_20261015/queue.json',
               dict(cursor=1,active={'j':dict(job_id='j',species_id='s',token='t')},skipped=[]))

    def test_complete_recovery_copies_original_receipt_unchanged(self):
        result=recover(self.root,'j')
        self.assertEqual(result['quantum_jobs_launched'],0)
        self.assertEqual(sha(self.original/'receipt.json'),sha(self.stage/'receipt.json'))
        self.assertEqual(readj(self.root/'operations/goat_bulk_20261015/queue.json')['active'],{})

    def test_partial_copy_retained_before_replacement(self):
        (self.original/'goat.out').write_text('partial earlier copy')
        result=recover(self.root,'j')
        retained=Path(result['retained_conflicting_partials'][0]['retained_path'])
        self.assertEqual(retained.read_text(),'partial earlier copy')
        self.assertEqual(sha(self.original/'goat.out'),sha(self.stage/'goat.out'))

    def test_running_receipt_never_promoted(self):
        self.receipt['status']='RUNNING';writej(self.stage/'receipt.json',self.receipt)
        with self.assertRaises(ValueError):recover(self.root,'j')
        self.assertFalse((self.original/'receipt.json').exists())

    def test_changed_scratch_blocks_recovery(self):
        (self.stage/'goat.out').write_text('unbound replacement')
        with self.assertRaises(ValueError):recover(self.root,'j')
        self.assertFalse((self.original/'goat.out').exists())

    def test_active_producer_flock_blocks_recovery(self):
        key='goat_'+hashlib.sha256(b'j').hexdigest()
        with single_execution(self.root,key):
            with self.assertRaises(DeferredWork):recover(self.root,'j')
        self.assertFalse((self.original/'receipt.json').exists())

    def test_failed_receipt_remains_failed(self):
        self.receipt['status']='FAILED';writej(self.stage/'receipt.json',self.receipt)
        recover(self.root,'j')
        self.assertEqual(readj(self.original/'receipt.json')['status'],'FAILED')


if __name__=='__main__':unittest.main()
