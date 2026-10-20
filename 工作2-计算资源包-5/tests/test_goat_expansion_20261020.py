import tempfile
import unittest
from pathlib import Path
from clusterflow.common import writej, sha, readj
from clusterflow.goat_delegation import resolve, FILES
from scripts.expand_goat_20261020 import resource_input
from tests.test_goat_bulk import BulkTests
from scripts.goat_bulk_worker import verify_scope


class ResourceTests(unittest.TestCase):
    def test_only_nprocs_changes(self):
        text='! XTB GOAT\n%pal nprocs 8 end\n%maxcore 2000\n* xyzfile 0 1 start.xyz\n'
        self.assertEqual(resource_input(text),text.replace('nprocs 8','nprocs 16'))

    def test_wrong_or_duplicate_process_fields_rejected(self):
        for text in ('%pal nprocs 4 end','nprocs 8\nnprocs 8','! XTB GOAT'):
            with self.assertRaises(ValueError):resource_input(text)


class CohortScopeTests(BulkTests):
    def test_sixteen_cpu_scope_requires_real_matching_inputs(self):
        self.policy.update(cpus_per_job=16,slots_per_node=2)
        for job in self.policy['job_ids']:
            (self.root/self.indexed[job]['job_directory']/'goat.inp').write_text('%pal nprocs 16 end')
        path=self.root/'scope.json';writej(path,self.policy)
        verify_scope(self.root,path,sha(path))
        (self.root/self.indexed[self.policy['job_ids'][0]]['job_directory']/'goat.inp').write_text('%pal nprocs 8 end')
        with self.assertRaises(ValueError):verify_scope(self.root,path,sha(path))


class DelegationTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.cohort=self.root/'operations/cohort'
        self.original=self.root/'jobs/j1';self.target=self.cohort/'jobs/j1'
        self.original.mkdir(parents=True);self.target.mkdir(parents=True)
        for name in FILES:
            text='%pal nprocs 8 end\n! XTB GOAT' if name=='goat.inp' else 'immutable '+name
            (self.original/name).write_text(text)
            (self.target/name).write_text(resource_input(text) if name=='goat.inp' else text)
        self.row=dict(job_id='j1',job_directory='jobs/j1',species_id='s1',initial_release=True,
                      input_sha256=sha(self.original/'goat.inp'))
        writej(self.root/'registry/jobs.json',[self.row])
        writej(self.cohort/'registry/jobs.json',[dict(self.row,input_sha256=sha(self.target/'goat.inp'))])
        self.record=dict(kind='RESOURCE_COHORT_DELEGATION',token='t',job_id='j1',cpus=16,
                         cohort=str(self.cohort),canonical_manifest_sha256=sha(self.root/'registry/jobs.json'),
                         cohort_manifest_sha256=sha(self.cohort/'registry/jobs.json'),
                         canonical_inputs={n:sha(self.original/n) for n in FILES})
        self.save()

    def save(self):
        writej(self.original/'delegated_goat.json',self.record)
        writej(self.original/'.goat-run.lock',self.record)

    def test_resource_only_delegation_valid(self):
        self.assertEqual(resolve(self.root,self.row)[0],self.cohort)

    def test_geometry_change_rejected(self):
        (self.target/'start.xyz').write_text('different geometry')
        with self.assertRaises(ValueError):resolve(self.root,self.row)

    def test_chemical_input_change_rejected(self):
        (self.target/'goat.inp').write_text('%pal nprocs 16 end\n! HF')
        with self.assertRaises(ValueError):resolve(self.root,self.row)

    def test_other_active_owner_rejected(self):
        writej(self.original/'.goat-run.lock',{'token':'other'})
        with self.assertRaises(ValueError):resolve(self.root,self.row)

    def test_chained_delegation_rejected(self):
        writej(self.target/'delegated_goat.json',{})
        with self.assertRaises(ValueError):resolve(self.root,self.row)

    def test_manifest_change_rejected(self):
        writej(self.cohort/'registry/jobs.json',[])
        with self.assertRaises(ValueError):resolve(self.root,self.row)

    def test_external_cohort_rejected(self):
        self.record['cohort']='/tmp/unrelated';self.save()
        with self.assertRaises(ValueError):resolve(self.root,self.row)


if __name__=='__main__':unittest.main()
