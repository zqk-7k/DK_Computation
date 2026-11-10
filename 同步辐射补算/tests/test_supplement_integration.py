"""Integration and failure-boundary tests; fixtures are never DFT/FEFF results."""
from pathlib import Path
import json,shutil,sys,tempfile,unittest
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
sys.path.insert(0,str(ROOT/'tests'))
import supplement as flow
import prepare_inputs
import export_supplement
from build_candidates import build,load_cif
from qe_tasks import prepare,collect_one
from test_candidates_qe import output_fixture

class Integration(unittest.TestCase):
    def test_full_model_pipeline_contract_and_tamper_rejection(self):
        with tempfile.TemporaryDirectory(prefix='SUPPLEMENT_INTEGRATION_TEST_ONLY_') as temp:
            root=Path(temp);shutil.copytree(ROOT/'structures_original',root/'structures_original')
            shutil.copytree(ROOT/'docs',root/'docs')
            c=flow.read(ROOT/'config.json');c['dft']['command']=['TEST_ONLY_NOT_EXECUTED']
            config=root/'config.json';flow.write(config,c)
            models=build(config);jobs=prepare(config,allow_missing=True)
            self.assertEqual(len(models),12)
            # Parser fixtures only: do not run QE or label any package input optimized.
            for job in jobs:
                s,_=load_cif(root/job['candidate_path']/'structure.cif')
                (root/job['path']/'qe.out').write_text(output_fixture(s),encoding='utf-8')
                result=collect_one(root,job,{'returncode':0,'fixture_only':True})
                self.assertEqual(result['status'],'optimized')
                self.assertTrue(flow.dft_done(job,root))
            selected=flow.coordination(config)
            self.assertEqual(len(selected),12)
            fc=flow.make_feff_config(config,selected)
            # Actual FEFF input preparation against parser-fixture geometries.
            self.assertEqual(prepare_inputs.main(['--config',str(fc),'--profile','supplement']),0)
            feff=flow.read(root/'manifests/supplement_jobs.json')
            self.assertEqual(len(feff),384)
            self.assertEqual(len({j['model_id'] for j in feff}),12)
            available={mode:[{'model':j['model_id'],'edge':j['edge']} for j in feff if j['mode']==mode] for mode in ['xanes','exafs']}
            _,comparisons=export_supplement.validate_comparisons(c,{'models_and_edges':available},root)
            self.assertEqual(len(comparisons),8)
            p=root/jobs[0]['path']/'qe.out';p.write_text(p.read_text()+'tampered',encoding='utf-8')
            self.assertFalse(flow.dft_done(jobs[0],root))

    def test_array_requests_MPI_processes_and_singleton_PBS(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);c=flow.read(ROOT/'config.json');config=root/'config.json';flow.write(config,c)
            script=flow.render_dft_array(root,config,c,12,'slurm','test')
            s=script.read_text(encoding='utf-8');self.assertIn('--ntasks=16',s);self.assertIn('--array=0-11%2',s)
            self.assertIn('SLURM_ARRAY_TASK_ID',s)
            script=flow.render_dft_array(root,config,c,1,'pbs','test2')
            s=script.read_text(encoding='utf-8');self.assertIn('mpiprocs=16',s);self.assertNotIn('#PBS -J',s)

    def test_MPI_probe_rejects_unverified_serial_executable(self):
        from probe_qe import probe
        with tempfile.TemporaryDirectory() as temp:
            fake=Path(temp)/'pw.x';fake.write_text('TEST_ONLY')
            class Response:
                stdout='Program PWSCF v.7.4.1\n Serial version\n';stderr=''
            with patch('probe_qe.subprocess.run',return_value=Response()):
                with self.assertRaisesRegex(RuntimeError,'MPI'):probe(fake,16)

if __name__=='__main__':unittest.main()
