"""Orchestration tests only. Mock spectra are temporary and never final data."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import workflow


FAKE_ENGINE=r'''# TEST ONLY: does not run FEFF or calculate a physical spectrum.
import json,math
from pathlib import Path
j=json.loads(Path('job.json').read_text())
e0={'W_L3':10207.,'Co_K':7709.,'Fe_K':7112.}[j['edge']]
if j['mode']=='xanes':
 with open('xmu.dat','w') as f:
  f.write('# TEST ONLY fixture; omega e k mu mu0 chi\n')
  for i in range(401):
   e=-20+i*.3;k=math.sqrt(max(e,0)/3.80998212)
   mu0=.5+math.atan(e/1.5)/math.pi
   extra=.5*math.exp(-((e-4)/2.5)**2)
   f.write(f'{e0+e} {e} {k} {mu0+extra} {mu0} {extra}\n')
else:
 with open('chi.dat','w') as f:
  f.write('# TEST ONLY fixture; k chi\n')
  for i in range(281):
   k=i*.05;chi=.15*math.sin(3.2*k)*math.exp(-.01*k*k)/max(k,1)
   f.write(f'{k} {chi}\n')
 Path('paths.dat').write_text('# TEST ONLY path listing\n1 2 1.9\n')
 with open('feff0001.dat','w') as f:
  f.write('# TEST ONLY path fixture\n')
  for i in range(1,100):f.write(f'{i*.1} 0 1 0 1 10 1\n')
'''


class CoordinatorTests(unittest.TestCase):
    def test_complete_requires_current_output(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);job=root/'job';job.mkdir()
            (job/'feff.inp').write_text('XANES\n')
            out=job/'xmu.dat';out.write_text('valid test payload')
            status={'status':'complete','input_sha256':workflow.sha(job/'feff.inp'),
                    'validation':{'files':[{'name':'xmu.dat','sha256':workflow.sha(out)}]}}
            workflow.write_json(job/'run_result.json',status)
            entry={'path':'job','mode':'xanes'}
            self.assertEqual(workflow.job_status(entry,root),'complete')
            out.write_text('changed')
            self.assertEqual(workflow.job_status(entry,root),'pending')
            out.unlink()
            self.assertEqual(workflow.job_status(entry,root),'pending')

    def test_pbs_single_job_and_slurm_array(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);config={'scheduler':{'max_concurrent':4}}
            script,_=workflow.render_array(root,root/'config.json',config,'pilot',[{'path':'a'}],'pbs','one')
            text=script.read_text(encoding='utf-8')
            self.assertNotIn('#PBS -J',text)
            self.assertIn('#PBS -S /bin/bash',text)
            self.assertIn('--index 0',text)
            script,_=workflow.render_array(root,root/'config.json',config,'pilot',[{'path':'a'}]*6,'slurm','six')
            self.assertIn('--array=0-5%4',script.read_text(encoding='utf-8'))

    def test_unknown_scheduler_state_does_not_prepare_or_submit(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);c=root/'config.json'
            workflow.write_json(c,{'scheduler':{'type':'slurm'}})
            workflow.write_json(root/'structure_checksums.json',[])
            workflow.write_json(root/'submissions/pilot/current_submission.json',{'job_id':'123','scheduler':'slurm'})
            with patch.object(workflow,'scheduler_active',return_value=None),patch.object(workflow,'prepare') as prep:
                with self.assertRaisesRegex(RuntimeError,'Cannot determine'):
                    workflow.run_profile(c,'pilot')
                prep.assert_not_called()

    def test_pbs_finished_state(self):
        response=subprocess.CompletedProcess([],0,'job_state = F\n','')
        with patch.object(workflow.subprocess,'run',return_value=response):
            self.assertFalse(workflow.scheduler_active('pbs','123.server'))

    def test_damaged_result_json_is_pending(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);(root/'run_result.json').write_text('{broken')
            self.assertEqual(workflow.job_status({'path':'.'},root),'pending')

    def test_end_to_end_temporary_mock_engine(self):
        # This validates coordinator hand-offs, not FEFF physics or Linux build.
        with tempfile.TemporaryDirectory(prefix='TEST_ONLY_XAS_COORDINATOR_') as t:
            root=Path(t);(root/'structures').mkdir()
            source=ROOT/'structures/Top-O1.cif'
            shutil.copy2(source,root/'structures/Top-O1.cif')
            workflow.write_json(root/'structure_checksums.json',[{'file':'structures/Top-O1.cif','sha256':workflow.sha(source)}])
            engine=root/'TEST_ONLY_engine.py';engine.write_text(FAKE_ENGINE,encoding='utf-8')
            c=copy.deepcopy(workflow.read_json(ROOT/'config.json'))
            c['software']['feff_command']=[sys.executable,str(engine)]
            c['scheduler'].update(type='local',max_concurrent=2)
            c['profiles']['core']={'models':['Top-O1.cif'],'edges':['W_L3','Co_K','Fe_K'],
                                   'site_selection':'representative_nearest_w'}
            config=root/'config.json';workflow.write_json(config,c)
            command=[sys.executable,str(ROOT/'scripts/workflow.py'),'execute','--config',str(config),
                     '--profile','core','--scheduler','local','--allow-local']
            result=subprocess.run(command,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=300)
            self.assertEqual(result.returncode,0,result.stdout[-4000:]+'\n'+result.stderr[-4000:])
            completion=workflow.read_json(root/'validation/workflow_completion.json')
            self.assertEqual(completion['status'],'computed_and_postprocessed')
            self.assertTrue(list((root/'results').rglob('*.pdf')))
            entries=workflow.read_json(root/'manifests/core_jobs.json')
            self.assertTrue(all(workflow.job_status(j,root)=='complete' for j in entries))
            # Resume keeps already completed engine output, verified by timestamps.
            before={p: p.stat().st_mtime_ns for p in (root/'jobs').rglob('chi.dat')}
            resumed=subprocess.run(command,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=300)
            self.assertEqual(resumed.returncode,0,resumed.stderr[-3000:])
            self.assertEqual(before,{p:p.stat().st_mtime_ns for p in before})


if __name__=='__main__':unittest.main(verbosity=2)
