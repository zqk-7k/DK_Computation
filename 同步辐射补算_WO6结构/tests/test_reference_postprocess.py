"""Software fixtures live only in temporary directories; never physical results."""
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest
import numpy as np

PACKAGE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("postprocess_reference", PACKAGE/"scripts/postprocess_reference.py")
pp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pp)


class PostprocessContractTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="XAS_SOFTWARE_FIXTURE_")
        self.root = Path(self.temp.name)
        self.old_root = pp.ROOT
        pp.ROOT = self.root
        shutil.copytree(PACKAGE/"references", self.root/"references")
        shutil.copytree(PACKAGE/"jobs_reference", self.root/"jobs_reference")
        self.xanes = self.root/"jobs_reference/full/xanes/xmu.dat"
        self.exafs = self.root/"jobs_reference/full/exafs/chi.dat"
        # Deliberately analytic fixture, not a model prediction or a measured spectrum.
        energy = np.arange(10080., 10320.1, .5)
        mu = 1/(1+np.exp(-(energy-10207)/3))
        xm = np.column_stack([energy, energy-10207, np.sqrt(np.maximum(energy-10207,0)/3.81), mu, mu*.99, mu*.01])
        k = np.arange(281)*.05
        ch = np.column_stack([k, np.sin(4*k)*np.exp(-k/5)])
        for path, data in [(self.xanes, xm),(self.exafs,ch)]:
            np.savetxt(path, data, header="TEST FIXTURE ONLY - not a FEFF calculation")
            job = json.loads((path.parent/"job.json").read_text())
            log = path.parent/"log1.dat"
            log.write_text("TEST FIXTURE ONLY\nConvergence reached in 8 iterations.\n")
            result = {"status":"complete", "mode":job["mode"],
                "source_sha256":job["source_sha256"], "input_sha256":job["feff_input_sha256"],
                "validation":{"scf_convergence":{"status":"confirmed","evidence":[{"file":"log1.dat","line":2}]},
                    "files":[{"name":path.name,"sha256":pp.digest(path),"bytes":path.stat().st_size},
                             {"name":log.name,"sha256":pp.digest(log),"bytes":log.stat().st_size}]}}
            (path.parent/"run_result.json").write_text(json.dumps(result))

    def tearDown(self):
        pp.ROOT = self.old_root
        self.temp.cleanup()

    def test_successful_source_gate(self):
        result = pp.validate_pair(self.xanes, self.exafs)
        self.assertEqual(result[0]["job"]["atom_count"],53)

    def test_changed_output_rejected(self):
        with self.exafs.open("a") as f:
            f.write("# changed\n")
        with self.assertRaisesRegex(ValueError,"output is untracked or changed"):
            pp.validate_pair(self.xanes,self.exafs)

    def test_mixed_profile_rejected(self):
        jobpath = self.exafs.parent/"job.json"
        job = json.loads(jobpath.read_text()); job["profile"]="small"
        jobpath.write_text(json.dumps(job))
        with self.assertRaisesRegex(ValueError,"different profiles"):
            pp.validate_pair(self.xanes,self.exafs)

    def test_unfinished_rejected(self):
        path = self.exafs.parent/"run_result.json"
        data=json.loads(path.read_text()); data["status"]="failed"
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError,"successful"):
            pp.validate_pair(self.xanes,self.exafs)

    def test_scf_unverified_rejected(self):
        path = self.exafs.parent/"run_result.json"
        data=json.loads(path.read_text()); data["validation"]["scf_convergence"]["status"]="unverified"
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError,"self-consistency"):
            pp.validate_pair(self.xanes,self.exafs)

    def test_missing_energy_window_rejected(self):
        data=pp.numeric_table(self.xanes,6)[100:]
        np.savetxt(self.xanes,data,header="TEST FIXTURE ONLY")
        path=self.xanes.parent/"run_result.json"
        result=json.loads(path.read_text()); result["validation"]["files"][0]["sha256"]=pp.digest(self.xanes)
        path.write_text(json.dumps(result))
        with self.assertRaisesRegex(ValueError,"does not cover"):
            pp.process(self.xanes,self.exafs,self.root/"results")

    def test_larch_and_export_contract(self):
        out=pp.process(self.xanes,self.exafs,self.root/"test_export")
        ft=np.loadtxt(out/"txt/W_L3_FT_k3.txt")
        np.testing.assert_allclose(ft[:,3],np.sqrt(ft[:,1]**2+ft[:,2]**2),rtol=1e-7,atol=1e-10)
        wt=np.loadtxt(out/"txt/W_L3_WT_k2_long.txt")
        self.assertTrue(np.isfinite(wt).all())
        self.assertTrue((out/"figures/W_L3_reference.pdf").is_file())
        self.assertTrue((out/"figures/W_L3_reference.svg").is_file())
        self.assertEqual(json.loads((out/"provenance.json").read_text())["energy_shifts"],[])


if __name__=="__main__":
    unittest.main(verbosity=2)
