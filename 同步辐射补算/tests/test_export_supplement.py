"""Isolated exporter checks; analytical fixtures are never scientific results."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import export_supplement as export


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def fixture(root, narrow=False):
    result = root / "TEST_ONLY_RESULTS"
    models = ["TEST_ONLY_hydroxylated", "TEST_ONLY_deprotonated", "TEST_ONLY_EG_bound"]
    entries = []
    for model_index, model in enumerate(models):
        for edge, threshold in [("W_L3", 10206), ("Co_K", 7710), ("Fe_K", 7115)]:
            directory = result / "means" / model / edge
            (directory / "xanes").mkdir(parents=True)
            (directory / "exafs").mkdir()
            x = np.arange(threshold-120, threshold+121, 1.0)
            if narrow and edge == "W_L3":
                x = np.arange(10200, 10280, 1.0)
            mu = 0.5 + 0.2*np.tanh((x-threshold)/5) + 0.03*model_index*np.exp(-((x-threshold-6)/8)**2)
            np.savetxt(directory / "xanes/mean_mu.csv", np.column_stack([x, mu, np.ones(len(x)), mu-1]),
                       delimiter=",", header="photon_energy_eV,mean_mu,mean_mu0,mean_native_chi_or_scaled_chi", comments="")
            k = np.arange(0,14.001,0.05)
            chi = np.sin(k*(3.1+0.1*model_index))*np.exp(-k/5)*.02
            chi_table=np.column_stack([k,chi,k*chi,k*k*chi,k**3*chi])
            headers="k_A_inv,chi_signed,k_chi_A_inv,k2_chi_A_inv2,k3_chi_A_inv3"
            for name,mask in [("chi_full_grid.csv",np.ones(len(k),bool)),("chi_common_3_to_12.csv",(k>=3)&(k<=12))]:
                np.savetxt(directory/"exafs"/name,chi_table[mask],delimiter=",",header=headers,comments="")
            r=np.arange(0,6.001,.04)
            value=np.exp(-((r-(1.72+.04*model_index))/.24)**2)*(1+0.2j)
            for weight in (1,2,3):
                np.savetxt(directory/"exafs"/f"FT_kweight{weight}.csv",
                           np.column_stack([r,value.real*weight,value.imag*weight,np.abs(value)*weight]),
                           delimiter=",",header="R_A_phase_uncorrected,real,imag,magnitude",comments="")
            kk,rr=np.meshgrid(k,r)
            wt=np.exp(-((rr-1.72)/.35)**2)*np.exp(-((kk-6)/2.2)**2)*(1+.15j)*(model_index+1)
            for part,array in [("REAL",wt.real),("IMAG",wt.imag),("MAG",np.abs(wt))]:
                np.savetxt(directory/"exafs"/f"WT_Cauchy_{part}.csv",np.column_stack([r,array]),delimiter=",",
                           header="R_A_phase_uncorrected_vs_k_A_inv,"+",".join(map(str,k)),comments="")
            np.savez_compressed(directory/"exafs/complex_transforms.npz",WT_k_A_inv=k,WT_R_A=r,WT_Cauchy_complex=wt)
            entries.append({"model":model,"edge":edge})
    report={"status":"complete","partial":False,"is_full_element_average":True,
            "fixture":"TEST ONLY program fixture; not a FEFF calculation",
            "parameters":{},"models_and_edges":{"xanes":entries,"exafs":entries}}
    write_json(result/"postprocess_report.json",report)
    inventory=[{"path":p.relative_to(result).as_posix(),"sha256":export.sha(p),"bytes":p.stat().st_size}
               for p in result.rglob('*') if p.is_file()]
    write_json(result/"file_inventory.json",inventory)
    config={"export":{"dpi":90,"w_energy_window_eV":[10100,10300],"comparisons":[
        {"id":"TEST_ONLY_comparison","title":"TEST ONLY exporter verification",
         "series":[{"model":models[0],"label":"Hydroxylated"},
                   {"model":models[1],"label":"Deprotonated (no EG)"},
                   {"model":models[2],"label":"EG bound at W"}]}]}}
    config_file=root/"config.json";write_json(config_file,config)
    return config_file,result,config,report


class ExportTests(unittest.TestCase):
    def test_complete_export_values_and_images(self):
        with tempfile.TemporaryDirectory(prefix="XAS_EXPORT_TEST_ONLY_") as temp:
            config,result,_,_=fixture(Path(temp))
            output=export.run(config,result)
            summary=export.read_json(output/"export_report.json")
            self.assertEqual(summary["status"],"complete")
            self.assertEqual(len(list((output/"figures").glob('*'))),16)
            for extension in export.FORMATS:
                self.assertEqual(len(list((output/"figures").glob('*.'+extension))),4)
            for table in summary["txt_tables"]:
                if table["source"].endswith('.csv'):
                    _,a=export.load_csv(result/table["source"])
                    b=np.loadtxt(output/table["output"])
                    np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-11)
            p=output/"txt/comparisons/TEST_ONLY_comparison/W_L3_XANES.txt"
            table=np.loadtxt(p)
            self.assertEqual(table[0,0],10100)
            self.assertEqual(table[-1,0],10300)
            peaks=summary["comparisons"][0]["W_first_shell_region_peaks"]
            self.assertGreater(peaks[0]["R_peak_A_phase_uncorrected"],1.6)
            self.assertNotEqual(peaks[0]["R_peak_A_phase_uncorrected"],1.5)

    def test_narrow_w_window_is_rejected_without_output(self):
        with tempfile.TemporaryDirectory(prefix="XAS_EXPORT_TEST_ONLY_") as temp:
            config,result,_,_=fixture(Path(temp),narrow=True)
            with self.assertRaisesRegex(export.DataError,"no padding or extrapolation"):
                export.run(config,result)
            self.assertEqual(list(result.glob('FINAL_EXPORT_*')),[])

    def test_changed_data_and_partial_report_are_rejected(self):
        with tempfile.TemporaryDirectory(prefix="XAS_EXPORT_TEST_ONLY_") as temp:
            config,result,_,report=fixture(Path(temp))
            path=next(result.rglob('mean_mu.csv'))
            path.write_text(path.read_text()+"\n",encoding='utf-8')
            with self.assertRaisesRegex(export.DataError,"changed or missing"):
                export.validate_source(result)
            report['status']='PARTIAL_DIAGNOSTIC';write_json(result/'postprocess_report.json',report)
            with self.assertRaisesRegex(export.DataError,"partial diagnostics"):
                export.validate_source(result)

    def test_voltage_labels_and_missing_reference_are_rejected(self):
        with tempfile.TemporaryDirectory(prefix="XAS_EXPORT_TEST_ONLY_") as temp:
            path,result,config,report=fixture(Path(temp))
            series=config['export']['comparisons'][0]['series']
            series[0]['label']='OCP'
            with self.assertRaisesRegex(export.DataError,"structural state labels"):
                export.validate_comparisons(config,report,path.parent)
            series[0]['label']='Hydroxylated'
            config['export']['energy_reference_shift_eV']={'W_L3':4.0}
            with self.assertRaisesRegex(export.DataError,"require energy_reference_evidence"):
                export.validate_comparisons(config,report,path.parent)

    def test_missing_state_is_not_replaced(self):
        with tempfile.TemporaryDirectory(prefix="XAS_EXPORT_TEST_ONLY_") as temp:
            path,result,config,report=fixture(Path(temp))
            config['export']['comparisons'][0]['series'][2]['model']='missing_model'
            with self.assertRaisesRegex(export.DataError,"Missing xanes"):
                export.validate_comparisons(config,report,path.parent)


if __name__=='__main__':
    unittest.main()
