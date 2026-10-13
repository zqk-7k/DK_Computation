"""Opt-in tiny ORCA parser/property smoke; NOT one of the 64 product calculations."""
from pathlib import Path
import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import json
import subprocess
import sys
import time
from rdkit import Chem

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from clusterflow.common import config,writej,sha
from clusterflow.inputs import render_input
from clusterflow.orcaio import parse_output

def model():
    # H2(18O) analyte + D2O environment. This is a six-atom engineering test,
    # not a DMSO-cluster validation or experimental structure.
    m=Chem.AddHs(Chem.MolFromSmiles('[18OH2]'))
    water=Chem.AddHs(Chem.MolFromSmiles('O'))
    for a in water.GetAtoms():
        if a.GetSymbol()=='H':a.SetIsotope(2)
    m=Chem.CombineMols(m,water)
    xyz=[(0,0,0),(.96,0,0),(-.24,.93,0),(2.86,0,0),(3.25,.86,0),(3.25,-.43,.745)]
    c=Chem.Conformer(6);c.Set3D(True)
    for i,p in enumerate(xyz):c.SetAtomPosition(i,p)
    m.AddConformer(c)
    rows=[]
    for a in m.GetAtoms():
        i=a.GetIdx();iso=a.GetIsotope()
        rows.append({'index_zero_based':i,'index_one_based':i+1,'element':a.GetSymbol(),'isotope':iso,
                     'isotope_mass_u':17.99915961286 if iso==18 else 2.01410177812 if iso==2 else None,
                     'component':'analyte' if i<3 else 'D2O_smoke_environment',
                     'environment_atom':i>=3,'source_atom_index_zero_based':i if i<3 else None,
                     'include_in_analyte_nmr':i in [1,2],'exchangeable_H':a.GetSymbol()=='H',
                     'formal_charge':0,'neighbours_zero_based':[n.GetIdx() for n in a.GetNeighbors()]})
    return m,rows

def main(args):
    if not args.execute:raise ValueError('--execute required; no automatic QM')
    base=Path(args.workdir).resolve()
    if not str(base).isascii():raise ValueError('Use an ASCII scratch path for ORCA smoke')
    if base.exists():raise ValueError('Smoke destination exists; refusing overwrite')
    base.mkdir(parents=True)
    m,amap=model();p=config(ROOT)
    def one(stage):
        folder=base/stage;folder.mkdir()
        local=copy.deepcopy(p);local[stage]['nprocs']=1;local[stage]['maxcore_mb']=1000
        text,pairs=render_input(m,stage,local,'DMSO',0,1,amap)
        (folder/'job.inp').write_text(text,encoding='utf-8')
        writej(folder/'atom_map.json',amap)
        start=time.monotonic()
        with (folder/'job.out').open('wb') as out,(folder/'job.stderr').open('wb') as err:
            done=subprocess.run([str(Path(args.orca).resolve()),'job.inp'],cwd=folder,stdout=out,stderr=err,timeout=240)
        prop=folder/'job.property.txt'
        parsed=parse_output((folder/'job.out').read_text(encoding='utf-8',errors='replace'),prop.read_text(encoding='utf-8',errors='replace') if prop.exists() else None)
        actual=[r['atom_index_zero_based'] for r in parsed['shielding']] if stage=='shielding' else [sorted([r['atom_i_zero_based'],r['atom_j_zero_based']]) for r in parsed['j_pairs']]
        expected=[1,2] if stage=='shielding' else pairs
        report={'stage':stage,'scope':'SIX_ATOM_WATER_DIMER_ENGINEERING_TEST_NOT_TARGET_PRODUCT',
                'returncode':done.returncode,'normal_termination':parsed['normal_termination'],
                'expected_analyte_nuclei_or_pairs':expected,'actual':actual,
                'passed':done.returncode==0 and parsed['normal_termination'] and sorted(actual)==sorted(expected),
                'elapsed_seconds':time.monotonic()-start,'input_sha256':sha(folder/'job.inp'),
                'output_sha256':sha(folder/'job.out'),'18O_and_D_inline_mass_annotations':True,
                'does_not_validate_GOAT_DMSO_thermochemistry_or_accuracy':True}
        writej(folder/'smoke_report.json',report)
        return report
    with ThreadPoolExecutor(max_workers=2) as pool:reports=list(pool.map(one,['shielding','j_PBE_pcJ3']))
    writej(base/'summary.json',reports)
    print(json.dumps(reports,indent=2))
    if not all(r['passed'] for r in reports):raise SystemExit(2)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--execute',action='store_true');p.add_argument('--orca',required=True);p.add_argument('--workdir',required=True)
    main(p.parse_args())
