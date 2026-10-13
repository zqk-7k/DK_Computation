"""Package generated inputs and read-only audit receipts; no QM launched."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile
from rdkit import Chem

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import micro64 as mc
import goat_job

def main(args):
    jobs=mc.readj(ROOT/'registry/jobs.json')
    source=mc.readj(ROOT/'registry/source_provenance.json')
    checks=[]
    for record in source['reference_files']:
        blob=subprocess.run(['git','-C',str(args.reference_clone),'show',f"{mc.SHA}:51计算包/{record['path']}"],capture_output=True,check=True).stdout
        remote_sha=hashlib.sha256(blob).hexdigest()
        if remote_sha!=record['sha256'] or mc.digest(ROOT/'reference_51'/record['path'])!=remote_sha:
            raise RuntimeError('Reference blob differs from pinned Git commit: '+record['path'])
        checks.append({'path':record['path'],'sha256':remote_sha,'pinned_git_blob_match':True})
    # Generated exchange formats; retain complete metadata in individual jobs.
    for filename,selected in [('structures_all_scenarios_start01.sdf',[j for j in jobs if j['start_index']==1]),
                              ('structures_all_938_starts.sdf',jobs)]:
        parts=[]
        for j in selected:
            mol=Chem.MolFromMolFile(str(ROOT/j['job_directory']/'cluster.mol'),removeHs=False)
            mol.SetProp('_Name',j['job_id'])
            block=Chem.MolToMolBlock(mol)
            parts.append(block+'\n>  <STATUS>\nINITIAL_COORDINATES_NOT_QM_OPTIMIZED\n\n'+
                         '>  <JOB_DIRECTORY>\n'+j['job_directory']+'\n\n'+
                         '>  <COMPOSITION>\n'+f"{j['DMSO_d6_count']} DMSO-d6 + {j['H2O_count']} H2O"+'\n\n$$$$\n')
        mc.put(ROOT/filename,''.join(parts))
    run=subprocess.run([sys.executable,'-B','-X','utf8','-m','unittest','discover','-s','tests','-v'],cwd=ROOT,capture_output=True,text=True,encoding='utf-8')
    # Receipts are new artifacts, not scientific results.
    mc.put(ROOT/'reports/unit_tests.txt',run.stdout+run.stderr)
    if run.returncode:
        raise RuntimeError('Engineering tests failed')
    preflight=goat_job.preflight(Path(args.orca),Path(args.xtb))
    mc.putj(ROOT/'reports/local_preflight.json',preflight)
    # Confirm the original v1 files covered by its delivered manifest unchanged.
    v1=Path(args.source)
    manifest=mc.readj(v1/'reports/file_manifest.json')
    if isinstance(manifest,dict):
        entries=manifest.get('files',manifest)
    else:
        entries=manifest
    preserved=[]
    if isinstance(entries,list):
        for e in entries:
            name=e.get('path',e.get('relative_path'))
            expected=e.get('sha256')
            if name and expected:
                actual=mc.digest(v1/name)
                if actual.lower()!=expected.lower():
                    raise RuntimeError('Source v1 delivered file changed: '+name)
                preserved.append(name)
    if not preserved:
        raise RuntimeError('Cannot interpret original v1 manifest; inspect schema')
    independent=mc.readj(ROOT/'reports/independent_audit.json')
    if independent.get('error_count') != 0 or independent.get('scope',{}).get('jobs_passing') != len(jobs):
        raise RuntimeError('Independent all-input audit failed or has incomplete coverage')
    delivery={'package':'ORCA_64_18O_EXPLICIT_51REF_v2','scope':mc.readj(ROOT/'reports/preparation.json'),
              'reference_files_verified_against_pinned_git':checks,'original_v1_hashes_verified_unchanged':len(preserved),
              'unit_test_exit_code':run.returncode,'independent_audit_report':'reports/independent_audit.json',
              'independent_audit_jobs_passed':independent['scope']['jobs_passing'],
              'independent_audit_error_count':independent['error_count'],
              'local_static_preflight_status':preflight['status'],
              'GOAT_DFT_NMR_production_executed':False,'fully_automated_downstream_pipeline':False,
              'notice':'INITIAL_MICROSOLVATION_MODELS_NOT_EXPERIMENTAL_DATA_NOT_QM_OPTIMIZED'}
    mc.putj(ROOT/'reports/delivery.json',delivery)
    filelist=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts or p.name=='delivery_file_manifest.json':continue
        filelist.append({'path':p.relative_to(ROOT).as_posix(),'sha256':mc.digest(p),'bytes':p.stat().st_size})
    mc.putj(ROOT/'reports/delivery_file_manifest.json',filelist)
    archive=ROOT.parent/(ROOT.name+'.zip')
    if archive.exists():raise RuntimeError('Archive exists; refusing to overwrite')
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for e in filelist+ [{'path':'reports/delivery_file_manifest.json'}]:
            z.write(ROOT/e['path'],arcname=ROOT.name+'/'+e['path'])
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:raise RuntimeError('ZIP CRC check failed')
    mc.putj(ROOT.parent/(ROOT.name+'_archive_receipt.json'),{'archive':str(archive),'sha256':mc.digest(archive),
            'bytes':archive.stat().st_size,'CRC_test':'PASS','included_files':len(filelist)+1})
    print(json.dumps(delivery,ensure_ascii=False,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-clone',required=True)
    p.add_argument('--source',required=True)
    p.add_argument('--orca',required=True)
    p.add_argument('--xtb',required=True)
    main(p.parse_args())
