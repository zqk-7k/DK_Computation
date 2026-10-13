"""Audit and archive this delivery. No production or smoke QM is launched."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from clusterflow.common import config, hash_protocol

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

def main(args):
    p=config(ROOT)
    if ROOT.name!=p['package_name']:raise ValueError('Package name mismatch')
    archive=ROOT.parent/(ROOT.name+'.zip')
    receipt=ROOT.parent/(ROOT.name+'-压缩包校验.json')
    if archive.exists() or receipt.exists():raise ValueError('Delivery archive already exists; use a new numbered version')
    original=Path(args.source_v2).resolve()
    old_manifest=read(original/'reports/delivery_file_manifest.json')
    changed=[r['path'] for r in old_manifest if not (original/r['path']).is_file() or sha(original/r['path'])!=r['sha256']]
    if changed:raise ValueError('Original v2 changed: '+str(changed[:10]))
    # Only immutable data/model/reference inputs are expected byte-identical in v3.
    preserved=[r for r in old_manifest if r['path'].startswith(('01_goat_explicit/','solute_models/','reference_51/','registry/')) or r['path'].endswith('.sdf')]
    bad=[r['path'] for r in preserved if not (ROOT/r['path']).is_file() or sha(ROOT/r['path'])!=r['sha256']]
    if bad:raise ValueError('Initial inputs changed: '+str(bad[:10]))
    smoke=Path(args.smoke_source).resolve()
    records=read(smoke/'summary.json')
    if len(records)!=2 or not all(r['passed'] for r in records):raise ValueError('Smoke evidence incomplete')
    dest=ROOT/'validation/real_smoke'
    for r in records:
        s=smoke/r['stage'];d=dest/r['stage'];d.mkdir(parents=True,exist_ok=True)
        if sha(s/'job.inp')!=r['input_sha256'] or sha(s/'job.out')!=r['output_sha256']:raise ValueError('Smoke input/output changed')
        for name in ['job.inp','job.out','job.stderr','job.property.txt','job.bibtex','atom_map.json','smoke_report.json']:
            target=d/name
            if target.exists() and sha(target)!=sha(s/name):raise ValueError('Smoke destination differs')
            if not target.exists():shutil.copy2(s/name,target)
    write(dest/'summary.json',records)
    write(dest/'scope.json',{'scope':'TWO_SIX_ATOM_WATER_DIMER_ENGINEERING_TESTS_ONLY',
          'target_product_calculations':0,'GOAT_calculations':0,'scientific_accuracy_validation':False,
          'not_validated':['DMSO cluster stability','sampling convergence','target thermochemistry','experimental NMR accuracy'],
          'large_restart_files_omitted_from_delivery':True})
    audit_path=ROOT/'reports/initial_938_audit_v3.json'
    if args.reuse_initial_audit:
        report=read(audit_path)
        if report['auditor_script_sha256']!=sha(ROOT/'scripts/independent_audit.py') or report['inputs']['jobs_manifest_sha256']!=sha(ROOT/'registry/jobs.json'):
            raise ValueError('Cached static audit inputs changed')
    else:
        spec=importlib.util.spec_from_file_location('initial_audit',ROOT/'scripts/independent_audit.py')
        audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
        report=audit.analyse(ROOT)
        write(audit_path,report)
    if report['error_count']:raise ValueError('Initial 938 audit failed')
    if args.prepare_only:
        print(json.dumps({'static_audit':str(audit_path),'jobs_passing':report['scope']['jobs_passing'],'errors':report['error_count']},ensure_ascii=False))
        return
    review=read(ROOT/'reports/workflow_independent_review.json')
    stale=[name for name,digest in review['code_snapshot_sha256'].items() if sha(ROOT/name)!=digest]
    if stale:raise ValueError('Independent code snapshot stale: '+str(stale))
    run=subprocess.run([sys.executable,'-B','-X','utf8','-m','unittest','discover','-s','tests','-v'],cwd=ROOT,capture_output=True,text=True,encoding='utf-8')
    (ROOT/'reports/unit_tests_v3.txt').write_text(run.stdout+run.stderr,encoding='utf-8')
    if run.returncode:raise ValueError('Unit tests failed: reports/unit_tests_v3.txt')
    status=subprocess.run([sys.executable,'-B','-X','utf8','workflow.py','status'],cwd=ROOT,capture_output=True,text=True,encoding='utf-8')
    if status.returncode:raise ValueError('Status command failed')
    write(ROOT/'reports/status_v3.json',json.loads(status.stdout))
    info={'package_name':ROOT.name,'generated_utc':datetime.now(timezone.utc).isoformat(),
          'protocol_hash':hash_protocol(p),'reference_commit':p['reference_commit'],
          'original_v2_manifest_entries_unchanged':len(old_manifest),'v3_immutable_entries_preserved':len(preserved),
          'initial_jobs_audited':report['scope']['jobs_fully_audited'],'initial_errors':report['error_count'],
          'unit_tests_exit_code':run.returncode,'small_real_engineering_tests':2,
          'production_QM_jobs_launched':0,'target_product_spectra_calculated':False,
          'production_scientific_validation':'NOT_COMPLETED',
          'missing_local_GOAT_dependency':'otool_xtb and associated runtime; user must configure compatible licensed/available software',
          'revDSD_active':False,'independent_scenarios':p['independent_scenarios'],
          'scenario_prerequisites':p['scenario_prerequisites']}
    write(ROOT/'reports/delivery_v3.json',info)
    files=sorted(x for x in ROOT.rglob('*') if x.is_file() and '__pycache__' not in x.parts and x.name!='delivery_file_manifest_v3.json')
    manifest=[{'path':x.relative_to(ROOT).as_posix(),'bytes':x.stat().st_size,'sha256':sha(x)} for x in files]
    write(ROOT/'reports/delivery_file_manifest_v3.json',manifest)
    files.append(ROOT/'reports/delivery_file_manifest_v3.json')
    with zipfile.ZipFile(archive,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for x in files:z.write(x,ROOT.name+'/'+x.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:raise ValueError('ZIP CRC error')
        for item in manifest:
            if hashlib.sha256(z.read(ROOT.name+'/'+item['path'])).hexdigest()!=item['sha256']:raise ValueError('ZIP content hash mismatch')
    result={'package_name':ROOT.name,'archive':str(archive),'archive_sha256':sha(archive),
            'archive_bytes':archive.stat().st_size,'archived_files':len(files),'CRC_and_contents_verified':True,
            'manifest_sha256':sha(ROOT/'reports/delivery_file_manifest_v3.json'),'scope':info}
    write(receipt,result)
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--source-v2',required=True)
    a.add_argument('--smoke-source',required=True)
    a.add_argument('--prepare-only',action='store_true')
    a.add_argument('--reuse-initial-audit',action='store_true')
    main(a.parse_args())
