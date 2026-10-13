"""No-QM tests, static audit, and non-overwriting, hash-verified delivery."""
from pathlib import Path
import argparse, hashlib, json, sys, unittest, zipfile
from datetime import datetime, timezone
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'scripts'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def code_hashes():
    return {p.relative_to(ROOT).as_posix():sha(p) for folder in ['clusterflow','scripts','tests','config','registry'] for p in sorted((ROOT/folder).rglob('*')) if p.is_file() and p.suffix in ('.py','.json')}
def tests():
    import numpy.testing
    suite=unittest.defaultTestLoader.discover(str(ROOT/'tests'))
    with (ROOT/'reports/engineering_tests_v5.log').open('w',encoding='utf-8') as stream:
        result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    report={'status':'PASS' if result.wasSuccessful() else 'FAIL','tests_run':result.testsRun,
        'failures':[(str(t),e) for t,e in result.failures], 'errors':[(str(t),e) for t,e in result.errors],
        'skipped':[(str(t),s) for t,s in result.skipped], 'source_sha256':code_hashes(),
        'notice':'Synthetic engineering tests, not actual quantum calculations or experimental validation.',
        'generated_utc':datetime.now(timezone.utc).isoformat()}
    save(ROOT/'reports/engineering_tests_v5.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},ensure_ascii=False,indent=2))
    if not result.wasSuccessful():raise SystemExit(1)
def audit():
    from independent_audit import analyse
    report=analyse(ROOT)
    save(ROOT/'reports/initial_cluster_audit_v5.json',report)
    print(json.dumps({k:report[k] for k in ('status','scope','error_count','errors')},ensure_ascii=False,indent=2))
    if report['error_count']:raise SystemExit(1)
def package():
    checks=[('engineering_tests_v5.json','PASS'),('initial_cluster_audit_v5.json','STATIC_INPUT_CHECKS_PASS_WITH_LIMITATIONS')]
    for name,status in checks:
        report=json.loads((ROOT/'reports'/name).read_text(encoding='utf-8'))
        if report.get('status')!=status:raise ValueError('Checks missing/failed: '+name)
        for relative,digest in report.get('source_sha256',{}).items():
            if sha(ROOT/relative)!=digest:raise ValueError('Source changed after testing: '+relative)
    target=ROOT.parent/(ROOT.name+'.zip')
    if target.exists():raise FileExistsError('Will not replace existing archive')
    manifest=ROOT/'reports/delivery_file_manifest_v5.json'
    files=sorted(p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p!=manifest and p.suffix!='.pyc' and not p.name.endswith('.inspect.ndjson'))
    records=[{'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files]
    save(manifest,{'package':ROOT.name,'status':'NO_ACTUAL_QM_RESULTS','manifest_excludes_itself':True,'files':records,'file_count':len(records)})
    with zipfile.ZipFile(target,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in files+[manifest]:z.write(p,(Path(ROOT.name)/p.relative_to(ROOT)).as_posix())
    with zipfile.ZipFile(target) as z:
        if z.testzip():raise ValueError('ZIP CRC failure')
        for record in records:
            if hashlib.sha256(z.read(ROOT.name+'/'+record['path'])).hexdigest()!=record['sha256']:
                raise ValueError('ZIP file hash mismatch')
    print(json.dumps({'zip':str(target),'sha256':sha(target),'bytes':target.stat().st_size,'files':len(records)+1,'CRC_and_hashes':'PASS'},ensure_ascii=False,indent=2))
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tests',action='store_true');parser.add_argument('--audit',action='store_true');parser.add_argument('--package',action='store_true')
    args=parser.parse_args()
    if args.tests:tests()
    if args.audit:audit()
    if args.package:package()
