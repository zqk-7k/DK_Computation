"""Run no-QM engineering tests and seal a non-overwriting v4 delivery archive."""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import json
import sys
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def save(path, value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

def tests():
    # Import numpy platform helpers before test fixtures prohibit child processes.
    import numpy.testing
    suite=unittest.defaultTestLoader.discover(str(ROOT/'tests'))
    with (ROOT/'reports/engineering_tests_v4.log').open('w',encoding='utf-8') as stream:
        result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    report={'status':'PASS' if result.wasSuccessful() else 'FAIL','tests_run':result.testsRun,
            'failures':[(str(t),e) for t,e in result.failures], 'errors':[(str(t),e) for t,e in result.errors],
            'skipped':[(str(t),s) for t,s in result.skipped],
            'notice':'Engineering and synthetic-fixture tests only; not actual QM or experimental validation.',
            'generated_utc':datetime.now(timezone.utc).isoformat()}
    save(ROOT/'reports/engineering_tests_v4.json',report)
    print(json.dumps(report,ensure_ascii=False,indent=2))
    if not result.wasSuccessful():raise SystemExit(1)

def package():
    for name,status in [('engineering_tests_v4.json','PASS'),('initial_cluster_audit_v4.json','STATIC_INPUT_CHECKS_PASS_WITH_LIMITATIONS')]:
        report=json.loads((ROOT/'reports'/name).read_text(encoding='utf-8'))
        if report.get('status')!=status:raise ValueError('Required checks missing or failed: '+name)
    target=ROOT.parent/(ROOT.name+'.zip')
    if target.exists():raise FileExistsError('Will not overwrite existing delivery archive: '+str(target))
    manifest_path=ROOT/'reports/delivery_file_manifest_v4.json'
    files=sorted(p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts
                 and p!=manifest_path and p.suffix!='.pyc' and not p.name.endswith('.inspect.ndjson'))
    records=[{'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files]
    save(manifest_path,{'package':ROOT.name,'status':'CODE_AND_INITIAL_STRUCTURES_NO_QM',
        'manifest_excludes_itself':True,'files':records,'file_count':len(records)})
    with zipfile.ZipFile(target,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in files+[manifest_path]:z.write(p,(Path(ROOT.name)/p.relative_to(ROOT)).as_posix())
    with zipfile.ZipFile(target) as z:
        bad=z.testzip()
        if bad:raise ValueError('Archive CRC failed: '+bad)
        for record in records:
            if hashlib.sha256(z.read(ROOT.name+'/'+record['path'])).hexdigest()!=record['sha256']:
                raise ValueError('Archive content hash mismatch: '+record['path'])
    print(json.dumps({'zip':str(target),'bytes':target.stat().st_size,'sha256':sha(target),'files':len(records)+1,'CRC_and_all_file_hashes':'PASS'},ensure_ascii=False,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tests',action='store_true');ap.add_argument('--package',action='store_true')
    args=ap.parse_args()
    if args.tests:tests()
    if args.package:package()
