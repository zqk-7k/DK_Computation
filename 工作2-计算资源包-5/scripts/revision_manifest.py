"""Build a revision inventory from git-visible package files, excluding itself.

Generated reports only, no model or protocol rewriting. Run after tests and
before commit. The original v5 delivery manifest remains historical evidence.
"""
import json
from pathlib import Path
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from clusterflow.common import sha,writej

def main():
    repository=ROOT.parent
    destination=ROOT/'reports/delivery_file_manifest_scnet_20261014.json'
    result=subprocess.run(['git','ls-files','-z','--cached','--others','--exclude-standard','--',ROOT.name],
                          cwd=repository,capture_output=True,check=True)
    paths=sorted({repository/name.decode('utf-8') for name in result.stdout.split(b'\0') if name})
    files=[]
    for path in paths:
        if path==destination:continue
        if not path.is_file() or path.is_symlink():raise ValueError('Unexpected nonregular delivery file: '+str(path))
        files.append({'path':str(path.relative_to(ROOT)),'bytes':path.stat().st_size,'sha256':sha(path)})
    writej(destination,{'revision':'SCNET_ORCA610_20261014','status':'SOURCE_INVENTORY_NOT_QUANTUM_RESULT',
                        'self_excluded':str(destination.relative_to(ROOT)),'files':files})
    print(json.dumps({'manifest':str(destination),'files':len(files)},ensure_ascii=False))

if __name__=='__main__':main()
