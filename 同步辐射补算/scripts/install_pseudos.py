#!/usr/bin/env python3
"""Download fixed official SSSP release and verify each required UPF."""
from pathlib import Path
import argparse,hashlib,json,tarfile,urllib.request,time

BASE='https://archive.materialscloud.org/records/rcyfm-68h65/files/'
NAME='SSSP_1.3.0_PBE_efficiency.tar.gz'
DIGEST='a58f1b3373f330179fd0832c48bb9a52'
ELEMENTS=('H','C','O','Co','Fe','W')

def digest(p,kind='md5'):
    h=hashlib.new(kind)
    with Path(p).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1048576),b''):h.update(chunk)
    return h.hexdigest()

def install(root):
    meta=json.loads((root/'docs/SSSP_1.3.0_PBE_efficiency.json').read_text(encoding='utf-8'))
    target=root/'software/pseudo';target.mkdir(parents=True,exist_ok=True)
    if not all((target/meta[e]['filename']).is_file() and digest(target/meta[e]['filename'])==meta[e]['md5'] for e in ELEMENTS):
        archive=root/'software/downloads'/NAME;archive.parent.mkdir(parents=True,exist_ok=True)
        if not archive.exists() or digest(archive)!=DIGEST:
            temporary=archive.with_suffix('.download')
            url=BASE+NAME+'?download=1&cache='+str(int(time.time()))
            with urllib.request.urlopen(url,timeout=120) as response,temporary.open('wb') as out:
                while chunk:=response.read(1048576):out.write(chunk)
            if digest(temporary)!=DIGEST:raise RuntimeError('SSSP archive digest mismatch; file retained for diagnosis')
            temporary.replace(archive)
        with tarfile.open(archive,'r:gz') as tar:
            members=tar.getmembers()
            for element in ELEMENTS:
                item=meta[element];found=[m for m in members if Path(m.name).name==item['filename'] and m.isfile()]
                if len(found)!=1:raise RuntimeError('Ambiguous/missing UPF '+element)
                data=tar.extractfile(found[0]).read()
                if hashlib.md5(data).hexdigest()!=item['md5']:raise RuntimeError('UPF checksum mismatch '+element)
                (target/item['filename']).write_bytes(data)
    report={'source':BASE+NAME,'archive_md5':DIGEST,'elements':{}}
    for e in ELEMENTS:
        p=target/meta[e]['filename'];assert digest(p)==meta[e]['md5']
        report['elements'][e]={**meta[e],'sha256':digest(p,'sha256')}
    (target/'provenance.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Verified official SSSP UPFs: '+', '.join(ELEMENTS))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',default=str(Path(__file__).resolve().parents[1]));a=p.parse_args();install(Path(a.root).resolve())
