"""Shared immutable identity, configuration and confined-path helpers."""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import uuid

def _io_path(path):
    """Use Windows extended paths for IO only, without changing stored identities."""
    path=Path(path)
    if os.name!='nt':return path
    name=os.path.abspath(str(path))
    if name.startswith('\\\\?\\'):return Path(name)
    if name.startswith('\\\\'):return Path('\\\\?\\UNC\\'+name[2:])
    return Path('\\\\?\\'+name)

def readj(path):
    return json.loads(_io_path(path).read_text(encoding='utf-8-sig'))

def writej(path,obj):
    path=_io_path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name('.write-'+uuid.uuid4().hex+'.tmp')
    try:
        temporary.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        os.replace(temporary,path)
    finally:
        if temporary.exists():temporary.unlink()

def sha(path):
    value=hashlib.sha256()
    with _io_path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):value.update(b)
    return value.hexdigest()

def hash_protocol(config):
    return hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()[:16]

protocol_hash=hash_protocol
sha256_file=sha
read_json=readj
write_json=writej

def config(root):
    return readj(Path(root)/'config/protocol.json')

def inside(root,name):
    root=Path(root).resolve(); relative=Path(name)
    path=(root/relative).resolve()
    if relative.is_absolute() or path==root or not path.is_relative_to(root):
        raise ValueError('Path must be relative and confined within package')
    return path
