#!/usr/bin/env python3
"""Package the source/input bundle; exclude Python caches and calculation output."""
import hashlib
from pathlib import Path
import zipfile

root = Path(__file__).resolve().parents[1]
target = root.with_suffix(".zip")
if target.exists():
    raise FileExistsError(f"Preserve earlier delivery: {target}")
files = sorted(p for p in root.rglob("*") if p.is_file() and not any(
    part in {"__pycache__", ".pytest_cache", ".venv", "results"} for part in p.relative_to(root).parts))
manifest = root/"PACKAGE_CONTENTS_SHA256.txt"
manifest.write_text("\n".join(hashlib.sha256(p.read_bytes()).hexdigest()+"  "+p.relative_to(root).as_posix()
                             for p in files if p != manifest)+"\n", encoding="utf-8")
if manifest not in files:
    files.append(manifest)
with zipfile.ZipFile(target,"x",compression=zipfile.ZIP_DEFLATED) as archive:
    for path in files:
        archive.write(path, arcname=root.name+"/"+path.relative_to(root).as_posix())
print(f"Created {target.name}: {len(files)} files; {target.stat().st_size} bytes")
print("SHA256", hashlib.sha256(target.read_bytes()).hexdigest())
