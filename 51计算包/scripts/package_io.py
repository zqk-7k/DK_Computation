#!/usr/bin/env python3
"""Small, dependency-free I/O primitives shared by the HPC tools."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any, Iterable


TERMINAL_STATES = frozenset({"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT"})
ACTIVE_STATES = frozenset({"RELEASED", "PENDING", "RUNNING", "COMPLETING"})
KNOWN_STATES = TERMINAL_STATES | ACTIVE_STATES | {"HELD"}


def atomic_write_text(path: Path, text: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_name, path)
    finally:
        try:
            os.unlink(tmp_name)
        except FileNotFoundError:
            pass


def atomic_write_json(path: Path, value: Any) -> None:
    atomic_write_text(Path(path), json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n")


def load_json(path: Path) -> Any:
    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def safe_relative(path: Path, root: Path) -> str:
    return Path(path).resolve().relative_to(Path(root).resolve()).as_posix()


def iter_package_files(root: Path, excluded_roots: Iterable[str] = ()):
    root = Path(root).resolve()
    excluded = {str(x).replace("\\", "/").strip("/") for x in excluded_roots}
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        rel = safe_relative(path, root)
        if any(rel == x or rel.startswith(x + "/") for x in excluded):
            continue
        yield path, rel
