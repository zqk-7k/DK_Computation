#!/usr/bin/env python3
"""Run one FEFF job, preserve logs, and verify mode-specific numerical outputs.

Only Python's standard library is required. An external command is an argv list,
never a shell string. Successful completion is run_result.json status=complete.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parents[1]
PINNED_COMMIT = "d861cf8b64342e2d63a0aeb62c8b78faca7942c3"


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path: Path, value: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def read_cards(path: Path) -> dict[str, list[str]]:
    cards = {}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(("*", "#", "!")):
            continue
        fields = stripped.split()
        if fields[0].upper() == "END":
            break
        if re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", fields[0]):
            cards[fields[0].upper()] = fields[1:]
    return cards


def numeric_table(path: Path, columns: int, min_rows: int = 10) -> tuple[list[list[float]], dict]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing/empty output: {path.name}")
    rows = []
    # FEFF writes descriptive headers, then a dashed separator/column headers.
    # Header lines beginning with # are skipped; only fully numerical rows count.
    for number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        line = line.strip()
        if not line or line.startswith(("#", "*", "!")):
            continue
        fields = line.split()
        try:
            values = [float(s.replace("D", "E").replace("d", "e")) for s in fields]
        except ValueError:
            # A numeric data line with overflow stars or NaN text must not vanish.
            if rows and fields and re.match(r"^[+-]?(?:\d|\.\d)", fields[0]) and len(fields) >= columns:
                raise ValueError(f"Malformed numeric output in {path.name}:{number}")
            continue
        if len(values) < columns:
            continue
        if not all(math.isfinite(v) for v in values):
            raise ValueError(f"Non-finite value in {path.name}:{number}")
        rows.append(values)
    if len(rows) < min_rows:
        raise ValueError(f"{path.name}: only {len(rows)} numerical rows; require {min_rows}")
    x = [r[0] for r in rows]
    if any(b <= a for a, b in zip(x, x[1:])):
        raise ValueError(f"{path.name}: first column must increase strictly")
    return rows, {"name": path.name, "rows": len(rows), "columns": len(rows[0]),
                  "x_min": min(x), "x_max": max(x), "bytes": path.stat().st_size,
                  "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def check_outputs(job: Path, mode: str, started: float) -> dict:
    error_file = job / ".feff.error"
    if error_file.exists() and error_file.stat().st_size:
        raise ValueError("FEFF left a nonempty .feff.error: " + error_file.read_text(errors="replace")[:1000])
    details = {"mode": mode, "files": [], "checks": []}
    main = job / ("xmu.dat" if mode == "xanes" else "chi.dat")
    rows, summary = numeric_table(main, 6 if mode == "xanes" else 2)
    if main.stat().st_mtime < started - 2:
        raise ValueError(f"Stale output timestamp: {main.name}")
    details["files"].append(summary)
    y = [r[3] if mode == "xanes" else r[1] for r in rows]
    if max(abs(v) for v in y) <= 1e-20 or max(y) - min(y) <= 1e-20:
        raise ValueError(f"{main.name}: spectrum is zero/constant")
    details["checks"].append("finite, increasing grid, nonconstant spectrum")
    if mode == "exafs":
        paths = job / "paths.dat"
        if not paths.is_file() or paths.stat().st_size == 0 or paths.stat().st_mtime < started - 2:
            raise ValueError("EXAFS requires fresh nonempty paths.dat")
        details["files"].append({"name": paths.name, "bytes": paths.stat().st_size,
                                 "sha256": hashlib.sha256(paths.read_bytes()).hexdigest()})
        files = sorted(job.glob("feff[0-9][0-9][0-9][0-9].dat"))
        if not files:
            raise ValueError("No feffNNNN.dat paths: set six-entry PRINT with final value >= 3")
        for path in files:
            if path.stat().st_mtime < started - 2:
                raise ValueError(f"Stale path output: {path.name}")
            _, info = numeric_table(path, 7, min_rows=5)
            details["files"].append(info)
        details["path_count"] = len(files)
        details["checks"].append("paths.dat and all feffNNNN.dat numerical tables validated")
    return details


def archive_previous_outputs(job: Path, run_id: str) -> str | None:
    # Only known generated outputs/logs are archived; feff.inp and job.json stay.
    exact = {"xmu.dat", "chi.dat", "paths.dat", "feff.bin", ".feff.error",
             "stdout.log", "stderr.log", "run_result.json"}
    candidates = [p for p in job.iterdir() if p.is_file() and
                  (p.name in exact or re.fullmatch(r"feff\d{4}\.dat", p.name))]
    if not candidates:
        return None
    archive = job / "previous_runs" / run_id
    archive.mkdir(parents=True, exist_ok=False)
    for path in candidates:
        shutil.move(str(path), str(archive / path.name))
    return str(archive.relative_to(job))


def resolve_command(config: dict) -> list[str]:
    command = config.get("software", {}).get("feff_command", [])
    if command in (None, []):
        wrapper = ROOT / "software" / "feff10" / "bin" / "run_feff10.sh"
        if not wrapper.is_file():
            raise ValueError(f"FEFF is not installed. Run scripts/install_feff10.sh; missing {wrapper}")
        return ["bash", str(wrapper)]
    if not isinstance(command, list) or not command or any(not isinstance(x, str) or not x for x in command):
        raise ValueError("software.feff_command must be a nonempty string array or [] (never a shell string)")
    expanded = [s.replace("{root}", str(ROOT)) for s in command]
    if ("/" in expanded[0] or "\\" in expanded[0]) and not Path(expanded[0]).is_absolute():
        expanded[0] = str((ROOT / expanded[0]).resolve())
    return expanded


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    job = args.job.resolve()
    if not job.is_dir():
        parser.error(f"Job directory does not exist: {job}")
    lock = job / ".workflow_run.lock"
    try:
        fd = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(fd, "w") as handle:
            handle.write(f"pid={os.getpid()}\nhost={os.uname().nodename if hasattr(os, 'uname') else 'unknown'}\nstart={utc_now()}\n")
    except FileExistsError:
        print(f"Job locked: {lock}; inspect the active process before removing a stale lock.", file=sys.stderr)
        return 2
    run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
    result = {"status": "running", "run_id": run_id, "job_directory": str(job),
              "started_utc": utc_now(), "return_code": None, "engine": "FEFF", "validation": None}
    proc = None
    started = time.time()
    def interrupt(signum, frame):
        raise KeyboardInterrupt(f"Interrupted by signal {signum}")
    old_handlers = {}
    for sig in (signal.SIGTERM, signal.SIGINT):
        old_handlers[sig] = signal.signal(sig, interrupt)
    try:
        config = json.loads(args.config.read_text(encoding="utf-8-sig"))
        metadata = json.loads((job / "job.json").read_text(encoding="utf-8-sig"))
        mode = str(metadata.get("mode", "")).lower()
        if mode not in {"xanes", "exafs"}:
            raise ValueError("job.json must specify mode='xanes' or mode='exafs'")
        inp = job / "feff.inp"
        result.update(mode=mode, input_sha256=hashlib.sha256(inp.read_bytes()).hexdigest())
        cards = read_cards(inp)
        if mode.upper() not in cards:
            raise ValueError(f"job.json mode={mode}, but feff.inp lacks {mode.upper()}")
        if mode == "exafs":
            values = cards.get("PRINT", [])
            if len(values) < 6 or int(values[5]) < 3:
                raise ValueError("EXAFS requires PRINT with six values and final value >=3, e.g. PRINT 1 0 0 0 0 3")
        command = resolve_command(config)
        result["command"] = command
        install_manifest = ROOT / "software" / "feff10" / "install_complete.json"
        if not config.get("software", {}).get("feff_command") and install_manifest.exists():
            result["software_installation"] = json.loads(install_manifest.read_text())
        else:
            result["software_installation"] = {"external_command": True, "version": "recorded in engine stdout/header; verify before production"}
        result["previous_outputs_archive"] = archive_previous_outputs(job, run_id)
        write_json(job / "run_result.json", result)
        env = os.environ.copy()
        env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        software = config.get("software", {})
        timeout = software.get("feff_timeout_seconds", software.get("timeout_seconds", 0))
        timeout = float(timeout) if timeout else None
        if timeout is not None and (timeout <= 0 or not math.isfinite(timeout)):
            raise ValueError("software.feff_timeout_seconds must be positive, or 0 for no runner timeout")
        started = time.time()
        with (job / "stdout.log").open("w", encoding="utf-8") as out, (job / "stderr.log").open("w", encoding="utf-8") as err:
            proc = subprocess.Popen(command, cwd=job, stdout=out, stderr=err, env=env,
                                    start_new_session=(os.name == "posix"))
            result["return_code"] = proc.wait(timeout=timeout)
        if result["return_code"] != 0:
            raise RuntimeError(f"FEFF command exited with code {result['return_code']}; inspect stdout.log and stderr.log")
        result["validation"] = check_outputs(job, mode, started)
        result["status"] = "complete"
        return_value = 0
    except (Exception, KeyboardInterrupt) as exc:
        result.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        if proc is not None and proc.poll() is None:
            try:
                if os.name == "posix":
                    os.killpg(proc.pid, signal.SIGTERM)
                else:
                    proc.terminate()
            except ProcessLookupError:
                pass
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                try:
                    if os.name == "posix":
                        os.killpg(proc.pid, signal.SIGKILL)
                    else:
                        proc.kill()
                except ProcessLookupError:
                    pass
                proc.wait()
            result["return_code"] = proc.returncode
        print(result["error"], file=sys.stderr)
        return_value = 1
    finally:
        duration = round(time.time() - started, 3)
        result.update(finished_utc=utc_now(), duration_seconds=duration, elapsed_seconds=duration)
        write_json(job / "run_result.json", result)
        lock.unlink(missing_ok=True)
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)
    print(json.dumps({"job": str(job), "status": result["status"], "result": str(job / "run_result.json")}))
    return return_value


if __name__ == "__main__":
    sys.exit(main())
