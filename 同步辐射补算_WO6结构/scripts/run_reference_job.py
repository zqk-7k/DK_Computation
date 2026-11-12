#!/usr/bin/env python3
"""Run one finite-reference FEFF task using an existing complete FEFF wrapper."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import shutil
import subprocess
import time
import uuid
from datetime import datetime, timezone

from prepare_reference_feff import ROOT, read_config, sha, write_json


class SCFValidationError(ValueError):
    """Numeric outputs exist, but the self-consistent potential is not accepted."""

    def __init__(self, message, status, validation):
        super().__init__(message)
        self.status = status
        self.validation = validation


def verify_scf_convergence(directory, started=None):
    """Recognize FEFF10 POT's actual success/failure markers, not file existence.

    Markers verified against src/POT/potsub.f90, commit
    d861cf8b64342e2d63a0aeb62c8b78faca7942c3, lines 498-502/547/567.
    They also occur in the earlier real FEFF10 log1.dat and stdout.log.
    """
    paths, evidence, problems = [], [], []
    for name in ("log1.dat", "stdout.log", "convergence.scf"):
        path = directory/name
        if not path.is_file() or not path.stat().st_size:
            continue
        if started is not None and path.stat().st_mtime < started-2:
            raise ValueError("Stale SCF convergence evidence: "+name)
        paths.append(path)
        lines = path.read_text(errors="replace").splitlines()
        if name == "convergence.scf":
            iterations = []
            for number, line in enumerate(lines, 1):
                fields = line.split()
                if len(fields) != 5 or not re.fullmatch(r"\d+", fields[0]):
                    continue
                try:
                    values = [float(x.replace("D", "E").replace("d", "e")) for x in fields]
                except ValueError:
                    continue
                if all(math.isfinite(x) for x in values) and values[0] > 0 and fields[-1] in {"0", "1"}:
                    iterations.append((number, line.strip(), fields[-1]))
            if iterations:
                number, line, passed = iterations[-1]
                item = {"file": name, "line": number, "text": line}
                if passed == "1":
                    evidence.append(item)
                else:
                    problems.append(item)
        else:
            for number, line in enumerate(lines, 1):
                if re.search(r"\bConvergence reached in\s+[1-9]\d*\s+iterations\.", line):
                    evidence.append({"file": name, "line": number, "text": line.strip()})
                if re.search(r"Convergence not reached|ran out of iterations", line, re.IGNORECASE):
                    problems.append({"file": name, "line": number, "text": line.strip()})
    result = {
        "status": "not_converged" if problems else ("confirmed" if evidence else "unverified"),
        "evidence": evidence, "problems": problems,
        "marker_source": "https://github.com/times-software/feff10/blob/d861cf8b64342e2d63a0aeb62c8b78faca7942c3/src/POT/potsub.f90",
    }
    return result, paths


def numeric_table(path, columns, minimum_rows=10):
    path = Path(path)
    if not path.is_file() or not path.stat().st_size:
        raise ValueError("Missing/empty output: "+str(path))
    rows = []
    for lineno, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        fields = line.split()
        if not fields or fields[0].startswith(("#", "*", "!")):
            continue
        try:
            row = [float(x.replace("D", "E").replace("d", "e")) for x in fields]
        except ValueError:
            if rows and re.match(r"^[+-]?(?:\d|\.\d)", fields[0]) and len(fields) >= columns:
                raise ValueError(f"Malformed numeric output {path.name}:{lineno}")
            continue
        if len(row) < columns:
            continue
        if not all(math.isfinite(x) for x in row):
            raise ValueError("Nonfinite output: "+path.name)
        rows.append(row)
    if len(rows) < minimum_rows or any(b[0] <= a[0] for a, b in zip(rows, rows[1:])):
        raise ValueError("Insufficient/nonincreasing data: "+path.name)
    return rows


def validate_outputs(directory, job, started=None):
    error = directory/".feff.error"
    if error.is_file() and error.stat().st_size:
        raise ValueError("FEFF left a nonempty .feff.error; inspect the native log")
    mode = job["mode"]
    main = directory/("xmu.dat" if mode == "xanes" else "chi.dat")
    rows = numeric_table(main, 6 if mode == "xanes" else 2)
    column = 3 if mode == "xanes" else 1
    yy = [r[column] for r in rows]
    if max(yy)-min(yy) <= 1e-20:
        raise ValueError("Constant/zero spectrum")
    paths = [main]
    report = {"mode": mode, "row_count": len(rows), "x_range": [rows[0][0], rows[-1][0]], "files": []}
    if mode == "xanes":
        low, high = job["absolute_window_eV"]
        if rows[0][0] > low+1e-5 or rows[-1][0] < high-1e-5:
            raise ValueError(f"Native xmu.dat does not cover {low:g}..{high:g} eV; rerun an expanded relative EGRID, do not pad or shift the curve")
        step = job["relative_egrid_eV"][2]
        if max(b[0]-a[0] for a, b in zip(rows, rows[1:])) > step*1.05+1e-5:
            raise ValueError("Native XANES grid has unexpected gaps")
        report["absolute_energy_coverage"] = "PASS: native xmu.dat column 1"
    else:
        if rows[0][0] > 2.5 or rows[-1][0] < 12.5:
            raise ValueError("EXAFS does not support the intended 3..12 A^-1 tapered transform")
        path_table = directory/"paths.dat"
        if not path_table.is_file() or not path_table.stat().st_size:
            raise ValueError("EXAFS has no paths.dat")
        paths.append(path_table)
        scattering = sorted(directory.glob("feff[0-9][0-9][0-9][0-9].dat"))
        if not scattering:
            raise ValueError("EXAFS has no scattering paths; PRINT final field must be >=3")
        for p in scattering:
            numeric_table(p, 7, minimum_rows=5)
        paths += scattering
        report["scattering_path_count"] = len(scattering)
    scf, scf_paths = verify_scf_convergence(directory, started)
    paths += scf_paths
    report["scf_convergence"] = scf
    for p in paths:
        if started is not None and p.stat().st_mtime < started-2:
            raise ValueError("Stale FEFF output: "+p.name)
        report["files"].append({"name": p.name, "sha256": sha(p), "bytes": p.stat().st_size})
    if scf["status"] != "confirmed":
        status = ("outputs_validated_scf_unverified" if scf["status"] == "unverified"
                  else "outputs_validated_scf_not_converged")
        raise SCFValidationError(
            "Numeric spectra validated, but SCF convergence is "+scf["status"]+
            "; inspect native POT logs and rerun or obtain manual verification. This task is not complete.",
            status, report)
    return report


def archive_previous_run(directory, run_id):
    """Preserve every non-input file/directory; never reuse a failed native file."""
    directory = directory.resolve()
    keep = {"feff.inp", "job.json", ".run.lock", "previous_runs"}
    prior = [p for p in directory.iterdir() if p.name not in keep]
    if not prior:
        return
    archived = directory/"previous_runs"/run_id
    if not archived.resolve().is_relative_to(directory):
        raise ValueError("Archive directory escapes the job directory")
    # Resolve all paths before moving anything, including possible links/junctions.
    for path in prior:
        if not path.resolve().is_relative_to(directory):
            raise ValueError("Refusing to archive a path outside the job directory: "+path.name)
        if not (archived/path.name).resolve().is_relative_to(directory):
            raise ValueError("Archive destination escapes the job directory: "+path.name)
    archived.mkdir(parents=True, exist_ok=False)
    for path in prior:
        shutil.move(str(path), str(archived/path.name))


def command_for(config, root):
    command = config.get("software", {}).get("feff_command", [])
    if command:
        if not isinstance(command, list) or not all(isinstance(x, str) and x for x in command):
            raise ValueError("software.feff_command must be an argv list")
        return [x.replace("{root}", str(root)) for x in command]
    wrapper = os.environ.get("FEFF_RUNNER") or shutil.which("run_feff10.sh") or shutil.which("feff10.sh")
    if not wrapper:
        raise ValueError("Set software.feff_command or FEFF_RUNNER to an existing complete FEFF10 wrapper; no DFT or software installation is required by this runner")
    return ["bash", str(Path(wrapper).expanduser().resolve())]


def complete_current(directory, job, command):
    p = directory/"run_result.json"
    if not p.is_file():
        return False
    try:
        result = json.loads(p.read_text(encoding="utf-8"))
        if result.get("status") != "complete" or result.get("input_sha256") != sha(directory/"feff.inp") or result.get("command") != command:
            return False
        files = result.get("validation", {}).get("files", [])
        if not files or any(not (directory/f["name"]).is_file() or sha(directory/f["name"]) != f["sha256"] for f in files):
            return False
        validate_outputs(directory, job)
        return True
    except (ValueError, KeyError, OSError):
        return False


def run(config_path=None, task_index=0, validate_only=False):
    config_path, config = read_config(config_path)
    root = config_path.parent
    jobs = json.loads((root/"jobs_reference/manifest.json").read_text(encoding="utf-8"))
    if not 0 <= task_index < len(jobs):
        raise ValueError("Task index out of range")
    job = jobs[task_index]
    directory = (root/job["path"]).resolve()
    if not directory.is_relative_to(root.resolve()):
        raise ValueError("Job directory escapes project")
    if sha(root/job["source_xyz"]) != job["source_sha256"] or sha(directory/"feff.inp") != job["feff_input_sha256"]:
        raise ValueError("Reference source/input changed; regenerate explicitly in a new project")
    if validate_only:
        print(json.dumps(validate_outputs(directory, job), indent=2))
        return
    command = command_for(config, root)
    if complete_current(directory, job, command):
        print("Already complete with verified current hashes:", job["id"])
        return
    lock = directory/".run.lock"
    fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    os.write(fd, f"pid={os.getpid()}\n".encode()); os.close(fd)
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:6]
    started = time.time()
    result = {"status": "running", "id": job["id"], "mode": job["mode"], "run_id": run_id,
              "command": command, "input_sha256": sha(directory/"feff.inp"), "source_sha256": job["source_sha256"],
              "started_utc": datetime.now(timezone.utc).isoformat(), "engine": "external complete FEFF wrapper"}
    process = None
    try:
        archive_previous_run(directory, run_id)
        write_json(directory/"run_result.json", result)
        timeout = float(config.get("software", {}).get("feff_timeout_seconds", 252000))
        if not math.isfinite(timeout) or timeout <= 0:
            raise ValueError("feff_timeout_seconds must be finite and positive")
        with (directory/"stdout.log").open("wb") as stdout, (directory/"stderr.log").open("wb") as stderr:
            process = subprocess.Popen(command, cwd=directory, stdout=stdout, stderr=stderr,
                                       start_new_session=(os.name == "posix"))
            process.wait(timeout=timeout)
        result["return_code"] = process.returncode
        if process.returncode:
            raise RuntimeError(f"FEFF wrapper returned {process.returncode}; see stdout.log and stderr.log")
        result["validation"] = validate_outputs(directory, job, started)
        result["status"] = "complete"
    except BaseException as exc:
        result["status"] = exc.status if isinstance(exc, SCFValidationError) else "failed"
        if isinstance(exc, SCFValidationError):
            result["validation"] = exc.validation
        result["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        if process is not None and process.poll() is None:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGTERM)
            else:
                process.terminate()
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                if os.name == "posix":
                    os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
                process.wait(timeout=10)
        result["duration_seconds"] = time.time()-started
        result["finished_utc"] = datetime.now(timezone.utc).isoformat()
        write_json(directory/"run_result.json", result)
        lock.unlink(missing_ok=True)
    print("Validated FEFF output:", job["id"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--task-index", type=int, required=True)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    run(args.config, args.task_index, args.validate_only)


if __name__ == "__main__":
    main()
