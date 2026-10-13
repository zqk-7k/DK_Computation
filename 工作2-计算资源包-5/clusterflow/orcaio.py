"""Small, auditable ORCA 6.1 output parser and single-job runner.

Parsing patterns were informed by the user's ORCA_52_NMR_HPC_portable_v1
process_nmr_predictions.py and verified against pre-existing ORCA 6.1 outputs.
No legacy regression coefficients, target-specific assignments or absolute
project paths are imported. This module does not validate a chemical structure.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import time
import uuid
from datetime import datetime, timezone

FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def _number(value: str) -> float:
    value = float(value.replace("D", "E").replace("d", "e"))
    if not math.isfinite(value):
        raise ValueError("Non-finite ORCA numeric value")
    return value


def _final_run(text: str) -> str:
    # Appended runs must never inherit a previous run's success or properties.
    starts = list(re.finditer(r"^\s*\*\s*O\s+R\s+C\s+A\s*\*\s*$", text, re.M))
    return text[starts[-1].start():] if starts else text


def _last_value(text: str, marker: str, pattern: str) -> float | None:
    """Read the final marker, never fall back when its value is truncated."""
    locations = list(re.finditer(marker, text, re.I | re.M))
    if not locations:
        return None
    match = re.match(pattern, text[locations[-1].start():], re.I)
    return _number(match.group(1)) if match else None


def _frequency_block(text: str, warnings: list[str]) -> list[float]:
    starts = list(re.finditer(r"^\s*VIBRATIONAL FREQUENCIES\s*$", text, re.M))
    if not starts:
        return []
    block = text[starts[-1].end():]
    # A truncated final block is not replaced by an earlier complete one.
    end = re.search(r"^\s*(?:NORMAL MODES|THERMOCHEMISTRY AT|IR SPECTRUM)\b", block, re.M)
    if end is None:
        warnings.append("Final frequency block has no closing section; frequencies withheld.")
        return []
    block = block[:end.start()]
    rows = re.findall(rf"^\s*(\d+):\s*({FLOAT})\s+cm(?:\*\*|\^)?-1\b[^\n]*", block, re.M)
    if not rows or [int(i) for i, _ in rows] != list(range(len(rows))):
        warnings.append("Final frequency block is empty or non-contiguous; frequencies withheld.")
        return []
    return [_number(v) for _, v in rows]


def _shielding_block(text: str, warnings: list[str]) -> list[dict]:
    markers = list(re.finditer(r"CHEMICAL SHIELDING SUMMARY\s*\(ppm\)", text))
    if not markers:
        return []
    block = text[markers[-1].end():]
    rows: list[dict] = []
    seen: set[int] = set()
    started = False
    closed = False
    pattern = re.compile(rf"^\s*(\d+)\s+([A-Z][a-z]?)\s+({FLOAT})\s+({FLOAT})\s*$")
    for line in block.splitlines():
        match = pattern.match(line)
        if match:
            started = True
            idx = int(match.group(1))
            if idx in seen:
                warnings.append("Duplicate shielding atom indices; final shielding block withheld.")
                return []
            seen.add(idx)
            rows.append({"atom_index_zero_based": idx, "element": match.group(2),
                         "isotropic_ppm": _number(match.group(3))})
        elif started:
            if line.strip() and re.match(r"^\s*\d+\s+[A-Z]", line):
                warnings.append("Truncated shielding row; final shielding block withheld.")
                return []
            closed = True
            break
    if not rows or not closed:
        warnings.append("Final shielding block is incomplete; shielding withheld.")
        return []
    return rows


def _property_array(section: str, name: str, integer: bool = False) -> list[list]:
    declaration = re.search(rf"&{re.escape(name)}\s+\[[^\n]*?&Dim\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)[^\n]*\][^\n]*\n", section)
    if declaration is None:
        raise ValueError(f"Missing property array {name}")
    n_rows, n_cols = map(int, declaration.groups())
    values = {}
    for line in section[declaration.end():].splitlines():
        stripped = line.strip()
        if stripped.startswith(("&", "$")):
            break
        parts = stripped.split()
        # Column-heading rows contain only n_cols entries, not n_cols + 1.
        if len(parts) != n_cols + 1 or not re.fullmatch(r"\d+", parts[0]):
            continue
        index = int(parts[0])
        if not 0 <= index < n_rows or index in values:
            raise ValueError(f"Invalid/duplicate row index in {name}")
        values[index] = [int(v) if integer else _number(v) for v in parts[1:]]
    if len(values) != n_rows:
        raise ValueError(f"Incomplete property array {name}")
    return [values[i] for i in range(n_rows)]


def _j_pairs(text: str | None, warnings: list[str]) -> list[dict]:
    if not text:
        return []
    marker = "$SCF_Spin_Spin_Coupling"
    start = text.rfind(marker)
    if start < 0:
        return []
    end = text.find("$End", start)
    if end < 0:
        warnings.append("Final spin-spin property section is incomplete; J withheld.")
        return []
    section = text[start:end]
    try:
        info = _property_array(section, "pairsInfo", integer=True)
        couplings = _property_array(section, "pairsTotalSSCIso")
        if len(info) != len(couplings):
            raise ValueError("J property array lengths differ")
        rows, seen = [], set()
        for atom, coupling in zip(info, couplings):
            if len(atom) != 4 or len(coupling) != 1:
                raise ValueError("Unexpected J property dimensions")
            i, j = sorted((atom[0], atom[2]))
            if i < 0 or i == j or (i, j) in seen:
                raise ValueError("Invalid or duplicate J atom pair")
            seen.add((i, j))
            rows.append({"atom_i_zero_based": i, "atom_j_zero_based": j,
                         "J_Hz": coupling[0]})
        return rows
    except (ValueError, OverflowError) as error:
        warnings.append(f"J parsing failed: {error}; J withheld.")
        return []


def parse_output(output_text: str, property_text: str | None = None) -> dict:
    """Parse final ORCA run. Missing data are None/[]; no inferred numbers.

    Frequencies include printed translation/rotation zeros. Shielding and J
    atom indices are zero-based as reported in output/property files, not the
    one-based indices of an input Nuclei list. Stage-specific acceptance and
    validation of the property/output association remain the caller's job.
    """
    text = _final_run(output_text)
    warnings: list[str] = []
    normal = bool(re.search(r"^\s*\**\s*ORCA TERMINATED NORMALLY\s*\**\s*$", text, re.M))
    opt_events = list(re.finditer(r"THE OPTIMIZATION HAS CONVERGED|OPTIMIZATION (?:DID NOT|HAS NOT) CONVERGE", text, re.I))
    # An OPTIMIZATION RUN DONE marker alone is not evidence of convergence.
    opt = bool(opt_events and "HAS CONVERGED" in opt_events[-1].group().upper())
    thermo_headers = list(re.finditer(rf"THERMOCHEMISTRY\s+AT\s+({FLOAT})\s*K", text, re.I))
    thermo = text[thermo_headers[-1].start():] if thermo_headers else ""
    temperature = _number(thermo_headers[-1].group(1)) if thermo_headers else None
    # A later frequency section invalidates earlier thermochemistry.
    freq_location = text.rfind("VIBRATIONAL FREQUENCIES")
    if thermo_headers and freq_location > thermo_headers[-1].start():
        thermo, temperature = "", None
    gibbs = _last_value(thermo, r"Final\s+Gibbs\s+free\s+energy", rf"Final\s+Gibbs\s+free\s+energy\s+\.\.\.\s*({FLOAT})\s+Eh")
    qrrho = re.findall(r"Quasi\s+RRHO\s+\.\.\.\s+(True|False)", thermo, re.I)
    cutoffs = re.findall(rf"reference\s+frequency\s+of\s+({FLOAT})\s+cm-1", thermo, re.I)
    return {
        "normal_termination": normal,
        "optimization_converged": opt,
        "final_energy_Eh": _last_value(text, r"FINAL\s+SINGLE\s+POINT\s+ENERGY", rf"FINAL\s+SINGLE\s+POINT\s+ENERGY\s+({FLOAT})"),
        "gibbs_free_energy_Eh": gibbs,
        "temperature_K": temperature,
        "frequencies_cm1": _frequency_block(text, warnings),
        "shielding": _shielding_block(text, warnings),
        "j_pairs": _j_pairs(property_text, warnings),
        "qrrho_enabled": qrrho[-1].lower() == "true" if qrrho else None,
        "qrrho_reference_frequency_cm1": _number(cutoffs[-1]) if cutoffs else None,
        "warnings": warnings,
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _input_hashes(directory: Path) -> dict[str, str]:
    paths = {"job.inp": directory / "job.inp", "job.json": directory / "job.json"}
    # Inline-coordinate jobs still use these artefacts for graph/atom-order
    # validation downstream, so they are immutable inputs too.
    for name in ("input.xyz", "input_structure.mol", "atom_map.json"):
        if (directory / name).is_file():
            paths[name] = directory / name
    text = paths["job.inp"].read_text(encoding="utf-8")
    references = re.findall(r'^\s*\*\s*xyzfile\s+[-+]?\d+\s+\d+\s+(.+?)\s*$', text, re.M | re.I)
    references += re.findall(r'^\s*%moinp\s+(.+?)\s*$', text, re.M | re.I)
    for raw in references:
        name = raw.split("#", 1)[0].strip().strip('"\'')
        path = (directory / name).resolve()
        # Portable jobs are self-contained. Do not read arbitrary outside files.
        try:
            relative = path.relative_to(directory)
        except ValueError as error:
            raise ValueError(f"Input reference escapes job directory: {name}") from error
        if not path.is_file():
            raise FileNotFoundError(f"Missing referenced input: {path}")
        paths[relative.as_posix()] = path
    return {name: _sha256(path) for name, path in sorted(paths.items())}


def _archive_existing(directory: Path, input_names: set[str], attempt: str) -> None:
    files = [p for p in directory.iterdir() if p.is_file() and
             (p.name.startswith("job.") or p.name in {"receipt.json", "propertyjob.property.txt"})
             and p.name not in input_names]
    if not files:
        return
    archive = directory / "attempts" / attempt
    archive.mkdir(parents=True, exist_ok=False)
    for path in files:
        path.rename(archive / path.name)


def _terminate_process(process: subprocess.Popen) -> None:
    """Stop only the child/session created by this runner, not unrelated jobs."""
    try:
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGTERM)
        else:
            process.terminate()
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        try:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGKILL)
            else:
                process.kill()
        except ProcessLookupError:
            pass
        process.wait()


def run_job(job_dir: Path, orca_executable: Path, timeout_seconds: float | None = None,
            resume: bool = True) -> dict:
    """Run exactly one local job via the serial ORCA driver, without shell/MPI.

    Same-input successful runs resume only if all recorded output hashes match.
    Same-input failures are NOT retried with resume=True; an explicit
    resume=False is required. Existing products/logs are archived, not deleted.
    Scheduler allocation and scientific acceptance are independent checks.
    """
    directory = Path(job_dir).resolve(strict=True)
    executable = Path(orca_executable).resolve(strict=True)
    if not directory.is_dir() or not executable.is_file():
        raise ValueError("job_dir must be a directory and ORCA must be an executable file")
    if timeout_seconds is not None and timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive")
    lock = directory / ".orca-run.lock"
    token = uuid.uuid4().hex
    lock_payload = {"token": token, "pid": os.getpid(), "started_utc": _utc()}
    # Never remove an existing/stale lock automatically: ownership is unknown.
    with lock.open("x", encoding="utf-8") as handle:
        json.dump(lock_payload, handle)
    try:
        inputs = _input_hashes(directory)
        metadata = json.loads((directory / "job.json").read_text(encoding="utf-8"))
        if not isinstance(metadata, dict) or not metadata.get("id"):
            raise ValueError("job.json must be an object with a nonempty id")
        for key in ("nprocs", "maxcore_mb"):
            if key in metadata and (not isinstance(metadata[key], int) or metadata[key] < 1):
                raise ValueError(f"Invalid {key} in job.json")
        executable_hash = _sha256(executable)
        receipt_path = directory / "receipt.json"
        previous = None
        if receipt_path.exists():
            previous = json.loads(receipt_path.read_text(encoding="utf-8"))
        matching = bool(previous and previous.get("input_sha256") == inputs and
                        previous.get("orca_executable_sha256") == executable_hash)
        if resume and matching:
            if previous.get("status") == "completed":
                outputs = previous.get("output_sha256", {})
                intact = bool(outputs) and all((directory / name).resolve().is_relative_to(directory)
                         and (directory / name).is_file() and _sha256(directory / name) == value
                         for name, value in outputs.items())
                if intact:
                    return {**previous, "skipped": True, "resume_reason": "matching_inputs_and_outputs"}
                raise ValueError("Completed-job output hashes changed; inspect and use resume=False explicitly")
            if previous.get("status") == "failed":
                return {**previous, "skipped": True, "resume_reason": "prior_failure_requires_resume_false"}
        attempt = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S") + "_" + token[:10]
        _archive_existing(directory, set(inputs), attempt)
        receipt = {"schema_version": "1.0", "id": metadata["id"], "stage": metadata.get("stage"),
                   "attempt": attempt, "started_utc": _utc(), "status": "running", "skipped": False,
                   "input_sha256": inputs, "orca_executable": str(executable),
                   "orca_executable_sha256": executable_hash, "command": [str(executable), "job.inp"]}
        _atomic_json(receipt_path, receipt)
        started = time.monotonic()
        returncode, failure, timed_out = None, None, False
        interruption = None
        process = None
        try:
            with (directory / "job.out").open("wb") as stdout, (directory / "job.stderr").open("wb") as stderr:
                process = subprocess.Popen([str(executable), "job.inp"], cwd=directory,
                                           stdout=stdout, stderr=stderr, shell=False,
                                           start_new_session=(os.name == "posix"))
                try:
                    returncode = process.wait(timeout=timeout_seconds)
                except subprocess.TimeoutExpired:
                    timed_out = True
                    failure = "timeout"
                    _terminate_process(process)
                    returncode = process.returncode
                except KeyboardInterrupt as error:
                    interruption = error
                    failure = "keyboard_interrupt"
                    _terminate_process(process)
                    returncode = process.returncode
        except OSError as error:
            failure = f"launch_error: {error}"
        property_path = next((directory / name for name in ("job.property.txt", "propertyjob.property.txt")
                              if (directory / name).is_file()), None)
        parsed = parse_output((directory / "job.out").read_text(encoding="utf-8", errors="replace")
                              if (directory / "job.out").exists() else "",
                              property_path.read_text(encoding="utf-8", errors="replace") if property_path else None)
        success = returncode == 0 and parsed["normal_termination"] and not timed_out and failure is None
        if not success and failure is None:
            failure = "nonzero_exit" if returncode != 0 else "missing_normal_termination"
        for name in metadata.get("required_outputs", []):
            path = (directory / name).resolve()
            if not path.is_relative_to(directory) or not path.is_file() or path.stat().st_size == 0:
                success, failure = False, f"missing_required_output: {name}"
        outputs = {p.name: _sha256(p) for p in directory.iterdir() if p.is_file()
                   and (p.name.startswith("job.") or p.name == "propertyjob.property.txt")
                   and p.name not in inputs}
        receipt.update({"status": "completed" if success else "failed", "finished_utc": _utc(),
                        "elapsed_seconds": time.monotonic() - started, "returncode": returncode,
                        "failure_reason": failure, "timed_out": timed_out,
                        "output_sha256": outputs, "parsed": parsed,
                        "scientific_acceptance": "NOT_EVALUATED_BY_RUNNER"})
        _atomic_json(receipt_path, receipt)
        if interruption is not None:
            raise interruption
        return receipt
    finally:
        # Only the process that wrote this token may remove this particular lock.
        try:
            if json.loads(lock.read_text(encoding="utf-8")).get("token") == token:
                lock.unlink()
        except (FileNotFoundError, ValueError, OSError):
            pass
