"""Read-only preflight and explicit, one-job GOAT launcher. No batch submission.

Preflight never starts ORCA, xTB or MPI. Binary markers are not a runtime smoke
test. A normal run still awaits independent chemical/ensemble acceptance.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
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
import uuid

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
REQUIRED_VERSION = "6.1.0"
INPUT_HASH_FIELDS = {
    "goat.inp": "input_sha256",
    "start.xyz": "xyz_sha256",
    "atom_map.json": "atom_map_sha256",
    "cluster.mol": "cluster_mol_sha256",
}
XTB_INSTALL_DOC = "https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/semiempirical.html"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, data: dict) -> None:
    temporary = path.with_name("." + path.name + "." + uuid.uuid4().hex + ".tmp")
    try:
        temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def embedded_versions(path: Path) -> list[str]:
    """Inspect the executable without executing it, including chunk boundaries."""
    found: set[str] = set()
    tail = b""
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            data = tail + block
            found.update(v.decode("ascii") for v in re.findall(
                rb"Program Version[ \t]+(\d+\.\d+\.\d+)(?![\d.])", data))
            tail = data[-128:]
    return sorted(found)


def executable_file(path: Path) -> bool:
    return path.is_file() and (os.name == "nt" or os.access(path, os.X_OK))


def preflight(orca: Path, xtb: Path, nprocs: int = 8) -> dict:
    orca, xtb = Path(orca).resolve(), Path(xtb).resolve()
    from clusterflow.runtime import wrapper_preflight
    wrapper = wrapper_preflight(orca, xtb)
    if wrapper is not None:
        return {**wrapper, 'orca_executable': str(orca), 'xtb_executable': str(xtb),
                'orca_sha256': digest(orca) if orca.is_file() else None,
                'xtb_sha256': digest(xtb) if xtb.is_file() else None,
                'quantum_jobs_launched': 0, 'nprocs': nprocs}
    reasons: list[str] = []
    versions = embedded_versions(orca) if orca.is_file() else []
    if not executable_file(orca):
        reasons.append("ORCA executable missing or not executable")
    if versions != [REQUIRED_VERSION]:
        reasons.append("Exact unique embedded ORCA Program Version 6.1.0 not verified")
    if not executable_file(xtb):
        reasons.append("External xTB executable missing or not executable")
    installed_xtb = orca.parent / ("otool_xtb.exe" if os.name == "nt" else "otool_xtb")
    if xtb != installed_xtb.resolve():
        reasons.append("--xtb must identify otool_xtb beside the ORCA driver; no binaries are copied automatically")
    if os.name == "nt" and not (orca.parent / "libiomp5md.dll").is_file():
        reasons.append("Windows xTB runtime libiomp5md.dll missing beside ORCA")
    mpi = shutil.which("mpiexec") or shutil.which("mpirun")
    if nprocs > 1 and not mpi:
        reasons.append("MPI launcher unavailable in this process PATH")
    if not isinstance(nprocs, int) or nprocs < 1:
        reasons.append("Invalid MPI process count")
    return {
        "status": "PREFLIGHT_STATIC_PASS" if not reasons else "PREFLIGHT_BLOCKED",
        "passed": not reasons, "reasons": reasons,
        "orca_executable": str(orca), "orca_sha256": digest(orca) if orca.is_file() else None,
        "orca_required_version": REQUIRED_VERSION, "embedded_versions": versions,
        "version_evidence": "READ_ONLY_BINARY_MARKER_NOT_RUNTIME_VERIFICATION",
        "xtb_executable": str(xtb), "xtb_sha256": digest(xtb) if xtb.is_file() else None,
        "xtb_version_status": "NOT_EXECUTED_NOT_VERIFIED_REQUIRE_6.7.1_OR_LATER",
        "xtb_install_documentation": XTB_INSTALL_DOC,
        "mpi_launcher": mpi, "mpi_build_compatibility": "SITE_VALIDATION_REQUIRED",
        "nprocs": nprocs, "runtime_goat_smoke": "UNVALIDATED_NOT_EXECUTED",
        "quantum_jobs_launched": 0,
    }


def xyz_frames(path: Path) -> list[list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    frames: list[list[str]] = []
    offset = 0
    while offset < len(lines):
        if not lines[offset].strip():
            offset += 1
            continue
        try:
            count = int(lines[offset])
        except ValueError as exc:
            raise ValueError("Invalid XYZ atom count") from exc
        if count < 1 or offset + count + 2 > len(lines):
            raise ValueError("Truncated or empty XYZ frame")
        symbols = []
        for line in lines[offset + 2:offset + count + 2]:
            fields = line.split()
            if len(fields) != 4 or not re.fullmatch(r"[A-Z][a-z]?", fields[0]):
                raise ValueError("XYZ requires exactly element/x/y/z, without isotope mass tokens")
            if not all(math.isfinite(float(v)) for v in fields[1:]):
                raise ValueError("XYZ has nonfinite coordinates")
            symbols.append(fields[0])
        frames.append(symbols)
        offset += count + 2
    if not frames:
        raise ValueError("No XYZ frames")
    return frames


def controller_released(root: Path, job_id: str, species_id: str) -> bool:
    """Held starts need a current, source-bound controller release, never a glob."""
    settings_path = root / 'config/automation.json'
    plan_path = root / 'registry/batch_plan_v5.json'
    if not settings_path.is_file() or not plan_path.is_file():
        return False
    if read_json(settings_path).get('automatic_handoff',{}).get('enabled') is not True:
        return False
    names = ('registry/batch_plan_v5.json','registry/jobs.json','registry/species.json',
             'config/protocol.json','config/automation.json')
    identity = {name:digest(root/name) for name in names}
    plans = {r['species_id']:r for r in read_json(plan_path)['species']}
    if species_id not in plans:
        return False
    for state_path in (root/'automation').glob('*/state.json'):
        if not state_path.resolve().is_relative_to(root.resolve()):
            continue
        if any(p.is_symlink() or (hasattr(p,'is_junction') and p.is_junction())
               for p in (state_path,state_path.parent,state_path.parent.parent)):
            continue
        state = read_json(state_path)
        target = state.get('species_id')
        if target not in plans or target != state_path.parent.name:
            continue
        if species_id not in (target,'TMS_DMSO') or state.get('identity_sha256') != identity:
            continue
        released = state.get('released',{}).get(species_id,[])
        pool = plans[species_id]['goat_start_ids']
        if isinstance(released,list) and released == pool[:len(released)] and job_id in released:
            return True
    return False


def staged_release_valid(root, row):
    proof_path = root / 'staged_authorization.json'
    if not proof_path.is_file(): return False
    proof = read_json(proof_path)
    original = Path(proof.get('original_root', '')).resolve()
    if original == root.resolve() or proof.get('job_id') != row.get('job_id'): return False
    manifest = original / 'registry/jobs.json'
    if not manifest.is_file() or digest(manifest) != proof.get('manifest_sha256'): return False
    originals = [r for r in read_json(manifest) if r.get('job_id') == row.get('job_id')]
    if len(originals) != 1: return False
    expected = dict(originals[0]); expected['job_directory'] = row['job_directory']
    return (expected == row and controller_released(original, row['job_id'], row.get('species_id')))


def validate_job(root: Path, job_id: str, allow_own_lock: bool = False) -> dict:
    root = Path(root).resolve()
    manifest_path = root / "registry" / "jobs.json"
    rows = read_json(manifest_path)
    if not isinstance(rows, list):
        raise ValueError("registry/jobs.json must be a list")
    matches = [r for r in rows if r.get("job_id") == job_id]
    if len(matches) != 1:
        raise ValueError("Job ID missing or duplicated")
    row = matches[0]
    if (row.get("initial_release") is not True
            and not controller_released(root,job_id,row.get('species_id'))
            and not staged_release_valid(root, row)):
        raise ValueError("Held job: initial_release must be true; no implicit wave release")
    if row.get("status") != "PREPARED_NOT_RUN":
        raise ValueError("Job is not PREPARED_NOT_RUN")
    relative = Path(row.get("job_directory", ""))
    directory = (root / relative).resolve()
    if relative.is_absolute() or directory == root or not directory.is_relative_to(root):
        raise ValueError("Job directory must be confined within this package")
    if not directory.is_dir():
        raise ValueError("Job directory missing")
    for entry in directory.iterdir():
        lower = entry.name.lower()
        forbidden = (lower.startswith("goat.") and lower != "goat.inp") or lower in {
            "receipt.json", "run_receipt.json", ".goat-run.lock"}
        forbidden = forbidden or lower.endswith((".out", ".gbw", ".hess", ".stderr", ".finalensemble.xyz"))
        if lower == ".goat-run.lock" and allow_own_lock:
            forbidden = False
        if forbidden:
            raise ValueError("Existing output/receipt/lock requires review; refusing overwrite: " + entry.name)
    hashes = {}
    for name, field in INPUT_HASH_FIELDS.items():
        path = directory / name
        expected = row.get(field)
        if not isinstance(expected, str) or not re.fullmatch(r"[A-Fa-f0-9]{64}", expected):
            raise ValueError("Missing/invalid manifest hash: " + field)
        if not path.is_file() or not path.resolve().is_relative_to(root) or digest(path) != expected.lower():
            raise ValueError("Input hash mismatch: " + name)
        hashes[name] = expected.lower()
    provenance = directory / "provenance.json"
    if not provenance.is_file():
        raise ValueError("Missing provenance.json")
    hashes["provenance.json"] = digest(provenance)
    if row.get("provenance_sha256") and row["provenance_sha256"].lower() != hashes["provenance.json"]:
        raise ValueError("Input hash mismatch: provenance.json")
    inp = (directory / "goat.inp").read_text(encoding="utf-8")
    references = re.findall(r"(?im)^\s*\*\s*xyzfile\s+[-+]?\d+\s+\d+\s+([^\r\n]+)", inp)
    if references != ["start.xyz"] or len(re.findall(r"(?im)^\s*\*\s*xyz", inp)) != 1:
        raise ValueError("GOAT input must reference exactly the local start.xyz")
    if re.search(r"(?im)^\s*%moinp|\$new_job|%compound|\bNative-", inp):
        raise ValueError("Unsupported external reference, chained job or native-xTB substitution")
    if not re.search(r"(?im)^!.*\bGOAT(?:-ENTROPY)?\b", inp):
        raise ValueError("Input is not a GOAT job")
    ranks = re.findall(r"(?i)\bnprocs\s+(\d+)", inp)
    maxcore = re.findall(r"(?i)%maxcore\s+(\d+)", inp)
    if len(ranks) != 1 or int(ranks[0]) < 1 or len(maxcore) != 1 or int(maxcore[0]) < 1:
        raise ValueError("Input needs one explicit nprocs and maxcore")
    frames = xyz_frames(directory / "start.xyz")
    if len(frames) != 1:
        raise ValueError("start.xyz must contain one frame")
    return {"row": row, "directory": directory, "input_hashes": hashes,
            "manifest_sha256": digest(manifest_path), "symbols": frames[0],
            "nprocs": int(ranks[0]), "maxcore_mb": int(maxcore[0])}


def check_allocation(nprocs: int, maxcore_mb: int) -> None:
    cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if cpus and int(cpus) < nprocs:
        raise ValueError("SLURM CPUs smaller than input nprocs")
    if int(os.environ.get("SLURM_NTASKS", "1")) != 1:
        raise ValueError("Use one serial ORCA driver task, never multiple drivers")
    memory = os.environ.get("WORK2_SLOT_MEMORY_MB") or os.environ.get("SLURM_MEM_PER_NODE")
    if not memory and os.environ.get("SLURM_MEM_PER_CPU"):
        memory = str(int(os.environ["SLURM_MEM_PER_CPU"]) * int(cpus or nprocs))
    if memory and nprocs * maxcore_mb > 0.75 * int(memory):
        raise ValueError("maxcore*nprocs exceeds 75% of SLURM memory")


def stop_child(process: subprocess.Popen) -> None:
    if process.poll() is not None:
        return
    if os.name == "posix":
        os.killpg(process.pid, signal.SIGTERM)
    else:
        # Restrict termination to this runner's still-live child process tree.
        subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGKILL)
        else:
            process.kill()
        process.wait()


def run_one(root: Path, job_id: str, orca: Path, xtb: Path,
            acknowledge_models: bool = False, dry_run: bool = False,
            timeout: float | None = None) -> dict:
    if not acknowledge_models:
        raise ValueError("Supply --acknowledge-models; this does not imply experimental validation")
    if timeout is not None and (not math.isfinite(timeout) or timeout <= 0):
        raise ValueError("Timeout must be finite and positive")
    job = validate_job(root, job_id)
    check_allocation(job["nprocs"], job["maxcore_mb"])
    checks = preflight(orca, xtb, job["nprocs"])
    if not checks["passed"]:
        raise ValueError("Preflight blocked: " + "; ".join(checks["reasons"]))
    command = [checks["orca_executable"], "goat.inp"]
    plan = {"job_id": job_id, "command": command, "job_directory": str(job["directory"]),
            "preflight": checks, "scientific_acceptance": "NOT_EVALUATED",
            "manifest_sha256": job["manifest_sha256"], "input_sha256": job["input_hashes"]}
    if dry_run:
        return {**plan, "status": "DRY_RUN_NO_EXECUTION", "quantum_jobs_launched": 0}
    directory = job["directory"]
    lock = directory / ".goat-run.lock"
    token = uuid.uuid4().hex
    with lock.open("x", encoding="utf-8") as stream:
        from clusterflow.runtime import owner_record
        json.dump(owner_record(token), stream)
    receipt_path = directory / "receipt.json"
    receipt = {**owner_record(token), **plan, "status": "RUNNING", "started_utc": utc(), "quantum_jobs_launched": 0}
    try:
        # Recheck after obtaining exclusive ownership, before making any output.
        fresh = validate_job(root, job_id, allow_own_lock=True)
        if fresh["input_hashes"] != job["input_hashes"] or fresh["manifest_sha256"] != job["manifest_sha256"]:
            raise ValueError("Inputs/manifest changed during preflight")
        for name, field in (("orca_executable", "orca_sha256"), ("xtb_executable", "xtb_sha256")):
            if digest(Path(checks[name])) != checks[field]:
                raise ValueError("Executable changed during preflight: " + name)
        atomic_json(receipt_path, receipt)
        env = os.environ.copy()
        env["PATH"] = str(Path(checks["orca_executable"]).parent) + os.pathsep + env.get("PATH", "")
        env["OMP_NUM_THREADS"] = "1"
        errors: list[str] = []
        process = None
        returncode = None
        try:
            from clusterflow.runtime import execution_environment, termination_handoff
            with execution_environment(orca, directory), termination_handoff(), (directory / "goat.out").open("xb") as stdout, (directory / "goat.stderr").open("xb") as stderr:
                env.update(os.environ)
                env['PATH'] = str(Path(checks['orca_executable']).parent) + os.pathsep + env.get('PATH','')
                env['OMP_NUM_THREADS'] = '1'
                process = subprocess.Popen(command, cwd=directory, env=env, stdout=stdout, stderr=stderr,
                                           shell=False, start_new_session=(os.name == "posix"))
                receipt["quantum_jobs_launched"] = 1
                returncode = process.wait(timeout=timeout)
        except (OSError, subprocess.TimeoutExpired, KeyboardInterrupt) as exc:
            errors.append(type(exc).__name__ + ": " + str(exc))
            if process is not None:
                stop_child(process)
                returncode = process.returncode
        output = (directory / "goat.out").read_text(encoding="utf-8", errors="replace") if (directory / "goat.out").exists() else ""
        if returncode != 0:
            errors.append("ORCA nonzero or absent exit code")
        versions = re.findall(r"Program Version\s+(\d+\.\d+\.\d+)", output)
        if not versions or set(versions) != {REQUIRED_VERSION}:
            errors.append("Actual ORCA output version mismatch or missing")
        if not re.search(r"(?m)^\s*\**\s*ORCA TERMINATED NORMALLY\s*\**\s*$", output):
            errors.append("Normal termination missing")
        try:
            ensemble = xyz_frames(directory / "goat.finalensemble.xyz")
            if any(symbols != job["symbols"] for symbols in ensemble):
                errors.append("Final ensemble atom count/order mismatch")
        except (OSError, ValueError) as exc:
            errors.append("Final ensemble missing/invalid: " + str(exc))
        for name, expected in job["input_hashes"].items():
            if not (directory / name).is_file() or digest(directory / name) != expected:
                errors.append("Input changed during execution: " + name)
        outputs = {p.name: digest(p) for p in directory.iterdir() if p.is_file()
                   and p.name.startswith("goat.") and p.name != "goat.inp"}
        receipt.update({"status": "FAILED" if errors else "COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW",
                        "finished_utc": utc(), "returncode": returncode,
                        "errors": errors, "output_sha256": outputs,
                        "scientific_acceptance": "NOT_EVALUATED_IDENTITY_STEREOCHEMISTRY_MOTIFS_AND_SAMPLING_REQUIRED"})
        atomic_json(receipt_path, receipt)
        return receipt
    except Exception as exc:
        # An unexpected post-launch I/O/parser failure must not leave a false
        # success. Before launch, merely report the exception without outputs.
        if receipt_path.exists():
            receipt.update({"status": "FAILED", "finished_utc": utc(),
                            "errors": [type(exc).__name__ + ": " + str(exc)],
                            "scientific_acceptance": "NOT_EVALUATED"})
            atomic_json(receipt_path, receipt)
        raise
    finally:
        if lock.exists() and read_json(lock).get("token") == token:
            lock.unlink()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("preflight", "run"):
        command = sub.add_parser(name)
        command.add_argument("--orca", type=Path, required=True)
        command.add_argument("--xtb", type=Path, required=True)
        if name == "run":
            command.add_argument("--job-id", required=True)
            command.add_argument("--acknowledge-models", action="store_true")
            command.add_argument("--dry-run", action="store_true")
            command.add_argument("--timeout", type=float)
    args = parser.parse_args(argv)
    try:
        if args.command == "preflight":
            result = preflight(args.orca, args.xtb)
            code = 0 if result["passed"] else 2
        else:
            result = run_one(ROOT, args.job_id, args.orca, args.xtb,
                             args.acknowledge_models, args.dry_run, args.timeout)
            code = 3 if result["status"] == "FAILED" else 0
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return code
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print("STOP: " + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
