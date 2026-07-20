#!/usr/bin/env python3
"""Shared, dependency-light helpers for the implicit-solvent 24 h workflow."""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import re
import shutil
import signal
import subprocess
import threading
import time
import uuid
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

from rdkit import Chem


ROOT = Path(__file__).resolve().parents[1]
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
_ACTIVE_PROCESS_LOCK = threading.Lock()
_ACTIVE_PROCESSES: dict[int, subprocess.Popen[str]] = {}


def runtime_environment_identity() -> dict[str, str]:
    """Return the numerical/post-processing environment that affects outputs."""

    packages = ("numpy", "scipy", "matplotlib", "psutil", "rdkit", "openpyxl")
    identity = {
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
    }
    for package in packages:
        try:
            identity[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            identity[package] = "NOT_INSTALLED"
    return identity


def utc_stamp() -> str:
    return time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + f".{uuid.uuid4().hex}.tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_banner_csv(path: Path) -> list[dict[str, str]]:
    """Read a CSV whose provenance banner is encoded as leading # lines."""

    data_lines = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        data_lines.extend(line for line in handle if not line.startswith("#"))
    return list(csv.DictReader(data_lines))


def artifact_manifest_issues(
    manifest_path: Path,
    *,
    workflow_root: Path,
    generated_root: Path,
) -> list[str]:
    """Verify every entry in an internal SHA256 manifest.

    Generated outputs can be relocated with their result folder; in that case
    the suffix below the generated-root folder name is resolved under the
    current generated root. Workflow files are resolved under the package root.
    """

    if not manifest_path.is_file():
        return ["manifest_missing"]
    metadata: dict[str, str] = {}
    data_lines: list[str] = []
    with manifest_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for line in handle:
            if line.startswith("#"):
                text = line[1:].strip()
                if ":" in text:
                    key, value = text.split(":", 1)
                    metadata[key.strip()] = value.strip()
            else:
                data_lines.append(line)
    rows = list(csv.DictReader(data_lines))
    if not rows:
        return ["manifest_empty"]
    issues: list[str] = []
    try:
        declared_count = int(metadata.get("file_count", ""))
    except ValueError:
        issues.append("manifest_file_count_missing_or_invalid")
    else:
        if declared_count != len(rows):
            issues.append(f"manifest_file_count_mismatch:{declared_count}!={len(rows)}")
    seen: set[str] = set()
    listed_generated: set[Path] = set()
    for index, row in enumerate(rows, start=1):
        role = row.get("role", "")
        raw = row.get("path", "")
        normalized = Path(raw.replace("\\", "/"))
        candidate = normalized if normalized.is_absolute() else workflow_root / normalized
        if role == "generated_output":
            parts = normalized.parts
            indices = [i for i, part in enumerate(parts) if part == generated_root.name]
            if indices:
                candidate = generated_root.joinpath(*parts[indices[-1] + 1 :])
            listed_generated.add(candidate.resolve())
        key = f"{role}:{raw}"
        if key in seen:
            issues.append(f"manifest_duplicate:{index}:{raw}")
            continue
        seen.add(key)
        if not candidate.is_file():
            issues.append(f"manifest_file_missing:{index}:{raw}")
            continue
        try:
            expected_size = int(row.get("size_bytes", ""))
        except ValueError:
            issues.append(f"manifest_size_invalid:{index}:{raw}")
            continue
        if candidate.stat().st_size != expected_size:
            issues.append(f"manifest_size_mismatch:{index}:{raw}")
        if sha256(candidate) != row.get("sha256", "").upper():
            issues.append(f"manifest_hash_mismatch:{index}:{raw}")
    actual_generated = {
        path.resolve()
        for path in generated_root.rglob("*")
        if path.is_file()
        and path.resolve() != manifest_path.resolve()
        and not path.name.endswith(".tmp")
    }
    missing_entries = sorted(str(path) for path in actual_generated - listed_generated)
    stale_entries = sorted(str(path) for path in listed_generated - actual_generated)
    if missing_entries:
        issues.append(f"manifest_generated_files_unlisted:{len(missing_entries)}")
    if stale_entries:
        issues.append(f"manifest_generated_entries_stale:{len(stale_entries)}")
    return issues


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + f".{uuid.uuid4().hex}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def read_xyz(path: Path) -> tuple[str, list[tuple[str, float, float, float]]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        raise ValueError(f"Malformed XYZ file: {path}")
    count = int(lines[0])
    atoms: list[tuple[str, float, float, float]] = []
    for line in lines[2 : 2 + count]:
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"Malformed XYZ atom line in {path}: {line!r}")
        atoms.append((fields[0], float(fields[1]), float(fields[2]), float(fields[3])))
    if len(atoms) != count:
        raise ValueError(f"XYZ count mismatch in {path}: {len(atoms)} != {count}")
    return lines[1], atoms


def write_xyz(path: Path, atoms: Sequence[tuple[str, float, float, float]], comment: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(len(atoms)), comment]
    lines.extend(f"{element:<2s} {x: .10f} {y: .10f} {z: .10f}" for element, x, y, z in atoms)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def xyz_elements(path: Path) -> list[str]:
    return [atom[0] for atom in read_xyz(path)[1]]


def physical_cpu_ids(cap: int) -> list[int]:
    """Return one schedulable logical CPU ID for each physical core.

    Slurm's ``--hint=nomultithread`` normally exposes one thread per core.  The
    topology check below also behaves correctly during a direct Linux run.
    """

    affinity = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else list(range(os.cpu_count() or 1))
    selected: list[int] = []
    seen: set[tuple[str, str]] = set()
    sys_cpu = Path("/sys/devices/system/cpu")
    if sys_cpu.is_dir():
        for cpu in affinity:
            topology = sys_cpu / f"cpu{cpu}" / "topology"
            try:
                package = (topology / "physical_package_id").read_text().strip()
                core = (topology / "core_id").read_text().strip()
            except OSError:
                package, core = "unknown", str(cpu)
            key = (package, core)
            if key in seen:
                continue
            seen.add(key)
            selected.append(cpu)
    else:
        selected = affinity
    if not selected:
        raise RuntimeError("No CPU is available in the current process affinity mask")
    return selected[: min(cap, len(selected))]


def process_environment(orca_path: Path) -> dict[str, str]:
    env = dict(os.environ)
    orca_dir = str(orca_path.resolve().parent)
    env["PATH"] = orca_dir + os.pathsep + env.get("PATH", "")
    env["LD_LIBRARY_PATH"] = orca_dir + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    for key in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        env[key] = "1"
    env["OMP_DYNAMIC"] = "FALSE"
    return env


def _register_process(process: subprocess.Popen[str]) -> None:
    with _ACTIVE_PROCESS_LOCK:
        _ACTIVE_PROCESSES[process.pid] = process


def _unregister_process(process: subprocess.Popen[str]) -> None:
    with _ACTIVE_PROCESS_LOCK:
        _ACTIVE_PROCESSES.pop(process.pid, None)


def terminate_active_processes(grace_seconds: float = 15.0) -> list[int]:
    """Terminate all currently registered ORCA/srun process groups."""

    with _ACTIVE_PROCESS_LOCK:
        processes = list(_ACTIVE_PROCESSES.values())
    terminated: list[int] = []
    for process in processes:
        if process.poll() is not None:
            continue
        terminated.append(process.pid)
        try:
            if os.name != "nt":
                os.killpg(process.pid, signal.SIGTERM)
            else:
                process.terminate()
        except ProcessLookupError:
            pass
    deadline = time.monotonic() + max(0.0, grace_seconds)
    for process in processes:
        if process.poll() is not None:
            continue
        remaining = max(0.0, deadline - time.monotonic())
        try:
            process.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            try:
                if os.name != "nt":
                    os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
            except ProcessLookupError:
                pass
    return terminated


def normal_termination(output_path: Path) -> bool:
    if not output_path.is_file():
        return False
    tail = output_path.read_text(encoding="utf-8", errors="replace")[-120000:]
    return "ORCA TERMINATED NORMALLY" in tail


def parse_shielding_summary(output_path: Path) -> list[dict[str, Any]]:
    text = output_path.read_text(encoding="utf-8", errors="replace")
    marker = "CHEMICAL SHIELDING SUMMARY (ppm)"
    if marker not in text:
        raise ValueError(f"No chemical shielding summary in {output_path}")
    block = text.rsplit(marker, 1)[1]
    rows: list[dict[str, Any]] = []
    started = False
    for line in block.splitlines():
        match = re.match(
            rf"^\s*(\d+)\s+([A-Z][a-z]?)\s+({FLOAT})\s+({FLOAT})\s*$",
            line,
        )
        if match:
            started = True
            rows.append(
                {
                    "orca_atom_index_0based": int(match.group(1)),
                    "element": match.group(2),
                    "isotropic_shielding_ppm": float(match.group(3)),
                    "anisotropy_ppm": float(match.group(4)),
                }
            )
        elif started and not line.strip():
            break
    if not rows:
        raise ValueError(f"Empty chemical shielding summary in {output_path}")
    return rows


def find_final_geometry(run_dir: Path, basename: str) -> Path | None:
    preferred = [
        run_dir / f"{basename}.xyz",
        run_dir / f"{basename}.globalminimum.xyz",
        run_dir / f"{basename}.finalensemble.xyz",
    ]
    for path in preferred:
        if path.is_file():
            return path
    candidates = sorted(run_dir.glob("*.xyz"), key=lambda path: path.stat().st_mtime, reverse=True)
    return candidates[0] if candidates else None


def geometry_same_composition(reference: Path, candidate: Path) -> bool:
    return xyz_elements(reference) == xyz_elements(candidate)


COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "I": 1.39,
}


def _distance(a: tuple[str, float, float, float], b: tuple[str, float, float, float]) -> float:
    return math.sqrt((a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2)


def _signed_tetra_volume(
    atoms: Sequence[tuple[str, float, float, float]], neighbor_indices0: Sequence[int]
) -> float:
    if len(neighbor_indices0) != 4:
        return 0.0
    points = [atoms[index][1:] for index in neighbor_indices0]
    base = points[3]
    vectors = [[point[k] - base[k] for k in range(3)] for point in points[:3]]
    return (
        vectors[0][0] * (vectors[1][1] * vectors[2][2] - vectors[1][2] * vectors[2][1])
        - vectors[0][1] * (vectors[1][0] * vectors[2][2] - vectors[1][2] * vectors[2][0])
        + vectors[0][2] * (vectors[1][0] * vectors[2][1] - vectors[1][1] * vectors[2][0])
    )


def frozen_geometry_issues(
    reference: Path,
    candidate: Path,
    bonds_path: Path,
    model_smiles: str,
) -> list[str]:
    """Check atom order, covalent topology, proton attachment and chirality.

    This deliberately uses conservative distance cutoffs.  It is a fail-safe,
    not a replacement for inspecting chemically difficult optimized structures.
    """

    _, reference_atoms = read_xyz(reference)
    _, candidate_atoms = read_xyz(candidate)
    if [atom[0] for atom in reference_atoms] != [atom[0] for atom in candidate_atoms]:
        return ["atom_order_or_composition_changed"]
    expected_bonds: dict[tuple[int, int], dict[str, Any]] = {}
    for row in read_csv(bonds_path):
        key = tuple(sorted((int(row["atom_i_1based"]) - 1, int(row["atom_j_1based"]) - 1)))
        expected_bonds[key] = {
            "order": float(row.get("bond_order", "1") or 1),
            "aromatic": str(row.get("is_aromatic", "false")).lower() == "true",
        }
    issues: list[str] = []
    for i, atom_i in enumerate(candidate_atoms):
        for j in range(i + 1, len(candidate_atoms)):
            atom_j = candidate_atoms[j]
            distance = _distance(atom_i, atom_j)
            radii_sum = COVALENT_RADII.get(atom_i[0], 0.8) + COVALENT_RADII.get(atom_j[0], 0.8)
            key = (i, j)
            if key in expected_bonds:
                bond = expected_bonds[key]
                order = float(bond["order"])
                if atom_i[0] == "H" or atom_j[0] == "H":
                    upper = 1.22 * radii_sum + 0.08
                elif bond["aromatic"]:
                    upper = 1.18 * radii_sum + 0.10
                elif order >= 2.5:
                    upper = 1.12 * radii_sum + 0.08
                elif order >= 1.5:
                    upper = 1.15 * radii_sum + 0.08
                else:
                    upper = 1.25 * radii_sum + 0.10
                if distance > upper:
                    issues.append(f"expected_bond_broken:{i + 1}-{j + 1}:{distance:.3f}A")
                elif distance < 0.60 * radii_sum:
                    issues.append(f"expected_bond_too_short:{i + 1}-{j + 1}:{distance:.3f}A")
            elif atom_i[0] != "H" or atom_j[0] != "H":
                # 1.12*(r_cov sum) captures a newly formed X-H or heavy-atom
                # covalent bond without treating ordinary hydrogen bonds as bonds.
                if distance < 1.10 * radii_sum:
                    issues.append(f"possible_new_bond:{i + 1}-{j + 1}:{distance:.3f}A")

    mol = Chem.AddHs(Chem.MolFromSmiles(model_smiles))
    if mol is None or mol.GetNumAtoms() != len(candidate_atoms):
        issues.append("smiles_atom_mapping_unavailable")
        return issues
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    for center_index0, _cip in Chem.FindMolChiralCenters(
        mol, includeUnassigned=False, includeCIP=True
    ):
        neighbors = sorted(neighbor.GetIdx() for neighbor in mol.GetAtomWithIdx(center_index0).GetNeighbors())
        if len(neighbors) != 4:
            continue
        reference_volume = _signed_tetra_volume(reference_atoms, neighbors)
        candidate_volume = _signed_tetra_volume(candidate_atoms, neighbors)
        if abs(reference_volume) > 0.02 and abs(candidate_volume) > 0.02:
            if reference_volume * candidate_volume < 0:
                issues.append(f"tetrahedral_stereochemistry_inverted:{center_index0 + 1}")
        elif abs(candidate_volume) <= 0.02:
            issues.append(f"tetrahedral_center_nearly_planar:{center_index0 + 1}")
    return issues


def orca_header(keywords: str, nprocs: int, maxcore_mb: int, slowconv: bool = False) -> list[str]:
    suffix = " SlowConv" if slowconv and "SlowConv" not in keywords else ""
    return [
        f"! {keywords}{suffix}",
        "",
        "%pal",
        f"  nprocs {nprocs}",
        "end",
        "",
        f"%maxcore {maxcore_mb}",
        "",
        "%scf",
        "  MaxIter 500",
        "end",
        "",
    ]


def build_opt_input(
    keywords: str,
    charge: int,
    multiplicity: int,
    xyz_name: str,
    nprocs: int,
    maxcore_mb: int,
    goat: dict[str, Any] | None = None,
) -> str:
    lines = orca_header(keywords, nprocs, maxcore_mb)
    if goat is not None:
        lines.extend(
            [
                "%goat",
                f"  NWorkers {goat['nworkers']}",
                f"  MinGlobalIter {goat['min_global_iter']}",
                f"  MaxGlobalIter {goat['max_global_iter']}",
                f"  MaxEn {goat['max_en_kcal_mol']}",
                f"  RandomSeed {str(goat['random_seed']).lower()}",
                "  AutoWall true",
                "  Align true",
                "  FreezeCisTrans true",
                "  KeepWorkerData false",
                "end",
                "",
            ]
        )
    lines.append(f"* xyzfile {charge} {multiplicity} {xyz_name}")
    return "\n".join(lines) + "\n"


def build_shield_input(
    keywords: str,
    charge: int,
    multiplicity: int,
    xyz_name: str,
    nprocs: int,
    maxcore_mb: int,
    elements: Sequence[str],
    iodine_override: Sequence[str] | None = None,
    slowconv: bool = False,
) -> str:
    lines = orca_header(keywords, nprocs, maxcore_mb, slowconv=slowconv)
    if iodine_override:
        lines.extend(list(iodine_override) + [""])
    lines.append(f"* xyzfile {charge} {multiplicity} {xyz_name}")
    lines.extend(["", "%eprnmr", "  Ori = GIAO"])
    if "H" in elements:
        lines.append("  Nuclei = ALL H { SHIFT }")
    if "C" in elements:
        lines.append("  Nuclei = ALL C { SHIFT }")
    lines.extend(["end", ""])
    return "\n".join(lines)


@dataclass
class RunRecord:
    job_id: str
    attempt: int
    status: str
    start_utc: str
    end_utc: str
    elapsed_seconds: float
    return_code: int | None
    timed_out: bool
    normal_termination: bool
    input_sha256: str
    geometry_sha256: str
    cpu_ids: list[int]
    launch_mode: str
    scratch_directory: str
    result_directory: str
    output_file: str
    error: str


def _next_attempt(job_root: Path) -> int:
    attempts = []
    for path in job_root.glob("attempt_*"):
        match = re.fullmatch(r"attempt_(\d{3})", path.name)
        if match:
            attempts.append(int(match.group(1)))
    return max(attempts, default=0) + 1


def run_orca_job(
    *,
    job_id: str,
    input_text: str,
    geometry_path: Path,
    orca_path: Path,
    scratch_root: Path,
    raw_root: Path,
    cpu_ids: Sequence[int],
    timeout_seconds: float,
    basename: str,
) -> tuple[RunRecord, Path]:
    """Run ORCA directly, pin it to physical CPUs, and atomically collect files."""

    if not orca_path.is_file():
        raise FileNotFoundError(f"ORCA executable not found: {orca_path}")
    job_root = raw_root / job_id
    job_root.mkdir(parents=True, exist_ok=True)
    attempt = _next_attempt(job_root)
    result_dir = job_root / f"attempt_{attempt:03d}"
    scratch = scratch_root / f"{job_id}_{attempt:03d}_{uuid.uuid4().hex[:10]}"
    scratch.mkdir(parents=True, exist_ok=False)
    input_path = scratch / f"{basename}.inp"
    geometry_copy = scratch / geometry_path.name
    input_path.write_text(input_text, encoding="utf-8")
    shutil.copy2(geometry_path, geometry_copy)
    output_path = scratch / f"{basename}.out"
    start = time.monotonic()
    start_utc = utc_stamp()
    return_code: int | None = None
    timed_out = False
    error = ""
    command = [str(orca_path.resolve()), input_path.name]
    launch_mode = "DIRECT"
    if (
        os.name != "nt"
        and os.environ.get("ORCA_USE_SLURM_STEPS") == "1"
        and os.environ.get("SLURM_JOB_ID")
        and shutil.which("srun")
    ):
        command = [
            "srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1",
            f"--cpus-per-task={len(cpu_ids)}", "--cpu-bind=cores",
            "--kill-on-bad-exit=1", *command,
        ]
        launch_mode = "SLURM_EXCLUSIVE_JOB_STEP"
    elif os.name != "nt" and shutil.which("taskset"):
        command = ["taskset", "-c", ",".join(str(cpu) for cpu in cpu_ids), *command]
        launch_mode = "TASKSET_DIRECT"
    process: subprocess.Popen[str] | None = None
    try:
        with output_path.open("w", encoding="utf-8", errors="replace") as output_handle:
            process = subprocess.Popen(
                command,
                cwd=scratch,
                stdout=output_handle,
                stderr=subprocess.STDOUT,
                text=True,
                env=process_environment(orca_path),
                start_new_session=(os.name != "nt"),
            )
            _register_process(process)
            try:
                return_code = process.wait(timeout=max(1.0, timeout_seconds))
            except subprocess.TimeoutExpired:
                timed_out = True
                if os.name != "nt":
                    os.killpg(process.pid, signal.SIGTERM)
                else:
                    process.terminate()
                try:
                    process.wait(timeout=60)
                except subprocess.TimeoutExpired:
                    if os.name != "nt":
                        os.killpg(process.pid, signal.SIGKILL)
                    else:
                        process.kill()
                    process.wait(timeout=30)
                return_code = process.returncode
            finally:
                _unregister_process(process)
    except Exception as exc:  # collection still proceeds
        error = f"{type(exc).__name__}: {exc}"
        if process is not None and process.poll() is None:
            try:
                process.kill()
            except Exception:
                pass

    result_tmp = job_root / f".attempt_{attempt:03d}_{uuid.uuid4().hex}.tmp"
    result_tmp.mkdir(parents=True, exist_ok=False)
    keep_suffixes = {
        ".inp", ".out", ".xyz", ".property.txt", ".gbw", ".hess", ".engrad",
        ".trj", ".opt", ".loc", ".cis", ".json", ".txt",
    }
    for path in scratch.iterdir():
        if not path.is_file():
            continue
        if path.suffix.lower() in keep_suffixes or "ensemble" in path.name or "globalminimum" in path.name:
            try:
                shutil.copy2(path, result_tmp / path.name)
            except OSError as exc:
                error = (error + "; " if error else "") + f"copy {path.name}: {exc}"
    os.replace(result_tmp, result_dir)
    collected_output = result_dir / output_path.name
    normal = normal_termination(collected_output)
    status = "COMPLETE" if return_code == 0 and normal and not timed_out else "TIMEOUT" if timed_out else "FAILED"
    elapsed = time.monotonic() - start
    record = RunRecord(
        job_id=job_id,
        attempt=attempt,
        status=status,
        start_utc=start_utc,
        end_utc=utc_stamp(),
        elapsed_seconds=elapsed,
        return_code=return_code,
        timed_out=timed_out,
        normal_termination=normal,
        input_sha256=hashlib.sha256(input_text.encode("utf-8")).hexdigest().upper(),
        geometry_sha256=sha256(geometry_path),
        cpu_ids=list(cpu_ids),
        launch_mode=launch_mode,
        scratch_directory=str(scratch),
        result_directory=str(result_dir),
        output_file=str(collected_output),
        error=error,
    )
    write_json(result_dir / "run_record.json", asdict(record))
    try:
        shutil.rmtree(scratch)
    except OSError:
        pass
    return record, result_dir


class CorePool:
    """A deterministic physical-core token pool for concurrent ORCA jobs."""

    def __init__(self, cpu_ids: Sequence[int]):
        self._free = sorted(cpu_ids)

    @property
    def available(self) -> int:
        return len(self._free)

    @property
    def total(self) -> int:
        return len(self._free)

    def acquire(self, count: int) -> list[int] | None:
        if len(self._free) < count:
            return None
        allocated = self._free[:count]
        del self._free[:count]
        return allocated

    def release(self, cpu_ids: Sequence[int]) -> None:
        overlap = set(cpu_ids) & set(self._free)
        if overlap:
            raise RuntimeError(f"CPU token double release: {sorted(overlap)}")
        self._free.extend(cpu_ids)
        self._free.sort()


def boltzmann_weight(energy_kcal: float, temperature_k: float = 298.15) -> float:
    return math.exp(-energy_kcal / (0.00198720425864083 * temperature_k))
