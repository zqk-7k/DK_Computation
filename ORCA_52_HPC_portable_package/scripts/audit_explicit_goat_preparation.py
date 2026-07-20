"""Audit the immutable explicit-solvent GOAT preparation bundle."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT / "data" / "GOAT_explicit_job_manifest.csv"
REPORT = PROJECT / "data" / "GOAT_explicit_preparation_audit.json"
ISSUES = PROJECT / "data" / "GOAT_explicit_preparation_issues.csv"
ISOTOPE_POLICY = PROJECT / "data" / "ORCA_6_1_1_isotope_mass_policy.json"
GOAT_XYZ_MASS_POLICY = "PLAIN_H_FOUR_COLUMN_EXTERNAL_XYZ"
DOWNSTREAM_D_MASS_POLICY = "INLINE_XYZ_M_EQUALS_FROM_SELECTED_GOAT_GEOMETRY"
DEUTERIUM_SYNTAX = "H {x_A} {y_A} {z_A} M=2.01410177812"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def xyz(path: Path) -> tuple[list[tuple[str, float, float, float]], list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    issues: list[str] = []
    try:
        count = int(lines[0].strip())
    except Exception:
        return [], ["invalid XYZ atom count"]
    if len(lines) < count + 2:
        return [], ["truncated XYZ"]
    atoms: list[tuple[str, float, float, float]] = []
    for line in lines[2 : count + 2]:
        parts = line.split()
        if len(parts) != 4:
            issues.append("external XYZ atom row is not exactly four fields")
            continue
        try:
            atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
        except ValueError:
            issues.append("non-numeric external XYZ coordinate")
    return atoms, issues


def main() -> int:
    manifest = rows(MANIFEST)
    issues: list[dict[str, str]] = []
    policy: dict[str, object] = {}
    if not ISOTOPE_POLICY.is_file():
        issues.append({"job_id": "PROJECT", "issue": "missing isotope mass policy"})
    else:
        try:
            policy = json.loads(ISOTOPE_POLICY.read_text(encoding="utf-8"))
        except Exception as exc:
            issues.append({"job_id": "PROJECT", "issue": f"invalid isotope policy JSON: {exc}"})
        else:
            goat_policy = policy.get("goat_external_xyz", {})
            freq_policy = policy.get("downstream_opt_freq", {})
            if not isinstance(goat_policy, dict) or goat_policy.get("policy") != GOAT_XYZ_MASS_POLICY:
                issues.append({"job_id": "PROJECT", "issue": "wrong GOAT isotope mass policy"})
            if not isinstance(freq_policy, dict) or freq_policy.get("policy") != DOWNSTREAM_D_MASS_POLICY:
                issues.append({"job_id": "PROJECT", "issue": "wrong downstream isotope mass policy"})
            elif freq_policy.get("deuterium_atom_syntax") != DEUTERIUM_SYNTAX:
                issues.append({"job_id": "PROJECT", "issue": "wrong downstream inline D syntax"})
            docs = policy.get("official_documentation", [])
            urls = {item.get("url", "") for item in docs if isinstance(item, dict)} if isinstance(docs, list) else set()
            if not urls or any(not url.startswith("https://www.faccts.de/docs/orca/6.1/") for url in urls):
                issues.append({"job_id": "PROJECT", "issue": "isotope policy lacks ORCA 6.1 official documentation URLs"})
    scenario_jobs: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    job_ids: set[str] = set()
    xyz_hashes: set[str] = set()
    minimum_distance = math.inf
    for row in manifest:
        job_id = row["job_id"]
        if job_id in job_ids:
            issues.append({"job_id": job_id, "issue": "duplicate job_id"})
        job_ids.add(job_id)
        scenario_jobs[(row["figure_id"], row["scenario_id"])].append(row)
        directory = Path(row["job_directory"])
        required = (
            "goat.inp",
            "start.xyz",
            "atom_map.csv",
            "isotope_map.csv",
            "structure_validation.json",
            "INPUT_BUNDLE_SHA256.txt",
        )
        for name in required:
            if not (directory / name).is_file():
                issues.append({"job_id": job_id, "issue": f"missing {name}"})
        if not all((directory / name).is_file() for name in required):
            continue
        if sha256(directory / "goat.inp") != row["goat_input_sha256"].upper():
            issues.append({"job_id": job_id, "issue": "goat.inp SHA-256 mismatch"})
        calculated_xyz_hash = sha256(directory / "start.xyz")
        if calculated_xyz_hash != row["start_xyz_sha256"].upper():
            issues.append({"job_id": job_id, "issue": "start.xyz SHA-256 mismatch"})
        if calculated_xyz_hash in xyz_hashes:
            issues.append({"job_id": job_id, "issue": "duplicate coordinate hash"})
        xyz_hashes.add(calculated_xyz_hash)
        input_text = (directory / "goat.inp").read_text(encoding="utf-8")
        active = re.sub(r"(?m)#.*$", "", input_text)
        for pattern, label in (
            (r"(?is)%pal\s+.*?nprocs\s+4.*?end", "nprocs 4"),
            (r"(?is)%goat\s+.*?NWorkers\s+4.*?end", "NWorkers 4"),
            (r"(?is)%goat\s+.*?RandomSeed\s+false.*?end", "RandomSeed false"),
            (r"(?im)^\s*%maxcore\s+1200", "%maxcore 1200"),
        ):
            if not re.search(pattern, active):
                issues.append({"job_id": job_id, "issue": f"input lacks {label}"})
        if "<" in active or "NOT_PREPARED" in active:
            issues.append({"job_id": job_id, "issue": "unresolved input placeholder"})
        atom_rows = rows(directory / "atom_map.csv")
        isotope_rows = rows(directory / "isotope_map.csv")
        atoms, xyz_issues = xyz(directory / "start.xyz")
        for problem in xyz_issues:
            issues.append({"job_id": job_id, "issue": problem})
        if len(atom_rows) != len(atoms):
            issues.append({"job_id": job_id, "issue": "atom-map/XYZ count mismatch"})
        isotope_indices = {int(r["atom_index_1based"]): float(r["mass_u"]) for r in isotope_rows}
        for isotope_row in isotope_rows:
            index = int(isotope_row["atom_index_1based"])
            mass = float(isotope_row["mass_u"])
            if index < 1 or index > len(atoms) or not (2.013 < mass < 2.015):
                issues.append({"job_id": job_id, "issue": "invalid isotope-map entry"})
                continue
            if atoms[index - 1][0] != "H" or isotope_row.get("element") != "H":
                issues.append({"job_id": job_id, "issue": "deuterium isotope-map entry is not H"})
            if isotope_row.get("goat_xyz_mass_annotation") != "NONE_PLAIN_H":
                issues.append({"job_id": job_id, "issue": "wrong GOAT isotope annotation policy"})
            if isotope_row.get("downstream_coordinate_source") != "REVIEWED_GOAT_SELECTED_GEOMETRY_ONLY":
                issues.append({"job_id": job_id, "issue": "wrong downstream coordinate source"})
            expected = f'H {{x_A}} {{y_A}} {{z_A}} M={isotope_row["mass_u"]}'
            if isotope_row.get("downstream_inline_mass_syntax") != expected:
                issues.append({"job_id": job_id, "issue": "wrong downstream inline isotope syntax"})
            if isotope_row.get("orca_6_1_manual_verified") != "true":
                issues.append({"job_id": job_id, "issue": "isotope syntax not marked manual-verified"})
        if re.search(r"(?im)^\s*\w+\s+[-+0-9.]\S*\s+[-+0-9.]\S*\s+[-+0-9.]\S*\s+M\s*=", (directory / "start.xyz").read_text(encoding="utf-8")):
            issues.append({"job_id": job_id, "issue": "external GOAT XYZ contains forbidden M= token"})
        validation = json.loads((directory / "structure_validation.json").read_text(encoding="utf-8"))
        mass_policy = validation.get("isotope_mass_policy", {})
        if mass_policy.get("goat_external_xyz") != GOAT_XYZ_MASS_POLICY:
            issues.append({"job_id": job_id, "issue": "validation lacks GOAT mass policy"})
        if mass_policy.get("downstream_opt_freq") != DOWNSTREAM_D_MASS_POLICY:
            issues.append({"job_id": job_id, "issue": "validation lacks downstream D mass policy"})
        if row.get("goat_xyz_mass_policy") != GOAT_XYZ_MASS_POLICY:
            issues.append({"job_id": job_id, "issue": "manifest lacks GOAT mass policy"})
        if row.get("downstream_deuterium_mass_policy") != DOWNSTREAM_D_MASS_POLICY:
            issues.append({"job_id": job_id, "issue": "manifest lacks downstream D mass policy"})
        if ISOTOPE_POLICY.is_file() and row.get("isotope_policy_sha256", "").upper() != sha256(ISOTOPE_POLICY):
            issues.append({"job_id": job_id, "issue": "manifest isotope policy hash mismatch"})
        for i, left in enumerate(atoms):
            for right in atoms[i + 1 :]:
                value = math.dist(left[1:4], right[1:4])
                minimum_distance = min(minimum_distance, value)
                if value < 0.70:
                    issues.append({"job_id": job_id, "issue": f"severe atomic overlap {value:.3f} A"})
                    break

    for (figure_id, scenario_id), group in scenario_jobs.items():
        if len(group) < 10:
            issues.append({"job_id": f"{figure_id}/{scenario_id}", "issue": "fewer than ten starts"})
        if len({r["start_xyz_sha256"] for r in group}) != len(group):
            issues.append({"job_id": f"{figure_id}/{scenario_id}", "issue": "non-unique starts within scenario"})

    temporary_files = [str(path) for path in PROJECT.rglob(".*.tmp")]
    if temporary_files:
        for path in temporary_files:
            issues.append({"job_id": "PROJECT", "issue": f"temporary file remains: {path}"})
    status_counts = Counter(row["primary_solvent"] for row in manifest)
    report = {
        "notice": "CALCULATED MODELS - NOT EXPERIMENTAL DATA",
        "manifest": str(MANIFEST),
        "manifest_sha256": sha256(MANIFEST),
        "jobs": len(manifest),
        "compounds": len({row["figure_id"] for row in manifest}),
        "scenarios": len(scenario_jobs),
        "unique_job_ids": len(job_ids),
        "unique_coordinate_hashes": len(xyz_hashes),
        "solvent_job_counts": dict(status_counts),
        "minimum_interatomic_distance_A": round(minimum_distance, 6),
        "temporary_files": temporary_files,
        "isotope_mass_policy_file": str(ISOTOPE_POLICY),
        "isotope_mass_policy_sha256": sha256(ISOTOPE_POLICY) if ISOTOPE_POLICY.is_file() else "",
        "goat_xyz_mass_policy": GOAT_XYZ_MASS_POLICY,
        "downstream_deuterium_mass_policy": DOWNSTREAM_D_MASS_POLICY,
        "issue_count": len(issues),
        "status": "PASS" if not issues else "FAIL",
    }
    REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with ISSUES.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("job_id", "issue"))
        writer.writeheader()
        writer.writerows(issues)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
