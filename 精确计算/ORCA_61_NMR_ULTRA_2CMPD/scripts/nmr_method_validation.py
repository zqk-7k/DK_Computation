"""Fail-closed validation for ORCA NMR method identity and MP2 density policy."""

from __future__ import annotations

import re
from pathlib import Path


class MethodValidationError(RuntimeError):
    pass


CORRELATED_METHOD_TOKENS = {
    "production_revdsd_pcsseg3": "REVDSDPBEP86D42021",
    "benchmark_ri_mp2_pcsseg3": "RIMP2",
    "benchmark_dlpno_mp2_pcsseg3": "DLPNOMP2",
}


def compact(value: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "", value.upper())


def method_family(method_id: str, fingerprint: str) -> str:
    identifier = method_id.strip().casefold()
    normalized = compact(fingerprint)
    if "vpt2" in identifier or "VPT2" in fingerprint.upper():
        return "vpt2_diagnostic"
    if "revdsd" in identifier or "REVDSDPBEP86D42021" in normalized:
        return "production_revdsd_pcsseg3"
    if "dlpno" in identifier or "DLPNOMP2" in normalized:
        return "benchmark_dlpno_mp2_pcsseg3"
    if "ri_mp2" in identifier or "RIMP2" in normalized:
        return "benchmark_ri_mp2_pcsseg3"
    if "pbe0" in identifier or "PBE0" in normalized:
        return "crosscheck_pbe0_pcsseg3"
    if "autoci" in identifier or "CCSDT" in normalized:
        return "diagnostic_autoci_ccsdt_pcsseg4_gas"
    raise MethodValidationError(
        f"Unsupported or ambiguous NMR method identity: {method_id!r} / {fingerprint!r}"
    )


def validate_nmr_method_evidence(
    input_path: Path,
    output_path: Path,
    method_id: str,
    fingerprint: str,
    *,
    allow_vpt2_diagnostic: bool = False,
    allow_autoci_diagnostic: bool = False,
) -> dict[str, str]:
    """Validate exact method identity and relaxed MP2 density in input and output.

    ORCA prints the parsed input in its output, so requiring the policy markers in
    both files also detects output/input mismatches and incomplete copied results.
    """
    if not input_path.is_file() or not output_path.is_file():
        raise MethodValidationError(f"Missing method evidence: {input_path} / {output_path}")
    input_text = input_path.read_text(encoding="utf-8", errors="replace")
    output_text = output_path.read_text(encoding="utf-8", errors="replace")
    family = method_family(method_id, fingerprint)
    if family == "vpt2_diagnostic" and not allow_vpt2_diagnostic:
        raise MethodValidationError(
            "VPT2 shielding/vibrational corrections are diagnostic-only and cannot enter the main delta scale"
        )
    if family == "diagnostic_autoci_ccsdt_pcsseg4_gas" and not allow_autoci_diagnostic:
        raise MethodValidationError(
            "Gas-phase AutoCI-CCSD(T)/pcSseg-4 is diagnostic-only and cannot enter "
            "the production main-delta or TMS approval path"
        )
    normalized_fingerprint = compact(fingerprint)
    input_compact = compact(input_text)
    output_compact = compact(output_text)
    if "PCSSEG3" not in normalized_fingerprint and family not in {
        "diagnostic_autoci_ccsdt_pcsseg4_gas",
        "vpt2_diagnostic",
    }:
        raise MethodValidationError(f"Method fingerprint lacks pcSseg-3: {fingerprint}")
    if family in CORRELATED_METHOD_TOKENS:
        method_token = CORRELATED_METHOD_TOKENS[family]
        if method_token not in normalized_fingerprint:
            raise MethodValidationError(
                f"Fingerprint does not identify {family}: {fingerprint}"
            )
        if method_token not in input_compact or method_token not in output_compact:
            raise MethodValidationError(
                f"Input/output method identity does not match {family}: {input_path} / {output_path}"
            )
        relaxed_pattern = re.compile(
            r"%mp2\b(?:(?!\bend\b).)*\bdensity\s+relaxed\b",
            flags=re.I | re.S,
        )
        if relaxed_pattern.search(input_text) is None or relaxed_pattern.search(output_text) is None:
            raise MethodValidationError(
                f"{family} NMR requires an explicitly relaxed MP2 density in input and output"
            )
        if (
            "PCSSEG3" not in input_compact
            or "PCSSEG3" not in output_compact
        ):
            raise MethodValidationError(f"{family} input/output lacks pcSseg-3 identity")
        if (
            "NOFROZENCORE" not in normalized_fingerprint
            or "NOFROZENCORE" not in input_compact
            or "NOFROZENCORE" not in output_compact
        ):
            raise MethodValidationError(f"{family} requires NoFrozenCore")
    elif family == "crosscheck_pbe0_pcsseg3":
        if "PBE0" not in input_compact or "PBE0" not in output_compact:
            raise MethodValidationError("PBE0 cross-check method identifier is missing")
    if family != "vpt2_diagnostic":
        if "NMR" not in input_compact:
            raise MethodValidationError("Shielding input lacks the ORCA NMR request")
        # ORCA's NMR request uses GIAO by default; a literal ``GIAO`` input
        # keyword is therefore not required.  The immutable method fingerprint
        # must still name the gauge convention explicitly.
        if "GIAO" not in normalized_fingerprint:
            raise MethodValidationError("Method fingerprint does not lock the GIAO convention")
        if "CHEMICALSHIELDINGSUMMARY" not in output_compact:
            raise MethodValidationError("Output lacks a chemical shielding summary")
    return {
        "method_family": family,
        "method_id": method_id,
        "method_fingerprint": fingerprint,
        "relaxed_mp2_density": (
            "verified_in_input_and_output"
            if family in CORRELATED_METHOD_TOKENS
            else "not_applicable"
        ),
        "vpt2_policy": (
            "diagnostic_only" if family == "vpt2_diagnostic" else "not_vpt2"
        ),
        "autoci_policy": (
            "diagnostic_only"
            if family == "diagnostic_autoci_ccsdt_pcsseg4_gas"
            else "not_autoci"
        ),
    }
