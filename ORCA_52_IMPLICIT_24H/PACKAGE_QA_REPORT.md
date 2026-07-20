# Package QA report

Date: 2026-07-20

Status: **READY FOR LINUX HPC PILOT; NOT YET A 52-COMPOUND PRODUCTION RESULT**

Completed checks:

- 52/52 ordered Figure 3 species present.
- 52/52 solute-only coordinate sets match the approved model SMILES element order.
- No explicit DMSO, water or counterion atom remains in an NMR electronic model.
- F3_016/F3_017 charges are +1; F3_021 is +2; F3_035/F3_044 are -1.
- F3_035/F3_043/F3_044 contain no carbon and are marked ^13C N/A.
- D2O exchange policy leaves no persistent analyte ^1H for F3_035/F3_043/F3_044.
- Offline end-to-end fixture passed: input generation, shift referencing, empirical-J rendering, full-point CSV/TXT, PNG/SVG, Excel and final audit.
- Formula-HRMS engine completed all 52 compounds and 122 candidate-ion envelopes in the local QA run.
- Real ORCA 6.1.1 smoke test passed for F3_001: xTB/ALPB optimization followed by PBE0/pcSseg-1/CPCM H/C shielding; 4/4 requested H/C nuclei parsed.
- Real ORCA 6.1.1 iodine-branch smoke test passed for F3_006: the pcSseg H/C calculation with I `def2-TZVP/def2-ECP` completed; 13/13 requested H/C nuclei parsed.
- Python syntax compilation and the offline QA suite passed after the final code changes.
- Final evidence-chain audit checks raw ORCA output, retained input/geometry hashes, run UUID, workflow fingerprint, reference calculation, peak-area conservation and signal-centroid conservation.
- Exact Python, NumPy, SciPy, Matplotlib, psutil, RDKit and openpyxl versions are recorded and included in the workflow fingerprint; a changed environment cannot silently resume the same run.
- The final SHA-256 manifest is generated only after audit logs and `pipeline_final.json` exist; it excludes only itself.

Production checks that must occur on the target server:

- ORCA Linux build and its required OpenMPI version load correctly.
- Slurm exposes at least 92 physical cores under `--hint=nomultithread`.
- Node-local NVMe has at least 300 GiB free.
- Concurrent ORCA job steps are isolated correctly by the site's Slurm/OpenMPI configuration; run a representative pilot before production.
- A representative pilot confirms actual P90 wall times before relying on the Tier 2 coverage estimate.
- `FINAL_AUDIT.json` is `PASS` after the production run.

The 24 h design enforces deadline-aware stopping and baseline-first scheduling. It cannot guarantee that every molecule reaches Tier 2, or guarantee 52/52 success in the presence of SCF, node or filesystem failures before a real server pilot. Any fallback is retained with its true method label.
