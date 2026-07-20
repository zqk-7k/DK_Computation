# ORCA Figure 3 HPC portable package

This folder is the portable calculation handoff for 52 Figure 3 compounds. It
contains all 1,120 prepared GOAT input bundles, the locked structure/condition
registries, downstream job builders, reference-model templates, NMR/HRMS
post-processing code, tests, scheduler templates and integrity metadata.

**Scientific boundary:** every result from this package is calculated. A final
spectrum with added solvent, water, noise, line width or baseline remains a
synthetic acquisition appearance and must not be described as measured data.

## Package contents

- `01_goat_explicit/` — 1,120 immutable six-file GOAT input bundles.
- `data/GOAT_explicit_job_manifest.csv` — canonical package-relative manifest.
- `registry/` — 52 explicit structures, InChIKey/SMILES registry and visual
  reference.
- `config/` — solvent/speciation choices, calculation methods and acquisition
  appearance policy.
- `jobs/` — staged approval, GOAT selection and minimum/Boltzmann gate schemas.
- `reference_jobs/` — six calculated reference-model templates; still draft.
- `scripts/` — structure preparation, gate drafting, DFT/property builders,
  runners and output processors.
- `postprocess/` — production NMR/HRMS parsing, tabulation and plotting layer.
- `hpc/` — site configuration, Slurm/PBS templates, submission/status helpers,
  examples and scheduler tests.
- `tests/` — small source fixtures and Python contract tests only.
- `package_metadata/` — portability audit, inventory and SHA-256 manifest.

ORCA, xTB and MPI are deliberately not included.

## Before upload

From the package root, verify the archive:

```bash
python3 package_metadata/audit_portable_package.py
```

The required result is `status: PASS`, 1,120 GOAT jobs, 52 compounds, 112
fixed scenarios, zero Windows absolute runtime paths and zero calculation
outputs/caches/binaries.

After transfer and before editing the site configuration, verify every shipped
file with:

```bash
sha256sum -c package_metadata/PACKAGE_SHA256SUMS.txt
```

## On the supercomputer

1. Upload the whole folder without flattening its directory structure.
2. Load the site ORCA, xTB, Python and scheduler modules. ORCA 6.1.1 is the
   locked target version; record any version substitution before calculation.
3. Copy the example site configuration under `hpc/config/`, then set the actual
   account, partition/queue, wall time, memory, cores, scratch path and ORCA/xTB
   commands. Never hard-code credentials in the package.
4. Run the scheduler smoke tests in `hpc/tests/` and submit the documented small
   pilot subset first.
5. Submit GOAT arrays using the canonical manifest:
   `data/GOAT_explicit_job_manifest.csv`.
6. After GOAT completion, download or inspect the ensemble files and complete
   the GOAT ensemble-selection gate before building DFT jobs.
7. Continue in the exact stage order documented in
   `docs/WORKFLOW_AND_MANUAL_GATES.md`.

Scheduler-specific commands and site profiles are documented under
`hpc/README.md`. The scheduler wrapper must resolve every manifest path from the
package root; it must not depend on the original Windows directory.

On Linux, use `hpc/bin/` and `hpc/submit/` as the production execution layer.
The `scripts/run_*_serial.py` and older local runners are retained only for
workflow provenance and local compatibility; they do not replace the tested
SLURM/PBS task runner.

## Integrity and restart

Each GOAT bundle has its own `INPUT_BUNDLE_SHA256.txt`. The manifest also stores
the expected input, coordinate and bundle hashes. The complete upload checksum
list is `package_metadata/PACKAGE_SHA256SUMS.txt`.

Do not edit an input in place after submission. If a scientifically justified
change is required, create a new scenario/job identifier, regenerate hashes and
record the change. Scheduler outputs, logs and scratch files should be written
outside the immutable input tree.

## Output expectations

The final result tree should keep these layers separate:

- raw ORCA/xTB outputs and restart files;
- selected conformers plus manual gate decisions;
- optimized minima, frequencies, QRRHO values and Boltzmann tables;
- raw shielding and J-coupling tables;
- clean calculated NMR spectra;
- separately labelled synthetic acquisition appearance;
- formula-level exact-mass and isotope-envelope predictions;
- TXT/CSV data, figures, provenance JSON and SHA-256 records.

See `NOTICE_CALCULATED_ONLY.txt` before using any output in a manuscript or
Supplementary Information.
