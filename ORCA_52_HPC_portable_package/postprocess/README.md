# Calculated spectra post-processing package

This directory is the portable, Linux-compatible post-processing layer for the
Figure 3 calculation workflow.  It **does not launch ORCA**.  It consumes
completed ORCA shielding/J outputs and keeps the following data layers separate:

1. `P0_calculated_sticks`: ORCA-derived calculated transitions.
2. `P1_clean_calculated_profile`: clean broadened analyte-only model.
3. `P2_synthetic_appearance`: synthetic FID/processing appearance with named
   solvent, water, noise and baseline components.  This is not acquired data.

The HRMS branch uses an explicitly declared ion formula, exact isotope masses
and natural abundances.  ORCA does not predict LC-MS fragmentation, retention
time, ion-source abundance, or an experimental `found` value; none is generated.

## Figure contract

- **Conclusion:** every displayed line can be traced to either a calculated
  analyte transition or a separately named synthetic acquisition component.
- **Evidence chain:** atom shielding/J output → P0 sticks → P1 clean model →
  optional P2 synthetic appearance.
- **Archetype:** compact quantitative spectrum/audit panel.
- **Exports:** editable SVG plus PNG preview; CSV and tab-delimited TXT are the
  source data.
- **Main review risks:** calculated reference scale, conformer weighting,
  exchangeable-proton treatment, approximate first-order fallback for large spin
  systems, and accidental presentation of P2 as experimental.

## Linux environment

From the project root:

```bash
conda create -n orca-post python=3.11 numpy matplotlib rdkit -c conda-forge
conda activate orca-post
python postprocess/production_spectra.py \
  run --config postprocess/examples/benzamide_demo_config.json
```

All paths in the example JSON are relative to the JSON file.  A transferred
package therefore remains runnable if the project directory topology is kept.
For an HPC result directory, pass an explicit Linux path:

```bash
python postprocess/production_spectra.py run \
  --config postprocess/examples/F3_001.json \
  --output-root /scratch/$USER/C_N_results/postprocessed
```

The P2 generator uses empirical DMSO calibration tables shipped in
`real_spectrum_calibration/`. D2O outputs remain explicitly labelled as
assumed when no local D2O calibration exists.

## Required JSON fields

```json
{
  "compound_id": "F3_001",
  "display_name": "compound name",
  "solvent": "DMSO-d6",
  "nuclei": ["1H", "13C"],
  "atom_map": "relative/path/atom_map.csv",
  "shielding_output": "relative/path/shield.property.txt",
  "j_output": "relative/path/j.property.txt",
  "reference_shielding_ppm": {"1H": 31.0, "13C": 186.0},
  "hrms_ions": [
    {"label": "[M+H]+", "ion_formula": "C7H8NO", "charge": 1}
  ]
}
```

For a Boltzmann ensemble, replace the three single-conformer fields with:

```json
"conformers": [
  {
    "conformer_id": "conf_01",
    "weight": 0.70,
    "atom_map": "...",
    "shielding_output": "...",
    "j_output": "..."
  },
  {
    "conformer_id": "conf_02",
    "weight": 0.30,
    "atom_map": "...",
    "shielding_output": "...",
    "j_output": "..."
  }
]
```

Weights are normalized, but the supplied sum is preserved in metadata.  All
conformers must have the same atom order.  `reference_shielding_ppm` must come
from a reference calculation at the declared compatible theory level; this
program does not fit shifts to an experimental spectrum.

## Outputs

Each compound receives atom-level and transition-level CSV/TXT, metadata JSON,
P0/P1 SVG/PNG, optional P2 complex FID and component ledger, plus formula-HRMS
sticks/profile data.  `_aggregate_workbook_ready` contains normalized CSV/TXT
tables that can later be imported into one workbook.

Run the parser/benzamide QA tests with:

```bash
python -m unittest discover -s postprocess/tests -v
```
