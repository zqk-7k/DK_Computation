# Post-processing contract

Post-processing must never change scientific method identity or conceal disagreement. Its first job is to convert approved ORCA outputs into traceable numerical tables. Visual realism is a second, optional layer.

## Required first layer: CLEAN_CALCULATED

Inputs must have passed condition/species, composition, structure, QRRHO, method, reference and coverage gates. Produce:

- per-atom shielding and assignment table;
- per-conformer and population-weighted shifts;
- method/scenario uncertainty table;
- selected J branch and coupling table;
- ideal transition/line list and spectrum CSV;
- a plot labelled `CALCULATED — NOT EXPERIMENTAL`;
- final provenance JSON and hashes.

PBE0 and revDSD are exported as separate spectra. PBE/pcJ-3 and PBE0/pcJ-2 are exported separately until the J-method gate selects one. Different fixed compositions are separate files and panels unless an external population model has been approved.

## Optional second layer: SYNTHETIC_APPEARANCE

Copy a clean calculated spectrum and apply only a frozen profile based on `synthetic_appearance_profile.template.json`. The output filename, plot title and metadata must contain `SYNTHETIC_APPEARANCE` and `NOT EXPERIMENTAL`.

Allowed display operations include assumed Lorentzian/Gaussian/Voigt width, known solvent/HOD/internal-standard peaks, a parameterized baseline and fixed-seed noise. Each parameter requires a source category:

- `USER_SPECIFIED`;
- `INDEPENDENT_EXPERIMENT_LIBRARY`;
- `DISPLAY_ASSUMPTION`.

Do not derive atom-specific shifts, peak heights or widths from the target experimental spectrum and then call the result a prediction. If the target spectrum is used for calibration, the result must be labelled `EXPERIMENT_FITTED_RECONSTRUCTION`, outside the prediction output.

## Exchangeable protons

NH/OH assignments remain in tables. Their line widths and observability require an exchange model or independent experimental parameters. In D2O they may be absent after H/D exchange. Static J values touching exchangeable H are retained but are excluded from coherent fine structure by default.

## Carbon-free compounds

For `F3_035`, `F3_043` and `F3_044`, emit a machine-readable `13C: NOT_APPLICABLE` receipt. Do not render a blank `13C` plot as if it were a calculated spectrum.
