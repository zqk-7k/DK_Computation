# Working rules

- Preserve source structures and completed outputs. Use new directories for changed settings.
- This package starts with a finite crystallographic PW12O40 reference, not a supported single-atom model.
- Never periodically replicate the reference CIF storage box. Use the finite XYZ coordinates.
- The seven-atom W1O6 fragment is a geometry diagnostic. Do not optimize it as a stable molecule or publish its quantitative XANES as the catalyst.
- Run FEFF/QE on compute nodes using the real scheduler. Reuse authorized existing installations; do not install global software or request credentials.
- Do not impose spectral shifts, FT peak positions, oxidation states, active-site conclusions, or potential labels.
- Keep actual xmu.dat, chi.dat, FEFF inputs, logs, checksums and all failed/converged attempts. Do not pad missing energy ranges.
- Use the same local-relaxation mask across paired catalyst models. Neighboring Co/Fe must not be frozen to manufacture unchanged spectra.
- Calculation completion requires actual native outputs; validation fixtures are not spectra and must not enter the delivered results directory.
