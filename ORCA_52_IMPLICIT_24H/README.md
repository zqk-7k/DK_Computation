# ORCA 52-compound implicit-solvent NMR/HRMS workflow

This is a deadline-aware, whole-node workflow for the 52 Figure 3 structures.
It uses implicit solvent only. It produces calculated ^1H/^13C NMR predictions,
full numeric profiles, figures, workbooks, and formula-level HRMS isotope
envelopes. It does **not** produce experimental data.

The 24 h profile is baseline-first: DFT shielding coverage is attempted before
optional higher-level refinements. Proton splittings are a labelled first-order
empirical geometry/topology model, not ab initio ORCA J couplings. One selected
conformer is used per compound; no Boltzmann ensemble is claimed.

Read [README_中文.md](README_中文.md) before running. The production entry point is:

```bash
bash hpc/submit_24h.sh /absolute/path/to/orca
```

ORCA must be invoked by its full path. The workflow never starts the ORCA driver
through an outer `mpirun`.
