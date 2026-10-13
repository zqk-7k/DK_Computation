# ORCA 51-compound NMR PRACTICAL_HIGH_ACCURACY workflow

This folder is a production-oriented, auditable calculation plan for the 51 Figure 3 compounds other than benzamide (`F3_005`). It is based on the benzamide pilot, but it does **not** claim that one pilot validates all chemical classes, nor that any method is guaranteed to reproduce an experimental spectrum exactly.

The default mode is `PRACTICAL_HIGH_ACCURACY`. `ULTRA` is a diagnostic escalation branch for small/light systems or unresolved method disagreement. The package deliberately separates:

- `CLEAN_CALCULATED`: chemical shifts, couplings, assignments and ideal calculated line shapes;
- `SYNTHETIC_APPEARANCE`: optional, explicitly labelled display modelling for line width, solvent peaks, baseline and noise.

Static ORCA calculations do not predict instrument noise, phasing, baseline, receiver response, concentration-dependent intensity, exchange broadening or all experimental line widths from first principles. Synthetic appearance is therefore never presented as experimental data.

## Scientific workflow

1. **Lock identity and conditions.** Confirm structure, stereochemistry, charge/protonation, salt form, isotopic solvent, temperature, concentration, pH/pD where relevant, reference convention and requested nuclei. Unconfirmed fields remain visibly provisional.
2. **Sample fixed-composition scenarios separately.** GOAT explores each microsolvation composition independently. A continuum-only cluster, a cluster with two explicit DMSO molecules, and a cluster with two DMSO plus one water molecule are separate models and are never mixed by a Boltzmann average.
3. **Run adaptive GOAT waves.** The usual `2 + 3 + 5` successful structures are minimum first-pass targets for the three scenario positions, not a convergence guarantee. Extra candidates are kept unreleased and are started in waves only when the coverage gate requires them.
4. **Validate structures before energy work.** Enforce graph identity, atom mapping, stereochemistry, charge, multiplicity, fixed cluster composition and solvent-motif integrity. Deduplicate with cluster-aware geometry and energy rules.
5. **Two-level optimization/frequency treatment.** Use `r2SCAN-3c` as the broad prescreen. Reoptimize the selected high-weight/motif-complete set with `PBE0-D4/def2-TZVPP`, `VeryTightOpt`, `VeryTightSCF` and `DefGrid3`. Small rigid systems and `ULTRA` may receive full PBE0 coverage.
6. **QRRHO and conformer populations.** Apply one frozen QRRHO policy and one temperature per locked condition. Average only conformers with the same molecular composition, charge and multiplicity.
7. **Baseline shieldings.** Calculate the approved conformer set with `PBE0/pcSseg-3`. This is the complete baseline branch.
8. **Independent revDSD validation.** Use relaxed-density `revDSD-PBEP86-D4/2021/pcSseg-3` on representatives covering at least 95% of the PBE0 population in each fixed-composition branch. Do not splice raw PBE0 and revDSD shieldings into one spectrum. `NoFrozenCore` is restricted to a gated small/light `ULTRA` branch; it is not the default for all 51 compounds.
9. **Reference and calibration.** Reference each electronic-structure branch with a same-method, same-basis, same-solvent and same-version reference calculation. Any empirical calibration is fitted on an independent, predeclared calibration set and then frozen; the target spectrum is not used to fit itself.
10. **J and spin simulation.** Keep `PBE/pcJ-3` and `PBE0/pcJ-2` as separate method branches and select them only after the cross-compound J-method gate. Do not average the two J methods. Exchangeable-proton J values may be tabulated, but are excluded from coherent fine structure by default.
11. **Render two output layers.** Generate the clean calculated spectrum first. Only then may a separate synthetic-appearance profile add assumed line widths, solvent/residual-water peaks, baseline and noise.

The machine-readable policy is in `config/`. Detailed rationale and gates are in `docs/`.

### Downstream-builder status

The GOAT inputs, wave controller and HPC runner are execution-ready after a real site
profile is supplied. Downstream jobs cannot be pre-generated before GOAT structures,
reviewed conformer gates and QRRHO evidence exist. The included audited compatibility
builder can construct the direct GOAT -> PBE0 route, but the new default
`r2SCAN-3c -> PBE0` hand-off remains deliberately **FAIL_CLOSED** until its
hash-bound screening gate is commissioned and fixture-tested. This is a visible safety
boundary, not a completed-spectrum claim; see `docs/DOWNSTREAM_BUILDERS.md`.

## Adaptive GOAT stopping rule

A branch may stop at its numerical target only when all of the following also hold:

- important low-energy basins have been rediscovered or independently supported;
- predeclared solute/solvent hydrogen-bond and orientation motifs are covered;
- no important new basin appears in the most recent completed wave;
- the last two completed waves preserve peak order and weighted shifts within the engineering thresholds: nonexchangeable `1H` 0.03–0.05 ppm, `13C` 0.3–0.5 ppm and exchangeable `NH` 0.2 ppm;
- identity, stereochemistry, composition and duplicate gates pass.

Do not cancel running jobs merely because the target count has been reached. Only unreleased/held candidates may be suppressed. A wave is evaluated after every released job reaches a terminal state so that fast jobs are not preferentially selected.

## Special fail-closed routes

- `F3_006` contains iodine. Production property work remains blocked until the X2C/ZORA or ECP choice, basis coverage, auxiliary bases and picture-change treatment have passed a version-specific ORCA 6.1.1 input test.
- F-, P- and Cl-containing compounds require an explicit basis/auxiliary-basis preflight.
- `F3_035`, `F3_043` and `F3_044` contain no carbon; the workflow must not fabricate a `13C` result.
- In D2O, rapidly exchanging analyte hydrogens may be absent or strongly broadened. This is an experimental exchange issue, not a missing static shielding calculation.
- `F3_016` and `F3_017` full-salt and cation-only states are separate species. They are never averaged together.
- The working states for `F3_019`, `F3_020`, `F3_021`, and the identity of `F3_023` require sample-owner confirmation before a result can be labelled production.

See `docs/SPECIAL_CASES_AND_FAIL_CLOSED_GATES.md` for the complete routing rules.

## Server execution order

Before upload, run the package audit and verify the checksum inventory. Upload the entire folder without credentials or private keys. On the cluster, use the site-specific profile under `hpc/` and execute stages in this order:

1. condition/species lock audit;
2. input-bundle checksum and ORCA/OpenMPI/basis preflight;
3. GOAT wave 1, structural validation and coverage decision;
4. additional GOAT waves only when requested by the adaptive controller;
5. r2SCAN-3c optimization/frequency prescreen;
6. PBE0-D4 production reoptimization/frequency work;
7. QRRHO, deduplication and fixed-composition population freeze;
8. PBE0 shielding baseline and same-method reference;
9. gated revDSD validation and same-method reference;
10. separate J branches and method gate;
11. spin simulation, `CLEAN_CALCULATED` export, then optional `SYNTHETIC_APPEARANCE` export;
12. final provenance, logs, hashes and failure receipts.

Use the full path to the serial ORCA driver. Do not start the ORCA driver itself with `mpirun`. Stage inputs to node-local scratch and copy results back atomically. Resource numbers are site profiles, not universal constants.

## Scope of explicit solvent

Two DMSO molecules and one water molecule are modelling choices for local interactions. One explicit water does **not** establish the real water content, activity, exchange rate or aggregation state of a sample. Results from different explicit-solvent stoichiometries are reported as separate sensitivity branches unless an independently justified population model exists.

## Package map

- `registry/`: authoritative Figure 3 structure registry and source artwork.
- `config/`: machine-readable scientific policies and compound/condition tables.
- `01_goat_explicit/`: generated fixed-composition candidate inputs.
- `gates/`: gate receipts and decisions.
- `jobs/`, `hpc/`, `scripts/`: generated jobs, cluster profiles and workflow controls.
- `templates/`: condition, gate, calibration and provenance templates.
- `postprocess/`: output-layer rules and appearance templates.
- `docs/`: workflow rationale, validation limits, routing and literature basis.
- `package_metadata/`: inventory, checksums and audit receipts.

## Status and interpretation

This package is a **calculation protocol**, not 51 completed spectra. All unrun stages, provisional conditions and fail-closed gates must remain visible. Agreement for rigid aromatic nonexchangeable protons may be good after reference/calibration and conformer control, but larger flexible molecules, ions, salts, exchangeable NH/OH and heavy-element compounds have greater uncertainty. Peak heights, baseline and noise are display assumptions unless experimental acquisition and lineshape parameters are supplied.

## Primary documentation and literature

- ORCA 6.1 NMR tutorial: <https://www.faccts.de/docs/orca/6.1/tutorials/spec/NMR.html>
- ORCA 6.1 GOAT manual: <https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/goat.html>
- ORCA 6.1 parallel execution guidance: <https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/parallel.html>
- Schattenberg and Kaupp, high-level shielding benchmark: <https://doi.org/10.1021/acs.jctc.1c00919>
- Liang *et al.*, pcSseg basis-set assessment: <https://doi.org/10.1021/acs.jctc.2c00933>
- Yesiltepe *et al.*, automated conformer-averaged NMR workflows: <https://doi.org/10.1186/s13321-018-0305-8>
- DELTA50 experimental shift benchmark: <https://doi.org/10.3390/molecules28062449>
- Grimme, quasi-rigid-rotor/harmonic-oscillator thermochemistry: <https://doi.org/10.1002/chem.201200497>

The LVK candidate-pair/sky/time analysis from the earlier report is scientifically unrelated to molecular NMR calculation. It is intentionally not used by this package and cannot affect conformer, method or spectrum gates.
