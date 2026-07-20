# Figure 3 structure audit for ORCA preparation

## Outcome

- All 52 structures parse and sanitize in RDKit.
- F3_001–F3_052 are retained in numerical order in the delivered SDF.
- Five previously ambiguous entries now have explicit modelling reference states: F3_016, F3_017, F3_019, F3_021 and F3_023.
- These choices remove 2D-structure ambiguity for geometry setup. They do **not** prove the salt, protonation state, enantiomer or tautomer present in an experimental sample.
- ORCA NMR work must still define solvent, dielectric treatment, ion-pair placement, hydration and conformer populations.

## Finalized reference states

1. **F3_016 — L-serinamide hydrochloride**: S alpha centre; explicit `NH3+` and `Cl−`. This follows the hydrochloride identity stated in the author Supplementary Information.
2. **F3_017 — L-leucinamide hydrochloride**: S alpha centre; explicit `NH3+` and `Cl−`. This follows the hydrochloride identity stated in the author Supplementary Information.
3. **F3_019 — L-asparagine zwitterion**: S alpha centre with explicit `NH3+`/`COO−`. This is the common near-neutral aqueous reference state; the original product drawing did not itself prove L/D.
4. **F3_021 — L-argininamide dihydrochloride**: S alpha centre; alpha-ammonium plus guanidinium dication and two explicit chlorides. The dihydrochloride is a common isolated form, but the author source did not document the actual counterion state.
5. **F3_023 — 3-(2-hydroxybenzamido)-L-tyrosine zwitterion**: the secondary-amide connectivity in the source drawing is retained, the alpha centre is defined as S, and the amino acid is shown as a zwitterion. The earlier “nitrotyrosine” wording is incompatible with the drawing because this structure has no nitro group. This identity is a standardized author-selected modelling reference and still requires sample-specific confirmation.

## Other definite corrections retained

1. **F3_006**: 2-iodobenzamide, matching the supplied original Figure 3.
2. **F3_014–F3_018, F3_020 and F3_021**: source-supported alpha stereochemistry restored where available.
3. **F3_017**: leucinamide, not valinamide.
4. **F3_029**: N-methylpyridine-3-carboxamide; not an N-methylpyridinium structure.

## Machine-readability notes

- Use `Figure3_52_structures_corrected.sdf` or the CDX as the identity source for ORCA setup, not optical recognition of the PNG.
- The registry records formula, isomeric SMILES, InChI, InChIKey, charge, stereochemical status and a modelling caveat for every entry.
- F3_022 cholamide retains all encoded stereocentres.
- The earlier open-database search matrix predates several corrections and is not an identity authority; any characterization search should be repeated with the final InChIKeys.
