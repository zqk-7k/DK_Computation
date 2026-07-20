# Complete calculation workflow and manual gates

This upload package separates highly parallel calculations from scientific
decisions that cannot be automated safely.

## Stage order

1. **Condition and structure lock.** Use the 52-entry registry, solvent table,
   charge, multiplicity, stereochemistry, pD assumptions, counterions and
   explicit-solvent scenarios already stored in `registry/`, `config/` and
   `data/`.
2. **GOAT conformer search.** Submit the 1,120 independent jobs in
   `data/GOAT_explicit_job_manifest.csv`. Each of the 112 fixed scenarios has
   ten independently embedded starting structures.
3. **GOAT ensemble selection gate.** Do not automatically pass all structures
   to DFT. Verify normal termination, connectivity, stereochemistry, cluster
   composition, proton transfer, counterion placement and duplicate
   conformers. Record accepted coordinate SHA-256 values with
   `scripts/draft_goat_selection_gate.py` and the gate template in `jobs/`.
4. **DFT optimization and frequency.** Build jobs only from accepted GOAT
   structures. The locked method is r2SCAN-3c with the declared continuum,
   TightOpt and frequency analysis. Apply D masses during frequency jobs using
   the recorded inline `M=2.01410177812` policy.
5. **Minimum/QRRHO/Boltzmann gate.** Reject wrong connectivity, unintended
   proton transfer, escaped solvent clusters and unresolved imaginary modes.
   Boltzmann-average only conformers with identical composition, charge and
   protonation state. Never mix pD/speciation scenarios.
6. **NMR property jobs.** Calculate GIAO shieldings at PBE0/pcSseg-2 and J
   couplings only after the representative-method gate. F3_006 requires the
   dedicated iodine basis/ECP and relativistic sensitivity gate.
7. **Reference models.** The six models in `reference_jobs/` are supplied as
   `DRAFT_NOT_APPROVED`. They require their own optimization, frequency and
   shielding calculations before calculated shifts can be referenced.
8. **Post-processing.** Keep raw shieldings/J values, clean calculated spectra
   and synthetic acquisition appearance in separate layers. Export atom-level
   tables, transition/stick tables, full XY data, complex FID where requested,
   figures, method metadata and hashes.
9. **HRMS.** Treat formula-based exact-mass/isotope calculations as a separate
   branch. ORCA does not generate an experimental LC-HRMS fragmentation trace.
10. **Final QA.** Run the scientific output audit, retain every input/output
    hash, and label all deliverables as calculated or synthetic.

## Why the two manual gates remain

The author's approval fixes the planned species and solvent assumptions. It
cannot pre-approve calculation-dependent facts such as whether a GOAT conformer
has changed connectivity or whether a DFT stationary point has an imaginary
frequency. These decisions must be made after the corresponding outputs exist.

## Recommended parallelism

Use scheduler arrays across independent jobs, while keeping one ORCA job within
the CPU and memory limits declared in `hpc/config/`. Do not launch all 1,120
jobs on a login node. Start with the smoke-test subset, inspect outputs, and
then release arrays in bounded batches according to site policy.
