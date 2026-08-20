# Downstream builders: current scope and fail-closed boundary

These scripts were copied without scientific-method edits from the previously audited
`ORCA_61_NMR_ULTRA_2CMPD` package so that their provenance and gate semantics remain
visible. They are **builders only**: no script in this set launches ORCA.

## Included scripts

- `workflow_schema.py`: hashes, path confinement, XYZ/map and completed-state checks.
- `draft_goat_selection_gate.py`: creates a review draft from completed GOAT evidence.
- `build_dft_jobs.py`: consumes a fully reviewed GOAT gate and builds the **direct
  PBE0-D4/def2-TZVPP opt/freq route** used by the source package.
- `draft_minimum_boltzmann_gate.py`: creates a QRRHO/weight review draft from completed
  opt/freq evidence.
- `build_property_jobs.py`: builds separate PBE0/pcSseg-3 and revDSD/pcSseg-3 shielding
  branches plus gated J jobs; raw method branches are never averaged.
- `build_j_benchmark_jobs.py`: an explicit fail-closed placeholder. The source
  two-compound benchmark hard-coded an excluded compound and is not reused.

## Important: r2SCAN two-layer route is not silently implemented here

The new 51-compound scientific policy specifies r2SCAN-3c broad prescreening followed
by selected PBE0-D4 production opt/freq. The copied `build_dft_jobs.py` goes directly
from an approved GOAT selection to PBE0-D4; it does **not** implement the intermediate
r2SCAN/QRRHO selection. Therefore:

1. the direct builder may be used only as a documented compatibility/ULTRA route;
2. it must not be described as the default two-layer workflow;
3. the default 51-compound production route remains **FAIL_CLOSED** until a separate
   r2SCAN builder and a hash-bound r2SCAN-to-PBE0 gate are commissioned and tested.

## Gate and configuration compatibility

The copied scripts intentionally retain their original strict schemas. They expect
source-package artefacts such as completed runner states, approved selection gates,
`config/nucleus_conditions.csv`, a compatible DFT resource-policy schema, workflow
profile bindings, approved same-method reference outputs, a J-method gate, and the
F3_006 iodine/relativistic receipt. The new package does not fabricate these approvals.
If a required artefact or hash is absent, the builder exits rather than guessing.

Before production use, adapt these scripts in a reviewed change set to the 51-compound
manifest fields and `hpc/resource_profiles.json`; then add fixture-based tests for the
r2SCAN -> PBE0 hand-off, D2O isotope maps, fixed-composition isolation, F3_006 and
F3_023 fail-closed gates, and carbon-free F3_035/F3_043/F3_044 outputs.

## Syntax checks

From the package root:

```bash
python -m py_compile scripts/workflow_schema.py scripts/draft_goat_selection_gate.py \
  scripts/draft_minimum_boltzmann_gate.py scripts/build_dft_jobs.py \
  scripts/build_property_jobs.py scripts/build_j_benchmark_jobs.py
python scripts/build_dft_jobs.py --help
python scripts/build_property_jobs.py --help
```

`--help` is safe: it imports and validates command-line definitions only; it neither
creates jobs nor starts calculations.
