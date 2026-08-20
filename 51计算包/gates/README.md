# Scientific and operational gates

`adaptive_goat_controller.py` does not infer chemical convergence from a count alone. After every released wave reaches a terminal state, write one assessment at:

`gates/goat_assessments/<compound_id>/<scenario_id>.json`

Start from `goat_assessment.template.json` and validate it against `package_metadata/schemas/goat_gate_assessment.schema.json`. `evaluated_job_ids` must include every released candidate, including failed candidates, so a quick-finishing subset cannot bias the stop decision.

The stop gate requires the target number of accepted results, intact identity/composition/stereochemistry, low-energy basin rediscovery, required microsolvation motif coverage, and stable weighted shifts/peak order. Default numerical limits are 0.05 ppm for non-exchangeable 1H, 0.5 ppm for 13C and 0.2 ppm for NH. These are workflow convergence checks, not claims of experimental accuracy.

When the stop gate passes, the controller leaves all unused candidates `HELD`; it never cancels running work. `stop_request.json` and `stop_final_receipt.json` document that behavior.
