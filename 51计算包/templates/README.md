# Templates

These templates define auditable records; they are not completed results. Copy a template into the relevant compound/stage directory, replace every `null`/placeholder, then hash and freeze it. A production gate must reject records that still contain required provisional fields.

- `condition_species_lock.template.json`: structure, sample condition and requested nuclei.
- `fixed_composition_scenario.template.json`: one microsolvation branch; different compositions get different files.
- `goat_stop_receipt.template.json`: wave-complete stop/escalation evidence.
- `method_branch_receipt.template.json`: method/version/reference and coverage receipt.
- `final_provenance.template.json`: final output traceability.

Never store passwords, SSH keys, API tokens or private cluster credentials in these records.
