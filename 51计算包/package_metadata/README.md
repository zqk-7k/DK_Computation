# Package metadata

The JSON Schemas define the state, assessment and receipt contracts used by the portable runners. `audit_report.json` and `SHA256SUMS` are generated artifacts; regenerate them after the rest of the package is finalised:

```bash
python scripts/audit_package.py --package-root . --full-hash
python scripts/write_package_checksums.py --package-root .
python scripts/write_package_checksums.py --package-root . --verify
```

The static audit does not run ORCA and must not be described as computational validation.
