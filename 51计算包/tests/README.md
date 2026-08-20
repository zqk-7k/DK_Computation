# Tests

Run the dependency-free control-plane tests with:

```bash
python -m unittest discover -s tests -v
```

They test controller safety, failure classification, profile memory bounds, input/profile agreement and fail-closed Slurm rendering. They do not run ORCA and are not validation of predicted spectra.
