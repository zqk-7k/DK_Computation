# HPC execution contract

1. Copy `site_config.example.json` outside the repository, fill in the site-specific absolute ORCA 6.1.1 driver path, modules, and node-local scratch policy, then protect the file from accidental edits.
2. Validate profiles and inputs with `python scripts/audit_package.py --package-root . --site-config /path/site.json --live-site`.
3. Render one Slurm script per job with `scripts/render_slurm.py`. `scheduler_cpus` and `orca_nprocs` are separate by design; the latter must also match `%pal nprocs` in the ORCA input.
4. Submit held GOAT candidates only through `adaptive_goat_controller.py`. It is a dry-run unless `--apply` is given. It never cancels `RUNNING` or `COMPLETING` jobs.
5. Run `sync_scheduler_states.py` before controller evaluation. It is also dry-run by default and records scheduler-side failures that occurred before the ORCA runner could write its own terminal state.
6. The runner invokes ORCA by its absolute path without outer `mpirun`/`srun`, uses node-local scratch, writes a heartbeat, performs hash-verified stage-out into a temporary destination, then atomically publishes an attempt directory and receipt.

The resource profiles are starting points, not promises that every site or molecule will fit. Memory failures may move up the declared resource ladder; SCF, numerical, identity, and imaginary-frequency failures require cause-specific review rather than a blind memory retry.
