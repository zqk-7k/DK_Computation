# Server start commands

These commands assume one exclusive Slurm node with 104 physical CPU cores,
754 GiB RAM, node-local NVMe, Ubuntu 22.04, and the ORCA 6.1.1 Linux shared
OpenMPI build already installed.

## One-time setup on the login node

```bash
git clone git@github.com:zqk-7k/DK_Computation.git
cd DK_Computation/ORCA_52_IMPLICIT_24H
bash hpc/bootstrap_environment.sh
```

If the cluster does not define `SLURM_TMPDIR`, set the administrator-confirmed
node-local NVMe mount before submission:

```bash
export NVME_SCRATCH_ROOT=/local_nvme/$USER
```

## Submit the 23 h 55 min production allocation

```bash
bash hpc/submit_24h.sh /opt/orca_6_1_1/orca \
  --partition=cpu \
  --account=YOUR_ACCOUNT \
  --qos=normal
```

Remove or change the partition/account/QOS flags to match the site. The first
argument must be the absolute ORCA driver path. Do not prepend `mpirun`.

## Monitor

```bash
squeue -u "$USER"
tail -f logs/slurm-JOBID.out
bash hpc/status.sh results/production_JOBID
```

The run is accepted only when:

- `pipeline_final.json` has `status: COMPLETE`;
- `FINAL_AUDIT.json` has `status: PASS`;
- `FINAL_FILE_MANIFEST_SHA256.csv` exists.

Queue wait is not part of the 24 h node-allocation clock. All NMR and HRMS
outputs are calculated predictions and are labelled as non-experimental.
