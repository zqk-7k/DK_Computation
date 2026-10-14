"""Explicit runtime-version verification on a compute node; starts no QM job."""
import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from clusterflow.common import sha, writej
from clusterflow.runtime import site_config, SITE, check_storage, owner_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', required=True, type=Path)
    parser.add_argument('--execute-version-check', action='store_true', required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('Run inside a current compute-node allocation, after loading site environment')
    if args.record.exists(): raise ValueError('Choose a new runtime evidence path; never overwrite')
    site = site_config()
    outputs = {}
    for name, executable in [('wrapper', site['orca_wrapper']), ('xtb', site['xtb'])]:
        result = subprocess.run([executable, '--version'], capture_output=True, text=True, timeout=180, check=False)
        output = result.stdout + result.stderr
        if result.returncode != 0: raise ValueError(name + ' version command failed: ' + output[-2000:])
        outputs[name] = output
    if not re.search(r'ORCA\s+6\.1\.0(?:\D|$)', outputs['wrapper']):
        raise ValueError('Actual wrapper did not confirm ORCA 6.1.0')
    if not re.search(r'\b6\.7\.\d+\b', outputs['xtb']):
        raise ValueError('Historical xTB 6.7.x route not verified; inspect rather than silently substitute')
    scratch = Path(site['shared_scratch'])
    scratch.mkdir(parents=True, exist_ok=True)
    storage = check_storage(scratch, site['disk']['admit_free_gib'])
    record = {**owner_record('explicit-runtime-version-check'), 'version': '6.1.0',
              'wrapper_sha256': sha(site['orca_wrapper']), 'container_sha256': sha(site['container']),
              'xtb_sha256': sha(site['xtb']), 'site_sha256': sha(os.environ.get('WORK2_SITE_CONFIG', SITE)),
              'version_outputs': outputs, 'storage': storage, 'quantum_jobs_launched': 0,
              'pilot_status': 'NOT_RUN', 'notice': 'Version verification is NOT GOAT/DFT/J validation'}
    writej(args.record, record)
    print(args.record.resolve())


if __name__ == '__main__': main()
