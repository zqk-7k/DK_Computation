"""Compute-node MPI launcher probe; runs hostname processes, never quantum work."""
import os
from pathlib import Path
import socket
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from clusterflow.common import writej
from clusterflow.runtime import site_config

def main(nprocs=8):
    if not os.environ.get('SLURM_JOB_ID'):raise ValueError('Compute allocation required')
    if nprocs not in (8,16):raise ValueError('Unsupported MPI probe size')
    site=site_config();env={**os.environ,**site['environment']}
    command=['/public/software/apps/apptainer/1.2.4/bin/singularity','exec',site['container'],
             '/bin/sh','-c',f'mpirun --bind-to none -np {nprocs} /bin/hostname']
    result=subprocess.run(command,env=env,capture_output=True,text=True,timeout=120)
    hosts=result.stdout.splitlines();local=socket.gethostname().split('.')[0]
    passed=result.returncode==0 and len(hosts)==nprocs and all(h.split('.')[0]==local for h in hosts)
    record={'status':'PASS' if passed else 'FAIL','hostname':socket.gethostname(),
            'slurm_job_id':os.environ['SLURM_JOB_ID'],'command':command,
            'stdout':result.stdout,'stderr':result.stderr,'returncode':result.returncode,
            'quantum_jobs_launched':0,'nprocs':nprocs,
            'notice':'Local MPI launcher probe only; not ORCA/ensemble acceptance.'}
    path=ROOT/'operations/mpi_probes'/(os.environ['SLURM_JOB_ID']+'_'+socket.gethostname()+'.json')
    writej(path,record)
    print(str(path),flush=True)
    if not passed:raise RuntimeError('MPI must launch the requested processes on the allocated local node: '+str(record))

if __name__=='__main__':main()
