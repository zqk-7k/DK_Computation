#!/usr/bin/env python3
"""Resume-aware Linux HPC coordinator. No FEFF work runs on a login node implicitly."""
from __future__ import annotations
import argparse
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time

SCRIPT_ROOT=Path(__file__).resolve().parent.parent


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write_json(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,ensure_ascii=False),encoding='utf-8')
    temporary.replace(path)


def stamp(): return datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@contextmanager
def coordinator_lock(root):
    lock_path=Path(root)/'.workflow_coordinator.lock'
    with lock_path.open('a+b') as handle:
        handle.seek(0)
        if handle.read(1)==b'':handle.write(b'0');handle.flush()
        handle.seek(0)
        try:
            if os.name=='posix':
                import fcntl
                fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
            else:
                import msvcrt
                msvcrt.locking(handle.fileno(),msvcrt.LK_NBLCK,1)
        except (OSError,BlockingIOError) as exc:
            raise RuntimeError('Another workflow coordinator holds this project lock. Resume that process rather than submitting duplicates.') from exc
        try:yield
        finally:
            handle.seek(0)
            if os.name=='posix':fcntl.flock(handle,fcntl.LOCK_UN)
            else:msvcrt.locking(handle.fileno(),msvcrt.LK_UNLCK,1)


def checked_run(command, **kwargs):
    print('+',shlex.join(map(str,command)),flush=True)
    return subprocess.run(list(map(str,command)),check=True,**kwargs)


def scheduler_kind(requested, config):
    requested=requested if requested!='auto' else config.get('scheduler',{}).get('type','auto')
    if requested!='auto': return requested
    found=[name for name,exe in [('slurm','sbatch'),('pbs','qsub')] if shutil.which(exe)]
    if len(found)==1:return found[0]
    if len(found)>1:raise RuntimeError('Both sbatch and qsub exist. Set scheduler.type to the actual scheduler.')
    raise RuntimeError('No supported scheduler found. Load its module or use --scheduler local --allow-local INSIDE an allocated compute session. No implicit login-node execution.')


def job_dir(job, root):
    path=Path(job['path'])
    path=path if path.is_absolute() else root/path
    path=path.resolve()
    if not path.is_relative_to(root.resolve()):raise ValueError('Job escaped workflow root')
    return path


def job_status(job,root):
    path=job_dir(job,root);r=path/'run_result.json'
    if not r.exists():return 'pending'
    try:data=read_json(r)
    except (ValueError,OSError):return 'pending'
    if not (path/'feff.inp').exists() or data.get('input_sha256')!=sha(path/'feff.inp'):return 'pending'
    status=data.get('status','pending')
    if status=='complete':
        files=(data.get('validation') or {}).get('files',[])
        if not files:return 'pending'
        for record in files:
            target=path/record['name']
            if not target.is_file() or sha(target)!=record.get('sha256'):return 'pending'
        if job.get('mode')=='exafs' and not (path/'paths.dat').is_file():return 'pending'
    return status


def manifest_path(root, config, profile):
    return root/config.get('paths',{}).get('manifests_dir','manifests')/f'{profile}_jobs.json'


def prepare(config_path,profile):
    checked_run([sys.executable,SCRIPT_ROOT/'scripts'/'prepare_inputs.py','--config',config_path,'--profile',profile])


def _directive(value,name):
    value=str(value)
    if '\n' in value or '\r' in value:raise ValueError(f'Newline in {name}')
    return value


def render_array(root,config_path,config,profile,tasks,scheduler,run_id):
    directory=root/'submissions'/profile/run_id;directory.mkdir(parents=True,exist_ok=True)
    task_file=directory/'tasks.json';write_json(task_file,tasks)
    settings=config.get('scheduler',{})
    cpus=int(settings.get('cpus_per_task',1));mem=int(settings.get('memory_gb',8))
    maximum=int(settings.get('max_concurrent',4));wall=_directive(settings.get('walltime','04:00:00'),'walltime')
    if min(cpus,mem,maximum)<1:raise ValueError('CPU, memory and concurrency must be positive')
    log_dir=root/'logs'/profile;log_dir.mkdir(parents=True,exist_ok=True)
    lines=['#!/usr/bin/env bash']
    if scheduler=='slurm':
        lines += [f'#SBATCH --job-name=xas_{profile}',f'#SBATCH --array=0-{len(tasks)-1}%{maximum}',
                  '#SBATCH --nodes=1','#SBATCH --ntasks=1',f'#SBATCH --cpus-per-task={cpus}',
                  f'#SBATCH --mem={mem}G',f'#SBATCH --time={wall}',
                  f'#SBATCH --output={shlex.quote(str(log_dir/"%A_%a.out"))}',
                  f'#SBATCH --error={shlex.quote(str(log_dir/"%A_%a.err"))}']
        for field,flag in [('partition','partition'),('account','account')]:
            if settings.get(field):lines.append(f'#SBATCH --{flag}={_directive(settings[field],field)}')
    elif scheduler=='pbs':
        flavor=settings.get('pbs_flavor','openpbs')
        lines += [f'#PBS -N xas_{profile}','#PBS -S /bin/bash',f'#PBS -l walltime={wall}', '#PBS -j oe',
                  f'#PBS -o {shlex.quote(str(log_dir))}']
        if flavor=='torque':
            if len(tasks)>1:lines.append(f'#PBS -t 0-{len(tasks)-1}%{maximum}')
            lines.append(f'#PBS -l nodes=1:ppn={cpus},mem={mem}gb')
        else:
            if len(tasks)>1:lines.append(f'#PBS -J 0-{len(tasks)-1}%{maximum}')
            lines.append(f'#PBS -l select=1:ncpus={cpus}:mem={mem}gb')
        if settings.get('queue'):lines.append(f'#PBS -q {_directive(settings["queue"],"queue")}')
        if settings.get('account'):lines.append(f'#PBS -A {_directive(settings["account"],"account")}')
    else:raise ValueError(scheduler)
    lines += ['set -euo pipefail',f'cd {shlex.quote(str(root))}',
              'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1']
    lines += [str(command) for command in settings.get('module_commands',[])]
    command=[sys.executable,str(SCRIPT_ROOT/'scripts'/'run_array_task.py'),
             '--task-list',str(task_file),'--config',str(config_path)]
    if scheduler=='pbs' and len(tasks)==1:command += ['--index','0']
    lines.append(shlex.join(command))
    script=directory/'array.sh';script.write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
    return script,task_file


def scheduler_active(kind,job_id):
    if kind=='slurm':
        result=subprocess.run(['squeue','-h','-j',job_id,'-o','%T'],capture_output=True,text=True)
        if result.returncode:
            if 'invalid job id' in result.stderr.lower():return False
            return None
        return bool(result.stdout.strip())
    result=subprocess.run(['qstat','-f',job_id],capture_output=True,text=True)
    if result.returncode==0:
        match=re.search(r'job_state\s*=\s*([A-Za-z])',result.stdout)
        if match:return match.group(1).upper() not in {'F','C'}
        return None
    if 'unknown job' in (result.stdout+result.stderr).lower():return False
    return None


def run_profile(config_path,profile,requested='auto',allow_local=False,dry_run=False):
    config_path=Path(config_path).resolve();root=config_path.parent;config=read_json(config_path)
    kind=scheduler_kind(requested,config)
    state_file=root/'submissions'/profile/'current_submission.json'
    previous=read_json(state_file) if state_file.exists() else {}
    contract={'feff':config.get('feff',{}),'profile':config.get('profiles',{}).get(profile),
              'software':config.get('software',{}),'paths':config.get('paths',{}),
              'structures':[(r['file'],sha(root/r['file'])) for r in read_json(root/'structure_checksums.json')]}
    contract_sha=hashlib.sha256(json.dumps(contract,sort_keys=True).encode()).hexdigest()
    active=None
    if kind!='local' and previous.get('job_id') and previous.get('scheduler')==kind:
        active=scheduler_active(kind,previous['job_id'])
        if active is None:
            raise RuntimeError('Cannot determine previous scheduler job state. Restore scheduler query access before resubmitting; duplicate work is not submitted.')
        if active is True and previous.get('contract_sha256')!=contract_sha:
            raise RuntimeError('A running array uses different inputs/settings. Inspect it before modifying this profile; no inputs were rewritten.')
    prepare(config_path,profile)
    jobs=read_json(manifest_path(root,config,profile))
    if isinstance(jobs,dict):jobs=jobs['jobs']
    pending=[j for j in jobs if job_status(j,root)!='complete']
    if not pending:
        print(f'{profile}: all {len(jobs)} jobs already complete');return
    if kind=='local':
        if not allow_local:raise RuntimeError('Local FEFF execution requires --allow-local in an allocated compute session.')
        if dry_run:
            print(f'DRY RUN: {len(pending)} local tasks');return
        for job in pending:job['feff_input_sha256']=sha(job_dir(job,root)/'feff.inp')
        task_file=root/'submissions'/profile/(stamp()+'_local_tasks.json');write_json(task_file,pending)
        def task(i):
            return subprocess.run([sys.executable,str(SCRIPT_ROOT/'scripts'/'run_array_task.py'),
                                   '--task-list',str(task_file),'--config',str(config_path),'--index',str(i)]).returncode
        with ThreadPoolExecutor(max_workers=int(config.get('scheduler',{}).get('max_concurrent',4))) as pool:
            codes=list(pool.map(task,range(len(pending))))
        if any(codes):raise RuntimeError('Local task failure. Inspect job run_result.json and logs, then resume.')
        return
    signature=hashlib.sha256((''.join(j['id']+sha(job_dir(j,root)/'feff.inp') for j in jobs)).encode()).hexdigest()
    if active is True and previous.get('signature')!=signature:
        raise RuntimeError('Existing array is active with different inputs. Do not overwrite running jobs; inspect/cancel it before changing configuration.')
    if active is True and previous.get('signature')==signature:
        job_id=previous['job_id'];print('Resuming monitor of',job_id,flush=True)
        if dry_run:return
    else:
        run_id=stamp()
        for job in pending:job['feff_input_sha256']=sha(job_dir(job,root)/'feff.inp')
        snapshot=root/f'config_snapshot_{profile}_{run_id}.json';write_json(snapshot,config)
        script,task_file=render_array(root,snapshot,config,profile,pending,kind,run_id)
        if dry_run:
            print('DRY RUN array:',script);return
        command=['sbatch','--parsable',str(script)] if kind=='slurm' else ['qsub',str(script)]
        result=subprocess.run(command,capture_output=True,text=True)
        (script.parent/'submission_stdout.txt').write_text(result.stdout,encoding='utf-8')
        (script.parent/'submission_stderr.txt').write_text(result.stderr,encoding='utf-8')
        if result.returncode:raise RuntimeError(f'Submission failed; inspect {script.parent}/submission_stderr.txt and correct queue/account settings')
        job_id=result.stdout.strip().split(';')[0].splitlines()[-1]
        if not re.fullmatch(r'[0-9][A-Za-z0-9_.\[\]-]*',job_id):raise RuntimeError('Unrecognized scheduler job ID: '+job_id)
        write_json(state_file,dict(scheduler=kind,job_id=job_id,signature=signature,contract_sha256=contract_sha,
                                   config_snapshot=str(snapshot),task_file=str(task_file),submitted_utc=stamp()))
        print('Submitted',kind,job_id,'tasks',len(pending),flush=True)
    poll=max(5,int(config.get('scheduler',{}).get('poll_seconds',30)))
    deadline=time.monotonic()+3600*float(config.get('scheduler',{}).get('max_wait_hours',168))
    inactive_since=None;query_failures=0
    while time.monotonic()<deadline:
        states=[job_status(j,root) for j in jobs]
        print(f'{profile}: complete={states.count("complete")}/{len(states)}, failed={states.count("failed")}',flush=True)
        if all(s=='complete' for s in states):return
        active=scheduler_active(kind,job_id)
        if active is None:
            query_failures+=1
            if query_failures>=3:raise RuntimeError('Scheduler status query failed three times. Tasks were not cancelled; restore access and resume monitoring.')
        else:query_failures=0
        if active is False:
            inactive_since=inactive_since or time.monotonic()
            if time.monotonic()-inactive_since>60:
                raise RuntimeError('Array has ended with missing/failed outputs. Inspect scheduler and per-job logs; rerun after correcting the cause. Completed tasks are retained.')
        else:inactive_since=None
        time.sleep(poll)
    raise RuntimeError('Monitoring timeout; submitted tasks were NOT cancelled. Resume with the same command.')


def compare_profiles(config_path,left,right):
    import numpy as np
    config_path=Path(config_path).resolve();root=config_path.parent;c=read_json(config_path)
    def load(profile):
        entries=read_json(manifest_path(root,c,profile))
        if isinstance(entries,dict):entries=entries['jobs']
        return {(j['model_id'],j['edge'],j['absorber_index'],j['mode']):j for j in entries}
    a,b=load(left),load(right)
    if set(a)!=set(b):raise ValueError('Convergence profiles must have identical absorbing sites and modes')
    records=[]
    for key in a:
        mode=key[-1]
        from run_feff_job import numeric_table
        filename='xmu.dat' if mode=='xanes' else 'chi.dat'
        aa=np.asarray(numeric_table(job_dir(a[key],root)/filename,6 if mode=='xanes' else 2)[0])
        bb=np.asarray(numeric_table(job_dir(b[key],root)/filename,6 if mode=='xanes' else 2)[0])
        if mode=='xanes':
            lo=max(aa[:,0].min(),bb[:,0].min());hi=min(aa[:,0].max(),bb[:,0].max())
            grid=np.linspace(lo,hi,600)
            ya=np.interp(grid,aa[:,0],aa[:,3]);yb=np.interp(grid,bb[:,0],bb[:,3])
        else:
            grid=np.arange(c['convergence'].get('kmin',3),c['convergence'].get('kmax',12)+.01,.05)
            if any(v[:,0].min()>grid[0] or v[:,0].max()<grid[-1]-.002 for v in [aa,bb]):
                raise ValueError('Convergence chi(k) does not cover comparison range')
            ya=np.interp(grid,aa[:,0],aa[:,1])*grid**2;yb=np.interp(grid,bb[:,0],bb[:,1])*grid**2
        if not np.isfinite(ya).all() or not np.isfinite(yb).all():raise ValueError('Nonfinite convergence data')
        rms=float(np.sqrt(np.mean((ya-yb)**2))/max(np.sqrt(np.mean(yb**2)),1e-12))
        tolerance=float(c['convergence'].get(f'relative_rms_tolerance_{mode}',.03))
        records.append(dict(model=key[0],edge=key[1],site=key[2],mode=mode,relative_rms=rms,tolerance=tolerance,passed=rms<=tolerance))
    report=dict(left=left,right=right,passed=all(x['passed'] for x in records),comparisons=records,
                scope='Numerical radius screen on representative sites only; physical model/charge/spin and every production geometry still require review')
    write_json(root/'validation'/f'convergence_{left}_vs_{right}.json',report)
    return report


def doctor(config_path):
    config_path=Path(config_path).resolve();root=config_path.parent
    report={'python':sys.version,'platform':sys.platform,'commands':{},'scheduler_queries':{},'input_hashes_ok':True}
    for cmd in ['sbatch','squeue','qsub','qstat','gfortran','ifort','ifx','make','curl','git']:
        report['commands'][cmd]=shutil.which(cmd)
    for command in [['sinfo','-h','-o','%P|%a|%c|%m|%l'],['qstat','-Q']]:
        if shutil.which(command[0]):
            try:
                res=subprocess.run(command,capture_output=True,text=True,timeout=20)
                report['scheduler_queries'][command[0]]={'returncode':res.returncode,'stdout':res.stdout[:20000],'stderr':res.stderr[:2000]}
            except subprocess.TimeoutExpired:report['scheduler_queries'][command[0]]={'timeout':True}
    for record in read_json(root/'structure_checksums.json'):
        if sha(root/record['file'])!=record['sha256']:report['input_hashes_ok']=False
    write_json(root/'validation'/'environment_report.json',report)
    print(json.dumps(report,indent=2))
    if not report['input_hashes_ok']:raise RuntimeError('Original CIF checksums changed')


def execute(args):
    path=Path(args.config).resolve();root=path.parent;c=read_json(path)
    profiles=c.get('convergence',{}).get('profiles',['pilot','convergence_1','convergence_2'])
    previous=None;chosen=None
    for profile in profiles:
        run_profile(path,profile,args.scheduler,args.allow_local,args.dry_run)
        if args.dry_run:
            print('Dry run stops before numerical convergence and production submission.');return
        if previous is not None:
            report=compare_profiles(path,previous,profile)
            if report['passed']:chosen=profile;break
        previous=profile
    if chosen is None:raise RuntimeError('Radius screening not converged. Inspect validation/convergence_*.json; extend justified parameters/resources before production.')
    production=json.loads(json.dumps(c))
    production['feff'].update(c['profiles'][chosen].get('feff',{}))
    production['numerical_screen']={'selected_profile':chosen,'screened_at':'Top-O1 representative W/Co/Fe sites','time_utc':stamp()}
    production_path=root/'production_config.json';write_json(production_path,production)
    run_profile(production_path,args.profile,args.scheduler,args.allow_local)
    checked_run([sys.executable,SCRIPT_ROOT/'scripts'/'postprocess.py','--config',production_path,'--profile',args.profile])
    write_json(root/'validation'/'workflow_completion.json',{
        'status':'computed_and_postprocessed','profile':args.profile,'production_config_sha256':sha(production_path),
        'completed_utc':stamp(),'data_kind':'calculated_static_geometry_XAS',
        'remaining_scientific_limits':['DFT charge/U/magnetism metadata not supplied','No applied-potential assignment','Finite-radius screen uses selected sites, not every possible model assumption']})
    checked_run([sys.executable,SCRIPT_ROOT/'scripts'/'package_results.py','--config',production_path,'--profile',args.profile])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['doctor','prepare','run','execute','postprocess','compare'])
    p.add_argument('--config',default=str(SCRIPT_ROOT/'config.json'))
    p.add_argument('--profile',default='core')
    p.add_argument('--scheduler',choices=['auto','slurm','pbs','local'],default='auto')
    p.add_argument('--allow-local',action='store_true')
    p.add_argument('--dry-run',action='store_true')
    p.add_argument('--left',default='pilot');p.add_argument('--right',default='convergence_1')
    a=p.parse_args()
    if a.action in {'prepare','run','execute'}:
        with coordinator_lock(Path(a.config).resolve().parent):
            if a.action=='prepare':prepare(a.config,a.profile)
            elif a.action=='run':run_profile(a.config,a.profile,a.scheduler,a.allow_local,a.dry_run)
            else:execute(a)
    elif a.action=='doctor':doctor(a.config)
    elif a.action=='compare':print(json.dumps(compare_profiles(a.config,a.left,a.right),indent=2))
    else:checked_run([sys.executable,SCRIPT_ROOT/'scripts'/'postprocess.py','--config',a.config,'--profile',a.profile])


if __name__=='__main__':
    try:main()
    except (RuntimeError,ValueError,subprocess.CalledProcessError) as exc:
        print('WORKFLOW STOPPED:',exc,file=sys.stderr);raise SystemExit(2)
