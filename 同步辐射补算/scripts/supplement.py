#!/usr/bin/env python3
"""Candidate relaxation -> coordination audit -> real FEFF -> TXT and figures."""
from __future__ import annotations
import argparse,hashlib,json,os,re,shlex,shutil,subprocess,sys,time,zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from datetime import datetime,timezone
import workflow as engine

ROOT=Path(__file__).resolve().parents[1]
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n',encoding='utf-8');tmp.replace(p)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def stamp():return datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
def call(script,*args):subprocess.run([sys.executable,str(ROOT/'scripts'/script),*map(str,args)],check=True)
def safe(root,name):
    p=(root/name).resolve()
    if not p.is_relative_to(root.resolve()):raise ValueError('Path leaves project: '+str(name))
    return p

def doctor(config):
    c=read(config);root=config.parent
    commands={n:shutil.which(n) for n in ['sbatch','squeue','qsub','qstat','pw.x','mpirun','mpiexec','mpif90','mpiifort','ifort','gfortran','cmake','git']}
    report={'time_utc':stamp(),'commands':commands,'original_structures':{},'notes':['No heavy computation is performed by doctor','Read platform module/queue policy before execution']}
    for e in read(root/'original_checksums.json'):
        report['original_structures'][e['file']]=sha(root/e['file'])==e['sha256']
    for cmd in [['sinfo','-h','-o','%P|%a|%c|%m|%l'],['qstat','-Q']]:
        if shutil.which(cmd[0]):
            x=subprocess.run(cmd,capture_output=True,text=True,timeout=30)
            report[cmd[0]]={'returncode':x.returncode,'stdout':x.stdout[:12000],'stderr':x.stderr[:2000]}
    write(root/'validation/environment.json',report)
    print(json.dumps(report,indent=2))
    if not all(report['original_structures'].values()):raise ValueError('Original CIF hashes changed')

def runtime_config(config):
    root=config.parent;c=read(config)
    if not c['dft'].get('command'):
        ptr=root/'software/qe_executable.txt'
        exe=ptr.read_text(encoding='utf-8').strip() if ptr.exists() else shutil.which('pw.x')
        if not exe or not Path(exe).is_file():raise RuntimeError('Install/load Quantum ESPRESSO first')
        ranks=int(c['dft_scheduler'].get('ranks',16))
        from probe_qe import probe
        write(root/'software/qe_probe.json',probe(exe,ranks))
        launcher=shutil.which('mpirun') or shutil.which('mpiexec')
        if not launcher and ranks>1:raise RuntimeError('Load the MPI runtime matching pw.x; no serial fallback')
        c['dft']['command']=([launcher,'-np',str(ranks)] if ranks>1 else [])+[exe,'-in','qe.in']
    path=root/'runtime_config.json';write(path,c);return path

def dft_done(job,root):
    directory=safe(root,job['path']);p=directory/'result.json'
    if not p.is_file():return False
    try:
        r=read(p)
        return (r.get('status')=='optimized' and r.get('input_sha256')==sha(directory/'qe.in')
                and r.get('settings_sha256')==job['settings_sha256']
                and r.get('output_sha256')==sha(directory/'qe.out')
                and r.get('relaxed_sha256')==sha(safe(root,r['relaxed_cif'])))
    except (ValueError,KeyError,OSError):return False

def render_dft_array(root,config,c,number,kind,runid):
    settings=c['dft_scheduler'];common=c['scheduler']
    ranks=int(settings['ranks']);mem=int(settings['memory_gb']);parallel=int(settings['max_concurrent'])
    if min(ranks,mem,parallel)<1:raise ValueError('Invalid DFT resource request')
    wall=engine._directive(settings['walltime'],'walltime');logs=root/'logs/dft';logs.mkdir(parents=True,exist_ok=True)
    lines=['#!/usr/bin/env bash']
    if kind=='slurm':
        lines += [f'#SBATCH --job-name=xas_relax',f'#SBATCH --array=0-{number-1}%{parallel}',
                  '#SBATCH --nodes=1',f'#SBATCH --ntasks={ranks}','#SBATCH --cpus-per-task=1',
                  f'#SBATCH --mem={mem}G',f'#SBATCH --time={wall}',f'#SBATCH --output={shlex.quote(str(logs/"%A_%a.out"))}',f'#SBATCH --error={shlex.quote(str(logs/"%A_%a.err"))}']
        for name in ['partition','account']:
            if common.get(name):lines.append(f'#SBATCH --{name}={engine._directive(common[name],name)}')
        index='${SLURM_ARRAY_TASK_ID}'
    else:
        flavor=common.get('pbs_flavor','openpbs');flag='-t' if flavor=='torque' else '-J'
        lines += ['#PBS -N xas_relax','#PBS -S /bin/bash',f'#PBS -l walltime={wall}','#PBS -j oe',f'#PBS -o {shlex.quote(str(logs))}']
        if number>1:lines.append(f'#PBS {flag} 0-{number-1}%{parallel}')
        lines.append(f'#PBS -l nodes=1:ppn={ranks},mem={mem}gb' if flavor=='torque' else f'#PBS -l select=1:ncpus={ranks}:mpiprocs={ranks}:mem={mem}gb')
        for name,flag in [('queue','-q'),('account','-A')]:
            if common.get(name):lines.append(f'#PBS {flag} {engine._directive(common[name],name)}')
        index='${PBS_ARRAY_INDEX:-${PBS_ARRAYID:-0}}'
    lines += ['set -euo pipefail',f'cd {shlex.quote(str(root))}','export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1']
    lines += list(common.get('module_commands',[]))
    cmd=[sys.executable,str(ROOT/'scripts/qe_tasks.py'),'run','--config',str(config),'--task']
    lines.append(shlex.join(cmd)+' "'+index+'"')
    p=root/'submissions'/f'dft_{runid}.sh';p.parent.mkdir(exist_ok=True);p.write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n');return p

def run_dft(config,kind,allow_local=False):
    root=config.parent;c=read(config);jobs=read(root/'dft/jobs.json')
    if not jobs:raise RuntimeError('No DFT jobs')
    state_path=root/'submissions/dft_current.json'
    contract=hashlib.sha256(json.dumps({'jobs':jobs,'dft':c['dft']},sort_keys=True).encode()).hexdigest()
    old=read(state_path) if state_path.exists() else {};active=False
    if old.get('job_id') and old.get('scheduler')==kind and kind!='local':
        active=engine.scheduler_active(kind,old['job_id'])
        if active is None:raise RuntimeError('Cannot query existing DFT array; do not duplicate it')
        if active and old.get('contract')!=contract:raise RuntimeError('DFT settings changed while array active')
    if all(dft_done(j,root) for j in jobs):print('All DFT tasks complete');return
    if kind=='local':
        if not allow_local:raise RuntimeError('Local mode requires --allow-local INSIDE allocated compute session')
        def task(i):
            if not dft_done(jobs[i],root):call('qe_tasks.py','run','--config',config,'--task',i)
        with ThreadPoolExecutor(max_workers=int(c['dft_scheduler']['max_concurrent'])) as pool:list(pool.map(task,range(len(jobs))))
        if not all(dft_done(j,root) for j in jobs):raise RuntimeError('Incomplete local DFT results')
        return
    if active:jobid=old['job_id']
    else:
        runid=stamp();snapshot=root/f'dft_config_{runid}.json';write(snapshot,c)
        script=render_dft_array(root,snapshot,c,len(jobs),kind,runid)
        result=subprocess.run((['sbatch','--parsable'] if kind=='slurm' else ['qsub'])+[str(script)],capture_output=True,text=True)
        write(root/'submissions'/f'dft_submit_{runid}.json',{'returncode':result.returncode,'stdout':result.stdout,'stderr':result.stderr})
        if result.returncode:raise RuntimeError('DFT submission failed; inspect submissions/dft_submit_*.json')
        jobid=result.stdout.strip().split(';')[0]
        if not re.fullmatch(r'[0-9][A-Za-z0-9_.\[\]-]*',jobid):raise RuntimeError('Unrecognized scheduler ID')
        write(state_path,{'scheduler':kind,'job_id':jobid,'contract':contract,'snapshot':snapshot.name})
    print('Monitoring DFT array',jobid,flush=True)
    deadline=time.monotonic()+3600*float(c['scheduler']['max_wait_hours']);inactive=None;errors=0
    while time.monotonic()<deadline:
        n=sum(dft_done(j,root) for j in jobs);print(f'DFT: {n}/{len(jobs)} validated',flush=True)
        if n==len(jobs):return
        active=engine.scheduler_active(kind,jobid)
        errors=errors+1 if active is None else 0
        if errors>=3:raise RuntimeError('Scheduler query failed; existing tasks retained. Restore access and resume')
        if active is False:
            inactive=inactive or time.monotonic()
            if time.monotonic()-inactive>60:raise RuntimeError('DFT array ended with missing/failed relaxations; inspect outputs and resolve before FEFF')
        else:inactive=None
        time.sleep(max(5,int(c['scheduler'].get('poll_seconds',30))))
    raise RuntimeError('Monitor timed out; jobs were not cancelled. Resume same command')

def coordination(config):
    from pymatgen.core import Structure
    root=config.parent;c=read(config);settings=c['coordination'];jobs=read(root/'dft/jobs.json')
    records=[];summaries=[];selected=root/'structures_relaxed';selected.mkdir(exist_ok=True)
    for job in jobs:
        if not dft_done(job,root):raise RuntimeError('Relaxation incomplete: '+job['id'])
        result=read(safe(root,job['path'])/'result.json');source=safe(root,result['relaxed_cif']);s=Structure.from_file(source)
        summaries.append((job['id'],result['metrics']['final_energy_Ry'],result['metrics']['final_energy_eV'],result['metrics']['max_atomic_force_eV_A'],len(s)))
        tungsten=[i for i,x in enumerate(s) if x.specie.symbol=='W']
        if not tungsten:raise RuntimeError('Relaxed model lacks W: '+job['id'])
        for i in tungsten:
            near=sorted(s.get_neighbors(s[i],4.0),key=lambda n:n.nn_distance)
            values=[{'element':n.specie.symbol,'distance_A':float(n.nn_distance),'site_index':int(n.index)} for n in near]
            oxy=[x for x in values if x['element']=='O' and x['distance_A']<=settings['oxygen_cutoff_A']]
            direct=[x for x in values if x['element'] not in {'O','H'} and x['distance_A']<=settings['nonoxygen_cutoff_A']]
            passed=len(oxy)>=settings['minimum_oxygen_neighbors'] and not direct
            records.append({'model':job['id'],'w_site':i,'passed_all_O_screen':passed,'oxygen_distances_A':[x['distance_A'] for x in oxy],
                            'nonoxygen_heavy_neighbors':direct,'all_neighbors_to_4A':values,
                            'cutoff_sensitivity':{str(cut):{e:sum(x['element']==e and x['distance_A']<=cut for x in values) for e in ['O','Co','Fe','W','C','H']} for cut in settings['inspect_cutoffs_A']}})
        target=selected/(job['id']+'.cif')
        if target.exists() and sha(target)!=sha(source):raise RuntimeError('A selected structure changed; preserve old FEFF project before replacing it')
        if not target.exists():shutil.copy2(source,target)
    report={'records':records,'all_passed':all(r['passed_all_O_screen'] for r in records),
            'interpretation':'Distance-based first-shell screen, not a bond-order measurement. W-O FT peak position is not imposed.'}
    write(root/'validation/relaxed_coordination.json',report)
    with (root/'validation/DFT_optimization_summary.txt').open('w',encoding='utf-8',newline='\n') as f:
        f.write('# Total energies cannot be directly ranked across different stoichiometries without reservoir/reference corrections.\nmodel\tenergy_Ry\tenergy_eV\tmax_force_eV_per_A\tn_atoms\n')
        for row in summaries:f.write('\t'.join(map(str,row))+'\n')
    with (root/'validation/relaxed_W_O_distances.txt').open('w',encoding='utf-8',newline='\n') as f:
        f.write('# Calculated relaxed geometry; distances in angstrom\nmodel\tW_site\tneighbor_element\tdistance_A\tneighbor_index\n')
        for record in records:
            for n in record['all_neighbors_to_4A']:f.write(f"{record['model']}\t{record['w_site']}\t{n['element']}\t{n['distance_A']:.8f}\t{n['site_index']}\n")
    if not report['all_passed']:
        if settings.get('failure_policy','report_and_calculate')=='stop':
            raise RuntimeError('Some candidates fail the all-O screen. Preserve/report them; do not force distances. See validation/relaxed_coordination.json')
        print('Some candidates did not retain the target all-O shell. They remain in the comparison as explicitly audited alternative structures; no distance correction is applied.',flush=True)
    return [j['id']+'.cif' for j in jobs]

def make_feff_config(config,models):
    c=read(config);root=config.parent
    chosen=[n for n in models if n in ['W_sub_Fe__hydroxylated__none.cif','W_sub_Fe__Co_deprotonated__EG_at_W.cif']]
    if len(chosen)!=2:raise RuntimeError('Missing paired representative structures for radius screen')
    f={'project':c['project'],'data_kind':c['data_kind'],'paths':{'structures_dir':'structures_relaxed','jobs_dir':'feff_jobs','manifests_dir':'manifests'},
       'software':c['software'],'feff':c['feff'],'scheduler':{**c['scheduler'],**c['feff_scheduler']},'postprocess':c['postprocess'],
       'profiles':{'screen_small':{'models':chosen,'edges':['W_L3','Co_K','Fe_K'],'site_selection':'representative_nearest_w'},
                   'screen_large':{'models':chosen,'edges':['W_L3','Co_K','Fe_K'],'site_selection':'representative_nearest_w','feff':{'cluster_radius_angstrom':12,'scf_radius_angstrom':7,'fms_radius_angstrom':10}},
                   'supplement':{'models':models,'edges':['W_L3','Co_K','Fe_K'],'site_selection':'all'}},
       'convergence':{'relative_rms_tolerance_xanes':0.03,'relative_rms_tolerance_exafs':0.03,'kmin':3,'kmax':12}}
    # Original engine requires a checksum list relative to project root.
    write(root/'structure_checksums.json',[{'file':'structures_relaxed/'+n,'sha256':sha(root/'structures_relaxed'/n),'bytes':(root/'structures_relaxed'/n).stat().st_size} for n in models])
    path=root/'feff_config.json';write(path,f);return path

def package_results(root,export_dir):
    jobs=read(root/'dft/jobs.json');feffjobs=read(root/'manifests/supplement_jobs.json')
    if not all(dft_done(j,root) for j in jobs):raise RuntimeError('Cannot package incomplete DFT')
    if not all(engine.job_status(j,root)=='complete' for j in feffjobs):raise RuntimeError('Cannot package incomplete FEFF')
    files=[]
    for folder in ['structures_original','candidates','structures_relaxed','scripts','docs','validation','manifests','tests']:
        files += [p for p in (root/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc']
    for directory in [safe(root,j['path']) for j in jobs]+[engine.job_dir(j,root) for j in feffjobs]:
        files += [p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in {'.in','.inp','.out','.dat','.json','.log','.cif','.txt'}]
    # Retain all processed numeric data, including complex FT/WT, plus requested exports.
    result_dir=export_dir.parent
    files += [p for p in result_dir.rglob('*') if p.is_file()]
    files += [p for p in root.iterdir() if p.is_file() and p.suffix in {'.json','.md','.txt','.sh'}]
    for p in [root/'software/python_environment.lock.txt',root/'software/pseudo/provenance.json',root/'software/qe_source_commit.txt',root/'software/qe_probe.json',root/'software/feff10/install_complete.json']:
        if p.exists():files.append(p)
    files=sorted(set(files));out=root/'deliverables'/('CALCULATED_W_O_SUPPLEMENT_'+stamp()+'.zip');out.parent.mkdir(exist_ok=True)
    manifest=[{'file':p.relative_to(root).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files]
    with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in files:z.write(p,p.relative_to(root).as_posix())
        z.writestr('RESULTS_SHA256.json',json.dumps(manifest,indent=2))
    with zipfile.ZipFile(out) as z:
        if z.testzip() is not None:raise RuntimeError('Result ZIP integrity failed')
    print('FINAL RESULT PACKAGE:',out);return out

def execute(config,requested,allow_local=False,prepare_only=False):
    root=config.parent;c=read(config)
    call('build_candidates.py','--config',config)
    if prepare_only:
        print('Candidate preparation finished. No DFT/FEFF result claimed.');return
    kind=engine.scheduler_kind(requested,c)
    runtime=runtime_config(config)
    call('qe_tasks.py','prepare','--config',runtime)
    run_dft(runtime,kind,allow_local)
    call('qe_tasks.py','collect','--config',runtime)
    models=coordination(runtime)
    feff=make_feff_config(runtime,models)
    for profile in ['screen_small','screen_large']:engine.run_profile(feff,profile,kind,allow_local)
    comparison=engine.compare_profiles(feff,'screen_small','screen_large')
    if not comparison['passed']:raise RuntimeError('New-model FEFF radius screen failed. Review/extend radii; do not relax tolerance to obtain desired curves')
    f=read(feff);f['feff'].update(f['profiles']['screen_large']['feff']);production=root/'feff_production_config.json';write(production,f)
    engine.run_profile(production,'supplement',kind,allow_local)
    call('postprocess.py','--config',production,'--profile','supplement')
    results=sorted((root/c['postprocess']['output_dir']/'supplement').glob('COMPLETE_*'))
    if not results:raise RuntimeError('No complete postprocessing output')
    result=results[-1]
    call('export_supplement.py','--config',config,'--result-dir',result)
    exports=sorted(result.glob('FINAL_EXPORT_*'))
    if not exports:raise RuntimeError('No requested TXT/image export')
    if read(exports[-1]/'export_report.json').get('status')!='complete':raise RuntimeError('Final TXT/image export incomplete')
    write(root/'validation/workflow_completion.json',{'status':'computed_relaxed_and_exported','time_utc':stamp(),'result_dir':result.relative_to(root).as_posix(),
       'condition_mapping':'not_established; structural model labels only',
       'all_O_coordination_hypothesis_retained':read(root/'validation/relaxed_coordination.json')['all_passed'],
       'active_site_conclusion':'Not predetermined. Element-specific spectral change alone is not proof of an exclusive active site.'})
    package_results(root,exports[-1])

def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['doctor','prepare','execute','coordination','package'])
    p.add_argument('--config',default=str(ROOT/'config.json'));p.add_argument('--scheduler',default='auto',choices=['auto','slurm','pbs','local']);p.add_argument('--allow-local',action='store_true');p.add_argument('--export-dir')
    a=p.parse_args();config=Path(a.config).resolve()
    if a.action=='doctor':doctor(config)
    elif a.action=='coordination':coordination(config)
    elif a.action=='package':package_results(config.parent,Path(a.export_dir).resolve())
    else:
        with engine.coordinator_lock(config.parent):execute(config,a.scheduler,a.allow_local,a.action=='prepare')

if __name__=='__main__':
    try:main()
    except (ValueError,RuntimeError,subprocess.CalledProcessError) as exc:
        print('SUPPLEMENT STOPPED:',exc,file=sys.stderr);raise SystemExit(2)
