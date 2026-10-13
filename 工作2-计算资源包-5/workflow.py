"""工作2-计算资源包-（5）: target final stable minima before NMR handoff."""
from pathlib import Path
import argparse
import json
import sys
import math
from clusterflow.common import config,hash_protocol,inside,readj,writej,sha

ROOT=Path(__file__).resolve().parent
DEFAULT_SCENARIO='DMSO_02_W01_trace_water_sensitivity'
STAGES=['prescreen','optfreq','shielding','j_PBE_pcJ3','j_PBE0_pcJ2']

def display_spectrum_settings(settings, linewidth_h=None, linewidth_c=None):
    """Display-only overrides; never mutate the scientific protocol or QM hash."""
    result=dict(settings['spectrum'])
    for key,value in [('proton_linewidth_Hz',linewidth_h),('carbon_linewidth_Hz',linewidth_c)]:
        if value is not None:
            if not math.isfinite(value) or value<=0:
                raise ValueError('Linewidth must be finite and positive in Hz')
            result[key]=float(value)
    return result

def summary():
    p=config(ROOT); jobs=readj(ROOT/'registry/jobs.json')
    ss=readj(ROOT/'registry/species.json')
    states={f.parent.name:readj(f).get('last_status','STATE_PREPARED') for f in (ROOT/'automation').glob('*/state.json')}
    downstream=list((ROOT/'jobs').glob('*/*/*/*/*/job.json')) if (ROOT/'jobs').exists() else []
    downstream+=list((ROOT/'automation').glob('*/wave_*/workspace/jobs/*/*/*/*/*/job.json'))
    return {'package':p['package_name'],'protocol_hash':hash_protocol(p),
            'products':len({s['compound_id'] for s in ss}-{'TMS'}),
            'models_including_TMS':len(ss),'GOAT_initial_inputs':len(jobs),
            'scenario_prerequisites':p['scenario_prerequisites'],
            'example_scenario':DEFAULT_SCENARIO,
            'models':[{'species_id':s['species_id'],'compound_id':s['compound_id']} for s in ss],
            'downstream_jobs_prepared':len(downstream),'automatic_controller_states':states,
            'stable_target_definition':'FINAL_DFT_UNIQUE_VERIFIED_MINIMA_NOT_STARTS',
            'scientific_convergence_approved':False,
            'status':'AUTOMATION_STATES_PRESENT_INSPECT_PER_SPECIES' if states else 'CODE_AVAILABLE_REAL_SAMPLING_DFT_AND_BENCHMARK_NOT_RUN'}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    sub=ap.add_subparsers(dest='command',required=True)
    sub.add_parser('status')
    sub.add_parser('auto-plan',help='Read the stable-conformer targets; does not run or prepare QM')
    for name in ['auto-step','auto-run']:
        q=sub.add_parser(name,help='Automatic handoff; no QM unless --execute --acknowledge-models')
        choice=q.add_mutually_exclusive_group(required=True)
        choice.add_argument('--species');choice.add_argument('--all',action='store_true')
        q.add_argument('--execute',action='store_true');q.add_argument('--acknowledge-models',action='store_true')
        q.add_argument('--orca');q.add_argument('--xtb');q.add_argument('--scratch-root')
        q.add_argument('--max-steps',type=int,default=10000)
        q.add_argument('--poll-seconds',type=float,default=15.0)
    q=sub.add_parser('preflight');q.add_argument('--orca',required=True)
    q=sub.add_parser('draft-review');q.add_argument('--job-id',required=True);q.add_argument('--out',required=True)
    q=sub.add_parser('import-goat');q.add_argument('--job-id',required=True);q.add_argument('--review',required=True)
    for name in ['prepare','advance']:
        q=sub.add_parser(name);q.add_argument('--species',required=True);q.add_argument('--scenario',default=DEFAULT_SCENARIO)
        if name=='advance':q.add_argument('--from-stage',choices=['prescreen','optfreq'],required=True)
    q=sub.add_parser('manifest');q.add_argument('--stage',choices=STAGES,required=True);q.add_argument('--species');q.add_argument('--scenario',default=DEFAULT_SCENARIO)
    q=sub.add_parser('goat-manifest');q.add_argument('--species');q.add_argument('--scenario',default=DEFAULT_SCENARIO)
    q=sub.add_parser('inspect');q.add_argument('--job',required=True)
    q=sub.add_parser('diagnose');q.add_argument('--job',required=True)
    q=sub.add_parser('draft-issue-review');q.add_argument('--job',required=True);q.add_argument('--out',required=True)
    q=sub.add_parser('multiplet');q.add_argument('--panel',required=True)
    for name in ['run','run-task']:
        q=sub.add_parser(name);q.add_argument('--orca',required=True);q.add_argument('--acknowledge-models',action='store_true');q.add_argument('--dry-run',action='store_true');q.add_argument('--retry',action='store_true');q.add_argument('--timeout',type=float);q.add_argument('--scratch-root')
        if name=='run':q.add_argument('--job',required=True)
        else:q.add_argument('--manifest',required=True);q.add_argument('--index',type=int,required=True)
    q=sub.add_parser('benchmark');q.add_argument('--panel',required=True);q.add_argument('--targets',nargs='+',required=True)
    q=sub.add_parser('collect');q.add_argument('--species',required=True);q.add_argument('--scenario',default=DEFAULT_SCENARIO)
    q.add_argument('--j-method',choices=['PBE_pcJ3','PBE0_pcJ2']);q.add_argument('--exploratory',action='store_true');q.add_argument('--benchmark-gate');q.add_argument('--convergence-gate');q.add_argument('--no-plots',action='store_true')
    q.add_argument('--linewidth-h',type=float,help='Display-only 1H Lorentzian FWHM in Hz; does not change QM protocol')
    q.add_argument('--linewidth-c',type=float,help='Display-only 13C Lorentzian FWHM in Hz; does not change QM protocol')
    a=ap.parse_args()
    if a.command=='status':r=summary()
    elif a.command=='auto-plan':
        r=readj(ROOT/'registry/batch_plan_v5.json')
    elif a.command in ['auto-step','auto-run']:
        from clusterflow.autoflow import step
        import time
        if a.execute and not a.acknowledge_models:
            raise ValueError('Use --acknowledge-models before executing unvalidated solvent/structure models')
        if a.max_steps<1 or not math.isfinite(a.poll_seconds) or not 0<a.poll_seconds<=60:
            raise ValueError('max-steps must be positive; poll-seconds must be within (0,60]')
        rows=readj(ROOT/'registry/batch_plan_v5.json')['species']
        active=[x['species_id'] for x in rows if x['compound_id']!='TMS'] if a.all else [a.species]
        results={}; executed=0
        limit=a.max_steps if a.command=='auto-run' and a.execute else len(active)
        for index in range(limit):
            if not active:break
            sid=active.pop(0)
            try:
                result=step(ROOT,sid,execute=a.execute,orca=a.orca,xtb=a.xtb,scratch_root=a.scratch_root)
            except (ValueError,RuntimeError,FileNotFoundError) as exc:
                result={'status':'AUTOMATION_ERROR_REVIEW_REQUIRED','species_id':sid,'reason':str(exc),'quantum_jobs_launched':0}
            results[sid]=result
            executed+=result.get('quantum_jobs_launched',0)
            print(json.dumps(result,ensure_ascii=False,default=str),flush=True)
            status=result['status']
            stop=any(word in status for word in ('REVIEW','BLOCKED','ERROR'))
            if not stop and a.command=='auto-run' and a.execute:
                active.append(sid)
                if not result.get('quantum_jobs_launched'):
                    time.sleep(a.poll_seconds)
        r={'status':'AUTOMATIC_HANDOFF_SESSION_STOPPED','reports':results,'remaining_species':active,
           'quantum_jobs_launched':executed,'scientific_convergence_approved':False,
           'notice':'Review terminal statuses; a stable-count target is not an accuracy or sampling-convergence certificate.'}
    elif a.command=='preflight':
        from clusterflow.runner import static_preflight
        r=static_preflight(ROOT,a.orca)
    elif a.command=='goat-manifest':
        p=config(ROOT)
        if a.scenario not in p['independent_scenarios']:raise ValueError('Unknown scenario')
        rows=[j for j in readj(ROOT/'registry/jobs.json') if j['scenario_id']==a.scenario and j['initial_release'] and (not a.species or j['species_id']==a.species)]
        if not rows:raise ValueError('No matched GOAT starts')
        target=ROOT/'manifests'/('goat_'+a.scenario+'_'+hash_protocol(rows)+'.json')
        writej(target,rows);r={'manifest':str(target),'jobs':len(rows),'submitted':False,'continuous_only_prerequisite':False}
    elif a.command in ['draft-review','import-goat','prepare','advance','manifest']:
        from clusterflow import pipeline
        if a.command=='draft-review':
            r=pipeline.draft_goat_review(ROOT,a.job_id);target=inside(ROOT,a.out)
            if target.exists():raise ValueError('Review draft path exists; will not overwrite')
            writej(target,r);r={'draft':str(target),'approved':False}
        elif a.command=='import-goat':r=pipeline.import_goat(ROOT,a.job_id,Path(a.review))
        elif a.command=='prepare':r=pipeline.prepare_prescreen(ROOT,a.species,a.scenario)
        elif a.command=='advance':r=pipeline.advance(ROOT,a.species,a.scenario,a.from_stage)
        else:
            rows=pipeline.manifest(ROOT,a.stage,a.species,a.scenario)
            target=ROOT/'manifests'/(a.stage+'_'+hash_protocol(rows)+'.json')
            writej(target,rows);r={'manifest':str(target),'jobs':len(rows),'submitted':False}
    elif a.command=='inspect':
        from clusterflow.gates import inspect_job
        r=inspect_job(inside(ROOT,a.job),config(ROOT));r.pop('mol',None)
    elif a.command=='diagnose':
        from clusterflow.issue_review import diagnose_job
        r=diagnose_job(inside(ROOT,a.job),config(ROOT))
    elif a.command=='draft-issue-review':
        from clusterflow.issue_review import draft_scientific_review
        target=inside(ROOT,a.out)
        if target.exists():raise ValueError('Review draft exists; refusing overwrite')
        r=draft_scientific_review(inside(ROOT,a.job),config(ROOT));writej(target,r)
        r={'draft':str(target),'approved':False}
    elif a.command=='multiplet':
        from clusterflow.multiplet import validate_multiplets
        r=validate_multiplets(ROOT,Path(a.panel))
    elif a.command in ['run','run-task']:
        from clusterflow.runner import run_one
        if a.command=='run-task':
            rows=readj(inside(ROOT,a.manifest))
            if not 0<=a.index<len(rows):raise ValueError('Array index outside manifest')
            row=rows[a.index];relative=row.get('job_dir',row.get('job_directory'))
            if not relative:raise ValueError('Manifest lacks job directory')
        else:relative=a.job
        r=run_one(ROOT,relative,a.orca,acknowledge=a.acknowledge_models,dry_run=a.dry_run,timeout=a.timeout,retry=a.retry,scratch_root=a.scratch_root)
    elif a.command=='benchmark':
        from clusterflow.benchmark import run_benchmark
        r=run_benchmark(ROOT,Path(a.panel),a.targets)
    else:
        from clusterflow.collect import collect_species
        from clusterflow.policy import resolve_j_method
        display_settings=display_spectrum_settings(config(ROOT),a.linewidth_h,a.linewidth_c)
        a.j_method=resolve_j_method(ROOT,a.species,a.j_method)
        r=collect_species(ROOT,a.species,a.scenario,a.j_method,exploratory=a.exploratory,benchmark_gate=Path(a.benchmark_gate) if a.benchmark_gate else None,convergence_gate=Path(a.convergence_gate) if a.convergence_gate else None)
        if not a.no_plots:
            from clusterflow.spectra import make_spectrum,export_spectrum
            from datetime import datetime,timezone
            import uuid
            rows=r.get('spectrum_atom_rows',r['atom_rows'])
            folder=ROOT/'results'/hash_protocol(config(ROOT))/a.species/a.scenario/'spectra'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')+'_'+uuid.uuid4().hex[:8])
            status='EXPLORATORY CALCULATED' if a.exploratory else 'CALCULATED'
            plots={}
            for nuc in ['1H','13C']:
                spec=make_spectrum(rows,r['j_pairs'],nuc,display_settings)
                writej(folder/(nuc+'_model.json'),{k:v for k,v in spec.items() if k not in ['x_ppm','y']})
                plots[nuc]=export_spectrum(spec,folder/nuc,a.species+' | '+a.scenario+' | '+status)
            r['spectra']=plots
            r['spectrum_display_settings']=display_settings
            writej(folder/'display_settings.json',{'settings':display_settings,'prediction_path':r['result_path'],'prediction_sha256':sha(inside(ROOT,r['result_path'])),'display_only_no_QM_rerun':True})
    print(json.dumps(r,ensure_ascii=False,indent=2,default=lambda x:x.tolist() if hasattr(x,'tolist') else str(x)))

if __name__=='__main__':
    try:main()
    except (ValueError,RuntimeError,FileNotFoundError) as exc:
        print('STOP: '+str(exc),file=sys.stderr);sys.exit(2)
