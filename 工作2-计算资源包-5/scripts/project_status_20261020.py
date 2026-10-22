"""Read-only account allocation and canonical/delegated GOAT evidence audit."""
import collections,datetime,json,re,shutil,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from clusterflow.common import readj,writej,sha
from clusterflow.observation import observe


def optional(path):
 try:return readj(path)
 except FileNotFoundError:return None


def audit(root):
 root=Path(root).resolve();now=time.time();site=readj(root/'config/site_scnet.json')
 observation_errors=[]
 def observed(path, missing_ok=False):return observe(path,observation_errors,missing_ok=missing_ok)
 command=['squeue','-A',site['slurm']['account'],'-r','-h','-o','%i|%j|%T|%C|%D|%Z|%S|%e']
 output=subprocess.check_output(command,text=True,timeout=60);alloc=[];totals=collections.Counter()
 for line in output.splitlines():
  f=line.split('|');project=f[5]==str(root) or f[5].startswith(str(root)+'/');r=dict(job_id=f[0],name=f[1],state=f[2],cpus=int(f[3]),nodes=int(f[4]),workdir=f[5],start=f[6],end=f[7],project=project);alloc.append(r)
  kind='project' if project else 'other'
  if f[2] in ('RUNNING','COMPLETING','CONFIGURING'):
   totals[kind+'_allocated_cpus']+=r['cpus'];totals['account_allocated_cpus']+=r['cpus']
  elif f[2]=='PENDING':totals[kind+'_pending_cpus']+=r['cpus']
 roots=[root]+[root/'operations'/n for n in ('expansion_20261020','recovery_cohort_20261020') if (root/'operations'/n/'registry/jobs.json').exists()]
 live=[];stale=[];outcomes=collections.Counter();initial=collections.Counter();failed=[];delegated=0;unique=set();registries={str(r):readj(r/'registry/jobs.json') for r in roots};hashes={str(r):sha(r/'registry/jobs.json') for r in roots}
 for workspace in roots:
  for row in registries[str(workspace)]:
   d=workspace/row['job_directory'];job=row['job_id']
   delegation=optional(d/'delegated_goat.json')
   if delegation:
    c=delegation.get('cohort');assert c in registries and delegation['cohort_manifest_sha256']==hashes[c]
    delegated+=1;continue
   if job in unique:raise ValueError('Duplicate independent owner in audited cohorts: '+job)
   unique.add(job);before_errors=len(observation_errors);receipt=observed(d/'receipt.json',True);owner=observed(d/'.goat-run.lock',True);observed_state='PENDING_OR_NOT_ATTEMPTED'
   if owner and owner.get('portable_stage_root'):
    stage=Path(owner['portable_stage_root'])/'job';h=observed(stage/'runtime_heartbeat.json',True);age=max(0,time.time()-(h or {}).get('updated_epoch',0))
    inp=(d/'goat.inp').read_text();cpus=int(re.findall(r'(?i)\bnprocs\s+(\d+)',inp)[0]);p=stage/'goat.out';size=p.stat().st_size if p.exists() else 0
    r=dict(job_id=job,species_id=row['species_id'],compound_id=row['compound_id'],cpus=cpus,hostname=owner.get('hostname'),slurm_job_id=owner.get('slurm_job_id'),heartbeat_age_seconds=round(age,1),cohort=str(workspace),output_bytes=size,stage=str(stage))
    (live if age<180 else stale).append(r)
    observed_state='RUNNING_FRESH_HEARTBEAT' if age<180 else 'OWNER_WITHOUT_FRESH_HEARTBEAT'
   if receipt and receipt.get('status')!='RUNNING':
    status=receipt['status']
    if (d/'technical_recovery_audit.json').exists():status=readj(d/'technical_recovery_audit.json')['status']
    outcomes[status]+=1
    observed_state=status
    if status=='FAILED':
     disposition=optional(d/'attempt_disposition.json')
     if disposition and disposition.get('source_receipt_sha256')!=sha(d/'receipt.json'):raise ValueError('Failed-seed disposition binding changed')
     failed.append(dict(job_id=job,errors=receipt.get('errors',[]),disposition=disposition))
   if len(observation_errors)>before_errors and observed_state=='PENDING_OR_NOT_ATTEMPTED':observed_state='UNREADABLE_RUNTIME_RECORD'
   if row.get('initial_release'):initial[observed_state]+=1
 pools=[]
 for workspace in roots:
  for name in ('operations/goat_bulk_20261015','operations/pilot_20261020','operations/retry_20261022'):
   b=workspace/name;scope=optional(b/'scope.json')
   if not scope:continue
   q=optional(b/'queue.json') or dict(active={},cursor=0,skipped=[]);nodes=[v for p in (b/'nodes').glob('*.json') if (v:=observed(p)) is not None]
   pools.append(dict(workspace=str(workspace),base=str(b),scope_tasks=len(scope['job_ids']),cursor=q['cursor'],unclaimed=len(scope['job_ids'])-q['cursor'],claims=len(q['active']),skipped=len(q['skipped']),submission=optional(b/'submission.json'),nodes_rejected=[n for n in nodes if n.get('status')=='NODE_REJECTED_BEFORE_QM'],fresh_worker_nodes=sum(n.get('status')=='RUNNING' and now-n.get('updated_epoch',0)<180 for n in nodes),results=dict(collections.Counter((observed(p) or {'result':{'status':'UNREADABLE_RECORD'}})['result']['status'] for p in (b/'results').glob('*.json')))))
 risks=[]
 for a in alloc:
  if not a['project'] or a['state']!='RUNNING':continue
  start=datetime.datetime.fromisoformat(a['start']);end=datetime.datetime.fromisoformat(a['end']);duration=(end-start).total_seconds()/3600
  a['actual_walltime_hours']=round(duration,2);a['remaining_hours']=round((end.timestamp()-now)/3600,2)
  if duration<719:risks.append(dict(job_id=a['job_id'],actual_walltime_hours=a['actual_walltime_hours'],remaining_hours=a['remaining_hours'],notice='Actual deadline shorter than requested 30 days; inspect before admitting more work.'))
 report=dict(status='LIVE_PROJECT_AND_ACCOUNT_AUDIT',checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),cpu=dict(totals),allocations=alloc,fresh_goat_jobs=len(live),fresh_goat_configured_cpus=sum(r['cpus'] for r in live),active_compounds=len({r['compound_id'] for r in live}-{'TMS'}),live=live,stale=stale,terminal_outcomes=dict(outcomes),unrecovered_failed=failed,delegated_tasks=delegated,pools=pools,walltime_risks=risks,shared_free_gib=round(shutil.disk_usage(site['shared_scratch']).free/2**30,2),notice='Allocated CPUs, configured live-QM CPUs and utilization are distinct. Heartbeats do not establish 100% CPU utilization or scientific acceptance.')
 report['initial_start_states']=dict(initial)
 report['observation_errors']=observation_errors
 report['monitor_health']='DEGRADED_RECORDS_REPORTED' if observation_errors else 'OK'
 writej(root/'operations/project_status_20261020.json',report);return report

if __name__=='__main__':print(json.dumps(audit(ROOT),ensure_ascii=False,indent=2))
