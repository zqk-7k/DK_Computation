"""Replace only a never-started pending gang with independent single-node jobs."""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from clusterflow.common import readj,writej,sha
from scripts.packed_worker import queue_lock


def no_started_work(cohort,base,queue_lines):
 lines=[l.split('|') for l in queue_lines.splitlines() if l.strip()]
 if not lines or any(l[1]!='PENDING' for l in lines):raise ValueError('Every old array element must still be PENDING')
 q=readj(base/'queue.json') if (base/'queue.json').exists() else {'active':{},'cursor':0}
 if q['active'] or q['cursor'] or any((base/'results').glob('*.json')):raise ValueError('Pool already consumed work')
 rows={r['job_id']:r for r in readj(cohort/'registry/jobs.json')}
 for j in readj(base/'scope.json')['job_ids']:
  d=cohort/rows[j]['job_directory']
  if any((d/n).exists() for n in ('.goat-run.lock','receipt.json','portable_stageout.json')) or any((d/'stage_intents').glob('*.json')):raise ValueError('A computation already has execution evidence')
 return lines


def repack(root,label):
 root=Path(root).resolve();cohort=root/'operations'/('expansion_20261020' if label=='main' else 'recovery_cohort_20261020');base=cohort/'operations/goat_bulk_20261015';ledger=base/'submission.json';site=readj(cohort/'config/site_scnet.json')
 with queue_lock(root/'operations/submission_lock'):
  old=readj(ledger)
  if old.get('layout')=='SINGLE_NODE_ARRAY':return old
  if old['status']!='SUBMITTED_NOT_YET_VERIFIED_RUNNING':raise ValueError('Reconcile original submission state first')
  jid=old['job_id'];cmd=['squeue','-r','-h','-j',jid,'-o','%i|%T']
  before=subprocess.check_output(cmd,text=True,timeout=60);no_started_work(cohort,base,before)
  nodes=64 if label=='main' else 5
  newcmd=[s for s in old['command'] if not s.startswith(('--nodes=','--array=')) and not s.endswith('.sbatch')]+['--nodes=1','--array=0-'+str(nodes-1)+'%'+str(nodes)]
  if '--time-min=720:00:00' not in newcmd:raise ValueError('Full actual 30-day minimum required')
  if old['cpus']!=nodes*32:raise ValueError('CPU total mismatch')
  count=subprocess.check_output(['squeue','-A',site['slurm']['account'],'-r','-h','-o','%i'],text=True,timeout=60)
  if len(count.splitlines())-len(before.splitlines())+nodes>=190:raise ValueError('Queue headroom insufficient')
  launcher=cohort/'scripts/scnet_goat_expansion.sbatch'
  subprocess.run(newcmd+['--test-only',str(launcher)],capture_output=True,text=True,check=True,timeout=60)
  history=base/'submission_history';history.mkdir(exist_ok=True)
  if (history/('repack_'+jid+'.json')).exists():raise ValueError('Repacking already attempted; reconcile before further actions')
  writej(history/('repack_'+jid+'.json'),dict(status='INTENT',original=old,command=newcmd,queue_before=before,at_epoch=time.time()))
  subprocess.run(['scontrol','hold',jid],capture_output=True,text=True,check=True,timeout=60)
  try:
   held=subprocess.check_output(cmd,text=True,timeout=60);no_started_work(cohort,base,held)
  except Exception:
   subprocess.run(['scontrol','release',jid],capture_output=True,text=True,timeout=60);raise
  subprocess.run(['scancel',jid],capture_output=True,text=True,check=True,timeout=60)
  for _ in range(12):
   remaining=subprocess.check_output(cmd,text=True,timeout=60)
   if not remaining.strip():break
   time.sleep(2)
  else:raise ValueError('Old pending allocation not conclusively gone; new submission not attempted')
  writej(history/('cancelled_'+jid+'.json'),dict(original=old,status='CANCELLED_BEFORE_ANY_QM',queue_before=before,queue_after=remaining,at_epoch=time.time()))
  scope=base/'scope.json';env=dict(os.environ,WORK2_PACKAGE_ROOT=str(cohort),WORK2_PYTHON=site['python'],WORK2_PACKED_PROFILE='large_memory',WORK2_RUNTIME_EVIDENCE=str(root/'operations/runtime_124440267.json'),WORK2_GOAT_BULK_BASE=str(base),WORK2_GOAT_BULK_SCOPE=str(scope),WORK2_GOAT_BULK_SCOPE_SHA256=sha(scope))
  result=dict(status='SUBMISSION_INTENT',layout='SINGLE_NODE_ARRAY',command=newcmd,at_epoch=time.time(),nodes=nodes,cpus=nodes*32,replaced_pending_job=jid,scope_sha256=sha(scope),scientific_inputs_unchanged=True)
  writej(ledger,result)
  try:
   p=subprocess.run(newcmd+[str(launcher)],capture_output=True,text=True,check=True,timeout=60,env=env);newid=p.stdout.strip().split(';')[0]
   if not newid.isdigit():raise ValueError('Unknown sbatch output')
  except Exception as e:writej(ledger,dict(result,status='SUBMISSION_UNCERTAIN',error=repr(e)));raise
  result.update(status='SUBMITTED_NOT_YET_VERIFIED_RUNNING',job_id=newid);writej(ledger,result);return result

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('label',choices=['main','recovery']);a=p.parse_args();print(json.dumps(repack(ROOT,a.label),ensure_ascii=False,indent=2))
