"""A single bounded recovery wave; retain failed output and prove orphan exit."""
import argparse,contextlib,hashlib,json,os,shutil,subprocess,sys,time,uuid
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from clusterflow.common import readj,writej,sha,inside
from clusterflow.scheduling import single_execution
from scripts.goat_job import validate_job, INPUT_HASH_FIELDS
from scripts.goat_portable import IMMUTABLE_FILES
from scripts.packed_worker import queue_lock
from scripts.submit_goat_bulk import IDENTITIES
from scripts.expand_goat_20261020 import resource_input
from scripts.recover_goat_stageout_20261020 import validate_terminal
RELATIVE='operations/recovery_cohort_20261020'
NAMES=[('F2_E21_source','04'),('F2_E21_source','07'),('F2_E22_face_2','02'),('F2_S24_source','02'),('F2_S24_source','04'),('F2_S24_source','06'),('F2_E14_source','05'),('F2_E13_cis','03'),('F4_D04_rep_R','02'),('F4_D04_rep_R','03')]
JOBS=[n+'__DMSO__DMSO_02_W01_trace_water_sensitivity__s'+s for n,s in NAMES]
PROBE="""import os,json,sys,time
from pathlib import Path
t=json.loads(sys.argv[1]);matches=[];failures=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit() or int(p.name)==os.getpid():continue
 try:
  if p.stat().st_uid!=os.getuid():continue
  cmd=(p/'cmdline').read_bytes().replace(bytes([0]),b' ').decode(errors='replace');cwd=str((p/'cwd').resolve())
  if int(p.name)==t['pid'] or any(x and (cwd.startswith(x) or x in cmd) for x in (t['stage'],t.get('private_runtime'))):matches.append({'pid':int(p.name),'cwd':cwd,'command':cmd[:200]})
 except FileNotFoundError:pass
 except OSError as e:failures.append(str(e))
print(json.dumps({'matches':matches,'failures':failures,'at_epoch':time.time()}))
"""

def verify_original_inputs(root,row):
 d=inside(root,row['job_directory'])
 if row.get('initial_release') is not True:raise ValueError('Unreleased original start')
 for name,field in dict(INPUT_HASH_FIELDS,**{'provenance.json':'provenance_sha256'}).items():
  if (d/name).is_symlink() or sha(d/name)!=row[field].lower():raise ValueError('Original input identity changed: '+name)

def inspect(root,job,row,queue):
 d=inside(root,row['job_directory']);verify_original_inputs(root,row)
 if (d/'delegated_goat.json').exists():raise ValueError('Already delegated: '+job)
 intents=list((d/'stage_intents').glob('*.json'))
 if not intents:
  if job not in JOBS[-2:] or any((d/n).exists() for n in ('.goat-run.lock','receipt.json','portable_stageout.json')):raise ValueError('Unexpected fresh start')
  validate_job(root,job)
  return dict(reason='FRESH_INDEPENDENT_SEED_AFTER_EMPTY_ENSEMBLE',retry_number=0)
 if len(intents)!=1:raise ValueError('Unexpected repeated attempts')
 intent=readj(intents[0]);stage=Path(intent['stage']);receipt=readj(stage/'receipt.json')
 if receipt.get('status')=='FAILED':
  validate_terminal(root,job,intent)
  if receipt.get('errors')==['ORCA nonzero or absent exit code']:raise ValueError('Recover valid output instead of retry')
  return dict(reason='ONE_RETRY_AFTER_TERMINAL_DAMAGED_OR_MISSING_OUTPUT',retry_number=1,receipt_sha256=sha(stage/'receipt.json'),stage=str(stage),original_status='FAILED')
 if receipt.get('status')!='RUNNING' or job not in JOBS[6:8]:raise ValueError('Unknown producer state')
 h=readj(stage/'runtime_heartbeat.json');claim=queue['active'][job]
 if time.time()-h['updated_epoch']<3600:raise ValueError('Fresh producer heartbeat')
 if claim['hostname']!='i15r4n18' or str(claim['slurm_job_id'])!='124489624':raise ValueError('Changed orphan allocation')
 payload=dict(pid=h['pid'],stage=str(stage),private_runtime=h.get('private_runtime'))
 cmd=['srun','--overlap','--jobid='+str(claim['slurm_job_id']),'--nodes=1','--ntasks=1','--cpus-per-task=1','--mem=256M','--nodelist='+claim['hostname'],'/public/home/qikai/.conda/envs/orca-nmr-py310/bin/python','-c',PROBE,json.dumps(payload)]
 result=subprocess.run(cmd,capture_output=True,text=True,timeout=90,check=True);proof=json.loads(result.stdout)
 if proof['matches'] or proof['failures']:raise ValueError('Orphan producer not conclusively absent: '+str(proof))
 return dict(reason='ONE_RETRY_AFTER_PROVEN_ORPHAN_EXIT',retry_number=1,probe=proof,heartbeat_sha256=sha(stage/'runtime_heartbeat.json'),receipt_sha256=sha(stage/'receipt.json'),stage=str(stage),original_status='RUNNING_INTERRUPTED_WITHOUT_TERMINAL_RECEIPT')

def prepare(root):
 root=Path(root).resolve();dest=root/RELATIVE;base=root/'operations/goat_bulk_20261015'
 if dest.exists():raise ValueError('Recovery cohort already exists; reconcile instead of preparing twice')
 rows={r['job_id']:r for r in readj(root/'registry/jobs.json')}
 with queue_lock(root/'operations/submission_lock'),contextlib.ExitStack() as stack:
  for job in JOBS:stack.enter_context(single_execution(root,'goat_'+hashlib.sha256(job.encode()).hexdigest()))
  queue=readj(base/'queue.json');proofs={j:inspect(root,j,rows[j],queue) for j in JOBS}
  dest.mkdir(parents=True);writej(dest/'PREPARATION_INTENT.json',dict(at_epoch=time.time(),proofs=proofs,authorization='User requested handling current failures and new calculations on 2026-10-20'))
  for name in ('scripts','clusterflow','config'):shutil.copytree(root/name,dest/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
  (dest/'registry').mkdir()
  for name in ('species.json','batch_plan_v5.json'):shutil.copy2(root/'registry'/name,dest/'registry'/name)
  new=[]
  for job in JOBS:
   row=rows[job];source=inside(root,row['job_directory']);target=dest/row['job_directory'];target.mkdir(parents=True)
   for name in IMMUTABLE_FILES:shutil.copy2(source/name,target/name)
   (target/'goat.inp').write_text(resource_input((source/'goat.inp').read_text()));new.append(dict(row,input_sha256=sha(target/'goat.inp')))
  writej(dest/'registry/jobs.json',new)
  for job in JOBS:validate_job(dest,job)
  with queue_lock(base):
   queue=readj(base/'queue.json')
   for job in JOBS:
    source=inside(root,rows[job]['job_directory']);proof=proofs[job]
    if 'heartbeat_sha256' in proof and sha(Path(proof['stage'])/'runtime_heartbeat.json')!=proof['heartbeat_sha256']:raise ValueError('Orphan heartbeat changed')
    token=uuid.uuid4().hex;archive=source/'recovery_history'/('bounded_retry_20261020_'+token);archive.mkdir(parents=True);writej(archive/'proof.json',proof)
    lock=source/'.goat-run.lock'
    if lock.exists():lock.rename(archive/'original_owner.json')
    record=dict(kind='RESOURCE_COHORT_DELEGATION',token=token,job_id=job,canonical_manifest_sha256=sha(root/'registry/jobs.json'),canonical_inputs={n:sha(source/n) for n in IMMUTABLE_FILES},cohort=str(dest),cohort_manifest_sha256=sha(dest/'registry/jobs.json'),cpus=16,created_epoch=time.time(),reason=proof['reason'],retry_number=proof['retry_number'])
    with (source/'.goat-run.lock').open('x') as f:json.dump(record,f)
    with (source/'delegated_goat.json').open('x') as f:json.dump(record,f)
    claim=queue['active'].pop(job,None)
    if claim:writej(archive/'old_active_claim.json',claim)
    cp=base/'claims'/(job+'.json')
    if cp.exists():cp.rename(archive/'old_pool_claim.json')
    result=base/'results'/(job+'.json')
    if result.exists():writej(archive/'old_pool_result.json',readj(result))
    writej(result,dict(job_id=job,result=dict(status='BOUNDED_RECOVERY_DELEGATED',cohort=str(dest),reason=proof['reason']),finished_epoch=time.time()))
   writej(base/'queue.json',queue)
  scope=dict(schema=1,stage='goat',job_ids=JOBS,operations_directory='operations/goat_bulk_20261015',cpus_per_job=16,slots_per_node=2,growth_gib_per_job=8,free_floor_gib=128,legacy_active_jobs=460,identity_sha256={n:sha(dest/n) for n in IDENTITIES},created_epoch=time.time(),canonical_root=str(root),notice='One bounded retry; all historical raw results retained. Two fresh independent lansoprazole seeds.')
  writej(dest/'operations/goat_bulk_20261015/scope.json',scope)
  (dest/'operations/packed').mkdir(parents=True);(dest/'operations/packed/STOP_NEW_CLAIMS').touch();(dest/'operations/ROLLOUT_PAUSED').touch()
  report=dict(status='PREPARED_NOT_SUBMITTED',tasks=len(JOBS),retry_tasks=8,fresh_lansoprazole_seeds=2,proofs=proofs,at_epoch=time.time());writej(dest/'PREPARED.json',report);return report

def submit(root):
 root=Path(root);dest=root/RELATIVE;base=dest/'operations/goat_bulk_20261015';scope=base/'scope.json';readj(dest/'PREPARED.json');site=readj(dest/'config/site_scnet.json')
 with queue_lock(root/'operations/submission_lock'):
  ledger=base/'submission.json'
  if ledger.exists():return readj(ledger)
  cmd=['sbatch','--parsable','--account='+site['slurm']['account'],'--partition='+site['slurm']['partition'],'--qos='+site['slurm']['qos'],'--time=720:00:00','--time-min=720:00:00','--nodes=1','--array=0-4%5','--ntasks-per-node=1','--cpus-per-task=32','--mem=114000M','--job-name=work2-64-recover','--chdir='+str(dest),'--exclude=a10r2n13,k18r2n[10-17]','--output='+str(base/'slurm-%j.out'),'--error='+str(base/'slurm-%j.err'),str(dest/'scripts/scnet_goat_expansion.sbatch')]
  env=dict(os.environ,WORK2_PACKAGE_ROOT=str(dest),WORK2_PYTHON=site['python'],WORK2_PACKED_PROFILE='large_memory',WORK2_RUNTIME_EVIDENCE=str(root/'operations/runtime_124440267.json'),WORK2_GOAT_BULK_BASE=str(base),WORK2_GOAT_BULK_SCOPE=str(scope),WORK2_GOAT_BULK_SCOPE_SHA256=sha(scope))
  subprocess.run(cmd[:-1]+['--test-only',cmd[-1]],capture_output=True,text=True,check=True,timeout=60,env=env)
  record=dict(status='SUBMISSION_INTENT',command=cmd,at_epoch=time.time(),cpus=160,nodes=5);writej(ledger,record)
  try:
   r=subprocess.run(cmd,capture_output=True,text=True,check=True,timeout=60,env=env);jid=r.stdout.strip().split(';')[0]
   if not jid.isdigit():raise ValueError('Unrecognized job ID')
  except Exception as e:writej(ledger,dict(record,status='SUBMISSION_UNCERTAIN',error=repr(e)));raise
  record.update(status='SUBMITTED_NOT_YET_VERIFIED_RUNNING',job_id=jid);writej(ledger,record);return record

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','submit']);args=p.parse_args();print(json.dumps(prepare(ROOT) if args.action=='prepare' else submit(ROOT),ensure_ascii=False,indent=2))
