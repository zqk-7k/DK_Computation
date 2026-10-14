"""Durable multi-process, per-node packed controller for this 64-product package.

Execute ONLY inside the supplied Slurm launcher. Finite parents may be chained
with afterany; persisted work is revisited on the next parent. Never cancel a
healthy calculation, fabricate a receipt, clear unknown locks, or fill the QOS
with one Slurm job per quantum calculation.
"""
import argparse
from contextlib import contextmanager
from datetime import datetime
import fcntl
import multiprocessing
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import time
import uuid
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from clusterflow.common import readj,writej
from clusterflow.runtime import site_config,owner_record,check_storage,wrapper_preflight,process_start
from clusterflow.operations import scheduler_terminal,retire_terminal_lock,coverage,recover_stageout,retry_dft_infrastructure,diagnosis


@contextmanager
def queue_lock(base):
    base.mkdir(parents=True,exist_ok=True)
    with (base/'queue.mutex').open('a') as stream:
        fcntl.flock(stream,fcntl.LOCK_EX)
        try:yield
        finally:fcntl.flock(stream,fcntl.LOCK_UN)


def end_time():
    value=os.environ.get('SLURM_JOB_END_TIME')
    if value and value.isdigit():return int(value)
    result=subprocess.run(['scontrol','show','job',os.environ['SLURM_JOB_ID'],'-o'],
                          capture_output=True,text=True,timeout=20,check=True)
    for part in result.stdout.split():
        if part.startswith('EndTime='):
            return datetime.fromisoformat(part.split('=',1)[1]).timestamp()
    raise ValueError('Cannot verify actual scheduler EndTime; refuse blind long-running claim')


def claim(base,site,species):
    with queue_lock(base):
        claims=base/'claims';claims.mkdir(exist_ok=True)
        active=[]
        for path in claims.glob('*.json'):
            record=readj(path)
            local_live=(record.get('hostname')==socket.gethostname() and record.get('process_start')
                        and process_start(record.get('pid'))==record['process_start'])
            if scheduler_terminal(record) and not local_live:
                archive=base/'claim_history';archive.mkdir(exist_ok=True)
                path.rename(archive/(path.stem+'_'+uuid.uuid4().hex+'.json'))
            else:active.append(record)
        snapshot=check_storage(site['shared_scratch'],site['disk']['admit_free_gib'])
        budget=site['disk']['unmeasured_task_growth_gib']
        if snapshot['free_bytes']/1024**3 < site['disk']['admit_free_gib']+sum(r['growth_gib'] for r in active)+budget:
            return None
        cursor_path=base/'cursor.json'
        cursor=readj(cursor_path).get('index',0) if cursor_path.is_file() else 0
        busy={r['species_id'] for r in active}
        for offset in range(len(species)):
            index=(cursor+offset)%len(species);sid=species[index]
            if sid in busy:continue
            blocked=base/'blocked'/(sid+'.json')
            if blocked.is_file():
                reason=readj(blocked)
                if not reason.get('retry_on_different_node') or reason.get('hostname')==socket.gethostname():continue
            directory=ROOT/'automation'/sid
            lock=directory/'.autoflow.lock'
            if lock.exists() and not retire_terminal_lock(lock):continue
            state=readj(directory/'state.json') if (directory/'state.json').is_file() else {}
            if state.get('phase')=='REVIEW':continue
            record={**owner_record(uuid.uuid4().hex),'species_id':sid,'growth_gib':budget}
            path=claims/(sid+'.json')
            # Outer flock + exclusive write protects against duplicate claims.
            with path.open('x') as stream:
                import json
                json.dump(record,stream)
            writej(cursor_path,{'index':index+1})
            return path,record
        return None


def worker(slot,slots,deadline,max_steps):
    site=site_config();base=ROOT/'operations/packed'
    os.environ['SLURM_NTASKS']='1'  # One serial driver per isolated worker, never mpirun the driver.
    os.environ['SLURM_CPUS_PER_TASK']=str(site['packed']['slot_allocated_cpus'])
    os.environ['WORK2_SLOT_MEMORY_MB']=str(site['packed']['slot_memory_mb'])
    if hasattr(os,'sched_getaffinity'):
        available=sorted(os.sched_getaffinity(0))
        size=len(available)//slots
        if size<site['packed']['slot_allocated_cpus']:raise ValueError('Insufficient isolated CPU allocation')
        os.sched_setaffinity(0,available[slot*size:(slot+1)*size])
    from clusterflow.autoflow import step
    species=[r['species_id'] for r in readj(ROOT/'registry/batch_plan_v5.json')['species'] if r['compound_id']!='TMS']
    idle=0
    for _ in range(max_steps):
        if os.environ.get('WORK2_DRAIN_REQUESTED')=='1' or time.time()>deadline-site['slurm']['drain_before_end_seconds']:break
        if (base/'STOP_NEW_CLAIMS').exists() or (base/('DRAIN_'+os.environ['SLURM_JOB_ID'])).exists():break
        if (base/'quarantine'/(socket.gethostname()+'.json')).exists():break
        check_storage(site['local_runtime_root'],site['disk']['local_admit_gib'])
        check_storage('/dev/shm',site['disk']['shm_min_gib'],minimum_inodes=100)
        selected=claim(base,site,species)
        if selected is None:
            idle+=1
            writej(base/f'idle_{socket.gethostname()}_{slot}.json',{'time':time.time(),'reason':'NO_ADMISSIBLE_WORK_OR_SPACE','coverage':coverage(ROOT)})
            if idle>=10:break  # Finite idle parent; successor/monitor reads durable coverage.
            time.sleep(30);continue
        idle=0;path,record=selected
        try:
            reconcile_species(record['species_id'],site)
            result=step(ROOT,record['species_id'],execute=True,orca=site['orca_wrapper'],
                        xtb=site['xtb'],scratch_root=site['shared_scratch'])
        except Exception as exc:
            result={'status':'RECOVERY_REQUIRED','reason':type(exc).__name__+': '+str(exc)}
        finally:
            with queue_lock(base):
                if path.is_file() and readj(path).get('token')==record['token']:
                    archive=base/'claim_history';archive.mkdir(exist_ok=True)
                    path.rename(archive/(path.stem+'_'+record['token']+'.json'))
        writej(base/'latest'/(record['species_id']+'.json'),
               {**record,'returned_epoch':time.time(),'result':result})
        if any(word in result.get('status','') for word in ('REVIEW','RECOVERY_REQUIRED','ERROR','BLOCKED')):
            repair=ROOT/'automation'/record['species_id']/'recovery_status.json'
            findings=readj(repair).get('findings',[]) if repair.is_file() else []
            different=any(v.get('category')=='NODE_READ_ONLY' for v in findings)
            writej(base/'blocked'/(record['species_id']+'.json'), {
                **record,'result':result,'retry_on_different_node':different,
                'notice':'Durable repair request; not silently dropped or endlessly polled. After targeted repair archive this block with operations unblock.'})
        if not result.get('quantum_jobs_launched'):time.sleep(5)


def reconcile_species(sid,site):
    """Recover completed stage-out first, then at most one diagnosed node retry.

    Never rerun a completed quantum job just because its parent died. All other
    problems remain durable repair requests rather than disappearing silently.
    """
    directories=[ROOT/r['job_directory'] for r in readj(ROOT/'registry/jobs.json') if r['species_id']==sid]
    base=ROOT/'automation'/sid
    directories.extend(p.parent for p in base.glob('wave_*/workspace/jobs/*/*/*/*/*/job.json'))
    findings=[]
    for directory in directories:
        for intent in sorted((directory/'stage_intents').glob('*.json')):
            if (directory/'receipt.json').is_file():break
            record=readj(intent)
            if not scheduler_terminal(record):continue
            try:
                result=recover_stageout(ROOT,intent,execute=True)
                findings.append(result)
            except (ValueError,OSError,KeyError) as exc:
                findings.append({'intent':str(intent),'status':'REPAIR_REQUIRED','reason':str(exc)})
        receipt_path=directory/'receipt.json'
        if not receipt_path.is_file():continue
        receipt=readj(receipt_path)
        if receipt.get('status') not in ('failed','FAILED'):continue
        report=diagnosis(directory);findings.append(report)
        if report['category']=='NODE_READ_ONLY' and (directory/'job.json').is_file():
            workspace=next((p for p in directory.parents if (p/'config/protocol.json').is_file()),None)
            if workspace is not None:
                result=retry_dft_infrastructure(workspace,str(directory.relative_to(workspace)),
                    site['orca_wrapper'],site['shared_scratch'],execute=True)
                findings.append({'bounded_node_retry':result})
                break
    writej(base/'recovery_status.json',{'time':time.time(),'findings':findings,
           'notice':'No scientific acceptance is fabricated; unknown failures require targeted repair.'})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execute',action='store_true',required=True)
    p.add_argument('--max-steps',type=int,default=10000)
    a=p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):raise ValueError('Compute allocation required')
    site=site_config();check=wrapper_preflight(site['orca_wrapper'],site['xtb'])
    if not check or not check['passed']:raise ValueError('Current runtime evidence required: '+str(check))
    slots=site['packed']['slots_per_node']
    if slots*site['packed']['slot_memory_mb']>int(os.environ.get('SLURM_MEM_PER_NODE','0')):
        raise ValueError('Packed slot reservations exceed actual node allocation')
    deadline=end_time()
    if (ROOT/'operations/packed/quarantine'/(socket.gethostname()+'.json')).exists():
        raise ValueError('Node quarantined; no new quantum work until explicit health revalidation')
    children=[multiprocessing.Process(target=worker,args=(slot,slots,deadline,a.max_steps)) for slot in range(slots)]
    for child in children:child.start()
    while any(child.is_alive() for child in children):
        states=[{'pid':child.pid,'alive':child.is_alive(),'exitcode':child.exitcode} for child in children]
        writej(ROOT/'operations/packed/workers'/(socket.gethostname()+'_'+os.environ['SLURM_JOB_ID']+'.json'),
               {'updated_epoch':time.time(),'workers':states,'live_slots':sum(r['alive'] for r in states)})
        if any(child.exitcode not in (None,0) for child in children):
            writej(ROOT/'operations/packed/quarantine'/(socket.gethostname()+'.json'),
                   {'reason':'Worker health/runtime failure; inspect before new admissions','workers':states,'time':time.time()})
        for child in children:child.join(timeout=1)
        time.sleep(5)
    coverage(ROOT,write=True)
    # A zero Slurm parent code alone is never reported as completed chemistry.
    if any(child.exitcode for child in children):raise SystemExit(2)


if __name__=='__main__':main()
