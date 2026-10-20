"""Prepare and submit an explicitly commissioned, source-bound 16-CPU cohort.

Canonical inputs and receipts are never overwritten. Node-failure retries keep
all previous scratch; a durable delegation and common execution lock prevent
the original pool from duplicating them.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from clusterflow.common import readj, writej, sha, inside
from scripts.packed_worker import queue_lock
from scripts.goat_job import INPUT_HASH_FIELDS, validate_job
from scripts.goat_portable import IMMUTABLE_FILES
from scripts.submit_goat_bulk import IDENTITIES
from clusterflow.scheduling import single_execution

RELATIVE = 'operations/expansion_20261020'


def resource_input(text):
    matches = re.findall(r'(?i)\bnprocs\s+(\d+)', text)
    if matches != ['8']:
        raise ValueError('Only an exact 8 to 16 process override is allowed')
    return re.sub(r'(?i)(\bnprocs\s+)8\b', r'\g<1>16', text, count=1)


def accounting():
    result = subprocess.run(
        ['sacct', '-X', '-n', '-P', '-S', '2026-10-15', '-j', '124489624',
         '-o', 'JobIDRaw,JobID,State,ExitCode,End,NodeList'],
        text=True, capture_output=True, check=True, timeout=60)
    return [line.split('|') for line in result.stdout.splitlines() if line.strip()]


def prepare(root):
    root = Path(root).resolve()
    base = root/'operations/goat_bulk_20261015'
    destination = root/RELATIVE
    if destination.exists():
        raise ValueError('Cohort exists: reconcile its ledger; never repeat preparation')
    source_rows = {r['job_id']: r for r in readj(root/'registry/jobs.json')}
    original_scope = readj(base/'scope.json')
    scheduler = accounting()
    dead = [r for r in scheduler if r[0] == '124489626'
            and r[1] == '124489624_1' and r[2] == 'NODE_FAIL']
    if len(dead) != 1:
        raise ValueError('Exact failed allocation identity not confirmed')
    excluded = ['a10r2n13'] + ['k18r2n'+str(i) for i in range(10,18)]
    with queue_lock(root/'operations/submission_lock'), queue_lock(base):
        queue = readj(base/'queue.json')
        pending = []
        for job in original_scope['job_ids'][queue['cursor']:]:
            row = source_rows[job]
            directory = inside(root, row['job_directory'])
            if any((directory/n).exists() for n in
                   ('.goat-run.lock','receipt.json','portable_stageout.json','delegated_goat.json')):
                continue
            if any((directory/'stage_intents').glob('*.json')):
                continue
            validate_job(root, job)
            pending.append(job)
        interrupted = []
        for job, claim in queue['active'].items():
            if str(claim.get('slurm_job_id')) != '124489626':
                continue
            directory = inside(root, source_rows[job]['job_directory'])
            owner = readj(directory/'.goat-run.lock')
            if str(owner.get('slurm_job_id')) != '124489626':
                raise ValueError('Interrupted owner mismatch')
            stage = Path(owner['portable_stage_root'])/'job'
            receipt = readj(stage/'receipt.json')
            if receipt.get('status') != 'RUNNING':
                raise ValueError('Recover terminal result before considering a retry')
            interrupted.append(job)
        if len(pending) < 128 or len(interrupted) != 28:
            raise ValueError(f'Review changed scope: {len(pending)} fresh, {len(interrupted)} interrupted')
        selected = interrupted + pending
        destination.mkdir(parents=True)
        writej(destination/'PREPARATION_INTENT.json',
               dict(created_epoch=time.time(), canonical_root=str(root),
                    fresh=pending, interrupted=interrupted, scheduler=scheduler,
                    source_identities={n:sha(root/n) for n in IDENTITIES},
                    authorization='User requested another 2000 cores and recovery on 2026-10-20',
                    notice='Node-failure work restarts from the original seed; not an exact GOAT checkpoint resume.'))
        for name in ('scripts','clusterflow','config'):
            shutil.copytree(root/name, destination/name,
                            ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        (destination/'registry').mkdir()
        for name in ('batch_plan_v5.json','species.json'):
            shutil.copy2(root/'registry'/name, destination/'registry'/name)
        rows = []
        for job in selected:
            old = source_rows[job]
            source = inside(root, old['job_directory'])
            target = destination/old['job_directory']
            target.mkdir(parents=True)
            bindings = {name:sha(source/name) for name in IMMUTABLE_FILES}
            for name, field in INPUT_HASH_FIELDS.items():
                if bindings[name] != old[field].lower():
                    raise ValueError('Canonical input identity changed: '+job)
            if bindings['provenance.json'] != old['provenance_sha256'].lower():
                raise ValueError('Canonical provenance identity changed')
            for name in IMMUTABLE_FILES:
                shutil.copy2(source/name, target/name)
            (target/'goat.inp').write_text(resource_input((source/'goat.inp').read_text()))
            row = dict(old, input_sha256=sha(target/'goat.inp'))
            rows.append(row)
            writej(target/'resource_override.json',
                   dict(canonical_root=str(root), canonical_job_id=job,
                        canonical_manifest_sha256=sha(root/'registry/jobs.json'),
                        canonical_inputs=bindings, cpus_before=8, cpus_after=16,
                        change='nprocs only; same geometry, seed, model and method'))
        writej(destination/'registry/jobs.json', rows)
        for row in rows:
            validate_job(destination, row['job_id'])
        # Fully stage and validate the cohort before taking canonical ownership.
        for job in selected:
            source = inside(root, source_rows[job]['job_directory'])
            digest = 'goat_' + hashlib.sha256(job.encode()).hexdigest()
            with single_execution(root, digest):
                token = uuid.uuid4().hex
                if job in interrupted:
                    lock = source/'.goat-run.lock'
                    owner = readj(lock)
                    if str(owner.get('slurm_job_id')) != '124489626':
                        raise ValueError('Interrupted owner changed')
                    archive = source/'recovery_history'/('node_fail_20261020_'+token)
                    archive.mkdir(parents=True)
                    lock.rename(archive/'original_goat_lock.json')
                    writej(archive/'claim.json', queue['active'][job])
                    writej(base/'results'/(job+'.json'),
                           dict(queue['active'][job], finished_epoch=time.time(),
                                result=dict(status='NODE_FAIL_DELEGATED_TO_RESOURCE_COHORT',
                                            cohort=str(destination), old_scratch_retained=True)))
                    del queue['active'][job]
                    claim_path=base/'claims'/(job+'.json')
                    if claim_path.exists():
                        claim_path.rename(archive/'pool_claim.json')
                elif (source/'.goat-run.lock').exists():
                    raise ValueError('Fresh job acquired an owner during preparation')
                record = dict(kind='RESOURCE_COHORT_DELEGATION', token=token,
                              canonical_manifest_sha256=sha(root/'registry/jobs.json'),
                              canonical_inputs={n:sha(source/n) for n in IMMUTABLE_FILES},
                              cohort=str(destination), cohort_manifest_sha256=sha(destination/'registry/jobs.json'),
                              job_id=job, cpus=16, created_epoch=time.time(),
                              reason='NODE_FAIL_RETRY' if job in interrupted else 'FIRST_EXECUTION',
                              retry_number=1 if job in interrupted else 0)
                with (source/'.goat-run.lock').open('x') as stream:
                    json.dump(record,stream)
                with (source/'delegated_goat.json').open('x') as stream:
                    json.dump(record,stream,ensure_ascii=False,indent=2)
        writej(base/'queue.json',queue)
        # Pilot uses two fresh production starts; it does not duplicate a test seed.
        pilot = pending[:2]
        main = [job for job in selected if job not in pilot]
        identities = {name:sha(destination/name) for name in IDENTITIES}
        for label, ids, directory in (
                ('pilot',pilot,'operations/pilot_20261020'),
                ('main',main,'operations/goat_bulk_20261015')):
            scope = dict(schema=1,stage='goat',job_ids=ids,operations_directory=directory,
                         cpus_per_job=16,slots_per_node=2,growth_gib_per_job=8,free_floor_gib=128,
                         legacy_active_jobs=300,identity_sha256=identities,
                         created_epoch=time.time(),excluded_nodes=excluded,
                         canonical_root=str(root),notice='Released initial starts only; resource-only override; retained interrupted attempts.')
            writej(destination/directory/'scope.json',scope)
        (destination/'operations/packed').mkdir(parents=True,exist_ok=True)
        (destination/'operations/packed/STOP_NEW_CLAIMS').touch()
        (destination/'operations/ROLLOUT_PAUSED').touch()
        runtime=root/'operations/runtime_124440267.json'
        report=dict(status='PREPARED_NOT_SUBMITTED',cohort=str(destination),
                    fresh=len(pending),interrupted=len(interrupted),pilot_ids=pilot,
                    main_tasks=len(main),runtime_evidence=str(runtime),
                    excluded_nodes=excluded,original_inputs_unchanged=True)
        writej(destination/'PREPARED.json',report)
        return report


def submit(root, label, execute=False):
    root=Path(root).resolve()
    cohort=root/RELATIVE
    preparation=readj(cohort/'PREPARED.json')
    directory='operations/pilot_20261020' if label=='pilot' else 'operations/goat_bulk_20261015'
    base=cohort/directory
    scope=base/'scope.json'
    with queue_lock(root/'operations/submission_lock'):
        journal=base/'submission.json'
        if journal.exists():
            return readj(journal)
        if label=='main':
            pilot=readj(cohort/'PILOT_ACCEPTED.json')
            if pilot.get('status')!='REAL_GOAT_16CPU_TWO_SLOT_START_VERIFIED':
                raise ValueError('Real pilot not accepted')
        site=readj(cohort/'config/site_scnet.json')
        count=1 if label=='pilot' else 64
        result=subprocess.run(['squeue','-r','-A',site['slurm']['account'],'-h','-o','%T|%C'],
                              capture_output=True,text=True,check=True,timeout=60)
        allocated=sum(int(r.split('|')[1]) for r in result.stdout.splitlines()
                      if r.split('|')[0] in ('RUNNING','COMPLETING','CONFIGURING'))
        if allocated+count*32>10000 or len(result.stdout.splitlines())+count>=190:
            raise ValueError('Account or queue headroom insufficient')
        if shutil.disk_usage(site['shared_scratch']).free/2**30 < 4096:
            raise ValueError('Need 4096 GiB free before expansion')
        command=['sbatch','--parsable','--account='+site['slurm']['account'],
                 '--partition='+site['slurm']['partition'],'--qos='+site['slurm']['qos'],
                 '--time=720:00:00','--time-min=720:00:00',
                 '--ntasks-per-node=1','--cpus-per-task=32','--mem=114000M',
                 '--nodes=1',
                 '--job-name=work2-64-'+('pilot16' if label=='pilot' else 'add2048'),
                 '--chdir='+str(cohort),'--exclude='+','.join(preparation['excluded_nodes']),
                 '--output='+str(base/'slurm-%A_%a.out'),'--error='+str(base/'slurm-%A_%a.err')]
        if label=='main':command.append('--array=0-63%64')
        launcher=cohort/'scripts/scnet_goat_expansion.sbatch'
        test=subprocess.run(command+['--test-only',str(launcher)],capture_output=True,text=True,check=True,timeout=60)
        report=dict(status='TEST_ONLY',nodes=count,cpus=count*32,label=label,created_epoch=time.time(),
                    command=command,preflight=test.stderr,scope_sha256=sha(scope),
                    account_allocated_cpus_before=allocated)
        writej(base/'preflight.json',report)
        if not execute:return report
        env=dict(os.environ,WORK2_PACKAGE_ROOT=str(cohort),WORK2_PYTHON=site['python'],
                 WORK2_PACKED_PROFILE='large_memory',WORK2_RUNTIME_EVIDENCE=preparation['runtime_evidence'],
                 WORK2_GOAT_BULK_BASE=str(base),WORK2_GOAT_BULK_SCOPE=str(scope),
                 WORK2_GOAT_BULK_SCOPE_SHA256=sha(scope))
        writej(journal,dict(report,status='SUBMISSION_INTENT'))
        try:
            launch=subprocess.run(command+[str(launcher)],env=env,capture_output=True,text=True,check=True,timeout=60)
            job=launch.stdout.strip().split(';')[0]
            if not job.isdigit():raise ValueError('Unrecognized sbatch result')
        except Exception as exc:
            writej(journal,dict(report,status='SUBMISSION_UNCERTAIN',error=str(exc)))
            raise
        report.update(status='SUBMITTED_NOT_YET_VERIFIED_RUNNING',job_id=job)
        writej(journal,report)
        return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('command',choices=('prepare','pilot','main'))
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    if args.command=='prepare':
        if not args.execute:raise ValueError('Preparation reserves canonical task ownership; use --execute')
        result=prepare(ROOT)
    else:result=submit(ROOT,args.command,args.execute)
    print(json.dumps(result,ensure_ascii=False,indent=2))
