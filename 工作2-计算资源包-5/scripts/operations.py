"""Read-only by default; explicit audit report/recovery interfaces."""
import argparse
import json
from pathlib import Path
import sys
import uuid
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from clusterflow.operations import coverage, diagnosis, retry_dft_infrastructure, recover_stageout,cleanup_plan

def main():
    p=argparse.ArgumentParser(description=__doc__)
    s=p.add_subparsers(dest='command',required=True)
    a=s.add_parser('coverage');a.add_argument('--write',action='store_true')
    a=s.add_parser('cleanup-plan');a.add_argument('--write',action='store_true')
    a=s.add_parser('diagnose');a.add_argument('--job',required=True,type=Path)
    a=s.add_parser('retry-node-failure');a.add_argument('--job',required=True)
    a.add_argument('--orca',required=True);a.add_argument('--scratch-root',required=True)
    a.add_argument('--execute',action='store_true')
    a=s.add_parser('recover-stageout');a.add_argument('--intent',required=True,type=Path)
    a.add_argument('--execute',action='store_true')
    a=s.add_parser('unblock');a.add_argument('--species',required=True)
    a.add_argument('--reason',required=True);a.add_argument('--execute',action='store_true')
    a=p.parse_args()
    if a.command=='coverage':r=coverage(ROOT,write=a.write)
    elif a.command=='cleanup-plan':r=cleanup_plan(ROOT,write=a.write)
    elif a.command=='diagnose':r=diagnosis(a.job)
    elif a.command=='recover-stageout':r=recover_stageout(ROOT,a.intent,execute=a.execute)
    elif a.command=='unblock':
        from clusterflow.common import readj,writej
        ids={v['species_id'] for v in readj(ROOT/'registry/batch_plan_v5.json')['species']}
        if a.species not in ids:raise ValueError('Unknown registered species')
        path=ROOT/'operations/packed/blocked'/(a.species+'.json')
        r={'status':'UNBLOCK_PLAN_NOT_SCIENTIFIC_APPROVAL','species':a.species,'reason':a.reason}
        if a.execute and path.exists():
            archive=ROOT/'operations/packed/block_history';archive.mkdir(exist_ok=True)
            target=archive/(a.species+'_'+uuid.uuid4().hex+'.json');path.rename(target)
            writej(target.with_suffix('.decision.json'),r)
            r['status']='REQUEUE_ALLOWED_EXISTING_SCIENTIFIC_GATES_UNCHANGED'
    else:r=retry_dft_infrastructure(ROOT,a.job,a.orca,a.scratch_root,execute=a.execute)
    print(json.dumps(r,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
