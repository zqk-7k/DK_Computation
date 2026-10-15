"""One bounded audit/expansion tick; called by the Vultr systemd timer."""
import argparse
import json
import os
from pathlib import Path
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from clusterflow.common import readj,writej
from clusterflow.priority import ranked_species
from clusterflow.stage_release import audit
from scripts.submit_cohort import submit


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--execute',action='store_true');a=p.parse_args()
    if (ROOT/'operations/goat_bulk_20261015/scope.json').exists():
        from scripts.goat_bulk_status import audit as bulk_audit
        report=bulk_audit(ROOT)
        writej(ROOT/'operations/rollout_latest.json',report)
        print(json.dumps(report,ensure_ascii=False,indent=2))
        return
    policy=readj(ROOT/'config/rollout_20261015.json')
    report={'updated_epoch':time.time(),'stage_release':audit(ROOT),'submissions':[]}
    ranked_species(ROOT,[r['species_id'] for r in readj(ROOT/'registry/batch_plan_v5.json')['species']],write=True)
    stopped=(ROOT/'operations/packed/STOP_NEW_CLAIMS').exists() or (ROOT/'operations/ROLLOUT_PAUSED').exists()
    if stopped: report['status']='PAUSED_BY_OPERATOR'
    else:
        os.environ['WORK2_RUNTIME_EVIDENCE']=str(ROOT/policy['canonical_runtime_evidence'])
        for cohort in policy['cohorts']:
            result=submit(ROOT,cohort['name'],ROOT/cohort['scope'],cohort['nodes'],
                          execute=a.execute,require_goat=cohort['require_valid_goat'])
            report['submissions'].append({'name':cohort['name'],**result})
        report['status']='AUDITED_NOT_SCIENTIFIC_ACCEPTANCE'
    writej(ROOT/'operations/rollout_latest.json',report)
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
