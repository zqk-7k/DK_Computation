"""Archive superseded v2 notes in this NEW copy, never touch the original v2."""
from pathlib import Path
import json
import shutil

ROOT=Path(__file__).resolve().parents[1]
assert ROOT.name=='工作2-计算资源包-（3）'
NAMES=['README_先读.txt','docs/流程差异与未完成事项.txt',
       'config/stage_methods_51_aligned.json','scripts/finish_delivery.py',
       'reports/delivery.json','reports/delivery_file_manifest.json',
       'reports/unit_tests.txt','reports/local_preflight.json','reports/independent_audit.json']

for name in NAMES:
    source=(ROOT/name).resolve();target=(ROOT/'history/v2'/name).resolve()
    if not source.is_relative_to(ROOT) or not target.is_relative_to(ROOT):raise ValueError('Escaping archive path')
    if not source.exists():continue
    if target.exists():raise ValueError('History target already exists')
    target.parent.mkdir(parents=True,exist_ok=True)
    source.rename(target)

manifest={'package_name':ROOT.name,'sequence_number':3,
          'replaces_for_current_work':'ORCA_64_18O_EXPLICIT_51REF_v2',
          'original_v1_v2_untouched':True,'shielding_method':'PBE0/pcSseg-3 only',
          'active_protocol':'config/protocol.json','original_reference':'reference_51/',
          'reference_warning':'Historical reference may contain excluded methods. It is not the active execution configuration.',
          'initial_cluster_data':'Same938 validated inputs as v2; only downstream code upgraded.'}
(ROOT/'PACKAGE.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(manifest,ensure_ascii=False,indent=2))
