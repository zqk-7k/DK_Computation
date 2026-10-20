"""Validate a resource-only GOAT delegation before using its scientific output."""
from pathlib import Path
import re
from .common import readj, sha, inside

FILES = ('goat.inp','start.xyz','cluster.mol','atom_map.json','isotope_map.json','provenance.json')


def resolve(root, row):
    root=Path(root).resolve()
    original=inside(root,row['job_directory'])
    record=readj(original/'delegated_goat.json')
    if record.get('kind')!='RESOURCE_COHORT_DELEGATION' or record.get('job_id')!=row['job_id']:
        raise ValueError('Invalid GOAT delegation identity')
    if record.get('canonical_manifest_sha256')!=sha(root/'registry/jobs.json'):
        raise ValueError('Delegation canonical manifest changed')
    cohort=Path(record['cohort']).resolve()
    if not cohort.is_relative_to(root/'operations') or cohort==root:
        raise ValueError('Delegated cohort escapes production operations')
    if sha(cohort/'registry/jobs.json')!=record.get('cohort_manifest_sha256'):
        raise ValueError('Delegated manifest changed')
    rows=[r for r in readj(cohort/'registry/jobs.json') if r.get('job_id')==row['job_id']]
    if len(rows)!=1:raise ValueError('Missing or duplicate delegated job')
    delegated=rows[0]
    if {k:v for k,v in row.items() if k!='input_sha256'} != {k:v for k,v in delegated.items() if k!='input_sha256'}:
        raise ValueError('Delegation changed chemical identity or release scope')
    target=inside(cohort,delegated['job_directory'])
    for name in FILES:
        if sha(original/name)!=record['canonical_inputs'].get(name):
            raise ValueError('Delegation source inputs changed')
        if name!='goat.inp' and sha(original/name)!=sha(target/name):
            raise ValueError('Delegation changed scientific input: '+name)
    original_input=(original/'goat.inp').read_text()
    if re.findall(r'(?i)\bnprocs\s+(\d+)',original_input)!=['8']:
        raise ValueError('Unexpected original GOAT resource layout')
    expected=re.sub(r'(?i)(\bnprocs\s+)8\b',r'\g<1>16',original_input,count=1)
    if record.get('cpus')!=16 or (target/'goat.inp').read_text()!=expected:
        raise ValueError('Delegation must change only nprocs 8 to 16')
    if sha(target/'goat.inp')!=delegated['input_sha256']:
        raise ValueError('Delegated input hash mismatch')
    lock=original/'.goat-run.lock'
    if lock.exists() and readj(lock)!=record:
        raise ValueError('A different canonical GOAT owner exists')
    if (target/'delegated_goat.json').exists():
        raise ValueError('Chained GOAT delegation not permitted')
    return cohort, record
