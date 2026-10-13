"""Per-product choices requested by the sample owner; not method validation."""
from .common import config, readj

def compound_for_species(root,species_id):
    rows=[s for s in readj(root/'registry/species.json') if s['species_id']==species_id]
    if len(rows)!=1:raise ValueError('Unknown or duplicate species identifier')
    return rows[0]['compound_id']

def j_methods_for(settings,compound_id):
    route=settings.get('j_routing',{})
    if compound_id in route.get('figure4_compounds',[]):
        return [route.get('figure4_method','PBE_pcJ3')]
    return route.get('other_candidates',settings.get('j_methods',['PBE_pcJ3','PBE0_pcJ2']))

def resolve_j_method(root,species_id,requested=None):
    settings=config(root)
    if not settings.get('j_routing'):
        return requested or 'PBE0_pcJ2'
    compound=compound_for_species(root,species_id)
    methods=j_methods_for(settings,compound)
    selected=requested or (methods[0] if len(methods)==1 else settings['j_routing']['other_exploratory_default'])
    if selected not in methods:raise ValueError('J method conflicts with declared compound routing: '+compound)
    return selected
