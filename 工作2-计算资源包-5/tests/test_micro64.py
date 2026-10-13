import copy
import importlib.util
import json
from pathlib import Path
import sys
import unittest
import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import micro64 as mc
import mass_handoff as mh

class MicroTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.jobs = mc.readj(ROOT/'registry/jobs.json')
        cls.primary = next(j for j in cls.jobs if j['compound_id']=='F2_E01' and j['DMSO_d6_count']==2 and j['H2O_count']==1)
        cls.directory = ROOT/cls.primary['job_directory']

    def test_64_products(self):
        self.assertEqual(len({j['compound_id'] for j in self.jobs}-{'TMS'}),64)

    def test_full_first_wave_counts(self):
        plan=mc.readj(ROOT/'registry/sampling_plan_v4.json')
        self.assertEqual(len(self.jobs),plan['summary']['total_starts'])
        self.assertEqual(len({(j['species_id'],j['scenario_id']) for j in self.jobs}),len(plan['species']))

    def test_cluster18o_and_D(self):
        m,a,x = mh.load_cluster(self.directory, self.directory/'start.xyz')
        self.assertEqual(sum(r['isotope']==18 for r in a),1)
        self.assertEqual(sum(r['isotope']==2 for r in a),12)
        self.assertEqual(len(a)-self.primary['solute_atom_count'],23)

    def test_mass_block(self):
        m,a,x = mh.load_cluster(self.directory,self.directory/'start.xyz')
        block=mh.coordinate_block(m,a,x)
        self.assertEqual(block.count('M=17.99915961286'),1)
        self.assertEqual(block.count('M=2.01410177812'),12)

    def test_corrupt_mass_rejected(self):
        m,a,x=mh.load_cluster(self.directory,self.directory/'start.xyz')
        next(r for r in a if r['isotope']==18)['isotope_mass_u']=16.0
        with self.assertRaises(ValueError):mh.coordinate_block(m,a,x)

    def test_charge_rejected(self):
        m,a,x=mh.load_cluster(self.directory,self.directory/'start.xyz')
        with self.assertRaises(ValueError):mh.coordinate_block(m,a,x,1)

    def test_solvent_excluded(self):
        m,a,x=mh.load_cluster(self.directory,self.directory/'start.xyz')
        masks=mh.analyte_masks(a)
        limit=self.primary['solute_atom_count']
        self.assertTrue(all(i<limit for v in masks['zero_based_analyte_nuclei'].values() for i in v))
        self.assertTrue(all(i<limit and j<limit for i,j in masks['HH_HF_CF_pairs']))

    def test_solvent_included_rejected(self):
        m,a,x=mh.load_cluster(self.directory,self.directory/'start.xyz')
        next(r for r in a if r['environment_atom'])['include_in_analyte_nmr']=True
        with self.assertRaises(ValueError):mh.analyte_masks(a)

    def test_exact_plain_xyz(self):
        lines=(self.directory/'start.xyz').read_text().splitlines()
        self.assertTrue(all(len(s.split())==4 for s in lines[2:]))
        self.assertNotIn('M=', '\n'.join(lines))

    def test_water_not_deuterated(self):
        a=mc.readj(self.directory/'atom_map.json')
        water=[r for r in a if r['component'].startswith('H2O')]
        self.assertEqual(len(water),3)
        self.assertTrue(all(r['isotope']==0 for r in water))

    def test_scenarios_not_pooled(self):
        scenarios=mc.readj(ROOT/'registry/scenarios.json')
        self.assertEqual(len({s['pooling_key'] for s in scenarios}),len(mc.readj(ROOT/'registry/sampling_plan_v4.json')['species']))
        self.assertFalse(mc.readj(ROOT/'config/microsolvation_policy.json')['cross_composition_Boltzmann_pooling'])

    def test_no_optional_model_silent_activation(self):
        self.assertTrue(all(j['charge']==0 and 'Methanol' not in j['species_id'] for j in self.jobs))
        self.assertTrue(mc.readj(ROOT/'registry/held_species.json'))

    def test_source_atom_order(self):
        m,a,x=mh.load_cluster(self.directory,self.directory/'start.xyz')
        src=Chem.MolFromMolBlock((ROOT/'solute_models'/self.primary['species_id']/'model.mol').read_text(encoding='utf-8'),removeHs=False)
        self.assertEqual([(r['element'],r['isotope']) for r in a[:src.GetNumAtoms()]],
                         [(at.GetSymbol(),at.GetIsotope()) for at in src.GetAtoms()])

    def test_only_declared_initial_model(self):
        report=mc.readj(ROOT/'reports/preparation.json')
        self.assertFalse(report['sampling_converged'])
        self.assertFalse(report['DFT_run'])
        self.assertFalse(report['ORCA_GOAT_run'])

if __name__=='__main__':unittest.main()
