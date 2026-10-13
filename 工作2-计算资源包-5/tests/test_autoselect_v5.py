"""Synthetic source-bound selection fixtures; no quantum software is executed."""
import copy
from pathlib import Path
import sys
import unittest
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import test_cluster_pipeline as fixture
from test_cluster_pipeline import SPECIES, SCENARIO
from clusterflow import autoselect as a, pipeline as p, geometry as g
from clusterflow.common import readj, writej, sha


class AutomaticSelectionTests(unittest.TestCase):
    setUp = fixture.ClusterPipelineTests.setUp
    completed_goat = fixture.ClusterPipelineTests.completed_goat
    completed_opt = fixture.ClusterPipelineTests.completed_opt

    def policy(self):
        writej(self.root / 'config/automation.json', {'automatic_handoff': {'enabled': True}, 'selection': {'rmsd_A': .10}})

    def ensemble(self, count=3, energies=None):
        self.policy()
        self.completed_goat()
        energies = energies or [-99.99, -100.0, -99.998][:count]
        self.mols = []
        texts = []
        for index, energy in enumerate(energies, 1):
            path = ROOT / '01_goat_explicit' / SPECIES / SCENARIO / f'start_{index:02d}' / 'cluster.mol'
            mol = Chem.MolFromMolBlock(path.read_text(encoding='utf-8'), removeHs=False)
            self.assertEqual(g.validate_geometry(self.mol, mol, self.amap, self.cfg), [])
            self.mols.append(mol)
            texts.append(g.xyz_text(mol, f'Energy = {energy:.12f} Eh'))
        (self.goat / 'goat.finalensemble.xyz').write_text(''.join(texts), encoding='utf-8')
        self.rebind()

    def rebind(self):
        receipt = readj(self.goat / 'receipt.json')
        receipt['output_sha256'] = {name: sha(self.goat / name) for name in ('goat.out', 'goat.finalensemble.xyz')}
        writej(self.goat / 'receipt.json', receipt)

    def test_source_without_completed_receipt_is_not_candidate(self):
        self.policy()
        with self.assertRaisesRegex(ValueError, 'receipt'):
            a.screen_batch(self.root, [self.row['job_id']], 3)

    def test_explicit_energies_rank_not_frame_order_and_keep_reserve(self):
        self.ensemble()
        report = a.screen_batch(self.root, [self.row['job_id']], 2)
        self.assertTrue(report['target_reached'])
        self.assertEqual(report['selected'][0]['frame_index_one_based'], 2)
        self.assertEqual(report['selected'][0]['energy_Eh'], -100.0)
        self.assertEqual(len(report['reserve']), 1)
        self.assertFalse(report['human_review_completed'])
        self.assertFalse(report['sampling_convergence_established'])

    def test_real_completed_batch_required_even_if_enough_first_source_candidates(self):
        self.ensemble()
        receipt = readj(self.goat / 'receipt.json')
        receipt['status'] = 'RUNNING'
        writej(self.goat / 'receipt.json', receipt)
        with self.assertRaisesRegex(ValueError, 'not completed'):
            a.screen_batch(self.root, [self.row['job_id']], 1)

    def test_explicit_energy_ambiguity_missing_and_nonfinite_rejected(self):
        frame = g.xyz_text(self.mol, 'Energy = -123 Eh -124 Eh')
        for comment in ('Energy = -123 Eh -124 Eh', 'Energy unknown', 'NaN Eh', 'inf Eh', '-123'):
            with self.subTest(comment=comment), self.assertRaises(ValueError):
                a.parse_frame_energies(g.parse_xyz_frames(g.xyz_text(self.mol, comment)), '')

    def test_native_zero_based_official_table_binds_bare_energy_units(self):
        energies = [-100.0, -99.999]
        frames = g.parse_xyz_frames(''.join(g.xyz_text(self.mol, str(value)) for value in energies))
        output = '# Final ensemble info #\n Conformer Energy Degen. % total % cumul.\n (kcal/mol)\n ---\n 0 0.000 1 70.0 70.0\n 1 0.628 1 30.0 100.0\n Lowest energy conformer : -100.000000 Eh\n'
        self.assertEqual(a.parse_frame_energies(frames, output), energies)
        with self.assertRaisesRegex(ValueError, 'disagrees'):
            a.parse_frame_energies(frames, output.replace('0.628', '2.000'))
        with self.assertRaisesRegex(ValueError, 'Incomplete'):
            a.parse_frame_energies(frames, output.replace(' 1 0.628 1 30.0 100.0\n', ''))

    def test_whole_solvent_permutation_deduplicated(self):
        self.ensemble(count=1)
        moved = Chem.Mol(self.mols[0])
        components = {}
        for row in self.amap:
            components.setdefault(row['component'], []).append(row['index_zero_based'])
        first, second = [ids for name, ids in components.items() if name.startswith('DMSO')]
        xyz = self.mols[0].GetConformer().GetPositions()
        for i, j in zip(first, second):
            moved.GetConformer().SetAtomPosition(i, xyz[j])
            moved.GetConformer().SetAtomPosition(j, xyz[i])
        (self.goat / 'goat.finalensemble.xyz').write_text(g.xyz_text(self.mols[0], '-100 Eh') + g.xyz_text(moved, '-99.999 Eh'), encoding='utf-8')
        self.rebind()
        report = a.screen_batch(self.root, [self.row['job_id']], 2)
        self.assertEqual(report['available_unique_count'], 1)
        self.assertEqual(len(report['duplicates']), 1)
        self.assertFalse(report['target_reached'])
        with self.assertRaisesRegex(ValueError, 'target not reached'):
            a.import_selection(self.root, report)
        self.assertFalse((self.root / 'imports').exists())

    def test_invalid_dissociated_geometry_excluded_not_counted(self):
        self.ensemble(count=1)
        moved = Chem.Mol(self.mols[0])
        for row in self.amap:
            if row['component'].startswith('H2O'):
                i = row['index_zero_based']
                pos = moved.GetConformer().GetAtomPosition(i)
                moved.GetConformer().SetAtomPosition(i, (pos.x + 100, pos.y, pos.z))
        (self.goat / 'goat.finalensemble.xyz').write_text(g.xyz_text(moved, '-100 Eh'), encoding='utf-8')
        self.rebind()
        report = a.screen_batch(self.root, [self.row['job_id']], 1)
        self.assertFalse(report['target_reached'])
        self.assertEqual(len(report['excluded']), 1)
        self.assertEqual(report['selected'], [])

    def test_automatic_import_provenance_and_pipeline_revalidation(self):
        self.ensemble()
        report = a.screen_batch(self.root, [self.row['job_id']], 2)
        result = a.import_selection(self.root, report)
        self.assertEqual(result['imported'], 2)
        index = p._validated_imports(self.root, SPECIES, SCENARIO)
        review = readj(self.root / index['sources'][0]['review'])
        self.assertEqual(review['status'], a.STATUS)
        self.assertNotEqual(review['status'], 'APPROVED_FOR_PRESCREEN')
        self.assertFalse(review['human_review_completed'])
        prepared = p.prepare_prescreen(self.root, SPECIES, SCENARIO)
        self.assertEqual(prepared['job_count'], 2)
        self.assertEqual(a.import_selection(self.root, report)['status'], 'AUTOMATIC_SELECTION_ALREADY_IMPORTED')

    def test_tampered_report_has_no_partial_import(self):
        self.ensemble()
        report = a.screen_batch(self.root, [self.row['job_id']], 2)
        report['selected'][0]['energy_Eh'] = -999
        with self.assertRaisesRegex(ValueError, 'tampered'):
            a.import_selection(self.root, report)
        self.assertFalse((self.root / 'imports').exists())

    def test_altered_policy_invalidates_prepared_import(self):
        self.ensemble()
        report = a.screen_batch(self.root, [self.row['job_id']], 2)
        a.import_selection(self.root, report)
        writej(self.root / 'config/automation.json', {'automatic_handoff': {'enabled': True}, 'selection': {'rmsd_A': .2}})
        with self.assertRaisesRegex(ValueError, 'policy changed'):
            p._validated_imports(self.root, SPECIES, SCENARIO)

    def test_previous_selection_preserved_or_replaced_by_duplicate_alias(self):
        self.ensemble()
        first = a.screen_batch(self.root, [self.row['job_id']], 1)
        second = a.screen_batch(self.root, [self.row['job_id']], 3, previous_selected_ids=[first['selected'][0]['conformer_id']])
        self.assertEqual(first['selected'][0]['conformer_id'], second['selected'][0]['conformer_id'])
        self.assertEqual(second['available_unique_count'], 3)

    def test_automatic_import_disabled_by_default(self):
        self.ensemble()
        writej(self.root / 'config/automation.json', {'automatic_handoff': {'enabled': False}})
        with self.assertRaisesRegex(ValueError, 'not explicitly enabled'):
            a.screen_batch(self.root, [self.row['job_id']], 1)

    def test_completed_receipt_with_active_lock_rejected(self):
        self.ensemble(count=1)
        (self.goat / '.goat-run.lock').write_text('unit-test-only', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'active runner lock'):
            a.screen_batch(self.root, [self.row['job_id']], 1)

    def test_extreme_high_energy_is_reserve_not_quota_filler(self):
        self.ensemble(energies=[-100, -99, -98])
        report = a.screen_batch(self.root, [self.row['job_id']], 3)
        self.assertEqual(report['energy_eligible_count'], 1)
        self.assertEqual(len(report['selected']), 1)
        self.assertEqual(len(report['reserve']), 2)
        self.assertFalse(report['target_reached'])

    def final_fixture(self):
        self.ensemble(energies=[-100, -99.999, -99.998])
        a.import_selection(self.root, a.screen_batch(self.root, [self.row['job_id']], 3))
        prepared = p.prepare_prescreen(self.root, SPECIES, SCENARIO)
        for relative in prepared['jobs']:
            self.completed_opt(self.root / relative)
        second = p.advance(self.root, SPECIES, SCENARIO, 'prescreen')
        for relative in second['prepared_jobs']:
            self.completed_opt(self.root / relative)
        return second['prepared_jobs']

    def test_final_dft_count_requires_actual_frequency_checked_minima(self):
        directories = self.final_fixture()
        result = p.final_stable_candidates(self.root, SPECIES, SCENARIO)
        self.assertEqual(result['eligible_count'], 3)
        self.assertTrue(result['all_prepared_optfreq_pass'])
        self.assertFalse(result['sampling_convergence_established'])
        self.completed_opt(self.root / directories[0], negative=True)
        changed = p.final_stable_candidates(self.root, SPECIES, SCENARIO)
        self.assertEqual(changed['eligible_count'], 2)
        self.assertTrue(changed['issues'])
        self.assertEqual(changed['human_reviewed_approximations'], [])
        self.assertFalse(changed['all_prepared_optfreq_pass'])

    def test_final_energy_window_configuration_controls_eligible_count(self):
        directories = self.final_fixture()
        for index, relative in enumerate(directories):
            directory = self.root / relative
            output = (directory / 'job.out').read_text(encoding='utf-8')
            (directory / 'job.out').write_text(output.replace('-123.456 Eh', f'{-123.456 + index * .001:.6f} Eh'), encoding='utf-8')
            receipt = readj(directory / 'receipt.json')
            receipt['output_sha256']['job.out'] = sha(directory / 'job.out')
            writej(directory / 'receipt.json', receipt)
        policy = readj(self.root / 'config/automation.json')
        policy['stable_minima'] = {'energy_window_kcal_mol': .8}
        writej(self.root / 'config/automation.json', policy)
        result = p.final_stable_candidates(self.root, SPECIES, SCENARIO)
        self.assertEqual(result['eligible_count'], 2)
        self.assertEqual(result['unique_stable_count'], 3)
        self.assertEqual(result['maximum_relative_G_kcal_mol'], .8)

    def test_public_advance_cannot_bypass_final_quota(self):
        self.final_fixture()
        writej(self.root / 'registry/batch_plan_v5.json', {'species': [{'species_id': SPECIES, 'scenario_id': SCENARIO, 'stable_conformer_target': 5}]})
        result = p.advance(self.root, SPECIES, SCENARIO, 'optfreq')
        self.assertEqual(result['status'], 'BLOCKED_FINAL_STABLE_TARGET')
        self.assertEqual(result['eligible_stable_count'], 3)
        self.assertEqual(result['prepared_jobs'], [])
        self.assertFalse(list((self.root / 'jobs').glob('*/' + SPECIES + '/' + SCENARIO + '/shielding/*')))

    def test_v5_same_geometry_different_energy_never_inflates_unique_count(self):
        p5 = copy.deepcopy(self.cfg)
        p5.setdefault('conformer_policy', {})['final_stable_target_plan'] = 'registry/batch_plan_v5.json'
        rows = [{'mol': self.mol, 'job': {'conformer_id': 'a'}, 'parsed': {'gibbs_free_energy_Eh': -100}},
                {'mol': self.mol, 'job': {'conformer_id': 'b'}, 'parsed': {'gibbs_free_energy_Eh': -99.999}}]
        unique, duplicates = p.deduplicate(rows, p5, self.amap)
        self.assertEqual(len(unique), 1)
        self.assertTrue(duplicates[0]['energy_review_recommended'])
        self.assertFalse(duplicates[0]['large_energy_review_required'])
        rows[1]['parsed']['gibbs_free_energy_Eh'] = -99.99
        unique, duplicates = p.deduplicate(rows, p5, self.amap)
        self.assertEqual(len(unique), 1)
        self.assertTrue(duplicates[0]['large_energy_review_required'])


if __name__ == '__main__':
    unittest.main()
