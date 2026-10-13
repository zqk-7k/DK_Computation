"""Synthetic, temporary pipeline fixtures only. Never executes ORCA or xTB."""
import copy
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import numpy.testing  # Load Windows platform helpers before the subprocess guard.
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from clusterflow import pipeline as p, geometry as g, inputs, gates
from clusterflow.common import config, hash_protocol, readj, writej, sha

SPECIES = 'F2_E01_rep_R__DMSO'
SCENARIO = 'DMSO_02_W01_trace_water_sensitivity'


class ClusterPipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='clusterflow_test_')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name) / '中文测试包'
        self.root.mkdir()
        self.cfg = readj(ROOT / 'config/protocol.json')
        writej(self.root / 'config/protocol.json', self.cfg)
        species = next(s for s in readj(ROOT / 'registry/species.json') if s['species_id'] == SPECIES)
        writej(self.root / 'registry/species.json', [species])
        source = ROOT / '01_goat_explicit' / SPECIES / SCENARIO / 'start_01'
        self.goat = self.root / '01_goat_explicit' / SPECIES / SCENARIO / 'start_01'
        self.goat.mkdir(parents=True)
        for name in ('goat.inp', 'start.xyz', 'cluster.mol', 'atom_map.json', 'provenance.json'):
            shutil.copyfile(source / name, self.goat / name)
        self.row = copy.deepcopy(next(r for r in readj(ROOT / 'registry/jobs.json')
                                      if r['species_id'] == SPECIES and r['scenario_id'] == SCENARIO and r['start_index'] == 1))
        writej(self.root / 'registry/jobs.json', [self.row])
        self.mol = Chem.MolFromMolBlock((source / 'cluster.mol').read_text(encoding='utf-8'), removeHs=False)
        self.amap = readj(source / 'atom_map.json')
        self.child = patch('subprocess.Popen', side_effect=AssertionError('Quantum/subprocess execution forbidden in pipeline tests'))
        self.child.start()
        self.addCleanup(self.child.stop)

    def completed_goat(self):
        # This is a fixture, not an actual completed scientific calculation.
        moved = Chem.Mol(self.mol)
        conformer = moved.GetConformer()
        for i, xyz in enumerate(conformer.GetPositions()):
            conformer.SetAtomPosition(i, xyz + [.2, -.1, .3])
        (self.goat / 'goat.finalensemble.xyz').write_text(g.xyz_text(moved, 'SYNTHETIC UNIT TEST ONLY'), encoding='utf-8')
        (self.goat / 'goat.out').write_text('SYNTHETIC UNIT TEST ONLY\nProgram Version 6.1.1\nORCA TERMINATED NORMALLY\n', encoding='utf-8')
        receipt = {'job_id': self.row['job_id'], 'status': 'COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW',
                   'input_sha256': {name: sha(self.goat / name) for name in ('goat.inp', 'start.xyz', 'atom_map.json', 'cluster.mol', 'provenance.json')},
                   'output_sha256': {name: sha(self.goat / name) for name in ('goat.out', 'goat.finalensemble.xyz')}}
        writej(self.goat / 'receipt.json', receipt)

    def approved_review(self):
        value = p.draft_goat_review(self.root, self.row['job_id'])
        value.update({'status': 'APPROVED_FOR_PRESCREEN', 'reviewer': 'UNIT_TEST_FIXTURE_NOT_HUMAN_APPROVAL',
                      'reviewed_utc': '2026-01-01T00:00:00Z', 'reason': 'Synthetic unit test only',
                      'selected_frames': [{'frame_index_one_based': 1, 'conformer_id': 'c000'}],
                      'checks': {key: True for key in p.REVIEW_CHECKS}})
        path = self.root / 'review.json'
        writej(path, value)
        return path

    def prepared(self):
        self.completed_goat()
        p.import_goat(self.root, self.row['job_id'], self.approved_review())
        return p.prepare_prescreen(self.root, SPECIES, SCENARIO)

    def completed_opt(self, directory, negative=False):
        directory = Path(directory)
        mol = Chem.MolFromMolBlock((directory / 'input_structure.mol').read_text(encoding='utf-8'), removeHs=False)
        (directory / 'job.xyz').write_text((directory / 'input.xyz').read_text(encoding='utf-8'), encoding='utf-8')
        lines = ['SYNTHETIC UNIT TEST ONLY', 'Program Version 6.1.1', 'THE OPTIMIZATION HAS CONVERGED', 'VIBRATIONAL FREQUENCIES']
        for i in range(3 * mol.GetNumAtoms()):
            value = -20 if negative and i == 6 else (0 if i < 6 else 100)
            lines.append(f'{i}: {value:.2f} cm**-1')
        lines += ['NORMAL MODES', 'THERMOCHEMISTRY AT 298.15 K', 'Quasi RRHO ... True',
                  'reference frequency of 100.0 cm-1', 'Final Gibbs free energy ... -123.456 Eh', 'ORCA TERMINATED NORMALLY']
        (directory / 'job.out').write_text('\n'.join(lines) + '\n', encoding='utf-8')
        writej(directory / 'receipt.json', {'status': 'completed',
               'input_sha256': {name: sha(directory / name) for name in ('job.inp', 'job.json', 'input.xyz', 'input_structure.mol', 'atom_map.json')},
               'output_sha256': {name: sha(directory / name) for name in ('job.out', 'job.xyz')}})

    def test_initial_seed_without_completed_goat_never_promoted(self):
        with self.assertRaisesRegex(ValueError, 'receipt'):
            p.draft_goat_review(self.root, self.row['job_id'])
        with self.assertRaises(FileNotFoundError):
            p.prepare_prescreen(self.root, SPECIES, SCENARIO)
        self.assertFalse((self.root / 'jobs').exists())

    def test_scenario_directly_imports_actual_frame_without_continuum(self):
        result = self.prepared()
        self.assertEqual(result['job_count'], 1)
        directory = self.root / result['jobs'][0]
        actual = g.parse_xyz((directory / 'input.xyz').read_text())[1]
        source = g.parse_xyz((self.goat / 'goat.finalensemble.xyz').read_text())[1]
        seed = g.parse_xyz((self.goat / 'start.xyz').read_text())[1]
        np.testing.assert_allclose(actual, source, rtol=0, atol=1e-8)
        self.assertFalse(np.allclose(actual, seed))
        self.assertFalse(list((self.root / 'jobs').glob('*/' + SPECIES + '/DMSO_00_W00_continuum')))

    def test_draft_review_never_counts_as_approval(self):
        self.completed_goat()
        review = p.draft_goat_review(self.root, self.row['job_id'])
        path = self.root / 'draft.json'
        writej(path, review)
        with self.assertRaisesRegex(ValueError, 'approved'):
            p.import_goat(self.root, self.row['job_id'], path)

    def test_modified_goat_output_or_review_binding_rejected(self):
        self.completed_goat()
        path = self.approved_review()
        review = readj(path)
        review['ensemble_sha256'] = '0' * 64
        writej(path, review)
        with self.assertRaisesRegex(ValueError, 'hashes'):
            p.import_goat(self.root, self.row['job_id'], path)
        (self.goat / 'goat.out').write_text('modified', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'binding'):
            p.draft_goat_review(self.root, self.row['job_id'])

    def test_no_import_after_downstream_prepared(self):
        self.prepared()
        with self.assertRaisesRegex(ValueError, 'before preparing'):
            p.import_goat(self.root, self.row['job_id'], self.root / 'review.json')

    def test_unrun_or_negative_frequency_blocks_advance(self):
        prepared = self.prepared()
        result = p.advance(self.root, SPECIES, SCENARIO, 'prescreen')
        self.assertTrue(result['status'].startswith('BLOCKED'))
        directory = self.root / prepared['jobs'][0]
        self.completed_opt(directory, negative=True)
        result = p.advance(self.root, SPECIES, SCENARIO, 'prescreen')
        self.assertTrue(result['status'].startswith('BLOCKED'))
        self.assertIn('imaginary_frequencies_require_mode_review', result['details'][0]['errors'])

    def test_full_temporary_stage_handoff_and_source_revalidation(self):
        prepared = self.prepared()
        self.completed_opt(self.root / prepared['jobs'][0])
        first = p.advance(self.root, SPECIES, SCENARIO, 'prescreen')
        self.assertEqual(first['status'], 'advanced')
        self.completed_opt(self.root / first['prepared_jobs'][0])
        second = p.advance(self.root, SPECIES, SCENARIO, 'optfreq')
        self.assertEqual(second['status'], 'advanced')
        self.assertEqual(set(second['prepared_jobs'][0]['j_jobs']), {'PBE_pcJ3', 'PBE0_pcJ2'})
        result = p.validate_ensemble_sources(self.root, SPECIES, SCENARIO)
        self.assertTrue(result['accepted'])
        ensemble = result['ensemble']
        self.assertFalse(ensemble['sampling_convergence_established'])
        self.assertEqual(ensemble['conformers'][0]['weight'], 1)
        shielding = self.root / ensemble['conformers'][0]['shielding_job']
        meta = readj(shielding / 'job.json')
        self.assertTrue(meta['method_fingerprint'])
        self.assertTrue(all(not self.amap[i]['environment_atom'] for i in meta['expected_nuclei_zero_based']))
        (shielding / 'input.xyz').write_text((shielding / 'input.xyz').read_text().replace('0.200', '1.200'), encoding='utf-8')
        # Stronger deliberate perturbation than string matching uncertain decimals.
        mol = g.with_xyz(self.mol, (shielding / 'input.xyz').read_text())
        xyz = mol.GetConformer().GetAtomPosition(0)
        mol.GetConformer().SetAtomPosition(0, (xyz.x + 1, xyz.y, xyz.z))
        (shielding / 'input.xyz').write_text(g.xyz_text(mol), encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'coordinates'):
            p.validate_ensemble_sources(self.root, SPECIES, SCENARIO)

    def test_sampling_review_boolean_does_not_establish_convergence(self):
        self.completed_goat()
        path = self.approved_review()
        value = readj(path)
        value['sampling_convergence_established'] = True
        writej(path, value)
        p.import_goat(self.root, self.row['job_id'], path)
        index = readj(p._index_path(self.root, SPECIES, SCENARIO))
        self.assertFalse(index['sampling_convergence_established'])
        self.assertTrue(index['sources'][0]['sampling_review_claim'])

    def test_solvent_component_permutation_is_same_basin(self):
        changed = Chem.Mol(self.mol)
        groups = {}
        for row in self.amap:
            groups.setdefault(row['component'], []).append(row['index_zero_based'])
        first, second = [indices for component, indices in groups.items() if component.startswith('DMSO')]
        original = self.mol.GetConformer().GetPositions()
        for i, j in zip(first, second):
            changed.GetConformer().SetAtomPosition(i, original[j])
            changed.GetConformer().SetAtomPosition(j, original[i])
        self.assertLess(p.cluster_rms(self.mol, changed, self.amap), 1e-6)

    def test_distinct_solvent_motif_is_not_solute_only_duplicate(self):
        changed = Chem.Mol(self.mol)
        first_component = next(r['component'] for r in self.amap if r['environment_atom'])
        for row in self.amap:
            if row['component'] == first_component:
                i = row['index_zero_based']
                pos = changed.GetConformer().GetAtomPosition(i)
                changed.GetConformer().SetAtomPosition(i, (pos.x + 1, pos.y, pos.z))
        self.assertGreater(p.cluster_rms(self.mol, changed, self.amap), .1)


if __name__ == '__main__':
    unittest.main()
