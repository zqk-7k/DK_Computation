"""Controller tests use temporary inputs/mocks only; never execute QM programs."""
import copy
import json
from pathlib import Path
import tempfile
import shutil
import unittest
from unittest.mock import patch

from clusterflow import autoflow as a
from clusterflow.common import readj, writej

ROOT = Path(__file__).resolve().parents[1]
SID = 'F2_E01_rep_R__DMSO'


class AutoflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='af5_')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cfg = readj(ROOT / 'config/protocol.json')
        writej(self.root / 'config/protocol.json', self.cfg)
        writej(self.root / 'config/automation.json', {'automatic_handoff': {'enabled': True}})
        self.plans = []
        self.jobs = []
        for sid in (SID, 'TMS_DMSO'):
            ids = []
            for number in range(1, 13):
                job_id = sid + '_' + str(number)
                ids.append(job_id)
                relative = f'goat/{sid}/{number:02d}'
                (self.root / relative).mkdir(parents=True)
                self.jobs.append({'job_id': job_id, 'species_id': sid, 'scenario_id': a.SCENARIO, 'job_directory': relative})
            self.plans.append({'species_id': sid, 'compound_id': 'TMS' if sid == 'TMS_DMSO' else 'F2_E01',
                               'stable_conformer_target': 3, 'initial_candidate_batch_size': 6,
                               'initial_goat_starts': 7, 'addition_goat_starts': 4, 'goat_start_ids': ids})
        writej(self.root / 'registry/batch_plan_v5.json', {'species': self.plans})
        writej(self.root / 'registry/jobs.json', self.jobs)
        writej(self.root / 'registry/species.json', [{'species_id': p['species_id']} for p in self.plans])

    def test_default_no_execution_releases_whole_initial_batch_idempotently(self):
        with patch.object(a, '_execute_goat', side_effect=AssertionError('No QM allowed')):
            first = a.step(self.root, SID)
            second = a.step(self.root, SID)
        self.assertEqual(first['status'], 'WAITING_WHOLE_GOAT_BATCH')
        self.assertEqual(first['released_goat_ids'], second['released_goat_ids'])
        self.assertEqual([len(v) for v in first['released_goat_ids'].values()], [7, 7])
        self.assertEqual(first['quantum_jobs_launched'], 0)
        self.assertFalse(first['scientific_convergence_approved'])
        self.assertFalse(list((self.root / 'automation').glob('*/wave_*/workspace')))

    def test_explicit_execution_requires_paths(self):
        with self.assertRaisesRegex(ValueError, 'requires ORCA'):
            a.step(self.root, SID, execute=True)

    def test_explicit_step_runs_at_most_one_job(self):
        with patch.object(a, '_execute_goat', return_value={'status': 'SYNTHETIC_TEST_ONLY'}) as run:
            result = a.step(self.root, SID, execute=True, orca='dummy', xtb='dummy', scratch_root='dummy')
        self.assertEqual(run.call_count, 1)
        self.assertEqual(result['quantum_jobs_launched'], 1)
        self.assertEqual(len(result['released_goat_ids'][SID]), 7)

    def test_whole_batch_running_prevents_early_selection(self):
        (self.root / self.jobs[0]['job_directory'] / '.goat-run.lock').write_text('{}')
        with patch('clusterflow.autoselect.screen_batch', side_effect=AssertionError('Must wait whole batch')):
            result = a.step(self.root, SID)
        self.assertEqual(result['status'], 'WAITING_WHOLE_GOAT_BATCH')
        self.assertEqual(result['job_statuses'][SID][self.jobs[0]['job_id']]['status'], 'RUNNING')

    def test_lock_precedes_apparently_complete_receipt(self):
        row = self.jobs[0]
        directory = self.root / row['job_directory']
        writej(directory / 'receipt.json', {'status': 'completed'})
        (directory / '.goat-run.lock').write_text('{}')
        with patch.object(a.pipeline, '_goat_source', side_effect=AssertionError('Must not read active outputs')):
            result = a._job_status(self.root, row, goat=True)
        self.assertFalse(result['terminal'])

    def test_failed_job_is_never_complete_or_automatically_retried(self):
        row = self.jobs[0]
        writej(self.root / row['job_directory'] / 'receipt.json', {'status': 'FAILED'})
        self.assertEqual(a._job_status(self.root, row, goat=True)['status'], 'FAILED_REVIEW')

    def test_output_without_receipt_requires_review(self):
        row = self.jobs[0]
        (self.root / row['job_directory'] / 'goat.out').write_text('SYNTHETIC INCOMPLETE OUTPUT')
        self.assertEqual(a._job_status(self.root, row, goat=True)['status'], 'INCOMPLETE_REVIEW')

    def test_changed_protocol_cannot_resume_old_state(self):
        a.step(self.root, SID)
        changed = copy.deepcopy(self.cfg)
        changed['temperature_K'] = 273.15
        writej(self.root / 'config/protocol.json', changed)
        with self.assertRaisesRegex(ValueError, 'changed'):
            a.step(self.root, SID)

    def test_exclusive_controller_lock(self):
        folder = self.root / 'automation' / SID
        folder.mkdir(parents=True)
        (folder / '.autoflow.lock').write_text('{}')
        with self.assertRaises(FileExistsError):
            a.step(self.root, SID)

    def test_stable_shortage_prevents_any_property_preparation(self):
        plans = {SID: self.plans[0]}
        stable = {'eligible_count': 2, 'all_prepared_optfreq_pass': True, 'issues': []}
        with patch.object(a.pipeline, '_validated_imports'), patch.object(a.pipeline, 'prepare_prescreen'), \
                patch.object(a.pipeline, 'manifest', return_value=[{'job_dir': 'dummy'}]), \
                patch.object(a, '_job_status', return_value={'terminal': True, 'status': 'COMPLETE'}), \
                patch.object(a.pipeline, 'advance', return_value={'status': 'advanced'}) as advance, \
                patch.object(a.pipeline, 'final_stable_candidates', return_value=stable):
            result = a._stage_step(self.root, plans, execute=False, orca=None, scratch_root=None)
        self.assertEqual(result['status'], 'DFT_STABLE_MINIMA_SHORTAGE')
        self.assertEqual([call.args[-1] for call in advance.call_args_list], ['prescreen'])

    def test_approximate_minima_never_qualify_for_automatic_nmr(self):
        stable = {'eligible_count': 5, 'all_prepared_optfreq_pass': False,
                  'issues': ['human-reviewed negative mode is not certified stable']}
        with patch.object(a.pipeline, '_validated_imports'), patch.object(a.pipeline, 'prepare_prescreen'), \
                patch.object(a.pipeline, 'manifest', return_value=[{'job_dir': 'dummy'}]), \
                patch.object(a, '_job_status', return_value={'terminal': True, 'status': 'COMPLETE'}), \
                patch.object(a.pipeline, 'advance', return_value={'status': 'advanced'}) as advance, \
                patch.object(a.pipeline, 'final_stable_candidates', return_value=stable):
            result = a._stage_step(self.root, {SID: self.plans[0]}, execute=False, orca=None, scratch_root=None)
        self.assertEqual(result['status'], 'AWAITING_SCIENTIFIC_REVIEW')
        self.assertEqual(advance.call_count, 1)

    def test_reserve_used_before_more_goat_sources(self):
        state = {'waves': [{'selections': {SID: {'available_unique_count': 30, 'energy_eligible_count': 10,
                                               'selected': [1, 2, 3, 4, 5, 6]}}}],
                 'targets': {SID: 6}, 'released': {SID: self.plans[0]['goat_start_ids'][:7]}}
        changes = a._expand_candidates(state, {SID: self.plans[0]}, [SID])
        self.assertEqual(changes[SID]['kind'], 'EXISTING_GOAT_RESERVE')
        self.assertEqual(state['targets'][SID], 10)
        self.assertEqual(len(state['released'][SID]), 7)

    def test_high_energy_reserve_does_not_fill_quota(self):
        state = {'waves': [{'selections': {SID: {'available_unique_count': 100, 'energy_eligible_count': 6,
                                               'selected': [1, 2, 3, 4, 5, 6]}}}],
                 'targets': {SID: 6}, 'released': {SID: list(self.plans[0]['goat_start_ids'])}}
        self.assertEqual(a._expand_candidates(state, {SID: self.plans[0]}, [SID]), {})

    def test_finite_additions_and_no_duplicate_release(self):
        state = {'released': {SID: self.plans[0]['goat_start_ids'][:7]}}
        plans = {SID: self.plans[0]}
        self.assertEqual(len(a._release(state, plans)), 4)
        self.assertEqual(len(a._release(state, plans)), 1)
        self.assertEqual(a._release(state, plans), [])
        self.assertEqual(len(set(state['released'][SID])), 12)

    def test_advisory_nmr_comparison_excludes_exchangeable_h(self):
        before = {'metadata': {k: 'x' for k in ('species_id', 'scenario_id', 'temperature_K', 'j_method', 'protocol_hash')},
                  'atom_rows': [{'atom_index_zero_based': 0, 'nucleus': '1H', 'shift_ppm': 1., 'exchangeable_H': False},
                                {'atom_index_zero_based': 1, 'nucleus': '13C', 'shift_ppm': 10.},
                                {'atom_index_zero_based': 2, 'nucleus': '1H', 'shift_ppm': 5., 'exchangeable_H': True}],
                  'j_pairs': []}
        after = copy.deepcopy(before)
        after['atom_rows'][0]['shift_ppm'] += .02
        after['atom_rows'][1]['shift_ppm'] += .2
        after['atom_rows'][2]['shift_ppm'] += 2.
        result = a.compare_predictions(before, after)
        self.assertTrue(result['shifts_within_reference_limits'])
        self.assertFalse(result['scientific_convergence_approved'])
        after['atom_rows'][0]['shift_ppm'] += 1.
        self.assertFalse(a.compare_predictions(before, after)['shifts_within_reference_limits'])

    def test_confined_paths_and_bound_copy(self):
        with self.assertRaises(ValueError):
            a._safe(self.root, '../escape')
        source = self.root / 'source'
        dest = self.root / 'copy'
        source.write_text('original')
        a._copy_bound(source, dest)
        source.write_text('changed')
        with self.assertRaisesRegex(ValueError, 'replace'):
            a._copy_bound(source, dest)

    def test_real_synthetic_handoff_optfreq_collapse_blocks_nmr(self):
        """Actual parser/import/DFT-input flow, but ONLY fabricated unit-test outputs."""
        from rdkit import Chem
        from clusterflow import geometry
        from clusterflow.common import sha
        import test_cluster_pipeline as fixture
        source_rows = readj(ROOT / 'registry/jobs.json')
        species_rows = readj(ROOT / 'registry/species.json')
        jobs, plans = [], []
        for sid in (SID, 'TMS_DMSO'):
            row = copy.deepcopy(next(r for r in source_rows if r['species_id'] == sid and r['start_index'] == 1))
            jobs.append(row)
            folder = self.root / row['job_directory']
            folder.mkdir(parents=True)
            original = ROOT / row['job_directory']
            for name in ('goat.inp', 'start.xyz', 'cluster.mol', 'atom_map.json', 'provenance.json'):
                shutil.copyfile(original / name, folder / name)
            texts = []
            for index in range(1, 4):
                model = Chem.MolFromMolBlock((original.parent / f'start_{index:02d}' / 'cluster.mol').read_text(encoding='utf-8'), removeHs=False)
                texts.append(geometry.xyz_text(model, f'SYNTHETIC UNIT TEST Energy = {-100 + index * .0001:.8f} Eh'))
            (folder / 'goat.finalensemble.xyz').write_text(''.join(texts), encoding='utf-8')
            (folder / 'goat.out').write_text('SYNTHETIC UNIT TEST ONLY\nProgram Version 6.1.1\nORCA TERMINATED NORMALLY\n')
            writej(folder / 'receipt.json', {'status': 'completed', 'job_id': row['job_id'],
                   'input_sha256': {n: sha(folder / n) for n in ('goat.inp', 'start.xyz', 'cluster.mol', 'atom_map.json', 'provenance.json')},
                   'output_sha256': {n: sha(folder / n) for n in ('goat.out', 'goat.finalensemble.xyz')}})
            plans.append({'species_id': sid, 'stable_conformer_target': 3, 'initial_candidate_batch_size': 6,
                          'initial_goat_starts': 1, 'addition_goat_starts': 1, 'goat_start_ids': [row['job_id']]})
        writej(self.root / 'registry/jobs.json', jobs)
        writej(self.root / 'registry/species.json', [s for s in species_rows if s['species_id'] in (SID, 'TMS_DMSO')])
        writej(self.root / 'registry/batch_plan_v5.json', {'species': plans})
        with patch('subprocess.Popen', side_effect=AssertionError('No real subprocess allowed')):
            first = a.step(self.root, SID)
            self.assertEqual(first['status'], 'WAITING_DFT_BATCH')
            self.assertEqual(len(first['jobs']), 3)
            workspace = Path(first['workspace'])
            same = a.step(self.root, SID)
            self.assertEqual(same['jobs'], first['jobs'])
            for _ in range(6):
                result = a.step(self.root, SID)
                if result['status'] != 'WAITING_DFT_BATCH':
                    break
                first_geometry = None
                for row in result['jobs']:
                    directory = workspace / row['job_dir']
                    fixture.ClusterPipelineTests.completed_opt(self, directory)
                    if row['stage'] == 'optfreq':
                        # Simulate three actual optimizations falling into one basin.
                        if first_geometry is None:
                            first_geometry = (directory / 'job.xyz').read_bytes()
                        (directory / 'job.xyz').write_bytes(first_geometry)
                        receipt = readj(directory / 'receipt.json')
                        receipt['output_sha256']['job.xyz'] = sha(directory / 'job.xyz')
                        writej(directory / 'receipt.json', receipt)
            result = a.step(self.root, SID)
        self.assertEqual(result['status'], 'AWAITING_SCIENTIFIC_REVIEW')
        self.assertIn('target unmet', result['reason'])
        self.assertEqual(result['stable_reports'][SID]['eligible_count'], 1)
        self.assertFalse(a.pipeline.manifest(workspace, 'shielding'))
        self.assertFalse(a.pipeline.manifest(workspace, 'j_PBE_pcJ3'))
        # Frozen source alteration is never silently accepted on resume.
        source_output = self.root / jobs[0]['job_directory'] / 'goat.out'
        source_output.write_text('tampered')
        with self.assertRaisesRegex(ValueError, 'binding changed'):
            a.step(self.root, SID)


if __name__ == '__main__':
    unittest.main()
