"""Independent v5 adversarial software audit; all ORCA output is SYNTHETIC.

Temporary engineering fixtures are not actual QM or experimental evidence.
No quantum executable, shell subprocess, or cluster submission may run here.
"""
from __future__ import annotations

import copy
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).parent))
from clusterflow import autoselect, geometry, pipeline, autoflow, gates, convergence
from clusterflow.common import readj, writej, sha, config, _io_path
import test_cluster_pipeline as fixture_library
from test_v4_integration import approve_negative_fixture
import test_independent_workflow as historical_fixtures

SCENARIO = 'DMSO_02_W01_trace_water_sensitivity'
SPECIES = fixture_library.SPECIES


class IndependentNativeGoatFormat(unittest.TestCase):
    """Format examples follow the public ORCA 6.1 tutorial, not real results."""

    @staticmethod
    def frames(comments):
        return geometry.parse_xyz_frames(''.join(f'1\n{c}\nH 0 0 0\n' for c in comments))

    def test_documented_zero_index_multiline_native_energy_table(self):
        out = '''SYNTHETIC FORMAT FIXTURE ONLY
        # Final ensemble info #
        Conformer     Energy     Degen.   % total   % cumul.
                      (kcal/mol)
        -------------------------------------------------
              0       0.000        1        80.0       80.0
              1       0.628        1        20.0      100.0
        Lowest energy conformer : -10.000000 Eh
        ORCA TERMINATED NORMALLY
        '''
        values = autoselect.parse_frame_energies(self.frames(['-10.0', '-9.999']), out)
        self.assertEqual(values, [-10.0, -9.999])

    def test_native_table_cannot_relabel_kj_as_kcal(self):
        out = '''# Final ensemble info #
        Conformer Energy Degen. % total % cumul.
                  (kJ/mol)
        0 0.000 1 80.0 80.0
        1 0.628 1 20.0 100.0
        Lowest energy conformer : -10.000000 Eh
        '''
        with self.assertRaises(ValueError):
            autoselect.parse_frame_energies(self.frames(['-10.0', '-9.999']), out)

    def test_native_table_must_corroborate_every_frame(self):
        out = '''# Final ensemble info #
        Conformer Energy Degen. % total % cumul.
                  (kcal/mol)
        0 0.000 1 80.0 80.0
        Lowest energy conformer : -10.000000 Eh
        '''
        with self.assertRaises(ValueError):
            autoselect.parse_frame_energies(self.frames(['-10.0', '-9.999']), out)

    def test_explicit_units_sortable_but_nonfinite_or_unitless_unbound_rejected(self):
        self.assertEqual(autoselect.parse_frame_energies(self.frames(['-10.1 Eh', '-10.2 hartree']), ''),
                         [-10.1, -10.2])
        for comment in ['NaN Eh', 'inf Eh', '-10.0', '0 kcal/mol', '-10 Eh and -20 Eh']:
            with self.subTest(comment=comment), self.assertRaises(ValueError):
                autoselect.parse_frame_energies(self.frames([comment]), '')


class IndependentAutomaticSelection(unittest.TestCase):
    def setUp(self):
        self.fixture = fixture_library.ClusterPipelineTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        automation_path = ROOT / 'config/automation.json'
        settings = readj(automation_path) if automation_path.is_file() else {
            'automatic_handoff': {'enabled': True},
            'selection': {'rmsd_A': .1, 'representative_energy_window_kcal_mol': 6.0}}
        writej(self.root / 'config/automation.json', settings)
        self.rows = [self.fixture.row]

    def real_initial_model(self, start):
        path = ROOT / '01_goat_explicit' / SPECIES / SCENARIO / f'start_{start:02d}' / 'cluster.mol'
        return Chem.MolFromMolBlock(path.read_text(encoding='utf-8'), removeHs=False)

    def add_source(self, start):
        row = copy.deepcopy(next(r for r in readj(ROOT / 'registry/jobs.json')
                                  if r['species_id'] == SPECIES and r['start_index'] == start))
        source, target = ROOT / row['job_directory'], self.root / row['job_directory']
        target.mkdir(parents=True)
        for name in ('goat.inp', 'start.xyz', 'cluster.mol', 'atom_map.json', 'provenance.json'):
            shutil.copyfile(source / name, target / name)
        self.rows.append(row)
        writej(self.root / 'registry/jobs.json', self.rows)
        return row

    def completed(self, row, molecules, energies):
        directory = self.root / row['job_directory']
        (directory / 'goat.finalensemble.xyz').write_text(''.join(
            geometry.xyz_text(m, f'{e:.12f} Eh SYNTHETIC UNIT TEST ONLY')
            for m, e in zip(molecules, energies)), encoding='utf-8')
        (directory / 'goat.out').write_text(
            'SYNTHETIC UNIT TEST ONLY\nProgram Version 6.1.0\nORCA TERMINATED NORMALLY\n', encoding='utf-8')
        writej(directory / 'receipt.json', {
            'job_id': row['job_id'], 'status': 'COMPLETED_NORMAL_AWAITING_SCIENTIFIC_REVIEW',
            'input_sha256': {name: sha(directory / name) for name in
                            ('goat.inp', 'start.xyz', 'atom_map.json', 'cluster.mol', 'provenance.json')},
            'output_sha256': {name: sha(directory / name) for name in ('goat.out', 'goat.finalensemble.xyz')}})

    def test_source_bindings_include_duplicate_only_completed_start(self):
        second = self.add_source(2)
        model = self.real_initial_model(1)
        self.completed(self.rows[0], [model], [-100.0])
        self.completed(second, [model], [-100.0])
        report = autoselect.screen_batch(self.root, [r['job_id'] for r in self.rows], 1)
        self.assertEqual(report['available_unique_count'], 1)
        self.assertEqual(len(report['source_bindings']), 2)
        self.assertEqual(report['rediscovery'][0]['job_count'], 2)
        self.assertFalse(report['human_review_completed'])
        self.assertFalse(report['sampling_convergence_established'])
        self.assertFalse(report['selected'][0]['DFT_minimum_verified'])
        imported = autoselect.import_selection(self.root, report)
        self.assertEqual(imported['quantum_jobs_launched'], 0)
        pipeline.prepare_prescreen(self.root, SPECIES, SCENARIO)
        # Even the unselected duplicate-only source remains immutable evidence.
        out = self.root / second['job_directory'] / 'goat.out'
        out.write_text(out.read_text(encoding='utf-8') + 'tampered', encoding='utf-8')
        with self.assertRaises(ValueError):
            pipeline._validated_imports(self.root, SPECIES, SCENARIO)

    def test_completed_receipt_with_live_runner_lock_is_not_a_finished_batch(self):
        self.completed(self.rows[0], [self.real_initial_model(1)], [-100.0])
        lock = self.root / self.rows[0]['job_directory'] / '.goat-run.lock'
        lock.write_text('{"status":"SYNTHETIC_RUNNING"}', encoding='utf-8')
        with self.assertRaises(ValueError):
            autoselect.screen_batch(self.root, [self.rows[0]['job_id']], 1)
        self.assertFalse((self.root / 'imports').exists())

    def test_seed_only_never_becomes_completed_candidate(self):
        with self.assertRaises(ValueError):
            autoselect.screen_batch(self.root, [self.rows[0]['job_id']], 1)
        self.assertFalse((self.root / 'imports').exists())

    def test_resealed_selection_cannot_change_energy_or_claim_human_approval(self):
        self.completed(self.rows[0], [self.real_initial_model(1)], [-100.0])
        original = autoselect.screen_batch(self.root, [self.rows[0]['job_id']], 1)
        for field, value in [('human_review_completed', True), ('sampling_convergence_established', True)]:
            changed = copy.deepcopy(original)
            changed[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                autoselect.import_selection(self.root, changed)
        changed = copy.deepcopy(original)
        changed['selected'][0]['energy_Eh'] -= 1
        with self.assertRaises(ValueError):
            autoselect.import_selection(self.root, changed)
        self.assertFalse((self.root / 'imports').exists())

    def test_high_energy_independent_structure_is_not_forced_to_fill_quota(self):
        models = [self.real_initial_model(i) for i in (1, 2, 3)]
        self.completed(self.rows[0], models, [-100, -99.999, -99.8])
        report = autoselect.screen_batch(self.root, [self.rows[0]['job_id']], 3)
        self.assertEqual(report['available_unique_count'], 3)
        self.assertEqual(report['energy_eligible_count'], 2)
        self.assertFalse(report['target_reached'])
        self.assertEqual(len(report['selected']), 2)
        self.assertEqual(len(report['reserve']), 1)
        self.assertGreater(report['reserve'][0]['relative_energy_kcal_mol'], 100)
        with self.assertRaises(ValueError):
            autoselect.import_selection(self.root, report)


class IndependentFinalStableMinimumGate(IndependentAutomaticSelection):
    """Real pipeline checks on fabricated temporary outputs, not scientific proof."""

    def prepared_final_stage(self, number=3):
        models = [self.real_initial_model(i) for i in range(1, number + 1)]
        self.completed(self.rows[0], models, [-100 + .001 * i for i in range(number)])
        screened = autoselect.screen_batch(self.root, [self.rows[0]['job_id']], number)
        self.assertEqual(len(screened['selected']), number)
        autoselect.import_selection(self.root, screened)
        batch = pipeline.prepare_prescreen(self.root, SPECIES, SCENARIO)
        for path in batch['jobs']:
            self.fixture.completed_opt(self.root / path)
        advanced = pipeline.advance(self.root, SPECIES, SCENARIO, 'prescreen')
        self.assertEqual(advanced['status'], 'advanced')
        self.assertEqual(len(advanced['prepared_jobs']), number)
        directories = [self.root / path for path in advanced['prepared_jobs']]
        for directory in directories:
            self.fixture.completed_opt(directory)
        return directories

    def edit_completed_output(self, directory, *, model=None, gibbs=None, negative=False):
        self.fixture.completed_opt(directory, negative=negative)
        if model is not None:
            (directory / 'job.xyz').write_text(geometry.xyz_text(model, 'SYNTHETIC FINAL DFT OUTPUT'), encoding='utf-8')
        if gibbs is not None:
            path = directory / 'job.out'
            path.write_text(path.read_text(encoding='utf-8').replace('-123.456 Eh', f'{gibbs} Eh'), encoding='utf-8')
        receipt = readj(directory / 'receipt.json')
        receipt['output_sha256'] = {name: sha(directory / name) for name in ('job.out', 'job.xyz')}
        writej(directory / 'receipt.json', receipt)

    def assert_no_nmr_jobs(self):
        for stage in ('shielding', 'j_PBE_pcJ3', 'j_PBE0_pcJ2'):
            self.assertEqual(pipeline.manifest(self.root, stage, SPECIES, SCENARIO), [])

    def test_goat_count_is_not_final_dft_minimum_count(self):
        directories = self.prepared_final_stage()
        collapsed = Chem.MolFromMolBlock((directories[0] / 'input_structure.mol').read_text(encoding='utf-8'), removeHs=False)
        self.edit_completed_output(directories[1], model=collapsed)
        report = pipeline.final_stable_candidates(self.root, SPECIES, SCENARIO)
        self.assertTrue(report['all_prepared_optfreq_pass'])
        self.assertEqual(report['eligible_count'], 2)
        self.assertEqual(len(report['duplicates']), 1)
        self.assertFalse(report['sampling_convergence_established'])
        stage = autoflow._stage_step(self.root, {SPECIES: {'stable_conformer_target': 3}},
                                    execute=False, orca=None, scratch_root=None)
        self.assertEqual(stage['status'], 'DFT_STABLE_MINIMA_SHORTAGE')
        self.assert_no_nmr_jobs()

    def test_reviewed_imaginary_mode_approximation_never_counts_as_stable(self):
        directories = self.prepared_final_stage()
        self.edit_completed_output(directories[2], negative=True)
        approve_negative_fixture(directories[2], config(self.root))
        self.assertTrue(gates.inspect_job(directories[2], config(self.root))['accepted'])
        report = pipeline.final_stable_candidates(self.root, SPECIES, SCENARIO)
        self.assertEqual(report['eligible_count'], 2)
        self.assertFalse(report['all_prepared_optfreq_pass'])
        self.assertTrue(any('not_verified_stable' in error for item in report['issues'] for error in item['errors']))
        self.assert_no_nmr_jobs()

    def test_nondefault_final_energy_window_controls_eligible_not_report_only(self):
        directories = self.prepared_final_stage()
        self.edit_completed_output(directories[2], gibbs=-123.450)
        baseline = pipeline.final_stable_candidates(self.root, SPECIES, SCENARIO)
        self.assertEqual(baseline['eligible_count'], 3)
        policy = readj(self.root / 'config/automation.json')
        policy.setdefault('stable_minima', {})['energy_window_kcal_mol'] = 1.0
        writej(self.root / 'config/automation.json', policy)
        result = pipeline.final_stable_candidates(self.root, SPECIES, SCENARIO)
        self.assertEqual(result['maximum_relative_G_kcal_mol'], 1.0)
        self.assertEqual(result['eligible_count'], 2)
        self.assertEqual(result['unique_stable_count'], 3)
        self.assert_no_nmr_jobs()

    def test_final_source_edit_invalidates_stability_count(self):
        directories = self.prepared_final_stage()
        self.assertEqual(pipeline.final_stable_candidates(self.root, SPECIES, SCENARIO)['eligible_count'], 3)
        path = directories[0] / 'job.out'
        path.write_text(path.read_text(encoding='utf-8') + 'tampered after receipt', encoding='utf-8')
        result = pipeline.final_stable_candidates(self.root, SPECIES, SCENARIO)
        self.assertFalse(result['all_prepared_optfreq_pass'])
        self.assertEqual(result['eligible_count'], 2)
        self.assert_no_nmr_jobs()

    def test_direct_advance_cannot_bypass_stable_minimum_target(self):
        directories = self.prepared_final_stage()
        model = Chem.MolFromMolBlock((directories[0] / 'input_structure.mol').read_text(encoding='utf-8'), removeHs=False)
        self.edit_completed_output(directories[1], model=model)
        writej(self.root / 'registry/batch_plan_v5.json', {'species': [{
            'species_id': SPECIES, 'scenario_id': SCENARIO, 'stable_conformer_target': 3}]})
        result = pipeline.advance(self.root, SPECIES, SCENARIO, 'optfreq')
        self.assertEqual(result['status'], 'BLOCKED_FINAL_STABLE_TARGET')
        self.assertEqual(result['eligible_stable_count'], 2)
        self.assertEqual(result['prepared_jobs'], [])
        self.assert_no_nmr_jobs()

    def test_sufficient_good_minima_do_not_silently_admit_reviewed_negative_mode(self):
        directories = self.prepared_final_stage(4)
        self.edit_completed_output(directories[3], negative=True)
        approve_negative_fixture(directories[3], config(self.root))
        self.assertEqual(pipeline.final_stable_candidates(self.root, SPECIES, SCENARIO)['eligible_count'], 3)
        result = autoflow._stage_step(self.root, {SPECIES: {'stable_conformer_target': 3}},
                                     execute=False, orca=None, scratch_root=None)
        self.assertEqual(result['status'], 'AWAITING_SCIENTIFIC_REVIEW')
        self.assert_no_nmr_jobs()

    def test_same_geometry_different_energy_cannot_inflate_stable_count(self):
        directories = self.prepared_final_stage()
        model = Chem.MolFromMolBlock((directories[0] / 'input_structure.mol').read_text(encoding='utf-8'), removeHs=False)
        self.edit_completed_output(directories[1], model=model, gibbs=-123.455)
        report = pipeline.final_stable_candidates(self.root, SPECIES, SCENARIO)
        self.assertLessEqual(report['eligible_count'], 2)
        self.assert_no_nmr_jobs()


class IndependentV5PlanContract(unittest.TestCase):
    def test_targets_are_final_minima_and_reference_is_separate(self):
        plan = readj(ROOT / 'registry/batch_plan_v5.json')
        rows = plan['species']
        self.assertEqual(len(plan['products']), 64)
        self.assertEqual(len(rows), 67)
        self.assertEqual(sum(r['stable_conformer_target'] for r in rows if r['compound_id'] != 'TMS'), 526)
        self.assertEqual(next(r for r in rows if r['compound_id'] == 'TMS')['stable_conformer_target'], 3)
        for row in rows:
            self.assertIn(row['stable_conformer_target'], (3, 5, 10, 15))
            self.assertGreaterEqual(row['initial_candidate_batch_size'], row['stable_conformer_target'])
            self.assertEqual((row['explicit_DMSO_d6'], row['explicit_H2O'], row['temperature_K']), (2, 1, 298.15))
            self.assertFalse(row['sampling_converged'])
        self.assertFalse(config(ROOT)['execution']['quantum_jobs_started_during_package_build'])

    def test_shift_change_advisory_never_approves_convergence(self):
        before = {'metadata': {'species_id': 'S', 'scenario_id': SCENARIO, 'temperature_K': 298.15,
                              'j_method': 'PBE_pcJ3', 'protocol_hash': 'synthetic'},
                  'atom_rows': [{'atom_index_zero_based': 0, 'nucleus': '1H', 'shift_ppm': 1.0,
                                 'exchangeable_H': False},
                                {'atom_index_zero_based': 1, 'nucleus': '13C', 'shift_ppm': 10.0,
                                 'exchangeable_H': False}], 'j_pairs': []}
        after = copy.deepcopy(before)
        after['atom_rows'][0]['shift_ppm'] += .4
        after['atom_rows'][1]['shift_ppm'] += 2
        result = autoflow.compare_predictions(before, after)
        self.assertEqual(result['status'], 'ADVISORY_ONLY')
        self.assertFalse(result['shifts_within_reference_limits'])
        self.assertFalse(result['scientific_convergence_approved'])
        small = autoflow.compare_predictions(before, copy.deepcopy(before))
        self.assertTrue(small['shifts_within_reference_limits'])
        self.assertFalse(small['scientific_convergence_approved'])

    def test_legacy_gate_and_v5_advisory_policy_are_explicitly_distinct(self):
        # Numeric-policy unit test ONLY. Actual parser/source checks are exercised
        # without mocks in IndependentTwoWaveEndToEnd below.
        with tempfile.TemporaryDirectory(prefix='v5advice_') as directory:
            root = Path(directory)
            previous = historical_fixtures.IndependentBenchmarkAndConvergence()
            gate = previous.convergence_fixture(root)
            result_path = root / gate['waves'][-1]['result']['path']
            result = readj(result_path)
            for role, key in [('target', 'source'), ('reference', 'reference_source')]:
                path = root / f'{role}_final.json'
                writej(path, result['metadata']['ensemble_snapshot'][role])
                result['metadata'][key] = {'ensemble_path': path.name, 'ensemble_sha256': sha(path)}
            writej(result_path, result)
            gate['waves'][-1]['result']['sha256'] = sha(result_path)
            writej(root / 'gate.json', gate)
            settings = config(ROOT)
            atoms = {0: {'element': 'H', 'exchangeable_H': False}, 1: {'element': 'C'}}
            for advisory in (True, False):
                settings['conformer_policy']['convergence_thresholds_are_advisory'] = advisory
                writej(root / 'config/protocol.json', settings)
                with patch.object(convergence, '_completed_sampling', return_value={'TARGET': 1, 'TMS_DMSO': 1}), \
                     patch.object(convergence, '_recompute_snapshot', side_effect=[({0: 1., 1: 10.}, atoms),
                                                                                ({0: 1.2, 1: 11.}, atoms)]):
                    if advisory:
                        report = convergence.validate_gate(root, 'gate.json', 'TARGET', SCENARIO, 'PBE0_pcJ2')
                        self.assertTrue(report['numerical_references_are_advisory'])
                        self.assertEqual(set(report['advisory_exceedances']), {'nonexchangeable_1H_ppm', '13C_ppm'})
                    else:
                        with self.assertRaisesRegex(ValueError, 'stability gate failed'):
                            convergence.validate_gate(root, 'gate.json', 'TARGET', SCENARIO, 'PBE0_pcJ2')

    def test_held_release_requires_current_identity_and_prefix_pool(self):
        from scripts.goat_job import controller_released
        with tempfile.TemporaryDirectory(prefix='v5held_') as directory:
            root = Path(directory)
            names = ('registry/batch_plan_v5.json', 'registry/jobs.json', 'registry/species.json',
                     'config/protocol.json', 'config/automation.json')
            writej(root / names[0], {'species': [{'species_id': SPECIES, 'goat_start_ids': ['job1', 'job2']}]})
            writej(root / names[1], [])
            writej(root / names[2], [])
            writej(root / names[3], config(ROOT))
            writej(root / names[4], {'automatic_handoff': {'enabled': True}})
            state = {'species_id': SPECIES, 'identity_sha256': {n: sha(root / n) for n in names},
                     'released': {SPECIES: ['job1']}}
            path = root / 'automation' / SPECIES / 'state.json'
            writej(path, state)
            self.assertTrue(controller_released(root, 'job1', SPECIES))
            self.assertFalse(controller_released(root, 'job2', SPECIES))
            state['released'][SPECIES] = ['job2']
            writej(path, state)
            self.assertFalse(controller_released(root, 'job2', SPECIES))
            state['released'][SPECIES] = ['job1', 'job2']
            state['identity_sha256']['config/protocol.json'] = '0' * 64
            writej(path, state)
            self.assertFalse(controller_released(root, 'job2', SPECIES))


class IndependentTwoWaveEndToEnd(unittest.TestCase):
    """Unmocked scientific/parser gates, SYNTHETIC output files, no QM execution."""
    short_deployment = True

    def setUp(self):
        self.base = IndependentAutomaticSelection()
        # Test in a package-like short deployment root. Long Windows TEMP roots
        # have a separate IO-path regression; never assert universal deep-path
        # support for third-party plotting/native chemistry executables.
        factory = tempfile.TemporaryDirectory
        def short_temp(*args, **kwargs):
            kwargs['dir'] = str(ROOT.parents[1])
            kwargs['prefix'] = 'v5qa_'
            return factory(*args, **kwargs)
        if self.short_deployment:
            with patch.object(fixture_library.tempfile, 'TemporaryDirectory', side_effect=short_temp):
                self.base.setUp()
        else:
            self.base.setUp()
        self.addCleanup(self.base.doCleanups)
        self.root = self.base.root
        source_species = readj(ROOT / 'registry/species.json')
        source_jobs = readj(ROOT / 'registry/jobs.json')
        self.ids = [SPECIES, 'TMS_DMSO']
        writej(self.root / 'registry/species.json', [r for r in source_species if r['species_id'] in self.ids])
        rows = []
        for sid in self.ids:
            for start in (1, 2):
                row = copy.deepcopy(next(r for r in source_jobs if r['species_id'] == sid and r['start_index'] == start))
                target = self.root / row['job_directory']
                if not target.exists():
                    target.mkdir(parents=True)
                    for name in ('goat.inp', 'start.xyz', 'cluster.mol', 'atom_map.json', 'provenance.json'):
                        shutil.copyfile(ROOT / row['job_directory'] / name, target / name)
                rows.append(row)
        self.rows = rows
        writej(self.root / 'registry/jobs.json', rows)
        plan = readj(ROOT / 'registry/batch_plan_v5.json')
        plan['species'] = [r for r in plan['species'] if r['species_id'] in self.ids]
        for row in plan['species']:
            row['goat_start_ids'] = [r['job_id'] for r in rows if r['species_id'] == row['species_id']]
            row['initial_goat_starts'] = 1
            row['addition_goat_starts'] = 1
            row['stable_conformer_target'] = 3
            row['initial_candidate_batch_size'] = 3
        plan['engineering_test_notice'] = 'SYNTHETIC TEMPORARY CONTROLLER FIXTURE; NOT PRODUCTION PLAN'
        writej(self.root / 'registry/batch_plan_v5.json', plan)

    def completed_source(self, sid, start):
        row = next(r for r in self.rows if r['species_id'] == sid and r['start_index'] == start)
        models = []
        for number in range(1 + (start - 1) * 3, 4 + (start - 1) * 3):
            path = ROOT / '01_goat_explicit' / sid / SCENARIO / f'start_{number:02d}' / 'cluster.mol'
            models.append(Chem.MolFromMolBlock(path.read_text(encoding='utf-8'), removeHs=False))
        self.base.completed(row, models, [-100 + .001 * (i + (start - 1) * 3) for i in range(3)])

    def completed_property(self, workspace, row):
        directory = workspace / row['job_dir']
        meta = readj(directory / 'job.json')
        amap = readj(directory / 'atom_map.json')
        output = ['SYNTHETIC END-TO-END UNIT TEST ONLY', 'Program Version 6.1.0']
        outputs = ['job.out']
        if meta['stage'] == 'shielding':
            output.append('CHEMICAL SHIELDING SUMMARY (ppm)')
            for idx in meta['expected_nuclei_zero_based']:
                symbol = amap[idx]['element']
                sigma = {'H': 31.0, 'C': 190.0}[symbol]
                if meta['species_id'] != 'TMS_DMSO':
                    sigma -= (1 + .01 * idx) if symbol == 'H' else (20 + .1 * idx)
                output.append(f'{idx} {symbol} {sigma:.10f} 0.0')
            output.append('END OF SYNTHETIC SHIELDING SECTION')
        else:
            pairs = meta['expected_j_pairs_zero_based']
            prop = ['$SCF_Spin_Spin_Coupling', f'&pairsInfo [ &Dim ({len(pairs)}, 4) ]']
            prop += [f'{k} {i} 1 {j} 1' for k, (i, j) in enumerate(pairs)]
            prop += [f'&pairsTotalSSCIso [ &Dim ({len(pairs)}, 1) ]']
            prop += [f'{k} 7.0' for k in range(len(pairs))]
            prop += ['$End']
            (directory / 'job.property.txt').write_text('\n'.join(prop) + '\n', encoding='utf-8')
            outputs.append('job.property.txt')
        output += ['ORCA TERMINATED NORMALLY']
        (directory / 'job.out').write_text('\n'.join(output) + '\n', encoding='utf-8')
        writej(directory / 'receipt.json', {
            'status': 'completed',
            'input_sha256': {n: sha(directory / n) for n in
                            ('job.inp', 'job.json', 'input.xyz', 'input_structure.mol', 'atom_map.json')},
            'output_sha256': {n: sha(directory / n) for n in outputs}})
        checked = gates.inspect_job(directory, config(workspace))
        self.assertTrue(checked['accepted'], checked['errors'])

    def drive_one_wave(self):
        for _ in range(15):
            result = autoflow.step(self.root, SPECIES, execute=False)
            self.assertEqual(result['quantum_jobs_launched'], 0)
            if result['status'] in ('WAITING_DFT_BATCH', 'WAITING_PROPERTY_BATCH'):
                workspace = Path(result['workspace'])
                for row in result['jobs']:
                    directory = workspace / row['job_dir']
                    if (directory / 'receipt.json').is_file():
                        continue
                    if row['stage'] in ('prescreen', 'optfreq'):
                        self.base.fixture.completed_opt(directory)
                    else:
                        self.completed_property(workspace, row)
                continue
            return result
        self.fail('Synthetic controller did not finish its bounded stage sequence')

    def test_two_waves_matching_tms_automatic_handoff_no_false_validation(self):
        self.completed_source(SPECIES, 1)
        waiting = autoflow.step(self.root, SPECIES, execute=False)
        self.assertEqual(waiting['status'], 'WAITING_WHOLE_GOAT_BATCH')
        self.assertEqual(list((self.root / 'automation').glob('*/wave_*')), [])
        self.completed_source('TMS_DMSO', 1)
        first = self.drive_one_wave()
        self.assertEqual(first['status'], 'FIRST_OR_EXPANDED_NMR_WAVE_COMPLETE_NEXT_BATCH_RELEASED')
        state_path = self.root / 'automation' / SPECIES / 'state.json'
        state = readj(state_path)
        self.assertEqual(len(state['waves']), 1)
        self.assertFalse(state['scientific_convergence_approved'])
        first_prediction = self.root / state['waves'][0]['workspace'] / state['waves'][0]['prediction']['path']
        frozen_digest = sha(first_prediction)
        self.completed_source(SPECIES, 2)
        self.completed_source('TMS_DMSO', 2)
        second = self.drive_one_wave()
        self.assertEqual(second['status'], 'AWAITING_SCIENTIFIC_REVIEW')
        state = readj(state_path)
        self.assertEqual(sum('prediction' in w for w in state['waves']), 2)
        self.assertFalse(state['scientific_convergence_approved'])
        self.assertEqual(sha(first_prediction), frozen_digest)
        self.assertTrue(state['advice']['shifts_within_reference_limits'])
        self.assertFalse(state['advice']['scientific_convergence_approved'])
        for wave in state['waves']:
            prediction = readj(self.root / wave['workspace'] / wave['prediction']['path'])
            reference = prediction['metadata']['ensemble_snapshot']['reference']
            self.assertEqual(reference['species_id'], 'TMS_DMSO')
            self.assertEqual(reference['scenario_id'], SCENARIO)
            self.assertEqual(reference['temperature_K'], 298.15)
            self.assertFalse(reference['sampling_convergence_established'])
        # Corruption of an old wave remains detectable after later waves exist.
        _io_path(first_prediction).write_text(_io_path(first_prediction).read_text(encoding='utf-8') + ' ', encoding='utf-8')
        with self.assertRaises(ValueError):
            autoflow.step(self.root, SPECIES, execute=False)


if __name__ == '__main__':
    unittest.main()
