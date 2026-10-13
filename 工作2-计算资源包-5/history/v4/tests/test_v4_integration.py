"""Independent v4 integration audit. All output fixtures are SYNTHETIC.

Never launches quantum software, changes production inputs, or certifies
experimental structures, convergence, or NMR agreement. Temporary test files
are deleted by unittest cleanup.
"""
from __future__ import annotations

import contextlib
import copy
import importlib.util
import io
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).parent))
from clusterflow import pipeline, gates, issue_review, policy, collect
from clusterflow.common import config, readj, writej, sha, hash_protocol
from clusterflow.spectra import make_spectrum
import test_cluster_pipeline as pf
import test_multiplet as mf

SCENARIO = 'DMSO_02_W01_trace_water_sensitivity'
DRUGS = {'F4_D01', 'F4_D02', 'F4_D03', 'F4_D04', 'F4_D05', 'F4_D07'}


def load_cli():
    spec = importlib.util.spec_from_file_location('v4_independent_cli', ROOT / 'workflow.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def approve_negative_fixture(directory, settings):
    evidence = directory / 'INDEPENDENT_SYNTHETIC_MODE_REVIEW.txt'
    evidence.write_text('SYNTHETIC UNIT TEST ONLY. NOT HUMAN SCIENTIFIC APPROVAL.', encoding='utf-8')
    value = issue_review.draft_scientific_review(directory, settings)
    value.update(status='ACCEPT_AS_REVIEWED_APPROXIMATION',
                 reviewer='Independent synthetic fixture reviewer',
                 reviewed_utc='2026-10-14T00:00:00Z',
                 rationale='Synthetic engineering fixture only.',
                 limitations='No physical minimum or real review demonstrated.',
                 checks={k: True for k in issue_review.REVIEW_CHECKS},
                 mode_review_evidence=[{'path': evidence.name, 'sha256': sha(evidence),
                                       'description': 'Synthetic fixture evidence, not experimental.'}])
    path = directory / 'scientific_review.json'
    writej(path, value)
    return path


class V4IndependentConfiguration(unittest.TestCase):
    def setUp(self):
        self.settings = config(ROOT)
        self.species = readj(ROOT / 'registry/species.json')

    def test_package_name_and_no_revdsd_active_stage(self):
        self.assertEqual(self.settings['package_name'], '工作2-计算资源包-（4）')
        self.assertFalse(any('revdsd' in key.lower() for key in self.settings))
        self.assertNotIn('revDSD', str(self.settings.get('j_methods')))

    def test_all_six_drugs_use_only_j3_and_reject_j2_request(self):
        rows = [s for s in self.species if s['compound_id'] in DRUGS]
        self.assertEqual({s['compound_id'] for s in rows}, DRUGS)
        self.assertEqual(len(rows), 6)
        for row in rows:
            self.assertEqual(policy.j_methods_for(self.settings, row['compound_id']), ['PBE_pcJ3'])
            self.assertEqual(policy.resolve_j_method(ROOT, row['species_id']), 'PBE_pcJ3')
            with self.assertRaisesRegex(ValueError, 'routing'):
                policy.resolve_j_method(ROOT, row['species_id'], 'PBE0_pcJ2')
            with self.assertRaisesRegex(ValueError, 'routing'):
                collect.collect_species(ROOT, row['species_id'], SCENARIO,
                                        'PBE0_pcJ2', exploratory=True, write=False)

    def test_non_drug_still_has_two_benchmark_candidates(self):
        for row in self.species:
            if row['compound_id'] not in DRUGS | {'TMS'}:
                self.assertEqual(set(policy.j_methods_for(self.settings, row['compound_id'])),
                                 {'PBE_pcJ3', 'PBE0_pcJ2'})

    def test_model_independent_no_continuum_prerequisite(self):
        self.assertEqual(self.settings['scenario_prerequisites'], {})
        self.assertEqual(self.settings['example_scenario'], SCENARIO)
        self.assertEqual(self.settings['scientific_issue_policy'], 'user_review')
        self.assertFalse(self.settings['execution']['automatic_cluster_submission'])
        self.assertFalse(self.settings['execution']['quantum_jobs_started_during_package_build'])

    def test_linewidths_are_h1_c2_and_adjustable_without_changing_sticks(self):
        for nucleus, symbol, key, expected in [('1H', 'H', 'proton_linewidth_Hz', 1),
                                               ('13C', 'C', 'carbon_linewidth_Hz', 2)]:
            atoms = [{'atom_index_zero_based': 0, 'element': symbol, 'shift_ppm': 3.0}]
            before = make_spectrum(atoms, [], nucleus, self.settings)
            settings = copy.deepcopy(self.settings)
            settings['spectrum'][key] = 3.25
            after = make_spectrum(atoms, [], nucleus, settings)
            self.assertEqual(before['linewidth_FWHM_Hz'], expected)
            self.assertEqual(after['linewidth_FWHM_Hz'], 3.25)
            self.assertEqual(before['sticks'], after['sticks'])

    def test_current_cli_status_never_says_calculations_finished(self):
        cli = load_cli()
        status = cli.summary()
        self.assertEqual(status['products'], 64)
        self.assertEqual(status['downstream_jobs_prepared'], 0)
        self.assertIn('NOT_RUN', status['status'])
        self.assertEqual(cli.DEFAULT_SCENARIO, SCENARIO)

    def test_display_overrides_preserve_protocol_and_reject_invalid_widths(self):
        cli = load_cli()
        original = copy.deepcopy(self.settings)
        previous_hash = hash_protocol(self.settings)
        display = cli.display_spectrum_settings(self.settings, 1.6, 2.7)
        self.assertEqual(display['proton_linewidth_Hz'], 1.6)
        self.assertEqual(display['carbon_linewidth_Hz'], 2.7)
        self.assertEqual(self.settings, original)
        self.assertEqual(hash_protocol(self.settings), previous_hash)
        for value in (0, -1, float('nan'), float('inf')):
            with self.assertRaisesRegex(ValueError, 'finite and positive'):
                cli.display_spectrum_settings(self.settings, value, None)

    def test_collect_cli_display_override_never_changes_protocol_or_executes_qm(self):
        cli = load_cli()
        protocol_sha = sha(ROOT / 'config/protocol.json')
        rows = [{'atom_index_zero_based': 0, 'element': 'H', 'shift_ppm': 2},
                {'atom_index_zero_based': 1, 'element': 'C', 'shift_ppm': 30}]
        synthetic = {'atom_rows': rows, 'j_pairs': [],
                     'result_path': 'SYNTHETIC_TEST_NO_FILE_prediction.json'}
        captured = []
        widths = []
        def render(*args):
            result = make_spectrum(*args)
            widths.append(result['linewidth_FWHM_Hz'])
            return result
        with patch.object(sys, 'argv', ['workflow.py', 'collect', '--species', pf.SPECIES,
                            '--exploratory', '--linewidth-h', '1.6', '--linewidth-c', '2.7']), \
                patch('clusterflow.collect.collect_species', return_value=synthetic), \
                patch('clusterflow.spectra.make_spectrum', side_effect=render), \
                patch('clusterflow.spectra.export_spectrum', return_value={'SYNTHETIC_NO_FILE': True}), \
                patch.object(cli, 'sha', return_value='SYNTHETIC_HASH_NOT_A_REAL_PREDICTION'), \
                patch.object(cli, 'writej', side_effect=lambda path, data: captured.append((path, data))), \
                patch('subprocess.Popen', side_effect=AssertionError('Quantum execution forbidden')), \
                contextlib.redirect_stdout(io.StringIO()):
            cli.main()
        self.assertEqual(widths, [1.6, 2.7])
        self.assertEqual(sha(ROOT / 'config/protocol.json'), protocol_sha)
        provenance = next(data for path, data in captured if path.name == 'display_settings.json')
        self.assertTrue(provenance['display_only_no_QM_rerun'])
        self.assertEqual(provenance['prediction_path'], synthetic['result_path'])

    def test_plan_is_1020_starts_not_1020_proven_minima(self):
        plan = readj(ROOT / 'registry/sampling_plan_v4.json')
        self.assertEqual(plan['summary']['total_starts'], 1020)
        self.assertEqual(plan['summary']['product_starts'], 1012)
        self.assertEqual(plan['summary']['TMS_starts'], 8)
        self.assertEqual(len(plan['products']), 64)
        for row in plan['species']:
            self.assertTrue(row['n_starts_is_not_n_final_conformers'])
            self.assertTrue(row['advisory_no_quota'])
            self.assertFalse(row['sampling_converged'])
            self.assertEqual((row['explicit_DMSO_d6'], row['explicit_H2O']), (2, 1))
            self.assertEqual(row['outer_sampling_model'], 'ALPB(DMSO)')
            self.assertEqual(row['outer_DFT_model'], 'CPCM(DMSO)')

    def test_main_registry_release_matches_plan_after_generation(self):
        plan = readj(ROOT / 'registry/sampling_plan_v4.json')
        jobs = readj(ROOT / 'registry/jobs.json')
        released = [j for j in jobs if j['initial_release']]
        self.assertEqual(len(released), 1020)
        self.assertEqual({j['scenario_id'] for j in released}, {SCENARIO})
        self.assertEqual(len({j['job_id'] for j in jobs}), len(jobs))
        for row in plan['species']:
            actual = [j for j in released if j['species_id'] == row['species_id']]
            self.assertEqual(len(actual), row['n_starts'])
            self.assertEqual({j['start_index'] for j in actual}, set(range(1, row['n_starts'] + 1)))
        self.assertTrue(all(j['status'] == 'PREPARED_NOT_RUN' for j in jobs))

    def test_cli_goat_manifest_defaults_to_1020_main_starts_without_launch(self):
        cli = load_cli()
        captured = []
        with patch.object(sys, 'argv', ['workflow.py', 'goat-manifest']), \
                patch.object(cli, 'writej', side_effect=lambda p, rows: captured.extend(rows)), \
                patch('subprocess.Popen', side_effect=AssertionError('No execution permitted')), \
                contextlib.redirect_stdout(io.StringIO()) as stream:
            cli.main()
        self.assertEqual(len(captured), 1020)
        self.assertEqual({r['scenario_id'] for r in captured}, {SCENARIO})
        self.assertIn('"submitted": false', stream.getvalue())


class V4IndependentHandoff(unittest.TestCase):
    def setup_fixture(self, species_id=pf.SPECIES):
        fixture = pf.ClusterPipelineTests()
        with patch.object(pf, 'SPECIES', species_id):
            fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        return fixture

    def prepared(self, fixture, species_id):
        fixture.completed_goat()
        pipeline.import_goat(fixture.root, fixture.row['job_id'], fixture.approved_review())
        return pipeline.prepare_prescreen(fixture.root, species_id, SCENARIO)

    def test_six_real_drug_structures_handoff_to_j3_only_synthetic_outputs(self):
        species = [s for s in readj(ROOT / 'registry/species.json') if s['compound_id'] in DRUGS]
        for row in species:
            with self.subTest(compound=row['compound_id']):
                sid = row['species_id']
                fixture = self.setup_fixture(sid)
                start = self.prepared(fixture, sid)
                fixture.completed_opt(fixture.root / start['jobs'][0])
                stage = pipeline.advance(fixture.root, sid, SCENARIO, 'prescreen')
                self.assertEqual(stage['status'], 'advanced')
                fixture.completed_opt(fixture.root / stage['prepared_jobs'][0])
                final = pipeline.advance(fixture.root, sid, SCENARIO, 'optfreq')
                self.assertEqual(final['status'], 'advanced')
                self.assertEqual(final['quantum_jobs_launched'], 0)
                self.assertEqual(set(final['prepared_jobs'][0]['j_jobs']), {'PBE_pcJ3'})
                manifest = pipeline.manifest(fixture.root, 'j_PBE_pcJ3', sid, SCENARIO)
                self.assertEqual(len(manifest), 1)
                self.assertEqual(pipeline.manifest(fixture.root, 'j_PBE0_pcJ2', sid, SCENARIO), [])
                inp = (fixture.root / manifest[0]['job_dir'] / 'job.inp').read_text(encoding='utf-8')
                self.assertIn('PBE pcJ-3', inp)
                self.assertIn('CPCM(DMSO)', inp)
                verified = pipeline.validate_ensemble_sources(fixture.root, sid, SCENARIO)
                self.assertTrue(verified['accepted'])
                self.assertFalse(verified['ensemble']['sampling_convergence_established'])
                fixture.doCleanups()

    def reviewed_optfreq_ensemble(self):
        f = self.setup_fixture()
        sid = f.row['species_id']
        start = self.prepared(f, sid)
        f.completed_opt(f.root / start['jobs'][0])
        stage = pipeline.advance(f.root, sid, SCENARIO, 'prescreen')
        job = f.root / stage['prepared_jobs'][0]
        f.completed_opt(job, negative=True)
        review = approve_negative_fixture(job, f.cfg)
        final = pipeline.advance(f.root, sid, SCENARIO, 'optfreq')
        self.assertEqual(final['status'], 'advanced')
        return f, sid, job, review, final

    def test_negative_exception_is_retained_in_ensemble_and_property_sources(self):
        f, sid, job, review, result = self.reviewed_optfreq_ensemble()
        row = result['prepared_jobs'][0]
        snapshot = row['scientific_review']
        self.assertEqual(snapshot['sha256'], sha(review))
        self.assertEqual(snapshot['scientific_status'], 'HUMAN_REVIEWED_APPROXIMATION')
        self.assertIn('No physical minimum', snapshot['limitations'])
        for path in [row['shielding_job'], *row['j_jobs'].values()]:
            self.assertEqual(readj(f.root / path / 'job.json')['source']['scientific_review'], snapshot)
        self.assertTrue(pipeline.validate_ensemble_sources(f.root, sid, SCENARIO)['accepted'])

    def test_changed_approved_review_cannot_be_silently_reused_downstream(self):
        f, sid, job, review, _ = self.reviewed_optfreq_ensemble()
        value = readj(review)
        value['limitations'] += ' Modified after handoff.'
        writej(review, value)
        self.assertTrue(gates.inspect_job(job, f.cfg)['accepted'])
        with self.assertRaisesRegex(ValueError, 'scientific review changed'):
            pipeline.validate_ensemble_sources(f.root, sid, SCENARIO)

    def test_historical_convergence_snapshot_rejects_changed_human_review(self):
        from clusterflow import convergence
        f, sid, job, review, _ = self.reviewed_optfreq_ensemble()
        verified = pipeline.validate_ensemble_sources(f.root, sid, SCENARIO)
        ensemble = verified['ensemble']
        frozen = f.root / 'SYNTHETIC_HISTORICAL_ENSEMBLE.json'
        writej(frozen, ensemble)
        result = {'metadata': {'species_id': sid, 'scenario_id': SCENARIO,
                  'j_method': 'PBE_pcJ3', 'protocol_hash': hash_protocol(f.cfg),
                  'ensemble_snapshot': {'target': ensemble},
                  'frozen_ensembles': {'target': {'path': frozen.name, 'sha256': sha(frozen)}},
                  'source': {'atom_map_sha256': ensemble['atom_map_sha256']}}}
        changed = readj(review)
        changed['limitations'] += ' Changed after historical wave was frozen.'
        writej(review, changed)
        self.assertTrue(gates.inspect_job(job, f.cfg)['accepted'])
        with self.assertRaisesRegex(ValueError, '[Ss]cientific.review|review.*changed'):
            convergence._recompute_snapshot(f.root, result, sid, SCENARIO, 'PBE_pcJ3')

    def historical_fixture(self, fixture, species_id, ensemble):
        frozen = fixture.root / 'SYNTHETIC_HISTORICAL_ENSEMBLE.json'
        writej(frozen, ensemble)
        return {'metadata': {'species_id': species_id, 'scenario_id': SCENARIO,
                'j_method': 'PBE_pcJ3', 'protocol_hash': hash_protocol(fixture.cfg),
                'ensemble_snapshot': {'target': ensemble},
                'frozen_ensembles': {'target': {'path': frozen.name, 'sha256': sha(frozen)}},
                'source': {'atom_map_sha256': ensemble['atom_map_sha256']}}}

    def test_historical_convergence_rejects_resealed_stored_G_change(self):
        from clusterflow import convergence
        f, sid, job, review, _ = self.reviewed_optfreq_ensemble()
        ensemble = pipeline.validate_ensemble_sources(f.root, sid, SCENARIO)['ensemble']
        ensemble['conformers'][0]['G_Eh'] += 0.01
        result = self.historical_fixture(f, sid, ensemble)
        with self.assertRaisesRegex(ValueError, 'stored Gibbs energy'):
            convergence._recompute_snapshot(f.root, result, sid, SCENARIO, 'PBE_pcJ3')

    def test_historical_convergence_rejects_old_property_review_after_reapproval(self):
        from clusterflow import convergence
        f, sid, job, review, _ = self.reviewed_optfreq_ensemble()
        ensemble = pipeline.validate_ensemble_sources(f.root, sid, SCENARIO)['ensemble']
        value = readj(review)
        value['limitations'] += ' Synthetic reapproval after property preparation.'
        writej(review, value)
        checked = gates.inspect_job(job, f.cfg)
        self.assertTrue(checked['accepted'])
        # Re-seal an edited historical object: shielding still points to the
        # earlier review, which must not be concealed by a new wave JSON hash.
        ensemble['conformers'][0]['scientific_review'] = pipeline._review_snapshot(checked)
        result = self.historical_fixture(f, sid, ensemble)
        with self.assertRaisesRegex(ValueError, 'shielding source scientific review'):
            convergence._recompute_snapshot(f.root, result, sid, SCENARIO, 'PBE_pcJ3')

    def test_historical_convergence_rejects_property_input_parent_hash_change(self):
        from clusterflow import convergence
        f, sid, job, review, _ = self.reviewed_optfreq_ensemble()
        ensemble = pipeline.validate_ensemble_sources(f.root, sid, SCENARIO)['ensemble']
        prop = f.root / ensemble['conformers'][0]['shielding_job']
        meta = readj(prop / 'job.json')
        meta['source']['geometry_sha256'] = '0' * 64
        writej(prop / 'job.json', meta)
        result = self.historical_fixture(f, sid, ensemble)
        with self.assertRaisesRegex(ValueError, 'source hash changed'):
            convergence._recompute_snapshot(f.root, result, sid, SCENARIO, 'PBE_pcJ3')

    def test_large_negative_review_is_not_a_hidden_numeric_cutoff(self):
        f = self.setup_fixture()
        prepared = self.prepared(f, f.row['species_id'])
        job = f.root / prepared['jobs'][0]
        f.completed_opt(job, negative=True)
        output = job / 'job.out'
        output.write_text(output.read_text(encoding='utf-8').replace('-20.00', '-500.00'), encoding='utf-8')
        receipt = readj(job / 'receipt.json')
        receipt['output_sha256']['job.out'] = sha(output)
        writej(job / 'receipt.json', receipt)
        self.assertFalse(gates.inspect_job(job, f.cfg)['accepted'])
        approve_negative_fixture(job, f.cfg)
        checked = gates.inspect_job(job, f.cfg)
        self.assertTrue(checked['accepted'], checked['errors'])
        self.assertEqual(checked['scientific_status'], 'HUMAN_REVIEWED_APPROXIMATION')
        self.assertFalse(checked['scientific_review']['physical_minimum_certified'])

    def test_diagnose_and_review_draft_cli_execute_no_qm_and_cannot_approve(self):
        f = self.setup_fixture()
        prepared = self.prepared(f, f.row['species_id'])
        cli = load_cli()
        cli.ROOT = f.root
        with patch.object(sys, 'argv', ['workflow.py', 'diagnose', '--job', prepared['jobs'][0]]), \
                contextlib.redirect_stdout(io.StringIO()) as stream:
            cli.main()
        self.assertIn('ADVISORY_ONLY_NOT_SCIENTIFIC_APPROVAL', stream.getvalue())
        self.assertIn('not_run_or_output_missing', stream.getvalue())


class V4IndependentMultiplet(unittest.TestCase):
    def setUp(self):
        self.fixture = mf.MultipletWorkflowTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)

    def test_experimental_field_cannot_be_overridden_by_nested_display_settings(self):
        self.fixture.panel['spectrum_settings'] = {'spectrum': {'field_MHz_H': 800}}
        result = self.fixture.run_panel()
        self.assertEqual(result['model']['field_MHz'], self.fixture.panel['experiment']['field_MHz'])

    def test_numeric_disagreement_is_advisory_not_autoapproval_or_deletion(self):
        self.fixture.panel['regions'][0]['lines'][0]['ppm'] += .05
        result = self.fixture.run_panel()
        self.assertEqual(result['status'], 'METRICS_COMPLETE_AWAITING_HUMAN_REVIEW')
        self.assertEqual(result['regions'][0]['recommendation'], 'REVIEW_REQUIRED')
        self.assertTrue(result['regions'][0]['metric_flags'])
        self.assertNotIn('approved', result)

    def test_missing_experiment_does_not_create_completed_validation(self):
        self.fixture.panel['experimental_artifacts'] = []
        with self.assertRaisesRegex(ValueError, 'Experimental source artifacts'):
            self.fixture.run_panel()


if __name__ == '__main__':
    unittest.main()
