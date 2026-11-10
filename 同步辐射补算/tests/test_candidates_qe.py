"""Software tests only; no fixture is a DFT calculation or scientific result."""
from pathlib import Path
import hashlib
import json
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from build_candidates import build, load_cif, read_json, sha
from qe_tasks import collect_one, dipole_clearance, dipole_region, parse_output, prepare, settings


def output_fixture(structure, force=.0001, electronic=True, ionic=True, done=True):
    lines = ['TEST ONLY PROGRAM PARSER FIXTURE; not an electronic-structure calculation']
    lines += ['convergence has been achieved in 18 iterations' if electronic else 'convergence NOT achieved after 250 iterations']
    lines += ['!    total energy = -100.123456 Ry', 'Forces acting on atoms (cartesian axes, Ry/au):']
    for i in range(len(structure)):
        lines.append(f'atom {i+1:5d} type 1 force = {force:.9f} 0.000000000 0.000000000')
    lines += ['Total force = 0.0009']
    if ionic:
        lines += ['bfgs converged in 5 scf cycles and 4 bfgs steps', 'End of BFGS Geometry Optimization']
    lines += ['Begin final coordinates', 'ATOMIC_POSITIONS (angstrom)']
    for site in structure:
        lines.append(site.specie.symbol+' '+' '.join(f'{v:.12f}' for v in site.coords))
    lines += ['End final coordinates']
    if done:
        lines += ['JOB DONE.']
    return '\n'.join(lines)+'\n'


class CandidateQETests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='XAS_CANDIDATE_TEST_ONLY_')
        cls.root = Path(cls.temp.name)
        (cls.root/'structures_original').mkdir()
        source = ROOT/'structures_original/CoFe-LDH.cif'
        if not source.exists():
            source = ROOT.parent/'XAS_CIF_HPC_WORKFLOW_20261029/structures/CoFe-LDH.cif'
        shutil.copy2(source, cls.root/'structures_original/CoFe-LDH.cif')
        pseudo_dir = cls.root/'software/pseudo'; pseudo_dir.mkdir(parents=True)
        pseudos = {}
        for element in ['H', 'C', 'O', 'Co', 'Fe', 'W']:
            filename = element+'.TEST_ONLY.UPF'
            content = ('TEST ONLY NON-PHYSICAL PSEUDOPOTENTIAL PLACEHOLDER '+element).encode()
            (pseudo_dir/filename).write_bytes(content)
            pseudos[element] = dict(file=filename, sha256=hashlib.sha256(content).hexdigest())
        cls.config = cls.root/'config.json'
        cls.config.write_text(json.dumps({'dft': {'pseudopotentials': pseudos}}), encoding='utf-8')
        cls.models = build(cls.config)
        cls.jobs = prepare(cls.config)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_12_crossed_models_and_six_oxygen_initial_shell(self):
        self.assertEqual(len(self.models), 12)
        factors = {(j['anchor'], j['background'], j['adsorbate']) for j in self.models}
        self.assertEqual(len(factors), 12)
        for job in self.models:
            metadata = read_json(self.root/job['path']/'metadata.json')
            self.assertEqual(metadata['status'], 'candidate_not_DFT_optimized')
            self.assertTrue(metadata['initial_coordination']['all_heavy_neighbours_are_O'])
            self.assertEqual(metadata['initial_coordination']['oxygen_count'], 6)
            self.assertTrue(metadata['initial_clash_check']['passed'])

    def test_atom_identity_round_trip_and_explicit_stoichiometry(self):
        for job in self.models:
            metadata = read_json(self.root/job['path']/'metadata.json')
            structure, _ = load_cif(self.root/job['path']/'structure.cif')
            self.assertEqual([a['id'] for a in metadata['atom_map']], [s.label for s in structure])
            self.assertEqual([a['element'] for a in metadata['atom_map']], [s.specie.symbol for s in structure])
            expected_h = 32 - int(job['background'] == 'Co_deprotonated') + 4*int(job['adsorbate'] != 'none')
            self.assertEqual(metadata['composition']['H'], expected_h)
            self.assertEqual(metadata['composition']['O'], 32+int(job['adsorbate'] != 'none'))
            self.assertEqual(metadata['composition'].get('C', 0), 2*int(job['adsorbate'] != 'none'))
            if metadata['ligand_exchange']:
                self.assertTrue(metadata['ligand_exchange']['retained_original_O'])

    def test_qe_input_is_spin_polarized_and_no_atomic_constraints(self):
        text = (self.root/self.jobs[1]['path']/'qe.in').read_text()
        self.assertIn('ecutwfc=90', text)
        self.assertIn('ecutrho=1080', text)
        self.assertIn('nspin=2', text)
        self.assertIn('dipfield=.true.', text)
        self.assertIn('eamp=0.0', text)
        self.assertNotIn('HUBBARD', text)
        atom_lines = text.split('ATOMIC_POSITIONS angstrom\n')[1].split('K_POINTS')[0].splitlines()
        self.assertTrue(all(len(line.split()) == 4 for line in atom_lines))

    def test_qe_parser_requires_all_convergence_conditions(self):
        structure, _ = load_cif(self.root/self.models[0]['path']/'structure.cif')
        dft = settings(read_json(self.config))
        xyz, metrics = parse_output(output_fixture(structure), structure, dft)
        self.assertEqual(xyz.shape, (len(structure), 3))
        self.assertEqual(metrics['failures'], [])
        for kwargs in [dict(done=False), dict(ionic=False), dict(electronic=False), dict(force=.01)]:
            _, bad = parse_output(output_fixture(structure, **kwargs), structure, dft)
            self.assertTrue(bad['failures'], kwargs)
        _, bad = parse_output(output_fixture(structure)+'convergence NOT achieved after 250 iterations\n', structure, dft)
        self.assertTrue(bad['failures'])

    def test_collect_retains_atom_identity_and_provenance(self):
        job = self.jobs[1]
        directory = self.root/job['path']
        structure, _ = load_cif(self.root/job['candidate_path']/'structure.cif')
        (directory/'qe.out').write_text(output_fixture(structure), encoding='utf-8')
        result = collect_one(self.root, job, {'returncode': 0, 'fixture_only': True})
        self.assertEqual(result['status'], 'optimized')
        relaxed, _ = load_cif(self.root/result['relaxed_cif'])
        self.assertEqual([s.label for s in relaxed], [s.label for s in structure])
        self.assertEqual(result['relaxed_coordination']['oxygen_count'], 6)
        self.assertEqual(result['relaxed_sha256'], sha(self.root/result['relaxed_cif']))
        self.assertTrue(result['fixture_only'])

    def test_changed_settings_refuse_stale_input_reuse(self):
        config = read_json(self.config)
        config['dft']['ecutwfc_Ry'] = 95
        alternate = self.root/'changed.json'; alternate.write_text(json.dumps(config), encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'overwrite|changed'):
            prepare(alternate)

    def test_pseudo_digest_change_rejected(self):
        path = self.root/'software/pseudo/W.TEST_ONLY.UPF'; original = path.read_bytes()
        try:
            path.write_bytes(b'changed test only pseudo')
            with self.assertRaisesRegex(ValueError, 'digest mismatch'):
                prepare(self.config)
        finally:
            path.write_bytes(original)

    def test_dipole_window_stays_in_empty_vacuum_for_all_candidates(self):
        for job in self.jobs:
            structure, _ = load_cif(self.root/job['candidate_path']/'structure.cif')
            region = job['dipole_region']
            self.assertTrue(dipole_clearance(structure, structure.cart_coords, region)['passed'])
            moved = structure.cart_coords.copy()
            moved[0, 2] = (region['emaxpos']+.5*region['eopreg'])*region['cell_c_A']
            self.assertFalse(dipole_clearance(structure, moved, region)['passed'])
        structure, _ = load_cif(self.root/self.models[0]['path']/'structure.cif')
        dft = settings(read_json(self.config))
        dft['dipole_vacuum_margin_A'] = 20.
        with self.assertRaisesRegex(ValueError, 'vacuum region'):
            dipole_region(structure, dft)


if __name__ == '__main__':
    unittest.main()
