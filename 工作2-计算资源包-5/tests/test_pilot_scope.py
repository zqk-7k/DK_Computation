import json
from pathlib import Path
import tempfile
import unittest
from scripts.packed_worker import selected_species

ROOT=Path(__file__).resolve().parents[1]

class PilotScopeTests(unittest.TestCase):
    def test_pilot_includes_reference_and_four_representatives(self):
        selected=selected_species(ROOT,ROOT/'config/pilot_scope_20261014.json')
        self.assertEqual(len(selected),5)
        self.assertIn('TMS_DMSO',selected)
        self.assertIn('F4_D04_rep_R__DMSO',selected)
        self.assertNotIn('TMS_DMSO',selected_species(ROOT))

    def test_invalid_scope_never_falls_back_to_all_products(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'scope.json'
            for ids in ([],['UNKNOWN'],['TMS_DMSO','TMS_DMSO'],'TMS_DMSO'):
                path.write_text(json.dumps({'species_ids':ids}))
                with self.assertRaises(ValueError):selected_species(ROOT,path)
