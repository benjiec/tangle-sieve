import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

from sieve.alphafold_pdockq2_patches import analyze, input_models, model_confidence, write_summary
from tests.sieve.test_alphafold_manifest import archive_bytes
from tests.sieve.test_alphafold_pdockq2_patches import CIF, DATA
from tests.scripts.helpers import load_script


class ConfidenceTests(unittest.TestCase):
    def test_all_input_types(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            folder = root / 'fold'
            folder.mkdir()
            members = {}
            for number in (0, 2):
                members.update({f'fold_model_{number}.cif': CIF,
                                f'fold_full_data_{number}.json': json.dumps(DATA),
                                f'fold_summary_confidences_{number}.json': json.dumps({'iptm': 0.2 + number / 10, 'ptm': 0.8})})
            for name, text in members.items():
                (folder / name).write_text(text)
            archive = root / 'fold.zip'
            archive.write_bytes(archive_bytes(members))
            outer = root / 'outer.zip'
            outer.write_bytes(archive_bytes({'inner.zip': archive_bytes(members)}))
            for source, data, expected in [(folder, None, [0, 2]), (archive, None, [0, 2]),
                                           (f'{archive}::.', None, [0, 2]),
                                           (f'{outer}::inner.zip!/.', None, [0, 2]),
                                           (folder/'fold_model_0.cif', folder/'fold_full_data_0.json', [0])]:
                with self.subTest(source=source):
                    scores = {}
                    models = list(input_models(source, data, model_scores=scores))
                    self.assertEqual([m[1] for m in models], expected)
                    for src, number, _, _ in models:
                        self.assertEqual(scores[(src, number)], {'iptm': 0.2 + number / 10, 'ptm': 0.8})
            renamed = folder/'renamed.cif'
            renamed.write_text(CIF)
            scores = {}
            list(input_models(renamed, folder/'fold_full_data_0.json', model_scores=scores))
            self.assertEqual(scores[(str(renamed), '')]['ptm'], 0.8)

    def test_invalid_and_missing_scores(self):
        for value in (-1, 1.1, float('nan'), float('inf'), True, '0.5', []):
            with self.subTest(value=value), self.assertRaises(ValueError):
                model_confidence({'iptm': value})
        for value in ([], None):
            with self.assertRaises(ValueError):
                model_confidence(value)
        self.assertEqual(model_confidence({'iptm': None, 'ptm': 0}), {'iptm': None, 'ptm': 0})

    def test_summary_without_interfaces(self):
        out = io.StringIO()
        write_summary([], {('fold', 0): 'A', ('fold', 1): 'A'}, out,
                      model_scores={('fold', 0): {'iptm': 0, 'ptm': 1}})
        self.assertIn('Model 0 | ipTM 0.000000 | pTM 1.000000', out.getvalue())
        self.assertIn('Model 1 | ipTM n/a (unavailable) | pTM n/a (unavailable)', out.getvalue())

    def test_fold_hash_consistency_and_missing_metadata(self):
        from tests.sieve.test_alphafold_manifest import cif
        from sieve.alphafold_manifest import fold_id
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            for number in (0, 1):
                (root/f'fold_model_{number}.cif').write_text(cif(('A', 'A')) + CIF[CIF.index('loop_'):])
                (root/f'fold_full_data_{number}.json').write_text(json.dumps(DATA))
            (root/'fold_job_request.json').write_text('[{"modelSeeds":["42"]}]')
            ids = {}
            analyze(root, fold_ids=ids)
            self.assertEqual(ids, {str(root): fold_id(42, ['A', 'A'])})
            (root/'fold_model_1.cif').write_text(cif(('A', 'C')) + CIF[CIF.index('loop_'):])
            with self.assertRaisesRegex(ValueError, 'composition differs'):
                analyze(root, fold_ids={})
            (root/'fold_model_1.cif').write_text(CIF)
            ids = {}
            analyze(root, fold_ids=ids)
            self.assertIsNone(ids[str(root)])

    def test_cli_structure_opt_in(self):
        script = load_script(Path(__file__).resolve().parents[2] / 'scripts' / 'alphafold-pdockq2-patches.py')
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'fold_model_0.cif').write_text(CIF)
            (root/'fold_full_data_0.json').write_text(json.dumps(DATA))
            summary = root/'fold_summary_confidences_0.json'
            summary.write_text(json.dumps({'iptm': 0.71, 'ptm': 0.82}))
            args = [str(root), '--no-dssp', '--min-residues-per-chain', '1', '--min-contacts', '1']
            for structure in (False, True):
                out = io.StringIO()
                with contextlib.redirect_stdout(out):
                    self.assertEqual(script.main(args + (['--structure'] if structure else [])), 0)
                self.assertEqual('Chain secondary structure' in out.getvalue(), structure)
                self.assertIn('Model 0 | ipTM 0.710000 | pTM 0.820000', out.getvalue())
                self.assertIn('patch sequence:', out.getvalue())
            summary.unlink()
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                self.assertEqual(script.main(args), 0)
            self.assertIn('ipTM n/a (unavailable) | pTM n/a (unavailable)', out.getvalue())
