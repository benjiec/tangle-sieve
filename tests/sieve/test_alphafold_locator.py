import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

from sieve.alphafold_files import fold_locator
from sieve.alphafold_manifest import run
from sieve.alphafold_pdockq2_patches import input_models
from tests.sieve.test_alphafold_manifest import archive_bytes, cif
from tests.sieve.test_alphafold_pdockq2_patches import CIF, DATA
from tests.scripts.helpers import load_script


class LocatorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.path = self.root / 'batch with spaces.zip'
        self.members = {'wanted/fold_model_0.cif': CIF,
                        'wanted/fold_full_data_0.json': json.dumps(DATA),
                        'other/fold_model_0.cif': 'invalid unrelated model',
                        'broken.zip': 'invalid unrelated archive'}

    def test_exact_selection_and_nested_root(self):
        for prefix in ('wanted', 'inner.zip!/wanted', 'inner.zip!/.'):
            with self.subTest(prefix=prefix):
                members = self.members
                if prefix.endswith('!/.'):
                    members = {Path(k).name: v for k, v in self.members.items() if k.startswith('wanted/')}
                self.path.write_bytes(archive_bytes({'inner.zip': archive_bytes(members)})
                                      if '!/' in prefix else archive_bytes(members))
                locator = fold_locator(str(self.path), prefix)
                models = list(input_models(locator))
                self.assertEqual(models, [(locator, 0, CIF, DATA)])

    def test_missing_fold_data_and_duplicate_model(self):
        self.path.write_bytes(archive_bytes(self.members))
        with self.assertRaisesRegex(ValueError, 'no model CIF'):
            list(input_models(f'{self.path}::missing'))
        with self.assertRaisesRegex(ValueError, 'nested archive'):
            list(input_models(f'{self.path}::absent.zip!/wanted'))
        with self.assertRaisesRegex(ValueError, 'only valid with CIF'):
            list(input_models(f'{self.path}::wanted', 'data.json'))
        for locator in (f'{self.path}::', f'{self.path}::inner.zip!/'):
            with self.assertRaises(ValueError):
                list(input_models(locator))
        del self.members['wanted/fold_full_data_0.json']
        self.path.write_bytes(archive_bytes(self.members))
        with self.assertRaisesRegex(ValueError, 'exactly one archive member'):
            list(input_models(f'{self.path}::wanted'))
        self.members['wanted/other_model_0.cif'] = CIF
        self.path.write_bytes(archive_bytes(self.members))
        with self.assertRaisesRegex(ValueError, 'unique model numbers'):
            list(input_models(f'{self.path}::wanted'))

    def test_duplicate_data_rejected(self):
        self.path.write_bytes(archive_bytes(self.members))
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            with zipfile.ZipFile(self.path, 'a') as archive:
                archive.writestr('wanted/fold_full_data_0.json', json.dumps(DATA))
        with self.assertRaisesRegex(ValueError, 'exactly one archive member'):
            list(input_models(f'{self.path}::wanted'))

    def test_manifest_to_patch_cli(self):
        script = load_script(Path(__file__).resolve().parents[2] / 'scripts' / 'alphafold-pdockq2-patches.py')
        # Entity metadata for the two alanine chains in the patch fixture.
        combined = cif(('A', 'A')) + CIF[CIF.index('loop_'):]
        members = {'fold_model_0.cif': combined, 'fold_full_data_0.json': json.dumps(DATA)}
        folder = self.root / 'extracted'
        folder.mkdir()
        for name, text in members.items():
            (folder / name).write_text(text)
        self.path.write_bytes(archive_bytes({'inner.zip': archive_bytes(members)}))
        fasta = self.root / 'candidates.faa'
        fasta.write_text('>alanine full header\nA\n')
        for source in (folder, self.path):
            manifest = self.root / 'manifest.tsv'
            self.assertEqual(run(source, fasta, manifest, io.StringIO()), 0)
            with manifest.open() as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            locator = next(r['fold_locator'] for r in rows if r['source_path'] == str(source))
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                self.assertEqual(script.main([locator, '--no-dssp', '--min-residues-per-chain', '1', '--min-contacts', '1']), 0)
            self.assertIn('Model 0', output.getvalue())
