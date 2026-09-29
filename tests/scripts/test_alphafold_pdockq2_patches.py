import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest

from tests.scripts.helpers import load_script
from tests.sieve.test_alphafold_pdockq2_patches import CIF, DATA

SCRIPT = load_script(str(Path(__file__).resolve().parents[2] / 'scripts' / 'alphafold-pdockq2-patches.py'))


class PatchCommandTests(unittest.TestCase):
    def test_cli_success_and_failure_preserves_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cif, data, output = root/'model.cif', root/'data.json', root/'patches.tsv'
            cif.write_text(CIF)
            data.write_text(json.dumps(DATA))
            args = [str(cif), '--full-data', str(data), '--output', str(output),
                    '--min-residues-per-chain', '1', '--min-contacts', '1']
            self.assertEqual(SCRIPT.main(args), 0)
            before = output.read_text()
            rows = list(csv.DictReader(io.StringIO(before), delimiter='\t'))
            self.assertEqual([r['scope'] for r in rows], ['full', 'patch'])
            self.assertEqual(rows[1]['interface residues 1'], '1')
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                SCRIPT.main(args + ['--chains', 'A', 'X'])
            self.assertEqual(output.read_text(), before)
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                SCRIPT.main([str(cif), '--full-data', str(data), '--output', str(cif)])
            self.assertEqual(cif.read_text(), CIF)

    def test_help_mentions_specification(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output), self.assertRaises(SystemExit) as error:
            SCRIPT.main(['--help'])
        self.assertEqual(error.exception.code, 0)
        self.assertIn('alphafold-pdockq2-patches.md', output.getvalue())


class SummaryTests(unittest.TestCase):
    def row(self, model, a, b, score, scope='patch'):
        return {'source': 'example.zip', 'model': model, 'scope': scope,
                'chain 1': a, 'chain 2': b, 'pDockQ2 max': score,
                'pDockQ2 1 to 2': score, 'pDockQ2 2 to 1': score,
                'contact distance min': 5.0, 'contact distance max': 7.0,
                'contact distance mean': 6.0, 'contact distance median': 6.0,
                'residue count 1': 1, 'residue count 2': 1,
                **{f'nearest contact distance {stat} {side}': 5.0
                   for side in (1, 2) for stat in ('min', 'max', 'median', 'mean')},
                'contact count': 10, 'graph residues': f'{a}:1,{b}:1',
                'patch id': f'p{model}', 'status': 'leaf',
                'region 1': '1-3', 'region 2': '60-65',
                'interface residues 1': str(model * 10 + 1),
                'interface residues 2': str(model * 10 + 2)}

    def test_top_five_last_chain_excludes_other_pairs_and_reference(self):
        rows = [self.row(m, 'A', 'C', m / 10) for m in range(7)]
        rows += [self.row(0, 'A', 'B', 1.0), self.row(0, 'B', 'C', 1.0, 'full')]
        targets = {('example.zip', m): 'C' for m in range(7)}
        output = io.StringIO()
        SCRIPT.write_summary(rows, targets, output)
        text = output.getvalue().split('Top 5 patches', 1)[1]
        self.assertEqual(text.count('Region 1'), 5)
        self.assertIn('1. Model 6', text)
        self.assertIn('min 5.000 | max 7.000 | median 6.000 | average 6.000', text)
        self.assertNotIn('Model 0 |', text)
        self.assertNotIn('A–B', text)
        self.assertNotIn('1.000000', text)

    def test_last_chain_uses_cif_order_not_lexical_order(self):
        from sieve.alphafold_pdockq2_patches import analyze
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cif, data = root/'model.cif', root/'data.json'
            cif.write_text(CIF.replace('ALA A 1', 'ALA Z 1').replace('90 1 A 1', '90 1 Z 1'))
            data.write_text(json.dumps(dict(DATA, token_chain_ids=['Z', 'B'])))
            targets = {}
            analyze(cif, data, summary_targets=targets)
            self.assertEqual(list(targets.values()), ['B'])

    def test_no_matches_and_separate_sources(self):
        output = io.StringIO()
        SCRIPT.write_summary([], {('first.zip', 0): 'C', ('second.zip', 0): 'Z'}, output)
        self.assertEqual(output.getvalue().count('No qualifying patches'), 2)

    def test_overlap_threshold_containment_and_chain_identity(self):
        from sieve.alphafold_pdockq2_patches import distinct_patches
        def patch(a, b):
            row = self.row(0, 'A', 'B', 0.5)
            row.update({'interface residues 1': ','.join(map(str, a)),
                        'interface residues 2': ','.join(map(str, b))})
            return row
        first = patch(range(1, 6), range(1, 6))
        exactly_ten = patch([1, 10, 11, 12, 13], range(10, 15))
        twenty = patch([1, 2, 20, 21, 22], range(20, 25))
        contained = patch([1], [1])
        self.assertEqual(distinct_patches([first, first, contained, twenty, exactly_ten]), [first, exactly_ten])
        self.assertEqual(len(distinct_patches([patch([1], [2]), patch([2], [1])])), 2)

    def test_stdout_only_creates_no_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cif, data = root/'model.cif', root/'data.json'
            cif.write_text(CIF)
            data.write_text(json.dumps(DATA))
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                self.assertEqual(SCRIPT.main([str(cif), '--full-data', str(data),
                                             '--min-residues-per-chain', '1', '--min-contacts', '1']), 0)
            self.assertIn('Model unnumbered', output.getvalue())
            self.assertEqual(set(root.iterdir()), {cif, data})

    def test_full_pairs_precede_last_chain_patches(self):
        rows = [self.row(0, 'A', 'B', 0.9, 'full'),
                self.row(0, 'A', 'C', 0.6, 'full'), self.row(0, 'A', 'C', 0.8)]
        output = io.StringIO()
        SCRIPT.write_summary(rows, {('example.zip', 0): 'C'}, output)
        full, patches = output.getvalue().split('Top 5 patches', 1)
        self.assertIn('A–B', full)
        self.assertIn('A→B 0.900000', full)
        self.assertIn('A–C', full)
        self.assertNotIn('All interchain contact pairs', full)
        self.assertNotIn('Nearest contact', full)
        self.assertIn('Contacts: 10', full)
        self.assertNotIn('A–B', patches)
        self.assertIn('max pDockQ2 0.800000', patches)

    def test_full_no_contacts_prints_zero_counts(self):
        row = self.row(0, 'A', 'B', 0, 'full')
        for key in row:
            if key.startswith(('contact distance', 'nearest contact distance')):
                row[key] = ''
        row.update({'contact count': 0, 'residue count 1': 0, 'residue count 2': 0})
        output = io.StringIO()
        SCRIPT.write_summary([row], {('example.zip', 0): 'B'}, output)
        self.assertIn('Contacts: 0 | residues A: 0, B: 0', output.getvalue())
        self.assertIn('No qualifying patches', output.getvalue())
