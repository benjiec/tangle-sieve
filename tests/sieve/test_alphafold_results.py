import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import warnings
import zipfile

from sieve.alphafold_results import run
from tests.scripts.helpers import load_script


def zip_bytes(members):
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, 'w') as archive:
        for name, value in members.items():
            archive.writestr(name, value)
    return stream.getvalue()


class ResultsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.manifest = self.root/'manifest.tsv'
        self.out, self.err = io.StringIO(), io.StringIO()

    def manifest_rows(self, rows):
        with self.manifest.open('w', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t')
            writer.writerow(['fold_locator', 'fold_id', 'sequence_id'])
            writer.writerows(rows)

    def directory(self, name, members):
        folder = self.root/name
        folder.mkdir()
        for member, value in members.items():
            (folder/member).write_text(value)
        return str(folder)

    def run_results(self):
        return run(self.manifest, 'b911edee', self.out, self.err)

    def test_filter_order_and_no_coordinate_or_full_data_reads(self):
        folder = self.directory('second', {
            'fold_model_10.cif': 'not a CIF', 'fold_model_2.cif': 'not a CIF',
            'fold_summary_confidences_10.json': '{"iptm":0.8,"ptm":0.9}',
            'fold_summary_confidences_2.json': '{"iptm":0,"ptm":1}',
        })
        first = self.directory('first', {'fold_model_0.cif': 'invalid'})
        self.manifest_rows([(folder, '12345678', 'b911edee'),
                            ('/does/not/exist', '87654321', 'aaaaaaaa'),
                            (first, '87654321', 'b911edee')])
        original = Path.read_text
        def read(path, *args, **kwargs):
            self.assertTrue(path.name.endswith('.json'))
            return original(path, *args, **kwargs)
        with patch.object(Path, 'read_text', read):
            self.assertEqual(self.run_results(), 0)
        self.assertEqual(self.out.getvalue(), f'{folder}\t12345678\t2\t0\t1\n'
                         f'{folder}\t12345678\t10\t0.8\t0.9\n'
                         f'{first}\t87654321\t0\tNA\tNA\n')
        self.assertEqual(self.err.getvalue(), '')

    def test_archives_and_nested_selection(self):
        members = {'wanted/fold_model_0.cif': 'invalid',
                   'wanted/fold_summary_confidences_0.json': '{"iptm":null,"ptm":0.7}',
                   'other/fold_model_0.cif': 'invalid', 'bad.zip': 'invalid'}
        for nested in (False, True):
            with self.subTest(nested=nested):
                archive = self.root/'batch.zip'
                archive.write_bytes(zip_bytes({'inner.zip': zip_bytes(members)}) if nested else zip_bytes(members))
                locator = f'{archive}::' + ('inner.zip!/' if nested else '') + 'wanted'
                self.manifest_rows([(locator, '12345678', 'b911edee')])
                self.out = io.StringIO()
                self.assertEqual(self.run_results(), 0)
                self.assertEqual(self.out.getvalue(), f'{locator}\t12345678\t0\tNA\t0.7\n')

    def test_missing_field(self):
        folder = self.directory('one', {'fold_model_0.cif': '',
                                      'fold_summary_confidences_0.json': '{"iptm":0.1}'})
        self.manifest_rows([(folder, '12345678', 'b911edee')])
        self.assertEqual(self.run_results(), 0)
        self.assertTrue(self.out.getvalue().endswith('\t0.1\tNA\n'))

    def test_model_and_fold_errors_continue(self):
        members = {}
        for number, summary in enumerate(['broken', '{"ptm":2}', '{"iptm":true}', '[]', '{}']):
            members[f'fold_model_{number}.cif'] = ''
            members[f'fold_summary_confidences_{number}.json'] = summary
        folder = self.directory('bad-models', members)
        good = self.directory('good', {'fold_model_0.cif': ''})
        self.manifest_rows([('/missing', '12345678', 'b911edee'),
                            (folder, '12345678', 'b911edee'), (good, '87654321', 'b911edee')])
        self.assertEqual(self.run_results(), 1)
        self.assertEqual(len(self.err.getvalue().splitlines()), 5)
        self.assertEqual(self.out.getvalue(), f'{folder}\t12345678\t4\tNA\tNA\n'
                         f'{good}\t87654321\t0\tNA\tNA\n')

    def test_duplicate_models_and_summary_members(self):
        folder = self.directory('duplicate', {'a_model_0.cif': '', 'b_model_0.cif': ''})
        archive = self.root/'duplicate.zip'
        archive.write_bytes(zip_bytes({'fold_model_0.cif': '', 'fold_summary_confidences_0.json': '{}'}))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            with zipfile.ZipFile(archive, 'a') as stream:
                stream.writestr('fold_summary_confidences_0.json', '{}')
        self.manifest_rows([(folder, '12345678', 'b911edee'), (f'{archive}::.', '87654321', 'b911edee')])
        self.assertEqual(self.run_results(), 1)
        self.assertEqual(self.out.getvalue(), '')
        self.assertIn('duplicate model numbers', self.err.getvalue())
        self.assertIn('duplicate summary file', self.err.getvalue())

    def test_no_match_and_header_only(self):
        for rows in ([], [('/missing', '12345678', 'aaaaaaaa')]):
            self.manifest_rows(rows)
            self.assertEqual(self.run_results(), 0)
        self.assertEqual(self.out.getvalue(), '')
        self.assertEqual(self.err.getvalue(), '')

    def test_invalid_manifest_fails_before_output(self):
        for text in ('', 'wrong\n', 'fold_locator\tfold_id\tsequence_id\nloc\t12345678\n',
                     'fold_locator\tfold_id\tsequence_id\nloc\tbad\tb911edee\n'):
            self.manifest.write_text(text)
            with self.assertRaises(ValueError):
                self.run_results()
            self.assertEqual(self.out.getvalue(), '')
        self.manifest_rows([])
        with self.assertRaises(ValueError):
            run(self.manifest, 'invalid')

    def test_cli_only_data_stdout(self):
        script = load_script(Path(__file__).resolve().parents[2]/'scripts'/'alphafold-list-results.py')
        folder = self.directory('with spaces', {'fold_model_0.cif': ''})
        self.manifest_rows([(folder, '12345678', 'b911edee')])
        with contextlib.redirect_stdout(self.out), contextlib.redirect_stderr(self.err):
            self.assertEqual(script.main([str(self.manifest), 'b911edee']), 0)
        self.assertEqual(self.out.getvalue(), f'{folder}\t12345678\t0\tNA\tNA\n')

    def test_highest_cli_modes_across_folds(self):
        script = load_script(Path(__file__).resolve().parents[2]/'scripts'/'alphafold-list-results.py')
        first = self.directory('first', {'fold_model_0.cif': '',
                                      'fold_summary_confidences_0.json': '{"iptm":0.9,"ptm":0.1}'})
        second = self.directory('second', {'fold_model_10.cif': '', 'fold_model_2.cif': '',
                                         'fold_summary_confidences_10.json': '{"iptm":0.8,"ptm":0.8}',
                                         'fold_summary_confidences_2.json': '{"iptm":0.2,"ptm":0.99}'})
        self.manifest_rows([(first, '12345678', 'b911edee'), (second, '87654321', 'b911edee')])
        cases = [([], [(first, '0'), (second, '2'), (second, '10')]),
                 (['--highest'], [(second, '10')]),
                 (['--highest', 'ipTM8_pTM2'], [(second, '10')]),
                 (['--highest', 'ipTM'], [(first, '0')]),
                 (['--highest', 'pTM'], [(second, '2')])]
        for flags, expected in cases:
            with self.subTest(flags=flags):
                out, err = io.StringIO(), io.StringIO()
                with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                    self.assertEqual(script.main([str(self.manifest), 'b911edee', *flags]), 0)
                rows = list(csv.reader(io.StringIO(out.getvalue()), delimiter='\t'))
                self.assertEqual([(r[0], r[2]) for r in rows], expected)
                self.assertTrue(all(len(r) == 5 for r in rows))
                self.assertEqual(err.getvalue(), '')
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as caught:
            script.main([str(self.manifest), 'b911edee', '--highest', 'invalid'])
        self.assertEqual(caught.exception.code, 2)

    def test_highest_missing_values_zero_and_no_eligible_rows(self):
        folder = self.directory('missing', {'fold_model_0.cif': '', 'fold_model_1.cif': '',
                                          'fold_summary_confidences_0.json': '{"iptm":0}',
                                          'fold_summary_confidences_1.json': '{"ptm":0}'})
        self.manifest_rows([(folder, '12345678', 'b911edee')])
        for metric, expected in [('ipTM', '\t0\t0\tNA\n'), ('pTM', '\t1\tNA\t0\n'), ('ipTM8_pTM2', '')]:
            out = io.StringIO()
            self.assertEqual(run(self.manifest, 'b911edee', out, self.err, highest=metric), 0)
            if expected:
                self.assertTrue(out.getvalue().endswith(expected))
            else:
                self.assertEqual(out.getvalue(), '')
        self.manifest_rows([])
        out = io.StringIO()
        self.assertEqual(run(self.manifest, 'b911edee', out, self.err, highest='ipTM'), 0)
        self.assertEqual(out.getvalue(), '')

    def test_highest_ties_preserve_manifest_then_numeric_order(self):
        # Both weighted scores are exactly 0.24; decimal scoring avoids float tie noise.
        first = self.directory('z-first', {'fold_model_10.cif': '', 'fold_model_2.cif': '',
                                         'fold_summary_confidences_10.json': '{"iptm":0.1,"ptm":0.8}',
                                         'fold_summary_confidences_2.json': '{"iptm":0.2,"ptm":0.4}'})
        second = self.directory('a-second', {'fold_model_0.cif': '',
                                           'fold_summary_confidences_0.json': '{"iptm":0.2,"ptm":0.4}'})
        self.manifest_rows([(first, '12345678', 'b911edee'), (second, '87654321', 'b911edee')])
        self.assertEqual(run(self.manifest, 'b911edee', self.out, self.err, highest='ipTM8_pTM2'), 0)
        self.assertEqual(self.out.getvalue(), f'{first}\t12345678\t2\t0.2\t0.4\n')

    def test_highest_reports_errors_and_returns_best_valid_row(self):
        folder = self.directory('valid', {'fold_model_0.cif': '',
                                         'fold_summary_confidences_0.json': '{"iptm":0.6,"ptm":0.7}'})
        self.manifest_rows([(folder, '12345678', 'b911edee'), ('/missing', '87654321', 'b911edee')])
        self.assertEqual(run(self.manifest, 'b911edee', self.out, self.err, highest='ipTM'), 1)
        self.assertEqual(self.out.getvalue(), f'{folder}\t12345678\t0\t0.6\t0.7\n')
        self.assertIn('ERROR:', self.err.getvalue())
