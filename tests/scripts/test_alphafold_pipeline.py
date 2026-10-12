import contextlib
import csv
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from sieve.alphafold_result_format import ResultRecord, write_result, read_patch_inputs
from tests.scripts.helpers import load_script
from tests.sieve.test_alphafold_pdockq2_patches import CIF, DATA
from tests.sieve.test_alphafold_results import zip_bytes

SCRIPTS = Path(__file__).resolve().parents[2]/'scripts'
SCRIPT = load_script(SCRIPTS/'alphafold-pdockq2-patches.py')


def record(locator, model):
    stream = io.StringIO()
    write_result(csv.writer(stream, delimiter='\t', lineterminator='\n'),
                 ResultRecord(locator, '12345678', model, None, 0.99))
    return stream.getvalue()


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.folder = self.root/'fold with spaces'
        self.folder.mkdir()
        self.members = {'fold_model_0.cif': CIF, 'fold_full_data_0.json': json.dumps(DATA),
                        'fold_summary_confidences_0.json': '{"iptm":0.2,"ptm":0.3}',
                        'fold_model_1.cif': 'invalid CIF',
                        'fold_summary_confidences_1.json': 'invalid JSON'}
        for name, text in self.members.items():
            (self.folder/name).write_text(text)
        self.args = ['-', '--no-dssp', '--min-residues-per-chain', '1', '--min-contacts', '1']

    def call(self, stdin, extra=()):
        out = io.StringIO()
        with patch('sys.stdin', io.StringIO(stdin)), contextlib.redirect_stdout(out):
            self.assertEqual(SCRIPT.main(self.args + list(extra)), 0)
        return out.getvalue()

    def test_selected_model_skips_other_files_directory_and_nested_zip(self):
        archive = self.root/'outer.zip'
        archive.write_bytes(zip_bytes({'inner.zip': zip_bytes(self.members)}))
        for locator in (str(self.folder), f'{archive}::inner.zip!/.'):
            with self.subTest(locator=locator):
                output = self.call(record(locator, 0))
                self.assertIn('Model 0 | ipTM 0.200000 | pTM 0.300000', output)
                self.assertNotIn('Model 1', output)
                self.assertNotIn('pTM 0.990000', output)

    def test_locator_only_matches_direct_and_multiple_records_separate(self):
        (self.folder/'fold_model_1.cif').write_text(CIF)
        (self.folder/'fold_full_data_1.json').write_text(json.dumps(DATA))
        (self.folder/'fold_summary_confidences_1.json').write_text('{}')
        direct = io.StringIO()
        with contextlib.redirect_stdout(direct):
            SCRIPT.main([str(self.folder), *self.args[1:]])
        self.assertEqual(self.call(str(self.folder)+'\n'), direct.getvalue())
        out = self.call(record(str(self.folder), 1) + record(str(self.folder), 0))
        self.assertEqual(out.count('Fold ID:'), 2)
        self.assertLess(out.index('Model 1'), out.index('Model 0'))

    def test_empty_stdin_no_output_or_file_change(self):
        output = self.root/'out.tsv'
        output.write_text('original')
        self.assertEqual(self.call('\n\n', ['--output', str(output)]), '')
        self.assertEqual(output.read_text(), 'original')

    def test_invalid_input_or_missing_model_preserves_output(self):
        output = self.root/'out.tsv'
        output.write_text('original')
        for text in ('locator\tbad\n', record(str(self.folder), 99),
                     record(str(self.folder), 0) + 'bad\trow\n',
                     record(str(self.folder), 0).replace('\t0\t', '\t-1\t')):
            with self.subTest(text=text):
                out, err = io.StringIO(), io.StringIO()
                with patch('sys.stdin', io.StringIO(text)), contextlib.redirect_stdout(out), contextlib.redirect_stderr(err), self.assertRaises(SystemExit):
                    SCRIPT.main(self.args + ['--output', str(output)])
                self.assertEqual(out.getvalue(), '')
                self.assertIn('stdin record', err.getvalue())
                self.assertEqual(output.read_text(), 'original')

    def test_shared_format_quotes_and_validation(self):
        locator = '/path/with\ttab and "quote"\nand newline'
        parsed = list(read_patch_inputs(io.StringIO(record(locator, 12))))
        self.assertEqual(parsed[0].locator, locator)
        self.assertEqual(parsed[0].model_id, 12)
        for text in ('a\t12345678\t0\tnan\tNA\n', '"unterminated',
                     'a\t12345678\t0\t1.1\tNA\n'):
            with self.assertRaises(ValueError):
                list(read_patch_inputs(io.StringIO(text)))

    def test_output_cannot_overwrite_piped_input(self):
        with patch('sys.stdin', io.StringIO(str(self.folder)+'\n')), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            SCRIPT.main(self.args + ['--output', str(self.folder)])
        self.assertTrue(self.folder.is_dir())

    def test_actual_list_results_pipe_selects_highest_only(self):
        # Invalid CIF and absent full-data on model 1 must never be read.
        (self.folder/'fold_summary_confidences_1.json').write_text('{"iptm":0.1,"ptm":0.1}')
        manifest = self.root/'manifest.tsv'
        with manifest.open('w', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t')
            writer.writerow(['fold_locator', 'fold_id', 'sequence_id'])
            writer.writerow([str(self.folder), '12345678', 'b911edee'])
        producer = subprocess.Popen([sys.executable, str(SCRIPTS/'alphafold-list-results.py'),
                                     str(manifest), 'b911edee', '--highest'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            consumer = subprocess.run([sys.executable, str(SCRIPTS/'alphafold-pdockq2-patches.py'),
                                       *self.args], stdin=producer.stdout, capture_output=True, text=True, timeout=30)
            producer.stdout.close()
            _, errors = producer.communicate(timeout=30)
            self.assertEqual(producer.returncode, 0, errors)
            self.assertEqual(consumer.returncode, 0, consumer.stderr)
            self.assertIn('Model 0 | ipTM 0.200000 | pTM 0.300000', consumer.stdout)
            self.assertNotIn('Model 1', consumer.stdout)
        finally:
            if producer.poll() is None:
                producer.kill()
                producer.communicate()
