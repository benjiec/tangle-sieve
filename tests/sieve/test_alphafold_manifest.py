import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from sieve.alphafold_files import discover, model_number
from sieve.alphafold_manifest import chain_sequences, fold_id, read_candidates, run


def cif(sequences=('AC', 'AC', 'DE')):
    names = dict(A='ALA', C='CYS', D='ASP', E='GLU')
    text = 'data_test\nloop_\n_entity.id\n_entity.type\n'
    text += ''.join(f'{i} polymer\n' for i in range(1, len(sequences)+1))
    text += '#\nloop_\n_struct_asym.id\n_struct_asym.entity_id\n'
    text += ''.join(f'{chr(64+i)} {i}\n' for i in range(1, len(sequences)+1))
    text += '#\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n'
    text += ''.join(f'{i} {n} {names[a]}\n' for i, seq in enumerate(sequences, 1) for n, a in enumerate(seq, 1))
    return text+'#\n'


def archive_bytes(members):
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, 'w') as archive:
        for name, text in members.items():
            archive.writestr(name, text)
    return stream.getvalue()


class ManifestTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.inputs = self.root / 'folds'
        self.inputs.mkdir()
        self.fasta = self.root / 'candidates.faa'
        self.fasta.write_text('>first full header\nAC\n>second full header\nDE\n')
        self.output = self.root / 'manifest.tsv'
        self.log = io.StringIO()

    def fold(self, name='one', sequences=('AC', 'AC', 'DE')):
        folder = self.inputs / name
        folder.mkdir(parents=True, exist_ok=True)
        for i in range(5):
            (folder / f'fold_model_{i}.cif').write_text(cif(sequences))
        return folder

    def run_manifest(self):
        return run(self.inputs, self.fasta, self.output, self.log)

    def rows(self):
        with self.output.open() as stream:
            return list(csv.DictReader(stream, delimiter='\t'))

    def test_hash(self):
        seqs = ['DE', 'AC', 'AC']
        expected = hashlib.sha1(b'["AC","AC","DE"]').hexdigest()[-8:]
        self.assertEqual(fold_id(seqs), expected)
        self.assertEqual(fold_id(seqs), fold_id(reversed(seqs)))
        self.assertNotEqual(fold_id(['A', 'CD']), fold_id(['AC', 'D']))
        self.assertNotEqual(fold_id(['AC']), fold_id(['AC', 'AC']))

    def test_fasta_errors_collect_and_prevent_scan(self):
        self.fasta.write_text('stray\n>one\nac\n>two\nAC\n>\n>bad\nA*\n')
        with patch('sieve.alphafold_manifest.discover', side_effect=AssertionError):
            self.assertEqual(self.run_manifest(), 1)
        for message in ['duplicate sequence', 'empty FASTA header', 'empty sequence', 'invalid sequence', 'before first']:
            self.assertIn(message, self.log.getvalue())
        self.assertFalse(self.output.exists())

    def test_empty_fasta(self):
        self.fasta.write_text('')
        self.assertEqual(self.run_manifest(), 1)
        self.assertIn('no FASTA records', self.log.getvalue())

    def test_full_header_and_normalization(self):
        self.fasta.write_text('>  id description\twith spaces  \n a c\n')
        errors = []
        self.assertEqual(read_candidates(self.fasta, errors.append), {'AC': '  id description\twith spaces  '})
        self.assertEqual(errors, [])

    def test_copies_and_cache_preserved_without_parse(self):
        self.fold()
        self.assertEqual(self.run_manifest(), 0)
        rows = self.rows()
        self.assertEqual([r['copies'] for r in rows], ['2', '1'])
        self.assertEqual(rows[0]['chain_ids'], 'A,B')
        self.assertEqual(rows[0]['description'], 'first full header')
        self.fasta.write_text('>new full header\nAC\n>second\nDE\n>added\nAA\n')
        with patch('sieve.alphafold_manifest.identify', side_effect=AssertionError):
            self.assertEqual(self.run_manifest(), 0)
            self.assertEqual(self.run_manifest(), 0)
        self.assertEqual(self.rows(), rows)

    def test_legacy_cache_migration(self):
        self.fold()
        self.assertEqual(self.run_manifest(), 0)
        rows = self.rows()
        self.assertEqual(set(rows[0]), {'source_path', 'fold_path', 'fold_id',
                                      'description', 'chain_ids', 'copies'})
        for extra in ({'sequence': 'AC', 'source_fingerprint': 'old'},
                      {'sequence_sha256': 'old', 'source_fingerprint': 'old'}):
            with self.subTest(columns=extra):
                legacy = [dict(row, **extra) for row in rows]
                with self.output.open('w', newline='') as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(legacy[0]), delimiter='\t')
                    writer.writeheader()
                    writer.writerows(legacy)
                with patch('sieve.alphafold_manifest.identify', side_effect=AssertionError):
                    self.assertEqual(self.run_manifest(), 0)
                self.assertEqual(self.rows(), rows)

    def test_fail_continue_then_retry(self):
        self.fold('bad', ('AA',))
        self.fold('good')
        self.assertEqual(self.run_manifest(), 1)
        self.assertIn('ERROR: AA not found in ', self.log.getvalue())
        self.assertEqual(len(self.rows()), 2)
        self.fasta.write_text(self.fasta.read_text()+'>missing\nAA\n')
        self.assertEqual(self.run_manifest(), 0)
        self.assertEqual(len(self.rows()), 3)
        self.assertIn('cached: 1', self.log.getvalue())

    def test_missing_sequences_each_start_error_line(self):
        folder = self.fold('missing', ('AA', 'CC', 'AA'))
        self.assertEqual(self.run_manifest(), 1)
        lines = [line for line in self.log.getvalue().splitlines() if line.startswith('ERROR:')]
        self.assertEqual(lines, [
            f'ERROR: {sequence} not found in {folder.resolve()}::. chain {chain}'
            for chain, sequence in [('A', 'AA'), ('B', 'CC'), ('C', 'AA')]
        ])
        self.assertFalse(self.output.exists())

    def test_changed_fold_stays_cached_until_rows_deleted(self):
        folder = self.fold()
        self.assertEqual(self.run_manifest(), 0)
        original = self.output.read_bytes()
        (folder / 'fold_model_1.cif').write_text(cif(('AA',)))
        with patch.object(Path, 'read_text', side_effect=AssertionError('cached model read')):
            self.assertEqual(self.run_manifest(), 0)
        self.assertEqual(self.output.read_bytes(), original)
        self.output.write_text(original.decode().splitlines()[0] + '\n')
        self.assertEqual(self.run_manifest(), 1)
        self.assertIn('differ across models', self.log.getvalue())

    def test_removed_fasta_sequence_leaves_cache_unchanged(self):
        self.fold()
        self.assertEqual(self.run_manifest(), 0)
        original = self.output.read_bytes()
        self.fasta.write_text('>first\nAC\n')
        self.assertEqual(self.run_manifest(), 0)
        self.assertEqual(self.output.read_bytes(), original)

    def test_invalid_fasta_leaves_legacy_cache_untouched(self):
        self.fold()
        self.assertEqual(self.run_manifest(), 0)
        rows = [dict(row, source_fingerprint='old') for row in self.rows()]
        with self.output.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        original = self.output.read_bytes()
        self.fasta.write_text('>first\nAC\n>duplicate\nAC\n')
        self.assertEqual(self.run_manifest(), 1)
        self.assertEqual(self.output.read_bytes(), original)

    def test_nested_zip_and_duplicates_and_templates(self):
        self.fold()
        inner = archive_bytes({'deep/fold_model_0.cif': cif(), 'deep/templates/hit.cif': 'bad'})
        (self.inputs / 'outer.zip').write_bytes(archive_bytes({'inside.zip': inner, 'other/fold_model_0.cif': cif()}))
        (self.inputs / '._junk.zip').write_text('ignored')
        self.assertEqual(self.run_manifest(), 0)
        rows = self.rows()
        self.assertEqual(len(rows), 6)
        self.assertEqual(len({r['fold_id'] for r in rows}), 1)
        self.assertEqual({r['fold_path'] for r in rows}, {'.', 'inside.zip!/deep', 'other'})

    def test_bad_archives_continue(self):
        self.fold()
        (self.inputs / 'bad.zip').write_text('bad')
        (self.inputs / 'nested.zip').write_bytes(archive_bytes({'bad.zip': b'bad', 'ok/fold_model_0.cif': cif()}))
        self.assertEqual(self.run_manifest(), 1)
        self.assertEqual(len(self.rows()), 4)

    def test_duplicate_model_number(self):
        folder = self.fold()
        (folder / 'other_model_0.cif').write_text(cif())
        self.assertEqual(self.run_manifest(), 1)
        self.assertIn('duplicate model numbers', self.log.getvalue())

    def test_unknown_residue_and_nonpolymer_rejected(self):
        for text in [cif().replace('ALA', 'UNK'), cif().replace('1 polymer', '1 non-polymer')]:
            with self.assertRaises(ValueError):
                chain_sequences(text)

    def test_invalid_cache_and_input_overwrite(self):
        self.output.write_text('wrong\n')
        with self.assertRaises(ValueError):
            self.run_manifest()
        with self.assertRaises(ValueError):
            run(self.inputs, self.fasta, self.fasta)

    def test_no_folds(self):
        self.assertEqual(self.run_manifest(), 1)
        self.assertIn('no folds found', self.log.getvalue())

    def test_unrecognized_fold_reports_error(self):
        self.fold('good')
        bad = self.inputs / 'bad'
        bad.mkdir()
        (bad / 'fold_job_request.json').write_text('{}')
        self.assertEqual(self.run_manifest(), 1)
        self.assertIn('no recognized model CIF', self.log.getvalue())
        self.assertEqual(len(self.rows()), 2)

    def test_checkpoint_survives_later_interruption(self):
        from sieve.alphafold_manifest import identify
        self.fold('a')
        self.fold('b')
        calls = []
        def interrupted(fold, candidates):
            calls.append(fold)
            if len(calls) == 2:
                raise KeyboardInterrupt()
            return identify(fold, candidates)
        with patch('sieve.alphafold_manifest.identify', side_effect=interrupted):
            with self.assertRaises(KeyboardInterrupt):
                self.run_manifest()
        self.assertEqual(len(self.rows()), 2)
        self.assertEqual(self.run_manifest(), 0)
        self.assertEqual(len(self.rows()), 4)

    def test_failed_atomic_replace_keeps_previous_manifest(self):
        self.fold('a')
        self.assertEqual(self.run_manifest(), 0)
        original = self.output.read_bytes()
        self.fold('b')
        with patch('sieve.alphafold_cache.os.replace', side_effect=OSError('disk failure')):
            with self.assertRaises(OSError):
                self.run_manifest()
        self.assertEqual(self.output.read_bytes(), original)

    def test_model_filter(self):
        self.assertEqual(model_number('dir/fold_model_4.cif'), 4)
        for name in ['templates/hit.cif', '__MACOSX/fold_model_0.cif', '._fold_model_0.cif']:
            self.assertIsNone(model_number(name))


if __name__ == '__main__':
    unittest.main()
