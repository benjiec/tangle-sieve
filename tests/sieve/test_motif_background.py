import os
import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio.SeqIO.FastaIO import SimpleFastaParser
from sieve.motif_background import genomic_background, valid_start_intervals


class TestMotifBackground(unittest.TestCase):
    def test_valid_starts_match_exhaustive_reference(self):
        rng = random.Random(12)
        for length in (1, 5, 10, 20):
            for sequence in ('', 'N' * 100, 'A' * 100, 'N' * 50 + 'ACGT' * 20,
                             ''.join(rng.choices('ACGTNn', k=100))):
                with self.subTest(length=length, sequence=sequence):
                    actual = [i for start, end in valid_start_intervals(sequence, length) for i in range(start, end)]
                    expected = [i for i in range(len(sequence) - length + 1)
                                if sequence[i:i + length].upper().count('N') <= length // 10]
                    self.assertEqual(actual, expected)

    def test_sampling_coordinates_quality_cache_and_invalidation(self):
        with tempfile.TemporaryDirectory() as directory, patch.dict(os.environ, {'XDG_CACHE_HOME': ''}):
            with patch.dict(os.environ, {'XDG_CACHE_HOME': directory}):
                source = Path(directory) / 'genome.fna'
                sequences = {'gapped': 'N' * 300 + 'ACGT' * 100, 'tiny': 'A', 'other': 'TGCA' * 100}
                source.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in sequences.items()))
                path = genomic_background(source, 20, nseq=100)
                with open(path) as fasta:
                    rows = list(SimpleFastaParser(fasta))
                self.assertEqual(len(rows), 100)
                seen = set()
                for title, sequence in rows:
                    contig, bounds = title.split()[1].split(':')
                    start, end = map(int, bounds.split('-'))
                    seen.add(contig)
                    self.assertEqual(sequence, sequences[contig][start:end])
                    self.assertEqual(len(sequence), 20)
                    self.assertLessEqual(sequence.count('N'), 2)
                self.assertEqual(seen, {'gapped', 'other'})
                original = Path(path).read_text()
                with patch('sieve.motif_background.valid_start_intervals', side_effect=AssertionError('cache miss')):
                    self.assertEqual(genomic_background(source, 20, nseq=100), path)
                self.assertNotEqual(genomic_background(source, 21, nseq=100), path)
                source.write_text(source.read_text() + '>new\n' + 'A' * 30 + '\n')
                self.assertNotEqual(genomic_background(source, 20, nseq=100), path)
                self.assertEqual(Path(path).read_text(), original)

    def test_no_valid_windows_and_invalid_parameters(self):
        with tempfile.TemporaryDirectory() as directory, patch.dict(os.environ, {'XDG_CACHE_HOME': directory}):
            source = Path(directory) / 'genome.fna'
            for sequence in ('NNNNNNNNNNNN', 'AAA', ''):
                source.write_text('>ctg\n' + sequence + '\n')
                with self.assertRaisesRegex(ValueError, 'No 10-base'):
                    genomic_background(source, 10)
            for length, count in [(0, 10), (True, 10), (10, 0), (10, 1.5)]:
                with self.assertRaises(ValueError):
                    genomic_background(source, length, count)
