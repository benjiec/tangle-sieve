import csv
import io
import os
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from tests.scripts.helpers import load_script


class TestAlignTranscriptsToLoci(unittest.TestCase):
    def setUp(self):
        repo = Path(__file__).resolve().parents[2]
        self.script = load_script(str(repo / "scripts/align-transcripts-to-loci.py"))
        self.aligner = self.script.make_aligner()

    def align(self, locus, transcript):
        return self.script.align_pair(self.aligner, locus, transcript, 2)

    def test_perfect_subsequence_ignores_locus_flanks(self):
        row = self.align("C" * 200 + "A" * 100 + "C" * 200, "A" * 100)
        self.assertEqual(row["score"], 100)
        self.assertEqual(row["transcript_coverage"], 1)
        self.assertEqual((row["locus_start"], row["locus_end"]), (200, 300))

    def test_partial_match(self):
        row = self.align("A" * 60, "A" * 100)
        self.assertEqual(row["score"], 60)
        self.assertEqual(row["transcript_coverage"], 0.6)

    def test_one_long_gap(self):
        row = self.align("A" * 50 + "C" * 40 + "A" * 50, "A" * 100)
        self.assertEqual(row["score"], 78)
        self.assertEqual((row["gap_count"], row["gap_bases"]), (1, 40))

    def test_ten_regions_and_local_fallback(self):
        # Repeated A blocks allow a better mismatch-based alignment for gap=1.
        for gap, score in [(1, 78), (5, 59.5), (10, 37), (20, 10)]:
            with self.subTest(gap=gap):
                row = self.align(("C" * gap).join(["A" * 10] * 10), "A" * 100)
                self.assertAlmostEqual(row["score"], score)
                self.assertEqual(row["gap_count"], 9 if gap in (5, 10) else 0)

    def test_reverse_complement_coordinates_in_original_transcript(self):
        row = self.align("ACGGA", "TCCGTTTT")
        self.assertEqual(row["strand"], "-")
        self.assertEqual((row["transcript_start"], row["transcript_end"]), (0, 5))
        self.assertEqual(row["raw_score"], 10)

    def test_no_positive_match(self):
        row = self.align("CCCC", "AAAA")
        self.assertEqual(row["score"], 0)
        self.assertEqual(row["strand"], ".")
        self.assertEqual(row["locus_start"], "")

    def test_ambiguity_is_not_rewarded(self):
        self.assertEqual(self.align("NNNN", "NNNN")["score"], 0)

    def test_mismatch_and_query_insertion(self):
        row = self.align("A" * 10 + "C" + "A" * 10, "A" * 10 + "G" + "A" * 10)
        self.assertEqual(row["raw_score"], 37)
        locus = "ACGATCGAGCTAGCATGCGA"
        row = self.align(locus, locus[:10] + "C" + locus[10:])
        self.assertEqual(row["raw_score"], 35)
        self.assertEqual(row["gap_bases"], 1)
        self.assertAlmostEqual(row["transcript_coverage"], 20 / 21)

    def test_forward_wins_orientation_tie(self):
        self.assertEqual(self.align("ACGT", "ACGT")["strand"], "+")

    def test_custom_scoring(self):
        aligner = self.script.make_aligner(4, -6, -10, -2)
        row = self.script.align_pair(aligner, "A" * 50 + "C" * 40 + "A" * 50, "A" * 100, 4)
        self.assertEqual(row["score"], 78)
        self.assertEqual(row["raw_score"], 312)

    def test_invalid_scoring(self):
        for values in [(0, -3, -5, -1), (2, 1, -5, -1), (2, -3, 1, -1),
                       (2, -3, -5, 1), (float("nan"), -3, -5, -1),
                       (2, -3, -5, float("-inf"))]:
            with self.subTest(values=values), self.assertRaises(ValueError):
                self.script.make_aligner(*values)

    def test_fasta_validation_and_normalization(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.fa"
            path.write_text(">rna description\nacgu\n")
            self.assertEqual(self.script.read_fasta(str(path)), [("rna", "ACGT")])
            for text in ["", ">empty\n", ">x\nACGT\n>x\nACGT\n", ">x\nAC-G\n"]:
                path.write_text(text)
                with self.subTest(text=text), self.assertRaises(ValueError):
                    self.script.read_fasta(str(path))

    def test_cli_all_pairs_order_ties_and_custom_parameters(self):
        with tempfile.TemporaryDirectory() as directory:
            loci, transcripts, output = [os.path.join(directory, name) for name in ("l.fa", "t.fa", "out.tsv")]
            Path(loci).write_text(">partial\nAAAA\n>exact\nAAAAAAAA\n>tied\nAAAAAAAA\n")
            Path(transcripts).write_text(">t1\nAAAAAAAA\n>t2\nCCCC\n")
            stdout = io.StringIO()
            with redirect_stdout(stdout):
                self.assertEqual(self.script.main([loci, transcripts, output, "--match", "3"]), 0)
            with open(output) as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(len(rows), 6)
            self.assertEqual([r["locus_id"] for r in rows[:3]], ["exact", "tied", "partial"])
            self.assertEqual([r["rank"] for r in rows[:3]], ["1", "1", "3"])
            self.assertEqual(float(rows[0]["score"]), 100)
            self.assertEqual(float(rows[0]["raw_score"]), 24)
            self.assertTrue(all(float(r["score"]) == 0 for r in rows[3:]))
            self.assertEqual(stdout.getvalue(),
                             "t1: exact, 100.0\n"
                             "t1: tied, 100.0\n"
                             "t2: partial, 0.0\n"
                             "t2: exact, 0.0\n"
                             "t2: tied, 0.0\n")
            Path(loci).write_text(">partial\nAAAA\n>exact\nAAAAAAAA\n")
            Path(transcripts).write_text(">t1\nAAAAAAAA\n")
            stdout = io.StringIO()
            with redirect_stdout(stdout):
                self.assertEqual(self.script.main([loci, transcripts, output]), 0)
            self.assertEqual(stdout.getvalue(),
                             "t1: exact, 100.0\n")
            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                self.script.main([loci, transcripts, output, "--match", "0"])


if __name__ == "__main__":
    unittest.main()
