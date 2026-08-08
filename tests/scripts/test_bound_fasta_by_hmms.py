import argparse
import io
import os
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sieve.protein import ProteinHMMAlignment
from tests.scripts.helpers import load_script


class TestBoundFastaByHmms(unittest.TestCase):

    def setUp(self):
        repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.script = load_script(os.path.join(repo, "scripts", "bound-fasta-by-hmms.py"))

    def hit(self, accession, model_length, hmm_start, hmm_end, sequence_start, sequence_end):
        return self.script.DomainHit(
            accession, "model", model_length, hmm_start, hmm_end, sequence_start, sequence_end,
        )

    class Alignment:
        def __init__(self, mappings):
            self.mappings = mappings

        def aa_hmm_pos_1b(self, position):
            return self.mappings.get(position)

    def bound(self, sequences, n_hits, c_hits, n_mappings=None, c_mappings=None,
              n_position=1, c_position=None, acc_desc=None, report=None):
        n_mappings = {} if n_mappings is None else n_mappings
        c_mappings = {} if c_mappings is None else c_mappings
        n_alignments = {
            accession: self.Alignment(n_mappings.get(accession, {n_position: (2, "B")}))
            for accession in sequences
        }
        default_c_positions = {
            hit.sequence_accession: hit.model_length if c_position is None else c_position
            for hit in c_hits
        }
        c_alignments = {
            accession: self.Alignment(c_mappings.get(
                accession,
                {default_c_positions.get(accession, 1): (8, "H")},
            ))
            for accession in sequences
        }
        with patch.object(
            self.script,
            "hmm_align_sequences",
            side_effect=[n_alignments, c_alignments],
        ):
            return self.script.align_and_bound(
                sequences, n_hits, c_hits, "n.hmm", "c.hmm",
                n_position=n_position, c_position=c_position,
                acc_desc=acc_desc, report=report,
            )

    def test_bounds_at_default_model_endpoints(self):
        sequences = {"X": "ABCDEFGHIJ"}
        bounded = self.bound(
            sequences,
            [self.hit("X", 5, 1, 4, 2, 5)],
            [self.hit("X", 6, 2, 6, 4, 8)],
        )
        self.assertEqual(bounded, {"X_bounded_2-8": "BCDEFGH"})

    def test_uses_alternative_model_positions(self):
        bounded = self.bound(
            {"X": "ABCDEFGHIJ"},
            [self.hit("X", 5, 1, 5, 2, 7)],
            [self.hit("X", 200, 10, 199, 5, 9)],
            n_mappings={"X": {2: (3, "C")}},
            c_mappings={"X": {197: (9, "I")}},
            n_position=2,
            c_position=197,
        )
        self.assertEqual(bounded, {"X_bounded_3-9": "CDEFGHI"})

    def test_crops_at_exact_coordinates_across_an_alignment_insertion(self):
        alignment = MultipleSeqAlignment([SeqRecord(Seq("ABCDE"), id="X")])
        alignment.column_annotations["reference_annotation"] = "xx.xx"
        mapped = ProteinHMMAlignment(alignment, "X")
        bounded = self.script.bounded_sequences(
            {"X": "ABCDE"},
            {"X": (2, 3)},
            {"X": mapped},
            {"X": mapped},
        )
        self.assertEqual(bounded, {"X_bounded_2-4": "BCD"})

    def test_appends_accession_description(self):
        bounded = self.bound(
            {"X": "ABCDEFGHIJ"},
            [self.hit("X", 5, 1, 5, 2, 5)],
            [self.hit("X", 6, 1, 6, 4, 8)],
            acc_desc="example.v1",
        )
        self.assertEqual(bounded, {"X_bounded_2-8|example.v1": "BCDEFGH"})

    def test_rejects_whitespace_in_accession_description(self):
        with self.assertRaisesRegex(argparse.ArgumentTypeError, "cannot contain whitespace"):
            self.script.accession_description("two words")

    def test_requires_requested_n_and_c_model_boundaries(self):
        cases = [
            ([self.hit("X", 5, 2, 5, 2, 5)], [self.hit("X", 6, 1, 6, 4, 8)]),
            ([self.hit("X", 5, 1, 5, 2, 5)], [self.hit("X", 6, 1, 5, 4, 8)]),
        ]
        for n_hits, c_hits in cases:
            messages = []
            bounded = self.bound(
                {"X": "ABCDEFGHIJ"}, n_hits, c_hits, report=messages.append,
            )
            self.assertEqual(bounded, {})
            self.assertEqual(len(messages), 1)

    def test_rejects_missing_multiple_reversed_and_out_of_range_hits(self):
        valid_n = self.hit("X", 5, 1, 5, 2, 5)
        valid_c = self.hit("X", 6, 1, 6, 4, 8)
        cases = [
            ([], [valid_c], "found 0"),
            ([valid_n, valid_n], [valid_c], "found 2"),
            ([valid_n], [], "found 0"),
            ([valid_n], [valid_c, valid_c], "found 2"),
        ]
        for n_hits, c_hits, expected in cases:
            messages = []
            bounded = self.bound(
                {"X": "ABCDEFGHIJ"}, n_hits, c_hits, report=messages.append,
            )
            self.assertEqual(bounded, {})
            self.assertIn(expected, messages[0])

    def test_rejects_invalid_positions_and_unknown_hit_accessions(self):
        with self.assertRaisesRegex(ValueError, "at least 1"):
            self.script.qualifying_hits({}, [], [], n_position=0)
        with self.assertRaisesRegex(ValueError, "unknown FASTA accession"):
            self.script.qualifying_hits({}, [self.hit("other", 5, 1, 5, 1, 5)], [])

    def test_rejects_deletions_and_invalid_mapped_boundaries(self):
        valid_n = [self.hit("X", 5, 1, 5, 2, 5)]
        valid_c = [self.hit("X", 6, 1, 6, 4, 8)]
        cases = [
            ({"X": {}}, {"X": {6: (8, "H")}}, "aligns to a deletion"),
            ({"X": {1: (2, "B")}}, {"X": {}}, "aligns to a deletion"),
            ({"X": {1: (9, "I")}}, {"X": {6: (8, "H")}}, "is after"),
            ({"X": {1: (0, "A")}}, {"X": {6: (8, "H")}}, "outside sequence"),
        ]
        for n_mappings, c_mappings, expected in cases:
            messages = []
            bounded = self.bound(
                {"X": "ABCDEFGHIJ"}, valid_n, valid_c,
                n_mappings=n_mappings, c_mappings=c_mappings, report=messages.append,
            )
            self.assertEqual(bounded, {})
            self.assertIn(expected, messages[0])

    def test_rejects_duplicate_fasta_accessions(self):
        with tempfile.NamedTemporaryFile("w", suffix=".faa") as f:
            f.write(">X\nAAA\n>X\nBBB\n")
            f.flush()
            with self.assertRaisesRegex(ValueError, "Duplicate FASTA accession: X"):
                self.script.read_unique_fasta(f.name)

    def test_search_uses_cut_ga_and_parses_model_length(self):
        def fake_run(hmm, fasta, domtblout, use_cut_ga):
            self.assertEqual((hmm, fasta, use_cut_ga), ("n.hmm", "input.faa", True))
            with open(domtblout, "w", encoding="utf-8") as f:
                f.write(
                    "X - 100 nmodel NMODEL 42 1e-10 50 0 1 1 1e-10 1e-10 50 0 "
                    "1 42 3 44 3 44 0.99 description\n"
                )
            return subprocess.CompletedProcess([], 0)

        with patch.object(self.script, "run_hmmsearch", side_effect=fake_run):
            hits = self.script.search("n.hmm", "input.faa")
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].model_name, "nmodel")
        self.assertEqual(
            (hits[0].model_length, hits[0].hmm_start, hits[0].hmm_end,
             hits[0].sequence_start, hits[0].sequence_end),
            (42, 1, 42, 3, 44),
        )

    def test_search_propagates_hmmsearch_failure(self):
        failure = subprocess.CalledProcessError(1, ["hmmsearch"])
        with patch.object(self.script, "run_hmmsearch", side_effect=failure):
            with self.assertRaises(subprocess.CalledProcessError):
                self.script.search("n.hmm", "input.faa")

    def test_aligns_qualifying_sequences_in_two_batches_and_propagates_failure(self):
        sequences = {"X": "ABCDEFGHIJ", "Y": "ABCDEFGHIJ"}
        n_hits = [self.hit("X", 5, 1, 5, 2, 5)]
        c_hits = [self.hit("X", 6, 1, 6, 4, 8)]
        n_alignment = self.Alignment({1: (2, "B")})
        c_alignment = self.Alignment({6: (8, "H")})
        with patch.object(
            self.script,
            "hmm_align_sequences",
            side_effect=[{"X": n_alignment}, {"X": c_alignment}],
        ) as align:
            bounded = self.script.align_and_bound(
                sequences, n_hits, c_hits, "n.hmm", "c.hmm", report=lambda message: None,
            )
        self.assertEqual(bounded, {"X_bounded_2-8": "BCDEFGH"})
        self.assertEqual(align.call_args_list[0].args, ("n.hmm", {"X": "ABCDEFGHIJ"}))
        self.assertEqual(align.call_args_list[1].args, ("c.hmm", {"X": "ABCDEFGHIJ"}))

        with patch.object(self.script, "hmm_align_sequences", side_effect=RuntimeError("bad alignment")):
            with self.assertRaisesRegex(RuntimeError, "bad alignment"):
                self.script.align_and_bound(sequences, n_hits, c_hits, "n.hmm", "c.hmm")

    def test_main_runs_each_search_once_and_writes_stdout(self):
        with tempfile.NamedTemporaryFile("w", suffix=".faa") as f:
            f.write(">X\nABCDEFGHIJ\n")
            f.flush()
            n_hit = self.hit("X", 5, 1, 5, 2, 5)
            c_hit = self.hit("X", 6, 1, 6, 4, 8)
            n_alignment = self.Alignment({1: (2, "B")})
            c_alignment = self.Alignment({6: (8, "H")})
            stdout = io.StringIO()
            with patch.object(self.script, "search", side_effect=[[n_hit], [c_hit]]) as search:
                with patch.object(
                    self.script,
                    "hmm_align_sequences",
                    side_effect=[{"X": n_alignment}, {"X": c_alignment}],
                ):
                    with patch("sys.stdout", stdout):
                        self.script.main([
                            "--n-hmm", "n.hmm", "--c-hmm", "c.hmm", "--fasta", f.name,
                        ])
            self.assertEqual(search.call_count, 2)
            self.assertEqual(stdout.getvalue(), ">X_bounded_2-8\nBCDEFGH\n")


if __name__ == "__main__":
    unittest.main()
