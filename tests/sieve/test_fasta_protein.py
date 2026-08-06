import unittest

from sieve.fasta_protein import FastaProtein


class TestFastaProtein(unittest.TestCase):

    def test_leader_sequence_candidates_include_full_sequence_and_early_methionines(self):
        protein = FastaProtein("p1", "AAMQQMQQQQQM")

        candidates = protein.leader_sequence_candidates()

        self.assertEqual(
            [(candidate.accession, candidate.start_label, candidate.start_aa_1b, candidate.sequence) for candidate in candidates],
            [
                ("p1", "", 1, "AAMQQMQQQQQM"),
                ("p1_with_leader_3_M", "3", 3, "MQQMQQQQQM"),
                ("p1_with_leader_6_M", "6", 6, "MQQQQQM"),
            ],
        )
        self.assertEqual([candidate.protein_start_aa_1b for candidate in candidates], [1, 3, 6])

    def test_leader_sequence_candidates_always_include_full_sequence_without_methionines(self):
        protein = FastaProtein("p1", "AAAA")

        candidates = protein.leader_sequence_candidates()

        self.assertEqual(
            [(candidate.accession, candidate.start_label, candidate.start_aa_1b, candidate.sequence) for candidate in candidates],
            [("p1", "", 1, "AAAA")],
        )

    def test_anchor_leader_candidates_stop_at_sequence_start_and_include_full_sequence(self):
        protein = FastaProtein("p1", "MAAMAAAAM")

        candidates = protein.leader_sequence_candidates_at_anchor(
            anchor_aa_1b=5,
            anchor_label="PF1",
            window_start=-5,
            window_end=2,
        )

        self.assertEqual(
            [(candidate.accession, candidate.start_label, candidate.start_aa_1b, candidate.sequence) for candidate in candidates],
            [
                ("p1", "", 1, "MAAMAAAAM"),
                ("p1_with_leader_u4_PF1_M", "u4_PF1", -4, "MAAMAAAAM"),
                ("p1_with_leader_u1_PF1_M", "u1_PF1", -1, "MAAAAM"),
            ],
        )
        self.assertEqual([candidate.protein_start_aa_1b for candidate in candidates], [1, 1, 4])

    def test_anchor_leader_candidates_include_downstream_anchor_window(self):
        protein = FastaProtein("p1", "AAAAMAA")

        candidates = protein.leader_sequence_candidates_at_anchor(
            anchor_aa_1b=4,
            anchor_label="PF1",
            window_start=-1,
            window_end=2,
        )

        self.assertEqual(
            [(candidate.start_label, candidate.start_aa_1b, candidate.sequence) for candidate in candidates],
            [
                ("", 1, "AAAAMAA"),
                ("2_PF1", 2, "MAA"),
            ],
        )
        self.assertEqual([candidate.protein_start_aa_1b for candidate in candidates], [1, 5])

    def test_genomic_locus_is_not_available(self):
        protein = FastaProtein("p1", "MAA")

        with self.assertRaisesRegex(ValueError, "FASTA-only protein p1"):
            protein.genomic_locus()

        with self.assertRaisesRegex(ValueError, "FASTA-only protein p1"):
            protein.genomic_locus_with_leader()
