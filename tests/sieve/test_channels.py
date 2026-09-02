import math
import os
import tempfile
import unittest

from sieve.channels import (
    CANONICAL_AMINO_ACIDS,
    EISENBERG_HYDROPHOBICITY,
    KYTE_DOOLITTLE_HYDROPATHY,
    TOP_IDP_DISORDER_PROPENSITY,
    ChannelBackground,
    CompositionBiasChannel,
    DipeptideFrequencyChannel,
    DisorderPropensityChannel,
    HydropathyChannel,
    HydrophobicMomentChannel,
    NetChargeChannel,
    SequenceEntropyChannel,
    ShortMotifChannel,
    iter_fasta_records,
)


class FastaFixture:

    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tempdir.cleanup()

    def write_fasta(self, text):
        path = os.path.join(self.tempdir.name, "proteome.faa")
        with open(path, "w", encoding="utf-8") as fasta:
            fasta.write(text)
        return path


class TestProteinChannels(unittest.TestCase):

    def test_scales_cover_exactly_the_canonical_amino_acids(self):
        for scale in (
            KYTE_DOOLITTLE_HYDROPATHY,
            EISENBERG_HYDROPHOBICITY,
            TOP_IDP_DISORDER_PROPENSITY,
        ):
            self.assertEqual(set(scale), CANONICAL_AMINO_ACIDS)

    def test_composition_uses_centered_windows_and_make_function_zscores(self):
        channel = CompositionBiasChannel("acidic", 1, "A")
        function = channel.make_function(ChannelBackground(0.5, 0.5))

        raw_values = channel.raw_values("aaac")
        values = function("aaac")

        for actual, expected in zip(raw_values, [1.0, 1.0, 2 / 3, 0.5]):
            self.assertAlmostEqual(actual, expected)
        for actual, expected in zip(values, [1.0, 1.0, 1 / 3, 0.0]):
            self.assertAlmostEqual(actual, expected)

    def test_net_charge_counts_only_k_r_d_e(self):
        channel = NetChargeChannel("charge", 1)

        self.assertEqual(channel.raw_values("KDA"), [0.0, 0.0, -0.5])

    def test_hydropathy_uses_kyte_doolittle_mean(self):
        channel = HydropathyChannel("hydropathy", 1)

        self.assertEqual(channel.raw_values("IV"), [4.35, 4.35])

    def test_hydrophobic_moment_is_normalized_vector_magnitude(self):
        scale = {amino_acid: 0.0 for amino_acid in CANONICAL_AMINO_ACIDS}
        scale["A"] = 1.0
        channel = HydrophobicMomentChannel(
            "moment", 1, angle=180.0, scale=scale
        )

        values = channel.raw_values("AA")

        self.assertAlmostEqual(values[0], 0.0, places=15)
        self.assertAlmostEqual(values[1], 0.0, places=15)

    def test_disorder_propensity_uses_top_idp_mean(self):
        channel = DisorderPropensityChannel("disorder", 1)

        values = channel.raw_values("PW")

        self.assertAlmostEqual(values[0], (0.987 - 0.884) / 2)
        self.assertAlmostEqual(values[1], (0.987 - 0.884) / 2)

    def test_entropy_is_unnormalized_shannon_entropy_in_bits(self):
        channel = SequenceEntropyChannel("entropy", 1)

        values = channel.raw_values("AAC")

        expected_middle = -(2 / 3) * math.log2(2 / 3) - (1 / 3) * math.log2(1 / 3)
        self.assertAlmostEqual(values[0], 0.0)
        self.assertAlmostEqual(values[1], expected_middle)
        self.assertAlmostEqual(values[2], 1.0)

    def test_dipeptide_frequency_counts_overlapping_pairs(self):
        channel = DipeptideFrequencyChannel("aa", 1, "AA")

        self.assertEqual(channel.raw_values("AACA"), [1.0, 0.5, 0.0, 0.0])

    def test_short_motif_marks_all_residues_in_overlapping_matches(self):
        channel = ShortMotifChannel("motif", "ACA")

        self.assertEqual(channel.raw_values("ACACA"), [1, 1, 1, 1, 1])
        self.assertIsNone(channel.compute_background("unused.faa"))

    def test_short_motif_leaves_residues_outside_matches_unmarked(self):
        channel = ShortMotifChannel("motif", "ACD")

        self.assertEqual(
            channel.make_function()("ACDQQACD"),
            [1, 1, 1, 0, 0, 1, 1, 1],
        )

    def test_empty_sequence_returns_empty_vector(self):
        self.assertEqual(HydropathyChannel("hydro", 2).raw_values(""), [])
        self.assertEqual(ShortMotifChannel("motif", "A").raw_values(""), [])

    def test_noncanonical_sequence_raises(self):
        channel = HydropathyChannel("hydro", 1)

        with self.assertRaisesRegex(ValueError, "noncanonical.*X"):
            channel.raw_values("AX")

    def test_invalid_background_raises(self):
        for sigma in (0.0, -1.0, float("nan"), float("inf")):
            with self.subTest(sigma=sigma):
                with self.assertRaisesRegex(ValueError, "bg_sigma"):
                    ChannelBackground(0.0, sigma)
        with self.assertRaisesRegex(TypeError, "requires a ChannelBackground"):
            NetChargeChannel("charge", 1).make_function()

    def test_invalid_short_name_and_radius_raise(self):
        with self.assertRaisesRegex(ValueError, "short_name"):
            NetChargeChannel("net-charge", 1)
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            SequenceEntropyChannel("entropy", -1)
        with self.assertRaisesRegex(TypeError, "integer"):
            SequenceEntropyChannel("entropy", 1.5)

    def test_dipeptide_requires_positive_radius_and_two_residue_sequence(self):
        with self.assertRaisesRegex(ValueError, "at least 1"):
            DipeptideFrequencyChannel("aa", 0, "AA")
        channel = DipeptideFrequencyChannel("aa", 1, "AA")
        with self.assertRaisesRegex(ValueError, "at least two"):
            channel.raw_values("A")

    def test_custom_scale_must_cover_all_canonical_amino_acids(self):
        with self.assertRaisesRegex(ValueError, "missing"):
            HydropathyChannel("hydro", 1, scale={"A": 1.0})

    def test_zero_length_motif_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "zero-length"):
            ShortMotifChannel("motif", "A*")

    def test_definition_contains_reusable_channel_arguments(self):
        channel = CompositionBiasChannel("acidic", radius=5, residues="ED")

        self.assertEqual(channel.definition(), {
            "short_name": "acidic",
            "type": "composition_bias",
            "arguments": {"radius": 5, "residues": "DE"},
        })


class TestProteinChannelBackgrounds(FastaFixture, unittest.TestCase):

    def test_background_weights_residues_not_proteins(self):
        path = self.write_fasta(">one\nA\n\n>two\nc\nc\n")

        background = CompositionBiasChannel("a", 0, "A").compute_background(path)

        self.assertAlmostEqual(background.bg_mu, 1 / 3)
        self.assertAlmostEqual(background.bg_sigma, math.sqrt(2) / 3)

    def test_each_normalized_channel_computes_a_background(self):
        path = self.write_fasta(
            ">proteome_part_1\nACDEFGHIKLMN\n"
            ">proteome_part_2\nPQRSTVWYACDE\n"
        )
        channels = [
            CompositionBiasChannel("composition", 1, "AG"),
            NetChargeChannel("charge", 1),
            HydropathyChannel("hydro", 1),
            HydrophobicMomentChannel("moment", 1),
            DisorderPropensityChannel("disorder", 1),
            SequenceEntropyChannel("entropy", 1),
            DipeptideFrequencyChannel("ac", 1, "AC"),
        ]

        for channel in channels:
            with self.subTest(channel=channel.short_name):
                background = channel.compute_background(path)
                self.assertTrue(math.isfinite(background.bg_mu))
                self.assertGreater(background.bg_sigma, 0.0)

    def test_computed_background_standardizes_its_source_values(self):
        path = self.write_fasta(">one\nA\n>two\nCC\n")
        channel = CompositionBiasChannel("a", 0, "A")
        function = channel.make_function(channel.compute_background(path))

        values = function("A") + function("CC")

        self.assertAlmostEqual(sum(values) / len(values), 0.0)
        self.assertAlmostEqual(
            math.sqrt(sum(value * value for value in values) / len(values)),
            1.0,
        )

    def test_dipeptide_background_skips_short_records(self):
        path = self.write_fasta(">short\nA\n>usable\nAACA\n")

        background = DipeptideFrequencyChannel("aa", 1, "AA").compute_background(path)

        self.assertAlmostEqual(background.bg_mu, 0.375)
        self.assertAlmostEqual(background.bg_sigma, math.sqrt(0.171875))

    def test_dipeptide_background_raises_if_no_record_is_usable(self):
        path = self.write_fasta(">one\nA\n>two\nC\n")

        with self.assertRaisesRegex(ValueError, "no usable positions"):
            DipeptideFrequencyChannel("aa", 1, "AA").compute_background(path)

    def test_zero_variance_background_raises(self):
        path = self.write_fasta(">one\nAAAA\n")

        with self.assertRaisesRegex(ValueError, "standard deviation"):
            CompositionBiasChannel("a", 1, "A").compute_background(path)

    def test_fasta_parser_uses_first_header_token_and_normalizes_case(self):
        path = self.write_fasta(">one description\na\nc\n")

        self.assertEqual(list(iter_fasta_records(path)), [("one", "AC")])

    def test_noncanonical_fasta_residue_raises(self):
        path = self.write_fasta(">one\nAAX\n")

        with self.assertRaisesRegex(ValueError, "noncanonical.*X"):
            HydropathyChannel("hydro", 1).compute_background(path)

    def test_fasta_parser_can_report_and_skip_noncanonical_records(self):
        path = self.write_fasta(">valid_1\nAC\n>invalid\nAX\n>valid_2\nDE\n")
        skipped = []

        records = list(iter_fasta_records(
            path,
            on_invalid_sequence=lambda sequence_id, error: skipped.append(
                (sequence_id, str(error))
            ),
        ))

        self.assertEqual(records, [("valid_1", "AC"), ("valid_2", "DE")])
        self.assertEqual(skipped, [(
            "invalid",
            "sequence contains noncanonical amino acids: X",
        )])

    def test_duplicate_fasta_sequence_id_raises(self):
        path = self.write_fasta(">one first\nA\n>one second\nC\n")

        with self.assertRaisesRegex(ValueError, "duplicate FASTA sequence ID"):
            list(iter_fasta_records(path))

    def test_malformed_fasta_raises(self):
        cases = {
            "no header": "ACD\n",
            "empty header": ">\nACD\n",
            "empty record": ">one\n>two\nACD\n",
            "no records": "\n\n",
        }
        for name, fasta_text in cases.items():
            with self.subTest(name=name):
                path = self.write_fasta(fasta_text)
                with self.assertRaises(ValueError):
                    HydropathyChannel("hydro", 1).compute_background(path)
