import math
import os
import tempfile
import unittest

from sieve.channels import (
    CANONICAL_AMINO_ACIDS,
    EISENBERG_HYDROPHOBICITY,
    KYTE_DOOLITTLE_HYDROPATHY,
    TOP_IDP_DISORDER_PROPENSITY,
    compute_composition_background,
    compute_dipeptide_background,
    compute_disorder_background,
    compute_entropy_background,
    compute_hydrophobic_moment_background,
    compute_hydropathy_background,
    compute_net_charge_background,
    mk_composition_bias,
    mk_dipeptide_frequency,
    mk_disorder_propensity,
    mk_hydrophobic_moment,
    mk_hydropathy,
    mk_net_charge,
    mk_sequence_entropy,
    mk_short_motif,
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

    def test_composition_bias_uses_centered_truncated_windows_and_zscores(self):
        channel = mk_composition_bias(1, "A", bg_mu=0.5, bg_sigma=0.5)

        values = channel("aaac")

        self.assertEqual(len(values), 4)
        for actual, expected in zip(values, [1.0, 1.0, 1.0 / 3.0, 0.0]):
            self.assertAlmostEqual(actual, expected)

    def test_net_charge_counts_only_k_r_d_e(self):
        channel = mk_net_charge(1, bg_mu=0.0, bg_sigma=1.0)

        values = channel("KDA")

        self.assertEqual(values, [0.0, 0.0, -0.5])

    def test_hydropathy_uses_kyte_doolittle_mean(self):
        channel = mk_hydropathy(1, bg_mu=0.0, bg_sigma=1.0)

        values = channel("IV")

        self.assertEqual(values, [4.35, 4.35])

    def test_hydrophobic_moment_is_normalized_vector_magnitude(self):
        scale = {amino_acid: 0.0 for amino_acid in CANONICAL_AMINO_ACIDS}
        scale["A"] = 1.0
        channel = mk_hydrophobic_moment(
            1,
            bg_mu=0.0,
            bg_sigma=1.0,
            angle=180.0,
            scale=scale,
        )

        values = channel("AA")

        self.assertAlmostEqual(values[0], 0.0, places=15)
        self.assertAlmostEqual(values[1], 0.0, places=15)

    def test_disorder_propensity_uses_top_idp_mean(self):
        channel = mk_disorder_propensity(1, bg_mu=0.0, bg_sigma=1.0)

        values = channel("PW")

        self.assertAlmostEqual(values[0], (0.987 - 0.884) / 2)
        self.assertAlmostEqual(values[1], (0.987 - 0.884) / 2)

    def test_entropy_is_unnormalized_shannon_entropy_in_bits(self):
        channel = mk_sequence_entropy(1, bg_mu=0.0, bg_sigma=1.0)

        values = channel("AAC")

        self.assertAlmostEqual(values[0], 0.0)
        self.assertAlmostEqual(values[1], -(2 / 3) * math.log2(2 / 3) - (1 / 3) * math.log2(1 / 3))
        self.assertAlmostEqual(values[2], 1.0)

    def test_dipeptide_frequency_counts_overlapping_pairs(self):
        channel = mk_dipeptide_frequency(1, "AA", bg_mu=0.0, bg_sigma=1.0)

        values = channel("AACA")

        self.assertEqual(values, [1.0, 0.5, 0.0, 0.0])

    def test_short_motif_marks_all_residues_in_overlapping_matches(self):
        channel = mk_short_motif("ACA")

        values = channel("ACACA")

        self.assertEqual(values, [1, 1, 1, 1, 1])

    def test_short_motif_leaves_residues_outside_matches_unmarked(self):
        channel = mk_short_motif("ACD")

        values = channel("ACDQQACD")

        self.assertEqual(values, [1, 1, 1, 0, 0, 1, 1, 1])

    def test_empty_sequence_returns_empty_vector(self):
        self.assertEqual(
            mk_hydropathy(2, bg_mu=0.0, bg_sigma=1.0)(""),
            [],
        )
        self.assertEqual(mk_short_motif("A")(""), [])

    def test_noncanonical_sequence_raises(self):
        channel = mk_hydropathy(1, bg_mu=0.0, bg_sigma=1.0)

        with self.assertRaisesRegex(ValueError, "noncanonical.*X"):
            channel("AX")

    def test_invalid_background_raises_when_generator_is_built(self):
        for sigma in (0.0, -1.0, float("nan"), float("inf")):
            with self.subTest(sigma=sigma):
                with self.assertRaisesRegex(ValueError, "bg_sigma"):
                    mk_net_charge(1, bg_mu=0.0, bg_sigma=sigma)

    def test_invalid_radius_raises(self):
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            mk_sequence_entropy(-1, bg_mu=0.0, bg_sigma=1.0)
        with self.assertRaisesRegex(TypeError, "integer"):
            mk_sequence_entropy(1.5, bg_mu=0.0, bg_sigma=1.0)

    def test_dipeptide_requires_positive_radius_and_two_residue_sequence(self):
        with self.assertRaisesRegex(ValueError, "at least 1"):
            mk_dipeptide_frequency(0, "AA", bg_mu=0.0, bg_sigma=1.0)
        channel = mk_dipeptide_frequency(1, "AA", bg_mu=0.0, bg_sigma=1.0)
        with self.assertRaisesRegex(ValueError, "at least two"):
            channel("A")

    def test_custom_scale_must_cover_all_canonical_amino_acids(self):
        with self.assertRaisesRegex(ValueError, "missing"):
            mk_hydropathy(1, bg_mu=0.0, bg_sigma=1.0, scale={"A": 1.0})

    def test_zero_length_motif_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "zero-length"):
            mk_short_motif("A*")


class TestProteinChannelBackgrounds(FastaFixture, unittest.TestCase):

    def test_composition_background_weights_residues_not_proteins(self):
        path = self.write_fasta(">one\nA\n\n>two\nc\nc\n")

        bg_mu, bg_sigma = compute_composition_background(path, 0, "A")

        self.assertAlmostEqual(bg_mu, 1 / 3)
        self.assertAlmostEqual(bg_sigma, math.sqrt(2) / 3)

    def test_each_normalized_channel_has_a_background_function(self):
        path = self.write_fasta(
            ">proteome_part_1\nACDEFGHIKLMN\n"
            ">proteome_part_2\nPQRSTVWYACDE\n"
        )
        calculations = [
            lambda: compute_composition_background(path, 1, "AG"),
            lambda: compute_net_charge_background(path, 1),
            lambda: compute_hydropathy_background(path, 1),
            lambda: compute_hydrophobic_moment_background(path, 1),
            lambda: compute_disorder_background(path, 1),
            lambda: compute_entropy_background(path, 1),
            lambda: compute_dipeptide_background(path, 1, "AC"),
        ]

        for calculation in calculations:
            with self.subTest(calculation=calculation):
                bg_mu, bg_sigma = calculation()
                self.assertTrue(math.isfinite(bg_mu))
                self.assertGreater(bg_sigma, 0.0)

    def test_computed_background_standardizes_its_source_values(self):
        path = self.write_fasta(">one\nA\n>two\nCC\n")
        bg_mu, bg_sigma = compute_composition_background(path, 0, "A")
        channel = mk_composition_bias(0, "A", bg_mu, bg_sigma)

        values = channel("A") + channel("CC")

        self.assertAlmostEqual(sum(values) / len(values), 0.0)
        self.assertAlmostEqual(
            math.sqrt(sum(value * value for value in values) / len(values)),
            1.0,
        )

    def test_dipeptide_background_skips_short_records(self):
        path = self.write_fasta(">short\nA\n>usable\nAACA\n")

        bg_mu, bg_sigma = compute_dipeptide_background(path, 1, "AA")

        self.assertAlmostEqual(bg_mu, 0.375)
        self.assertAlmostEqual(bg_sigma, math.sqrt(0.171875))

    def test_dipeptide_background_raises_if_no_record_is_usable(self):
        path = self.write_fasta(">one\nA\n>two\nC\n")

        with self.assertRaisesRegex(ValueError, "no usable positions"):
            compute_dipeptide_background(path, 1, "AA")

    def test_zero_variance_background_raises(self):
        path = self.write_fasta(">one\nAAAA\n")

        with self.assertRaisesRegex(ValueError, "standard deviation"):
            compute_composition_background(path, 1, "A")

    def test_noncanonical_fasta_residue_raises(self):
        path = self.write_fasta(">one\nAAX\n")

        with self.assertRaisesRegex(ValueError, "noncanonical.*X"):
            compute_hydropathy_background(path, 1)

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
                    compute_hydropathy_background(path, 1)
