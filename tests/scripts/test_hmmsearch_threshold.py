import csv
import os
import tempfile
import unittest

from tangle.sequence import write_fasta_from_dict

from tests.fixtures import DefaultsFixture
from tests.scripts.helpers import load_script


class TestHmmsearchThresholdScript(unittest.TestCase):

    def setUp(self):
        self.fx = DefaultsFixture(self)
        self.repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.script = load_script(os.path.join(self.repo, "scripts", "hmmsearch-threshold.py"))

    def tearDown(self):
        self.fx.cleanup()

    def write_tsv(self, path, headers, rows):
        with open(path, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def write_lines(self, path, lines):
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    def test_positive_spec_uses_literal_and_numeric_operators(self):
        with tempfile.TemporaryDirectory() as tmpd:
            spec = os.path.join(tmpd, "positive.tsv")
            self.write_lines(spec, [
                "column_regex\toperator\tvalue",
                "Literal\teq\t1",
                "Numeric\tnum_eq\t1",
                "Leader.call\\('mTP'\\)\tgt\t1",
            ])

            specs = self.script.read_positive_specs(spec)

        self.assertFalse(self.script.value_matches("001", "eq", "1"))
        self.assertTrue(self.script.value_matches("001", "num_eq", "1"))
        self.assertTrue(self.script.row_is_positive({
            "Literal": "001",
            "Numeric": "001",
            "Leader.call('mTP')": "0",
        }, specs))
        self.assertTrue(self.script.row_is_positive({
            "Literal": "0",
            "Numeric": "0",
            "Leader.call('mTP')": "2",
        }, specs))

    def test_keeps_best_scoring_sequence_per_protein_accession(self):
        entries = [
            {
                "sequence accession": "p1_with_leader_a",
                "protein accession": "p1",
                "bitscore": 40.0,
            },
            {
                "sequence accession": "p1_with_leader_b",
                "protein accession": "p1",
                "bitscore": 80.0,
            },
            {
                "sequence accession": "p2",
                "protein accession": "p2",
                "bitscore": 30.0,
            },
        ]

        kept = self.script.best_entry_per_protein(entries)

        self.assertEqual(
            [entry["sequence accession"] for entry in kept],
            ["p1_with_leader_b", "p2"],
        )

    def test_discovers_threshold_stats_and_marks_best_row(self):
        with tempfile.TemporaryDirectory() as tmpd:
            hmm = os.path.join(tmpd, "hmm.tsv")
            rules = os.path.join(tmpd, "rules.tsv")
            fasta = os.path.join(tmpd, "sequences.faa")
            spec = os.path.join(tmpd, "positive.tsv")
            output = os.path.join(tmpd, "thresholds.tsv")

            self.write_tsv(hmm, [
                "sequence accession",
                "HMM model",
                "domain e-value",
                "domain bitscore",
            ], [
                {"sequence accession": "p1_a", "HMM model": "m", "domain e-value": "1e-5", "domain bitscore": "40"},
                {"sequence accession": "p1_b", "HMM model": "m", "domain e-value": "1e-8", "domain bitscore": "80"},
                {"sequence accession": "p2", "HMM model": "m", "domain e-value": "1e-4", "domain bitscore": "30"},
                {"sequence accession": "p3", "HMM model": "m", "domain e-value": "1e-3", "domain bitscore": "20"},
            ])
            self.write_tsv(rules, [
                "protein accession",
                "sequence accession",
                "genome accession",
                "Pfam.matches('PF00001')",
            ], [
                {"protein accession": "p1", "sequence accession": "p1_a", "genome accession": "g1", "Pfam.matches('PF00001')": "true"},
                {"protein accession": "p1", "sequence accession": "p1_b", "genome accession": "g1", "Pfam.matches('PF00001')": "true"},
                {"protein accession": "p2", "sequence accession": "p2", "genome accession": "g2", "Pfam.matches('PF00001')": "false"},
                {"protein accession": "p3", "sequence accession": "p3", "genome accession": "g3", "Pfam.matches('PF00001')": "true"},
                {"protein accession": "p4", "sequence accession": "p4", "genome accession": "g4", "Pfam.matches('PF00001')": "false"},
            ])
            write_fasta_from_dict({
                "p1_a": "MA",
                "p1_b": "MMA",
                "p2": "MG",
                "p3": "MP",
                "p4": "MT",
            }, fasta)
            self.write_lines(spec, ["Pfam\\.matches.+\teq\ttrue"])

            rows = self.script.discover_threshold(hmm, rules, fasta, spec)
            self.script.write_threshold_stats(rows, output)

            with open(output, "r", encoding="utf-8", newline="") as f:
                output_rows = list(csv.DictReader(f, delimiter="\t"))

        self.assertEqual(
            [(row["threshold bitscore"], row["tp"], row["fp"], row["tn"], row["fn"], row["selected"]) for row in output_rows],
            [
                ("80", "1", "0", "2", "1", "true"),
                ("30", "1", "1", "1", "1", ""),
                ("20", "2", "1", "1", "0", ""),
                ("0", "2", "2", "0", "0", ""),
            ],
        )
        self.assertEqual(output_rows[0]["sensitivity"], "0.5")
        self.assertEqual(output_rows[0]["specificity"], "1")
        self.assertEqual(output_rows[0]["balanced accuracy"], "0.75")

    def test_taxon_filters_entries_before_computing_stats(self):
        self.fx.write_taxonomy_rows([
            {
                "Genome Accession": "g1",
                "Domain": "Eukaryota",
                "Phylum": "Cnidaria",
            },
            {
                "Genome Accession": "g2",
                "Domain": "Eukaryota",
                "Phylum": "Arthropoda",
            },
        ])
        with tempfile.TemporaryDirectory() as tmpd:
            hmm = os.path.join(tmpd, "hmm.tsv")
            rules = os.path.join(tmpd, "rules.tsv")
            fasta = os.path.join(tmpd, "sequences.faa")
            spec = os.path.join(tmpd, "positive.tsv")

            self.write_tsv(hmm, [
                "sequence accession",
                "domain bitscore",
            ], [
                {"sequence accession": "p1", "domain bitscore": "50"},
                {"sequence accession": "p2", "domain bitscore": "10"},
            ])
            self.write_tsv(rules, [
                "protein accession",
                "sequence accession",
                "genome accession",
                "Rule",
            ], [
                {"protein accession": "p1", "sequence accession": "p1", "genome accession": "g1", "Rule": "true"},
                {"protein accession": "p2", "sequence accession": "p2", "genome accession": "g2", "Rule": "false"},
            ])
            write_fasta_from_dict({"p1": "MA", "p2": "MG"}, fasta)
            self.write_lines(spec, ["Rule\teq\ttrue"])

            rows = self.script.discover_threshold(hmm, rules, fasta, spec, taxon="cnidaria")

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["tp"], 1)
        self.assertEqual(rows[0]["fp"], 0)
        self.assertEqual(rows[0]["tn"], 0)
        self.assertEqual(rows[0]["fn"], 0)

    def test_reports_missing_columns_and_invalid_numeric_values(self):
        with tempfile.TemporaryDirectory() as tmpd:
            spec = os.path.join(tmpd, "positive.tsv")
            self.write_lines(spec, ["Missing\teq\ttrue"])
            specs = self.script.read_positive_specs(spec)

        with self.assertRaisesRegex(ValueError, "did not match any"):
            self.script.row_is_positive({"Rule": "true"}, specs)
        with self.assertRaisesRegex(ValueError, "Expected numeric value"):
            self.script.value_matches("not-a-number", "gt", "1")


if __name__ == "__main__":
    unittest.main()
