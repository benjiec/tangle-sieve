import csv
import os
import sys
import tempfile
import unittest

from sieve.artifacts import rule_results_tsv, sequences_fasta
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

    def write_filter_module(self, tmpd, lines):
        module_path = os.path.join(tmpd, "positive_filters.py")
        self.write_lines(module_path, lines)
        return "positive_filters.is_positive"

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

    def test_duplicate_rule_result_sequence_accessions_use_first_representative(self):
        with tempfile.TemporaryDirectory() as tmpd:
            hmm = os.path.join(tmpd, "hmm.tsv")
            artifacts = os.path.join(tmpd, "artifacts")
            os.makedirs(artifacts)
            self.write_tsv(hmm, [
                "sequence accession",
                "domain bitscore",
            ], [
                {"sequence accession": "same_accession", "domain bitscore": "50"},
            ])
            self.write_tsv(rule_results_tsv(artifacts), [
                "protein accession",
                "sequence accession",
                "genome accession",
                "Rule",
            ], [
                {"protein accession": "same_accession", "sequence accession": "same_accession", "genome accession": "g1", "Rule": "true"},
                {"protein accession": "same_accession", "sequence accession": "same_accession", "genome accession": "g2", "Rule": "false"},
            ])
            write_fasta_from_dict({"same_accession": "MA"}, sequences_fasta(artifacts))
            filter_spec = self.write_filter_module(tmpd, [
                "from sieve.result_filters import Field",
                "is_positive = Field('Rule').eq('true')",
            ])

            sys.path.insert(0, tmpd)
            try:
                rows = self.script.discover_threshold(hmm, artifacts, filter_spec)
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("positive_filters", None)

        self.assertEqual(rows[0]["tp"], 1)
        self.assertEqual(rows[0]["fp"], 0)

    def test_discovers_threshold_stats_and_marks_best_row(self):
        with tempfile.TemporaryDirectory() as tmpd:
            hmm = os.path.join(tmpd, "hmm.tsv")
            artifacts = os.path.join(tmpd, "artifacts")
            os.makedirs(artifacts)
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
            self.write_tsv(rule_results_tsv(artifacts), [
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
            }, sequences_fasta(artifacts))
            filter_spec = self.write_filter_module(tmpd, [
                "from sieve.result_filters import FieldRegex",
                "is_positive = FieldRegex(r\"Pfam\\.matches.+\").any().eq(\"true\")",
            ])

            sys.path.insert(0, tmpd)
            try:
                rows = self.script.discover_threshold(hmm, artifacts, filter_spec)
                self.script.write_threshold_stats(rows, output)
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("positive_filters", None)

            with open(output, "r", encoding="utf-8", newline="") as f:
                lines = f.readlines()
            output_rows = list(csv.DictReader([
                line for line in lines
                if not line.startswith("# ")
            ], delimiter="\t"))
            comments = [line.rstrip("\n") for line in lines if line.startswith("# ")]

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
        self.assertEqual(comments, [
            "# selected threshold bitscore 80",
            "# false positives",
            "# false negatives",
            "# p3 20",
        ])

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
            artifacts = os.path.join(tmpd, "artifacts")
            os.makedirs(artifacts)

            self.write_tsv(hmm, [
                "sequence accession",
                "domain bitscore",
            ], [
                {"sequence accession": "p1", "domain bitscore": "50"},
                {"sequence accession": "p2", "domain bitscore": "10"},
            ])
            self.write_tsv(rule_results_tsv(artifacts), [
                "protein accession",
                "sequence accession",
                "genome accession",
                "Rule",
            ], [
                {"protein accession": "p1", "sequence accession": "p1", "genome accession": "g1", "Rule": "true"},
                {"protein accession": "p2", "sequence accession": "p2", "genome accession": "g2", "Rule": "false"},
            ])
            write_fasta_from_dict({"p1": "MA", "p2": "MG"}, sequences_fasta(artifacts))
            filter_spec = self.write_filter_module(tmpd, [
                "from sieve.result_filters import Field",
                "is_positive = Field('Rule').eq('true')",
            ])

            sys.path.insert(0, tmpd)
            try:
                rows = self.script.discover_threshold(hmm, artifacts, filter_spec, taxon="cnidaria")
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("positive_filters", None)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["tp"], 1)
        self.assertEqual(rows[0]["fp"], 0)
        self.assertEqual(rows[0]["tn"], 0)
        self.assertEqual(rows[0]["fn"], 0)

    def test_selected_threshold_error_details_include_false_positives(self):
        with tempfile.TemporaryDirectory() as tmpd:
            output = os.path.join(tmpd, "thresholds.tsv")
            rows = [{
                "threshold bitscore": 10.0,
                "tp": 1,
                "fp": 1,
                "tn": 0,
                "fn": 0,
                "sensitivity": 1.0,
                "specificity": 0.0,
                "balanced accuracy": 0.5,
                "selected": "true",
                "false positives": [{"sequence accession": "fp1", "bitscore": 12.5}],
                "false negatives": [],
            }]

            self.script.write_threshold_stats(rows, output)

            with open(output, "r", encoding="utf-8") as f:
                comments = [line.rstrip("\n") for line in f if line.startswith("# ")]

        self.assertEqual(comments, [
            "# selected threshold bitscore 10",
            "# false positives",
            "# fp1 12.5",
            "# false negatives",
        ])

    def test_main_uses_standard_artifact_paths(self):
        with tempfile.TemporaryDirectory() as tmpd:
            hmm = os.path.join(tmpd, "hmm.tsv")
            artifacts = os.path.join(tmpd, "artifacts")
            os.makedirs(artifacts)
            output = os.path.join(tmpd, "thresholds.tsv")
            self.write_tsv(hmm, [
                "sequence accession",
                "domain bitscore",
            ], [
                {"sequence accession": "p1", "domain bitscore": "50"},
            ])
            self.write_tsv(rule_results_tsv(artifacts), [
                "protein accession",
                "sequence accession",
                "genome accession",
                "Rule",
            ], [
                {"protein accession": "p1", "sequence accession": "p1", "genome accession": "g1", "Rule": "true"},
            ])
            write_fasta_from_dict({"p1": "MA"}, sequences_fasta(artifacts))
            filter_spec = self.write_filter_module(tmpd, [
                "from sieve.result_filters import Field",
                "is_positive = Field('Rule').eq('true')",
            ])

            sys.path.insert(0, tmpd)
            try:
                self.script.main([
                    "--hmmsearch-tsv", hmm,
                    "--artifacts-dir", artifacts,
                    "--positive-filter", filter_spec,
                    "--output", output,
                ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("positive_filters", None)

            with open(output, "r", encoding="utf-8", newline="") as f:
                output_rows = list(csv.DictReader(f, delimiter="\t"))

        self.assertEqual(output_rows[0]["threshold bitscore"], "50")
        self.assertEqual(output_rows[0]["tp"], "1")


if __name__ == "__main__":
    unittest.main()
