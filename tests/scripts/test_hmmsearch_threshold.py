import csv
import os
import tempfile
import unittest
from subprocess import CompletedProcess
from unittest.mock import patch

from sieve.artifacts import SEQUENCE_HEADERS, rule_results_tsv, sequences_fasta, sequences_tsv
from sieve.hmmsearch import DomtbloutHit
from tangle.sequence import write_fasta_from_dict

from tests.scripts.helpers import load_script


class TestHmmsearchThresholdScript(unittest.TestCase):

    def setUp(self):
        repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.script = load_script(os.path.join(repo, "scripts", "hmmsearch-threshold.py"))

    def write_artifacts(self, artifacts, sequences, rule_rows):
        os.makedirs(artifacts)
        write_fasta_from_dict(sequences, sequences_fasta(artifacts))
        with open(sequences_tsv(artifacts), "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=SEQUENCE_HEADERS, delimiter="\t")
            writer.writeheader()
            writer.writerows({
                "protein accession": row["protein accession"],
                "genome accession": row.get("genome accession", ""),
                "sequence accession": accession,
                "start label": "",
                "start aa 1b": "1",
            } for accession, row in zip(sequences, rule_rows))
        with open(rule_results_tsv(artifacts), "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["protein accession", "sequence accession", "genome accession", "pass all"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(rule_rows)

    def hit(self, accession, bitscore, domain_bitscore=None):
        if domain_bitscore is None:
            domain_bitscore = bitscore
        return DomtbloutHit(
            sequence_accession=accession,
            model_accession="model",
            model_name="model",
            full_evalue=1e-10,
            full_bitscore=bitscore,
            domain_evalue=1e-10,
            domain_bitscore=domain_bitscore,
            hmm_start=1,
            hmm_end=10,
            sequence_start=1,
            sequence_end=10,
        )

    def test_best_hmm_hits_keeps_highest_full_sequence_score(self):
        best = self.script.best_hmm_hits([
            self.hit("p1", 20, domain_bitscore=90),
            self.hit("p1", 40, domain_bitscore=10),
            self.hit("p2", 30, domain_bitscore=80),
        ])
        self.assertEqual(best["p1"]["bitscore"], 40)
        self.assertEqual(best["p2"]["bitscore"], 30)

    def test_join_labels_accessions_in_positive_fasta_as_positive(self):
        entries = self.script.joined_entries(
            [self.hit("p1_a", 40), self.hit("p1_b", 80)],
            {"p1_a": "MA"},
            {"p1_a": "MA", "p1_b": "MMA"},
        )
        rows = self.script.threshold_stats(entries)

        self.assertEqual(len(entries), 2)
        self.assertEqual([entry["positive"] for entry in entries], [True, False])
        self.assertEqual((rows[0]["tp"], rows[0]["fp"], rows[0]["fn"]), (0, 1, 1))

    def test_artifact_entries_still_use_pass_all_labels(self):
        with tempfile.TemporaryDirectory() as tmpd:
            artifacts = os.path.join(tmpd, "artifacts")
            rule_rows = [
                {"protein accession": "p1", "sequence accession": "p1_a", "genome accession": "g1", "pass all": "true"},
                {"protein accession": "p1", "sequence accession": "p1_b", "genome accession": "g1", "pass all": "false"},
            ]
            self.write_artifacts(artifacts, {"p1_a": "MA", "p1_b": "MMA"}, rule_rows)
            entries = self.script.joined_artifact_entries(
                [self.hit("p1_a", 40), self.hit("p1_b", 80)],
                rule_rows,
                self.script.read_sequence_rows(artifacts),
                sequences_fasta(artifacts),
            )

        self.assertEqual([entry["positive"] for entry in entries], [True, False])

    def test_discover_threshold_runs_hmmsearch_without_cut_ga_and_scores_missing_hits_zero(self):
        with tempfile.TemporaryDirectory() as tmpd:
            profile = os.path.join(tmpd, "profile.hmm")
            positives = os.path.join(tmpd, "positives.faa")
            all_fasta = os.path.join(tmpd, "all.faa")
            output = os.path.join(tmpd, "thresholds.tsv")
            with open(profile, "w", encoding="utf-8"):
                pass
            write_fasta_from_dict({"p1": "MA"}, positives)
            write_fasta_from_dict({"p1": "MA", "p2": "MG"}, all_fasta)
            commands = []
            searched_sequences = []

            def run(cmd, check, capture_output, text):
                commands.append(cmd)
                searched_sequences.append(
                    self.script.read_fasta_ignoring_duplicate_accessions(cmd[-1])
                )
                domtblout = cmd[cmd.index("--domtblout") + 1]
                line = " ".join([
                    "p1", "-", "100", "model", "MODEL", "100", "1e-10", "50", "0",
                    "1", "1", "1e-10", "1e-10", "5", "0", "1", "10", "1", "10",
                    "1", "10", "0.99",
                ])
                with open(domtblout, "w", encoding="utf-8") as f:
                    f.write(line + "\n")
                return CompletedProcess(cmd, 0, stdout="", stderr="")

            with patch("subprocess.run", side_effect=run):
                rows, domtblout = self.script.discover_threshold(
                    profile,
                    positive_fasta=positives,
                    all_fasta=all_fasta,
                )
                self.script.write_threshold_stats(rows, output, domtblout)

            self.assertNotIn("--cut_ga", commands[0])
            self.assertEqual(commands[0][-2], profile)
            self.assertEqual(searched_sequences[0], {"p1": "MA", "p2": "MG"})
            self.assertEqual([row["threshold bitscore"] for row in rows], [50.0, 0.0])
            with open(output, "r", encoding="utf-8") as f:
                text = f.read()
            self.assertIn("50\t1\t0\t1\t0", text)
            self.assertIn("# hmmsearch domtblout\n", text)
            self.assertIn("# p1 - 100 model MODEL 100 1e-10 50", text)

    def test_rejects_positive_accession_absent_from_all_fasta(self):
        with self.assertRaisesRegex(ValueError, "absent from all FASTA: p3"):
            self.script.joined_entries([], {"p1": "MA", "p3": "MX"}, {"p1": "MA", "p2": "MG"})

    def test_rejects_positive_sequence_mismatch(self):
        with self.assertRaisesRegex(ValueError, "sequences disagree: p1"):
            self.script.joined_entries([], {"p1": "MX"}, {"p1": "MA", "p2": "MG"})

    def test_requires_positive_and_negative_sequences(self):
        with self.assertRaisesRegex(ValueError, "Positive FASTA contains no sequences"):
            self.script.joined_entries([], {}, {"p1": "MA"})
        with self.assertRaisesRegex(ValueError, "no negative sequences"):
            self.script.joined_entries([], {"p1": "MA"}, {"p1": "MA"})

    def test_duplicate_accessions_keep_first_sequence(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "duplicate.faa")
            with open(fasta, "w", encoding="utf-8") as f:
                f.write(">p1 first\nMA\n>p2\nMG\n>p1 repeated\nMX\n")
            self.assertEqual(
                self.script.read_fasta_ignoring_duplicate_accessions(fasta),
                {"p1": "MA", "p2": "MG"},
            )

    def test_main_accepts_fasta_interface(self):
        with (
            patch.object(self.script, "discover_threshold", return_value=([], "")) as discover,
            patch.object(self.script, "write_threshold_stats") as write,
        ):
            result = self.script.main([
                "--hmm", "model.hmm",
                "--positive-fasta", "positive.faa",
                "--all-fasta", "all.faa",
                "--output", "thresholds.tsv",
            ])

        self.assertEqual(result, 0)
        discover.assert_called_once_with(
            "model.hmm",
            artifacts_dir=None,
            positive_fasta="positive.faa",
            all_fasta="all.faa",
        )
        write.assert_called_once_with([], "thresholds.tsv", "")

    def test_main_accepts_artifacts_interface(self):
        with (
            patch.object(self.script, "discover_threshold", return_value=([], "")) as discover,
            patch.object(self.script, "write_threshold_stats"),
        ):
            result = self.script.main([
                "--hmm", "model.hmm",
                "--artifacts-dir", "artifacts",
            ])

        self.assertEqual(result, 0)
        discover.assert_called_once_with(
            "model.hmm",
            artifacts_dir="artifacts",
            positive_fasta=None,
            all_fasta=None,
        )

    def test_main_rejects_mixed_or_incomplete_input_modes(self):
        invalid_arguments = [
            ["--hmm", "model.hmm"],
            ["--hmm", "model.hmm", "--positive-fasta", "positive.faa"],
            [
                "--hmm", "model.hmm",
                "--artifacts-dir", "artifacts",
                "--positive-fasta", "positive.faa",
                "--all-fasta", "all.faa",
            ],
        ]
        for arguments in invalid_arguments:
            with self.subTest(arguments=arguments), self.assertRaises(SystemExit):
                self.script.main(arguments)


if __name__ == "__main__":
    unittest.main()
