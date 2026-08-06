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

    def write_tsv(self, path, headers, rows):
        with open(path, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def write_artifacts(self, artifacts, sequences, rule_rows):
        os.makedirs(artifacts)
        write_fasta_from_dict(sequences, sequences_fasta(artifacts))
        self.write_tsv(sequences_tsv(artifacts), SEQUENCE_HEADERS, [
            {
                "protein accession": row["protein accession"],
                "genome accession": row.get("genome accession", ""),
                "sequence accession": accession,
                "start label": "",
                "start aa 1b": "1",
            }
            for accession, row in zip(sequences, rule_rows)
        ])
        self.write_tsv(
            rule_results_tsv(artifacts),
            ["protein accession", "sequence accession", "genome accession", "pass all"],
            rule_rows,
        )

    def hit(self, accession, bitscore):
        return DomtbloutHit(
            sequence_accession=accession,
            model_accession="model",
            model_name="model",
            full_evalue=1e-10,
            full_bitscore=bitscore,
            domain_evalue=1e-10,
            domain_bitscore=bitscore,
            hmm_start=1,
            hmm_end=10,
            sequence_start=1,
            sequence_end=10,
        )

    def test_best_hmm_hits_keeps_highest_domain_score(self):
        best = self.script.best_hmm_hits([
            self.hit("p1", 20), self.hit("p1", 40), self.hit("p2", 30),
        ])
        self.assertEqual(best["p1"]["bitscore"], 40)
        self.assertEqual(best["p2"]["bitscore"], 30)

    def test_join_labels_each_candidate_from_pass_all_without_protein_collapse(self):
        with tempfile.TemporaryDirectory() as tmpd:
            artifacts = os.path.join(tmpd, "artifacts")
            rule_rows = [
                {"protein accession": "p1", "sequence accession": "p1_a", "genome accession": "g1", "pass all": "true"},
                {"protein accession": "p1", "sequence accession": "p1_b", "genome accession": "g1", "pass all": "false"},
            ]
            self.write_artifacts(artifacts, {"p1_a": "MA", "p1_b": "MMA"}, rule_rows)
            entries = self.script.joined_entries(
                [self.hit("p1_a", 40), self.hit("p1_b", 80)],
                rule_rows,
                self.script.read_sequence_rows(artifacts),
                sequences_fasta(artifacts),
            )
            rows = self.script.threshold_stats(entries)

        self.assertEqual(len(entries), 2)
        self.assertEqual([entry["positive"] for entry in entries], [True, False])
        self.assertEqual((rows[0]["tp"], rows[0]["fp"], rows[0]["fn"]), (0, 1, 1))

    def test_discover_threshold_runs_hmmsearch_without_cut_ga_and_scores_missing_hits_zero(self):
        with tempfile.TemporaryDirectory() as tmpd:
            artifacts = os.path.join(tmpd, "artifacts")
            profile = os.path.join(tmpd, "profile.hmm")
            output = os.path.join(tmpd, "thresholds.tsv")
            with open(profile, "w", encoding="utf-8"):
                pass
            rule_rows = [
                {"protein accession": "p1", "sequence accession": "p1", "genome accession": "g1", "pass all": "true"},
                {"protein accession": "p2", "sequence accession": "p2", "genome accession": "g2", "pass all": "false"},
            ]
            self.write_artifacts(artifacts, {"p1": "MA", "p2": "MG"}, rule_rows)
            commands = []

            def run(cmd, check, capture_output, text):
                commands.append(cmd)
                domtblout = cmd[cmd.index("--domtblout") + 1]
                line = " ".join([
                    "p1", "-", "100", "model", "MODEL", "100", "1e-10", "50", "0",
                    "1", "1", "1e-10", "1e-10", "50", "0", "1", "10", "1", "10",
                    "1", "10", "0.99",
                ])
                with open(domtblout, "w", encoding="utf-8") as f:
                    f.write(line + "\n")
                return CompletedProcess(cmd, 0, stdout="", stderr="")

            with patch("subprocess.run", side_effect=run):
                rows = self.script.discover_threshold(profile, artifacts)
                self.script.write_threshold_stats(rows, output)

            self.assertNotIn("--cut_ga", commands[0])
            self.assertEqual(commands[0][-2:], [profile, sequences_fasta(artifacts)])
            self.assertEqual([row["threshold bitscore"] for row in rows], [50.0, 0.0])
            with open(output, "r", encoding="utf-8") as f:
                text = f.read()
            self.assertIn("50\t1\t0\t1\t0", text)

    def test_requires_one_rule_result_per_candidate(self):
        with tempfile.TemporaryDirectory() as tmpd:
            artifacts = os.path.join(tmpd, "artifacts")
            rule_rows = [
                {"protein accession": "p1", "sequence accession": "p1", "genome accession": "", "pass all": "true"},
            ]
            self.write_artifacts(artifacts, {"p1": "MA", "p2": "MG"}, rule_rows)
            with self.assertRaisesRegex(ValueError, "no metadata"):
                self.script.joined_entries(
                    [], rule_rows, self.script.read_sequence_rows(artifacts), sequences_fasta(artifacts),
                )


if __name__ == "__main__":
    unittest.main()
