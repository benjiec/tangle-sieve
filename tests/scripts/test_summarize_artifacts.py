import csv
import io
import os
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from tangle.sequence import write_fasta_from_dict

from tests.scripts.helpers import load_script


class TestSummarizeArtifactsScript(unittest.TestCase):

    def setUp(self):
        self.repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.script = load_script(os.path.join(self.repo, "scripts", "summarize-artifacts.py"))

    def write_artifacts(self, artifacts):
        os.makedirs(artifacts)
        with open(os.path.join(artifacts, "rule-results.tsv"), "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=[
                "rule result",
                "pass all",
                "contig accession",
                "sequence accession",
                "protein accession",
                "genome accession",
            ], delimiter="\t")
            writer.writeheader()
            writer.writerows([
                {"protein accession": "p1|description_one", "sequence accession": "seq1|description_one", "genome accession": "", "contig accession": "c1", "pass all": "true", "rule result": "x"},
                {"protein accession": "p2", "sequence accession": "seq2", "genome accession": "g2", "contig accession": "c2", "pass all": "false", "rule result": "y"},
                {"protein accession": "p3|description_ten", "sequence accession": "seq10|description_ten", "genome accession": "", "contig accession": "c3", "pass all": "true", "rule result": "z"},
            ])
        write_fasta_from_dict({"seq1|description_one": "MAAA", "seq2": "MBBB", "seq10|description_ten": "MCCC"}, os.path.join(artifacts, "sequences.faa"))

    def fake_hmmsearch(self, expected_hmm, expected_fasta, lines, returncode=0):
        def run(cmd, check, capture_output, text):
            self.assertEqual(cmd[0:2], ["hmmsearch", "--domtblout"])
            self.assertEqual(cmd[3:], [expected_hmm, expected_fasta])
            self.assertFalse(check)
            self.assertTrue(capture_output)
            self.assertTrue(text)
            if returncode == 0:
                with open(cmd[2], "w", encoding="utf-8") as f:
                    f.write("\n".join(lines))
            return subprocess.CompletedProcess(cmd, returncode, stdout="summary\n", stderr="bad hmm\n" if returncode else "")
        return run

    def test_writes_named_columns_and_highest_score_for_exact_sequence_accession(self):
        with tempfile.TemporaryDirectory() as tmpd:
            artifacts = os.path.join(tmpd, "artifacts")
            self.write_artifacts(artifacts)
            hmm = os.path.join(tmpd, "models", "example.hmm")
            output = os.path.join(tmpd, "summary.tsv")
            lines = [
                "seq1|description_one - 4 model - 10 1e-5 8 0 1 2 1e-4 1e-4 7.5 0 1 4 1 4 1 4 0.9 hit",
                "seq1|description_one - 4 model - 10 1e-5 9 0 2 2 1e-5 1e-5 8.25 0 5 8 1 4 1 4 0.9 hit",
                "seq10|description_ten - 4 model - 10 1e-5 7 0 1 1 1e-3 1e-3 6.0 0 1 4 1 4 1 4 0.9 hit",
            ]
            fake_run = self.fake_hmmsearch(hmm, os.path.join(artifacts, "sequences.faa"), lines)

            with patch.object(self.script.subprocess, "run", side_effect=fake_run):
                self.script.main(["--artifacts-dir", artifacts, "--category", "algae", "--hmm", hmm, "--output", output])

            with open(output, "r", encoding="utf-8", newline="") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(rows, [
                {"protein accession": "p1|description_one", "sequence accession": "seq1|description_one", "genome accession": "", "contig accession": "c1", "pass all": "true", "category": "algae", "hmm model": "example.hmm", "highest domain score": "8.25", "genome description": "description_one"},
                {"protein accession": "p2", "sequence accession": "seq2", "genome accession": "g2", "contig accession": "c2", "pass all": "false", "category": "algae", "hmm model": "example.hmm", "highest domain score": "", "genome description": ""},
                {"protein accession": "p3|description_ten", "sequence accession": "seq10|description_ten", "genome accession": "", "contig accession": "c3", "pass all": "true", "category": "algae", "hmm model": "example.hmm", "highest domain score": "6.0", "genome description": "description_ten"},
            ])

    def test_appends_without_repeating_header_and_rejects_wrong_header(self):
        with tempfile.TemporaryDirectory() as tmpd:
            output = os.path.join(tmpd, "summary.tsv")
            row = {column: column for column in self.script.RULE_RESULT_COLUMNS}
            self.script.write_summary([row], "one", "/a/first.hmm", {}, output)
            self.script.write_summary([row], "two", "/b/second.hmm", {}, output)
            with open(output, "r", encoding="utf-8", newline="") as f:
                lines = f.readlines()
            self.assertEqual(len(lines), 3)
            self.assertEqual(lines[0].rstrip("\r\n").split("\t"), self.script.SUMMARY_COLUMNS)

            bad_output = os.path.join(tmpd, "bad.tsv")
            with open(bad_output, "w", encoding="utf-8") as f:
                f.write("wrong\theader\n")
            with self.assertRaisesRegex(ValueError, "Unexpected columns"):
                self.script.write_summary([], "x", "model.hmm", {}, bad_output)

    def test_rejects_missing_required_rule_result_column(self):
        with tempfile.TemporaryDirectory() as tmpd:
            with open(os.path.join(tmpd, "rule-results.tsv"), "w", encoding="utf-8") as f:
                f.write("protein accession\tsequence accession\n")
            with self.assertRaisesRegex(ValueError, "genome accession"):
                self.script.read_rule_results(tmpd)

    def test_reports_hmmsearch_failure(self):
        with tempfile.TemporaryDirectory() as tmpd:
            artifacts = os.path.join(tmpd, "artifacts")
            self.write_artifacts(artifacts)
            hmm = os.path.join(tmpd, "example.hmm")
            fake_run = self.fake_hmmsearch(hmm, os.path.join(artifacts, "sequences.faa"), [], returncode=2)
            stderr = io.StringIO()
            with patch.object(self.script.subprocess, "run", side_effect=fake_run), patch("sys.stderr", stderr):
                with self.assertRaises(SystemExit) as cm:
                    self.script.highest_domain_scores(artifacts, hmm)
            self.assertEqual(cm.exception.code, 2)
            self.assertIn("hmmsearch failed:", stderr.getvalue())
            self.assertIn("bad hmm", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
