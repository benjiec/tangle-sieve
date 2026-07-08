import csv
import io
import os
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from tangle.sequence import write_fasta_from_dict

from tests.scripts.helpers import load_script


class TestHmmsearchSequencesScript(unittest.TestCase):

    def setUp(self):
        self.repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    def test_runs_hmmsearch_against_artifact_sequences_and_writes_domain_tsv(self):
        script = load_script(os.path.join(self.repo, "scripts", "hmmsearch-sequences.py"))

        def fake_run(cmd, check, capture_output, text):
            self.assertEqual(cmd[0], "hmmsearch")
            self.assertEqual(cmd[1], "--domtblout")
            self.assertEqual(cmd[3:], [hmm_file, os.path.join(artifacts, "sequences.faa")])
            self.assertFalse(check)
            self.assertTrue(capture_output)
            self.assertTrue(text)
            with open(cmd[2], "w", encoding="utf-8") as f:
                f.write("\n".join([
                    "# target name accession tlen query name accession qlen E-value score bias # of c-Evalue i-Evalue score bias hmm from hmm to ali from ali to env from env to acc description",
                    "seq1 - 120 modelA - 80 1e-40 140.0 0.0 1 1 2e-20 3e-22 75.5 0.1 4 55 10 61 8 63 0.95 first hit",
                    "seq2 - 90 modelA - 80 1e-10 70.0 0.0 1 2 5e-05 7e-06 42.0 0.0 1 30 20 49 18 50 0.80 second hit",
                    "",
                ]))
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

        with tempfile.TemporaryDirectory() as tmpd:
            artifacts = os.path.join(tmpd, "artifacts")
            os.makedirs(artifacts)
            write_fasta_from_dict({"seq1": "MSEQ", "seq2": "MPEP"}, os.path.join(artifacts, "sequences.faa"))
            hmm_file = os.path.join(tmpd, "model.hmm")
            output = os.path.join(tmpd, "hits.tsv")
            with open(hmm_file, "w", encoding="utf-8") as f:
                f.write("HMMER3/f\n")

            with patch.object(script.subprocess, "run", side_effect=fake_run):
                rows = script.run_hmmsearch(artifacts, hmm_file, output)

            self.assertEqual(rows, [
                {
                    "sequence accession": "seq1",
                    "HMM model": "modelA",
                    "domain e-value": "3e-22",
                    "domain bitscore": "75.5",
                    "query start": "10",
                    "query end": "61",
                    "hmm start": "4",
                    "hmm end": "55",
                },
                {
                    "sequence accession": "seq2",
                    "HMM model": "modelA",
                    "domain e-value": "7e-06",
                    "domain bitscore": "42.0",
                    "query start": "20",
                    "query end": "49",
                    "hmm start": "1",
                    "hmm end": "30",
                },
            ])
            with open(output, "r", encoding="utf-8", newline="") as f:
                tsv_rows = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(tsv_rows, rows)

    def test_prints_hmmsearch_stderr_on_failure(self):
        script = load_script(os.path.join(self.repo, "scripts", "hmmsearch-sequences.py"))

        def fake_run(cmd, check, capture_output, text):
            return subprocess.CompletedProcess(cmd, 1, stdout="summary\n", stderr="bad hmm\n")

        with tempfile.TemporaryDirectory() as tmpd:
            artifacts = os.path.join(tmpd, "artifacts")
            os.makedirs(artifacts)
            write_fasta_from_dict({"seq1": "MSEQ"}, os.path.join(artifacts, "sequences.faa"))
            output = os.path.join(tmpd, "hits.tsv")
            stderr = io.StringIO()

            with patch.object(script.subprocess, "run", side_effect=fake_run), patch("sys.stderr", stderr):
                with self.assertRaises(SystemExit) as cm:
                    script.run_hmmsearch(artifacts, "model.hmm", output)

            self.assertEqual(cm.exception.code, 1)
            self.assertIn("hmmsearch failed:", stderr.getvalue())
            self.assertIn("bad hmm", stderr.getvalue())
            self.assertIn("summary", stderr.getvalue())
