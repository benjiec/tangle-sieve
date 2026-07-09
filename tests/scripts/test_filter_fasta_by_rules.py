import csv
import os
import sys
import tempfile
import unittest
from subprocess import CompletedProcess
from unittest.mock import patch

from sieve.artifacts import rule_results_tsv, sequences_fasta
from tangle.sequence import read_fasta_as_dict, write_fasta_from_dict

from tests.scripts.helpers import load_script


def domtblout_line(sequence, model_name, model_accession, full_score, domain_score, ali_from=7):
    return " ".join([
        sequence,
        "-",
        "100",
        model_name,
        model_accession,
        "200",
        "1e-20",
        str(full_score),
        "0.0",
        "1",
        "1",
        "1e-10",
        "1e-10",
        str(domain_score),
        "0.0",
        "3",
        "50",
        str(ali_from),
        str(ali_from + 47),
        str(ali_from - 1),
        str(ali_from + 48),
        "0.98",
    ])


class TestFilterFastaByRulesScript(unittest.TestCase):

    def setUp(self):
        self.repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.script = load_script(os.path.join(self.repo, "scripts", "filter-fasta-by-rules.py"))

    def write_lines(self, path, lines):
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    def write_rule_module(self, tmpd, lines):
        path = os.path.join(tmpd, "fasta_rules.py")
        self.write_lines(path, lines)
        return "fasta_rules.mnsod_rule"

    def read_rows(self, path):
        with open(path, "r", encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    def test_runs_hmmsearches_and_writes_standard_artifacts(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "proteins.faa")
            pfam_hmm = os.path.join(tmpd, "pfam.hmm")
            ko_hmm = os.path.join(tmpd, "ko.hmm")
            thresholds = os.path.join(tmpd, "thresholds.tsv")
            artifacts = os.path.join(tmpd, "artifacts")
            self.write_lines(pfam_hmm, ["PFAM"])
            self.write_lines(ko_hmm, ["KO"])
            write_fasta_from_dict({
                "p1": "MAAAAA",
                "p2": "MBBBBB",
                "p3": "MCCCCC",
            }, fasta)
            self.write_lines(thresholds, [
                "model threshold score_type definition",
                "K04564 100 full manganese superoxide dismutase",
            ])
            rule_spec = self.write_rule_module(tmpd, [
                "from sieve.rules import KO, Pfam, Rules",
                "mnsod_rule = Rules(Pfam.matches('PF00081') & KO.matches('K04564'))",
                "",
            ])

            def fake_run(cmd, check, capture_output, text):
                domtblout = cmd[cmd.index("--domtblout") + 1]
                if "--cut_ga" in cmd:
                    self.write_lines(domtblout, [
                        domtblout_line("p1", "SOD_Fe_N", "PF00081.28", 200, 180),
                        domtblout_line("p2", "SOD_Fe_N", "PF00081.28", 200, 180),
                    ])
                else:
                    self.write_lines(domtblout, [
                        domtblout_line("p1", "K04564", "-", 150, 80),
                        domtblout_line("p2", "K04564", "-", 50, 200),
                    ])
                return CompletedProcess(cmd, 0, stdout="", stderr="")

            sys.path.insert(0, tmpd)
            try:
                with patch("subprocess.run", side_effect=fake_run) as run:
                    self.script.main([
                        "--fasta", fasta,
                        "-r", rule_spec,
                        "--artifacts-dir", artifacts,
                        "--pfam-hmm", pfam_hmm,
                        "--ko-hmm", ko_hmm,
                        "--ko-thresholds", thresholds,
                    ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("fasta_rules", None)

            calls = [call.args[0] for call in run.call_args_list]
            self.assertIn("--cut_ga", calls[0])
            self.assertNotIn("--cut_ga", calls[1])
            self.assertTrue(os.path.exists(os.path.join(artifacts, "pfam.domtblout")))
            self.assertTrue(os.path.exists(os.path.join(artifacts, "ko.domtblout")))
            rows = self.read_rows(rule_results_tsv(artifacts))
            self.assertEqual(
                [(row["protein accession"], row["pass all"]) for row in rows],
                [("p1", "true"), ("p2", "false"), ("p3", "false")],
            )
            self.assertEqual(read_fasta_as_dict(sequences_fasta(artifacts)), {
                "p1": "MAAAAA",
                "p2": "MBBBBB",
                "p3": "MCCCCC",
            })

    def test_ko_hmm_requires_threshold_file(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "proteins.faa")
            ko_hmm = os.path.join(tmpd, "ko.hmm")
            artifacts = os.path.join(tmpd, "artifacts")
            write_fasta_from_dict({"p1": "MAAAAA"}, fasta)
            self.write_lines(ko_hmm, ["KO"])
            rule_spec = self.write_rule_module(tmpd, [
                "from sieve.rules import KO, Rules",
                "mnsod_rule = Rules(KO.matches('K04564'))",
                "",
            ])

            sys.path.insert(0, tmpd)
            try:
                with self.assertRaisesRegex(ValueError, "--ko-thresholds is required"):
                    self.script.main([
                        "--fasta", fasta,
                        "-r", rule_spec,
                        "--artifacts-dir", artifacts,
                        "--ko-hmm", ko_hmm,
                    ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("fasta_rules", None)

    def test_genomic_only_rules_become_rule_errors_for_fasta_proteins(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "proteins.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            write_fasta_from_dict({"p1": "MAAAAA"}, fasta)
            rule_spec = self.write_rule_module(tmpd, [
                "from sieve.rules import Rule, Rules, RULE_TRUE",
                "class GenomicRule(Rule):",
                "    label = 'GenomicRule'",
                "    def evaluate(self, context):",
                "        context.protein.genomic_locus()",
                "        return RULE_TRUE",
                "mnsod_rule = Rules(GenomicRule())",
                "",
            ])

            sys.path.insert(0, tmpd)
            try:
                self.script.main([
                    "--fasta", fasta,
                    "-r", rule_spec,
                    "--artifacts-dir", artifacts,
                ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("fasta_rules", None)

            rows = self.read_rows(rule_results_tsv(artifacts))
            self.assertEqual(rows[0]["GenomicRule"], "error")
            self.assertEqual(rows[0]["pass all"], "error")
