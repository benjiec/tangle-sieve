import csv
import io
import os
import sys
import tempfile
import unittest
from unittest.mock import patch

from tangle.manifest import ManifestTable
from tangle.sequence import read_fasta_as_dict

from sieve.protein import CuratedProtein, SEQUENCE_SOURCE_NCBI
from tests.fixtures import DefaultsFixture
from tests.scripts.helpers import load_script


class TestFilterProteinByRulesScript(unittest.TestCase):

    def setUp(self):
        CuratedProtein.clear_cache()
        self.fx = DefaultsFixture(self)
        self.repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    def tearDown(self):
        CuratedProtein.clear_cache()
        self.fx.cleanup()

    def write_manifest_and_sequences(self):
        ManifestTable.write_tsv(str(self.fx.area_genomics / "sequences.tsv"), [
            dict(
                sequence_accession="p_true",
                sequence_database="g1",
                sequence_type="protein",
                sequence_source=SEQUENCE_SOURCE_NCBI,
            ),
            dict(
                sequence_accession="p_maybe",
                sequence_database="g1",
                sequence_type="protein",
                sequence_source=SEQUENCE_SOURCE_NCBI,
            ),
            dict(
                sequence_accession="p_false",
                sequence_database="g1",
                sequence_type="protein",
                sequence_source=SEQUENCE_SOURCE_NCBI,
            ),
        ])
        self.fx.write_ncbi_proteins("g1", {
            "p_true": "MT",
            "p_maybe": "MM",
            "p_false": "MF",
        })

    def write_rule_module(self, tmpd):
        module_path = os.path.join(tmpd, "constant_rules.py")
        with open(module_path, "w", encoding="utf-8") as f:
            f.write("\n".join([
                "from sieve.rules import Rules",
                "from tests.scripts.helpers import ConstantByProteinRule",
                "mnsod_rule = Rules(ConstantByProteinRule())",
                "",
            ]))

    def test_writes_rule_tsv_and_fasta(self):
        script = load_script(os.path.join(self.repo, "scripts", "filter-protein-by-rules.py"))
        self.write_manifest_and_sequences()

        with tempfile.TemporaryDirectory() as tmpd:
            self.write_rule_module(tmpd)
            artifacts = os.path.join(tmpd, "artifacts")
            fasta_output = os.path.join(tmpd, "res.faa")
            stdin = io.StringIO("p_true\tg1\np_maybe\tg1\np_false\tg1\np_true\tg1\n")
            sys.path.insert(0, tmpd)
            try:
                with patch("sys.stdin", stdin):
                    script.main([
                        "-r", "constant_rules.mnsod_rule",
                        "--artifacts-dir", artifacts,
                        "--fasta-output", fasta_output,
                    ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("constant_rules", None)

            self.assertTrue(os.path.exists(os.path.join(artifacts, "rule-results.tsv")))
            self.assertTrue(os.path.exists(os.path.join(artifacts, "genomes", "g1")))
            self.assertEqual(read_fasta_as_dict(fasta_output), {
                "p_true": "MT",
                "p_maybe": "MM",
            })

    def test_can_exclude_maybe_from_fasta(self):
        script = load_script(os.path.join(self.repo, "scripts", "filter-protein-by-rules.py"))
        self.write_manifest_and_sequences()

        with tempfile.TemporaryDirectory() as tmpd:
            self.write_rule_module(tmpd)
            fasta_output = os.path.join(tmpd, "res.faa")
            stdin = io.StringIO("p_true\tg1\np_maybe\tg1\n")
            sys.path.insert(0, tmpd)
            try:
                with patch("sys.stdin", stdin):
                    script.main([
                        "-r", "constant_rules.mnsod_rule",
                        "--fasta-output", fasta_output,
                        "--fasta-excludes-maybe",
                    ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("constant_rules", None)

            self.assertEqual(read_fasta_as_dict(fasta_output), {"p_true": "MT"})

    def test_ignores_input_proteins_missing_from_manifest(self):
        script = load_script(os.path.join(self.repo, "scripts", "filter-protein-by-rules.py"))
        self.write_manifest_and_sequences()

        with tempfile.TemporaryDirectory() as tmpd:
            self.write_rule_module(tmpd)
            artifacts = os.path.join(tmpd, "artifacts")
            stdin = io.StringIO("p_true\tg1\np_missing\tg1\n")
            stderr = io.StringIO()
            sys.path.insert(0, tmpd)
            try:
                with patch("sys.stdin", stdin), patch("sys.stderr", stderr):
                    script.main([
                        "-r", "constant_rules.mnsod_rule",
                        "--artifacts-dir", artifacts,
                    ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("constant_rules", None)

            self.assertIn("Ignoring p_missing\tg1", stderr.getvalue())
            with open(os.path.join(artifacts, "rule-results.tsv"), "r", encoding="utf-8", newline="") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(
                [(row["protein accession"], row["genome accession"]) for row in rows],
                [("p_true", "g1")],
            )
