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
            self.assertTrue(os.path.exists(os.path.join(artifacts, "genomic_locus_with_leader.tsv")))
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

    def test_writes_unfiltered_fasta_with_leader_for_all_filtered_inputs(self):
        script = load_script(os.path.join(self.repo, "scripts", "filter-protein-by-rules.py"))
        self.write_manifest_and_sequences()

        with tempfile.TemporaryDirectory() as tmpd:
            self.write_rule_module(tmpd)
            fasta_output = os.path.join(tmpd, "res.faa")
            unfiltered_fasta_output = os.path.join(tmpd, "unfiltered.faa")
            stdin = io.StringIO("p_true\tg1\np_maybe\tg1\np_false\tg1\n")
            sys.path.insert(0, tmpd)
            try:
                with patch("sys.stdin", stdin):
                    script.main([
                        "-r", "constant_rules.mnsod_rule",
                        "--fasta-output", fasta_output,
                        "--unfiltered-fasta-output", unfiltered_fasta_output,
                        "--fasta-excludes-maybe",
                    ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("constant_rules", None)

            self.assertEqual(read_fasta_as_dict(fasta_output), {"p_true": "MT"})
            self.assertEqual(read_fasta_as_dict(unfiltered_fasta_output), {
                "p_true": "MT",
                "p_maybe": "MM",
                "p_false": "MF",
            })

    def test_writes_genomic_locus_with_leader_artifact(self):
        script = load_script(os.path.join(self.repo, "scripts", "filter-protein-by-rules.py"))
        ManifestTable.write_tsv(str(self.fx.area_genomics / "sequences.tsv"), [
            dict(
                sequence_accession="p_locus",
                sequence_database="g1",
                sequence_type="protein",
                sequence_source=SEQUENCE_SOURCE_NCBI,
            ),
        ])
        self.fx.write_ncbi_proteins("g1", {"p_locus": "MGP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "AAACCCGGGTTTAAACCCGGGTTTAAA"})
        self.fx.write_gff("g1", "\n".join([
            "ctg1\tsrc\tmRNA\t4\t24\t.\t+\t.\tID=tx1",
            "ctg1\tsrc\tCDS\t10\t15\t.\t+\t0\tID=cds1;Parent=tx1;protein_id=p_locus",
            "ctg1\tsrc\tCDS\t19\t21\t.\t+\t0\tID=cds2;Parent=tx1;protein_id=p_locus",
            "ctg1\tsrc\tstart_codon\t10\t12\t.\t+\t0\tParent=tx1",
            "ctg1\tsrc\tstop_codon\t22\t24\t.\t+\t0\tParent=tx1",
            "",
        ]))

        with tempfile.TemporaryDirectory() as tmpd:
            self.write_rule_module(tmpd)
            artifacts = os.path.join(tmpd, "artifacts")
            stdin = io.StringIO("p_locus\tg1\n")
            sys.path.insert(0, tmpd)
            try:
                with patch("sys.stdin", stdin):
                    script.main([
                        "-r", "constant_rules.mnsod_rule",
                        "--artifacts-dir", artifacts,
                    ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("constant_rules", None)

            with open(
                os.path.join(artifacts, "genomic_locus_with_leader.tsv"),
                "r",
                encoding="utf-8",
                newline="",
            ) as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(
                [
                    (row["feature type"], row["feature index"], row["feature position 1b"])
                    for row in rows
                ],
                [
                    ("start", "1", "7"),
                    ("stop", "1", "19"),
                    ("dss", "1", "12"),
                    ("ass", "1", "16"),
                ],
            )
            for row in rows:
                self.assertEqual(row["protein accession"], "p_locus")
                self.assertEqual(row["contig accession"], "ctg1")
                self.assertEqual(row["locus start 1b"], "4")
                self.assertEqual(row["locus end 1b"], "24")
                self.assertEqual(row["strand"], "+")
                self.assertEqual(row["error"], "")
