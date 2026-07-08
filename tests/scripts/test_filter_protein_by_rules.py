import csv
import io
import os
import sys
import tempfile
import unittest
from unittest.mock import patch

from tangle.detected import DetectedTable
from tangle.manifest import ManifestTable
from tangle.sequence import read_fasta_as_dict

from sieve.protein import CuratedProtein, SEQUENCE_SOURCE_HMM_DETECTED, SEQUENCE_SOURCE_NCBI
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

    def detected_protein_row(self, protein_accession, genome_accession, contig_accession, q_start, q_end, t_start, t_end):
        return dict(
            detection_type="model",
            detection_method="hmm",
            batch="b1",
            query_accession=contig_accession,
            query_database=genome_accession,
            query_type="contig",
            target_accession=protein_accession,
            target_database=genome_accession,
            target_type="protein",
            target_model="HMM1",
            query_start=q_start,
            query_end=q_end,
            target_start=t_start,
            target_end=t_end,
            evalue=0.001,
            bitscore=10,
        )

    def write_rule_module(self, tmpd):
        module_path = os.path.join(tmpd, "constant_rules.py")
        with open(module_path, "w", encoding="utf-8") as f:
            f.write("\n".join([
                "from sieve.rules import Rules",
                "from tests.scripts.helpers import ConstantByProteinRule",
                "mnsod_rule = Rules(ConstantByProteinRule())",
                "",
            ]))

    def write_leader_rule_module(self, tmpd):
        module_path = os.path.join(tmpd, "leader_rules.py")
        with open(module_path, "w", encoding="utf-8") as f:
            f.write("\n".join([
                "from sieve.rules import Leader, Rules",
                "mnsod_rule = Rules(Leader().is_mTP())",
                "",
            ]))

    def write_leader_and_rule_module(self, tmpd):
        module_path = os.path.join(tmpd, "leader_and_rules.py")
        with open(module_path, "w", encoding="utf-8") as f:
            f.write("\n".join([
                "from sieve.rules import Leader, Rules",
                "from tests.scripts.helpers import ConstantByProteinRule",
                "mnsod_rule = Rules(Leader().is_mTP() & ConstantByProteinRule())",
                "",
            ]))

    def write_leader_or_rule_module(self, tmpd):
        module_path = os.path.join(tmpd, "leader_or_rules.py")
        with open(module_path, "w", encoding="utf-8") as f:
            f.write("\n".join([
                "from sieve.rules import Leader, Rules",
                "mnsod_rule = Rules(Leader().is_mTP() | Leader().is_SP())",
                "",
            ]))

    def write_pfam_anchor_leader_rule_module(self, tmpd):
        module_path = os.path.join(tmpd, "pfam_anchor_leader_rules.py")
        with open(module_path, "w", encoding="utf-8") as f:
            f.write("\n".join([
                "from sieve.rules import Leader, Rules",
                "mnsod_rule = Rules(Leader().upstreamOfPfam('PF00081').betweenAA(-5, 0).is_mTP())",
                "",
            ]))

    def test_writes_rule_tsv_and_sequences_fasta(self):
        script = load_script(os.path.join(self.repo, "scripts", "filter-protein-by-rules.py"))
        self.write_manifest_and_sequences()

        with tempfile.TemporaryDirectory() as tmpd:
            self.write_rule_module(tmpd)
            artifacts = os.path.join(tmpd, "artifacts")
            stdin = io.StringIO("p_true\tg1\np_maybe\tg1\np_false\tg1\np_true\tg1\n")
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

            self.assertTrue(os.path.exists(os.path.join(artifacts, "rule-results.tsv")))
            self.assertTrue(os.path.exists(os.path.join(artifacts, "genomes", "g1")))
            self.assertTrue(os.path.exists(os.path.join(artifacts, "genomic_locus_with_leader.tsv")))
            with open(os.path.join(artifacts, "rule-results.tsv"), "r", encoding="utf-8", newline="") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(
                [(row["protein accession"], row["sequence accession"]) for row in rows],
                [("p_true", "p_true"), ("p_maybe", "p_maybe"), ("p_false", "p_false")],
            )
            self.assertEqual(read_fasta_as_dict(os.path.join(artifacts, "sequences.faa")), {
                "p_true": "MT",
                "p_maybe": "MM",
                "p_false": "MF",
            })

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

    def test_preserves_rule_annotation_columns_in_recombined_results(self):
        script = load_script(os.path.join(self.repo, "scripts", "filter-protein-by-rules.py"))
        self.write_manifest_and_sequences()

        with tempfile.TemporaryDirectory() as tmpd:
            module_path = os.path.join(tmpd, "annotated_rules.py")
            with open(module_path, "w", encoding="utf-8") as f:
                f.write("\n".join([
                    "from sieve.rules import Rules",
                    "from tests.scripts.helpers import AnnotatedByProteinRule",
                    "mnsod_rule = Rules(AnnotatedByProteinRule())",
                    "",
                ]))
            artifacts = os.path.join(tmpd, "artifacts")
            stdin = io.StringIO("p_true\tg1\np_false\tg1\n")
            sys.path.insert(0, tmpd)
            try:
                with patch("sys.stdin", stdin):
                    script.main([
                        "-r", "annotated_rules.mnsod_rule",
                        "--artifacts-dir", artifacts,
                    ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("annotated_rules", None)

            with open(os.path.join(artifacts, "rule-results.tsv"), "r", encoding="utf-8", newline="") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(
                [(row["protein accession"], row["sequence accession"], row["contig accession"], row["Example.call"]) for row in rows],
                [("p_true", "p_true", "", "p_true_call"), ("p_false", "p_false", "", "p_false_call")],
            )

    def test_sequences_fasta_contains_all_filtered_inputs(self):
        script = load_script(os.path.join(self.repo, "scripts", "filter-protein-by-rules.py"))
        self.write_manifest_and_sequences()

        with tempfile.TemporaryDirectory() as tmpd:
            self.write_rule_module(tmpd)
            artifacts = os.path.join(tmpd, "artifacts")
            stdin = io.StringIO("p_true\tg1\np_maybe\tg1\np_false\tg1\n")
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

            self.assertEqual(read_fasta_as_dict(os.path.join(artifacts, "sequences.faa")), {
                "p_true": "MT",
                "p_maybe": "MM",
                "p_false": "MF",
            })

    def test_unfiltered_fasta_uses_rule_scoped_leader_accessions(self):
        script = load_script(os.path.join(self.repo, "scripts", "filter-protein-by-rules.py"))
        self.fx.write_manifest([
            dict(
                sequence_accession="p1",
                sequence_database="g1",
                sequence_type="protein",
                sequence_source=SEQUENCE_SOURCE_HMM_DETECTED,
            ),
        ])
        self.fx.write_detected_proteins("g1", {"p1": "AAAA"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "GCTGCTATGAAAGCTGCT"})
        self.fx.write_detected_rows("g1", [
            self.detected_protein_row("p1", "g1", "ctg1", 1, 6, 37, 38),
            self.detected_protein_row("p1", "g1", "ctg1", 13, 18, 39, 40),
        ])
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [
            dict(
                detection_type="sequence",
                detection_method="hmm",
                batch="b1",
                query_accession="p1",
                query_database="g1",
                query_type="protein",
                target_accession="PF00081.28",
                target_database="Pfam",
                target_type="protein",
                query_start=4,
                query_end=12,
                target_start=1,
                target_end=9,
                evalue=0.001,
                bitscore=10,
            ),
        ])

        def fake_run(cmd, check, capture_output, text):
            fasta_path = cmd[cmd.index("-v") + 1].split(":", 1)[0] + "/query.faa"
            ids = list(read_fasta_as_dict(fasta_path).keys())
            self.assertEqual(ids, ["p1_with_leader_u3_PF00081_anchor_M"])
            stdout = "# TargetP-2.0\np1_with_leader_u3_PF00081_anchor_M\tmTP\t0.1\t0.1\t0.8\t\n"
            return type("Completed", (), {"returncode": 0, "stdout": stdout, "stderr": ""})()

        with tempfile.TemporaryDirectory() as tmpd:
            self.write_pfam_anchor_leader_rule_module(tmpd)
            artifacts = os.path.join(tmpd, "artifacts")
            stdin = io.StringIO("p1\tg1\n")
            sys.path.insert(0, tmpd)
            try:
                with patch("sys.stdin", stdin), patch("sieve.rules.subprocess.run", side_effect=fake_run):
                    script.main([
                        "-r", "pfam_anchor_leader_rules.mnsod_rule",
                        "--artifacts-dir", artifacts,
                    ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("pfam_anchor_leader_rules", None)

            self.assertEqual(
                set(read_fasta_as_dict(os.path.join(artifacts, "sequences.faa")).keys()),
                {"p1_with_leader_u3_PF00081_anchor_M"},
            )

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
                reader = csv.DictReader(f, delimiter="\t")
                rows = list(reader)
            self.assertEqual(reader.fieldnames[:4], [
                "protein accession",
                "genome accession",
                "contig accession",
                "sequence accession",
            ])
            self.assertEqual(
                [
                    (row["feature type"], row["feature index"], row["feature position 1b"])
                    for row in rows
                ],
                [
                    ("start", "", "7"),
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

    def test_writes_hmm_detected_leader_candidate_starts_to_locus_artifact(self):
        script = load_script(os.path.join(self.repo, "scripts", "filter-protein-by-rules.py"))
        ManifestTable.write_tsv(str(self.fx.area_genomics / "sequences.tsv"), [
            dict(
                sequence_accession="p_locus",
                sequence_database="g1",
                sequence_type="protein",
                sequence_source=SEQUENCE_SOURCE_HMM_DETECTED,
            ),
        ])
        self.fx.write_detected_proteins("g1", {"p_locus": "KMMP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "NNNATGATGAAAATGATGCCC"})
        self.fx.write_detected_rows("g1", [
            self.detected_protein_row("p_locus", "g1", "ctg1", 10, 21, 1, 4),
        ])

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

        start_rows = [row for row in rows if row["feature type"] == "start"]
        self.assertEqual(
            [
                (
                    row["feature index"],
                    row["sequence accession"],
                    row["feature position 1b"],
                )
                for row in start_rows
            ],
            [
                ("", "p_locus_with_leader_u2_M", "1"),
                ("", "p_locus_with_leader_u1_M", "4"),
                ("", "p_locus_with_leader_2_M", "10"),
                ("", "p_locus_with_leader_3_M", "13"),
            ],
        )
        for row in rows:
            self.assertEqual(row["contig accession"], "ctg1")
            self.assertEqual(row["locus start 1b"], "4")
            self.assertEqual(row["locus end 1b"], "21")

    def test_hmm_detected_locus_artifact_has_no_start_for_original_sequence_fallback(self):
        script = load_script(os.path.join(self.repo, "scripts", "filter-protein-by-rules.py"))
        ManifestTable.write_tsv(str(self.fx.area_genomics / "sequences.tsv"), [
            dict(
                sequence_accession="p_locus",
                sequence_database="g1",
                sequence_type="protein",
                sequence_source=SEQUENCE_SOURCE_HMM_DETECTED,
            ),
        ])
        self.fx.write_detected_proteins("g1", {"p_locus": "RYDA"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "CGTTATGATGCT"})
        self.fx.write_detected_rows("g1", [
            self.detected_protein_row("p_locus", "g1", "ctg1", 1, 12, 1, 4),
        ])

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

        self.assertEqual([row["feature type"] for row in rows], [""])
        self.assertEqual(rows[0]["protein accession"], "p_locus")
        self.assertEqual(rows[0]["sequence accession"], "p_locus")
