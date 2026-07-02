import csv
import io
import os
import tempfile
import unittest
from subprocess import CompletedProcess
from unittest.mock import patch

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tangle.detected import DetectedTable
from sieve.protein import (
    CuratedProtein,
    ProteinHMMAlignment,
    SEQUENCE_SOURCE_NCBI,
)
from sieve.rules import (
    HMMAlignment,
    KO,
    Leader,
    Pfam,
    RULE_ERROR,
    RULE_FALSE,
    RULE_MAYBE,
    RULE_TRUE,
    Rules,
    TFMotifs,
    _edge_distance,
    _parse_gimme_scan_output,
    _parse_targetp_output,
)
from tests.test_protein import DefaultsFixture


class RulesFixture(DefaultsFixture):

    def manifest_row(self, protein_accession, genome_accession, source=SEQUENCE_SOURCE_NCBI):
        return dict(
            sequence_accession=protein_accession,
            sequence_database=genome_accession,
            sequence_type="protein",
            sequence_source=source,
        )

    def detected_row(self, protein_accession, genome_accession, target_accession, target_database):
        return dict(
            detection_type="sequence",
            detection_method="hmm",
            batch="b1",
            query_accession=protein_accession,
            query_database=genome_accession,
            query_type="protein",
            target_accession=target_accession,
            target_database=target_database,
            target_type="protein",
            query_start=1,
            query_end=10,
            target_start=1,
            target_end=10,
        )

    def write_protein_fixture(self, protein_accession, genome_accession, sequence="MGP"):
        self.write_manifest([self.manifest_row(protein_accession, genome_accession)])
        self.write_ncbi_proteins(genome_accession, {protein_accession: sequence})

    def write_three_exon_gene(self, protein_accession, genome_accession, strand="+"):
        self.write_manifest([self.manifest_row(protein_accession, genome_accession)])
        self.write_ncbi_proteins(genome_accession, {protein_accession: "MGP"})
        self.write_genomic_fasta(genome_accession, {"ctg1": "A" * 90})
        if strand == "+":
            cds_rows = [
                "ctg1\tsrc\tCDS\t1\t10\t.\t+\t0\tID=cds1;Parent=tx1;protein_id=%s" % protein_accession,
                "ctg1\tsrc\tCDS\t31\t40\t.\t+\t0\tID=cds2;Parent=tx1;protein_id=%s" % protein_accession,
                "ctg1\tsrc\tCDS\t71\t80\t.\t+\t0\tID=cds3;Parent=tx1;protein_id=%s" % protein_accession,
            ]
            mrna = "ctg1\tsrc\tmRNA\t1\t90\t.\t+\t.\tID=tx1"
        else:
            cds_rows = [
                "ctg1\tsrc\tCDS\t71\t80\t.\t-\t0\tID=cds1;Parent=tx1;protein_id=%s" % protein_accession,
                "ctg1\tsrc\tCDS\t31\t40\t.\t-\t0\tID=cds2;Parent=tx1;protein_id=%s" % protein_accession,
                "ctg1\tsrc\tCDS\t1\t10\t.\t-\t0\tID=cds3;Parent=tx1;protein_id=%s" % protein_accession,
            ]
            mrna = "ctg1\tsrc\tmRNA\t1\t90\t.\t-\t.\tID=tx1"
        self.write_gff(genome_accession, "\n".join([mrna] + cds_rows + [""]))

    def write_single_exon_gene(self, protein_accession, genome_accession):
        self.write_manifest([self.manifest_row(protein_accession, genome_accession)])
        self.write_ncbi_proteins(genome_accession, {protein_accession: "MGP"})
        self.write_genomic_fasta(genome_accession, {"ctg1": "A" * 30})
        self.write_gff(genome_accession, "\n".join([
            "ctg1\tsrc\tmRNA\t1\t30\t.\t+\t.\tID=tx1",
            "ctg1\tsrc\tCDS\t1\t30\t.\t+\t0\tID=cds1;Parent=tx1;protein_id=%s" % protein_accession,
            "",
        ]))


class TestRules(unittest.TestCase):

    def setUp(self):
        CuratedProtein.clear_cache()
        self.fx = RulesFixture(self)

    def tearDown(self):
        CuratedProtein.clear_cache()
        self.fx.cleanup()

    def read_tsv(self, path):
        with open(path, "r", encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    def test_pfam_matches_prefix_before_version_and_ko_matches_exactly(self):
        self.fx.write_protein_fixture("p1", "g1")
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [
            self.fx.detected_row("p1", "g1", "PF02777.24", "Pfam"),
        ])
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_ko_assigned.tsv"), [
            self.fx.detected_row("p1", "g1", "K04564", "KO"),
        ])

        with tempfile.TemporaryDirectory() as tmpd:
            out = os.path.join(tmpd, "rules.tsv")
            rows = Rules(Pfam.matches("PF02777") & KO.matches("K04564")).check([("p1", "g1")], out)

            self.assertEqual(rows[0]["pass all"], RULE_TRUE)
            self.assertEqual(self.read_tsv(out)[0]["Pfam.matches('PF02777')"], RULE_TRUE)

    def test_rules_check_traces_each_atomic_rule_concisely(self):
        self.fx.write_protein_fixture("p1", "g1")
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [
            self.fx.detected_row("p1", "g1", "PF02777.24", "Pfam"),
        ])
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_ko_assigned.tsv"), [
            self.fx.detected_row("p1", "g1", "K00001", "KO"),
        ])
        err = io.StringIO()

        with patch("sys.stderr", err):
            with tempfile.TemporaryDirectory() as tmpd:
                Rules(Pfam.matches("PF02777") & KO.matches("K04564")).check(
                    [("p1", "g1")],
                    os.path.join(tmpd, "rules.tsv"),
                )

        self.assertEqual(err.getvalue().splitlines(), [
            "[rules 1/2] Pfam.matches('PF02777'): 1 proteins",
            "[rules 1/2] done: true=1",
            "[rules 2/2] KO.matches('K04564'): 1 proteins",
            "[rules 2/2] done: false=1",
        ])

    def test_rules_check_trace_can_be_disabled(self):
        self.fx.write_protein_fixture("p1", "g1")
        err = io.StringIO()

        with patch("sys.stderr", err):
            with tempfile.TemporaryDirectory() as tmpd:
                Rules(Pfam.matches("PF02777")).check(
                    [("p1", "g1")],
                    os.path.join(tmpd, "rules.tsv"),
                    trace=False,
                )

        self.assertEqual(err.getvalue(), "")

    def test_rules_check_continues_with_error_values(self):
        self.fx.write_protein_fixture("p1", "g1")
        rule = Pfam.matches("PF00001") & Leader.is_mTP()

        with patch("sieve.rules.subprocess.run", side_effect=RuntimeError("no docker")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(rule).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["pass all"], RULE_FALSE)
        self.assertEqual(rows[0]["Leader.is_mTP()"], RULE_ERROR)

    def test_rules_check_writes_artifacts_for_batch_tools(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        targetp_output = "\n".join([
            "# TargetP-2.0",
            "seq0\tmTP\t0.1\t0.2\t0.7\tCS pos: 1-2. AA-AA. Pr: 0.1",
            "",
        ])
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "seq0 12 17 GM.5.0.Rel.0001 8.0 +",
            "seq0 25 30 GM.5.0.bZIP.0001 8.0 -",
            "",
        ])

        def fake_run(cmd, check, capture_output, text):
            if cmd[0] == "docker":
                return CompletedProcess(cmd, 0, stdout=targetp_output, stderr="targetp stderr")
            if cmd[0] == "gimme":
                return CompletedProcess(cmd, 0, stdout=gimme_output, stderr="gimme stderr")
            raise AssertionError(cmd)

        with patch("sieve.rules.subprocess.run", side_effect=fake_run):
            with tempfile.TemporaryDirectory() as tmpd:
                out = os.path.join(tmpd, "rules.tsv")
                artifacts = os.path.join(tmpd, "artifacts")
                rows = Rules(
                    Leader.is_mTP()
                    & TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron()
                ).check([("p1", "g1")], out, artifacts_dir=artifacts)

                self.assertEqual(rows[0]["pass all"], RULE_TRUE)
                artifact_dirs = os.listdir(artifacts)
                leader_dir = os.path.join(artifacts, next(d for d in artifact_dirs if d.startswith("Leader.is_mTP")))
                tf_dir = os.path.join(artifacts, next(d for d in artifact_dirs if d.startswith("TFMotifs.has_within")))

                self.assertTrue(os.path.exists(os.path.join(leader_dir, "query.faa")))
                self.assertTrue(os.path.exists(os.path.join(leader_dir, "command.txt")))
                with open(os.path.join(leader_dir, "stdout.txt"), "r", encoding="utf-8") as f:
                    self.assertEqual(f.read(), targetp_output)
                with open(os.path.join(leader_dir, "stderr.txt"), "r", encoding="utf-8") as f:
                    self.assertEqual(f.read(), "targetp stderr")

                self.assertTrue(os.path.exists(os.path.join(tf_dir, "locus.fna")))
                self.assertTrue(os.path.exists(os.path.join(tf_dir, "command.txt")))
                with open(os.path.join(tf_dir, "stdout.txt"), "r", encoding="utf-8") as f:
                    self.assertEqual(f.read(), gimme_output)
                with open(os.path.join(tf_dir, "stderr.txt"), "r", encoding="utf-8") as f:
                    self.assertEqual(f.read(), "gimme stderr")

    def test_rules_check_uses_absolute_artifact_paths_for_docker_mounts(self):
        self.fx.write_protein_fixture("p1", "g1")

        def fake_run(cmd, check, capture_output, text):
            mount_arg = cmd[cmd.index("-v") + 1]
            mounted_path = mount_arg.split(":", 1)[0]
            self.assertTrue(os.path.isabs(mounted_path))
            return CompletedProcess(cmd, 0, stdout="seq0\tmTP\t0.1\t0.2\t0.7\t\n", stderr="")

        with patch("sieve.rules.subprocess.run", side_effect=fake_run):
            with tempfile.TemporaryDirectory() as tmpd:
                old_cwd = os.getcwd()
                try:
                    os.chdir(tmpd)
                    rows = Rules(Leader.is_mTP()).check(
                        [("p1", "g1")],
                        "rules.tsv",
                        artifacts_dir="tmp",
                    )
                finally:
                    os.chdir(old_cwd)

        self.assertEqual(rows[0]["pass all"], RULE_TRUE)

    def test_or_and_and_composite_results_use_true_false_maybe_error(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "seq0 45 50 GM.5.0.Rel.0001 8.0 +",
            "",
        ])
        rule = (
            TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron(2)
            | Leader.is_mTP()
        ) & Pfam.matches("PF00001")

        def fake_run(cmd, check, capture_output, text):
            if cmd[0] == "gimme":
                return CompletedProcess(cmd, 0, stdout=gimme_output, stderr="")
            raise RuntimeError("targetp failure")

        with patch("sieve.rules.subprocess.run", side_effect=fake_run):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(rule).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["pass all"], RULE_FALSE)
        self.assertEqual(rows[0]["Leader.is_mTP()"], RULE_ERROR)
        self.assertEqual(
            rows[0]["TFMotifs.has_within(20, 'GM.5.0.Rel', 'GM.5.0.bZIP', min_score_threshold=8).in_intron(2)"],
            RULE_MAYBE,
        )

    def test_hmm_position_motif_and_coverage_rules_are_strict_about_gaps(self):
        self.fx.write_protein_fixture("p1", "g1")
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("ACD-E"), id="p1"),
        ])
        alignment.column_annotations["reference_annotation"] = "xxxxx"
        protein_alignment = ProteinHMMAlignment(alignment, "p1")

        with patch.object(CuratedProtein, "hmm_align", return_value=protein_alignment):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    HMMAlignment("/models/profile.hmm").is_at("ACD", 1)
                    & HMMAlignment("/models/profile.hmm").is_at("DE", 3)
                    & HMMAlignment("/models/profile.hmm").covers(1, 3)
                    & HMMAlignment("/models/profile.hmm").covers(1, 4)
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["HMMAlignment('profile.hmm').is_at('ACD', 1)"], RULE_TRUE)
        self.assertEqual(rows[0]["HMMAlignment('profile.hmm').is_at('DE', 3)"], RULE_FALSE)
        self.assertEqual(rows[0]["HMMAlignment('profile.hmm').covers(1, 3)"], RULE_TRUE)
        self.assertEqual(rows[0]["HMMAlignment('profile.hmm').covers(1, 4)"], RULE_FALSE)

    def test_leader_rules_batch_targetp_and_use_sequence_with_leader(self):
        self.fx.write_manifest([
            self.fx.manifest_row("p1", "g1"),
            self.fx.manifest_row("p2", "g2"),
        ])
        self.fx.write_ncbi_proteins("g1", {"p1": "MA"})
        self.fx.write_ncbi_proteins("g2", {"p2": "MG"})

        def fake_run(cmd, check, capture_output, text):
            self.assertEqual(cmd[0:5], ["docker", "run", "--rm", "--platform", "linux/amd64"])
            fasta_path = cmd[cmd.index("-v") + 1].split(":", 1)[0] + "/query.faa"
            with open(fasta_path, "r", encoding="utf-8") as f:
                fasta_text = f.read()
            self.assertIn(">seq0\nMA\n", fasta_text)
            self.assertIn(">seq1\nMG\n", fasta_text)
            return CompletedProcess(cmd, 0, stdout="\n".join([
                "# TargetP-2.0",
                "# ID\tPrediction\tnoTP\tSP\tmTP\tCS Position",
                "seq0\tmTP\t0.1\t0.2\t0.7\tCS pos: 1-2. AA-AA. Pr: 0.1",
                "seq1\tSP\t0.1\t0.8\t0.1\tCS pos: 1-2. AA-AA. Pr: 0.1",
                "",
            ]), stderr="")

        with patch("sieve.rules.subprocess.run", side_effect=fake_run):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(Leader.is_mTP()).check(
                    [("p1", "g1"), ("p2", "g2")],
                    os.path.join(tmpd, "rules.tsv"),
                )

        self.assertEqual([row["Leader.is_mTP()"] for row in rows], [RULE_TRUE, RULE_FALSE])

    def test_tf_motifs_pass_for_hits_within_intron_on_forward_gene_even_opposite_hit_strands(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "seq0 45 50 GM.5.0.Rel.0001 8.0 +",
            "seq0 70 75 GM.5.0.bZIP.0001 9.0 -",
            "seq0 55 60 GM.5.0.bZIP.0002 9.0 -",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron(2)
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["pass all"], RULE_TRUE)

    def test_tf_motifs_use_gene_direction_for_reverse_gene_intron_order(self):
        self.fx.write_three_exon_gene("p1", "g1", "-")
        locus = CuratedProtein("p1", "g1").genomic_locus_with_leader()
        self.assertEqual(locus.strand, -1)
        self.assertEqual(locus.cds_intervals_1b, [(11, 20), (51, 60), (81, 90)])

        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "seq0 62 66 GM.5.0.Rel.0001 8.0 -",
            "seq0 75 80 GM.5.0.bZIP.0001 8.0 +",
            "seq0 25 30 GM.5.0.Rel.0002 8.0 +",
            "seq0 35 39 GM.5.0.bZIP.0002 8.0 -",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron(2)
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["pass all"], RULE_TRUE)

    def test_tf_motifs_in_intron_without_number_matches_any_intron(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "seq0 12 17 GM.5.0.Rel.0001 8.0 +",
            "seq0 25 30 GM.5.0.bZIP.0001 8.0 -",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron()
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(
            rows[0]["TFMotifs.has_within(20, 'GM.5.0.Rel', 'GM.5.0.bZIP', min_score_threshold=8).in_intron()"],
            RULE_TRUE,
        )

    def test_tf_motifs_in_intron_without_number_is_false_for_single_exon_gene(self):
        self.fx.write_single_exon_gene("p1", "g1")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "seq0 2 6 GM.5.0.Rel.0001 8.0 +",
            "seq0 8 12 GM.5.0.bZIP.0001 8.0 -",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron()
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["pass all"], RULE_FALSE)

    def test_tf_motifs_return_maybe_when_no_pair_above_default_threshold(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "seq0 45 50 GM.5.0.Rel.0001 7.9 +",
            "seq0 55 60 GM.5.0.bZIP.0001 9.0 +",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron(2)
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["pass all"], RULE_MAYBE)

    def test_tf_motifs_return_maybe_when_hits_are_not_within_nearest_edge_distance(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "seq0 41 45 GM.5.0.Rel.0001 8.0 +",
            "seq0 67 70 GM.5.0.bZIP.0001 9.0 +",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron(2)
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["pass all"], RULE_MAYBE)

    def test_tf_motifs_without_scope_can_match_in_any_intron_or_exon(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "seq0 12 17 GM.5.0.Rel.0001 8.0 +",
            "seq0 25 30 GM.5.0.bZIP.0001 8.0 -",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP")
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(
            rows[0]["TFMotifs.has_within(20, 'GM.5.0.Rel', 'GM.5.0.bZIP', min_score_threshold=8)"],
            RULE_TRUE,
        )

    def test_tf_motifs_without_scope_ignores_hits_outside_exons_and_introns(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "seq0 82 85 GM.5.0.Rel.0001 8.0 +",
            "seq0 86 90 GM.5.0.bZIP.0001 8.0 -",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP")
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["pass all"], RULE_MAYBE)

    def test_tf_motifs_in_exon_matches_any_exon(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "seq0 32 35 GM.5.0.Rel.0001 8.0 +",
            "seq0 37 40 GM.5.0.bZIP.0001 8.0 -",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_exon()
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(
            rows[0]["TFMotifs.has_within(20, 'GM.5.0.Rel', 'GM.5.0.bZIP', min_score_threshold=8).in_exon()"],
            RULE_TRUE,
        )

    def test_tf_motifs_in_numbered_exon_uses_gene_direction(self):
        self.fx.write_three_exon_gene("p1", "g1", "-")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "seq0 52 55 GM.5.0.Rel.0001 8.0 -",
            "seq0 57 60 GM.5.0.bZIP.0001 8.0 +",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_exon(2)
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(
            rows[0]["TFMotifs.has_within(20, 'GM.5.0.Rel', 'GM.5.0.bZIP', min_score_threshold=8).in_exon(2)"],
            RULE_TRUE,
        )

    def test_tf_motifs_between_uses_locus_coordinates(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "seq0 82 85 GM.5.0.Rel.0001 8.0 +",
            "seq0 86 90 GM.5.0.bZIP.0001 8.0 -",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").between(90, 81)
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(
            rows[0]["TFMotifs.has_within(20, 'GM.5.0.Rel', 'GM.5.0.bZIP', min_score_threshold=8).between(90, 81)"],
            RULE_TRUE,
        )

    def test_tf_motifs_scope_methods_can_only_be_called_once(self):
        base = TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP")

        scoped_rules = [
            base.in_exon(),
            base.in_exon(2),
            base.in_intron(),
            base.in_intron(2),
            base.between(10, 20),
        ]

        for scoped_rule in scoped_rules:
            with self.assertRaises(ValueError):
                scoped_rule.in_exon()
            with self.assertRaises(ValueError):
                scoped_rule.in_intron()
            with self.assertRaises(ValueError):
                scoped_rule.between(10, 20)


class TestRuleParsingHelpers(unittest.TestCase):

    def test_parse_targetp_output(self):
        parsed = _parse_targetp_output("\n".join([
            "# TargetP-2.0",
            "# ID Prediction noTP SP mTP CS Position",
            "seq0\tmTP\t0.1\t0.2\t0.7\tCS pos: 1-2. AA-AA. Pr: 0.1",
            "seq1 noTP 0.9 0.1 0.0",
        ]))

        self.assertEqual(parsed, {"seq0": "mTP", "seq1": "noTP"})

    def test_parse_gimme_scan_output(self):
        parsed = _parse_gimme_scan_output("\n".join([
            "sequence start end feature score strand",
            "seq0 45 50 GM.5.0.Rel.0001 8.0 +",
            "seq0 55 60 GM.5.0.bZIP.0001 9.0 -",
        ]))

        self.assertEqual(len(parsed["seq0"]), 2)
        self.assertEqual(parsed["seq0"][1].feature, "GM.5.0.bZIP.0001")
        self.assertEqual(parsed["seq0"][1].strand, "-")

    def test_edge_distance_uses_nearest_edges_and_overlap(self):
        self.assertEqual(_edge_distance(10, 20, 25, 30), 5)
        self.assertEqual(_edge_distance(25, 30, 10, 20), 5)
        self.assertEqual(_edge_distance(10, 20, 20, 30), 0)
        self.assertEqual(_edge_distance(10, 20, 15, 25), 0)


if __name__ == "__main__":
    unittest.main()
