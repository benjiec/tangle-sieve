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
from sieve.fasta_protein import FastaProtein
from sieve.protein import (
    CuratedProtein,
    ProteinHMMAlignment,
    SEQUENCE_SOURCE_HMM_DETECTED,
    SEQUENCE_SOURCE_NCBI,
    hmm_align_sequences,
)
from sieve.rules import (
    HMMAlignment,
    KO,
    Leader,
    Pfam,
    RULE_ERROR,
    RULE_FALSE,
    RULE_MAYBE,
    RULE_TOO_FAR,
    RULE_TRUE,
    RULE_YES,
    RuleContext,
    Rules,
    TFMotifs,
    _edge_distance,
    _parse_gimme_scan_output,
    _parse_deeploc_csv,
    _parse_targetp_output,
    _targetp_call,
)
from tests.fixtures import DefaultsFixture


LEADER_MTP_LABEL = "Leader().betweenAA(-30, 3).is_mTP()"
LEADER_CALL_COLUMNS = ("Leader.call('noTP')", "Leader.call('mTP')", "Leader.call('SP')")
TARGETP_NO_TP_COLUMNS = ("80", "10", "10")
TARGETP_MTP_COLUMNS = ("10", "80", "10")
TARGETP_SP_COLUMNS = ("10", "10", "80")


class TestRuleContextHMMProfiles(unittest.TestCase):

    def test_batches_multiple_sequences_in_one_hmmalign_call(self):
        stockholm = "\n".join([
            "# STOCKHOLM 1.0",
            "candidate_1 AC-",
            "candidate_2 A-D",
            "#=GC RF xxx",
            "//",
            "",
        ])
        completed = CompletedProcess(["hmmalign"], 0, stdout=stockholm, stderr="")

        with patch("sieve.protein.subprocess.run", return_value=completed) as run:
            alignments = hmm_align_sequences(
                "profile.hmm",
                {"candidate_1": "AC", "candidate_2": "AD"},
            )

        run.assert_called_once()
        self.assertEqual(alignments["candidate_1"].aa_at_hmm_pos_1b(2), (2, "C"))
        self.assertIsNone(alignments["candidate_1"].aa_at_hmm_pos_1b(3))
        self.assertIsNone(alignments["candidate_2"].aa_at_hmm_pos_1b(2))
        self.assertEqual(alignments["candidate_2"].aa_at_hmm_pos_1b(3), (2, "D"))

    def test_hmmalign_failure_includes_stderr(self):
        completed = CompletedProcess(
            ["hmmalign"],
            1,
            stdout="",
            stderr="Failed to open HMM file profile.hmm",
        )

        with patch("sieve.protein.subprocess.run", return_value=completed):
            with self.assertRaisesRegex(RuntimeError, "Failed to open HMM file profile.hmm"):
                hmm_align_sequences("profile.hmm", {"candidate_1": "MA"})

    def test_resolves_bare_profile_to_registered_absolute_path(self):
        protein = FastaProtein("p1", "MA")
        profile = os.path.abspath(os.path.join("assets", "model.hmm"))
        alignment = object()

        with patch.object(protein, "hmm_align", return_value=alignment) as hmm_align:
            context = RuleContext(protein, hmm_profiles=[profile])
            self.assertIs(context.hmm_alignment("model.hmm"), alignment)

        hmm_align.assert_called_once_with(profile)

    def test_does_not_replace_profile_that_includes_a_directory(self):
        protein = FastaProtein("p1", "MA")
        registered = os.path.abspath(os.path.join("assets", "model.hmm"))
        explicit = os.path.join("other", "model.hmm")

        with patch.object(protein, "hmm_align", return_value=object()) as hmm_align:
            RuleContext(protein, hmm_profiles=[registered]).hmm_alignment(explicit)

        hmm_align.assert_called_once_with(explicit)

    def test_rejects_registered_profiles_with_ambiguous_basenames(self):
        protein = FastaProtein("p1", "MA")
        profiles = [
            os.path.join(os.sep, "models-a", "model.hmm"),
            os.path.join(os.sep, "models-b", "model.hmm"),
        ]

        with self.assertRaisesRegex(ValueError, "Ambiguous HMM profile basename 'model.hmm'"):
            RuleContext(protein, hmm_profiles=profiles)


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

    def write_protein_fixture(self, protein_accession, genome_accession, sequence="MGP"):
        self.write_manifest([self.manifest_row(protein_accession, genome_accession)])
        self.write_ncbi_proteins(genome_accession, {protein_accession: sequence})

    def write_three_exon_gene(self, protein_accession, genome_accession, strand="+", sequence="MGP"):
        self.write_manifest([self.manifest_row(protein_accession, genome_accession)])
        self.write_ncbi_proteins(genome_accession, {protein_accession: sequence})
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

    def leader_call_columns(self, row):
        return tuple(row[column] for column in LEADER_CALL_COLUMNS)

    def fake_targetp_with_expected_ids(self, expected_ids, predictions=None):
        if predictions is None:
            predictions = {}

        def fake_run(cmd, check, capture_output, text):
            fasta_path = cmd[cmd.index("-v") + 1].split(":", 1)[0] + "/query.faa"
            with open(fasta_path, "r", encoding="utf-8") as f:
                fasta_text = f.read()
            ids = [line[1:].strip() for line in fasta_text.splitlines() if line.startswith(">")]
            self.assertEqual(ids, expected_ids)
            output = ["# TargetP-2.0"]
            for sequence_id in ids:
                prediction = predictions.get(sequence_id, "mTP")
                probabilities = {
                    "noTP": {"noTP": 0.8, "SP": 0.1, "mTP": 0.1},
                    "SP": {"noTP": 0.1, "SP": 0.8, "mTP": 0.1},
                    "mTP": {"noTP": 0.1, "SP": 0.1, "mTP": 0.8},
                }[prediction]
                output.append(
                    f"{sequence_id}\t{prediction}\t{probabilities['noTP']}\t{probabilities['SP']}\t{probabilities['mTP']}\t"
                )
            output.append("")
            return CompletedProcess(cmd, 0, stdout="\n".join(output), stderr="")

        return fake_run

    def write_deeploc_csv(self, path, rows):
        fieldnames = [
            "Protein_ID",
            "Localizations",
            "Signals",
            "Membrane types",
            "Cytoplasm",
            "Endoplasmic reticulum",
            "Soluble",
        ]
        with open(path, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow({
                    "Protein_ID": "",
                    "Localizations": "",
                    "Signals": "",
                    "Membrane types": "",
                    "Cytoplasm": "",
                    "Endoplasmic reticulum": "",
                    "Soluble": "",
                } | row)

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

    def test_rules_check_includes_contig_accession_when_locus_is_available(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")

        with tempfile.TemporaryDirectory() as tmpd:
            out = os.path.join(tmpd, "rules.tsv")
            rows = Rules(Pfam.matches("PF00001")).check([("p1", "g1")], out)
            tsv_rows = self.read_tsv(out)

        self.assertEqual(rows[0]["contig accession"], "ctg1")
        self.assertEqual(tsv_rows[0]["contig accession"], "ctg1")

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
        rule = Pfam.matches("PF00001") & Leader().is_mTP()

        with patch("sieve.rules.subprocess.run", side_effect=RuntimeError("no docker")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(rule).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["pass all"], RULE_FALSE)
        self.assertEqual(rows[0][LEADER_MTP_LABEL], RULE_ERROR)

    def test_rules_check_writes_artifacts_for_batch_tools(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        targetp_output = "\n".join([
            "# TargetP-2.0",
            "p1\tmTP\t0.1\t0.2\t0.7\tCS pos: 1-2. AA-AA. Pr: 0.1",
            "",
        ])
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "p1_locus 12 17 GM.5.0.Rel.0001 8.0 +",
            "p1_locus 25 30 GM.5.0.bZIP.0001 8.0 -",
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
                    Leader().is_mTP()
                    & TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron()
                ).check([("p1", "g1")], out, artifacts_dir=artifacts)

                self.assertIn(RULE_TRUE, [row["pass all"] for row in rows])
                artifact_dirs = os.listdir(artifacts)
                leader_dir = os.path.join(artifacts, next(d for d in artifact_dirs if d.startswith("Leader_")))
                tf_dir = os.path.join(artifacts, next(d for d in artifact_dirs if d.startswith("TFMotifs.has_within")))

                self.assertTrue(os.path.exists(os.path.join(leader_dir, "query.faa")))
                self.assertTrue(os.path.exists(os.path.join(leader_dir, "command.txt")))
                with open(os.path.join(leader_dir, "stdout.txt"), "r", encoding="utf-8") as f:
                    self.assertEqual(f.read(), targetp_output)
                with open(os.path.join(leader_dir, "stderr.txt"), "r", encoding="utf-8") as f:
                    self.assertEqual(f.read(), "targetp stderr")

                self.assertTrue(os.path.exists(os.path.join(tf_dir, "locus.fna")))
                with open(os.path.join(tf_dir, "locus.fna"), "r", encoding="utf-8") as f:
                    self.assertIn(">p1_locus\n", f.read())
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
            return CompletedProcess(
                cmd,
                0,
                stdout="p1\tmTP\t0.1\t0.2\t0.7\t\n",
                stderr="",
            )

        with patch("sieve.rules.subprocess.run", side_effect=fake_run):
            with tempfile.TemporaryDirectory() as tmpd:
                old_cwd = os.getcwd()
                try:
                    os.chdir(tmpd)
                    rows = Rules(Leader().is_mTP()).check(
                        [("p1", "g1")],
                        "rules.tsv",
                        artifacts_dir="tmp",
                    )
                finally:
                    os.chdir(old_cwd)

        self.assertIn(RULE_TRUE, [row["pass all"] for row in rows])

    def test_or_and_and_composite_results_use_true_false_maybe_error(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "p1_locus 45 50 GM.5.0.Rel.0001 8.0 +",
            "",
        ])
        rule = (
            TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron(2)
            | Leader().is_mTP()
        ) & Pfam.matches("PF00001")

        def fake_run(cmd, check, capture_output, text):
            if cmd[0] == "gimme":
                return CompletedProcess(cmd, 0, stdout=gimme_output, stderr="")
            raise RuntimeError("targetp failure")

        with patch("sieve.rules.subprocess.run", side_effect=fake_run):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(rule).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["pass all"], RULE_FALSE)
        self.assertEqual(rows[0][LEADER_MTP_LABEL], RULE_ERROR)
        self.assertEqual(
            rows[0]["TFMotifs.has_within(20, 'GM.5.0.Rel', 'GM.5.0.bZIP', min_score_threshold=8).in_intron(2)"],
            "missing_GM.5.0.bZIP",
        )

    def test_hmm_position_motif_and_coverage_rules_are_strict_about_gaps(self):
        self.fx.write_protein_fixture("p1", "g1")
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("ACD-EFNGG"), id="p1"),
        ])
        alignment.column_annotations["reference_annotation"] = "xxxxxxxxx"
        protein_alignment = ProteinHMMAlignment(alignment, "p1")

        with patch.object(CuratedProtein, "hmm_align", return_value=protein_alignment):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    HMMAlignment("/models/profile.hmm").is_at("ACD", 1)
                    & HMMAlignment("/models/profile.hmm").is_at("DE", 3)
                    & HMMAlignment("/models/profile.hmm").covers(1, 3)
                    & HMMAlignment("/models/profile.hmm").covers(1, 4)
                    & HMMAlignment("/models/profile.hmm").covers(1, 5).between(1, 4)
                    & HMMAlignment("/models/profile.hmm").spans(1, 5).between(1, 4)
                    & HMMAlignment("/models/profile.hmm").spans(1, 5).between(2, 4)
                    & HMMAlignment("/models/profile.hmm").matches_regex("EFN[AGST]G", 5)
                    & HMMAlignment("/models/profile.hmm").matches_regex("EFN[AGST]G", 6)
                    & HMMAlignment("/models/profile.hmm").matches_regex("D.E", 3)
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["HMMAlignment('profile.hmm').is_at('ACD', 1)"], RULE_TRUE)
        self.assertEqual(rows[0]["HMMAlignment('profile.hmm').is_at('DE', 3)"], RULE_FALSE)
        self.assertEqual(rows[0]["HMMAlignment('profile.hmm').covers(1, 3)"], RULE_TRUE)
        self.assertEqual(rows[0]["HMMAlignment('profile.hmm').covers(1, 4)"], RULE_FALSE)
        self.assertEqual(rows[0]["HMMAlignment('profile.hmm').covers(1, 5).between(1, 4)"], RULE_FALSE)
        self.assertEqual(rows[0]["HMMAlignment('profile.hmm').spans(1, 5).between(1, 4)"], RULE_TRUE)
        self.assertEqual(rows[0]["HMMAlignment('profile.hmm').spans(1, 5).between(2, 4)"], RULE_FALSE)
        self.assertEqual(rows[0]["HMMAlignment('profile.hmm').matches_regex('EFN[AGST]G', 5)"], RULE_TRUE)
        self.assertEqual(rows[0]["HMMAlignment('profile.hmm').matches_regex('EFN[AGST]G', 6)"], RULE_FALSE)
        self.assertEqual(rows[0]["HMMAlignment('profile.hmm').matches_regex('D.E', 3)"], RULE_TRUE)

    def test_hmm_rules_align_leader_candidates_in_one_profile_batch(self):
        protein = FastaProtein("p1", "MAM")
        rule = (
            Leader().betweenAA(1, 3)
            & HMMAlignment("profile.hmm").spans(1, 2).between(1, 3)
        )

        def align_candidates(_profile, sequences_by_id):
            alignments = {}
            for sequence_id, sequence in sequences_by_id.items():
                alignment = MultipleSeqAlignment([
                    SeqRecord(Seq(sequence), id=sequence_id),
                ])
                alignment.column_annotations["reference_annotation"] = "x" * len(sequence)
                alignments[sequence_id] = ProteinHMMAlignment(alignment, sequence_id)
            return alignments

        with patch("sieve.rules.hmm_align_sequences", side_effect=align_candidates) as hmmalign:
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(rule).check_proteins([protein], os.path.join(tmpd, "rules.tsv"))

        hmmalign.assert_called_once()
        label = "HMMAlignment('profile.hmm').spans(1, 2).between(1, 3)"
        self.assertEqual(
            [(row["sequence accession"], row[label], row["pass all"]) for row in rows],
            [
                ("p1_with_leader_1_M", RULE_TRUE, RULE_TRUE),
                ("p1_with_leader_3_M", RULE_FALSE, RULE_FALSE),
            ],
        )

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
            self.assertIn(">p1\nMA\n", fasta_text)
            self.assertIn(">p2\nMG\n", fasta_text)
            return CompletedProcess(cmd, 0, stdout="\n".join([
                "# TargetP-2.0",
                "# ID\tPrediction\tnoTP\tSP\tmTP\tCS Position",
                "p1\tmTP\t0.1\t0.2\t0.7\tCS pos: 1-2. AA-AA. Pr: 0.1",
                "p2\tSP\t0.1\t0.8\t0.1\tCS pos: 1-2. AA-AA. Pr: 0.1",
                "",
            ]), stderr="")

        with patch("sieve.rules.subprocess.run", side_effect=fake_run):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(Leader().is_mTP()).check(
                    [("p1", "g1"), ("p2", "g2")],
                    os.path.join(tmpd, "rules.tsv"),
                )

        self.assertEqual(
            [(row["sequence accession"], row[LEADER_MTP_LABEL], self.leader_call_columns(row)) for row in rows],
            [
                ("p1", RULE_TRUE, ("10", "70", "20")),
                ("p2", RULE_FALSE, TARGETP_SP_COLUMNS),
            ],
        )

    def test_leader_rule_checks_all_alternative_start_candidates(self):
        self.fx.write_manifest([
            self.fx.manifest_row("p1", "g1", source=SEQUENCE_SOURCE_HMM_DETECTED),
        ])
        self.fx.write_detected_proteins("g1", {"p1": "KMMP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "NNNATGATGAAAATGATGCCC"})
        self.fx.write_detected_rows("g1", [
            self.fx.detected_protein_row("p1", "g1", "ctg1", 10, 21, 1, 4),
        ])

        def fake_run(cmd, check, capture_output, text):
            fasta_path = cmd[cmd.index("-v") + 1].split(":", 1)[0] + "/query.faa"
            with open(fasta_path, "r", encoding="utf-8") as f:
                fasta_text = f.read()
            self.assertIn(">p1_with_leader_u2_M\nMMKMMP\n", fasta_text)
            self.assertIn(">p1_with_leader_u1_M\nMKMMP\n", fasta_text)
            self.assertIn(">p1_with_leader_2_M\nMMP\n", fasta_text)
            self.assertIn(">p1_with_leader_3_M\nMP\n", fasta_text)
            return CompletedProcess(cmd, 0, stdout="\n".join([
                "# TargetP-2.0",
                "p1_with_leader_u2_M\tnoTP\t0.8\t0.1\t0.1\t",
                "p1_with_leader_u1_M\tSP\t0.1\t0.8\t0.1\t",
                "p1_with_leader_2_M\tmTP\t0.1\t0.1\t0.8\t",
                "p1_with_leader_3_M\tnoTP\t0.8\t0.1\t0.1\t",
                "",
            ]), stderr="")

        with patch("sieve.rules.subprocess.run", side_effect=fake_run):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(Leader().is_mTP()).check(
                    [("p1", "g1")],
                    os.path.join(tmpd, "rules.tsv"),
                )

        self.assertEqual(
            [(row["sequence accession"], row[LEADER_MTP_LABEL], self.leader_call_columns(row)) for row in rows],
            [
                ("p1_with_leader_u2_M", RULE_FALSE, TARGETP_NO_TP_COLUMNS),
                ("p1_with_leader_u1_M", RULE_FALSE, TARGETP_SP_COLUMNS),
                ("p1_with_leader_2_M", RULE_TRUE, TARGETP_MTP_COLUMNS),
                ("p1_with_leader_3_M", RULE_FALSE, TARGETP_NO_TP_COLUMNS),
            ],
        )

    def test_leader_rule_writes_low_confidence_targetp_probabilities_to_call_annotation(self):
        self.fx.write_manifest([
            self.fx.manifest_row("p1", "g1", source=SEQUENCE_SOURCE_HMM_DETECTED),
        ])
        self.fx.write_detected_proteins("g1", {"p1": "MP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "ATGCCC"})
        self.fx.write_detected_rows("g1", [
            self.fx.detected_protein_row("p1", "g1", "ctg1", 1, 6, 1, 2),
        ])

        targetp_output = "\n".join([
            "# TargetP-2.0",
            "# ID\tPrediction\tnoTP\tSP\tmTP\tCS Position",
            "p1_with_leader_1_M\tnoTP\t0.6\t0.3\t0.1\t",
            "",
        ])

        def fake_run(cmd, check, capture_output, text):
            return CompletedProcess(cmd, 0, stdout=targetp_output, stderr="")

        with patch("sieve.rules.subprocess.run", side_effect=fake_run):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(Leader().is_noTP()).check(
                    [("p1", "g1")],
                    os.path.join(tmpd, "rules.tsv"),
                )

        self.assertEqual(self.leader_call_columns(rows[0]), ("60", "10", "30"))

    def test_leader_rule_calls_targetp_on_original_sequence_without_any_m_start_candidates(self):
        self.fx.write_manifest([
            self.fx.manifest_row("p1", "g1"),
        ])
        self.fx.write_ncbi_proteins("g1", {"p1": "RYDA"})

        def fake_run(cmd, check, capture_output, text):
            fasta_path = cmd[cmd.index("-v") + 1].split(":", 1)[0] + "/query.faa"
            with open(fasta_path, "r", encoding="utf-8") as f:
                fasta_text = f.read()
            self.assertIn(">p1\nRYDA\n", fasta_text)
            return CompletedProcess(cmd, 0, stdout="p1\tnoTP\t0.8\t0.1\t0.1\t\n", stderr="")

        with patch("sieve.rules.subprocess.run", side_effect=fake_run):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(Leader().is_mTP()).check(
                    [("p1", "g1")],
                    os.path.join(tmpd, "rules.tsv"),
                )

        self.assertEqual(rows[0][LEADER_MTP_LABEL], RULE_FALSE)
        self.assertEqual(rows[0]["sequence accession"], "p1")
        self.assertEqual(self.leader_call_columns(rows[0]), TARGETP_NO_TP_COLUMNS)

    def test_leader_rule_writes_call_annotation_to_tsv(self):
        self.fx.write_manifest([
            self.fx.manifest_row("p1", "g1"),
        ])
        self.fx.write_ncbi_proteins("g1", {"p1": "MA"})

        def fake_run(cmd, check, capture_output, text):
            return CompletedProcess(cmd, 0, stdout="\n".join([
                "# TargetP-2.0",
                "# ID\tPrediction\tnoTP\tSP\tmTP\tCS Position",
                "p1\tnoTP\t0.8\t0.1\t0.1\t",
                "",
            ]), stderr="")

        with patch("sieve.rules.subprocess.run", side_effect=fake_run):
            with tempfile.TemporaryDirectory() as tmpd:
                out = os.path.join(tmpd, "rules.tsv")
                Rules(Leader().is_mTP()).check([("p1", "g1")], out)
                rows = self.read_tsv(out)

        self.assertEqual(
            [(row["sequence accession"], row[LEADER_MTP_LABEL], self.leader_call_columns(row)) for row in rows],
            [
                ("p1", RULE_FALSE, TARGETP_NO_TP_COLUMNS),
            ],
        )

    def test_leader_default_window_checks_candidates_between_minus_30_and_position_3(self):
        self.fx.write_manifest([
            self.fx.manifest_row("p1", "g1"),
        ])
        self.fx.write_ncbi_proteins("g1", {"p1": "MAMMAM"})

        with patch(
            "sieve.rules.subprocess.run",
            side_effect=self.fake_targetp_with_expected_ids([
                "p1",
                "p1_with_leader_3_M",
            ], {"p1": "noTP"}),
        ):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(Leader().is_mTP()).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(
            [(row["sequence accession"], row[LEADER_MTP_LABEL], self.leader_call_columns(row)) for row in rows],
            [
                ("p1", RULE_FALSE, TARGETP_NO_TP_COLUMNS),
                ("p1_with_leader_3_M", RULE_TRUE, TARGETP_MTP_COLUMNS),
            ],
        )

    def test_leader_between_aa_handles_zero_negative_and_positive_boundaries(self):
        self.fx.write_manifest([
            self.fx.manifest_row("p1", "g1", source=SEQUENCE_SOURCE_HMM_DETECTED),
        ])
        self.fx.write_detected_proteins("g1", {"p1": "MMMP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "NNNATGATGATGATGCCC"})
        self.fx.write_detected_rows("g1", [
            self.fx.detected_protein_row("p1", "g1", "ctg1", 7, 18, 1, 4),
        ])

        cases = [
            ((-1, 1), ["p1_with_leader_u1_M", "p1_with_leader_1_M"]),
            ((0, 1), ["p1_with_leader_1_M"]),
            ((-1, 0), ["p1_with_leader_u1_M"]),
        ]
        for (start, end), expected_ids in cases:
            with self.subTest(start=start, end=end):
                with patch(
                    "sieve.rules.subprocess.run",
                    side_effect=self.fake_targetp_with_expected_ids(expected_ids),
                ):
                    with tempfile.TemporaryDirectory() as tmpd:
                        rows = Rules(Leader().betweenAA(start, end).is_mTP()).check(
                            [("p1", "g1")],
                            os.path.join(tmpd, "rules.tsv"),
                        )

                self.assertEqual(rows[0][f"Leader().betweenAA({start}, {end}).is_mTP()"], RULE_TRUE)

    def test_leader_upstream_of_pfam_uses_target_start_to_infer_anchor(self):
        self.fx.write_manifest([
            self.fx.manifest_row("p1", "g1"),
        ])
        self.fx.write_ncbi_proteins("g1", {"p1": "AAMMAMAMMMGG"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "GCTGCTATGATGGCTATGGCTATGATGATGGGTGGT"})
        self.fx.write_gff("g1", "\n".join([
            "ctg1\tsrc\tmRNA\t1\t36\t.\t+\t.\tID=tx1",
            "ctg1\tsrc\tCDS\t1\t36\t.\t+\t0\tID=cds1;Parent=tx1;protein_id=p1",
            "",
        ]))
        row = self.fx.detected_row("p1", "g1", "PF00081.28", "Pfam")
        row.update(query_start=10, query_end=30, target_start=5, target_end=25)
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [row])

        with patch(
            "sieve.rules.subprocess.run",
            side_effect=self.fake_targetp_with_expected_ids([
                "p1_with_leader_u3_PF00081_M",
                "p1_with_leader_u2_PF00081_M",
            ]),
        ):
            with tempfile.TemporaryDirectory() as tmpd:
                rule = Leader().upstreamOfPfam("PF00081.28").betweenAA(-3, 0).is_mTP()
                rows = Rules(rule).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        label = "Leader().upstreamOfPfam('PF00081.28').betweenAA(-3, 0).is_mTP()"
        self.assertEqual(
            [(row["sequence accession"], row[label], self.leader_call_columns(row)) for row in rows],
            [
                ("p1_with_leader_u3_PF00081_M", RULE_TRUE, TARGETP_MTP_COLUMNS),
                ("p1_with_leader_u2_PF00081_M", RULE_TRUE, TARGETP_MTP_COLUMNS),
            ],
        )
        self.assertEqual(
            [
                candidate.accession
                for row in rows
                for candidate in Rules(rule).sequence_candidates_for_row(row)
            ],
            [
                "p1_with_leader_u3_PF00081_M",
                "p1_with_leader_u2_PF00081_M",
            ],
        )

    def test_bare_leader_discovers_all_methionines_in_inclusive_pfam_window(self):
        pfam_row = self.fx.detected_row("p1", "", "PF00081.28", "Pfam")
        pfam_row.update(query_start=6, query_end=15, target_start=1, target_end=10)
        protein = FastaProtein("p1", "MAMAAA", pfam_rows=[pfam_row])

        with patch("sieve.rules.subprocess.run") as run:
            with tempfile.TemporaryDirectory() as tmpd:
                rule = Leader().upstreamOfPfam("PF00081").betweenAA(-5, -3)
                rows = Rules(rule).check_proteins([protein], os.path.join(tmpd, "rules.tsv"))

        run.assert_not_called()
        label = "Leader().upstreamOfPfam('PF00081').betweenAA(-5, -3)"
        self.assertEqual(
            [(row["sequence accession"], row[label], row["pass all"]) for row in rows],
            [
                ("p1_with_leader_u5_PF00081_M", RULE_TRUE, RULE_TRUE),
                ("p1_with_leader_u3_PF00081_M", RULE_TRUE, RULE_TRUE),
            ],
        )

    def test_bare_leader_rejects_when_pfam_window_contains_no_methionine(self):
        pfam_row = self.fx.detected_row("p1", "", "PF00081.28", "Pfam")
        pfam_row.update(query_start=6, query_end=15, target_start=1, target_end=10)
        protein = FastaProtein("p1", "AAAAAA", pfam_rows=[pfam_row])

        with tempfile.TemporaryDirectory() as tmpd:
            rule = Leader().upstreamOfPfam("PF00081").betweenAA(-5, -1)
            rows = Rules(rule).check_proteins([protein], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows, [])

    def test_bare_leader_rejects_without_matching_pfam_hit(self):
        protein = FastaProtein("p1", "MAMAAA", pfam_rows=[])

        with tempfile.TemporaryDirectory() as tmpd:
            rule = Leader().upstreamOfPfam("PF00081").betweenAA(-5, -1)
            rows = Rules(rule).check_proteins([protein], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows, [])

    def test_leader_upstream_of_pfam_for_hmm_detected_ignores_hmm_profile_coordinates(self):
        self.fx.write_manifest([
            self.fx.manifest_row("p1", "g1", source=SEQUENCE_SOURCE_HMM_DETECTED),
        ])
        self.fx.write_detected_proteins("g1", {"p1": "AAA"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "ATGAAAGCTGCTGCT"})
        self.fx.write_detected_rows("g1", [
            self.fx.detected_protein_row("p1", "g1", "ctg1", 7, 15, 35, 37),
        ])
        row = self.fx.detected_row("p1", "g1", "PF00081.28", "Pfam")
        row.update(query_start=3, query_end=20, target_start=2, target_end=19)
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [row])

        with patch(
            "sieve.rules.subprocess.run",
            side_effect=self.fake_targetp_with_expected_ids([
                "p1_with_leader_u3_PF00081_M",
            ]),
        ):
            with tempfile.TemporaryDirectory() as tmpd:
                rule = Leader().upstreamOfPfam("PF00081").betweenAA(-5, 0).is_mTP()
                rows = Rules(rule).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        label = "Leader().upstreamOfPfam('PF00081').betweenAA(-5, 0).is_mTP()"
        self.assertEqual(
            [(row["sequence accession"], row[label], self.leader_call_columns(row)) for row in rows],
            [
                ("p1_with_leader_u3_PF00081_M", RULE_TRUE, TARGETP_MTP_COLUMNS),
            ],
        )

    def test_leader_upstream_of_pfam_for_hmm_detected_handles_domain_anchor_before_protein_start(self):
        self.fx.write_manifest([
            self.fx.manifest_row("p1", "g1", source=SEQUENCE_SOURCE_HMM_DETECTED),
        ])
        self.fx.write_detected_proteins("g1", {"p1": "AAA"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "ATGAAAGCTGCTGCT"})
        self.fx.write_detected_rows("g1", [
            self.fx.detected_protein_row("p1", "g1", "ctg1", 7, 15, 35, 37),
        ])
        row = self.fx.detected_row("p1", "g1", "PF00081.28", "Pfam")
        row.update(query_start=1, query_end=20, target_start=2, target_end=21)
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [row])

        with patch(
            "sieve.rules.subprocess.run",
            side_effect=self.fake_targetp_with_expected_ids([
                "p1_with_leader_u1_PF00081_M",
            ]),
        ):
            with tempfile.TemporaryDirectory() as tmpd:
                rule = Leader().upstreamOfPfam("PF00081").betweenAA(-5, 0).is_mTP()
                rows = Rules(rule).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["Leader().upstreamOfPfam('PF00081').betweenAA(-5, 0).is_mTP()"], RULE_TRUE)
        self.assertEqual(rows[0]["sequence accession"], "p1_with_leader_u1_PF00081_M")
        self.assertEqual(self.leader_call_columns(rows[0]), TARGETP_MTP_COLUMNS)

    def test_leader_upstream_of_pfam_uses_earliest_inferred_domain_start(self):
        self.fx.write_manifest([
            self.fx.manifest_row("p1", "g1"),
        ])
        self.fx.write_ncbi_proteins("g1", {"p1": "AAAAMMAMM"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "GCTGCTGCTGCTATGATGGCTATGATG"})
        self.fx.write_gff("g1", "\n".join([
            "ctg1\tsrc\tmRNA\t1\t27\t.\t+\t.\tID=tx1",
            "ctg1\tsrc\tCDS\t1\t27\t.\t+\t0\tID=cds1;Parent=tx1;protein_id=p1",
            "",
        ]))
        later = self.fx.detected_row("p1", "g1", "PF00081.28", "Pfam")
        later.update(query_start=9, query_end=20, target_start=1, target_end=12)
        earlier = self.fx.detected_row("p1", "g1", "PF00081.28", "Pfam")
        earlier.update(query_start=6, query_end=17, target_start=1, target_end=12)
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [later, earlier])

        with patch(
            "sieve.rules.subprocess.run",
            side_effect=self.fake_targetp_with_expected_ids([
                "p1_with_leader_u1_PF00081_M",
                "p1_with_leader_1_PF00081_M",
            ]),
        ):
            with tempfile.TemporaryDirectory() as tmpd:
                rule = Leader().upstreamOfPfam("PF00081").betweenAA(-1, 1).is_mTP()
                rows = Rules(rule).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        label = "Leader().upstreamOfPfam('PF00081').betweenAA(-1, 1).is_mTP()"
        self.assertEqual(
            [(row["sequence accession"], row[label], self.leader_call_columns(row)) for row in rows],
            [
                ("p1_with_leader_u1_PF00081_M", RULE_TRUE, TARGETP_MTP_COLUMNS),
                ("p1_with_leader_1_PF00081_M", RULE_TRUE, TARGETP_MTP_COLUMNS),
            ],
        )

    def test_leader_upstream_of_pfam_rejects_without_matching_pfam_hit(self):
        self.fx.write_manifest([
            self.fx.manifest_row("p1", "g1"),
        ])
        self.fx.write_ncbi_proteins("g1", {"p1": "MMA"})
        row = self.fx.detected_row("p1", "g1", "PF99999.1", "Pfam")
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [row])

        with patch("sieve.rules.subprocess.run") as run:
            with tempfile.TemporaryDirectory() as tmpd:
                rule = Leader().upstreamOfPfam("PF00081").betweenAA(-30, 0).is_mTP()
                rows = Rules(rule).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        run.assert_not_called()
        self.assertEqual(rows, [])

    def test_leader_upstream_of_pfam_fallback_original_must_be_in_anchor_window(self):
        row = self.fx.detected_row("p1", "", "PF00081.28", "Pfam")
        row.update(query_start=5, query_end=15, target_start=1, target_end=11)
        protein = FastaProtein("p1", "AAAAAA", pfam_rows=[row])

        with patch(
            "sieve.rules.subprocess.run",
            side_effect=self.fake_targetp_with_expected_ids(["p1"]),
        ):
            with tempfile.TemporaryDirectory() as tmpd:
                rule = Leader().upstreamOfPfam("PF00081").betweenAA(-5, -3).is_mTP()
                rows = Rules(rule).check_proteins([protein], os.path.join(tmpd, "rules.tsv"))

        label = "Leader().upstreamOfPfam('PF00081').betweenAA(-5, -3).is_mTP()"
        self.assertEqual(
            [(row["sequence accession"], row[label], self.leader_call_columns(row)) for row in rows],
            [("p1", RULE_TRUE, TARGETP_MTP_COLUMNS)],
        )

    def test_leader_upstream_of_pfam_rejects_fallback_original_outside_anchor_window(self):
        row = self.fx.detected_row("p1", "", "PF00081.28", "Pfam")
        row.update(query_start=5, query_end=15, target_start=1, target_end=11)
        protein = FastaProtein("p1", "AAAAAA", pfam_rows=[row])

        with patch("sieve.rules.subprocess.run") as run:
            with tempfile.TemporaryDirectory() as tmpd:
                rule = Leader().upstreamOfPfam("PF00081").betweenAA(-3, -1).is_mTP()
                rows = Rules(rule).check_proteins([protein], os.path.join(tmpd, "rules.tsv"))

        run.assert_not_called()
        self.assertEqual(rows, [])

    def test_leader_deeploc_mtp_uses_signal_and_score_annotations(self):
        protein = FastaProtein("p1", "MMT")
        with tempfile.TemporaryDirectory() as tmpd:
            deeploc_csv = os.path.join(tmpd, "deeploc.csv")
            self.write_deeploc_csv(deeploc_csv, [
                {
                    "Protein_ID": "p1",
                    "Localizations": "Cytoplasm",
                    "Signals": "",
                    "Cytoplasm": "0.7658",
                    "Endoplasmic reticulum": "0.2588",
                    "Soluble": "0.7442",
                },
                {
                    "Protein_ID": "p1_with_leader_2_M",
                    "Localizations": "Endoplasmic reticulum",
                    "Signals": "Mitochondrial transit peptide|Signal peptide",
                    "Cytoplasm": "0.1",
                    "Endoplasmic reticulum": "0.88",
                    "Soluble": "0.2",
                },
            ])
            rows = Rules(Leader().is_mTP(deeploc=True)).check_proteins(
                [protein],
                os.path.join(tmpd, "rules.tsv"),
                deeploc_csv=deeploc_csv,
            )

            output_rows = self.read_tsv(os.path.join(tmpd, "rules.tsv"))

        label = "Leader().betweenAA(-30, 3).is_mTP(deeploc=True)"
        self.assertEqual(
            [
                (
                    row["sequence accession"],
                    row[label],
                    row["Leader.call('mTP')"],
                    row["Leader.call('SP')"],
                    row["Leader.localization"],
                    row["Leader.call('Endoplasmic reticulum')"],
                )
                for row in rows
            ],
            [
                ("p1", RULE_FALSE, "0", "0", "Cytoplasm", "26"),
                ("p1_with_leader_2_M", RULE_TRUE, "100", "100", "Endoplasmic reticulum", "88"),
            ],
        )
        self.assertEqual(output_rows[1]["Leader.call('Soluble')"], "20")

    def test_leader_deeploc_sp_uses_signal_peptide(self):
        protein = FastaProtein("p1", "MMT")
        with tempfile.TemporaryDirectory() as tmpd:
            deeploc_csv = os.path.join(tmpd, "deeploc.csv")
            self.write_deeploc_csv(deeploc_csv, [
                {"Protein_ID": "p1", "Signals": "", "Localizations": "Cytoplasm"},
                {
                    "Protein_ID": "p1_with_leader_2_M",
                    "Signals": "Signal peptide|Peroxisomal targeting signal",
                    "Localizations": "Extracellular",
                },
            ])
            rows = Rules(Leader().is_SP(deeploc=True)).check_proteins(
                [protein],
                os.path.join(tmpd, "rules.tsv"),
                deeploc_csv=deeploc_csv,
            )

        label = "Leader().betweenAA(-30, 3).is_SP(deeploc=True)"
        self.assertEqual(
            [(row["sequence accession"], row[label], row["Leader.call('SP')"]) for row in rows],
            [("p1", RULE_FALSE, "0"), ("p1_with_leader_2_M", RULE_TRUE, "100")],
        )

    def test_leader_deeploc_localize_at_uses_exact_localization(self):
        protein = FastaProtein("p1", "MMT")
        with tempfile.TemporaryDirectory() as tmpd:
            deeploc_csv = os.path.join(tmpd, "deeploc.csv")
            self.write_deeploc_csv(deeploc_csv, [
                {"Protein_ID": "p1", "Localizations": "Endoplasmic reticulum", "Endoplasmic reticulum": "0.05"},
                {"Protein_ID": "p1_with_leader_2_M", "Localizations": "Cytoplasm", "Endoplasmic reticulum": "0.95"},
            ])
            rows = Rules(Leader().localize_at("Endoplasmic reticulum")).check_proteins(
                [protein],
                os.path.join(tmpd, "rules.tsv"),
                deeploc_csv=deeploc_csv,
            )

        label = "Leader().betweenAA(-30, 3).localize_at('Endoplasmic reticulum')"
        self.assertEqual(
            [
                (
                    row["sequence accession"],
                    row[label],
                    row["Leader.localization"],
                    row["Leader.call('Endoplasmic reticulum')"],
                )
                for row in rows
            ],
            [
                ("p1", RULE_TRUE, "Endoplasmic reticulum", "5"),
                ("p1_with_leader_2_M", RULE_FALSE, "Cytoplasm", "95"),
            ],
        )

    def test_leader_deeploc_rules_require_csv_when_evaluating(self):
        protein = FastaProtein("p1", "MMT")
        with tempfile.TemporaryDirectory() as tmpd:
            with self.assertRaisesRegex(ValueError, "--deeploc-csv is required"):
                Rules(Leader().is_mTP(deeploc=True)).check_proteins(
                    [protein],
                    os.path.join(tmpd, "rules.tsv"),
                )

    def test_leader_deeploc_missing_candidate_row_is_error(self):
        protein = FastaProtein("p1", "MMT")
        with tempfile.TemporaryDirectory() as tmpd:
            deeploc_csv = os.path.join(tmpd, "deeploc.csv")
            self.write_deeploc_csv(deeploc_csv, [
                {"Protein_ID": "p1", "Signals": "Mitochondrial transit peptide"},
            ])
            rows = Rules(Leader().is_mTP(deeploc=True)).check_proteins(
                [protein],
                os.path.join(tmpd, "rules.tsv"),
                deeploc_csv=deeploc_csv,
            )

        self.assertEqual(rows[0]["Leader().betweenAA(-30, 3).is_mTP(deeploc=True)"], RULE_ERROR)

    def test_leader_upstream_of_pfam_uses_ncbi_spliced_prefix_before_anchor(self):
        cases = [
            ("+", "ATGAAACCCCCCCCCCCGCTGCTGCT", [(1, 6), (16, 27)]),
            ("-", "AGCAGCAGCGGGGGGGGGTTTCAT", [(19, 24), (1, 12)]),
        ]
        for strand, contig_sequence, cds_intervals in cases:
            with self.subTest(strand=strand):
                CuratedProtein.clear_cache()
                self.fx.write_manifest([
                    self.fx.manifest_row("p1", "g1", source=SEQUENCE_SOURCE_NCBI),
                ])
                self.fx.write_ncbi_proteins("g1", {"p1": "MKAAAA"})
                self.fx.write_genomic_fasta("g1", {"ctg1": contig_sequence})
                cds_rows = [
                    f"ctg1\tsrc\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds{i};Parent=tx1;protein_id=p1"
                    for i, (start, end) in enumerate(cds_intervals, start=1)
                ]
                self.fx.write_gff("g1", "\n".join([
                    f"ctg1\tsrc\tmRNA\t1\t{len(contig_sequence)}\t.\t{strand}\t.\tID=tx1",
                    *cds_rows,
                    "",
                ]))
                row = self.fx.detected_row("p1", "g1", "PF00081.28", "Pfam")
                row.update(query_start=5, query_end=12, target_start=1, target_end=8)
                DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [row])

                with patch(
                    "sieve.rules.subprocess.run",
                    side_effect=self.fake_targetp_with_expected_ids([
                        "p1",
                    ], {"p1": "noTP"}),
                ):
                    with tempfile.TemporaryDirectory() as tmpd:
                        rule = Leader().upstreamOfPfam("PF00081").betweenAA(-5, 0).is_mTP()
                        rows = Rules(rule).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

                label = "Leader().upstreamOfPfam('PF00081').betweenAA(-5, 0).is_mTP()"
                self.assertEqual(
                    [(row["sequence accession"], row[label], self.leader_call_columns(row)) for row in rows],
                    [
                        ("p1", RULE_FALSE, TARGETP_NO_TP_COLUMNS),
                    ],
                )

    def test_leader_upstream_of_pfam_uses_hmm_detected_spliced_prefix_before_anchor(self):
        cases = [
            ("+", "ATGAAACCCCCCCCCCCGCTGCTGCT", [(1, 6), (16, 27)]),
            ("-", "AGCAGCAGCGGGGGGGGGTTTCAT", [(24, 19), (12, 1)]),
        ]
        for strand, contig_sequence, rowspecs in cases:
            with self.subTest(strand=strand):
                CuratedProtein.clear_cache()
                self.fx.write_manifest([
                    self.fx.manifest_row("p1", "g1", source=SEQUENCE_SOURCE_HMM_DETECTED),
                ])
                self.fx.write_detected_proteins("g1", {"p1": "MKAAAA"})
                self.fx.write_genomic_fasta("g1", {"ctg1": contig_sequence})
                self.fx.write_detected_rows("g1", [
                    self.fx.detected_protein_row("p1", "g1", "ctg1", rowspecs[0][0], rowspecs[0][1], 37, 38),
                    self.fx.detected_protein_row("p1", "g1", "ctg1", rowspecs[1][0], rowspecs[1][1], 39, 42),
                ])
                row = self.fx.detected_row("p1", "g1", "PF00081.28", "Pfam")
                row.update(query_start=5, query_end=12, target_start=1, target_end=8)
                DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [row])

                with patch(
                    "sieve.rules.subprocess.run",
                    side_effect=self.fake_targetp_with_expected_ids([
                        "p1_with_leader_u4_PF00081_M",
                    ]),
                ):
                    with tempfile.TemporaryDirectory() as tmpd:
                        rule = Leader().upstreamOfPfam("PF00081").betweenAA(-5, 0).is_mTP()
                        rows = Rules(rule).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

                self.assertEqual(rows[0]["Leader().upstreamOfPfam('PF00081').betweenAA(-5, 0).is_mTP()"], RULE_TRUE)
                self.assertEqual(rows[0]["sequence accession"], "p1_with_leader_u4_PF00081_M")
                self.assertEqual(self.leader_call_columns(rows[0]), TARGETP_MTP_COLUMNS)

    def test_leader_upstream_of_pfam_for_hmm_detected_also_uses_raw_upstream_anchor_context(self):
        coding_with_missing_middle = "GCTGCTATGAAAGCTGCT"
        cases = [
            ("+", coding_with_missing_middle, [(1, 6), (13, 18)]),
            ("-", str(Seq(coding_with_missing_middle).reverse_complement()), [(18, 13), (6, 1)]),
        ]
        for strand, contig_sequence, rowspecs in cases:
            with self.subTest(strand=strand):
                CuratedProtein.clear_cache()
                self.fx.write_manifest([
                    self.fx.manifest_row("p1", "g1", source=SEQUENCE_SOURCE_HMM_DETECTED),
                ])
                self.fx.write_detected_proteins("g1", {"p1": "AAAA"})
                self.fx.write_genomic_fasta("g1", {"ctg1": contig_sequence})
                self.fx.write_detected_rows("g1", [
                    self.fx.detected_protein_row("p1", "g1", "ctg1", rowspecs[0][0], rowspecs[0][1], 37, 38),
                    self.fx.detected_protein_row("p1", "g1", "ctg1", rowspecs[1][0], rowspecs[1][1], 39, 40),
                ])
                row = self.fx.detected_row("p1", "g1", "PF00081.28", "Pfam")
                row.update(query_start=4, query_end=12, target_start=1, target_end=9)
                DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [row])

                with patch(
                    "sieve.rules.subprocess.run",
                    side_effect=self.fake_targetp_with_expected_ids([
                        "p1_with_leader_u3_PF00081_anchor_M",
                    ]),
                ):
                    with tempfile.TemporaryDirectory() as tmpd:
                        rule = Leader().upstreamOfPfam("PF00081").betweenAA(-5, 0).is_mTP()
                        rows = Rules(rule).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

                self.assertEqual(rows[0]["Leader().upstreamOfPfam('PF00081').betweenAA(-5, 0).is_mTP()"], RULE_TRUE)
                self.assertEqual(rows[0]["sequence accession"], "p1_with_leader_u3_PF00081_anchor_M")
                self.assertEqual(self.leader_call_columns(rows[0]), TARGETP_MTP_COLUMNS)

    def test_leader_upstream_of_pfam_scoping_handles_ncbi_forward_and_reverse_gff_loci(self):
        coding = "GCTGCTATGATGGCTATG"
        cases = [
            ("+", coding, 1, 18),
            ("-", str(Seq(coding).reverse_complement()), 1, 18),
        ]
        for strand, contig_sequence, cds_start, cds_end in cases:
            with self.subTest(strand=strand):
                CuratedProtein.clear_cache()
                self.fx.write_manifest([
                    self.fx.manifest_row("p1", "g1", source=SEQUENCE_SOURCE_NCBI),
                ])
                self.fx.write_ncbi_proteins("g1", {"p1": "AAMMAM"})
                self.fx.write_genomic_fasta("g1", {"ctg1": contig_sequence})
                self.fx.write_gff("g1", "\n".join([
                    f"ctg1\tsrc\tmRNA\t{cds_start}\t{cds_end}\t.\t{strand}\t.\tID=tx1",
                    f"ctg1\tsrc\tCDS\t{cds_start}\t{cds_end}\t.\t{strand}\t0\tID=cds1;Parent=tx1;protein_id=p1",
                    "",
                ]))
                row = self.fx.detected_row("p1", "g1", "PF00081.28", "Pfam")
                row.update(query_start=4, query_end=10, target_start=1, target_end=7)
                DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [row])

                with patch(
                    "sieve.rules.subprocess.run",
                    side_effect=self.fake_targetp_with_expected_ids([
                        "p1_with_leader_u1_PF00081_M",
                        "p1_with_leader_1_PF00081_M",
                    ]),
                ):
                    with tempfile.TemporaryDirectory() as tmpd:
                        rule = Leader().upstreamOfPfam("PF00081").betweenAA(-1, 1).is_mTP()
                        rows = Rules(rule).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

                self.assertIn(RULE_TRUE, [row["pass all"] for row in rows])

    def test_leader_upstream_of_pfam_scoping_handles_hmm_detected_forward_and_reverse_loci(self):
        coding = "GCTGCTATGATGGCTATG"
        cases = [
            ("+", coding, 1, 18),
            ("-", str(Seq(coding).reverse_complement()), 18, 1),
        ]
        for strand, contig_sequence, query_start, query_end in cases:
            with self.subTest(strand=strand):
                CuratedProtein.clear_cache()
                self.fx.write_manifest([
                    self.fx.manifest_row("p1", "g1", source=SEQUENCE_SOURCE_HMM_DETECTED),
                ])
                self.fx.write_detected_proteins("g1", {"p1": "AAMMAM"})
                self.fx.write_genomic_fasta("g1", {"ctg1": contig_sequence})
                self.fx.write_detected_rows("g1", [
                    self.fx.detected_protein_row("p1", "g1", "ctg1", query_start, query_end, 1, 6),
                ])
                row = self.fx.detected_row("p1", "g1", "PF00081.28", "Pfam")
                row.update(query_start=4, query_end=10, target_start=1, target_end=7)
                DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [row])

                with patch(
                    "sieve.rules.subprocess.run",
                    side_effect=self.fake_targetp_with_expected_ids([
                        "p1_with_leader_u1_PF00081_M",
                        "p1_with_leader_1_PF00081_M",
                    ]),
                ):
                    with tempfile.TemporaryDirectory() as tmpd:
                        rule = Leader().upstreamOfPfam("PF00081").betweenAA(-1, 1).is_mTP()
                        rows = Rules(rule).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

                self.assertIn(RULE_TRUE, [row["pass all"] for row in rows])

    def test_tf_motifs_pass_for_hits_within_intron_on_forward_gene_even_opposite_hit_strands(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "p1_locus 45 50 GM.5.0.Rel.0001 8.0 +",
            "p1_locus 70 75 GM.5.0.bZIP.0001 9.0 -",
            "p1_locus 55 60 GM.5.0.bZIP.0002 9.0 -",
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
            "p1_locus 62 66 GM.5.0.Rel.0001 8.0 -",
            "p1_locus 75 80 GM.5.0.bZIP.0001 8.0 +",
            "p1_locus 25 30 GM.5.0.Rel.0002 8.0 +",
            "p1_locus 35 39 GM.5.0.bZIP.0002 8.0 -",
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
            "p1_locus 12 17 GM.5.0.Rel.0001 8.0 +",
            "p1_locus 25 30 GM.5.0.bZIP.0001 8.0 -",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron()
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(
            rows[0]["TFMotifs.has_within(20, 'GM.5.0.Rel', 'GM.5.0.bZIP', min_score_threshold=8).in_intron()"],
            RULE_YES,
        )

    def test_tf_motifs_in_intron_without_number_is_false_for_single_exon_gene(self):
        self.fx.write_single_exon_gene("p1", "g1")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "p1_locus 2 6 GM.5.0.Rel.0001 8.0 +",
            "p1_locus 8 12 GM.5.0.bZIP.0001 8.0 -",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron()
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["pass all"], RULE_FALSE)

    def test_tf_motifs_reports_missing_motif_when_no_pair_above_default_threshold(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "p1_locus 45 50 GM.5.0.Rel.0001 7.9 +",
            "p1_locus 55 60 GM.5.0.bZIP.0001 9.0 +",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron(2)
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["pass all"], RULE_FALSE)
        self.assertEqual(
            rows[0]["TFMotifs.has_within(20, 'GM.5.0.Rel', 'GM.5.0.bZIP', min_score_threshold=8).in_intron(2)"],
            "missing_GM.5.0.Rel",
        )

    def test_tf_motifs_returns_too_far_as_passing_when_hits_are_not_within_nearest_edge_distance(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "p1_locus 41 45 GM.5.0.Rel.0001 8.0 +",
            "p1_locus 67 70 GM.5.0.bZIP.0001 9.0 +",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron(2)
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["pass all"], RULE_TRUE)
        self.assertEqual(
            rows[0]["TFMotifs.has_within(20, 'GM.5.0.Rel', 'GM.5.0.bZIP', min_score_threshold=8).in_intron(2)"],
            RULE_TOO_FAR,
        )

    def test_tf_motifs_without_scope_can_match_in_any_intron_or_exon(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "p1_locus 12 17 GM.5.0.Rel.0001 8.0 +",
            "p1_locus 25 30 GM.5.0.bZIP.0001 8.0 -",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP")
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(
            rows[0]["TFMotifs.has_within(20, 'GM.5.0.Rel', 'GM.5.0.bZIP', min_score_threshold=8)"],
            RULE_YES,
        )

    def test_tf_motifs_without_scope_ignores_hits_outside_exons_and_introns(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "p1_locus 82 85 GM.5.0.Rel.0001 8.0 +",
            "p1_locus 86 90 GM.5.0.bZIP.0001 8.0 -",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP")
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(rows[0]["pass all"], RULE_FALSE)
        self.assertEqual(
            rows[0]["TFMotifs.has_within(20, 'GM.5.0.Rel', 'GM.5.0.bZIP', min_score_threshold=8)"],
            "missing_GM.5.0.Rel_and_GM.5.0.bZIP",
        )

    def test_tf_motifs_in_exon_matches_any_exon(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "p1_locus 32 35 GM.5.0.Rel.0001 8.0 +",
            "p1_locus 37 40 GM.5.0.bZIP.0001 8.0 -",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_exon()
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(
            rows[0]["TFMotifs.has_within(20, 'GM.5.0.Rel', 'GM.5.0.bZIP', min_score_threshold=8).in_exon()"],
            RULE_YES,
        )

    def test_tf_motifs_in_numbered_exon_uses_gene_direction(self):
        self.fx.write_three_exon_gene("p1", "g1", "-")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "p1_locus 52 55 GM.5.0.Rel.0001 8.0 -",
            "p1_locus 57 60 GM.5.0.bZIP.0001 8.0 +",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_exon(2)
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(
            rows[0]["TFMotifs.has_within(20, 'GM.5.0.Rel', 'GM.5.0.bZIP', min_score_threshold=8).in_exon(2)"],
            RULE_YES,
        )

    def test_tf_motifs_between_uses_locus_coordinates(self):
        self.fx.write_three_exon_gene("p1", "g1", "+")
        gimme_output = "\n".join([
            "sequence start end feature score strand",
            "p1_locus 82 85 GM.5.0.Rel.0001 8.0 +",
            "p1_locus 86 90 GM.5.0.bZIP.0001 8.0 -",
            "",
        ])

        with patch("sieve.rules.subprocess.run", return_value=CompletedProcess([], 0, stdout=gimme_output, stderr="")):
            with tempfile.TemporaryDirectory() as tmpd:
                rows = Rules(
                    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").between(90, 81)
                ).check([("p1", "g1")], os.path.join(tmpd, "rules.tsv"))

        self.assertEqual(
            rows[0]["TFMotifs.has_within(20, 'GM.5.0.Rel', 'GM.5.0.bZIP', min_score_threshold=8).between(90, 81)"],
            RULE_YES,
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

    def write_deeploc_csv(self, path, rows):
        fieldnames = [
            "Protein_ID",
            "Localizations",
            "Signals",
            "Membrane types",
            "Cytoplasm",
            "Endoplasmic reticulum",
        ]
        with open(path, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    def test_parse_deeploc_csv_maps_signals_and_scores(self):
        with tempfile.TemporaryDirectory() as tmpd:
            path = os.path.join(tmpd, "deeploc.csv")
            self.write_deeploc_csv(path, [
                {
                    "Protein_ID": "seq0",
                    "Localizations": "Endoplasmic reticulum",
                    "Signals": "Mitochondrial transit peptide",
                    "Membrane types": "Soluble",
                    "Cytoplasm": "0.7658",
                    "Endoplasmic reticulum": "0.2588",
                },
                {
                    "Protein_ID": "seq1",
                    "Localizations": "Cytoplasm",
                    "Signals": "",
                    "Membrane types": "Soluble",
                    "Cytoplasm": "0.5",
                    "Endoplasmic reticulum": "not-a-number",
                },
                {
                    "Protein_ID": "seq2",
                    "Localizations": "Extracellular",
                    "Signals": " Peroxisomal targeting signal | Signal peptide ",
                },
                {
                    "Protein_ID": "seq3",
                    "Localizations": "Extracellular",
                    "Signals": "Signal peptide|Mitochondrial transit peptide",
                },
            ])

            parsed = _parse_deeploc_csv(path)

        self.assertEqual(parsed["seq0"].prediction, "mTP")
        self.assertEqual(parsed["seq0"].localization, "Endoplasmic reticulum")
        self.assertEqual(parsed["seq0"].probability("mTP"), 1.0)
        self.assertEqual(parsed["seq0"].probability("SP"), 0.0)
        self.assertEqual(parsed["seq0"].probability("Cytoplasm"), 0.7658)
        self.assertEqual(parsed["seq0"].probability("Endoplasmic reticulum"), 0.2588)
        self.assertEqual(parsed["seq1"].prediction, "")
        self.assertEqual(parsed["seq1"].probability("mTP"), 0.0)
        self.assertEqual(parsed["seq1"].probability("SP"), 0.0)
        self.assertIsNone(parsed["seq1"].probability("Endoplasmic reticulum"))
        self.assertEqual(parsed["seq2"].prediction, "SP")
        self.assertEqual(parsed["seq2"].probability("mTP"), 0.0)
        self.assertEqual(parsed["seq2"].probability("SP"), 1.0)
        self.assertEqual(parsed["seq3"].prediction, "mTP")
        self.assertEqual(parsed["seq3"].probability("mTP"), 1.0)
        self.assertEqual(parsed["seq3"].probability("SP"), 1.0)

    def test_parse_targetp_output(self):
        parsed = _parse_targetp_output("\n".join([
            "# TargetP-2.0",
            "# ID Prediction noTP SP mTP CS Position",
            "seq0\tmTP\t0.1\t0.2\t0.7\tCS pos: 1-2. AA-AA. Pr: 0.1",
            "seq1 noTP 0.9 0.1 0.0",
        ]))

        self.assertEqual(parsed, {
            "seq0": _targetp_call("mTP", 0.1, 0.2, 0.7),
            "seq1": _targetp_call("noTP", 0.9, 0.1, 0.0),
        })

    def test_parse_gimme_scan_output(self):
        parsed = _parse_gimme_scan_output("\n".join([
            "sequence start end feature score strand",
            "p1_locus 45 50 GM.5.0.Rel.0001 8.0 +",
            "p1_locus 55 60 GM.5.0.bZIP.0001 9.0 -",
        ]))

        self.assertEqual(len(parsed["p1_locus"]), 2)
        self.assertEqual(parsed["p1_locus"][1].feature, "GM.5.0.bZIP.0001")
        self.assertEqual(parsed["p1_locus"][1].strand, "-")

    def test_edge_distance_uses_nearest_edges_and_overlap(self):
        self.assertEqual(_edge_distance(10, 20, 25, 30), 5)
        self.assertEqual(_edge_distance(25, 30, 10, 20), 5)
        self.assertEqual(_edge_distance(10, 20, 20, 30), 0)
        self.assertEqual(_edge_distance(10, 20, 15, 25), 0)


if __name__ == "__main__":
    unittest.main()
