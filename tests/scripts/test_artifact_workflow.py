import csv
import io
import os
import sys
import tempfile
import unittest
from subprocess import CompletedProcess
from unittest.mock import patch

from sieve.artifacts import input_fasta, sequences_fasta, sequences_tsv
from sieve.protein import CuratedProtein, SEQUENCE_SOURCE_NCBI
from tangle.sequence import read_fasta_as_dict, write_fasta_from_dict
from tangle.manifest import ManifestTable

from tests.fixtures import DefaultsFixture
from tests.scripts.helpers import load_script


def domtblout_line(sequence, model_name, model_accession, full_score, domain_score):
    return " ".join([
        sequence, "-", "100", model_name, model_accession, "200", "1e-20",
        str(full_score), "0.0", "1", "1", "1e-10", "1e-10", str(domain_score),
        "0.0", "3", "50", "7", "54", "6", "55", "0.98",
    ])


class TestArtifactWorkflow(unittest.TestCase):

    def setUp(self):
        self.repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.build_fasta = load_script(os.path.join(self.repo, "scripts", "build-fasta-artifacts.py"))
        self.check = load_script(os.path.join(self.repo, "scripts", "check-artifacts-by-rules.py"))

    def write_rule(self, directory, module, expression, imports="Leader, Rules"):
        path = os.path.join(directory, module + ".py")
        with open(path, "w", encoding="utf-8") as f:
            f.write(f"from sieve.rules import {imports}\n")
            f.write(f"rule = Rules({expression})\n")
        return module + ".rule"

    def read_tsv(self, path):
        with open(path, "r", encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    def test_fasta_builder_is_incremental_and_unions_candidates(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            write_fasta_from_dict({"p1": "MMT"}, fasta)
            first = self.write_rule(tmpd, "first_rules", "Leader().betweenAA(1, 1)")
            second = self.write_rule(tmpd, "second_rules", "Leader().betweenAA(2, 2)")
            sys.path.insert(0, tmpd)
            try:
                self.build_fasta.main([
                    "--fasta", fasta, "--rule", first, "--artifacts-dir", artifacts,
                ])
                self.build_fasta.main([
                    "--fasta", fasta, "--rule", second, "--artifacts-dir", artifacts,
                ])
                self.build_fasta.main([
                    "--fasta", fasta, "--rule", second, "--artifacts-dir", artifacts,
                ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("first_rules", None)
                sys.modules.pop("second_rules", None)

            self.assertEqual(read_fasta_as_dict(input_fasta(artifacts)), {"p1": "MMT"})
            self.assertEqual(read_fasta_as_dict(sequences_fasta(artifacts)), {
                "p1_with_leader_1_M": "MMT",
                "p1_with_leader_2_M": "MT",
            })
            self.assertEqual(
                [row["start aa 1b"] for row in self.read_tsv(sequences_tsv(artifacts))],
                ["1", "2"],
            )

    def test_fasta_builder_rejects_conflicting_input_accession(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            rule = self.write_rule(tmpd, "rules", "Leader().betweenAA(1, 2)")
            sys.path.insert(0, tmpd)
            try:
                write_fasta_from_dict({"p1": "MMT"}, fasta)
                self.build_fasta.main(["--fasta", fasta, "--rule", rule, "--artifacts-dir", artifacts])
                write_fasta_from_dict({"p1": "MGT"}, fasta)
                with self.assertRaisesRegex(ValueError, "Conflicting input sequence"):
                    self.build_fasta.main(["--fasta", fasta, "--rule", rule, "--artifacts-dir", artifacts])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("rules", None)

    def test_zero_candidate_original_keeps_metadata_for_stage_b(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            write_fasta_from_dict({"p1": "AAAA"}, fasta)
            rule = self.write_rule(tmpd, "rules", "Leader().betweenAA(1, 2)")
            sys.path.insert(0, tmpd)
            try:
                self.build_fasta.main(["--fasta", fasta, "--rule", rule, "--artifacts-dir", artifacts])
                self.check.main(["--rule", rule, "--artifacts-dir", artifacts])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("rules", None)

            manifest = self.read_tsv(sequences_tsv(artifacts))
            self.assertEqual(len(manifest), 1)
            self.assertEqual(manifest[0]["protein accession"], "p1")
            self.assertEqual(manifest[0]["sequence accession"], "")
            self.assertEqual(self.read_tsv(os.path.join(artifacts, "rule-results.tsv")), [])

    def test_stage_b_uses_recorded_candidates_without_leader_discovery(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            write_fasta_from_dict({"p1": "MMT"}, fasta)
            rule = self.write_rule(tmpd, "rules", "Leader().betweenAA(1, 2)")
            sys.path.insert(0, tmpd)
            try:
                self.build_fasta.main(["--fasta", fasta, "--rule", rule, "--artifacts-dir", artifacts])
                with patch(
                    "sieve.artifact_protein.ArtifactProtein.sequences_with_leader",
                    side_effect=AssertionError("leader discovery called"),
                ):
                    self.check.main(["--rule", rule, "--artifacts-dir", artifacts])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("rules", None)

            rows = self.read_tsv(os.path.join(artifacts, "rule-results.tsv"))
            self.assertEqual(
                [row["sequence accession"] for row in rows],
                ["p1_with_leader_1_M", "p1_with_leader_2_M"],
            )

    def test_stage_b_batch_aligns_recorded_candidates(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            profile = os.path.join(tmpd, "profile.hmm")
            open(profile, "w").close()
            write_fasta_from_dict({"p1": "MMT"}, fasta)
            build_rule = self.write_rule(tmpd, "build_rules", "Leader().betweenAA(1, 2)")
            check_rule = self.write_rule(
                tmpd,
                "check_rules",
                "HMMAlignment('profile.hmm').is_at('M', 1)",
                imports="HMMAlignment, Rules",
            )
            captured = {}

            def align(_profile, sequences):
                captured.update(sequences)
                return {
                    accession: unittest.mock.Mock(
                        aa_at_hmm_pos_1b=lambda position: (1, "M") if position == 1 else None,
                    )
                    for accession in sequences
                }

            sys.path.insert(0, tmpd)
            try:
                self.build_fasta.main([
                    "--fasta", fasta, "--rule", build_rule, "--artifacts-dir", artifacts,
                ])
                with patch("sieve.rules.hmm_align_sequences", side_effect=align):
                    self.check.main([
                        "--rule", check_rule,
                        "--artifacts-dir", artifacts,
                        "--hmm-dir", tmpd,
                    ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("build_rules", None)
                sys.modules.pop("check_rules", None)

            self.assertEqual(set(captured.values()), {"MMT", "MT"})
            self.assertTrue(all(row["pass all"] == "true" for row in self.read_tsv(
                os.path.join(artifacts, "rule-results.tsv")
            )))

    def test_stage_b_regenerates_pfam_and_ko_with_distinct_threshold_modes(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            pfam_hmm = os.path.join(tmpd, "pfam.hmm")
            ko_hmm = os.path.join(tmpd, "ko.hmm")
            thresholds = os.path.join(tmpd, "thresholds.tsv")
            for path in (pfam_hmm, ko_hmm):
                with open(path, "w", encoding="utf-8"):
                    pass
            with open(thresholds, "w", encoding="utf-8") as f:
                f.write("model threshold score_type\nK04564 100 full\n")
            write_fasta_from_dict({"p1": "MAAA"}, fasta)
            build_rule = self.write_rule(
                tmpd, "build_rules", "Pfam.matches('PF00081')", imports="Pfam, Rules",
            )
            check_rule = self.write_rule(
                tmpd,
                "check_rules",
                "Pfam.matches('PF00081') & KO.matches('K04564')",
                imports="KO, Pfam, Rules",
            )
            commands = []

            def run(cmd, check, capture_output, text):
                commands.append(cmd)
                domtblout = cmd[cmd.index("--domtblout") + 1]
                if "--cut_ga" in cmd:
                    line = domtblout_line("p1", "PFAM", "PF00081.1", 120, 110)
                else:
                    line = domtblout_line("p1", "K04564", "-", 120, 110)
                with open(domtblout, "w", encoding="utf-8") as f:
                    f.write(line + "\n")
                return CompletedProcess(cmd, 0, stdout="", stderr="")

            sys.path.insert(0, tmpd)
            try:
                self.build_fasta.main([
                    "--fasta", fasta, "--rule", build_rule, "--artifacts-dir", artifacts,
                ])
                with patch("subprocess.run", side_effect=run):
                    self.check.main([
                        "--rule", check_rule,
                        "--artifacts-dir", artifacts,
                        "--pfam-hmm", pfam_hmm,
                        "--ko-hmm", ko_hmm,
                        "--ko-thresholds", thresholds,
                    ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("build_rules", None)
                sys.modules.pop("check_rules", None)

            self.assertIn("--cut_ga", commands[0])
            self.assertNotIn("--cut_ga", commands[1])
            rows = self.read_tsv(os.path.join(artifacts, "rule-results.tsv"))
            self.assertEqual(rows[0]["pass all"], "true")


class TestCuratedArtifactBuilder(unittest.TestCase):

    def setUp(self):
        CuratedProtein.clear_cache()
        self.fx = DefaultsFixture(self)
        repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.build = load_script(os.path.join(repo, "scripts", "build-protein-artifacts.py"))

    def tearDown(self):
        CuratedProtein.clear_cache()
        self.fx.cleanup()

    def read_tsv(self, path):
        with open(path, "r", encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    def test_accepts_accession_genome_pipeline_and_writes_originals(self):
        ManifestTable.write_tsv(str(self.fx.area_genomics / "sequences.tsv"), [
            dict(
                sequence_accession="p1",
                sequence_database="g1",
                sequence_type="protein",
                sequence_source=SEQUENCE_SOURCE_NCBI,
            ),
        ])
        self.fx.write_ncbi_proteins("g1", {"p1": "MMT"})
        with tempfile.TemporaryDirectory() as tmpd:
            rule_path = os.path.join(tmpd, "rules.py")
            with open(rule_path, "w", encoding="utf-8") as f:
                f.write("from sieve.rules import Leader, Rules\n")
                f.write("rule = Rules(Leader().betweenAA(1, 2))\n")
            artifacts = os.path.join(tmpd, "artifacts")
            sys.path.insert(0, tmpd)
            try:
                with patch("sys.stdin", io.StringIO("p1\tg1\n")):
                    self.build.main(["--rule", "rules.rule", "--artifacts-dir", artifacts])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("rules", None)

            self.assertEqual(read_fasta_as_dict(input_fasta(artifacts)), {"p1": "MMT"})
            self.assertEqual(set(read_fasta_as_dict(sequences_fasta(artifacts))), {
                "p1_with_leader_1_M", "p1_with_leader_2_M",
            })
            self.assertTrue(os.path.exists(os.path.join(artifacts, "genomic_loci.tsv")))

    def test_same_accession_in_later_genome_is_ignored_when_sequence_matches(self):
        ManifestTable.write_tsv(str(self.fx.area_genomics / "sequences.tsv"), [
            dict(sequence_accession="p1", sequence_database="g1", sequence_type="protein", sequence_source=SEQUENCE_SOURCE_NCBI),
            dict(sequence_accession="p1", sequence_database="g2", sequence_type="protein", sequence_source=SEQUENCE_SOURCE_NCBI),
        ])
        self.fx.write_ncbi_proteins("g1", {"p1": "MMT"})
        self.fx.write_ncbi_proteins("g2", {"p1": "MMT"})
        with tempfile.TemporaryDirectory() as tmpd:
            rule_path = os.path.join(tmpd, "rules.py")
            with open(rule_path, "w", encoding="utf-8") as f:
                f.write("from sieve.rules import Leader, Rules\n")
                f.write("rule = Rules(Leader().betweenAA(1, 2))\n")
            artifacts = os.path.join(tmpd, "artifacts")
            sys.path.insert(0, tmpd)
            try:
                with (
                    patch("sys.stdin", io.StringIO("p1\tg1\np1\tg2\n")),
                    patch("sys.stderr", new_callable=io.StringIO) as stderr,
                ):
                    self.build.main(["--rule", "rules.rule", "--artifacts-dir", artifacts])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("rules", None)

            self.assertIn("Ignoring p1 from genome g2", stderr.getvalue())
            self.assertEqual(
                {row["genome accession"] for row in self.read_tsv(sequences_tsv(artifacts))},
                {"g1"},
            )
            self.assertEqual(
                {row["genome accession"] for row in self.read_tsv(os.path.join(artifacts, "genomic_loci.tsv"))},
                {"g1"},
            )

    def test_same_accession_in_later_genome_errors_when_sequence_differs(self):
        ManifestTable.write_tsv(str(self.fx.area_genomics / "sequences.tsv"), [
            dict(sequence_accession="p1", sequence_database="g1", sequence_type="protein", sequence_source=SEQUENCE_SOURCE_NCBI),
            dict(sequence_accession="p1", sequence_database="g2", sequence_type="protein", sequence_source=SEQUENCE_SOURCE_NCBI),
        ])
        self.fx.write_ncbi_proteins("g1", {"p1": "MMT"})
        self.fx.write_ncbi_proteins("g2", {"p1": "MGT"})
        with tempfile.TemporaryDirectory() as tmpd:
            rule_path = os.path.join(tmpd, "rules.py")
            with open(rule_path, "w", encoding="utf-8") as f:
                f.write("from sieve.rules import Leader, Rules\n")
                f.write("rule = Rules(Leader().betweenAA(1, 2))\n")
            artifacts = os.path.join(tmpd, "artifacts")
            sys.path.insert(0, tmpd)
            try:
                with patch("sys.stdin", io.StringIO("p1\tg1\np1\tg2\n")):
                    with self.assertRaisesRegex(ValueError, "Conflicting input sequence for accession p1"):
                        self.build.main(["--rule", "rules.rule", "--artifacts-dir", artifacts])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("rules", None)


if __name__ == "__main__":
    unittest.main()
