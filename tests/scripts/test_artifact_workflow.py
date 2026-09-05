import csv
import io
import os
import sys
import tempfile
import unittest
from subprocess import CompletedProcess
from unittest.mock import patch

from sieve.artifacts import input_fasta, sequences_fasta, sequences_tsv
from sieve.protein import CuratedProtein, LeaderSequenceCandidate, SEQUENCE_SOURCE_NCBI
from tangle.sequence import read_fasta_as_dict, write_fasta_from_dict
from tangle.manifest import ManifestTable
from tangle.detected import DetectedTable

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

    def test_fasta_builder_requires_accession_and_genome_description(self):
        invalid_identifiers = [
            "p1",
            "p1|description|extra",
            "|description",
            "p1|",
            "p1|description with spaces",
        ]
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            for identifier in invalid_identifiers:
                with self.subTest(identifier=identifier):
                    write_fasta_from_dict({identifier: "MAAA"}, fasta)
                    with self.assertRaisesRegex(ValueError, "FASTA identifier"):
                        self.build_fasta.read_fasta_sequences([fasta])

    def test_fasta_builder_combines_multiple_files_and_rejects_duplicate_identifiers(self):
        with tempfile.TemporaryDirectory() as tmpd:
            first_fasta = os.path.join(tmpd, "first.faa")
            second_fasta = os.path.join(tmpd, "second.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            profile = os.path.join(tmpd, "pfam.hmm")
            open(profile, "w").close()
            write_fasta_from_dict({"p1|first_genome": "MAAA"}, first_fasta)
            write_fasta_from_dict({"p2|second_genome": "MCCC"}, second_fasta)
            rule = self.write_rule(tmpd, "rules", "Leader().betweenAA(1, 1)")
            searched_sequences = {}

            def run_hmmsearch(_hmm, fasta, domtblout):
                searched_sequences.update(read_fasta_as_dict(fasta))
                open(domtblout, "w").close()

            sys.path.insert(0, tmpd)
            try:
                with patch.object(self.build_fasta, "run_hmmsearch", side_effect=run_hmmsearch):
                    self.build_fasta.main([
                        "--fasta", first_fasta,
                        "--fasta", second_fasta,
                        "--rule", rule,
                        "--artifacts-dir", artifacts,
                        "--pfam-hmm", profile,
                    ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("rules", None)

            self.assertEqual(read_fasta_as_dict(input_fasta(artifacts)), {
                "p1|first_genome": "MAAA",
                "p2|second_genome": "MCCC",
            })
            self.assertEqual(searched_sequences, {
                "p1|first_genome": "MAAA",
                "p2|second_genome": "MCCC",
            })
            self.assertEqual(read_fasta_as_dict(sequences_fasta(artifacts)), {
                "p1_with_leader_1_M|first_genome": "MAAA",
                "p2_with_leader_1_M|second_genome": "MCCC",
            })

            write_fasta_from_dict({"p1|first_genome": "MAAA"}, second_fasta)
            with self.assertRaisesRegex(ValueError, "Duplicate FASTA identifier: p1\\|first_genome"):
                self.build_fasta.read_fasta_sequences([first_fasta, second_fasta])

    def test_fasta_builder_is_incremental_and_unions_candidates(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            write_fasta_from_dict({"p1|genome_description": "MMT"}, fasta)
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

            self.assertEqual(read_fasta_as_dict(input_fasta(artifacts)), {"p1|genome_description": "MMT"})
            self.assertEqual(read_fasta_as_dict(sequences_fasta(artifacts)), {
                "p1_with_leader_1_M|genome_description": "MMT",
                "p1_with_leader_2_M|genome_description": "MT",
            })
            self.assertEqual(
                [row["start aa 1b"] for row in self.read_tsv(sequences_tsv(artifacts))],
                ["1", "2"],
            )
            self.assertEqual(
                [row["protein start aa 1b"] for row in self.read_tsv(sequences_tsv(artifacts))],
                ["1", "2"],
            )

    def test_fasta_builder_rejects_conflicting_input_accession(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            rule = self.write_rule(tmpd, "rules", "Leader().betweenAA(1, 2)")
            sys.path.insert(0, tmpd)
            try:
                write_fasta_from_dict({"p1|genome_description": "MMT"}, fasta)
                self.build_fasta.main(["--fasta", fasta, "--rule", rule, "--artifacts-dir", artifacts])
                write_fasta_from_dict({"p1|genome_description": "MGT"}, fasta)
                with self.assertRaisesRegex(ValueError, "Conflicting input sequence"):
                    self.build_fasta.main(["--fasta", fasta, "--rule", rule, "--artifacts-dir", artifacts])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("rules", None)

    def test_zero_candidate_original_keeps_metadata_for_stage_b(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            write_fasta_from_dict({"p1|genome_description": "AAAA"}, fasta)
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
            self.assertEqual(manifest[0]["protein accession"], "p1|genome_description")
            self.assertEqual(manifest[0]["sequence accession"], "")
            self.assertEqual(
                read_fasta_as_dict(input_fasta(artifacts)),
                {"p1|genome_description": "AAAA"},
            )
            self.assertEqual(read_fasta_as_dict(sequences_fasta(artifacts)), {})
            self.assertEqual(self.read_tsv(os.path.join(artifacts, "rule-results.tsv")), [])

    def test_stage_b_uses_recorded_candidates_without_leader_discovery(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            write_fasta_from_dict({"p1|genome_description": "MMT"}, fasta)
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
                ["p1_with_leader_1_M|genome_description", "p1_with_leader_2_M|genome_description"],
            )

    def test_check_artifacts_restores_deeploc_truncated_identifiers(self):
        with tempfile.TemporaryDirectory() as tmpd:
            source = os.path.join(tmpd, "deeploc.csv")
            destination = os.path.join(tmpd, "normalized.csv")
            with open(source, "w", encoding="utf-8", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=["Protein_ID", "Localizations", "Signals"])
                writer.writeheader()
                writer.writerows([
                    {"Protein_ID": "p1_with_leader_1_M", "Localizations": "Cytoplasm", "Signals": ""},
                    {"Protein_ID": "p2_with_leader_1_M_second_genome", "Localizations": "Mitochondrion", "Signals": "Mitochondrial transit peptide"},
                    {"Protein_ID": "unknown", "Localizations": "Nucleus", "Signals": ""},
                ])
            candidates = {
                ("p1|genome", ""): [LeaderSequenceCandidate(
                    accession="p1_with_leader_1_M|genome",
                    start_label="1",
                    start_aa_1b=1,
                    sequence="MAAA",
                )],
                ("p2|second_genome", ""): [LeaderSequenceCandidate(
                    accession="p2_with_leader_1_M|second_genome",
                    start_label="1",
                    start_aa_1b=1,
                    sequence="MCCC",
                )],
            }

            self.check.normalize_deeploc_ids(source, destination, candidates)

            with open(destination, "r", encoding="utf-8", newline="") as f:
                normalized = list(csv.DictReader(f))
            self.assertEqual(normalized[0]["Protein_ID"], "p1_with_leader_1_M|genome")
            self.assertEqual(normalized[0]["Localizations"], "Cytoplasm")
            self.assertEqual(normalized[1]["Protein_ID"], "p2_with_leader_1_M|second_genome")
            self.assertEqual(normalized[1]["Signals"], "Mitochondrial transit peptide")
            self.assertEqual(normalized[2]["Protein_ID"], "unknown")

    def test_check_artifacts_rejects_ambiguous_deeploc_truncated_identifier(self):
        with tempfile.TemporaryDirectory() as tmpd:
            source = os.path.join(tmpd, "deeploc.csv")
            destination = os.path.join(tmpd, "normalized.csv")
            with open(source, "w", encoding="utf-8", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=["Protein_ID", "Localizations", "Signals"])
                writer.writeheader()
            candidates = {
                ("p1|first", ""): [LeaderSequenceCandidate("candidate|first", "", 1, "MAAA")],
                ("p1|second", ""): [LeaderSequenceCandidate("candidate|second", "", 1, "MCCC")],
            }

            with self.assertRaisesRegex(ValueError, "ambiguous.*candidate"):
                self.check.normalize_deeploc_ids(source, destination, candidates)

    def test_stage_b_batch_aligns_recorded_candidates(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            profile = os.path.join(tmpd, "profile.hmm")
            open(profile, "w").close()
            write_fasta_from_dict({"p1|genome_description": "MMT"}, fasta)
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
            write_fasta_from_dict({"p1|genome_description": "MAAA"}, fasta)
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
                    line = domtblout_line("p1|genome_description", "PFAM", "PF00081.1", 120, 110)
                else:
                    line = domtblout_line("p1|genome_description", "K04564", "-", 120, 110)
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

    def test_fasta_builder_uses_unthresholded_ko_search_for_cterm_bound(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            ko_hmm = os.path.join(tmpd, "ko.hmm")
            with open(ko_hmm, "w", encoding="utf-8"):
                pass
            sequence = "M" + "A" * 99
            write_fasta_from_dict({"p1|genome_description": sequence}, fasta)
            unbounded_rule = self.write_rule(
                tmpd, "unbounded_rules", "KO.matches('K04564')", imports="KO, Rules",
            )
            bounded_rule = self.write_rule(
                tmpd,
                "bounded_rules",
                "KO.matches('K04564', bound_cterm=True)",
                imports="KO, Rules",
            )
            commands = []

            def run(cmd, check, capture_output, text):
                commands.append(cmd)
                domtblout = cmd[cmd.index("--domtblout") + 1]
                with open(domtblout, "w", encoding="utf-8") as f:
                    f.write(domtblout_line("p1|genome_description", "K04564", "-", 80, 70) + "\n")
                return CompletedProcess(cmd, 0, stdout="", stderr="")

            sys.path.insert(0, tmpd)
            try:
                self.build_fasta.main([
                    "--fasta", fasta,
                    "--rule", unbounded_rule,
                    "--artifacts-dir", artifacts,
                ])
                with patch("subprocess.run", side_effect=run):
                    self.build_fasta.main([
                        "--fasta", fasta,
                        "--rule", bounded_rule,
                        "--artifacts-dir", artifacts,
                        "--ko-hmm", ko_hmm,
                    ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("unbounded_rules", None)
                sys.modules.pop("bounded_rules", None)

            self.assertNotIn("--cut_ga", commands[0])
            candidates = read_fasta_as_dict(sequences_fasta(artifacts))
            self.assertEqual(candidates["p1|genome_description"], sequence)
            self.assertEqual(candidates["p1_to_K04564_54|genome_description"], sequence[:54])
            rows = self.read_tsv(sequences_tsv(artifacts))
            bounded = next(row for row in rows if row["sequence accession"] == "p1_to_K04564_54|genome_description")
            self.assertEqual(bounded["end aa 1b"], "54")

    def test_fasta_builder_requires_ko_hmm_for_cterm_bound(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            write_fasta_from_dict({"p1|genome_description": "MAAA"}, fasta)
            rule = self.write_rule(
                tmpd,
                "bounded_rules",
                "KO.matches('K04564', bound_cterm=True)",
                imports="KO, Rules",
            )
            sys.path.insert(0, tmpd)
            try:
                with self.assertRaisesRegex(ValueError, "--ko-hmm is required"):
                    self.build_fasta.main([
                        "--fasta", fasta, "--rule", rule, "--artifacts-dir", artifacts,
                    ])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("bounded_rules", None)

    def test_pfam_anchored_sequence_regex_requires_pfam_database(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            artifacts = os.path.join(tmpd, "artifacts")
            write_fasta_from_dict({"p1|genome_description": "MOTIF"}, fasta)
            rule = self.write_rule(
                tmpd,
                "anchored_rules",
                "Sequence.matches_regex('MOTIF').relativeToPfam('PF00001', 1, 5)",
                imports="Rules, Sequence",
            )
            sys.path.insert(0, tmpd)
            try:
                self.build_fasta.main([
                    "--fasta", fasta, "--rule", rule, "--artifacts-dir", artifacts,
                ])
                with self.assertRaisesRegex(ValueError, "--pfam-hmm is required"):
                    self.check.main(["--rule", rule, "--artifacts-dir", artifacts])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("anchored_rules", None)


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

    def test_zero_candidate_original_is_only_written_to_input_fasta(self):
        ManifestTable.write_tsv(str(self.fx.area_genomics / "sequences.tsv"), [
            dict(
                sequence_accession="p1",
                sequence_database="g1",
                sequence_type="protein",
                sequence_source=SEQUENCE_SOURCE_NCBI,
            ),
        ])
        self.fx.write_ncbi_proteins("g1", {"p1": "AAAA"})
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

            self.assertEqual(read_fasta_as_dict(input_fasta(artifacts)), {"p1": "AAAA"})
            self.assertEqual(read_fasta_as_dict(sequences_fasta(artifacts)), {})
            manifest = self.read_tsv(sequences_tsv(artifacts))
            self.assertEqual(len(manifest), 1)
            self.assertEqual(manifest[0]["protein accession"], "p1")
            self.assertEqual(manifest[0]["sequence accession"], "")

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

    def test_curated_builder_bounds_candidates_from_curated_ko_hits(self):
        ManifestTable.write_tsv(str(self.fx.area_genomics / "sequences.tsv"), [
            dict(sequence_accession="p1", sequence_database="g1", sequence_type="protein", sequence_source=SEQUENCE_SOURCE_NCBI),
        ])
        self.fx.write_ncbi_proteins("g1", {"p1": "MABCDEFGHIJ"})
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_ko_assigned.tsv"), [
            dict(
                detection_type="sequence", detection_method="hmm", batch="b1",
                query_accession="p1", query_database="g1", query_type="protein",
                target_accession="K04564", target_database="KO", target_type="protein",
                query_start=1, query_end=6, target_start=1, target_end=6,
            ),
        ])
        with tempfile.TemporaryDirectory() as tmpd:
            rule_path = os.path.join(tmpd, "rules.py")
            with open(rule_path, "w", encoding="utf-8") as f:
                f.write("from sieve.rules import KO, Rules\n")
                f.write("rule = Rules(KO.matches('K04564', bound_cterm=True))\n")
            artifacts = os.path.join(tmpd, "artifacts")
            sys.path.insert(0, tmpd)
            try:
                with patch("sys.stdin", io.StringIO("p1\tg1\n")):
                    self.build.main(["--rule", "rules.rule", "--artifacts-dir", artifacts])
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("rules", None)

            self.assertEqual(read_fasta_as_dict(sequences_fasta(artifacts)), {
                "p1_to_K04564_6": "MABCDE",
            })

if __name__ == "__main__":
    unittest.main()
