import io
import os
import tempfile
import unittest
from unittest.mock import patch

from tangle.detected import DetectedTable

from sieve.protein import CuratedProtein, SEQUENCE_SOURCE_NCBI
from tests.fixtures import DefaultsFixture
from tests.scripts.helpers import load_script


class TestKoFindMatchesScript(unittest.TestCase):

    def setUp(self):
        CuratedProtein.clear_cache()
        self.fx = DefaultsFixture(self)
        self.repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    def tearDown(self):
        CuratedProtein.clear_cache()
        self.fx.cleanup()

    def detected_row(
        self,
        protein_accession,
        genome_accession,
        ko_accession,
        evalue,
        query_start=1,
        query_end=10,
        custom_metric_name=None,
        custom_metric_value=None,
        bitscore=100,
        bitscore_threshold=None,
    ):
        return dict(
            detection_type="sequence",
            detection_method="hmm",
            batch="b1",
            query_accession=protein_accession,
            query_database=genome_accession,
            query_type="protein",
            target_accession=ko_accession,
            target_database="KO",
            target_type="protein",
            query_start=query_start,
            query_end=query_end,
            target_start=1,
            target_end=10,
            evalue=evalue,
            bitscore=bitscore,
            bitscore_threshold=bitscore_threshold,
            custom_metric_name=custom_metric_name,
            custom_metric_value=custom_metric_value,
        )

    def add_ncbi_proteins(self, genome_accession, sequences):
        self.fx.write_manifest([
            {
                "sequence_accession": accession,
                "sequence_database": genome_accession,
                "sequence_type": "protein",
                "sequence_source": SEQUENCE_SOURCE_NCBI,
            }
            for accession in sequences
        ])
        self.fx.write_ncbi_proteins(genome_accession, sequences)

    def test_filters_by_ko_and_max_evalue(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_ko_assigned.tsv"), [
            self.detected_row("p1", "g1", "K04564", 1e-20,
                              custom_metric_name="evalue-rank", custom_metric_value=1),
            self.detected_row("p1", "g1", "K04564", 1e-20,
                              custom_metric_name="evalue-rank", custom_metric_value=1),
            self.detected_row("p2", "g1", "K04564", 1e-2,
                              custom_metric_name="evalue-rank", custom_metric_value=1),
            self.detected_row("p3", "g2", "K00001", 1e-50,
                              custom_metric_name="evalue-rank", custom_metric_value=1),
        ])

        self.assertEqual(script.find_matches("K04564", max_evalue=1e-10), [("p1", "g1")])
        stdout = io.StringIO()
        with patch("sys.stdout", stdout):
            script.main(["K04564", "--max-evalue", "1e-10"])
        self.assertEqual(stdout.getvalue(), "p1\tg1\n")

    def test_filters_by_inclusive_match_start_and_end_positions(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_ko_assigned.tsv"), [
            self.detected_row("both_boundary", "g1", "K04564", 1e-20, 20, 40),
            self.detected_row("start_after", "g1", "K04564", 1e-20, 21, 30),
            self.detected_row("end_after", "g1", "K04564", 1e-20, 10, 41),
        ])

        self.assertEqual(
            script.find_matches(
                "K04564",
                match_starts_before=20,
                match_ends_before=40,
            ),
            [("both_boundary", "g1")],
        )

    def test_position_filters_must_be_satisfied_by_the_same_hit(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_ko_assigned.tsv"), [
            self.detected_row("different_hits", "g1", "K04564", 1e-20, 10, 100),
            self.detected_row("different_hits", "g1", "K04564", 1e-20, 50, 60),
            self.detected_row("one_matching_hit", "g1", "K04564", 1e-20, 15, 55),
            self.detected_row("one_matching_hit", "g1", "K04564", 1e-20, 80, 100),
        ])

        self.assertEqual(
            script.find_matches(
                "K04564",
                match_starts_before=20,
                match_ends_before=60,
            ),
            [("one_matching_hit", "g1")],
        )

    def test_main_forwards_position_filters(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        with patch.object(script, "find_matches", return_value=[]) as find_matches:
            script.main([
                "K04564",
                "--match-starts-before", "20",
                "--match-ends-before", "40",
                "--max-evalue-rank", "3",
            ])

        find_matches.assert_called_once_with(
            "K04564",
            None,
            taxon=None,
            match_starts_before=20,
            match_ends_before=40,
            max_evalue_rank=3.0,
        )

    def test_filters_by_inclusive_evalue_rank(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_ko_assigned.tsv"), [
            self.detected_row("rank2", "g1", "K00001", 1e-30),
            self.detected_row("rank2", "g1", "K04564", 1e-20),
            self.detected_row("rank3", "g1", "K00001", 1e-40),
            self.detected_row("rank3", "g1", "K00002", 1e-30),
            self.detected_row("rank3", "g1", "K04564", 1e-20),
            self.detected_row("rank4", "g1", "K00001", 1e-50),
            self.detected_row("rank4", "g1", "K00002", 1e-40),
            self.detected_row("rank4", "g1", "K00003", 1e-30),
            self.detected_row("rank4", "g1", "K04564", 1e-20),
        ])

        self.assertEqual(
            script.find_matches("K04564", max_evalue_rank=3),
            [("rank2", "g1"), ("rank3", "g1")],
        )

    def test_ranks_only_hits_meeting_an_available_bitscore_threshold(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_ko_assigned.tsv"), [
            self.detected_row(
                "eligible", "g1", "K00001", 1e-30,
                bitscore=99, bitscore_threshold=100,
            ),
            self.detected_row(
                "eligible", "g1", "K04564", 1e-20,
                bitscore=100, bitscore_threshold=100,
            ),
            self.detected_row(
                "below", "g1", "K04564", 1e-40,
                bitscore=99, bitscore_threshold=100,
            ),
            self.detected_row("no_threshold", "g1", "K04564", 1e-10),
        ])

        self.assertEqual(
            script.find_matches("K04564", max_evalue_rank=1),
            [("eligible", "g1"), ("no_threshold", "g1")],
        )

    def test_equal_evalues_share_the_lowest_rank(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_ko_assigned.tsv"), [
            self.detected_row("tied", "g1", "K00001", 1e-20),
            self.detected_row("tied", "g1", "K04564", 1e-20),
        ])

        self.assertEqual(
            script.find_matches("K04564", max_evalue_rank=1),
            [("tied", "g1")],
        )

    def test_main_defaults_max_evalue_rank_to_one(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        with patch.object(script, "find_matches", return_value=[]) as find_matches:
            script.main(["K04564"])

        find_matches.assert_called_once_with(
            "K04564",
            None,
            taxon=None,
            match_starts_before=None,
            match_ends_before=None,
            max_evalue_rank=1,
        )

    def test_position_filter_applies_only_to_rows_passing_evalue_rank(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_ko_assigned.tsv"), [
            self.detected_row(
                "different_hits", "g1", "K04564", 1e-10, query_start=10,
            ),
            self.detected_row(
                "different_hits", "g1", "K04564", 1e-20, query_start=50,
            ),
            self.detected_row(
                "same_hit", "g1", "K04564", 1e-20, query_start=20,
            ),
        ])

        self.assertEqual(
            script.find_matches(
                "K04564",
                match_starts_before=20,
                max_evalue_rank=1,
            ),
            [("same_hit", "g1")],
        )

    def test_filters_matches_by_taxon_at_any_rank(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_ko_assigned.tsv"), [
            self.detected_row("p1", "g1", "K04564", 1e-20,
                              custom_metric_name="evalue-rank", custom_metric_value=1),
            self.detected_row("p2", "g2", "K04564", 1e-20,
                              custom_metric_name="evalue-rank", custom_metric_value=1),
            self.detected_row("p3", "g3", "K04564", 1e-20,
                              custom_metric_name="evalue-rank", custom_metric_value=1),
        ])
        self.fx.write_taxonomy_rows([
            {
                "Genome Accession": "g1",
                "Domain": "Eukaryota",
                "Phylum": "Cnidaria",
            },
            {
                "Genome Accession": "g2",
                "Domain": "Eukaryota",
                "Phylum": "Arthropoda",
            },
        ])

        self.assertEqual(script.find_matches("K04564", taxon="cNiDaRiA"), [("p1", "g1")])

        stdout = io.StringIO()
        with patch("sys.stdout", stdout):
            script.main(["K04564", "--taxon", "cnidaria"])
        self.assertEqual(stdout.getvalue(), "p1\tg1\n")

    def test_taxon_filter_returns_no_matches_without_taxonomy_file(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_ko_assigned.tsv"), [
            self.detected_row("p1", "g1", "K04564", 1e-20),
        ])

        self.assertEqual(script.find_matches("K04564", taxon="Cnidaria"), [])

    def test_output_writes_full_sequences_as_fasta(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        self.add_ncbi_proteins("g1", {"p1": "MSEQONE", "p2": "MSEQTWO"})
        with tempfile.TemporaryDirectory() as tmpd:
            output = os.path.join(tmpd, "matches.faa")
            with patch.object(script, "find_matches", return_value=[("p2", "g1"), ("p1", "g1")]):
                stdout = io.StringIO()
                with patch("sys.stdout", stdout):
                    script.main(["K04564", "-o", output])

            self.assertEqual(stdout.getvalue(), "")
            with open(output, encoding="utf-8") as f:
                self.assertEqual(f.read(), ">p2\nMSEQTWO\n>p1\nMSEQONE\n")

    def test_output_writes_empty_fasta_when_there_are_no_matches(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        with tempfile.TemporaryDirectory() as tmpd:
            output = os.path.join(tmpd, "matches.faa")
            with patch.object(script, "find_matches", return_value=[]):
                script.main(["K04564", "-o", output])

            with open(output, encoding="utf-8") as f:
                self.assertEqual(f.read(), "")

    def test_output_ignores_matches_missing_from_manifest(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        with tempfile.TemporaryDirectory() as tmpd:
            output = os.path.join(tmpd, "matches.faa")
            stderr = io.StringIO()
            with (
                patch.object(script, "find_matches", return_value=[("missing", "g1")]),
                patch("sys.stderr", stderr),
            ):
                script.main(["K04564", "-o", output])

            with open(output, encoding="utf-8") as f:
                self.assertEqual(f.read(), "")
            self.assertEqual(
                stderr.getvalue(),
                "Ignoring missing\tg1: Cannot find protein missing from g1 in manifest\n",
            )

    def test_output_does_not_ignore_other_sequence_errors(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        self.add_ncbi_proteins("g1", {"p1": "MSEQ"})
        os.remove(self.fx.genome_dir("g1") / "protein.faa")
        with tempfile.TemporaryDirectory() as tmpd:
            output = os.path.join(tmpd, "matches.faa")
            with patch.object(script, "find_matches", return_value=[("p1", "g1")]):
                with self.assertRaises(FileNotFoundError):
                    script.main(["K04564", "-o", output])
