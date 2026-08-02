import io
import os
import unittest
from unittest.mock import patch

from tangle.detected import DetectedTable

from sieve.protein import CuratedProtein
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

    def detected_row(self, protein_accession, genome_accession, ko_accession, evalue):
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
            query_start=1,
            query_end=10,
            target_start=1,
            target_end=10,
            evalue=evalue,
        )

    def test_filters_by_ko_and_max_evalue(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_ko_assigned.tsv"), [
            self.detected_row("p1", "g1", "K04564", 1e-20),
            self.detected_row("p1", "g1", "K04564", 1e-20),
            self.detected_row("p2", "g1", "K04564", 1e-2),
            self.detected_row("p3", "g2", "K00001", 1e-50),
        ])

        self.assertEqual(script.find_matches("K04564", max_evalue=1e-10), [("p1", "g1")])
        stdout = io.StringIO()
        with patch("sys.stdout", stdout):
            script.main(["K04564", "--max-evalue", "1e-10"])
        self.assertEqual(stdout.getvalue(), "p1\tg1\n")

    def test_filters_matches_by_taxon_at_any_rank(self):
        script = load_script(os.path.join(self.repo, "scripts", "ko-find-matches.py"))
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_ko_assigned.tsv"), [
            self.detected_row("p1", "g1", "K04564", 1e-20),
            self.detected_row("p2", "g2", "K04564", 1e-20),
            self.detected_row("p3", "g3", "K04564", 1e-20),
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
