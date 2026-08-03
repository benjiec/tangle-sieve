import io
import os
import tempfile
import unittest
from unittest.mock import patch

from tangle.detected import DetectedTable

from sieve.protein import CuratedProtein, SEQUENCE_SOURCE_NCBI
from tests.fixtures import DefaultsFixture
from tests.scripts.helpers import load_script


class TestPfamFindMatchesScript(unittest.TestCase):

    def setUp(self):
        CuratedProtein.clear_cache()
        self.fx = DefaultsFixture(self)
        self.repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.script = load_script(os.path.join(self.repo, "scripts", "pfam-find-matches.py"))

    def tearDown(self):
        CuratedProtein.clear_cache()
        self.fx.cleanup()

    def detected_row(self, protein_accession, genome_accession, pfam_accession, evalue):
        return dict(
            detection_type="sequence",
            detection_method="hmm",
            batch="b1",
            query_accession=protein_accession,
            query_database=genome_accession,
            query_type="protein",
            target_accession=pfam_accession,
            target_database="Pfam",
            target_type="protein",
            query_start=1,
            query_end=10,
            target_start=1,
            target_end=10,
            evalue=evalue,
        )

    def test_matches_unversioned_and_versioned_accessions(self):
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [
            self.detected_row("p1", "g1", "PF00504", 1e-20),
            self.detected_row("p2", "g1", "PF00504.27", 1e-30),
            self.detected_row("p3", "g1", "PF005040.1", 1e-40),
            self.detected_row("p4", "g1", "PF00504.28", 1e-2),
        ])

        self.assertEqual(
            self.script.find_matches("PF00504", max_evalue=1e-10),
            [("p1", "g1"), ("p2", "g1")],
        )

    def test_filters_by_taxon_and_prints_tsv_without_output(self):
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [
            self.detected_row("p1", "g1", "PF00504.27", 1e-20),
            self.detected_row("p2", "g2", "PF00504.27", 1e-20),
        ])
        self.fx.write_taxonomy_rows([
            {"Genome Accession": "g1", "Domain": "Eukaryota", "Phylum": "Alveolata"},
            {"Genome Accession": "g2", "Domain": "Eukaryota", "Phylum": "Cnidaria"},
        ])
        stdout = io.StringIO()
        with patch("sys.stdout", stdout):
            self.script.main(["PF00504", "--taxon", "alveolata"])
        self.assertEqual(stdout.getvalue(), "p1\tg1\n")

    def test_output_writes_fasta_and_ignores_missing_manifest_entries(self):
        self.fx.write_manifest([{
            "sequence_accession": "p1",
            "sequence_database": "g1",
            "sequence_type": "protein",
            "sequence_source": SEQUENCE_SOURCE_NCBI,
        }])
        self.fx.write_ncbi_proteins("g1", {"p1": "MSEQ"})
        with tempfile.TemporaryDirectory() as tmpd:
            output = os.path.join(tmpd, "matches.faa")
            stderr = io.StringIO()
            with (
                patch.object(self.script, "find_matches", return_value=[("p1", "g1"), ("missing", "g1")]),
                patch("sys.stderr", stderr),
            ):
                self.script.main(["PF00504", "-o", output])

            with open(output, encoding="utf-8") as f:
                self.assertEqual(f.read(), ">p1\nMSEQ\n")
            self.assertIn("Ignoring missing\tg1:", stderr.getvalue())
