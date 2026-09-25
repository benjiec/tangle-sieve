import csv
import gzip
import io
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
import unittest
from unittest.mock import patch

from tests.fixtures import DefaultsFixture
from tests.scripts.helpers import load_script


class TestProteinNames(unittest.TestCase):
    def setUp(self):
        self.script = load_script(Path(__file__).resolve().parents[2] / "scripts" / "protein-names.py")
        self.fixture = DefaultsFixture(self)
        self.addCleanup(self.fixture.cleanup)
        self.addCleanup(self.script.CuratedProtein.clear_cache)
        self.input = self.fixture.root / "accessions.txt"

    def manifest(self, entries):
        self.fixture.write_manifest([
            dict(sequence_accession=accession, sequence_database=genome,
                 sequence_type=kind, sequence_source="ncbi", sequence_length=3)
            for accession, genome, kind in entries
        ])

    def fasta(self, genome, headers):
        path = self.fixture.genome_dir(genome) / "protein.faa"
        path.write_text("".join(f">{header}\nMA\nK\n" for header in headers))

    def run_script(self, content):
        self.input.write_text(content)
        output, errors = io.StringIO(), io.StringIO()
        with redirect_stdout(output), redirect_stderr(errors):
            status = self.script.main([str(self.input)])
        return status, output.getvalue(), errors.getvalue()

    def test_order_duplicates_multiple_genomes_and_one_scan_per_genome(self):
        self.manifest([("P1.1", "G1", "protein"), ("P2", "G1", "protein"),
                       ("P1.1", "G2", "protein"), ("P1.1", "G1", "protein"),
                       ("P2", "IGNORED", "contig")])
        self.fasta("G1", ["OTHER irrelevant", "P1.1 kinase [Species one]",
                          "P2 uncharacterized protein"])
        self.fasta("G2", ["P1.1 kinase [Species two]"])
        with patch.object(self.script, "read_names", wraps=self.script.read_names) as read:
            result = self.run_script(" P2 \n\nP1.1\nP2\n")
        self.assertEqual(read.call_count, 2)
        self.assertEqual(result, (0, "accession_id\tprotein_name\nP2\tuncharacterized protein\n"
                                 "P1.1\tkinase\nP2\tuncharacterized protein\n", ""))

    def test_headers_and_tsv_escaping(self):
        for description, expected in [
            ("enzyme [internal] subunit [Species]", "enzyme [internal] subunit"),
            ("enzyme [internal] subunit", "enzyme [internal] subunit"),
            ('enzyme\t"alpha" [Species]', 'enzyme\t"alpha"'),
            ("enzyme [unclosed", "enzyme [unclosed"),
        ]:
            with self.subTest(description=description):
                self.fasta("G1", [f"P1 {description}"])
                self.assertEqual(self.script.read_names("G1", {"P1"}), {"P1": expected})
        self.manifest([("P1", "G1", "protein")])
        self.fasta("G1", ['P1 enzyme\t"alpha" [Species]'])
        status, output, errors = self.run_script("P1\n")
        self.assertEqual((status, errors), (0, ""))
        self.assertEqual(list(csv.reader(io.StringIO(output), delimiter="\t"))[1],
                         ["P1", 'enzyme\t"alpha"'])

    def test_empty_input(self):
        self.assertEqual(self.run_script("\n  \n"), (0, "accession_id\tprotein_name\n", ""))

    def test_gzipped_fasta(self):
        self.manifest([("P1", "G1", "protein")])
        path = self.fixture.genome_dir("G1") / "protein.faa.gz"
        with gzip.open(path, "wt") as source:
            source.write(">P1 kinase [Species]\nMAK\n")
        self.assertEqual(self.run_script("P1"),
                         (0, "accession_id\tprotein_name\nP1\tkinase\n", ""))

    def test_missing_manifest(self):
        status, output, errors = self.run_script("P1\n")
        self.assertEqual((status, output), (1, ""))
        self.assertIn("sequences.tsv", errors)

    def test_missing_accession_file(self):
        with redirect_stdout(io.StringIO()) as output, redirect_stderr(io.StringIO()) as errors:
            self.assertEqual(self.script.main([str(self.input)]), 1)
        self.assertEqual(output.getvalue(), "")
        self.assertIn("accessions.txt", errors.getvalue())

    def test_lookup_failures_never_emit_partial_tsv(self):
        self.manifest([("P1", "G1", "protein"), ("P2", "G2", "protein"),
                       ("DNA", "G1", "contig")])
        self.fasta("G1", ["P1 kinase"])
        for content, expected in [("P1\nUNKNOWN", "manifest"), ("DNA", "manifest"),
                                  ("P1\nP2", "protein.faa")]:
            with self.subTest(content=content):
                status, output, errors = self.run_script(content)
                self.assertEqual((status, output), (1, ""))
                self.assertIn(expected, errors)
        for headers, expected in [(["P2.1 kinase"], "Cannot find proteins"),
                                  (["P2"], "Empty protein name"),
                                  (["P2 [Species]"], "Empty protein name"),
                                  (["P2 kinase", "P2 phosphatase"], "Conflicting")]:
            with self.subTest(headers=headers):
                self.fasta("G2", headers)
                status, output, errors = self.run_script("P1\nP2")
                self.assertEqual((status, output), (1, ""))
                self.assertIn(expected, errors)

    def test_conflicting_names_across_genomes(self):
        self.manifest([("P1", "G1", "protein"), ("P1", "G2", "protein")])
        self.fasta("G1", ["P1 kinase"])
        self.fasta("G2", ["P1 phosphatase"])
        status, output, errors = self.run_script("P1")
        self.assertEqual((status, output), (1, ""))
        self.assertIn("across genomes", errors)

    def test_missing_genome_and_quoted_accession(self):
        self.manifest([("P'1", "G1", "protein"), ("P2", "", "protein")])
        self.fasta("G1", ["P'1 kinase", "P'1 kinase"])
        self.assertEqual(self.run_script("P'1"),
                         (0, "accession_id\tprotein_name\nP'1\tkinase\n", ""))
        status, output, errors = self.run_script("P2")
        self.assertEqual((status, output), (1, ""))
        self.assertIn("Missing genome", errors)
