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

    def test_match_only_writes_every_matched_region_with_coordinates(self):
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [
            self.detected_row("p1", "g1", "PF00504.27", 1e-20) | {
                "query_start": 2,
                "query_end": 4,
            },
            self.detected_row("p1", "g1", "PF00504.27", 1e-20) | {
                "query_start": 6,
                "query_end": 8,
            },
        ])
        self.fx.write_manifest([{
            "sequence_accession": "p1",
            "sequence_database": "g1",
            "sequence_type": "protein",
            "sequence_source": SEQUENCE_SOURCE_NCBI,
        }])
        self.fx.write_ncbi_proteins("g1", {"p1": "ABCDEFGHIJ"})

        with tempfile.TemporaryDirectory() as tmpd:
            output = os.path.join(tmpd, "matches.faa")
            self.script.main(["PF00504", "--match-only", "-o", output])

            with open(output, encoding="utf-8") as f:
                self.assertEqual(
                    f.read(),
                    ">p1_PF00504_2_4\nBCD\n>p1_PF00504_6_8\nFGH\n",
                )

    def test_match_only_requires_output(self):
        with self.assertRaises(SystemExit):
            self.script.main(["PF00504", "--match-only"])

    def test_multiple_targets_require_all_passing_hits_in_same_database(self):
        hits = [
            ("complete", "g1", "PF07714.3", 1e-20),
            ("complete", "g1", "PF00069", 1e-20),
            ("complete", "g1", "PF00531.1", 1e-20),
            ("complete", "g1", "PF99999", 1e-20),
            ("complete", "g1", "PF07714.3", 1e-20),
            ("partial", "g1", "PF07714", 1e-20),
            ("partial", "g1", "PF00069", 1e-20),
            ("partial", "g1", "PF00531", 1e-2),
            ("split", "g1", "PF07714", 1e-20),
            ("split", "g2", "PF00069", 1e-20),
            ("split", "g2", "PF00531", 1e-20),
        ]
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [
            self.detected_row(*hit) | {"query_start": {"PF07714.3": 1, "PF07714": 1, "PF00069": 20, "PF00531.1": 40, "PF00531": 40, "PF99999": 60}[hit[2]], "query_end": {"PF07714.3": 10, "PF07714": 10, "PF00069": 30, "PF00531.1": 50, "PF00531": 50, "PF99999": 70}[hit[2]]} for hit in hits
        ])
        targets = ["PF07714.9", "PF00069", "PF00531", "PF07714"]
        self.assertEqual(self.script.find_matches(targets, max_evalue=1e-10), [("complete", "g1")])
        stdout = io.StringIO()
        with patch("sys.stdout", stdout):
            self.script.main([*targets, "--max-evalue", "1e-10"])
        self.assertEqual(stdout.getvalue(), "complete\tg1\n")
        self.assertEqual(self.script.find_matches([*targets, "PF12345"]), [])

    def test_exact_match_and_version_normalization(self):
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [
            self.detected_row("p1", "g1", "PF00504", 1e-20),
            self.detected_row("p2", "g1", "PF00504.27", 1e-20),
            self.detected_row("p3", "g1", "PF00504.28", 1e-20),
            self.detected_row("p3", "g1", "PF00504.27", 1e-20),
            self.detected_row("p4", "g1", "PF005040.27", 1e-20),
        ])
        self.assertEqual(self.script.find_matches("PF00504.99"), [("p1", "g1"), ("p2", "g1"), ("p3", "g1")])
        self.assertEqual(self.script.find_matches("PF00504", exact_match=True), [("p1", "g1")])
        self.assertEqual(self.script.find_matches(["PF00504.27", "PF00504.28"], exact_match=True, allow_overlap=True), [("p3", "g1")])
        stdout = io.StringIO()
        with patch("sys.stdout", stdout):
            self.script.main(["PF00504.27", "--exact-match"])
        self.assertEqual(stdout.getvalue(), "p2\tg1\np3\tg1\n")

    def test_empty_table(self):
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [])
        self.assertEqual(self.script.find_matches(["PF07714", "PF00069"]), [])

    def test_sequence_source_filters_manifest_by_accession_database_and_type(self):
        proteins = [("p1", "g1"), ("p1", "g2"), ("missing", "g1"),
                    ("blank", "g1"), ("gene", "g1")]
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [
            self.detected_row(p, g, "PF00504", 1e-20) for p, g in proteins
        ])
        self.fx.write_manifest([
            dict(sequence_accession=p, sequence_database=g, sequence_type=t, sequence_source=s)
            for p, g, t, s in [
                ("p1", "g1", "protein", "custom's source"),
                ("p1", "g2", "protein", "other"),
                ("blank", "g1", "protein", None),
                ("gene", "g1", "gene", "custom's source"),
            ]
        ])
        self.assertEqual(self.script.find_matches("PF00504", sequence_source="custom's source"), [("p1", "g1")])
        self.assertEqual(self.script.find_matches("PF00504", sequence_source="CUSTOM'S SOURCE"), [])
        with patch.object(self.script.Defaults, "area_sequence_manifest_tsv", side_effect=AssertionError("Unexpected manifest read")):
            self.assertEqual(len(self.script.find_matches("PF00504")), len(proteins))
        self.fx.write_manifest([])
        self.assertEqual(self.script.find_matches("PF00504", sequence_source="other"), [])

    def test_sequence_source_combines_with_filters_and_output_modes(self):
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [
            self.detected_row(p, g, target, evalue) | ({"query_start": 20, "query_end": 30} if target == "PF00069" else {})
            for p, g, evalue in [("keep", "g1", 1e-20), ("weak", "g1", 1),
                                 ("taxon", "g2", 1e-20), ("source", "g1", 1e-20)]
            for target in ["PF00504.27", "PF00069"]
        ])
        self.fx.write_manifest([
            dict(sequence_accession=p, sequence_database=g, sequence_type="protein", sequence_source=s)
            for p, g, s in [("keep", "g1", SEQUENCE_SOURCE_NCBI), ("weak", "g1", SEQUENCE_SOURCE_NCBI),
                            ("taxon", "g2", SEQUENCE_SOURCE_NCBI), ("source", "g1", "other")]
        ])
        self.fx.write_taxonomy_rows([
            {"Genome Accession": "g1", "Phylum": "Cnidaria"},
            {"Genome Accession": "g2", "Phylum": "Alveolata"},
        ])
        self.fx.write_ncbi_proteins("g1", {"keep": "ABCDEFGHIJKLM"})
        options = ["--sequence-source", SEQUENCE_SOURCE_NCBI, "--taxon", "cnidaria",
                   "--max-evalue", "1e-10", "--exact-match"]
        stdout = io.StringIO()
        with patch("sys.stdout", stdout):
            self.script.main(["PF00504.27", "PF00069", *options])
        self.assertEqual(stdout.getvalue(), "keep\tg1\n")
        with tempfile.TemporaryDirectory() as tmpd:
            output = os.path.join(tmpd, "matches.faa")
            self.script.main(["PF00504.27", "PF00069", *options, "-o", output])
            with open(output) as f:
                self.assertEqual(f.read(), ">keep\nABCDEFGHIJKLM\n")
            self.script.main(["PF00504.27", *options, "--match-only", "-o", output])
            with open(output) as f:
                self.assertEqual(f.read(), ">keep_PF00504.27_1_10\nABCDEFGHIJ\n")

    def test_invalid_arguments(self):
        for args in ([], ["PF00504", "PF00069", "--match-only", "-o", "unused.faa"],
                     ["PF00504", "PF00504", "--match-only"]):
            with self.subTest(args=args), patch("sys.stderr", io.StringIO()) as stderr:
                with self.assertRaises(SystemExit) as error:
                    self.script.main(args)
                self.assertEqual(error.exception.code, 2)
                if args:
                    self.assertIn("--match-only requires a single Pfam accession", stderr.getvalue())
        with self.assertRaises(ValueError):
            self.script.find_matches([])
        with self.assertRaises(ValueError):
            self.script.find_matches(["PF00504", "PF00069"], include_coordinates=True)

    def test_architecture_quantifiers_and_order(self):
        cases = [
            (["A", "B"], "AAB", True, True),
            (["A", "B"], "BAA", True, False),
            (["A", "B"], "BAA", False, True),
            (["A?", "B"], "AAB", True, False),
            (["A?", "B"], "AAB", False, False),
            (["A?", "B"], "B", True, True),
            (["A+", "B"], "B", False, False),
            (["A*", "B"], "B", True, True),
            (["A{3,5}"], "AA", False, False),
            (["A{3,5}"], "AAA", True, True),
            (["A{3,5}"], "AAAAA", False, True),
            (["A{3,5}"], "AAAAAA", True, False),
            (["A{2}"], "AA", False, True),
            (["A{0}"], "A", False, False),
            (["A?"], "", True, True),
            (["A", "B", "A"], "ABA", True, True),
            (["A", "B", "A"], "AAB", True, False),
            (["A?", "A", "B"], "AB", True, True),
            (["A*", "A{2}", "B"], "AAAB", True, True),
        ]
        for expressions, architecture, ordered, expected in cases:
            with self.subTest(expressions=expressions, architecture=architecture, ordered=ordered):
                hits = [(a, i * 10 + 1, i * 10 + 5) for i, a in enumerate(architecture)]
                self.assertEqual(self.script.architecture_matches(
                    hits, self.script.parse_targets(expressions), ordered=ordered), expected)

    def test_overlap_ties_and_duplicate_hits(self):
        terms = self.script.parse_targets(["A", "B"])
        for hits in ([("A", 1, 10), ("B", 10, 20)],
                     [("A", 1, 30), ("B", 10, 20)]):
            for ordered in (False, True):
                self.assertFalse(self.script.architecture_matches(hits, terms, ordered))
                self.assertTrue(self.script.architecture_matches(hits, terms, ordered, True))
        tied = [("A", 1, 10), ("B", 1, 20)]
        self.assertFalse(self.script.architecture_matches(tied, terms, True, True))
        self.assertTrue(self.script.architecture_matches(tied, terms, False, True))
        self.assertTrue(self.script.architecture_matches(
            [("A", 1, 10), ("A", 1, 10)], self.script.parse_targets(["A{1}"])))

    def test_ordered_cli_ignores_other_domains_and_checks_all_requested_hits(self):
        architectures = {"forward": "AXB", "reverse": "BXA", "excess": "AAXB", "absent": "X"}
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [
            self.detected_row(p, "g1", a + ".1", 1e-20) |
            {"query_start": i * 20 + 1, "query_end": i * 20 + 10}
            for p, architecture in architectures.items() for i, a in enumerate(architecture)
        ])
        stdout = io.StringIO()
        with patch("sys.stdout", stdout):
            self.script.main(["A?", "B", "--ordered"])
        self.assertEqual(stdout.getvalue(), "forward\tg1\n")
        self.assertEqual(self.script.find_matches(["A?", "B"]), [("forward", "g1"), ("reverse", "g1")])
        self.assertEqual(self.script.find_matches(["A{0}"], ordered=True), [("absent", "g1")])
        self.assertEqual(self.script.find_matches(["A.1?", "B.1"], exact_match=True, ordered=True), [("forward", "g1")])

    def test_allow_overlap_cli_and_reversed_coordinates(self):
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), [
            self.detected_row("p", "g", "A", 1e-20) | {"query_start": 10, "query_end": 1},
            self.detected_row("p", "g", "B", 1e-20) | {"query_start": 5, "query_end": 15},
        ])
        self.assertEqual(self.script.find_matches(["A", "B"]), [])
        stdout = io.StringIO()
        with patch("sys.stdout", stdout):
            self.script.main(["A", "B", "--ordered", "--allow-overlap"])
        self.assertEqual(stdout.getvalue(), "p\tg\n")

    def test_invalid_quantifiers_report_argument_errors(self):
        for expression in ("A{5,3}", "A{,3}", "A{3,}", "A{-1}", "A??", "A+?", "A{a}", "?", ""):
            with self.subTest(expression=expression), patch("sys.stderr", io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    self.script.main([expression])
                self.assertEqual(error.exception.code, 2)
