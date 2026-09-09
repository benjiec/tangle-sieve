import os
import tempfile
import unittest

from tangle.detected import DetectedTable

from tests.scripts.helpers import load_script


class TestPfamExtractRegionScript(unittest.TestCase):

    def setUp(self):
        repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.script = load_script(os.path.join(repo, "scripts", "pfam-extract-region.py"))

    def detected_row(self, accession, pfam, start, end):
        return dict(
            detection_type="sequence",
            detection_method="hmm",
            batch="b1",
            query_accession=accession,
            query_database="input",
            query_type="protein",
            target_accession=pfam,
            target_database="Pfam",
            target_type="protein",
            query_start=start,
            query_end=end,
        )

    def test_finds_exact_and_versioned_matches_and_removes_duplicate_rows(self):
        with tempfile.TemporaryDirectory() as tmpd:
            tsv = os.path.join(tmpd, "matches.tsv")
            row = self.detected_row("p1", "PF00023", 2, 4)
            DetectedTable.write_tsv(tsv, [
                row,
                row,
                self.detected_row("p1", "PF00023.35", 8, 6),
                self.detected_row("p1", "PF000230", 9, 10),
            ])
            self.assertEqual(
                self.script.find_match_regions(tsv, "PF00023"),
                {"p1": [(2, 4), (6, 8)]},
            )

    def test_extracts_first_start_through_last_end_at_inclusive_limits(self):
        reports = []
        extracted = self.script.extract_regions(
            {"four": "ABCDEFGHIJKLM", "nine": "ABCDEFGHIJKLM", "three": "ABCDE"},
            {
                "four": [(8, 9), (2, 3), (4, 5), (6, 7)],
                "nine": [(i, i) for i in range(1, 10)],
                "three": [(1, 1), (2, 2), (3, 3)],
            },
            "PF00023",
            4,
            9,
            report=reports.append,
        )
        self.assertEqual(extracted, {
            "four_PF00023_4_9": "BCDEFGHI",
            "nine_PF00023_4_9": "ABCDEFGHI",
        })
        self.assertEqual(
            reports,
            ["Ignoring three: expected 4-9 PF00023 matches, found 3"],
        )

    def test_rejects_too_many_missing_and_out_of_bounds_matches(self):
        reports = []
        extracted = self.script.extract_regions(
            {"ten": "ABCDEFGHIJ", "none": "ABC", "outside": "ABCDE"},
            {
                "ten": [(i, i) for i in range(1, 11)],
                "outside": [(1, 1), (2, 2), (3, 3), (4, 6)],
            },
            "PF00023",
            4,
            9,
            report=reports.append,
        )
        self.assertEqual(extracted, {})
        self.assertEqual(len(reports), 3)
        self.assertIn("found 10", reports[0])
        self.assertIn("found 0", reports[1])
        self.assertIn("outside sequence length 5", reports[2])

    def test_preserves_x_in_extracted_regions(self):
        reports = []
        matches = [(2, 2), (3, 3), (4, 4), (5, 5)]
        extracted = self.script.extract_regions(
            {
                "inside": "ABXDEF",
                "outside": "ABCDEX",
                "boundaries": "AXCDXF",
            },
            {
                "inside": matches,
                "outside": matches,
                "boundaries": matches,
            },
            "PF00023",
            4,
            9,
            report=reports.append,
        )
        self.assertEqual(extracted, {
            "inside_PF00023_4_9": "BXDE",
            "outside_PF00023_4_9": "BCDE",
            "boundaries_PF00023_4_9": "XCDX",
        })
        self.assertEqual(reports, [])

    def test_main_reads_and_writes_fasta(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            tsv = os.path.join(tmpd, "matches.tsv")
            output = os.path.join(tmpd, "output.faa")
            with open(fasta, "w", encoding="utf-8") as stream:
                stream.write(">p1 description\nABCDEFGHIJ\n")
            DetectedTable.write_tsv(tsv, [
                self.detected_row("p1", "PF00023", 2, 3),
                self.detected_row("p1", "PF00023", 4, 5),
                self.detected_row("p1", "PF00023", 6, 7),
                self.detected_row("p1", "PF00023", 8, 9),
            ])

            self.assertEqual(
                self.script.main([fasta, tsv, "PF00023", "4", "9", output]),
                0,
            )
            with open(output, encoding="utf-8") as stream:
                self.assertEqual(stream.read(), ">p1_PF00023_4_9\nBCDEFGHI\n")

    def test_rejects_duplicate_fasta_accessions(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            with open(fasta, "w", encoding="utf-8") as stream:
                stream.write(">p1 first\nABC\n>p1 second\nDEF\n")
            with self.assertRaisesRegex(ValueError, "Duplicate FASTA accession: p1"):
                self.script.read_unique_fasta(fasta)

    def test_rejects_invalid_limits(self):
        with self.assertRaises(ValueError):
            self.script.extract_regions({}, {}, "PF00023", 0, 9)
        with self.assertRaises(ValueError):
            self.script.extract_regions({}, {}, "PF00023", 10, 9)

        with self.assertRaises(SystemExit):
            self.script.main(["in.faa", "in.tsv", "PF00023", "9", "4", "out.faa"])


if __name__ == "__main__":
    unittest.main()
