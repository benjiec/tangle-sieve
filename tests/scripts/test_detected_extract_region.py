import contextlib
import io
import os
import tempfile
import unittest

from tangle.detected import DetectedTable

from tests.scripts.helpers import load_script


class TestDetectedExtractRegionScript(unittest.TestCase):

    def setUp(self):
        repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.script = load_script(os.path.join(repo, "scripts", "detected-extract-region.py"))

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

    def test_indices_and_offsets(self):
        cases = [
            ({}, "BCDEFGHI"),
            ({"left_index": 2}, "FGHI"),
            ({"right_index": 2}, "BCDEFG"),
            ({"left_index": 2, "right_index": 2}, "FG"),
            ({"left_offset": -1, "right_offset": 1}, "ABCDEFGHIJ"),
            ({"left_offset": 1, "right_offset": -1}, "CDEFGH"),
            ({"left_index": 2, "right_index": 2,
              "left_offset": -1, "right_offset": 1}, "EFGH"),
            ({"left_index": 2, "right_index": 2, "right_offset": -1}, "F"),
        ]
        for options, expected in cases:
            with self.subTest(options=options):
                reports = []
                self.assertEqual(self.script.extract_regions(
                    {"p1": "ABCDEFGHIJ"}, {"p1": [(8, 9), (2, 3), (6, 7)]},
                    "PF00023", 3, 3, report=reports.append, **options,
                ), {"p1_PF00023_3_3": expected})
                self.assertEqual(reports, [])

    def test_invalid_indices_and_adjusted_boundaries_are_skipped(self):
        cases = [
            ({side: value}, "domain indices")
            for side in ("left_index", "right_index") for value in (-1, 0, 4)
        ] + [
            ({"left_offset": -2}, "outside sequence length"),
            ({"left_offset": 9}, "outside sequence length"),
            ({"right_offset": 2}, "outside sequence length"),
            ({"right_offset": -9}, "outside sequence length"),
            ({"left_index": 3, "right_index": 1}, "exceeds end"),
            ({"left_offset": 4, "right_offset": -4}, "exceeds end"),
        ]
        for options, diagnostic in cases:
            with self.subTest(options=options):
                reports = []
                self.assertEqual(self.script.extract_regions(
                    {"p1": "ABCDEFGHIJ"}, {"p1": [(2, 3), (6, 7), (8, 9)]},
                    "PF00023", 3, 3, report=reports.append, **options,
                ), {})
                self.assertEqual(len(reports), 1)
                self.assertIn(diagnostic, reports[0])

    def test_nested_and_tied_domains_preserve_default_span(self):
        for options, expected in [({}, "BCDEFGHI"),
                                  ({"right_index": 3}, "BCDEF"),
                                  ({"right_index": 1}, "BCD")]:
            with self.subTest(options=options):
                self.assertEqual(self.script.extract_regions(
                    {"p1": "ABCDEFGHIJ"}, {"p1": [(5, 6), (2, 9), (2, 4)]},
                    "PF00023", 3, 3, **options,
                ), {"p1_PF00023_3_3": expected})

    def test_main_passes_indices_and_signed_offsets(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            tsv = os.path.join(tmpd, "matches.tsv")
            output = os.path.join(tmpd, "output.faa")
            with open(fasta, "w", encoding="utf-8") as stream:
                stream.write(">p1\nABCDEFGHIJ\n")
            DetectedTable.write_tsv(tsv, [
                self.detected_row("p1", "PF00023", 8, 9),
                self.detected_row("p1", "PF00023", 2, 3),
                self.detected_row("p1", "PF00023", 7, 6),
            ])
            self.assertEqual(self.script.main([
                fasta, tsv, "PF00023", "3", "3", output,
                "--left-index", "2", "--right-index", "2",
                "--left-offset", "-1", "--right-offset", "+1",
            ]), 0)
            with open(output, encoding="utf-8") as stream:
                self.assertEqual(stream.read(), ">p1_PF00023_3_3\nEFGH\n")

    def test_cli_requires_integer_options(self):
        for option in ("--left-index", "--right-index", "--left-offset", "--right-offset"):
            with self.subTest(option=option), self.assertRaises(SystemExit):
                self.script.main([
                    "in.faa", "in.tsv", "PF00023", "1", "3", "out.faa",
                    option, "1.5",
                ])

    def test_rejects_invalid_limits(self):
        with self.assertRaises(ValueError):
            self.script.extract_regions({}, {}, "PF00023", 0, 9)
        with self.assertRaises(ValueError):
            self.script.extract_regions({}, {}, "PF00023", 10, 9)

        with self.assertRaises(SystemExit):
            self.script.main(["in.faa", "in.tsv", "PF00023", "9", "4", "out.faa"])

    def test_domain_only_filters_counts_then_extracts_each_domain(self):
        reports = []
        result = self.script.extract_regions(
            {name: "ABCDEFGHIJ" for name in ("one", "two", "three", "none")},
            {"one": [(1, 1)], "two": [(8, 10), (2, 4)],
             "three": [(1, 2), (4, 5), (8, 9)]},
            "PF13676", 1, 2, domain_only=True, report=reports.append,
        )
        self.assertEqual(result, {
            "one_PF13676_1_2_domain_1": "A",
            "two_PF13676_1_2_domain_1": "BCD",
            "two_PF13676_1_2_domain_2": "HIJ",
        })
        self.assertEqual(len(reports), 2)
        self.assertIn("found 3", reports[0])
        self.assertIn("found 0", reports[1])

    def test_domain_only_offsets_and_individual_invalid_regions(self):
        cases = [
            (-1, 1, {1: "ABCDE", 2: "EFGHIJ"}, []),
            (1, -1, {1: "C", 2: "GH"}, []),
            (-2, 0, {2: "DEFGHI"}, ["domain 1", "outside sequence length"]),
            (0, 2, {1: "BCDEF"}, ["domain 2", "outside sequence length"]),
            (2, -1, {2: "H"}, ["domain 1", "exceeds end"]),
            (10, 0, {}, ["outside sequence length"]),
            (0, -10, {}, ["outside sequence length"]),
        ]
        for left, right, expected, diagnostics in cases:
            with self.subTest(left=left, right=right):
                reports = []
                result = self.script.extract_regions(
                    {"p": "ABCDEFGHIJ"}, {"p": [(2, 4), (6, 9)]},
                    "PF13676", 2, 2, domain_only=True,
                    left_offset=left, right_offset=right, report=reports.append,
                )
                self.assertEqual(result, {
                    f"p_PF13676_2_2_domain_{i}": seq for i, seq in expected.items()
                })
                self.assertEqual(len(reports), 2 - len(expected))
                for diagnostic in diagnostics:
                    self.assertIn(diagnostic, "\n".join(reports))

    def test_domain_only_keeps_overlapping_domains_separate(self):
        self.assertEqual(self.script.extract_regions(
            {"p": "ABCDEFGHIJ"}, {"p": [(5, 6), (2, 9), (2, 4)]},
            "PF13676", 3, 3, domain_only=True,
        ), {
            "p_PF13676_3_3_domain_1": "BCD",
            "p_PF13676_3_3_domain_2": "BCDEFGHI",
            "p_PF13676_3_3_domain_3": "EF",
        })

    def test_domain_only_rejects_indices_before_reading_inputs(self):
        for side in ("left", "right"):
            with self.subTest(side=side):
                with self.assertRaisesRegex(ValueError, "cannot be combined"):
                    self.script.extract_regions(
                        {}, {}, "PF13676", 1, 2, domain_only=True,
                        **{f"{side}_index": 1},
                    )
                stderr = io.StringIO()
                with contextlib.redirect_stderr(stderr), self.assertRaises(SystemExit) as error:
                    self.script.main([
                        "missing.faa", "missing.tsv", "PF13676", "1", "2", "out.faa",
                        "--domain-only", f"--{side}-index", "1",
                    ])
                self.assertEqual(error.exception.code, 2)
                self.assertIn("cannot be combined", stderr.getvalue())

    def test_main_domain_only_with_offsets(self):
        with tempfile.TemporaryDirectory() as tmpd:
            fasta = os.path.join(tmpd, "input.faa")
            tsv = os.path.join(tmpd, "matches.tsv")
            output = os.path.join(tmpd, "output.faa")
            with open(fasta, "w", encoding="utf-8") as stream:
                stream.write(">p\nABCDEFGHIJ\n")
            row = self.detected_row("p", "PF13676", 2, 4)
            DetectedTable.write_tsv(tsv, [
                self.detected_row("p", "PF13676.1", 9, 6), row, row,
                self.detected_row("p", "PF00023", 1, 10),
            ])
            self.assertEqual(self.script.main([
                fasta, tsv, "PF13676", "1", "2", output,
                "--domain-only", "--left-offset", "-1", "--right-offset", "+1",
            ]), 0)
            with open(output, encoding="utf-8") as stream:
                self.assertEqual(stream.read(),
                    ">p_PF13676_1_2_domain_1\nABCDE\n"
                    ">p_PF13676_1_2_domain_2\nEFGHIJ\n")


if __name__ == "__main__":
    unittest.main()
