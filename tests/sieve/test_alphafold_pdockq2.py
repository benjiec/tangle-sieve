import math
import unittest

from sieve.alphafold_pdockq2 import (
    Residue,
    collect_regions,
    parse_region,
    score_pair,
    score_model,
    score_pairs,
    validate_regions,
)


class RegionTests(unittest.TestCase):
    def test_parse_region(self):
        self.assertEqual(parse_region("A:20-40"), ("A", 20, 40))

    def test_rejects_invalid_region(self):
        for value in ("A", "A:0-2", "A:4-2", "A:x-2"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                parse_region(value)

    def test_collects_disjoint_regions_as_union(self):
        self.assertEqual(
            collect_regions(["A:20-40", "A:1-10", "B:2-3"]),
            {"A": [(1, 10), (20, 40)], "B": [(2, 3)]},
        )

    def test_rejects_overlapping_regions(self):
        with self.assertRaisesRegex(ValueError, "overlaps"):
            collect_regions(["A:1-10", "A:10-20"])

    def test_rejects_unknown_chain_and_absent_residue(self):
        residues = {"A": [Residue("A", 1, 0, 90.0, (0.0, 0.0, 0.0))]}
        with self.assertRaisesRegex(ValueError, "unknown chain"):
            validate_regions({"B": [(1, 1)]}, residues)
        with self.assertRaisesRegex(ValueError, "absent residue 2"):
            validate_regions({"A": [(1, 2)]}, residues)


class ScoreTests(unittest.TestCase):
    def setUp(self):
        self.residues = {
            "A": [
                Residue("A", 1, 0, 90.0, (0.0, 0.0, 0.0)),
                Residue("A", 2, 1, 50.0, (100.0, 0.0, 0.0)),
            ],
            "B": [
                Residue("B", 1, 2, 80.0, (0.0, 0.0, 7.0)),
                Residue("B", 2, 3, 40.0, (200.0, 0.0, 0.0)),
            ],
        }
        self.pae = [
            [0.0, 0.0, 2.0, 30.0],
            [0.0, 0.0, 30.0, 30.0],
            [20.0, 30.0, 0.0, 0.0],
            [30.0, 30.0, 0.0, 0.0],
        ]

    def test_scores_directional_pae(self):
        row = score_pair("A", "B", self.residues, self.pae, {}, cutoff=8.0)
        self.assertEqual(row["contact count"], 1)
        self.assertEqual(row["interface residues 1"], "1")
        self.assertEqual(row["interface residues 2"], "1")
        self.assertGreater(row["pDockQ2 1 to 2"], row["pDockQ2 2 to 1"])
        self.assertEqual(row["pDockQ2 max"], row["pDockQ2 1 to 2"])
        self.assertTrue(math.isclose(row["normalized PAE 1 to 2"], 1 / 1.04))

    def test_region_on_one_chain_leaves_other_unrestricted(self):
        row = score_pair(
            "A", "B", self.residues, self.pae, {"A": [(1, 1)]}, cutoff=8.0
        )
        self.assertEqual(row["region 1"], "1-1")
        self.assertEqual(row["region 2"], "all")
        self.assertEqual(row["contact count"], 1)

    def test_no_contacts_returns_zero(self):
        row = score_pair("A", "B", self.residues, self.pae, {}, cutoff=1.0)
        self.assertEqual(row["contact count"], 0)
        self.assertEqual(row["pDockQ2 max"], 0.0)

    def test_regions_on_both_chains_can_exclude_contacts(self):
        row = score_pair(
            "A",
            "B",
            self.residues,
            self.pae,
            {"A": [(2, 2)], "B": [(2, 2)]},
            cutoff=8.0,
        )
        self.assertEqual(row["contact count"], 0)

    def test_rejects_nonpositive_cutoff_before_parsing(self):
        with self.assertRaisesRegex(ValueError, "greater than zero"):
            score_model("", {}, cutoff=0)

    def test_full_rows_are_always_included_before_regional_rows(self):
        residues = dict(self.residues)
        residues["C"] = [Residue("C", 1, 4, 70.0, (300.0, 0.0, 0.0))]
        pae = [row + [30.0] for row in self.pae] + [[30.0] * 4 + [0.0]]
        rows = score_pairs(
            residues,
            pae,
            {"A": [(1, 1)], "B": [(1, 1)]},
            cutoff=8.0,
        )
        self.assertEqual(len(rows), 6)
        self.assertEqual([row["scope"] for row in rows], [
            "full", "regional", "full", "regional", "full", "regional"
        ])
        self.assertEqual(rows[0]["region 1"], "all")
        self.assertEqual(rows[1]["region 1"], "1-1")

    def test_pair_unaffected_by_regions_has_no_duplicate(self):
        residues = dict(self.residues)
        residues["C"] = [Residue("C", 1, 4, 70.0, (300.0, 0.0, 0.0))]
        pae = [row + [30.0] for row in self.pae] + [[30.0] * 4 + [0.0]]
        rows = score_pairs(residues, pae, {"A": [(1, 1)]}, cutoff=8.0)
        self.assertEqual(len(rows), 5)
        bc_rows = [row for row in rows if row["chain 1"] == "B" and row["chain 2"] == "C"]
        self.assertEqual(len(bc_rows), 1)
        self.assertEqual(bc_rows[0]["scope"], "full")


if __name__ == "__main__":
    unittest.main()
