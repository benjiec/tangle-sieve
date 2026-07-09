import os
import sys
import tempfile
import unittest

from sieve.result_filters import Field, FieldRegex, LeaderCall, load_result_filter


class TestResultFilters(unittest.TestCase):

    def test_field_exact_literal_and_numeric_comparisons_are_distinct(self):
        row = {"Literal": "001", "Numeric": "001"}

        self.assertFalse(Field("Literal").eq("1").matches(row))
        self.assertTrue(Field("Numeric").num_eq(1).matches(row))
        self.assertTrue(Field("Numeric").gt(0).matches(row))

    def test_field_regex_any_matches_any_matching_column(self):
        row = {
            "Pfam.matches('PF00001')": "false",
            "Pfam.matches('PF00002')": "true",
        }

        self.assertTrue(FieldRegex(r"Pfam\.matches.+").any().eq("true").matches(row))

    def test_field_regex_all_requires_all_matching_columns(self):
        row = {
            "HMMAlignment('a').is_at('H', 1)": "true",
            "HMMAlignment('a').is_at('H', 2)": "false",
        }

        self.assertFalse(FieldRegex(r"HMMAlignment.+").all().eq("true").matches(row))
        row["HMMAlignment('a').is_at('H', 2)"] = "true"
        self.assertTrue(FieldRegex(r"HMMAlignment.+").all().eq("true").matches(row))

    def test_field_regex_requires_explicit_all_or_any(self):
        with self.assertRaisesRegex(ValueError, "requires .all\\(\\) or .any\\(\\)"):
            FieldRegex(r"Pfam\.matches.+").eq("true").matches({"Pfam.matches('PF00001')": "true"})

    def test_composes_filters_with_python_grouping(self):
        a = Field("A").eq("true")
        b = Field("B").eq("true")
        c = Field("C").eq("true")
        grouped = (a | b) & c

        self.assertTrue(grouped.matches({"A": "false", "B": "true", "C": "true"}))
        self.assertFalse(grouped.matches({"A": "true", "B": "false", "C": "false"}))

    def test_inverts_filters_and_matches_value_regexes(self):
        row = {"description": "mitochondrial MnSOD candidate"}

        self.assertTrue(Field("description").matches(r"MnSOD").matches(row))
        self.assertTrue((~Field("description").matches(r"chloroplast")).matches(row))

    def test_leader_call_helper_targets_split_leader_column(self):
        row = {"Leader.call('mTP')": "80"}

        self.assertTrue(LeaderCall("mTP").gte(80).matches(row))

    def test_ge_and_le_alias_numeric_comparisons(self):
        row = {"score": "10"}

        self.assertTrue(Field("score").ge(10).matches(row))
        self.assertTrue(Field("score").le(10).matches(row))

    def test_validates_missing_columns_before_filtering_rows(self):
        filter_ = Field("Missing").eq("true") | Field("Rule").eq("true")

        with self.assertRaisesRegex(ValueError, "Missing"):
            filter_.validate_columns(["Rule"])

    def test_load_result_filter_imports_module_attribute(self):
        with tempfile.TemporaryDirectory() as tmpd:
            module_path = os.path.join(tmpd, "filters.py")
            with open(module_path, "w", encoding="utf-8") as f:
                f.write("\n".join([
                    "from sieve.result_filters import Field",
                    "is_positive = Field('Rule').eq('true')",
                    "",
                ]))

            sys.path.insert(0, tmpd)
            try:
                loaded = load_result_filter("filters.is_positive")
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("filters", None)

        self.assertTrue(loaded.matches({"Rule": "true"}))


if __name__ == "__main__":
    unittest.main()
