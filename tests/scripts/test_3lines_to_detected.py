import csv
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from tests.scripts.helpers import load_script


class TestThreeLinesToDetected(unittest.TestCase):
    def setUp(self):
        self.script = load_script(str(
            Path(__file__).resolve().parents[2] / "scripts" / "deeptm-to-detected.py"
        ))

    def rows(self, text):
        return list(self.script.detected_rows(
            self.script.parse_records(io.StringIO(text)), "test-batch"
        ))

    def test_coordinates_and_metadata(self):
        rows = self.rows(">NP_001056.1 | SP+TM\nABCDEFGHI\nSSSOOMMII\n")
        self.assertEqual([
            (r["target_accession"], r["query_start"], r["query_end"],
             r["target_start"], r["target_end"])
            for r in rows
        ], [
            ("signal", 1, 3, 1, 3),
            ("outside", 4, 5, 1, 2),
            ("membrane", 6, 7, 1, 2),
            ("inside", 8, 9, 1, 2),
        ])
        for row in rows:
            self.assertNotIn("target_model", row)
            self.assertEqual(row["query_accession"], "NP_001056.1")
            self.assertEqual(row["query_database"], "_")
            self.assertEqual(row["detection_type"], "sequence")
            self.assertEqual(row["detection_method"], "other")
            self.assertEqual(row["target_database"], "DeepTMHMM")
            self.assertEqual(row["batch"], "test-batch")

    def test_repeated_and_single_residue_regions(self):
        rows = self.rows(">p\nABCDE\nIMOMI\n")
        self.assertEqual([r["target_accession"] for r in rows],
                         ["inside", "membrane", "outside", "membrane", "inside"])
        self.assertEqual([(r["query_start"], r["query_end"]) for r in rows],
                         [(i, i) for i in range(1, 6)])

    def test_multiple_records_blank_separators_crlf_and_no_final_newline(self):
        rows = self.rows("\r\n>p|SP\r\nABC\r\nSSS\r\n\r\n>q\r\nD\r\nO")
        self.assertEqual([(r["query_accession"], r["query_start"], r["query_end"])
                          for r in rows], [("p", 1, 3), ("q", 1, 1)])

    def test_rejects_malformed_records(self):
        for text in (
            "ABC\nSSS\n", ">\nA\nS\n", "> | SP\nA\nS\n",
            ">p\n", ">p\nA\n", ">p\n\nS\n", ">p\nA\n\n",
            ">p\n>q\nS\n", ">p\nA\n>q\n", ">p\nAA\nS\n",
            ">p\nA\nX\n", ">p\nA\ns\n", ">p\nA A\nSSS\n",
            ">p\nAAA\nS S\n", ">p\nA\nS\nextra\n",
        ):
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.rows(text)

    def test_cli_writes_table_and_calls_unique_batch_once(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "input.3lines"
            output = Path(directory) / "output.tsv"
            source.write_text(">p\nAA\nSO\n>q\nA\nI\n")
            with patch.object(self.script, "unique_batch", return_value="batch-123") as batch:
                self.assertEqual(self.script.main([str(source), str(output)]), 0)
                batch.assert_called_once_with()
            with output.open() as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(len(rows), 3)
            self.assertTrue(all(r["target_model"] == "" for r in rows))
            self.assertEqual({r["batch"] for r in rows}, {"batch-123"})
            self.assertEqual(rows[-1]["query_start"], "1")
            self.assertEqual(rows[-1]["target_end"], "1")

    def test_empty_input_writes_header(self):
        self.assertEqual(self.rows("\n"), [])
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "input"
            output = Path(directory) / "output"
            source.write_text("")
            self.script.main([str(source), str(output)])
            with output.open() as stream:
                reader = csv.DictReader(stream, delimiter="\t")
                self.assertIn("query_accession", reader.fieldnames)
                self.assertEqual(list(reader), [])

    def test_invalid_later_record_preserves_output(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "input"
            output = Path(directory) / "output"
            source.write_text(">p\nA\nS\n>q\nAA\nS\n")
            output.write_text("existing output")
            with self.assertRaises(ValueError):
                self.script.main([str(source), str(output)])
            self.assertEqual(output.read_text(), "existing output")


if __name__ == "__main__":
    unittest.main()
