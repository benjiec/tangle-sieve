import csv
import os
import tempfile
import unittest

from tests.scripts.helpers import load_script


SCRIPT = load_script(os.path.join(
    os.path.dirname(__file__), "..", "..", "scripts", "alphafold-pdockq2-directory.py"
))


class AlphaFoldPDockQ2DirectoryTests(unittest.TestCase):
    def test_sha256_file_depends_on_contents(self):
        with tempfile.NamedTemporaryFile("wb") as stream:
            stream.write(b"alpha")
            stream.flush()
            first = SCRIPT.sha256_file(stream.name)
            stream.seek(0)
            stream.truncate()
            stream.write(b"beta")
            stream.flush()
            self.assertNotEqual(first, SCRIPT.sha256_file(stream.name))

    def test_cache_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "results.tsv")
            columns = [*SCRIPT.FIXED_COLUMNS, "A_B_pDockQ2_max"]
            rows = [{
                "zipfile": "fold.zip",
                "zipfile_sha256": "abc",
                "model": "0",
                "A_B_pDockQ2_max": "0.500000",
            }]
            SCRIPT.write_cache(path, columns, rows)
            self.assertEqual(SCRIPT.read_cache(path), (columns, rows))

    def test_cache_requires_fixed_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "results.tsv")
            with open(path, "w", encoding="utf-8", newline="") as stream:
                csv.writer(stream, delimiter="\t").writerow(["zipfile", "model"])
            with self.assertRaisesRegex(ValueError, "invalid header"):
                SCRIPT.read_cache(path)

    def test_cache_write_replaces_existing_file(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "results.tsv")
            with open(path, "w", encoding="utf-8") as stream:
                stream.write("old\n")
            columns = list(SCRIPT.FIXED_COLUMNS)
            SCRIPT.write_cache(path, columns, [{
                "zipfile": "fold.zip", "zipfile_sha256": "abc", "model": "0"
            }])
            read_columns, rows = SCRIPT.read_cache(path)
            self.assertEqual(read_columns, columns)
            self.assertEqual(rows[0]["zipfile_sha256"], "abc")


if __name__ == "__main__":
    unittest.main()
