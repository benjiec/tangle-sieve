import os
import tempfile
import unittest

from tests.scripts.helpers import load_script


class TestFastaFilterLengthScript(unittest.TestCase):

    def setUp(self):
        repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.script = load_script(os.path.join(repo, "scripts", "fasta-filter-length.py"))

    def run_script(self, fasta_text, minimum, maximum):
        with tempfile.TemporaryDirectory() as tmpd:
            input_fasta = os.path.join(tmpd, "input.faa")
            output_fasta = os.path.join(tmpd, "filtered.faa")
            with open(input_fasta, "w", encoding="utf-8") as stream:
                stream.write(fasta_text)
            self.assertEqual(
                self.script.main([
                    input_fasta,
                    str(minimum),
                    str(maximum),
                    output_fasta,
                ]),
                0,
            )
            with open(output_fasta, encoding="utf-8") as stream:
                return stream.read()

    def test_keeps_sequences_at_inclusive_limits(self):
        self.assertEqual(
            self.run_script(
                ">too_short\nAA\n"
                ">minimum full description\nAAA\n"
                ">middle\nAA\nAA\n"
                ">maximum\nAAAAA\n"
                ">too_long\nAAAAAA\n",
                3,
                5,
            ),
            ">minimum full description\nAAA\n>middle\nAAAA\n>maximum\nAAAAA\n",
        )

    def test_writes_empty_fasta_when_no_sequences_match(self):
        self.assertEqual(self.run_script(">short\nAA\n", 3, 5), "")

    def test_rejects_invalid_limits(self):
        with self.assertRaises(ValueError):
            self.script.filter_by_length("input.faa", 0, 5, "output.faa")
        with self.assertRaises(ValueError):
            self.script.filter_by_length("input.faa", 6, 5, "output.faa")
        with self.assertRaises(SystemExit):
            self.script.main(["input.faa", "6", "5", "output.faa"])
        with self.assertRaises(SystemExit):
            self.script.main(["input.faa", "0", "5", "output.faa"])


if __name__ == "__main__":
    unittest.main()
