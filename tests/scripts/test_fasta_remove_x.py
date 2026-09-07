import os
import tempfile
import unittest

from tests.scripts.helpers import load_script


class TestFastaRemoveXScript(unittest.TestCase):

    def setUp(self):
        repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.script = load_script(os.path.join(repo, "scripts", "fasta-remove-x.py"))

    def run_script(self, fasta_text):
        with tempfile.TemporaryDirectory() as tmpd:
            input_fasta = os.path.join(tmpd, "input.faa")
            output_fasta = os.path.join(tmpd, "without-x.faa")
            with open(input_fasta, "w", encoding="utf-8") as stream:
                stream.write(fasta_text)
            self.assertEqual(self.script.main([input_fasta, output_fasta]), 0)
            with open(output_fasta, encoding="utf-8") as stream:
                return stream.read()

    def test_removes_sequences_containing_x(self):
        self.assertEqual(
            self.run_script(
                ">keep description\nMSE\nQ\n"
                ">remove_start\nXAAA\n"
                ">remove_middle\nAAXA\n"
                ">remove_end\nAAAX\n"
                ">also_keep\nBBBB\n"
            ),
            ">keep description\nMSEQ\n>also_keep\nBBBB\n",
        )

    def test_lowercase_x_is_not_removed(self):
        self.assertEqual(
            self.run_script(">lowercase description\nAAxA\n"),
            ">lowercase description\nAAxA\n",
        )

    def test_writes_empty_fasta_for_empty_input(self):
        self.assertEqual(self.run_script(""), "")


if __name__ == "__main__":
    unittest.main()
