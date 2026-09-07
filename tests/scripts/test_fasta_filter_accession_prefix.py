import os
import tempfile
import unittest

from tests.scripts.helpers import load_script


class TestFastaFilterAccessionPrefixScript(unittest.TestCase):

    def setUp(self):
        repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.script = load_script(
            os.path.join(repo, "scripts", "fasta-filter-accession-prefix.py")
        )

    def run_script(self, big_text, filter_text):
        with tempfile.TemporaryDirectory() as tmpd:
            big_fasta = os.path.join(tmpd, "big.faa")
            filter_fasta = os.path.join(tmpd, "filter.faa")
            output_fasta = os.path.join(tmpd, "output.faa")
            with open(big_fasta, "w", encoding="utf-8") as stream:
                stream.write(big_text)
            with open(filter_fasta, "w", encoding="utf-8") as stream:
                stream.write(filter_text)

            self.assertEqual(
                self.script.main([big_fasta, filter_fasta, output_fasta]),
                0,
            )
            with open(output_fasta, encoding="utf-8") as stream:
                return stream.read()

    def test_retains_exact_and_prefix_accessions_in_big_fasta_order(self):
        self.assertEqual(
            self.run_script(
                ">ABC description retained\nMS\nEQ\n"
                ">EXACT full header\nAAAA\n"
                ">ABCD is not a prefix of ABC_suffix\nBBBB\n"
                ">MISSING\nCCCC\n",
                ">ABC_suffix filter description\nX\n"
                ">EXACT\nY\n",
            ),
            ">ABC description retained\nMSEQ\n>EXACT full header\nAAAA\n",
        )

    def test_does_not_match_substrings_or_filter_descriptions(self):
        self.assertEqual(
            self.run_script(
                ">BC\nAAA\n>DESCRIPTION\nBBB\n",
                ">ABC DESCRIPTION\nX\n",
            ),
            "",
        )

    def test_empty_filter_writes_empty_output(self):
        self.assertEqual(self.run_script(">ABC\nMSEQ\n", ""), "")


if __name__ == "__main__":
    unittest.main()
