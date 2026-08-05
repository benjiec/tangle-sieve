import os
import tempfile
import unittest

from tests.scripts.helpers import load_script


class TestFastaUniqueSequencesScript(unittest.TestCase):

    def setUp(self):
        repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.script = load_script(os.path.join(repo, "scripts", "fasta-unique-sequences.py"))

    def run_script(self, fasta_text):
        with tempfile.TemporaryDirectory() as tmpd:
            input_fasta = os.path.join(tmpd, "input.faa")
            output_fasta = os.path.join(tmpd, "unique.faa")
            with open(input_fasta, "w", encoding="utf-8") as f:
                f.write(fasta_text)
            self.script.main([input_fasta, output_fasta])
            with open(output_fasta, encoding="utf-8") as f:
                return f.read()

    def test_keeps_first_accession_for_each_unique_sequence_in_input_order(self):
        self.assertEqual(
            self.run_script(
                ">first description\nMSE\nQ\n"
                ">second\nAAAA\n"
                ">duplicate other description\nMSEQ\n"
                ">third\nBBBB\n"
                ">second_duplicate\nAAAA\n"
            ),
            ">first\nMSEQ\n>second\nAAAA\n>third\nBBBB\n",
        )

    def test_sequence_comparison_is_case_sensitive(self):
        self.assertEqual(
            self.run_script(">upper\nMSEQ\n>lower\nmseq\n"),
            ">upper\nMSEQ\n>lower\nmseq\n",
        )

    def test_writes_empty_fasta_for_empty_input(self):
        self.assertEqual(self.run_script(""), "")
