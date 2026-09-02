import csv
import json
import math
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from io import StringIO

from tests.scripts.helpers import load_script


SCRIPTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "scripts")
)


class TestComputeChannelsScripts(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.background_script = load_script(
            os.path.join(SCRIPTS_DIR, "compute-channel-backgrounds.py")
        )
        self.channel_script = load_script(
            os.path.join(SCRIPTS_DIR, "compute-protein-channels.py")
        )
        self.write(
            "example_channels.py",
            "\n".join([
                "from sieve.channel_set import ChannelSet",
                "from sieve.channels import CompositionBiasChannel, ShortMotifChannel",
                "CHANNELS = ChannelSet([",
                "    CompositionBiasChannel('alanine', radius=0, residues='A'),",
                "    ShortMotifChannel('ac_motif', pattern='AC'),",
                "])",
                "",
            ]),
        )
        sys.path.insert(0, self.tempdir.name)

    def tearDown(self):
        sys.path.remove(self.tempdir.name)
        sys.modules.pop("example_channels", None)
        self.tempdir.cleanup()

    def path(self, name):
        return os.path.join(self.tempdir.name, name)

    def write(self, name, content):
        path = self.path(name)
        with open(path, "w", encoding="utf-8") as output:
            output.write(content)
        return path

    def test_scripts_compute_backgrounds_then_position_wise_tsv(self):
        proteome_path = self.write(
            "proteome.faa",
            ">one\nA\n>invalid_background\nAX\n>two\nCC\n",
        )
        proteins_path = self.write(
            "proteins.faa",
            ">protein_1 description\nac\n"
            ">invalid_protein\nAX\n"
            ">protein_2\nC\n",
        )
        background_path = self.path("backgrounds.json")
        output_path = self.path("channels.tsv")

        background_stderr = StringIO()
        with redirect_stderr(background_stderr):
            result = self.background_script.main([
                proteome_path,
                "-c",
                "example_channels.CHANNELS",
                "-o",
                background_path,
            ])
        self.assertEqual(result, 0)
        channel_stderr = StringIO()
        with redirect_stderr(channel_stderr):
            result = self.channel_script.main([
                proteins_path,
                "-c",
                "example_channels.CHANNELS",
                "-b",
                background_path,
                "-o",
                output_path,
            ])
        self.assertEqual(result, 0)

        self.assertEqual(
            background_stderr.getvalue().splitlines(),
            [
                "Processing one",
                "Skipping invalid_background: sequence contains noncanonical amino acids: X",
                "Processing two",
            ],
        )
        self.assertEqual(
            channel_stderr.getvalue().splitlines(),
            [
                "Processing protein_1",
                "Skipping invalid_protein: sequence contains noncanonical amino acids: X",
                "Processing protein_2",
            ],
        )

        with open(background_path, encoding="utf-8") as background_file:
            backgrounds = json.load(background_file)
        self.assertEqual(backgrounds["channels"][0]["definition"]["short_name"], "alanine")
        with open(output_path, encoding="utf-8", newline="") as output:
            rows = list(csv.DictReader(output, delimiter="\t"))
        self.assertEqual(
            [(row["sequence_id"], row["position"], row["residue"]) for row in rows],
            [
                ("protein_1", "1", "A"),
                ("protein_1", "2", "C"),
                ("protein_2", "1", "C"),
            ],
        )
        self.assertAlmostEqual(float(rows[0]["alanine"]), math.sqrt(2))
        self.assertAlmostEqual(float(rows[1]["alanine"]), -1 / math.sqrt(2))
        self.assertEqual([row["ac_motif"] for row in rows], ["1", "1", "0"])

    def test_channel_tsv_rejects_duplicate_first_header_tokens(self):
        proteome_path = self.write("proteome.faa", ">one\nA\n>two\nCC\n")
        proteins_path = self.write(
            "proteins.faa",
            ">duplicate first\nA\n>duplicate second\nC\n",
        )
        background_path = self.path("backgrounds.json")
        with redirect_stderr(StringIO()):
            self.background_script.main([
                proteome_path,
                "-c",
                "example_channels.CHANNELS",
                "-o",
                background_path,
            ])

        with redirect_stderr(StringIO()):
            with self.assertRaisesRegex(ValueError, "duplicate FASTA sequence ID"):
                self.channel_script.main([
                    proteins_path,
                    "-c",
                    "example_channels.CHANNELS",
                    "-b",
                    background_path,
                    "-o",
                    self.path("channels.tsv"),
                ])
        self.assertFalse(os.path.exists(self.path("channels.tsv")))
