import json
import os
import sys
import tempfile
import unittest

from sieve.channel_set import ChannelSet, load_channel_set
from sieve.channels import (
    ChannelBackground,
    CompositionBiasChannel,
    NetChargeChannel,
    ShortMotifChannel,
)


class TestChannelSet(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tempdir.cleanup()

    def path(self, name):
        return os.path.join(self.tempdir.name, name)

    def write(self, name, content):
        path = self.path(name)
        with open(path, "w", encoding="utf-8") as output:
            output.write(content)
        return path

    def channel_set(self):
        return ChannelSet([
            CompositionBiasChannel("alanine", radius=0, residues="A"),
            ShortMotifChannel("ac_motif", pattern="AC"),
        ])

    def test_preserves_channel_order_and_rejects_invalid_members(self):
        channel_set = self.channel_set()

        self.assertEqual(
            [channel.short_name for channel in channel_set],
            ["alanine", "ac_motif"],
        )
        with self.assertRaisesRegex(ValueError, "at least one"):
            ChannelSet([])
        with self.assertRaisesRegex(TypeError, "Channel instances"):
            ChannelSet([object()])

    def test_rejects_duplicate_and_reserved_short_names(self):
        with self.assertRaisesRegex(ValueError, "duplicate"):
            ChannelSet([
                NetChargeChannel("charge", 1),
                NetChargeChannel("charge", 2),
            ])
        with self.assertRaisesRegex(ValueError, "reserved"):
            ChannelSet([NetChargeChannel("position", 1)])

    def test_compute_backgrounds_includes_none_for_motif(self):
        fasta_path = self.write("proteome.faa", ">one\nA\n>two\nCC\n")

        backgrounds = self.channel_set().compute_backgrounds(fasta_path)

        self.assertIsInstance(backgrounds["alanine"], ChannelBackground)
        self.assertIsNone(backgrounds["ac_motif"])

    def test_background_file_round_trip_and_functions(self):
        fasta_path = self.write("proteome.faa", ">one\nA\n>two\nCC\n")
        background_path = self.path("backgrounds.json")
        channel_set = self.channel_set()
        computed = channel_set.compute_backgrounds(fasta_path)

        channel_set.save_backgrounds(computed, background_path)
        loaded = channel_set.load_backgrounds(background_path)
        functions = channel_set.make_functions(loaded)

        self.assertEqual(set(functions), {"alanine", "ac_motif"})
        self.assertAlmostEqual(functions["alanine"]("A")[0], 2 ** 0.5)
        self.assertEqual(functions["ac_motif"]("AC"), [1, 1])
        with open(background_path, encoding="utf-8") as background_file:
            document = json.load(background_file)
        self.assertEqual(document["file_type"], "sieve-channel-backgrounds")
        self.assertIsNone(document["channels"][1]["background"])

    def test_background_file_rejects_different_channel_arguments(self):
        fasta_path = self.write("proteome.faa", ">one\nA\n>two\nCC\n")
        background_path = self.path("backgrounds.json")
        channel_set = self.channel_set()
        channel_set.save_backgrounds(
            channel_set.compute_backgrounds(fasta_path),
            background_path,
        )
        changed = ChannelSet([
            CompositionBiasChannel("alanine", radius=1, residues="A"),
            ShortMotifChannel("ac_motif", pattern="AC"),
        ])

        with self.assertRaisesRegex(ValueError, "does not match"):
            changed.load_backgrounds(background_path)

    def test_make_functions_requires_exact_background_keys(self):
        channel_set = self.channel_set()

        with self.assertRaisesRegex(ValueError, "missing"):
            channel_set.make_functions({"alanine": ChannelBackground(0, 1)})

    def test_load_channel_set_uses_module_attribute_path(self):
        self.write(
            "example_channels.py",
            "\n".join([
                "from sieve.channel_set import ChannelSet",
                "from sieve.channels import ShortMotifChannel",
                "CHANNELS = ChannelSet([ShortMotifChannel('motif', 'AC')])",
                "not_channels = object()",
                "",
            ]),
        )
        sys.path.insert(0, self.tempdir.name)
        try:
            loaded = load_channel_set("example_channels.CHANNELS")
            self.assertIsInstance(loaded, ChannelSet)
            with self.assertRaisesRegex(TypeError, "ChannelSet instance"):
                load_channel_set("example_channels.not_channels")
            with self.assertRaisesRegex(ValueError, "Cannot find"):
                load_channel_set("example_channels.missing")
        finally:
            sys.path.remove(self.tempdir.name)
            sys.modules.pop("example_channels", None)

    def test_load_channel_set_rejects_spec_without_attribute(self):
        with self.assertRaisesRegex(ValueError, "module.attribute"):
            load_channel_set("channels")
