from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

from sieve.region_accession import format_region_accession, parse_protein_accession


SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "merge-cluster-tsvs.py"
HEADER = "Cluster Repr\tCluster Member"


class TestRegionAccessions(unittest.TestCase):
    def test_round_trip(self):
        for prefix in ("p", "p_1", "_p_a_b_c_", "p_PF00531_1_1", "p_PF00531_1_1_domain_1"):
            for pfam in ("PF00531", "PF00531.10"):
                for domain in (None, 1, 12):
                    with self.subTest(prefix=prefix, pfam=pfam, domain=domain):
                        identifier = format_region_accession(prefix, pfam, 2, 12, domain)
                        self.assertEqual(parse_protein_accession(identifier), prefix)

    def test_invalid_suffixes(self):
        for identifier in (
            "p", "_PF00531_1_1", "p_PF0053_1_1", "p_XX00531_1_1",
            "p_PF00531.bad_1_1", "p_PF00531_1", "p_PF00531_01_1",
            "p_PF00531_0_1", "p_PF00531_2_1", "p_PF00531_1_1_domain_",
            "p_PF00531_1_1_domain_0", "p_PF00531_1_1_domain_2", "p_PF00531_1_1 ",
            "p x_PF00531_1_1",
        ):
            with self.subTest(identifier=identifier), self.assertRaises(ValueError):
                parse_protein_accession(identifier)


class TestMergeClusterTsvs(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.output = self.root / "output.tsv"

    def input(self, name, body):
        path = self.root / name
        path.write_bytes(body.encode("utf-8"))
        return path

    def merge(self, *inputs, output=None):
        return subprocess.run(
            [sys.executable, SCRIPT, "--output", output or self.output, *inputs],
            capture_output=True, text=True, check=False,
        )

    def test_combines_inputs_in_order(self):
        row = "repr\tp_a_b_c_PF00531_1_1"
        first = self.input("first.tsv", f"{HEADER}\r\n{row}\r\n{row}\r\n")
        empty = self.input("empty.tsv", HEADER)
        second = self.input("second.tsv", f"{HEADER}\nrepr\tother_PF00531_1_2_domain_2")
        result = self.merge(first, empty, second)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            self.output.read_text(),
            f"{HEADER}\tProtein Accession\n{row}\tp_a_b_c\n{row}\tp_a_b_c\n"
            "repr\tother_PF00531_1_2_domain_2\tother\n",
        )

    def test_invalid_later_input_leaves_output_intact(self):
        first = self.input("first.tsv", f"{HEADER}\nr\tp_PF00531_1_1\n")
        for body, line in (("wrong\n", 1), (f"{HEADER}\nr\n", 2), (f"{HEADER}\nr\tp\n", 2)):
            with self.subTest(body=body):
                later = self.input("later.tsv", body)
                self.output.write_text("existing output")
                result = self.merge(first, later)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(f"{later}:{line}:", result.stderr)
                self.assertEqual(self.output.read_text(), "existing output")
                self.assertEqual(set(self.root.iterdir()), {first, later, self.output})

    def test_rejects_input_as_output(self):
        source = self.input("source.tsv", HEADER)
        alias = self.root / "alias.tsv"
        alias.symlink_to(source)
        for output in (source, alias):
            with self.subTest(output=output):
                result = self.merge(source, output=output)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("output refers to input file", result.stderr)
                self.assertEqual(source.read_text(), HEADER)
