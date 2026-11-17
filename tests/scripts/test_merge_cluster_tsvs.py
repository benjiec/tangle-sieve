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

    def test_compares_sequences_by_cluster_and_member(self):
        a, b, c = (f"{name}_PF00531_1_1" for name in "abc")
        table = self.input("clusters.tsv", f"{HEADER}\n{a}\t{a}\n{a}\t{c}\n{b}\t{c}\n")
        first = self.input("first.fasta", f">{a}\n>{c}\nABCD\n>{a}\nAB\nCD\n")
        second = self.input("second.fasta", f">{b}\n>{b}\nABCd\n>{c}\nABCD\n")
        result = self.merge(table, "--cluster-fasta", first, second, first)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            self.output.read_text(),
            f"{HEADER}\tProtein Accession\tSame as Cluster Repr\n"
            f"{a}\t{a}\ta\ttrue\n{a}\t{c}\tc\ttrue\n{b}\t{c}\tc\tfalse\n",
        )

    def test_invalid_fasta_keeps_output_intact(self):
        a, b = (f"{name}_PF00531_1_1" for name in "ab")
        table = self.input("clusters.tsv", f"{HEADER}\n{a}\t{b}\n")
        valid = self.input("valid.fasta", f">{a}\n>{a}\nABC\n>{b}\nABC\n")
        for body, expected, other_files in (
            ("ABC\n", "broken.fasta:1:", []),
            (">\n", "broken.fasta:1:", []),
            (f">{a}\nABC\n", "before a cluster header", []),
            (f">{a}\n>{b}\nABC\n", "missing sequence for representative", []),
            (f">{a}\n>{a}\nABC\n", "clusters.tsv:2: missing FASTA sequences", []),
            (f">{a}\n>{a}\nABD\n>{b}\nABD\n", "conflicting sequences", [valid]),
            (f">{a}\n>{a}\nABC\n>{b}\nABC\n>{b}\nABD\n", "conflicting sequences", []),
            (f">{a}\n>{a}\nABC\n>{b}\nABC\n>{b}\n", "missing sequence for representative", []),
        ):
            with self.subTest(body=body):
                broken = self.input("broken.fasta", body)
                self.output.write_text("existing output")
                result = self.merge(table, "--cluster-fasta", *other_files, broken)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(expected, result.stderr)
                self.assertEqual(self.output.read_text(), "existing output")

    def test_rejects_fasta_as_output(self):
        table = self.input("clusters.tsv", HEADER)
        fasta = self.input("input.fasta", ">rep\n>rep\nABC\n")
        result = self.merge(table, "--cluster-fasta", fasta, output=fasta)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("output refers to input file", result.stderr)
        self.assertEqual(fasta.read_text(), ">rep\n>rep\nABC\n")
