import io
from pathlib import Path
import unittest
from contextlib import redirect_stdout

from sieve.protein import CuratedProtein
from tests.fixtures import DefaultsFixture
from tests.scripts.helpers import load_script


class TestProteinExons(unittest.TestCase):
    def setUp(self):
        self.script = load_script(Path(__file__).resolve().parents[2] / "scripts" / "protein-exons.py")
        self.fixture = DefaultsFixture(self)
        self.addCleanup(self.fixture.cleanup)
        self.addCleanup(CuratedProtein.clear_cache)

    def rows(self, lengths, strand="+", initial_phase=0):
        rows = []
        offset = -initial_phase
        position = 10
        for index, length in enumerate(lengths):
            rows.append(dict(seqid="NC_1", strand=strand, start=position,
                             end=position + length - 1,
                             phase=initial_phase if index == 0 else (-offset) % 3,
                             attrs={"ID": "cds-P1", "protein_id": "P1"}))
            offset += length
            position += length + 10
        if strand == "-":
            for row in rows:
                row["start"], row["end"] = 1000 - row["end"], 1000 - row["start"]
        return rows

    def test_boundary_offsets_and_short_exons(self):
        cases = [
            ([6, 9], ["MA", "KLT"]),
            ([7, 8], ["MA(K)", "(K)LT"]),
            ([8, 7], ["MA(K)", "(K)LT"]),
            ([7, 3, 5], ["MA(K)", "(K)(L)", "(L)T"]),
            ([7, 1, 7], ["MA(K)", "(K)", "(K)LT"]),
        ]
        for strand in ("+", "-"):
            for lengths, expected in cases:
                with self.subTest(strand=strand, lengths=lengths):
                    portions = self.script.exon_portions(self.rows(lengths, strand), "MAKLT")
                    self.assertEqual([text for _, text in portions], expected)

    def test_partial_start_stop_and_partial_end(self):
        for phase in (0, 1, 2):
            for trailing in range(4):
                with self.subTest(phase=phase, trailing=trailing):
                    portions = self.script.exon_portions(self.rows([6 + phase + trailing], initial_phase=phase), "MA*")
                    self.assertEqual(portions[0][1], "MA")
        portions = self.script.exon_portions(self.rows([6, 3]), "MA")
        self.assertEqual([text for _, text in portions], ["MA", ""])

    def test_invalid_annotations(self):
        mutations = [
            lambda rows: rows[1].update(phase=0),
            lambda rows: rows[1].update(seqid="other"),
            lambda rows: rows[1].update(strand="-"),
            lambda rows: rows[1].update(start=11),
            lambda rows: rows[0].update(start=0),
            lambda rows: rows[1]["attrs"].update(ID="other"),
        ]
        for mutate in mutations:
            rows = self.rows([7, 8])
            mutate(rows)
            with self.subTest(mutate=mutate), self.assertRaises(ValueError):
                self.script.exon_portions(rows, "MAKLT")
        for sequence in ("M", "MAKLTT"):
            with self.assertRaisesRegex(ValueError, "length"):
                self.script.exon_portions(self.rows([15]), sequence)

    def write_genome(self, genome, strand="+"):
        self.fixture.write_ncbi_proteins(genome, {"P1": "MAKLT"})
        lines = ["##gff-version 3"]
        for row in reversed(self.rows([7, 8], strand)):
            lines.append(f"NC_1\tNCBI\tCDS\t{row['start']}\t{row['end']}\t.\t{strand}\t{row['phase']}\tID=cds-P1;protein_id=P1")
        lines.append("NC_1\tNCBI\tCDS\t2000\t2002\t.\t+\t0\tprotein_id=OTHER")
        self.fixture.write_gff(genome, "\n".join(lines) + "\n##FASTA\n>NC_1\nAAA\n")

    def manifest(self, genomes):
        self.fixture.write_manifest([
            dict(sequence_accession="P1", sequence_database=genome,
                 sequence_type="protein", sequence_source="ncbi", sequence_length=5)
            for genome in genomes
        ])

    def test_multiple_genomes_stdout(self):
        self.manifest(["G2", "G1"])
        self.write_genome("G1")
        self.write_genome("G2", "-")
        output = io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(self.script.main(["P1"]), 0)
        self.assertEqual(output.getvalue(),
                         "genome G1\nexon 1, NC_1, 10, 16: MA(K)\n"
                         "exon 2, NC_1, 27, 34: (K)LT\n"
                         "genome G2\nexon 1, NC_1, 984, 990: MA(K)\n"
                         "exon 2, NC_1, 966, 973: (K)LT\n")

    def test_missing_manifest_entry_and_file(self):
        with self.assertRaises(FileNotFoundError):
            self.script.main(["P1"])
        self.manifest(["G1"])
        with self.assertRaisesRegex(ValueError, "manifest"):
            self.script.main(["OTHER"])
        with self.assertRaises(FileNotFoundError):
            self.script.main(["P1"])

    def test_missing_cds_invalid_phase_and_protein(self):
        self.manifest(["G1"])
        self.fixture.write_gff("G1", "##gff-version 3\n")
        with self.assertRaisesRegex(ValueError, "CDS rows"):
            self.script.main(["P1"])
        self.fixture.write_gff("G1", "NC_1\tNCBI\tCDS\t1\t15\t.\t+\t.\tprotein_id=P1\n")
        with self.assertRaisesRegex(ValueError, "phase"):
            self.script.main(["P1"])
        self.write_genome("G1")
        self.fixture.write_ncbi_proteins("G1", {"OTHER": "MAKLT"})
        with self.assertRaisesRegex(ValueError, "protein sequence"):
            self.script.main(["P1"])

    def test_no_partial_stdout_on_failure(self):
        self.manifest(["G1", "G2"])
        self.write_genome("G1")
        output = io.StringIO()
        with redirect_stdout(output), self.assertRaises(FileNotFoundError):
            self.script.main(["P1"])
        self.assertEqual(output.getvalue(), "")


if __name__ == "__main__":
    unittest.main()
