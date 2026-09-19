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

    def test_dna_gc_flag_combinations_both_strands(self):
        self.manifest(["G1", "G2"])
        for genome, strand, first, second in [
            ("G1", "+", "atggcta", "aactgact"),
            ("G2", "-", "tagccat", "agtcagtt"),
        ]:
            self.write_genome(genome, strand)
            contig = list("N" * 1000)
            for row, dna in zip(self.rows([7, 8], strand), [first, second]):
                contig[row["start"] - 1:row["end"]] = dna
            self.fixture.write_genomic_fasta(genome, {"NC_1": "".join(contig)})
        for flags, suffixes in [
            ([], ["", ""]),
            (["--dna"], [", ATGGCTA", ", AACTGACT"]),
            (["--gc"], [", 42.86%", ", 37.50%"]),
            (["--dna", "--gc"], [", ATGGCTA, 42.86%", ", AACTGACT, 37.50%"]),
        ]:
            with self.subTest(flags=flags):
                output = io.StringIO()
                with redirect_stdout(output):
                    self.script.main(["P1", *flags])
                lines = output.getvalue().splitlines()
                for offset in (0, 3):
                    self.assertTrue(lines[offset + 1].endswith(": MA(K)" + suffixes[0]))
                    self.assertTrue(lines[offset + 2].endswith(": (K)LT" + suffixes[1]))

    def test_gc_ambiguous_bases_and_extremes(self):
        self.manifest(["G1"])
        self.write_genome("G1")
        for dna, expected in [("gcnnnnn", "28.57%"), ("AAAAAAA", "0.00%"), ("GCGCGCG", "100.00%")]:
            self.fixture.write_genomic_fasta("G1", {"NC_1": "N" * 9 + dna + "N" * 20})
            output = io.StringIO()
            with redirect_stdout(output):
                self.script.main(["P1", "--gc"])
            self.assertTrue(output.getvalue().splitlines()[1].endswith(", " + expected))

    def test_dna_gc_require_genomic_file_contig_and_valid_bounds(self):
        self.manifest(["G1"])
        self.write_genome("G1")
        for flag in ("--dna", "--gc"):
            with self.subTest(flag=flag), self.assertRaises(FileNotFoundError):
                self.script.main(["P1", flag])
        for sequences, message in [({"OTHER": "N" * 40}, "contig sequence"), ({"NC_1": "N" * 15}, "exceed")]:
            self.fixture.write_genomic_fasta("G1", sequences)
            for flag in ("--dna", "--gc"):
                output = io.StringIO()
                with redirect_stdout(output), self.assertRaisesRegex(ValueError, message):
                    self.script.main(["P1", flag])
                self.assertEqual(output.getvalue(), "")

    def test_fasta_protein_suppresses_exons_and_ignores_gc(self):
        self.manifest(["G2", "G1"])
        self.write_genome("G1")
        self.write_genome("G2", "-")
        for flags in (["--fasta"], ["--fasta", "--gc"]):
            output = io.StringIO()
            with redirect_stdout(output):
                self.script.main(["P1", *flags])
            self.assertEqual(output.getvalue(), ">P1\nMAKLT\n>P1\nMAKLT\n")

    def test_fasta_dna_splices_exons_in_translation_order(self):
        self.manifest(["G1", "G2"])
        for genome, strand, first, second in [
            ("G1", "+", "atggcta", "aactgact"),
            ("G2", "-", "tagccat", "agtcagtt"),
        ]:
            self.write_genome(genome, strand)
            contig = list("N" * 1000)
            for row, dna in zip(self.rows([7, 8], strand), [first, second]):
                contig[row["start"] - 1:row["end"]] = dna
            self.fixture.write_genomic_fasta(genome, {"NC_1": "".join(contig)})
        for flags in (["--fasta", "--dna"], ["--fasta", "--dna", "--gc"]):
            output = io.StringIO()
            with redirect_stdout(output):
                self.script.main(["P1", *flags])
            self.assertEqual(
                output.getvalue(),
                ">P1_cds\nATGGCTAAACTGACT\n>P1_cds\nATGGCTAAACTGACT\n",
            )

    def test_fasta_dna_validates_genomic_input_before_stdout(self):
        self.manifest(["G1", "G2"])
        self.write_genome("G1")
        self.write_genome("G2")
        self.fixture.write_genomic_fasta("G1", {"NC_1": "N" * 1000})
        output = io.StringIO()
        with redirect_stdout(output), self.assertRaises(FileNotFoundError):
            self.script.main(["P1", "--fasta", "--dna"])
        self.assertEqual(output.getvalue(), "")

    def write_single_cds_genome(self, genome, accession, protein_sequence, dna):
        self.fixture.write_ncbi_proteins(genome, {accession: protein_sequence})
        self.fixture.write_gff(
            genome,
            f"##gff-version 3\nNC_1\tNCBI\tCDS\t1\t{len(dna)}\t.\t+\t0\t"
            f"ID=cds-{accession};protein_id={accession}\n",
        )
        self.fixture.write_genomic_fasta(genome, {"NC_1": dna})

    def test_multiple_accessions_preserve_order_in_all_output_modes(self):
        self.fixture.write_manifest([
            dict(sequence_accession="P1", sequence_database="G1", sequence_type="protein",
                 sequence_source="ncbi", sequence_length=2),
            dict(sequence_accession="P2", sequence_database="G2", sequence_type="protein",
                 sequence_source="ncbi", sequence_length=2),
        ])
        self.write_single_cds_genome("G1", "P1", "MA", "ATGGCT")
        self.write_single_cds_genome("G2", "P2", "GP", "GGTCCT")

        output = io.StringIO()
        with redirect_stdout(output):
            self.script.main(["P2", "P1"])
        self.assertEqual(
            output.getvalue(),
            "genome G2\nexon 1, NC_1, 1, 6: GP\n"
            "genome G1\nexon 1, NC_1, 1, 6: MA\n",
        )

        for flags, expected in [
            (["--fasta"], ">P2\nGP\n>P1\nMA\n"),
            (["--fasta", "--dna"], ">P2_cds\nGGTCCT\n>P1_cds\nATGGCT\n"),
        ]:
            output = io.StringIO()
            with redirect_stdout(output):
                self.script.main(["P2", "P1", *flags])
            self.assertEqual(output.getvalue(), expected)

    def test_multiple_accessions_repeat_duplicates(self):
        self.fixture.write_manifest([
            dict(sequence_accession="P1", sequence_database="G1", sequence_type="protein",
                 sequence_source="ncbi", sequence_length=2),
        ])
        self.write_single_cds_genome("G1", "P1", "MA", "ATGGCT")
        output = io.StringIO()
        with redirect_stdout(output):
            self.script.main(["P1", "P1", "--fasta"])
        self.assertEqual(output.getvalue(), ">P1\nMA\n>P1\nMA\n")

    def test_multiple_accessions_have_no_partial_output_on_failure(self):
        self.fixture.write_manifest([
            dict(sequence_accession="P1", sequence_database="G1", sequence_type="protein",
                 sequence_source="ncbi", sequence_length=2),
        ])
        self.write_single_cds_genome("G1", "P1", "MA", "ATGGCT")
        output = io.StringIO()
        with redirect_stdout(output), self.assertRaisesRegex(ValueError, "OTHER.*manifest"):
            self.script.main(["P1", "OTHER", "--fasta"])
        self.assertEqual(output.getvalue(), "")


if __name__ == "__main__":
    unittest.main()
