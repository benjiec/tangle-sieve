import unittest
from subprocess import CompletedProcess
from unittest.mock import patch

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tangle.detected import DetectedTable
from sieve import protein as protein_module
from sieve.protein import (
    CuratedProtein,
    ProteinHMMAlignment,
    SEQUENCE_SOURCE_HMM_DETECTED,
    SEQUENCE_SOURCE_NCBI,
)
from tests.fixtures import DefaultsFixture


class TestCuratedProtein(unittest.TestCase):

    def setUp(self):
        CuratedProtein.clear_cache()
        self.fx = DefaultsFixture(self)

    def tearDown(self):
        CuratedProtein.clear_cache()
        self.fx.cleanup()

    def manifest_row(self, protein_accession, genome_accession, source):
        return dict(
            sequence_accession=protein_accession,
            sequence_database=genome_accession,
            sequence_type="protein",
            sequence_source=source,
        )

    def detected_row(self, protein_accession, genome_accession, contig_accession, q_start, q_end, t_start, t_end):
        return dict(
            detection_type="model",
            detection_method="hmm",
            batch="b1",
            query_accession=contig_accession,
            query_database=genome_accession,
            query_type="contig",
            target_accession=protein_accession,
            target_database=genome_accession,
            target_type="protein",
            target_model="HMM1",
            query_start=q_start,
            query_end=q_end,
            target_start=t_start,
            target_end=t_end,
            evalue=0.001,
            bitscore=10,
        )

    def test_hmm_detected_sequence_and_locus(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_HMM_DETECTED)])
        self.fx.write_detected_proteins("g1", {"p1": "PGP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "AAACCCGGGTTTAAACCCGGG"})
        self.fx.write_detected_rows("g1", [
            self.detected_row("p1", "g1", "ctg1", 4, 9, 1, 2),
            self.detected_row("p1", "g1", "ctg1", 13, 18, 3, 4),
        ])

        protein = CuratedProtein("p1", "g1")
        locus = protein.genomic_locus()

        self.assertEqual(protein.sequence(), "PGP")
        self.assertEqual(protein.sequence_source, SEQUENCE_SOURCE_HMM_DETECTED)
        self.assertEqual(locus.sequence(), "CCCGGGTTTAAACCC")
        self.assertEqual(locus.cds_intervals_1b, [(1, 6), (10, 15)])
        self.assertEqual(locus.dss_positions_1b(), [6])
        self.assertEqual(locus.ass_positions_1b(), [10])

    def test_hmm_detected_leader_extends_to_upstream_methionine(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_HMM_DETECTED)])
        self.fx.write_detected_proteins("g1", {"p1": "PG"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "ATGAAACCCGGGTTT"})
        self.fx.write_detected_rows("g1", [
            self.detected_row("p1", "g1", "ctg1", 7, 12, 1, 2),
        ])

        protein = CuratedProtein("p1", "g1")
        locus = protein.genomic_locus_with_leader()

        self.assertEqual(protein.sequence_with_leader(), "MKPG")
        self.assertEqual(locus.start_1b, 1)
        self.assertEqual(locus.end_1b, 12)
        self.assertEqual(locus.sequence(), "ATGAAACCCGGG")

    def test_hmm_detected_leader_uses_last_upstream_methionine(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_HMM_DETECTED)])
        self.fx.write_detected_proteins("g1", {"p1": "PG"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "ATGAAAATGCCCGGG"})
        self.fx.write_detected_rows("g1", [
            self.detected_row("p1", "g1", "ctg1", 10, 15, 1, 2),
        ])

        protein = CuratedProtein("p1", "g1")
        locus = protein.genomic_locus_with_leader()

        self.assertEqual(protein.sequence_with_leader(), "MPG")
        self.assertEqual(locus.start_1b, 7)
        self.assertEqual(locus.end_1b, 15)
        self.assertEqual(locus.sequence(), "ATGCCCGGG")

    def test_hmm_detected_leader_sequence_candidates_include_upstream_and_first_ten_protein_positions(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_HMM_DETECTED)])
        self.fx.write_detected_proteins("g1", {"p1": "KMMP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "NNNATGATGAAAATGATGCCC"})
        self.fx.write_detected_rows("g1", [
            self.detected_row("p1", "g1", "ctg1", 10, 21, 1, 4),
        ])

        candidates = CuratedProtein("p1", "g1").leader_sequence_candidates()

        self.assertEqual(
            [
                (candidate.accession, candidate.start_label, candidate.start_aa_1b, candidate.sequence)
                for candidate in candidates
            ],
            [
                ("p1_with_leader_u2_M", "u2", -1, "MMKMMP"),
                ("p1_with_leader_u1_M", "u1", 0, "MKMMP"),
                ("p1_with_leader_2_M", "2", 2, "MMP"),
                ("p1_with_leader_3_M", "3", 3, "MP"),
            ],
        )

    def test_ncbi_leader_sequence_candidates_include_first_ten_protein_positions(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_NCBI)])
        self.fx.write_ncbi_proteins("g1", {"p1": "MMMGAAAAAMM"})

        candidates = CuratedProtein("p1", "g1").leader_sequence_candidates()

        self.assertEqual(
            [
                (candidate.accession, candidate.start_label, candidate.start_aa_1b, candidate.sequence)
                for candidate in candidates
            ],
            [
                ("p1_with_leader_1_M", "1", 1, "MMMGAAAAAMM"),
                ("p1_with_leader_2_M", "2", 2, "MMGAAAAAMM"),
                ("p1_with_leader_3_M", "3", 3, "MGAAAAAMM"),
                ("p1_with_leader_10_M", "10", 10, "MM"),
            ],
        )

    def test_ncbi_forward_leader_sequence_candidates_include_upstream_tail(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_NCBI)])
        self.fx.write_ncbi_proteins("g1", {"p1": "MP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "ATGATGAAAATGCCC"})
        self.fx.write_gff("g1", "\n".join([
            "ctg1\tsrc\tmRNA\t10\t15\t.\t+\t.\tID=tx1",
            "ctg1\tsrc\tCDS\t10\t15\t.\t+\t0\tID=cds1;Parent=tx1;protein_id=p1",
            "ctg1\tsrc\tstart_codon\t10\t12\t.\t+\t0\tParent=tx1",
            "",
        ]))

        candidates = CuratedProtein("p1", "g1").leader_sequence_candidates()

        self.assertEqual(
            [
                (candidate.accession, candidate.start_label, candidate.start_aa_1b, candidate.sequence)
                for candidate in candidates
            ],
            [
                ("p1_with_leader_u3_M", "u3", -2, "MMKMP"),
                ("p1_with_leader_u2_M", "u2", -1, "MKMP"),
                ("p1_with_leader_1_M", "1", 1, "MP"),
            ],
        )

    def test_ncbi_reverse_leader_sequence_candidates_include_upstream_tail(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_NCBI)])
        self.fx.write_ncbi_proteins("g1", {"p1": "MP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "A" * 12 + "GGGCATTTTCATCAT"})
        self.fx.write_gff("g1", "\n".join([
            "ctg1\tsrc\tmRNA\t13\t27\t.\t-\t.\tID=tx1",
            "ctg1\tsrc\tCDS\t13\t18\t.\t-\t0\tID=cds1;Parent=tx1;protein_id=p1",
            "ctg1\tsrc\tstart_codon\t16\t18\t.\t-\t0\tParent=tx1",
            "",
        ]))

        candidates = CuratedProtein("p1", "g1").leader_sequence_candidates()

        self.assertEqual(
            [
                (candidate.accession, candidate.start_label, candidate.start_aa_1b, candidate.sequence)
                for candidate in candidates
            ],
            [
                ("p1_with_leader_u3_M", "u3", -2, "MMKMP"),
                ("p1_with_leader_u2_M", "u2", -1, "MKMP"),
                ("p1_with_leader_1_M", "1", 1, "MP"),
            ],
        )

    def test_ncbi_leader_sequence_candidates_stop_at_upstream_stop_codon_boundary(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_NCBI)])
        self.fx.write_ncbi_proteins("g1", {"p1": "P"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "ATGTAAATGAAACCC"})
        self.fx.write_gff("g1", "\n".join([
            "ctg1\tsrc\tmRNA\t13\t15\t.\t+\t.\tID=tx1",
            "ctg1\tsrc\tCDS\t13\t15\t.\t+\t0\tID=cds1;Parent=tx1;protein_id=p1",
            "",
        ]))

        candidates = CuratedProtein("p1", "g1").leader_sequence_candidates()

        self.assertEqual(
            [
                (candidate.accession, candidate.start_label, candidate.start_aa_1b, candidate.sequence)
                for candidate in candidates
            ],
            [
                ("p1_with_leader_u2_M", "u2", -1, "MKP"),
            ],
        )

    def test_leader_sequence_candidates_fall_back_to_original_sequence_without_m_in_first_ten_positions_or_upstream_tail(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_NCBI)])
        self.fx.write_ncbi_proteins("g1", {"p1": "RYDAAAAAAAM"})

        candidates = CuratedProtein("p1", "g1").leader_sequence_candidates()

        self.assertEqual(
            [
                (candidate.accession, candidate.start_label, candidate.start_aa_1b, candidate.sequence)
                for candidate in candidates
            ],
            [
                ("p1", "", 1, "RYDAAAAAAAM"),
            ],
        )

    def test_hmm_detected_leader_stops_at_stop_codon_boundary(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_HMM_DETECTED)])
        self.fx.write_detected_proteins("g1", {"p1": "P"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "ATGTAAAAACCC"})
        self.fx.write_detected_rows("g1", [
            self.detected_row("p1", "g1", "ctg1", 10, 12, 1, 1),
        ])

        protein = CuratedProtein("p1", "g1")

        self.assertEqual(protein.sequence_with_leader(), "KP")
        self.assertEqual(protein.genomic_locus_with_leader().sequence(), "AAACCC")

    def test_hmm_detected_leader_trims_to_second_amino_acid_methionine(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_HMM_DETECTED)])
        self.fx.write_detected_proteins("g1", {"p1": "KMP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "AAAATGCCC"})
        self.fx.write_detected_rows("g1", [
            self.detected_row("p1", "g1", "ctg1", 1, 9, 1, 3),
        ])

        protein = CuratedProtein("p1", "g1")
        locus = protein.genomic_locus_with_leader()

        self.assertEqual(protein.sequence_with_leader(), "MP")
        self.assertEqual(locus.start_1b, 4)
        self.assertEqual(locus.end_1b, 9)
        self.assertEqual(locus.sequence(), "ATGCCC")
        self.assertEqual(locus.cds_intervals_1b, [(1, 6)])
        self.assertEqual(locus.start_codon_interval_1b, (1, 3))

    def test_hmm_detected_leader_trims_to_third_amino_acid_methionine(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_HMM_DETECTED)])
        self.fx.write_detected_proteins("g1", {"p1": "KKMP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "AAAAAAATGCCC"})
        self.fx.write_detected_rows("g1", [
            self.detected_row("p1", "g1", "ctg1", 1, 12, 1, 4),
        ])

        protein = CuratedProtein("p1", "g1")
        locus = protein.genomic_locus_with_leader()

        self.assertEqual(protein.sequence_with_leader(), "MP")
        self.assertEqual(locus.start_1b, 7)
        self.assertEqual(locus.end_1b, 12)
        self.assertEqual(locus.sequence(), "ATGCCC")
        self.assertEqual(locus.cds_intervals_1b, [(1, 6)])

    def test_hmm_detected_leader_trims_to_second_amino_acid_methionine_on_reverse_strand(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_HMM_DETECTED)])
        self.fx.write_detected_proteins("g1", {"p1": "KMP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "GGGCATTTT"})
        self.fx.write_detected_rows("g1", [
            self.detected_row("p1", "g1", "ctg1", 9, 1, 1, 3),
        ])

        protein = CuratedProtein("p1", "g1")
        locus = protein.genomic_locus_with_leader()

        self.assertEqual(protein.sequence_with_leader(), "MP")
        self.assertEqual(locus.start_1b, 6)
        self.assertEqual(locus.end_1b, 1)
        self.assertEqual(locus.sequence(), "ATGCCC")
        self.assertEqual(locus.cds_intervals_1b, [(1, 6)])

    def test_hmm_detected_leader_does_not_extend_existing_methionine(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_HMM_DETECTED)])
        self.fx.write_detected_proteins("g1", {"p1": "MP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "AAACCCGGG"})
        self.fx.write_detected_rows("g1", [
            self.detected_row("p1", "g1", "ctg1", 4, 9, 1, 2),
        ])

        protein = CuratedProtein("p1", "g1")

        self.assertEqual(protein.sequence_with_leader(), "MP")
        self.assertEqual(protein.genomic_locus_with_leader().start_1b, 4)

    def test_hmm_detected_leader_extends_reverse_strand(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_HMM_DETECTED)])
        self.fx.write_detected_proteins("g1", {"p1": "PG"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "AAAAAACCCGGGTTTCAT"})
        self.fx.write_detected_rows("g1", [
            self.detected_row("p1", "g1", "ctg1", 12, 7, 1, 2),
        ])

        protein = CuratedProtein("p1", "g1")
        locus = protein.genomic_locus_with_leader()

        self.assertEqual(protein.sequence_with_leader(), "MKPG")
        self.assertEqual(locus.start_1b, 18)
        self.assertEqual(locus.end_1b, 7)
        self.assertEqual(locus.sequence(), "ATGAAACCCGGG")

    def test_ncbi_forward_strand_sequence_and_locus_from_gff(self):
        self.fx.write_manifest([self.manifest_row("pncbi", "g1", SEQUENCE_SOURCE_NCBI)])
        self.fx.write_ncbi_proteins("g1", {"pncbi": "MGP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "AAACCCGGGTTTAAACCCGGGTTTAAA"})
        self.fx.write_gff("g1", "\n".join([
            "ctg1\tsrc\tgene\t4\t24\t.\t+\t.\tID=gene1",
            "ctg1\tsrc\tmRNA\t4\t24\t.\t+\t.\tID=tx1;Parent=gene1",
            "ctg1\tsrc\tCDS\t10\t15\t.\t+\t0\tID=cds1;Parent=tx1;protein_id=pncbi",
            "ctg1\tsrc\tCDS\t19\t21\t.\t+\t0\tID=cds2;Parent=tx1;protein_id=pncbi",
            "ctg1\tsrc\tstart_codon\t10\t12\t.\t+\t0\tParent=tx1",
            "ctg1\tsrc\tstop_codon\t22\t24\t.\t+\t0\tParent=tx1",
            "ctg1\tsrc\tgene\t25\t27\t.\t+\t.\tID=gene2",
            "",
        ]))

        protein = CuratedProtein("pncbi", "g1")
        locus = protein.genomic_locus()

        self.assertEqual(protein.sequence(), "MGP")
        self.assertEqual(locus.sequence(), "CCCGGGTTTAAACCCGGGTTT")
        self.assertEqual(locus.cds_intervals_1b, [(7, 12), (16, 18)])
        self.assertEqual(locus.start_codon_position_1b(), 7)
        self.assertEqual(locus.stop_codon_position_1b(), 19)
        self.assertEqual(locus.dss_positions_1b(), [12])
        self.assertEqual(locus.ass_positions_1b(), [16])

    def test_ncbi_locus_falls_back_to_two_pass_gff_scan_if_rows_not_grouped_by_gene(self):
        self.fx.write_manifest([self.manifest_row("pncbi", "g1", SEQUENCE_SOURCE_NCBI)])
        self.fx.write_ncbi_proteins("g1", {"pncbi": "MGP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "AAACCCGGGTTTAAACCCGGGTTTAAA"})
        self.fx.write_gff("g1", "\n".join([
            "ctg1\tsrc\tgene\t1\t3\t.\t+\t.\tID=gene_other",
            "ctg1\tsrc\tmRNA\t4\t24\t.\t+\t.\tID=tx1",
            "ctg1\tsrc\tCDS\t10\t15\t.\t+\t0\tID=cds1;Parent=tx1;protein_id=pncbi",
            "ctg1\tsrc\tgene\t25\t27\t.\t+\t.\tID=gene_late",
            "ctg1\tsrc\tCDS\t19\t21\t.\t+\t0\tID=cds2;Parent=tx1;protein_id=pncbi",
            "ctg1\tsrc\tstart_codon\t10\t12\t.\t+\t0\tParent=tx1",
            "ctg1\tsrc\tstop_codon\t22\t24\t.\t+\t0\tParent=tx1",
            "",
        ]))

        protein = CuratedProtein("pncbi", "g1")
        locus = protein.genomic_locus()

        self.assertEqual(locus.sequence(), "CCCGGGTTTAAACCCGGGTTT")
        self.assertEqual(locus.cds_intervals_1b, [(7, 12), (16, 18)])
        self.assertEqual(locus.start_codon_position_1b(), 7)
        self.assertEqual(locus.stop_codon_position_1b(), 19)

    def test_ncbi_reverse_strand_sequence_and_locus_from_gff(self):
        self.fx.write_manifest([self.manifest_row("pncbi", "g1", SEQUENCE_SOURCE_NCBI)])
        self.fx.write_ncbi_proteins("g1", {"pncbi": "MGP"})
        self.fx.write_genomic_fasta("g1", {"ctg1": "AAACCCGGGTTTAAACCCGGGTTTAAA"})
        self.fx.write_gff("g1", "\n".join([
            "ctg1\tsrc\tmRNA\t4\t24\t.\t-\t.\tID=tx1",
            "ctg1\tsrc\tCDS\t19\t21\t.\t-\t0\tID=cds1;Parent=tx1;protein_id=pncbi",
            "ctg1\tsrc\tCDS\t10\t15\t.\t-\t0\tID=cds2;Parent=tx1;protein_id=pncbi",
            "ctg1\tsrc\tstart_codon\t19\t21\t.\t-\t0\tParent=tx1",
            "ctg1\tsrc\tstop_codon\t4\t6\t.\t-\t0\tParent=tx1",
            "",
        ]))

        protein = CuratedProtein("pncbi", "g1")
        locus = protein.genomic_locus()

        self.assertEqual(protein.sequence(), "MGP")
        self.assertEqual(locus.start_1b, 24)
        self.assertEqual(locus.end_1b, 4)
        self.assertEqual(locus.strand, -1)
        self.assertEqual(locus.sequence(), "AAACCCGGGTTTAAACCCGGG")
        self.assertEqual(locus.cds_intervals_1b, [(4, 6), (10, 15)])
        self.assertEqual(locus.start_codon_position_1b(), 4)
        self.assertEqual(locus.stop_codon_position_1b(), 19)
        self.assertEqual(locus.dss_positions_1b(), [6])
        self.assertEqual(locus.ass_positions_1b(), [10])

    def test_detected_pfam_and_ko_return_rows_for_protein_query(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_HMM_DETECTED)])

        pfam_rows = [
            dict(detection_type="sequence", detection_method="hmm", batch="b1",
                 query_accession="p1", query_database="g1", query_type="protein",
                 target_accession="PF1", target_database="Pfam", target_type="protein",
                 query_start=1, query_end=10, target_start=1, target_end=10),
            dict(detection_type="sequence", detection_method="hmm", batch="b1",
                 query_accession="p2", query_database="g1", query_type="protein",
                 target_accession="PF2", target_database="Pfam", target_type="protein",
                 query_start=1, query_end=10, target_start=1, target_end=10),
        ]
        ko_rows = [
            dict(detection_type="sequence", detection_method="hmm", batch="b1",
                 query_accession="p1", query_database="g1", query_type="protein",
                 target_accession="K00001", target_database="KO", target_type="protein",
                 query_start=1, query_end=10, target_start=1, target_end=10),
        ]
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), pfam_rows)
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_ko_assigned.tsv"), ko_rows)

        protein = CuratedProtein("p1", "g1")

        self.assertEqual([row["target_accession"] for row in protein.detected_pfam()], ["PF1"])
        self.assertEqual([row["target_accession"] for row in protein.detected_ko()], ["K00001"])

    def test_reuses_loaded_duckdb_tables_across_curated_proteins(self):
        self.fx.write_manifest([
            self.manifest_row("p1", "g1", SEQUENCE_SOURCE_HMM_DETECTED),
            self.manifest_row("p2", "g1", SEQUENCE_SOURCE_NCBI),
        ])
        pfam_rows = [
            dict(detection_type="sequence", detection_method="hmm", batch="b1",
                 query_accession="p1", query_database="g1", query_type="protein",
                 target_accession="PF1", target_database="Pfam", target_type="protein",
                 query_start=1, query_end=10, target_start=1, target_end=10),
            dict(detection_type="sequence", detection_method="hmm", batch="b1",
                 query_accession="p2", query_database="g1", query_type="protein",
                 target_accession="PF2", target_database="Pfam", target_type="protein",
                 query_start=1, query_end=10, target_start=1, target_end=10),
        ]
        DetectedTable.write_tsv(str(self.fx.area_genomics / "protein_pfam.tsv"), pfam_rows)

        original_duckdb_load = protein_module.Schema.duckdb_load
        loaded_schemas = []

        def record_load(schema):
            loaded_schemas.append(schema.name)
            return original_duckdb_load(schema)

        with patch.object(protein_module.Schema, "duckdb_load", autospec=True, side_effect=record_load):
            self.assertEqual(CuratedProtein("p1", "g1").sequence_source, SEQUENCE_SOURCE_HMM_DETECTED)
            self.assertEqual(CuratedProtein("p2", "g1").sequence_source, SEQUENCE_SOURCE_NCBI)
            self.assertEqual([row["target_accession"] for row in CuratedProtein("p1", "g1").detected_pfam()], ["PF1"])
            self.assertEqual([row["target_accession"] for row in CuratedProtein("p2", "g1").detected_pfam()], ["PF2"])

        self.assertEqual(len(loaded_schemas), 2)

    def test_hmm_align_runs_hmmalign_and_returns_alignment(self):
        self.fx.write_manifest([self.manifest_row("p1", "g1", SEQUENCE_SOURCE_HMM_DETECTED)])
        self.fx.write_detected_proteins("g1", {"p1": "ACDE"})
        stockholm = "\n".join([
            "# STOCKHOLM 1.0",
            "p1 AC-DE",
            "#=GC RF xxxxx",
            "//",
            "",
        ])

        with patch("sieve.protein.subprocess.run") as mock_run:
            mock_run.return_value = CompletedProcess(
                args=["hmmalign", "profile.hmm", "protein.faa"],
                returncode=0,
                stdout=stockholm,
                stderr="",
            )

            alignment = CuratedProtein("p1", "g1").hmm_align("profile.hmm")

        cmd = mock_run.call_args.args[0]
        self.assertEqual(cmd[0:2], ["hmmalign", "profile.hmm"])
        self.assertTrue(cmd[2].endswith(".faa"))
        self.assertEqual(alignment.aa_hmm_pos_1b(1), (1, "A"))
        self.assertEqual(alignment.aa_hmm_pos_1b(2), (2, "C"))
        self.assertIsNone(alignment.aa_hmm_pos_1b(3))
        self.assertEqual(alignment.aa_hmm_pos_1b(4), (3, "D"))


class TestProteinHMMAlignment(unittest.TestCase):

    def test_aa_hmm_pos_1b_returns_sequence_position_and_letter(self):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("A-CDE"), id="p1"),
        ])
        alignment.column_annotations["reference_annotation"] = "xxx.x"

        protein_alignment = ProteinHMMAlignment(alignment, "p1")

        self.assertEqual(protein_alignment.aa_hmm_pos_1b(1), (1, "A"))
        self.assertIsNone(protein_alignment.aa_hmm_pos_1b(2))
        self.assertEqual(protein_alignment.aa_hmm_pos_1b(3), (2, "C"))
        self.assertEqual(protein_alignment.aa_hmm_pos_1b(4), (4, "E"))
