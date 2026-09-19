#!/usr/bin/env python3
"""Print GFF CDS coordinates and their portions of the annotated protein.

Coordinates are 1-based, inclusive. Parentheses identify residues whose codons
span CDS rows. Residues come from the manifest-selected protein FASTA, preserving
the annotated translation (including translation exceptions).
"""

import argparse
import sys
from urllib.parse import unquote

from Bio.Seq import Seq
from tangle.defaults import Defaults
from tangle.detected import DetectedTable
from tangle.manifest import ManifestTable
from sieve.protein import (
    CuratedProtein, _existing_path, _parse_gff_line, _rows_from_table, _sql_string,
)


def find_genomes(accession):
    path = _existing_path(Defaults.area_sequence_manifest_tsv())
    rows = _rows_from_table(ManifestTable, path, column_filters=[
        f"sequence_accession = {_sql_string(accession)}",
        "sequence_type = 'protein'",
    ])
    if not rows:
        raise ValueError(f"Cannot find protein {accession} in manifest")
    return sorted({row["sequence_database"] for row in rows})


def read_cds(protein):
    path = _existing_path(Defaults.ncbi_genome_gff(protein.genome_accession))
    rows = []
    with open(path, encoding="utf-8") as source:
        for line in source:
            if line.startswith("##FASTA"):
                break
            row = _parse_gff_line(line)
            if row is None or row["type"] != "CDS":
                continue
            row["attrs"] = {key: unquote(value) for key, value in row["attrs"].items()}
            if not protein._gff_row_matches_protein(row):
                continue
            phase = line.rstrip().split("\t")[7]
            if phase not in ("0", "1", "2"):
                raise ValueError(f"Invalid CDS phase for {protein.protein_accession}: {phase}")
            row["phase"] = int(phase)
            rows.append(row)
    if not rows:
        raise ValueError(f"Cannot find CDS rows for protein {protein.protein_accession}")
    return rows


def read_fragment_groups(accession, path):
    path = _existing_path(path)
    rows = _rows_from_table(DetectedTable, path, column_filters=[
        f"target_accession = {_sql_string(accession)}",
        "target_type = 'protein'",
    ])
    if not rows:
        raise ValueError(f"Cannot find protein {accession} in fragments TSV")
    groups = {}
    for row in rows:
        genome = row.get("target_database")
        if not genome:
            raise ValueError(f"Fragment for {accession} has no target_database")
        if row.get("query_database") != genome:
            raise ValueError(f"Fragment databases disagree for {accession}")
        groups.setdefault(genome, []).append(row)
    return [(genome, groups[genome]) for genome in sorted(groups)]


def fragment_portions(rows, genomic_sequences, accession):
    for row in rows:
        if row.get("target_start") is None or row.get("target_end") is None:
            raise ValueError(f"Fragment for {accession} has no target coordinates")
    rows = sorted(rows, key=lambda row: (row["target_start"], row["target_end"]))
    contigs = {row["query_accession"] for row in rows}
    strands = {1 if row["query_start"] <= row["query_end"] else -1 for row in rows}
    if len(contigs) != 1:
        raise ValueError(f"Fragments for {accession} span multiple contigs")
    if len(strands) != 1:
        raise ValueError(f"Fragments for {accession} span multiple strands")
    result = []
    previous_end = None
    for row in rows:
        target_start, target_end = row["target_start"], row["target_end"]
        if target_start < 1 or target_end < target_start or (previous_end is not None and target_start <= previous_end):
            raise ValueError(f"Fragments for {accession} have overlapping or invalid model coordinates")
        cds = {
            "seqid": row["query_accession"],
            "start": min(row["query_start"], row["query_end"]),
            "end": max(row["query_start"], row["query_end"]),
            "strand": "+" if row["query_start"] <= row["query_end"] else "-",
        }
        dna = exon_dna(cds, genomic_sequences)
        if len(dna) % 3:
            raise ValueError(f"Genomic fragment length is not divisible by three for {accession}")
        result.append((cds, str(Seq(dna).translate(table="Standard", to_stop=False))))
        previous_end = target_end
    return result


def collate_fragment_protein(fragment_rows, portions):
    rows = sorted(fragment_rows, key=lambda row: (row["target_start"], row["target_end"]))
    sequence = portions[0][1]
    for previous, row, (_cds, portion) in zip(rows, rows[1:], portions[1:]):
        sequence += "X" * max(0, row["target_start"] - previous["target_end"] - 1)
        sequence += portion
    return sequence


def exon_portions(rows, sequence):
    """Map CDS nucleotide intervals onto annotated amino acid positions."""
    if not rows:
        raise ValueError("No CDS rows")
    locations = {(row["seqid"], row["strand"]) for row in rows}
    if len(locations) != 1 or rows[0]["strand"] not in ("+", "-"):
        raise ValueError("CDS rows must share one contig and a defined strand")
    groups = {(row["attrs"].get("ID"), row["attrs"].get("Parent")) for row in rows}
    if len(groups) != 1:
        raise ValueError("Protein matches multiple CDS annotations; cannot assign exon order")
    rows = sorted(rows, key=lambda row: row["start"], reverse=rows[0]["strand"] == "-")
    genomic_order = sorted(rows, key=lambda row: row["start"])
    for index, row in enumerate(genomic_order):
        if row["start"] < 1 or row["end"] < row["start"]:
            raise ValueError("Invalid CDS coordinates")
        if index and row["start"] <= genomic_order[index - 1]["end"]:
            raise ValueError("Overlapping CDS rows cannot be mapped without resolving a frameshift")

    # Only the initial phase skips bases. Internal phase bases complete the
    # preceding codon and must remain in the concatenated coding sequence.
    offset = -rows[0]["phase"]
    intervals = []
    for index, row in enumerate(rows):
        if index and row["phase"] != (-offset) % 3:
            raise ValueError("CDS phase is inconsistent with preceding CDS lengths")
        end = offset + row["end"] - row["start"] + 1
        intervals.append((offset, end))
        offset = end
    sequence = sequence.rstrip("*")
    # NCBI CDS may include the terminal stop codon or a partial trailing codon.
    if not 0 <= offset - 3 * len(sequence) <= 3:
        raise ValueError("CDS length does not match annotated protein length")
    result = []
    for row, (start, end) in zip(rows, intervals):
        portion = []
        for aa in range(max(0, start // 3), min(len(sequence), (end + 2) // 3)):
            residue = sequence[aa]
            if start > aa * 3 or end < aa * 3 + 3:
                residue = f"({residue})"
            portion.append(residue)
        result.append((row, "".join(portion)))
    return result


def exon_dna(row, genomic_sequences):
    contig = genomic_sequences.get(row["seqid"])
    if contig is None:
        raise ValueError(f"Cannot find contig sequence {row['seqid']}")
    if row["end"] > len(contig):
        raise ValueError(f"CDS coordinates exceed contig sequence {row['seqid']}")
    dna = contig[row["start"] - 1:row["end"]].upper()
    if row["strand"] == "-":
        dna = str(Seq(dna).reverse_complement())
    return dna


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("protein_accessions", nargs="+")
    parser.add_argument("--dna", action="store_true", help="append the full CDS exon DNA in translation direction")
    parser.add_argument("--gc", action="store_true", help="append GC percentage (G+C divided by all exon bases, including ambiguous bases)")
    parser.add_argument("--fasta", action="store_true", help="print FASTA instead of enumerating CDS exons")
    parser.add_argument("--fragments-tsv", help="derive CDS fragments from this detected-fragments TSV instead of GFF")
    args = parser.parse_args(argv)
    lines = []
    try:
        for protein_accession in args.protein_accessions:
            if args.fragments_tsv:
                genome_rows = read_fragment_groups(protein_accession, args.fragments_tsv)
            else:
                genome_rows = [(genome, None) for genome in find_genomes(protein_accession)]
            for genome, fragment_rows in genome_rows:
                protein = CuratedProtein(protein_accession, genome)
                if fragment_rows is None:
                    protein.manifest_entry
                    rows = read_cds(protein)
                    protein_sequence = protein.sequence()
                    portions = exon_portions(rows, protein_sequence)
                else:
                    genomic_sequences = protein._genomic_sequences()
                    portions = fragment_portions(fragment_rows, genomic_sequences, protein_accession)
                    protein_sequence = collate_fragment_protein(fragment_rows, portions)
                if fragment_rows is None:
                    genomic_sequences = protein._genomic_sequences() if args.dna or (args.gc and not args.fasta) else None
                if args.fasta:
                    if args.dna:
                        sequence = "".join(exon_dna(row, genomic_sequences) for row, _ in portions)
                        lines.extend((f">{protein_accession}_cds", sequence))
                    else:
                        lines.extend((f">{protein_accession}", protein_sequence))
                    continue
                lines.append(f"genome {genome}")
                for number, (row, portion) in enumerate(portions, 1):
                    if genomic_sequences is not None:
                        dna = exon_dna(row, genomic_sequences)
                        if args.dna:
                            portion += f", {dna}"
                        if args.gc:
                            gc = 100 * (dna.count("G") + dna.count("C")) / len(dna)
                            portion += f", {gc:.2f}%"
                    lines.append(
                        f"exon {number}, {row['seqid']}, {row['start']}, {row['end']}: {portion}"
                    )
    finally:
        CuratedProtein.clear_cache()
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
