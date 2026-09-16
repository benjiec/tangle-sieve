#!/usr/bin/env python3
"""Print GFF CDS coordinates and their portions of the annotated protein.

Coordinates are 1-based, inclusive. Parentheses identify residues whose codons
span CDS rows. Residues come from the manifest-selected protein FASTA, preserving
the annotated translation (including translation exceptions).
"""

import argparse
import sys
from urllib.parse import unquote

from tangle.defaults import Defaults
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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("protein_accession")
    args = parser.parse_args(argv)
    lines = []
    try:
        for genome in find_genomes(args.protein_accession):
            protein = CuratedProtein(args.protein_accession, genome)
            rows = read_cds(protein)
            portions = exon_portions(rows, protein.sequence())
            lines.append(f"genome {genome}")
            for number, (row, portion) in enumerate(portions, 1):
                lines.append(
                    f"exon {number}, {row['seqid']}, {row['start']}, {row['end']}: {portion}"
                )
    finally:
        CuratedProtein.clear_cache()
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
