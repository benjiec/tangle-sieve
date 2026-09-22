#!/usr/bin/env python3
"""Rank locus CDS sequences by normalized local nucleotide alignment score."""

import argparse
import csv
import math
import os
import sys

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.Seq import Seq

from tangle import open_file_to_read, open_file_to_write


ALPHABET = "ACGTRYSWKMBDHVN"
FIELDS = [
    "transcript_id", "locus_id", "rank", "score", "raw_score", "strand",
    "transcript_length", "locus_length", "transcript_coverage",
    "identity", "gap_count", "gap_bases", "transcript_start",
    "transcript_end", "locus_start", "locus_end",
]


def make_aligner(match=2, mismatch=-3, gap_open=-5, gap_extend=-1):
    values = (match, mismatch, gap_open, gap_extend)
    if not all(math.isfinite(v) for v in values):
        raise ValueError("scoring parameters must be finite")
    if match <= 0 or any(v > 0 for v in values[1:]):
        raise ValueError("match must be positive; mismatch and gap scores must be nonpositive")
    matrix = substitution_matrices.Array(ALPHABET, dims=2)
    for a in ALPHABET:
        for b in ALPHABET:
            matrix[a, b] = match if a == b and a in "ACGT" else mismatch
    aligner = Align.PairwiseAligner(mode="local")
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    return aligner


def read_fasta(path):
    records = []
    seen = set()
    with open_file_to_read(path) as stream:
        for record in SeqIO.parse(stream, "fasta"):
            if record.id in seen:
                continue
            seen.add(record.id)
            sequence = str(record.seq).upper().replace("U", "T")
            if not sequence:
                raise ValueError(f"empty sequence in {path}: {record.id}")
            invalid = set(sequence) - set(ALPHABET)
            if invalid:
                raise ValueError(f"invalid nucleotide(s) in {record.id}: {''.join(sorted(invalid))}")
            records.append((record.id, sequence))
    if not records:
        raise ValueError(f"no FASTA records in {path}")
    return records


def align_pair(aligner, locus, transcript, match):
    candidates = []
    for strand, query in (("+", transcript), ("-", str(Seq(transcript).reverse_complement()))):
        alignments = aligner.align(locus, query)
        alignment = next(iter(alignments), None)
        row = dict(raw_score=0.0, score=0.0, strand=".",
                   transcript_coverage=0.0, identity=0.0, gap_count=0,
                   gap_bases=0, transcript_start="", transcript_end="",
                   locus_start="", locus_end="")
        if alignment is not None:
            coordinates = alignment.coordinates
            paired = matches = gap_count = gap_bases = 0
            for (a, b), (c, d) in zip(coordinates.T[:-1], coordinates.T[1:]):
                da, db = int(c - a), int(d - b)
                if da and db:
                    paired += db
                    matches += sum(x == y and x in "ACGT"
                                   for x, y in zip(locus[a:c], query[b:d]))
                else:
                    gap_count += 1
                    gap_bases += da + db
            start, end = map(int, coordinates[1, [0, -1]])
            if strand == "-":
                start, end = len(transcript) - end, len(transcript) - start
            row.update(raw_score=float(alignment.score),
                       score=100 * alignment.score / (match * len(transcript)),
                       strand=strand, transcript_coverage=paired / len(transcript),
                       identity=matches / paired if paired else 0.0,
                       gap_count=gap_count, gap_bases=gap_bases,
                       transcript_start=start, transcript_end=end,
                       locus_start=int(coordinates[0, 0]),
                       locus_end=int(coordinates[0, -1]))
        candidates.append(row)
    return max(candidates, key=ranking_key)


def ranking_key(row):
    return row["raw_score"], row["transcript_coverage"], -row["gap_bases"]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("loci_fasta", help="exon-only locus CDS nucleotide FASTA")
    parser.add_argument("transcripts_fasta", help="transcript nucleotide FASTA")
    parser.add_argument("output_tsv")
    parser.add_argument("--match", type=float, default=2)
    parser.add_argument("--mismatch", type=float, default=-3)
    parser.add_argument("--gap-open", type=float, default=-5)
    parser.add_argument("--gap-extend", type=float, default=-1,
                        help="gap length k scores gap-open + (k-1)*gap-extend")
    args = parser.parse_args(argv)
    try:
        aligner = make_aligner(args.match, args.mismatch, args.gap_open, args.gap_extend)
        loci = read_fasta(args.loci_fasta)
        transcripts = read_fasta(args.transcripts_fasta)
    except ValueError as error:
        parser.error(str(error))
    write_header = True
    if os.path.exists(args.output_tsv):
        with open_file_to_read(args.output_tsv) as existing:
            write_header = not existing.read(1)
    with open_file_to_write(args.output_tsv, "at") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        if write_header:
            writer.writeheader()
        for transcript_id, transcript in transcripts:
            rows = []
            for locus_id, locus in loci:
                row = align_pair(aligner, locus, transcript, args.match)
                row.update(transcript_id=transcript_id, locus_id=locus_id,
                           transcript_length=len(transcript), locus_length=len(locus))
                rows.append(row)
            rows.sort(key=ranking_key, reverse=True)
            previous = None
            rank = 0
            for position, row in enumerate(rows, 1):
                key = ranking_key(row)
                if key != previous:
                    rank = position
                previous = key
                row["rank"] = rank
                writer.writerow(row)
                if rank == 1:
                    print(f"{transcript_id}: {row['locus_id']}, {row['score']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
