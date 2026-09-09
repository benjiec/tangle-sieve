#!/usr/bin/env python3

import argparse
import sys

import duckdb

from tangle import open_file_to_read, open_file_to_write, unique_batch
from tangle.detected import DetectedTable
from tangle.models import CSVSource, Schema


def positive_int(value):
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("match count must be at least 1")
    return parsed


def _sql_string(value):
    return "'" + str(value).replace("'", "''") + "'"


def read_unique_fasta(path):
    sequences = {}
    accession = None
    parts = []

    def store():
        if accession is None:
            return
        if accession in sequences:
            raise ValueError(f"Duplicate FASTA accession: {accession}")
        sequences[accession] = "".join(parts)

    with open_file_to_read(path) as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                store()
                header = line[1:].strip()
                if not header:
                    raise ValueError("FASTA header has no accession")
                accession = header.split(None, 1)[0]
                parts = []
            else:
                if accession is None:
                    raise ValueError("FASTA sequence occurs before its first header")
                parts.append(line)
    store()
    return sequences


def find_match_regions(tsv, pfam_accession):
    accession = _sql_string(pfam_accession)
    version_prefix = _sql_string(f"{pfam_accession}.")
    filters = [
        f"(target_accession = {accession} OR starts_with(target_accession, {version_prefix}))"
    ]
    schema = Schema("__pfam_extract_region__" + unique_batch())
    schema.add_table(CSVSource(DetectedTable, tsv, load_filters=filters))
    schema.duckdb_load()
    try:
        rows = duckdb.execute(f"""
            SELECT DISTINCT query_accession, query_start, query_end
              FROM {schema.name}.{DetectedTable.name}
             ORDER BY query_accession, least(query_start, query_end),
                      greatest(query_start, query_end)
        """).fetchall()
    finally:
        schema.duckdb_drop()

    regions = {}
    for query_accession, query_start, query_end in rows:
        left = min(query_start, query_end)
        right = max(query_start, query_end)
        regions.setdefault(query_accession, []).append((left, right))
    return regions


def extract_regions(sequences, matches, pfam_accession, minimum, maximum, report=None):
    if minimum < 1 or maximum < 1:
        raise ValueError("match counts must be at least 1")
    if minimum > maximum:
        raise ValueError("minimum matches cannot exceed maximum matches")
    report = report or (lambda message: print(message, file=sys.stderr))
    extracted = {}

    for accession, sequence in sequences.items():
        sequence_matches = matches.get(accession, [])
        count = len(sequence_matches)
        if not minimum <= count <= maximum:
            report(
                f"Ignoring {accession}: expected {minimum}-{maximum} {pfam_accession} "
                f"matches, found {count}"
            )
            continue

        start = min(region[0] for region in sequence_matches)
        end = max(region[1] for region in sequence_matches)
        if start < 1 or end > len(sequence):
            report(
                f"Ignoring {accession}: match boundaries {start}-{end} are outside "
                f"sequence length {len(sequence)}"
            )
            continue

        output_accession = f"{accession}_{pfam_accession}_{minimum}_{maximum}"
        extracted[output_accession] = sequence[start - 1:end]

    return extracted


def write_fasta(sequences, output):
    with open_file_to_write(output, "wt") as stream:
        for accession, sequence in sequences.items():
            stream.write(f">{accession}\n{sequence}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract FASTA regions spanning a required number of Pfam matches."
    )
    parser.add_argument("fasta", help="input protein FASTA")
    parser.add_argument("tsv", help="detected-table TSV containing Pfam matches")
    parser.add_argument("pfam_accession", help="Pfam accession, for example PF00023")
    parser.add_argument("minimum", type=positive_int, help="minimum match count (inclusive)")
    parser.add_argument("maximum", type=positive_int, help="maximum match count (inclusive)")
    parser.add_argument("output", help="output FASTA")
    args = parser.parse_args(argv)

    if args.minimum > args.maximum:
        parser.error("minimum matches cannot exceed maximum matches")

    sequences = read_unique_fasta(args.fasta)
    matches = find_match_regions(args.tsv, args.pfam_accession)
    extracted = extract_regions(
        sequences,
        matches,
        args.pfam_accession,
        args.minimum,
        args.maximum,
    )
    write_fasta(extracted, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
