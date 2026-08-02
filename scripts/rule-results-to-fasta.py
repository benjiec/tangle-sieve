#!/usr/bin/env python3

import argparse
import csv
import sys

from tangle.sequence import read_fasta_as_dict, write_fasta_from_dict

from sieve.artifacts import rule_results_tsv, sequences_fasta
from sieve.taxonomy import normalized, read_taxonomy_rows, taxonomy_matches


def parse_rule_filter(value):
    if "=" not in value:
        raise ValueError(f"Expected rule filter in column=value form, got: {value}")
    column, expected = value.split("=", 1)
    column = column.strip()
    expected = expected.strip()
    if not column or not expected:
        raise ValueError(f"Expected rule filter in column=value form, got: {value}")
    return column, expected


def rule_filters(values):
    filters = [parse_rule_filter(value) for value in values]
    if not any(normalized(column) == "pass all" for column, _expected in filters):
        filters.insert(0, ("pass all", "true"))
    return filters


def read_rule_rows(artifacts_dir):
    with open(rule_results_tsv(artifacts_dir), "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def row_matches(row, filters, taxonomy_by_genome, taxon):
    for column, expected in filters:
        if column not in row:
            raise ValueError(f"Rule results are missing column: {column}")
        if row[column] != expected:
            return False
    if taxon is not None:
        taxonomy_row = taxonomy_by_genome.get(row.get("genome accession", ""))
        if taxonomy_row is None or not taxonomy_matches(taxonomy_row, taxon):
            return False
    return True


def select_sequences(artifacts_dir, filters, taxon=None):
    rows = read_rule_rows(artifacts_dir)
    sequences = read_fasta_as_dict(sequences_fasta(artifacts_dir))
    taxonomy_by_genome = read_taxonomy_rows() if taxon is not None else {}
    selected = {}
    for row in rows:
        if not row_matches(row, filters, taxonomy_by_genome, taxon):
            continue
        accession = row.get("sequence accession", "")
        if accession and accession not in selected:
            if accession not in sequences:
                raise ValueError(f"Cannot find sequence accession in sequences.faa: {accession}")
            selected[accession] = sequences[accession]
    return selected


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifacts-dir", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--rule", action="append", default=[])
    parser.add_argument("--taxon")
    args = parser.parse_args(argv)

    write_fasta_from_dict(
        select_sequences(args.artifacts_dir, rule_filters(args.rule), taxon=args.taxon),
        args.output,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
