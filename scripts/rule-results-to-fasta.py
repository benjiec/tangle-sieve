#!/usr/bin/env python3

import argparse
import csv
import os
import sys

from tangle.defaults import Defaults
from tangle.sequence import read_fasta_as_dict, write_fasta_from_dict


TAXONOMY_FIELDS = {
    "domain",
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
}


def parse_rule_filter(value):
    if "=" not in value:
        raise ValueError(f"Expected rule filter in column=value form, got: {value}")
    column, expected = value.split("=", 1)
    column = column.strip()
    expected = expected.strip()
    if not column or not expected:
        raise ValueError(f"Expected rule filter in column=value form, got: {value}")
    return column, expected


def normalized(value):
    return str(value).strip().lower()


def rule_filters(values):
    filters = [parse_rule_filter(value) for value in values]
    if not any(normalized(column) == "pass all" for column, _expected in filters):
        filters.insert(0, ("pass all", "true"))
    return filters


def read_rule_rows(artifacts_dir):
    with open(os.path.join(artifacts_dir, "rule-results.tsv"), "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_taxonomy_rows():
    taxonomy_tsv = Defaults.area_genome_taxon_tsv()
    if not os.path.exists(taxonomy_tsv):
        return {}
    with open(taxonomy_tsv, "r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return {
        genome_accession(row): row
        for row in rows
        if genome_accession(row)
    }


def genome_accession(row):
    for key in row:
        if normalized(key).replace("_", " ") in {"genome accession", "accession"}:
            return row[key]
    return ""


def taxonomy_matches(row, taxon):
    if taxon is None:
        return True
    expected = normalized(taxon)
    for key, value in row.items():
        if normalized(key) in TAXONOMY_FIELDS and normalized(value) == expected:
            return True
    return False


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
    sequences = read_fasta_as_dict(os.path.join(artifacts_dir, "sequences.faa"))
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
