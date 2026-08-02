#!/usr/bin/env python3

import argparse
import sys

import duckdb

from tangle.defaults import Defaults
from tangle.detected import DetectedTable
from tangle.models import CSVSource, Schema
from tangle import unique_batch

from sieve.taxonomy import read_taxonomy_rows, taxonomy_matches


def _sql_string(value):
    return "'" + str(value).replace("'", "''") + "'"


def find_matches(ko_accession, max_evalue=None, taxon=None):
    filters = [f"target_accession = {_sql_string(ko_accession)}"]
    if max_evalue is not None:
        filters.append(f"evalue <= {float(max_evalue)}")

    schema = Schema("__ko_find_matches__" + unique_batch())
    source = CSVSource(
        DetectedTable,
        Defaults.area_protein_ko_assigned_tsv(),
        load_filters=filters,
    )
    schema.add_table(source)
    schema.duckdb_load()
    try:
        query = f"""
            SELECT DISTINCT query_accession, query_database
              FROM {schema.name}.{DetectedTable.name}
             ORDER BY query_database, query_accession
        """
        matches = [
            (row["query_accession"], row["query_database"])
            for row in duckdb.execute(query).fetchdf().to_dict("records")
        ]
        if taxon is None:
            return matches
        taxonomy_by_genome = read_taxonomy_rows()
        return [
            match for match in matches
            if taxonomy_matches(taxonomy_by_genome.get(match[1], {}), taxon)
        ]
    finally:
        schema.duckdb_drop()


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("ko_accession")
    parser.add_argument("--max-evalue", type=float)
    parser.add_argument("--taxon")
    args = parser.parse_args(argv)

    for protein_accession, genome_accession in find_matches(
        args.ko_accession,
        args.max_evalue,
        taxon=args.taxon,
    ):
        print(f"{protein_accession}\t{genome_accession}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
