#!/usr/bin/env python3

import argparse
import sys

import duckdb

from tangle import open_file_to_write, unique_batch
from tangle.defaults import Defaults
from tangle.detected import DetectedTable
from tangle.models import CSVSource, Schema

from sieve.protein import CuratedProtein
from sieve.taxonomy import read_taxonomy_rows, taxonomy_matches


def _sql_string(value):
    return "'" + str(value).replace("'", "''") + "'"


def find_matches(
    ko_accession,
    max_evalue=None,
    taxon=None,
    match_starts_before=None,
    match_ends_before=None,
    max_evalue_rank=None,
):
    filters = [f"target_accession = {_sql_string(ko_accession)}"]
    if max_evalue is not None:
        filters.append(f"evalue <= {float(max_evalue)}")
    if match_starts_before is not None:
        filters.append(f"query_start <= {int(match_starts_before)}")
    if match_ends_before is not None:
        filters.append(f"query_end <= {int(match_ends_before)}")
    if max_evalue_rank is not None:
        filters.extend([
            "custom_metric_name = 'evalue-rank'",
            f"TRY_CAST(custom_metric_value AS DOUBLE) <= {float(max_evalue_rank)}",
        ])

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


def _is_missing_manifest_error(error):
    return str(error).startswith("Cannot find protein ") and str(error).endswith(" in manifest")


def write_matches_fasta(matches, output):
    with open_file_to_write(output, "wt") as f:
        for protein_accession, genome_accession in matches:
            protein = CuratedProtein(protein_accession, genome_accession)
            try:
                try:
                    sequence = protein.sequence()
                except ValueError as e:
                    if not _is_missing_manifest_error(e):
                        raise
                    print(
                        f"Ignoring {protein_accession}\t{genome_accession}: {e}",
                        file=sys.stderr,
                    )
                    continue
                f.write(f">{protein_accession}\n{sequence}\n")
            finally:
                CuratedProtein.clear_cache()


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("ko_accession")
    parser.add_argument("--max-evalue", type=float)
    parser.add_argument("--match-starts-before", type=int)
    parser.add_argument("--match-ends-before", type=int)
    parser.add_argument("--max-evalue-rank", type=float, default=1)
    parser.add_argument("--taxon")
    parser.add_argument("-o", "--output")
    args = parser.parse_args(argv)

    matches = find_matches(
        args.ko_accession,
        args.max_evalue,
        taxon=args.taxon,
        match_starts_before=args.match_starts_before,
        match_ends_before=args.match_ends_before,
        max_evalue_rank=args.max_evalue_rank,
    )
    if args.output is not None:
        write_matches_fasta(matches, args.output)
        return 0
    for protein_accession, genome_accession in matches:
        print(f"{protein_accession}\t{genome_accession}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
