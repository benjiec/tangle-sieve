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


def find_matches(pfam_accession, max_evalue=None, taxon=None, include_coordinates=False):
    accession = _sql_string(pfam_accession)
    version_prefix = _sql_string(f"{pfam_accession}.")
    filters = [
        f"(target_accession = {accession} OR starts_with(target_accession, {version_prefix}))"
    ]
    if max_evalue is not None:
        filters.append(f"evalue <= {float(max_evalue)}")

    schema = Schema("__pfam_find_matches__" + unique_batch())
    source = CSVSource(
        DetectedTable,
        Defaults.area_protein_pfam_tsv(),
        load_filters=filters,
    )
    schema.add_table(source)
    schema.duckdb_load()
    try:
        coordinate_columns = ", query_start, query_end" if include_coordinates else ""
        query = f"""
            SELECT DISTINCT query_accession, query_database{coordinate_columns}
              FROM {schema.name}.{DetectedTable.name}
             ORDER BY query_database, query_accession{coordinate_columns}
        """
        rows = duckdb.execute(query).fetchdf().to_dict("records")
        if include_coordinates:
            matches = [
                (
                    row["query_accession"],
                    row["query_database"],
                    row["query_start"],
                    row["query_end"],
                )
                for row in rows
            ]
        else:
            matches = [
                (row["query_accession"], row["query_database"])
                for row in rows
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


def write_matches_fasta(matches, output, target_accession=None, match_only=False):
    with open_file_to_write(output, "wt") as f:
        for match in matches:
            protein_accession, genome_accession = match[:2]
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
                output_accession = protein_accession
                if match_only:
                    query_start, query_end = match[2:]
                    left = min(query_start, query_end)
                    right = max(query_start, query_end)
                    sequence = sequence[left - 1:right]
                    output_accession = (
                        f"{protein_accession}_{target_accession}_{query_start}_{query_end}"
                    )
                f.write(f">{output_accession}\n{sequence}\n")
            finally:
                CuratedProtein.clear_cache()


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("pfam_accession")
    parser.add_argument("--max-evalue", type=float)
    parser.add_argument("--taxon")
    parser.add_argument("--match-only", action="store_true")
    parser.add_argument("-o", "--output")
    args = parser.parse_args(argv)

    if args.match_only and args.output is None:
        parser.error("--match-only requires --output")

    find_kwargs = {"taxon": args.taxon}
    if args.match_only:
        find_kwargs["include_coordinates"] = True
    matches = find_matches(args.pfam_accession, args.max_evalue, **find_kwargs)
    if args.output is not None:
        write_matches_fasta(
            matches,
            args.output,
            target_accession=args.pfam_accession,
            match_only=args.match_only,
        )
        return 0
    for protein_accession, genome_accession in matches:
        print(f"{protein_accession}\t{genome_accession}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
