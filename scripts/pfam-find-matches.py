#!/usr/bin/env python3

import argparse
from collections import Counter, defaultdict
import re
import sys

import duckdb

from tangle import open_file_to_write, unique_batch
from tangle.defaults import Defaults
from tangle.detected import DetectedTable
from tangle.manifest import ManifestTable
from tangle.models import CSVSource, Schema

from sieve.protein import CuratedProtein
from sieve.taxonomy import read_taxonomy_rows, taxonomy_matches


def _sql_string(value):
    return "'" + str(value).replace("'", "''") + "'"


def parse_targets(accessions, exact_match=False):
    """Parse whole-accession repetition bounds; plain accessions mean one or more."""
    terms = []
    for expression in accessions:
        match = re.fullmatch(r"([^*+?{}\s]+)([?*+]|\{[0-9]+(?:,[0-9]+)?\})?", expression)
        if match is None:
            raise ValueError(f"Invalid Pfam expression: {expression}")
        accession, quantifier = match.groups()
        if not exact_match:
            accession = accession.split(".", 1)[0]
        minimum, maximum = 1, None
        if quantifier == "?":
            minimum, maximum = 0, 1
        elif quantifier == "*":
            minimum = 0
        elif quantifier and quantifier.startswith("{"):
            bounds = [int(n) for n in quantifier[1:-1].split(",")]
            minimum, maximum = bounds[0], bounds[-1]
            if minimum > maximum:
                raise ValueError(f"Invalid repetition bounds: {expression}")
        terms.append((accession, minimum, maximum))
    return terms


def architecture_matches(hits, terms, ordered=False, allow_overlap=False):
    """Consume all distinct requested-domain hits, ignoring unrequested domains."""
    hits = sorted(set(hits), key=lambda hit: (hit[1], hit[2], hit[0]))
    if not allow_overlap and any(a[2] >= b[1] for a, b in zip(hits, hits[1:])):
        return False
    if not ordered:
        counts = Counter(hit[0] for hit in hits)
        bounds = {}
        for accession, minimum, maximum in dict.fromkeys(terms):
            lo, hi = bounds.get(accession, (0, 0))
            bounds[accession] = (lo + minimum, None if hi is None or maximum is None else hi + maximum)
        return all(counts[a] >= lo and (hi is None or counts[a] <= hi)
                   for a, (lo, hi) in bounds.items())
    # Different domains at the same start have no resolvable sequence order.
    if any(a[1] == b[1] and a[0] != b[0] for a, b in zip(hits, hits[1:])):
        return False
    positions = {0}
    for accession, minimum, maximum in terms:
        following = set()
        for start in positions:
            end = start
            while end < len(hits) and hits[end][0] == accession:
                end += 1
            limit = end if maximum is None else min(end, start + maximum)
            following.update(range(start + minimum, limit + 1))
        positions = following
    return len(hits) in positions


def find_matches(pfam_accession, max_evalue=None, taxon=None, include_coordinates=False,
                 exact_match=False, sequence_source=None, ordered=False, allow_overlap=False):
    accessions = [pfam_accession] if isinstance(pfam_accession, str) else list(pfam_accession)
    if not accessions:
        raise ValueError("At least one Pfam accession is required")
    if include_coordinates and len(accessions) > 1:
        raise ValueError("Coordinates require a single Pfam accession")
    terms = parse_targets(accessions, exact_match)
    targets = {term[0] for term in terms}
    target_column = "target_accession" if exact_match else "split_part(target_accession, '.', 1)"
    filters = []
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
        rows = duckdb.execute(f"""
            SELECT query_accession, query_database, {target_column} AS accession,
                   query_start, query_end
              FROM {schema.name}.{DetectedTable.name}
        """).fetchdf().to_dict("records")
        by_query = defaultdict(set)
        for row in rows:
            key = (row["query_accession"], row["query_database"])
            hits = by_query[key]  # Retain queries with zero requested-domain hits.
            if row["accession"] in targets:
                hits.add((row["accession"], min(row["query_start"], row["query_end"]),
                          max(row["query_start"], row["query_end"])))
        matches = []
        for key, hits in sorted(by_query.items(), key=lambda item: (item[0][1], item[0][0])):
            if architecture_matches(hits, terms, ordered, allow_overlap):
                if include_coordinates:
                    matches.extend((*key, start, end) for _, start, end in sorted(hits, key=lambda hit: (hit[1], hit[2])))
                else:
                    matches.append(key)
        if sequence_source is not None:
            manifest_rows = CSVSource(
                ManifestTable,
                Defaults.area_sequence_manifest_tsv(),
                load_filters=[
                    "sequence_type = 'protein'",
                    f"sequence_source = {_sql_string(sequence_source)}",
                ],
            ).values()
            eligible = {
                (row["sequence_accession"], row["sequence_database"])
                for row in manifest_rows
            }
            matches = [match for match in matches if match[:2] in eligible]
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
    parser.add_argument("pfam_accession", nargs="+", help="Pfam accessions that must all match")
    parser.add_argument("--ordered", action="store_true", help="Require the requested domain architecture in argument order")
    parser.add_argument("--allow-overlap", action="store_true", help="Allow overlapping requested-domain hits")
    parser.add_argument("--exact-match", action="store_true", help="Match accession versions exactly")
    parser.add_argument("--max-evalue", type=float)
    parser.add_argument("--taxon")
    parser.add_argument("--sequence-source", help="Require an exact sequence_source match in the protein manifest")
    parser.add_argument("--match-only", action="store_true")
    parser.add_argument("-o", "--output")
    args = parser.parse_args(argv)

    if args.match_only and len(args.pfam_accession) > 1:
        parser.error("--match-only requires a single Pfam accession")
    if args.match_only and args.output is None:
        parser.error("--match-only requires --output")

    try:
        parse_targets(args.pfam_accession, args.exact_match)
    except ValueError as error:
        parser.error(str(error))

    find_kwargs = {
        "taxon": args.taxon,
        "exact_match": args.exact_match,
        "sequence_source": args.sequence_source,
        "ordered": args.ordered,
        "allow_overlap": args.allow_overlap,
    }
    if args.match_only:
        find_kwargs["include_coordinates"] = True
    matches = find_matches(args.pfam_accession, args.max_evalue, **find_kwargs)
    if args.output is not None:
        write_matches_fasta(
            matches,
            args.output,
            target_accession=args.pfam_accession[0],
            match_only=args.match_only,
        )
        return 0
    for protein_accession, genome_accession in matches:
        print(f"{protein_accession}\t{genome_accession}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
