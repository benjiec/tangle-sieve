#!/usr/bin/env python3

import argparse
import sys

from sieve.artifact_builder import build_artifacts
from sieve.protein import CuratedProtein
from sieve.rule_loader import load_rules


def read_protein_keys(lines):
    keys = []
    seen = set()
    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            raise ValueError(f"Expected two tab-separated columns, got: {raw_line.rstrip()}")
        key = tuple(parts)
        if key not in seen:
            keys.append(key)
            seen.add(key)
    return keys


def curated_proteins(keys):
    proteins = []
    for protein_accession, genome_accession in keys:
        protein = CuratedProtein(protein_accession, genome_accession)
        try:
            protein.manifest_entry
        except ValueError as error:
            message = str(error)
            if message.startswith("Cannot find protein ") and message.endswith(" in manifest"):
                print(f"Ignoring {protein_accession}\t{genome_accession}: {error}", file=sys.stderr)
                continue
            raise
        proteins.append(protein)
    return proteins


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--rule", required=True)
    parser.add_argument("--artifacts-dir", required=True)
    args = parser.parse_args(argv)

    proteins = curated_proteins(read_protein_keys(sys.stdin))
    try:
        build_artifacts(proteins, load_rules(args.rule), args.artifacts_dir, write_loci=True)
    finally:
        CuratedProtein.clear_cache()
    return 0


if __name__ == "__main__":
    sys.exit(main())
