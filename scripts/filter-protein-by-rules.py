#!/usr/bin/env python3

import argparse
import os
import sys
import tempfile

from sieve.protein import CuratedProtein
from sieve.rule_loader import load_rules
from sieve.rules import RULE_MAYBE, RULE_TRUE
from tangle.sequence import write_fasta_from_dict


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
        key = (parts[0], parts[1])
        if key not in seen:
            keys.append(key)
            seen.add(key)
    return keys


def write_rule_fasta(rows, fasta_output, include_maybe=True):
    accepted = {RULE_TRUE}
    if include_maybe:
        accepted.add(RULE_MAYBE)

    fasta = {}
    for row in rows:
        if row["pass all"] not in accepted:
            continue
        protein_accession = row["protein accession"]
        genome_accession = row["genome accession"]
        fasta[protein_accession] = CuratedProtein(protein_accession, genome_accession).sequence()
    write_fasta_from_dict(fasta, fasta_output)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--rule", required=True)
    parser.add_argument("--artifacts-dir")
    parser.add_argument("--fasta-output")
    parser.add_argument("--fasta-excludes-maybe", action="store_true", default=False)
    args = parser.parse_args(argv)

    rules = load_rules(args.rule)
    protein_keys = read_protein_keys(sys.stdin)

    with tempfile.TemporaryDirectory() as tmpd:
        if args.artifacts_dir is not None:
            os.makedirs(args.artifacts_dir, exist_ok=True)
            output_tsv = os.path.join(args.artifacts_dir, "rule-results.tsv")
        else:
            output_tsv = os.path.join(tmpd, "rule-results.tsv")

        rows = rules.check(
            protein_keys,
            output_tsv,
            artifacts_dir=args.artifacts_dir,
        )

        if args.fasta_output is not None:
            write_rule_fasta(
                rows,
                args.fasta_output,
                include_maybe=not args.fasta_excludes_maybe,
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
