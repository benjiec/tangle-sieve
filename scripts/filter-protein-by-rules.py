#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys
import tempfile
from collections import defaultdict

from tangle import open_file_to_write
from sieve.protein import CuratedProtein
from sieve.rule_loader import load_rules
from sieve.rules import RULE_MAYBE, RULE_TRUE


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


def _safe_filename(value):
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
    return value[:120] or "genome"


def _is_missing_manifest_error(error):
    return str(error).startswith("Cannot find protein ") and str(error).endswith(" in manifest")


def filter_manifest_protein_keys(protein_keys):
    filtered = []
    for protein_accession, genome_accession in protein_keys:
        protein = CuratedProtein(protein_accession, genome_accession)
        try:
            protein.manifest_entry
        except ValueError as e:
            if _is_missing_manifest_error(e):
                print(f"Ignoring {protein_accession}\t{genome_accession}: {e}", file=sys.stderr)
                continue
            raise
        filtered.append((protein_accession, genome_accession))
    return filtered


def group_protein_keys_by_genome(protein_keys):
    grouped = defaultdict(list)
    for key in protein_keys:
        grouped[key[1]].append(key)
    return [
        (genome_accession, grouped[genome_accession])
        for genome_accession in sorted(grouped)
    ]


def write_rule_rows(rows, output_tsv, rules):
    headers = ["protein accession", "genome accession", "pass all"] + [
        rule.label for rule in rules.atomic_rules()
    ]
    with open(output_tsv, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def check_rules_by_genome(rules, protein_keys, output_tsv, artifacts_dir=None):
    rows = []
    with tempfile.TemporaryDirectory() as tmpd:
        for genome_accession, genome_keys in group_protein_keys_by_genome(protein_keys):
            group_tsv = os.path.join(tmpd, f"{_safe_filename(genome_accession)}.tsv")
            group_artifacts_dir = None
            if artifacts_dir is not None:
                group_artifacts_dir = os.path.join(
                    artifacts_dir,
                    "genomes",
                    _safe_filename(genome_accession),
                )
            try:
                rows.extend(rules.check(
                    genome_keys,
                    group_tsv,
                    artifacts_dir=group_artifacts_dir,
                ))
            finally:
                CuratedProtein.clear_cache()
    write_rule_rows(rows, output_tsv, rules)
    return rows


def write_rule_fasta(rows, fasta_output, include_maybe=True):
    accepted = {RULE_TRUE}
    if include_maybe:
        accepted.add(RULE_MAYBE)

    with open_file_to_write(fasta_output, "wt") as f:
        for row in rows:
            if row["pass all"] not in accepted:
                continue
            protein_accession = row["protein accession"]
            genome_accession = row["genome accession"]
            sequence = CuratedProtein(protein_accession, genome_accession).sequence()
            f.write(f">{protein_accession}\n{sequence}\n")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--rule", required=True)
    parser.add_argument("--artifacts-dir")
    parser.add_argument("--fasta-output")
    parser.add_argument("--fasta-excludes-maybe", action="store_true", default=False)
    args = parser.parse_args(argv)

    rules = load_rules(args.rule)
    protein_keys = filter_manifest_protein_keys(read_protein_keys(sys.stdin))

    with tempfile.TemporaryDirectory() as tmpd:
        if args.artifacts_dir is not None:
            os.makedirs(args.artifacts_dir, exist_ok=True)
            output_tsv = os.path.join(args.artifacts_dir, "rule-results.tsv")
        else:
            output_tsv = os.path.join(tmpd, "rule-results.tsv")

        rows = check_rules_by_genome(
            rules,
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
