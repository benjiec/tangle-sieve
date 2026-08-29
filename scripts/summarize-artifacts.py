#!/usr/bin/env python3

import argparse
import csv
import os
import subprocess
import sys
import tempfile

from sieve.artifacts import rule_results_tsv, sequences_fasta
from sieve.hmmsearch import parse_domtblout


RULE_RESULT_COLUMNS = [
    "protein accession",
    "sequence accession",
    "genome accession",
    "contig accession",
    "pass all",
]
SUMMARY_COLUMNS = RULE_RESULT_COLUMNS + [
    "category",
    "hmm model",
    "highest domain score",
    "genome description",
]


def read_rule_results(artifacts_dir):
    path = rule_results_tsv(artifacts_dir)
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = [column for column in RULE_RESULT_COLUMNS if column not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Missing required columns in {path}: {', '.join(missing)}")
        return list(reader)


def highest_domain_scores(artifacts_dir, hmm_file):
    with tempfile.TemporaryDirectory() as tmpd:
        domtblout = os.path.join(tmpd, "hmmsearch.domtblout")
        cmd = ["hmmsearch", "--domtblout", domtblout, hmm_file, sequences_fasta(artifacts_dir)]
        completed = subprocess.run(cmd, check=False, capture_output=True, text=True)
        if completed.returncode != 0:
            print("hmmsearch failed:", " ".join(cmd), file=sys.stderr)
            if completed.stderr:
                print(completed.stderr, end="" if completed.stderr.endswith("\n") else "\n", file=sys.stderr)
            if completed.stdout:
                print(completed.stdout, end="" if completed.stdout.endswith("\n") else "\n", file=sys.stderr)
            raise SystemExit(completed.returncode)

        scores = {}
        for hit in parse_domtblout(domtblout):
            previous = scores.get(hit.sequence_accession)
            if previous is None or hit.domain_bitscore > previous:
                scores[hit.sequence_accession] = hit.domain_bitscore
        return scores


def write_summary(rows, category, hmm_file, scores, output):
    output_exists = os.path.exists(output) and os.path.getsize(output) > 0
    if output_exists:
        with open(output, "r", encoding="utf-8", newline="") as f:
            fieldnames = csv.DictReader(f, delimiter="\t").fieldnames
        if fieldnames != SUMMARY_COLUMNS:
            raise ValueError(f"Unexpected columns in {output}: {fieldnames}")

    with open(output, "a", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_COLUMNS, delimiter="\t")
        if not output_exists:
            writer.writeheader()
        for row in rows:
            sequence_accession = row["sequence accession"]
            score = scores.get(sequence_accession)
            summary_row = {column: row[column] for column in RULE_RESULT_COLUMNS}
            summary_row.update({
                "category": category,
                "hmm model": os.path.basename(hmm_file),
                "highest domain score": "" if score is None else str(score),
                "genome description": (
                    sequence_accession.split("|", 1)[1]
                    if not row["genome accession"] and "|" in sequence_accession
                    else ""
                ),
            })
            writer.writerow(summary_row)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifacts-dir", required=True)
    parser.add_argument("--category", required=True)
    parser.add_argument("--hmm", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)

    rows = read_rule_results(args.artifacts_dir)
    scores = highest_domain_scores(args.artifacts_dir, args.hmm)
    write_summary(rows, args.category, args.hmm, scores, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
