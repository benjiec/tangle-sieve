#!/usr/bin/env python3

import argparse
import csv
import os
import subprocess
import sys
import tempfile

from sieve.artifacts import sequences_fasta


HEADERS = [
    "sequence accession",
    "HMM model",
    "domain e-value",
    "domain bitscore",
    "query start",
    "query end",
    "hmm start",
    "hmm end",
]


def parse_domtblout(path):
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(maxsplit=22)
            if len(parts) < 22:
                raise ValueError(f"Expected at least 22 domtblout columns, got: {raw_line.rstrip()}")
            rows.append({
                "sequence accession": parts[0],
                "HMM model": parts[3],
                "domain e-value": parts[12],
                "domain bitscore": parts[13],
                "query start": parts[17],
                "query end": parts[18],
                "hmm start": parts[15],
                "hmm end": parts[16],
            })
    return rows


def write_tsv(rows, output_tsv):
    with open(output_tsv, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=HEADERS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run_hmmsearch(artifacts_dir, hmm_file, output_tsv):
    sequences_faa = sequences_fasta(artifacts_dir)
    with tempfile.TemporaryDirectory() as tmpd:
        domtblout = os.path.join(tmpd, "hmmsearch.domtblout")
        cmd = ["hmmsearch", "--domtblout", domtblout, hmm_file, sequences_faa]
        completed = subprocess.run(cmd, check=False, capture_output=True, text=True)
        if completed.returncode != 0:
            print("hmmsearch failed:", " ".join(cmd), file=sys.stderr)
            if completed.stderr:
                print(completed.stderr, end="" if completed.stderr.endswith("\n") else "\n", file=sys.stderr)
            if completed.stdout:
                print(completed.stdout, end="" if completed.stdout.endswith("\n") else "\n", file=sys.stderr)
            raise SystemExit(completed.returncode)
        rows = parse_domtblout(domtblout)
    write_tsv(rows, output_tsv)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifacts-dir", required=True)
    parser.add_argument("--hmm", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args(argv)

    run_hmmsearch(args.artifacts_dir, args.hmm, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
