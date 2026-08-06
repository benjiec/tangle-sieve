#!/usr/bin/env python3

import argparse
import os
import sys
from collections import defaultdict

from tangle.sequence import read_fasta_as_dict

from sieve.artifact_builder import build_artifacts
from sieve.fasta_protein import FastaProtein
from sieve.hmmsearch import detected_rows_from_hits, parse_domtblout, run_hmmsearch
from sieve.rule_loader import load_rules


PFAM_DOMTBLOUT = "pfam.domtblout"
KO_DOMTBLOUT = "ko.domtblout"


def _rows_by_query(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["query_accession"]].append(row)
    return grouped


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("-r", "--rule", required=True)
    parser.add_argument("--artifacts-dir", required=True)
    parser.add_argument("--pfam-hmm")
    parser.add_argument("--ko-hmm")
    args = parser.parse_args(argv)

    rules = load_rules(args.rule)
    pfam_rows = []
    ko_rows = []
    if args.pfam_hmm is not None:
        os.makedirs(args.artifacts_dir, exist_ok=True)
        domtblout = os.path.join(args.artifacts_dir, PFAM_DOMTBLOUT)
        run_hmmsearch(args.pfam_hmm, args.fasta, domtblout)
        pfam_rows = detected_rows_from_hits(parse_domtblout(domtblout), "Pfam")
    if rules.cterm_bound_rules():
        if args.ko_hmm is None:
            raise ValueError("--ko-hmm is required for KO.matches(..., bound_cterm=True)")
        os.makedirs(args.artifacts_dir, exist_ok=True)
        domtblout = os.path.join(args.artifacts_dir, KO_DOMTBLOUT)
        run_hmmsearch(args.ko_hmm, args.fasta, domtblout, use_cut_ga=False)
        ko_rows = detected_rows_from_hits(parse_domtblout(domtblout), "KO")
    pfam_by_query = _rows_by_query(pfam_rows)
    ko_by_query = _rows_by_query(ko_rows)
    proteins = [
        FastaProtein(
            accession,
            sequence,
            pfam_rows=pfam_by_query.get(accession, []),
            ko_rows=ko_by_query.get(accession, []),
        )
        for accession, sequence in read_fasta_as_dict(args.fasta).items()
    ]
    build_artifacts(proteins, rules, args.artifacts_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
