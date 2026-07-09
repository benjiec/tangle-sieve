#!/usr/bin/env python3

import argparse
import os
import sys
from collections import defaultdict

from tangle.sequence import read_fasta_as_dict

from sieve.artifacts import rule_results_tsv, sequences_fasta
from sieve.fasta_protein import FastaProtein
from sieve.hmmsearch import (
    detected_rows_from_hits,
    parse_domtblout,
    read_ko_thresholds,
    run_hmmsearch,
)
from sieve.rule_artifacts import write_rule_fasta
from sieve.rule_loader import load_rules


PFAM_DOMTBLOUT = "pfam.domtblout"
KO_DOMTBLOUT = "ko.domtblout"


def _rows_by_query(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["query_accession"]].append(row)
    return grouped


def run_pfam_search(pfam_hmm, fasta, artifacts_dir):
    if pfam_hmm is None:
        return []
    domtblout = os.path.join(artifacts_dir, PFAM_DOMTBLOUT)
    run_hmmsearch(pfam_hmm, fasta, domtblout, use_cut_ga=True)
    return detected_rows_from_hits(parse_domtblout(domtblout), "Pfam")


def run_ko_search(ko_hmm, fasta, ko_thresholds, artifacts_dir):
    if ko_hmm is None:
        return []
    if ko_thresholds is None:
        raise ValueError("--ko-thresholds is required when --ko-hmm is provided")
    domtblout = os.path.join(artifacts_dir, KO_DOMTBLOUT)
    run_hmmsearch(ko_hmm, fasta, domtblout, use_cut_ga=False)
    return detected_rows_from_hits(
        parse_domtblout(domtblout),
        "KO",
        threshold_by_model=read_ko_thresholds(ko_thresholds),
    )


def build_fasta_proteins(fasta_path, pfam_rows, ko_rows):
    sequences = read_fasta_as_dict(fasta_path)
    pfam_by_query = _rows_by_query(pfam_rows)
    ko_by_query = _rows_by_query(ko_rows)
    return [
        FastaProtein(
            accession,
            sequence,
            pfam_rows=pfam_by_query.get(accession, []),
            ko_rows=ko_by_query.get(accession, []),
        )
        for accession, sequence in sequences.items()
    ]


def write_fasta_artifact(rows, proteins_by_accession, rules, output_fasta):
    def candidates_for_row(row):
        protein = proteins_by_accession[row["protein accession"]]
        return rules.scoped_sequence_candidates_for_protein_row(protein, row)
    write_rule_fasta(rows, output_fasta, candidates_for_row)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("-r", "--rule", required=True)
    parser.add_argument("--artifacts-dir", required=True)
    parser.add_argument("--pfam-hmm")
    parser.add_argument("--ko-hmm")
    parser.add_argument("--ko-thresholds")
    args = parser.parse_args(argv)

    os.makedirs(args.artifacts_dir, exist_ok=True)

    rules = load_rules(args.rule)
    pfam_rows = run_pfam_search(args.pfam_hmm, args.fasta, args.artifacts_dir)
    ko_rows = run_ko_search(args.ko_hmm, args.fasta, args.ko_thresholds, args.artifacts_dir)
    proteins = build_fasta_proteins(args.fasta, pfam_rows, ko_rows)
    proteins_by_accession = {
        protein.protein_accession: protein
        for protein in proteins
    }

    rows = rules.check_proteins(
        proteins,
        rule_results_tsv(args.artifacts_dir),
        artifacts_dir=args.artifacts_dir,
    )
    write_fasta_artifact(
        rows,
        proteins_by_accession,
        rules,
        sequences_fasta(args.artifacts_dir),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
