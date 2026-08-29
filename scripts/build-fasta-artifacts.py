#!/usr/bin/env python3

import argparse
import os
import sys
import tempfile
from collections import defaultdict

from tangle import open_file_to_read
from tangle.sequence import write_fasta_from_dict

from sieve.artifact_builder import build_artifacts
from sieve.fasta_protein import FastaProtein
from sieve.hmmsearch import detected_rows_from_hits, parse_domtblout, run_hmmsearch
from sieve.rule_loader import load_rules


PFAM_DOMTBLOUT = "pfam.domtblout"
KO_DOMTBLOUT = "ko.domtblout"


def _validate_identifier(identifier):
    if identifier.count("|") != 1:
        raise ValueError(
            f"FASTA identifier must have exactly one '|' in accession|genome_description format: {identifier}"
        )
    accession, genome_description = identifier.split("|")
    if not accession or not genome_description:
        raise ValueError(
            f"FASTA identifier must have nonempty accession and genome description: {identifier}"
        )
    if any(character.isspace() for character in identifier):
        raise ValueError(f"FASTA identifier cannot contain whitespace: {identifier}")


def read_fasta_sequences(paths):
    sequences = {}
    for path in paths:
        current_identifier = None
        sequence_parts = []
        with open_file_to_read(path) as f:
            for raw_line in f:
                line = raw_line.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if current_identifier is not None:
                        sequences[current_identifier] = "".join(sequence_parts)
                    current_identifier = line[1:].strip()
                    _validate_identifier(current_identifier)
                    if current_identifier in sequences:
                        raise ValueError(f"Duplicate FASTA identifier: {current_identifier}")
                    sequence_parts = []
                else:
                    sequence_parts.append(line.strip())
        if current_identifier is not None:
            if current_identifier in sequences:
                raise ValueError(f"Duplicate FASTA identifier: {current_identifier}")
            sequences[current_identifier] = "".join(sequence_parts)
    return sequences


def _rows_by_query(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["query_accession"]].append(row)
    return grouped


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, action="append")
    parser.add_argument("-r", "--rule", required=True)
    parser.add_argument("--artifacts-dir", required=True)
    parser.add_argument("--pfam-hmm")
    parser.add_argument("--ko-hmm")
    args = parser.parse_args(argv)

    sequences = read_fasta_sequences(args.fasta)
    rules = load_rules(args.rule)
    pfam_rows = []
    ko_rows = []
    with tempfile.TemporaryDirectory() as tmpd:
        combined_fasta = os.path.join(tmpd, "combined.faa")
        write_fasta_from_dict(sequences, combined_fasta)
        if args.pfam_hmm is not None:
            os.makedirs(args.artifacts_dir, exist_ok=True)
            domtblout = os.path.join(args.artifacts_dir, PFAM_DOMTBLOUT)
            run_hmmsearch(args.pfam_hmm, combined_fasta, domtblout)
            pfam_rows = detected_rows_from_hits(parse_domtblout(domtblout), "Pfam")
        if rules.cterm_bound_rules():
            if args.ko_hmm is None:
                raise ValueError("--ko-hmm is required for KO.matches(..., bound_cterm=True)")
            os.makedirs(args.artifacts_dir, exist_ok=True)
            domtblout = os.path.join(args.artifacts_dir, KO_DOMTBLOUT)
            run_hmmsearch(args.ko_hmm, combined_fasta, domtblout, use_cut_ga=False)
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
        for accession, sequence in sequences.items()
    ]
    build_artifacts(proteins, rules, args.artifacts_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
