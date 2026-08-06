#!/usr/bin/env python3

import argparse
import os
import sys
from collections import defaultdict

from tangle.sequence import read_fasta_as_dict

from sieve.artifact_protein import ArtifactProtein
from sieve.artifacts import input_fasta, read_sequence_rows, rule_results_tsv, sequences_fasta
from sieve.hmm_profiles import hmm_profiles_in_dirs
from sieve.hmmsearch import detected_rows_from_hits, parse_domtblout, read_ko_thresholds, run_hmmsearch
from sieve.protein import CuratedProtein, LeaderSequenceCandidate
from sieve.rule_loader import load_rules


PFAM_DOMTBLOUT = "pfam.domtblout"
KO_DOMTBLOUT = "ko.domtblout"


def _rows_by_query(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["query_accession"]].append(row)
    return grouped


def _uses_rule(rules, prefix):
    return any(rule.label.startswith(prefix) for rule in rules.atomic_rules())


def _run_detection_searches(args, rules):
    input_path = input_fasta(args.artifacts_dir)
    pfam_rows = []
    ko_rows = []
    if _uses_rule(rules, "Pfam."):
        if args.pfam_hmm is None:
            raise ValueError("--pfam-hmm is required for Pfam rules")
        path = os.path.join(args.artifacts_dir, PFAM_DOMTBLOUT)
        run_hmmsearch(args.pfam_hmm, input_path, path)
        pfam_rows = detected_rows_from_hits(parse_domtblout(path), "Pfam")
    if _uses_rule(rules, "KO."):
        if args.ko_hmm is None:
            raise ValueError("--ko-hmm is required for KO rules")
        if args.ko_thresholds is None:
            raise ValueError("--ko-thresholds is required for KO rules")
        path = os.path.join(args.artifacts_dir, KO_DOMTBLOUT)
        run_hmmsearch(args.ko_hmm, input_path, path, use_cut_ga=False)
        ko_rows = detected_rows_from_hits(
            parse_domtblout(path),
            "KO",
            threshold_by_model=read_ko_thresholds(args.ko_thresholds),
        )
    return pfam_rows, ko_rows


def load_artifact_proteins(artifacts_dir, pfam_rows, ko_rows):
    originals = read_fasta_as_dict(input_fasta(artifacts_dir))
    candidate_sequences = read_fasta_as_dict(sequences_fasta(artifacts_dir))
    manifest_rows = read_sequence_rows(artifacts_dir)
    pfam_by_query = _rows_by_query(pfam_rows)
    ko_by_query = _rows_by_query(ko_rows)

    metadata_by_protein = {}
    candidates_by_key = defaultdict(list)
    seen_candidates = set()
    for row in manifest_rows:
        protein_accession = row["protein accession"]
        genome_accession = row["genome accession"]
        sequence_accession = row["sequence accession"]
        if protein_accession not in originals:
            raise ValueError(f"Cannot find input protein sequence {protein_accession}")
        previous_genome = metadata_by_protein.setdefault(protein_accession, genome_accession)
        if previous_genome != genome_accession:
            raise ValueError(f"Protein accession occurs in multiple genomes: {protein_accession}")
        key = (protein_accession, genome_accession)
        candidates_by_key.setdefault(key, [])
        if not sequence_accession:
            continue
        if sequence_accession not in candidate_sequences:
            raise ValueError(f"Cannot find candidate sequence {sequence_accession}")
        if sequence_accession in seen_candidates:
            raise ValueError(f"Duplicate candidate accession in manifest: {sequence_accession}")
        seen_candidates.add(sequence_accession)
        candidates_by_key[key].append(LeaderSequenceCandidate(
            accession=sequence_accession,
            start_label=row["start label"],
            start_aa_1b=int(row["start aa 1b"]),
            sequence=candidate_sequences[sequence_accession],
        ))

    missing_metadata = set(originals) - set(metadata_by_protein)
    if missing_metadata:
        raise ValueError(f"Input proteins have no candidate metadata: {', '.join(sorted(missing_metadata)[:5])}")
    extra_sequences = set(candidate_sequences) - seen_candidates
    if extra_sequences:
        raise ValueError(f"Candidate sequences have no metadata: {', '.join(sorted(extra_sequences)[:5])}")

    proteins = [
        ArtifactProtein(
            accession,
            metadata_by_protein[accession],
            sequence,
            pfam_rows=pfam_by_query.get(accession, []),
            ko_rows=ko_by_query.get(accession, []),
        )
        for accession, sequence in originals.items()
    ]
    return proteins, candidates_by_key


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--rule", required=True)
    parser.add_argument("--artifacts-dir", required=True)
    parser.add_argument("--pfam-hmm")
    parser.add_argument("--ko-hmm")
    parser.add_argument("--ko-thresholds")
    parser.add_argument("--hmm-dir", action="append", default=[])
    parser.add_argument("--deeploc-csv")
    args = parser.parse_args(argv)

    rules = load_rules(args.rule)
    pfam_rows, ko_rows = _run_detection_searches(args, rules)
    proteins, candidates_by_key = load_artifact_proteins(args.artifacts_dir, pfam_rows, ko_rows)
    profiles = hmm_profiles_in_dirs(args.hmm_dir)
    profiles.extend(path for path in (args.pfam_hmm, args.ko_hmm) if path is not None)
    try:
        rules.check_proteins(
            proteins,
            rule_results_tsv(args.artifacts_dir),
            artifacts_dir=args.artifacts_dir,
            deeploc_csv=args.deeploc_csv,
            hmm_profiles=profiles,
            sequence_candidates_by_key=candidates_by_key,
        )
    finally:
        CuratedProtein.clear_cache()
    return 0


if __name__ == "__main__":
    sys.exit(main())
