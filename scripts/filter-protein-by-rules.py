#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys
import tempfile
from collections import defaultdict

from sieve.artifacts import rule_results_tsv, sequences_fasta
from sieve.protein import CuratedProtein, SEQUENCE_SOURCE_HMM_DETECTED
from sieve.rule_artifacts import write_rule_fasta, write_rule_rows
from sieve.rule_loader import load_rules


GENOMIC_LOCI_TSV = "genomic_loci.tsv"


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


def _group_candidate_entries_by_protein(candidate_entries):
    grouped = {}
    ordered_keys = []
    for entry in candidate_entries:
        row = entry["row"]
        key = (row["protein accession"], row["genome accession"])
        if key not in grouped:
            grouped[key] = {
                "row": row,
                "candidates": [],
                "seen_accessions": set(),
            }
            ordered_keys.append(key)
        group = grouped[key]
        for candidate in entry["candidates"]:
            if candidate.accession in group["seen_accessions"]:
                continue
            group["candidates"].append(candidate)
            group["seen_accessions"].add(candidate.accession)
    return [grouped[key] for key in ordered_keys]


def _protein_locus_for_candidates(protein, candidates):
    if protein.sequence_source == SEQUENCE_SOURCE_HMM_DETECTED:
        if any(candidate.start_label.endswith("_anchor") for candidate in candidates):
            return protein.genomic_locus()
        start_candidates = [
            candidate for candidate in candidates
            if candidate.start_label
        ]
        earliest_start_aa_1b = min(
            [candidate.start_aa_1b for candidate in start_candidates],
            default=1,
        )
        return protein._hmm_detected_locus(
            leader_prefix_len=max(0, 1 - earliest_start_aa_1b),
            start_trim_len=0,
        )
    return protein.genomic_locus_with_leader()


def _candidate_start_feature_position(protein, locus, candidate, earliest_start_aa_1b):
    if protein.sequence_source == SEQUENCE_SOURCE_HMM_DETECTED:
        if not candidate.start_label:
            return None
        return (candidate.start_aa_1b - earliest_start_aa_1b) * 3 + 1
    if not candidate.start_label:
        return locus.start_codon_position_1b()
    if candidate.start_aa_1b >= 1:
        try:
            codon_start, _codon_end = protein.protein_codon_interval_1b(candidate.start_aa_1b)
            return abs(codon_start - locus.start_1b) + 1
        except Exception:
            return None
    return None


def write_locus_artifact(candidate_entries, output_tsv):
    headers = [
        "protein accession",
        "genome accession",
        "contig accession",
        "sequence accession",
        "locus start 1b",
        "locus end 1b",
        "strand",
        "sequence length",
        "feature type",
        "feature index",
        "feature position 1b",
        "error",
    ]
    with open(output_tsv, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        grouped_entries = _group_candidate_entries_by_protein(candidate_entries)
        entries_by_genome = defaultdict(list)
        for entry in grouped_entries:
            entries_by_genome[entry["row"]["genome accession"]].append(entry)
        for genome_accession in sorted(entries_by_genome):
            try:
                for entry in entries_by_genome[genome_accession]:
                    protein_accession = entry["row"]["protein accession"]
                    candidates = entry["candidates"]
                    if not candidates:
                        continue
                    row = {
                        "protein accession": protein_accession,
                        "genome accession": genome_accession,
                    }
                    try:
                        protein = CuratedProtein(protein_accession, genome_accession)
                        locus = _protein_locus_for_candidates(protein, candidates)
                        start_candidates = [
                            candidate for candidate in candidates
                            if candidate.start_label
                        ]
                        earliest_start_aa_1b = min(
                            [candidate.start_aa_1b for candidate in start_candidates],
                            default=1,
                        )
                        base_row = row | {
                            "contig accession": locus.contig_accession,
                            "locus start 1b": locus.start_1b,
                            "locus end 1b": locus.end_1b,
                            "strand": "+" if locus.strand == 1 else "-",
                            "sequence length": len(locus.sequence()),
                            "error": "",
                        }
                        feature_rows = []
                        for candidate in candidates:
                            position = _candidate_start_feature_position(
                                protein,
                                locus,
                                candidate,
                                earliest_start_aa_1b,
                            )
                            if position is not None:
                                feature_rows.append(base_row | {
                                    "feature type": "start",
                                    "feature index": "",
                                    "sequence accession": candidate.accession,
                                    "feature position 1b": position,
                                })
                        stop_position = locus.stop_codon_position_1b()
                        if stop_position is not None:
                            feature_rows.append(base_row | {
                                "feature type": "stop",
                                "feature index": 1,
                                "sequence accession": "",
                                "feature position 1b": stop_position,
                            })
                        for i, position in enumerate(locus.dss_positions_1b(), start=1):
                            feature_rows.append(base_row | {
                                "feature type": "dss",
                                "feature index": i,
                                "sequence accession": "",
                                "feature position 1b": position,
                            })
                        for i, position in enumerate(locus.ass_positions_1b(), start=1):
                            feature_rows.append(base_row | {
                                "feature type": "ass",
                                "feature index": i,
                                "sequence accession": "",
                                "feature position 1b": position,
                            })
                        if feature_rows:
                            writer.writerows(feature_rows)
                        else:
                            writer.writerow(base_row | {
                                "feature type": "",
                                "feature index": "",
                                "sequence accession": candidates[0].accession,
                                "feature position 1b": "",
                            })
                    except Exception as e:
                        row["error"] = str(e)
                        writer.writerow(row)
            finally:
                CuratedProtein.clear_cache()


def _log_sequence_progress(index, total, row, message):
    protein_accession = row["protein accession"]
    genome_accession = row["genome accession"]
    print(
        f"[sequences {index}/{total}] {protein_accession} {genome_accession}: {message}",
        file=sys.stderr,
    )


def scoped_candidate_entries(rows, rules, trace=True):
    entries = []
    total = len(rows)
    for i, row in enumerate(rows, start=1):
        if trace:
            _log_sequence_progress(i, total, row, "discovering")
        try:
            candidates = rules.scoped_sequence_candidates_for_row(row)
        finally:
            CuratedProtein.clear_cache()
        if trace:
            _log_sequence_progress(i, total, row, f"done: candidates={len(candidates)}")
        entries.append({
            "row": row,
            "candidates": candidates,
        })
    return entries


def write_candidate_entries_fasta(candidate_entries, fasta_output):
    write_rule_fasta(
        candidate_entries,
        fasta_output,
        lambda entry: entry["candidates"],
    )


def sequence_only_rows(protein_keys):
    return [
        {
            "protein accession": protein_accession,
            "genome accession": genome_accession,
        }
        for protein_accession, genome_accession in protein_keys
    ]


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--rule", required=True)
    parser.add_argument("--artifacts-dir", required=True)
    parser.add_argument(
        "--sequences-only",
        action="store_true",
        help=f"write sequences.faa and {GENOMIC_LOCI_TSV} without evaluating rules",
    )
    args = parser.parse_args(argv)

    rules = load_rules(args.rule)
    protein_keys = filter_manifest_protein_keys(read_protein_keys(sys.stdin))

    os.makedirs(args.artifacts_dir, exist_ok=True)
    if args.sequences_only:
        candidate_entries = scoped_candidate_entries(sequence_only_rows(protein_keys), rules)
        write_locus_artifact(
            candidate_entries,
            os.path.join(args.artifacts_dir, GENOMIC_LOCI_TSV),
        )
        write_candidate_entries_fasta(
            candidate_entries,
            sequences_fasta(args.artifacts_dir),
        )
        return 0

    output_tsv = rule_results_tsv(args.artifacts_dir)
    rows = check_rules_by_genome(
        rules,
        protein_keys,
        output_tsv,
        artifacts_dir=args.artifacts_dir,
    )
    candidate_entries = scoped_candidate_entries(rows, rules)
    write_locus_artifact(
        candidate_entries,
        os.path.join(args.artifacts_dir, GENOMIC_LOCI_TSV),
    )
    write_candidate_entries_fasta(
        candidate_entries,
        sequences_fasta(args.artifacts_dir),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
