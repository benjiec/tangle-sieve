#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys
import tempfile
from collections import defaultdict

from tangle import open_file_to_write
from sieve.protein import CuratedProtein, SEQUENCE_SOURCE_HMM_DETECTED
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
    rules.write_rows(output_tsv, rows)


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


def write_locus_artifact(protein_keys, output_tsv):
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
        for genome_accession, genome_keys in group_protein_keys_by_genome(protein_keys):
            try:
                for protein_accession, _genome_accession in genome_keys:
                    row = {
                        "protein accession": protein_accession,
                        "genome accession": genome_accession,
                    }
                    try:
                        protein = CuratedProtein(protein_accession, genome_accession)
                        candidates = protein.sequences_with_leader()
                        if protein.sequence_source == SEQUENCE_SOURCE_HMM_DETECTED:
                            start_candidates = [
                                candidate for candidate in candidates
                                if candidate.start_label
                            ]
                            earliest_start_aa_1b = min(
                                [candidate.start_aa_1b for candidate in start_candidates],
                                default=1,
                            )
                            locus = protein._hmm_detected_locus(
                                leader_prefix_len=max(0, 1 - earliest_start_aa_1b),
                                start_trim_len=0,
                            )
                        else:
                            locus = protein.genomic_locus_with_leader()
                        base_row = row | {
                            "contig accession": locus.contig_accession,
                            "locus start 1b": locus.start_1b,
                            "locus end 1b": locus.end_1b,
                            "strand": "+" if locus.strand == 1 else "-",
                            "sequence length": len(locus.sequence()),
                            "error": "",
                        }
                        feature_rows = []
                        if protein.sequence_source == SEQUENCE_SOURCE_HMM_DETECTED:
                            for candidate in start_candidates:
                                feature_rows.append(base_row | {
                                    "feature type": "start",
                                    "feature index": "",
                                    "sequence accession": candidate.accession,
                                    "feature position 1b": (candidate.start_aa_1b - earliest_start_aa_1b) * 3 + 1,
                                })
                        else:
                            start_position = locus.start_codon_position_1b()
                            if start_position is not None:
                                feature_rows.append(base_row | {
                                    "feature type": "start",
                                    "feature index": "",
                                    "sequence accession": protein_accession,
                                    "feature position 1b": start_position,
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
                                "sequence accession": protein_accession,
                                "feature position 1b": "",
                            })
                    except Exception as e:
                        row["error"] = str(e)
                        writer.writerow(row)
            finally:
                CuratedProtein.clear_cache()


def write_unfiltered_rule_fasta(rows, fasta_output, rules):
    with open_file_to_write(fasta_output, "wt") as f:
        for row in rows:
            try:
                for candidate in rules.scoped_sequence_candidates_for_row(row):
                    f.write(f">{candidate.accession}\n{candidate.sequence}\n")
            finally:
                CuratedProtein.clear_cache()


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--rule", required=True)
    parser.add_argument("--artifacts-dir", required=True)
    args = parser.parse_args(argv)

    rules = load_rules(args.rule)
    protein_keys = filter_manifest_protein_keys(read_protein_keys(sys.stdin))

    os.makedirs(args.artifacts_dir, exist_ok=True)
    output_tsv = os.path.join(args.artifacts_dir, "rule-results.tsv")
    rows = check_rules_by_genome(
        rules,
        protein_keys,
        output_tsv,
        artifacts_dir=args.artifacts_dir,
    )
    write_locus_artifact(
        protein_keys,
        os.path.join(args.artifacts_dir, "genomic_locus_with_leader.tsv"),
    )
    write_unfiltered_rule_fasta(
        rows,
        os.path.join(args.artifacts_dir, "sequences.faa"),
        rules,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
