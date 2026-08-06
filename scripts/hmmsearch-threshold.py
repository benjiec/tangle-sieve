#!/usr/bin/env python3

import argparse
import csv
import math
import sys
import tempfile

from tangle.sequence import read_fasta_as_dict

from sieve.artifacts import read_sequence_rows, rule_results_tsv, sequences_fasta
from sieve.hmmsearch import parse_domtblout, run_hmmsearch


OUTPUT_HEADERS = [
    "threshold bitscore",
    "tp",
    "fp",
    "tn",
    "fn",
    "sensitivity",
    "specificity",
    "balanced accuracy",
    "selected",
]


def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def best_hmm_hits(hits):
    best = {}
    for hit in hits:
        accession = hit.sequence_accession
        bitscore = hit.domain_bitscore
        current = best.get(accession)
        if current is None or bitscore > current["bitscore"]:
            best[accession] = {"bitscore": bitscore, "hit": hit}
    return best


def joined_entries(hmm_hits, rule_rows, sequence_rows, fasta_path):
    sequences = read_fasta_as_dict(fasta_path)
    hits_by_sequence = best_hmm_hits(hmm_hits)
    manifest_by_sequence = {}
    for row in sequence_rows:
        accession = row["sequence accession"]
        if not accession:
            continue
        if accession in manifest_by_sequence:
            raise ValueError(f"Duplicate sequence accession in candidate metadata: {accession}")
        manifest_by_sequence[accession] = row
    manifest_accessions = set(manifest_by_sequence)
    if manifest_accessions != set(sequences):
        missing_fasta = sorted(manifest_accessions - set(sequences))
        missing_metadata = sorted(set(sequences) - manifest_accessions)
        if missing_fasta:
            raise ValueError(f"Candidate metadata refers to missing sequences: {', '.join(missing_fasta[:5])}")
        raise ValueError(f"Candidate sequences have no metadata: {', '.join(missing_metadata[:5])}")
    rules_by_sequence = {}
    for row in rule_rows:
        accession = row.get("sequence accession", "")
        if not accession:
            raise ValueError("Rule results row is missing sequence accession")
        if accession in rules_by_sequence:
            raise ValueError(f"Duplicate sequence accession in rule results: {accession}")
        if "pass all" not in row:
            raise ValueError("Rule results are missing pass all column")
        rules_by_sequence[accession] = row

    entries = []
    for accession in sequences:
        rule_row = rules_by_sequence.get(accession)
        if rule_row is None:
            raise ValueError(f"Cannot find FASTA accession in rule results: {accession}")
        manifest_row = manifest_by_sequence[accession]
        for column in ("protein accession", "genome accession"):
            if rule_row.get(column, "") != manifest_row.get(column, ""):
                raise ValueError(f"Candidate metadata and rule results disagree for {accession}: {column}")
        hit = hits_by_sequence.get(accession)
        entries.append({
            "sequence accession": accession,
            "protein accession": rule_row.get("protein accession", ""),
            "genome accession": rule_row.get("genome accession", ""),
            "bitscore": hit["bitscore"] if hit is not None else 0.0,
            "rule_row": rule_row,
            "hmm_hit": hit["hit"] if hit is not None else None,
            "positive": rule_row.get("pass all") == "true",
        })

    extra_hits = sorted(set(hits_by_sequence) - set(sequences))
    if extra_hits:
        raise ValueError(f"HMM search results contain accessions absent from FASTA: {', '.join(extra_hits[:5])}")
    extra_rules = sorted(set(rules_by_sequence) - set(sequences))
    if extra_rules:
        raise ValueError(f"Rule results contain accessions absent from FASTA: {', '.join(extra_rules[:5])}")
    return entries


def _ratio(numerator, denominator):
    if denominator == 0:
        return math.nan
    return numerator / denominator


def _format_metric(value):
    if math.isnan(value):
        return ""
    return f"{value:.6g}"


def threshold_stats(entries):
    if not entries:
        return []
    thresholds = sorted({entry["bitscore"] for entry in entries}, reverse=True)
    rows = []
    for threshold in thresholds:
        tp = fp = tn = fn = 0
        for entry in entries:
            predicted_positive = entry["bitscore"] >= threshold
            actual_positive = entry["positive"]
            if actual_positive and predicted_positive:
                tp += 1
            elif not actual_positive and predicted_positive:
                fp += 1
            elif not actual_positive and not predicted_positive:
                tn += 1
            else:
                fn += 1
        sensitivity = _ratio(tp, tp + fn)
        specificity = _ratio(tn, tn + fp)
        balanced_accuracy = _ratio(sensitivity + specificity, 2)
        false_positives = [
            entry for entry in entries
            if not entry["positive"] and entry["bitscore"] >= threshold
        ]
        false_negatives = [
            entry for entry in entries
            if entry["positive"] and entry["bitscore"] < threshold
        ]
        rows.append({
            "threshold bitscore": threshold,
            "tp": tp,
            "fp": fp,
            "tn": tn,
            "fn": fn,
            "sensitivity": sensitivity,
            "specificity": specificity,
            "balanced accuracy": balanced_accuracy,
            "selected": "",
            "false positives": false_positives,
            "false negatives": false_negatives,
        })
    best = best_threshold_row(rows)
    if best is not None:
        best["selected"] = "true"
    return rows


def best_threshold_row(rows):
    comparable = [
        row for row in rows
        if not math.isnan(row["balanced accuracy"])
    ]
    if not comparable:
        return None
    return max(comparable, key=lambda row: (
        row["balanced accuracy"],
        -math.inf if math.isnan(row["specificity"]) else row["specificity"],
        -math.inf if math.isnan(row["sensitivity"]) else row["sensitivity"],
        row["threshold bitscore"],
    ))


def write_threshold_stats(rows, output):
    stream = sys.stdout
    close = False
    if output is not None:
        stream = open(output, "w", encoding="utf-8", newline="")
        close = True
    try:
        writer = csv.DictWriter(stream, fieldnames=OUTPUT_HEADERS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "threshold bitscore": _format_metric(row["threshold bitscore"]),
                "tp": row["tp"],
                "fp": row["fp"],
                "tn": row["tn"],
                "fn": row["fn"],
                "sensitivity": _format_metric(row["sensitivity"]),
                "specificity": _format_metric(row["specificity"]),
                "balanced accuracy": _format_metric(row["balanced accuracy"]),
                "selected": row["selected"],
            })
        selected = selected_threshold_row(rows)
        if selected is not None:
            write_selected_error_details(stream, selected)
    finally:
        if close:
            stream.close()


def selected_threshold_row(rows):
    for row in rows:
        if row["selected"]:
            return row
    return None


def write_selected_error_details(stream, row):
    stream.write(f"# selected threshold bitscore {_format_metric(row['threshold bitscore'])}\n")
    stream.write("# false positives\n")
    for entry in row["false positives"]:
        stream.write(f"# {entry['sequence accession']} {_format_metric(entry['bitscore'])}\n")
    stream.write("# false negatives\n")
    for entry in row["false negatives"]:
        stream.write(f"# {entry['sequence accession']} {_format_metric(entry['bitscore'])}\n")


def discover_threshold(hmm_profile, artifacts_dir):
    with tempfile.TemporaryDirectory() as tmpd:
        domtblout = f"{tmpd}/hmmsearch.domtblout"
        run_hmmsearch(
            hmm_profile,
            sequences_fasta(artifacts_dir),
            domtblout,
            use_cut_ga=False,
        )
        hits = parse_domtblout(domtblout)
    entries = joined_entries(
        hits,
        read_tsv(rule_results_tsv(artifacts_dir)),
        read_sequence_rows(artifacts_dir),
        sequences_fasta(artifacts_dir),
    )
    return threshold_stats(entries)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--hmm", required=True)
    parser.add_argument("--artifacts-dir", required=True)
    parser.add_argument("-o", "--output")
    args = parser.parse_args(argv)

    rows = discover_threshold(
        args.hmm,
        args.artifacts_dir,
    )
    write_threshold_stats(rows, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
