#!/usr/bin/env python3

import argparse
import csv
import math
import os
import sys

from tangle.defaults import Defaults
from tangle.sequence import read_fasta_as_dict

from sieve.artifacts import rule_results_tsv, sequences_fasta
from sieve.result_filters import load_result_filter


TAXONOMY_FIELDS = {
    "domain",
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
}
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


def normalized(value):
    return str(value).strip().lower()


def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _numeric(value, context):
    try:
        return float(value)
    except (TypeError, ValueError) as e:
        raise ValueError(f"Expected numeric value for {context}, got: {value!r}") from e


def best_hmm_hits(hmm_rows):
    hits = {}
    for row in hmm_rows:
        accession = row.get("sequence accession", "")
        if not accession:
            raise ValueError("HMM search row is missing sequence accession")
        bitscore = _numeric(row.get("domain bitscore"), f"HMM bitscore for {accession}")
        current = hits.get(accession)
        if current is None or bitscore > current["bitscore"]:
            hits[accession] = {"bitscore": bitscore, "row": row}
    return hits


def joined_entries(hmm_rows, rule_rows, fasta_path):
    sequences = read_fasta_as_dict(fasta_path)
    hits_by_sequence = best_hmm_hits(hmm_rows)
    rules_by_sequence = {}
    for row in rule_rows:
        accession = row.get("sequence accession", "")
        if not accession:
            raise ValueError("Rule results row is missing sequence accession")
        rules_by_sequence.setdefault(accession, row)

    entries = []
    for accession in sequences:
        rule_row = rules_by_sequence.get(accession)
        if rule_row is None:
            raise ValueError(f"Cannot find FASTA accession in rule results: {accession}")
        hit = hits_by_sequence.get(accession)
        entries.append({
            "sequence accession": accession,
            "protein accession": rule_row.get("protein accession", ""),
            "genome accession": rule_row.get("genome accession", ""),
            "bitscore": hit["bitscore"] if hit is not None else 0.0,
            "rule_row": rule_row,
            "hmm_row": hit["row"] if hit is not None else {},
        })

    extra_hits = sorted(set(hits_by_sequence) - set(sequences))
    if extra_hits:
        raise ValueError(f"HMM search results contain accessions absent from FASTA: {', '.join(extra_hits[:5])}")
    return entries


def best_entry_per_protein(entries):
    best = {}
    for entry in entries:
        protein_accession = entry["protein accession"]
        if not protein_accession:
            raise ValueError(f"Rule results row is missing protein accession for {entry['sequence accession']}")
        current = best.get(protein_accession)
        if current is None or entry["bitscore"] > current["bitscore"]:
            best[protein_accession] = entry
    return [
        best[entry["protein accession"]]
        for entry in entries
        if best[entry["protein accession"]]["sequence accession"] == entry["sequence accession"]
    ]


def genome_accession(row):
    for key in row:
        if normalized(key).replace("_", " ") in {"genome accession", "accession"}:
            return row[key]
    return ""


def read_taxonomy_rows():
    taxonomy_tsv = Defaults.area_genome_taxon_tsv()
    if not os.path.exists(taxonomy_tsv):
        return {}
    with open(taxonomy_tsv, "r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return {
        genome_accession(row): row
        for row in rows
        if genome_accession(row)
    }


def taxonomy_matches(row, taxon):
    expected = normalized(taxon)
    for key, value in row.items():
        if normalized(key) in TAXONOMY_FIELDS and normalized(value) == expected:
            return True
    return False


def filter_by_taxon(entries, taxon):
    if taxon is None:
        return entries
    taxonomy_by_genome = read_taxonomy_rows()
    return [
        entry for entry in entries
        if taxonomy_matches(taxonomy_by_genome.get(entry["genome accession"], {}), taxon)
    ]


def label_entries(entries, positive_filter):
    if entries:
        positive_filter.validate_columns(entries[0]["rule_row"].keys())
    labeled = []
    for entry in entries:
        labeled.append(entry | {"positive": positive_filter.matches(entry["rule_row"])})
    return labeled


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


def discover_threshold(hmmsearch_tsv, artifacts_dir, positive_filter_spec, taxon=None):
    entries = joined_entries(
        read_tsv(hmmsearch_tsv),
        read_tsv(rule_results_tsv(artifacts_dir)),
        sequences_fasta(artifacts_dir),
    )
    entries = best_entry_per_protein(entries)
    entries = filter_by_taxon(entries, taxon)
    entries = label_entries(entries, load_result_filter(positive_filter_spec))
    return threshold_stats(entries)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--hmmsearch-tsv", required=True)
    parser.add_argument("--artifacts-dir", required=True)
    parser.add_argument("--positive-filter", required=True)
    parser.add_argument("--taxon")
    parser.add_argument("-o", "--output")
    args = parser.parse_args(argv)

    rows = discover_threshold(
        args.hmmsearch_tsv,
        args.artifacts_dir,
        args.positive_filter,
        taxon=args.taxon,
    )
    write_threshold_stats(rows, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
