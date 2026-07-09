#!/usr/bin/env python3

import argparse
import csv
import math
import os
import re
import sys

from tangle.defaults import Defaults
from tangle.sequence import read_fasta_as_dict


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
SPEC_OPERATORS = {"eq", "ne", "regex", "not_regex", "num_eq", "num_ne", "gt", "gte", "lt", "lte"}
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


def read_positive_specs(path):
    with open(path, "r", encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    rows = [
        (line_number, row)
        for line_number, row in enumerate(rows, start=1)
        if row and any(cell.strip() for cell in row)
    ]
    if rows and rows[0][1] == ["column_regex", "operator", "value"]:
        rows = rows[1:]
    specs = []
    for line_number, row in rows:
        if len(row) != 3:
            raise ValueError(f"Expected 3 columns in positive spec at line {line_number}, got {len(row)}")
        column_regex, operator, value = [cell.strip() for cell in row]
        if operator not in SPEC_OPERATORS:
            raise ValueError(f"Unknown positive spec operator at line {line_number}: {operator}")
        try:
            compiled = re.compile(column_regex)
        except re.error as e:
            raise ValueError(f"Invalid column regex at line {line_number}: {e}") from e
        specs.append((compiled, operator, value))
    if not specs:
        raise ValueError("Positive spec is empty")
    return specs


def _numeric(value, context):
    try:
        return float(value)
    except (TypeError, ValueError) as e:
        raise ValueError(f"Expected numeric value for {context}, got: {value!r}") from e


def value_matches(actual, operator, expected):
    actual = "" if actual is None else str(actual)
    if operator == "eq":
        return actual == expected
    if operator == "ne":
        return actual != expected
    if operator == "regex":
        return re.search(expected, actual) is not None
    if operator == "not_regex":
        return re.search(expected, actual) is None
    actual_number = _numeric(actual, "rule row")
    expected_number = _numeric(expected, "positive spec")
    if operator == "num_eq":
        return actual_number == expected_number
    if operator == "num_ne":
        return actual_number != expected_number
    if operator == "gt":
        return actual_number > expected_number
    if operator == "gte":
        return actual_number >= expected_number
    if operator == "lt":
        return actual_number < expected_number
    if operator == "lte":
        return actual_number <= expected_number
    raise ValueError(f"Unknown operator: {operator}")


def row_is_positive(row, specs):
    for column_regex, operator, expected in specs:
        matched_column = False
        for column, actual in row.items():
            if column_regex.fullmatch(column) is None:
                continue
            matched_column = True
            if value_matches(actual, operator, expected):
                return True
        if not matched_column:
            raise ValueError(f"Positive spec did not match any rule-results column: {column_regex.pattern}")
    return False


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
        if accession in rules_by_sequence:
            raise ValueError(f"Duplicate rule-results sequence accession: {accession}")
        rules_by_sequence[accession] = row

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


def label_entries(entries, specs):
    labeled = []
    for entry in entries:
        labeled.append(entry | {"positive": row_is_positive(entry["rule_row"], specs)})
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
    finally:
        if close:
            stream.close()


def discover_threshold(hmmsearch_tsv, rule_results_tsv, fasta, spec_tsv, taxon=None):
    entries = joined_entries(read_tsv(hmmsearch_tsv), read_tsv(rule_results_tsv), fasta)
    entries = best_entry_per_protein(entries)
    entries = filter_by_taxon(entries, taxon)
    entries = label_entries(entries, read_positive_specs(spec_tsv))
    return threshold_stats(entries)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--hmmsearch-tsv", required=True)
    parser.add_argument("--rule-results-tsv", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--positive-spec", required=True)
    parser.add_argument("--taxon")
    parser.add_argument("-o", "--output")
    args = parser.parse_args(argv)

    rows = discover_threshold(
        args.hmmsearch_tsv,
        args.rule_results_tsv,
        args.fasta,
        args.positive_spec,
        taxon=args.taxon,
    )
    write_threshold_stats(rows, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
