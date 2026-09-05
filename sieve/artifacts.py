import csv
import os
import tempfile

from tangle.sequence import read_fasta_as_dict


RULE_RESULTS_TSV = "rule-results.tsv"
SEQUENCES_FASTA = "sequences.faa"
INPUT_FASTA = "input.faa"
SEQUENCES_TSV = "sequences.tsv"

SEQUENCE_HEADERS = [
    "protein accession",
    "genome accession",
    "sequence accession",
    "start label",
    "start aa 1b",
    "protein start aa 1b",
    "end aa 1b",
]
LEGACY_SEQUENCE_HEADERS = SEQUENCE_HEADERS[:5]
PREVIOUS_SEQUENCE_HEADERS = SEQUENCE_HEADERS[:5] + ["end aa 1b"]


def rule_results_tsv(artifacts_dir):
    return os.path.join(artifacts_dir, RULE_RESULTS_TSV)


def sequences_fasta(artifacts_dir):
    return os.path.join(artifacts_dir, SEQUENCES_FASTA)


def input_fasta(artifacts_dir):
    return os.path.join(artifacts_dir, INPUT_FASTA)


def sequences_tsv(artifacts_dir):
    return os.path.join(artifacts_dir, SEQUENCES_TSV)


def _read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_sequence_rows(artifacts_dir):
    path = sequences_tsv(artifacts_dir)
    if not os.path.exists(path):
        return []
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames not in (SEQUENCE_HEADERS, PREVIOUS_SEQUENCE_HEADERS, LEGACY_SEQUENCE_HEADERS):
            raise ValueError(f"Unexpected columns in {path}: {reader.fieldnames}")
        rows = list(reader)
        if reader.fieldnames == LEGACY_SEQUENCE_HEADERS:
            for row in rows:
                row["end aa 1b"] = ""
        if reader.fieldnames != SEQUENCE_HEADERS:
            for row in rows:
                row["protein start aa 1b"] = ""
        return rows


def _atomic_write(path, write):
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=directory, prefix=".artifact-", suffix=".tmp")
    os.close(fd)
    try:
        write(temporary)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


def _write_fasta(path, sequences):
    def write(temporary):
        with open(temporary, "w", encoding="utf-8") as f:
            for accession, sequence in sequences.items():
                f.write(f">{accession}\n{sequence}\n")
    _atomic_write(path, write)


def _merge_sequences(existing_path, additions, label):
    existing = read_fasta_as_dict(existing_path) if os.path.exists(existing_path) else {}
    merged = dict(existing)
    for accession, sequence in additions.items():
        if accession in merged and merged[accession] != sequence:
            raise ValueError(f"Conflicting {label} sequence for accession {accession}")
        merged[accession] = sequence
    return merged


def merge_sequence_artifacts(artifacts_dir, originals, candidate_entries):
    """Validate and merge original proteins and discovered candidate sequences."""
    input_path = input_fasta(artifacts_dir)
    candidate_path = sequences_fasta(artifacts_dir)
    row_path = sequences_tsv(artifacts_dir)

    candidate_sequences = {}
    new_rows = []
    for entry in candidate_entries:
        protein_accession = entry["row"]["protein accession"]
        genome_accession = entry["row"].get("genome accession", "")
        if not entry["candidates"]:
            new_rows.append({
                "protein accession": protein_accession,
                "genome accession": genome_accession,
                "sequence accession": "",
                "start label": "",
                "start aa 1b": "",
                "protein start aa 1b": "",
                "end aa 1b": "",
            })
        for candidate in entry["candidates"]:
            previous = candidate_sequences.get(candidate.accession)
            if previous is not None and previous != candidate.sequence:
                raise ValueError(f"Conflicting candidate sequence for accession {candidate.accession}")
            candidate_sequences[candidate.accession] = candidate.sequence
            new_rows.append({
                "protein accession": protein_accession,
                "genome accession": genome_accession,
                "sequence accession": candidate.accession,
                "start label": candidate.start_label,
                "start aa 1b": str(candidate.start_aa_1b),
                "protein start aa 1b": (
                    "" if candidate.protein_start_aa_1b is None
                    else str(candidate.protein_start_aa_1b)
                ),
                "end aa 1b": "" if candidate.end_aa_1b is None else str(candidate.end_aa_1b),
            })

    merged_originals = _merge_sequences(input_path, originals, "input")
    merged_candidates = _merge_sequences(candidate_path, candidate_sequences, "candidate")
    existing_rows = read_sequence_rows(artifacts_dir)
    rows_by_identity = {
        row["sequence accession"] or (row["protein accession"], row["genome accession"], ""):
        row
        for row in existing_rows
    }
    for row in new_rows:
        accession = row["sequence accession"]
        identity = accession or (row["protein accession"], row["genome accession"], "")
        if identity in rows_by_identity and rows_by_identity[identity] != row:
            previous = rows_by_identity[identity]
            upgrade = dict(previous)
            if not upgrade.get("protein start aa 1b"):
                upgrade["protein start aa 1b"] = row["protein start aa 1b"]
            if upgrade != row:
                raise ValueError(f"Conflicting candidate metadata for accession {accession}")
        rows_by_identity[identity] = row

    genome_by_protein = {}
    for row in rows_by_identity.values():
        if row["protein accession"] not in merged_originals:
            raise ValueError(f"Candidate refers to missing input protein {row['protein accession']}")
        if row["sequence accession"] and row["sequence accession"] not in merged_candidates:
            raise ValueError(f"Candidate metadata refers to missing sequence {row['sequence accession']}")
        protein_accession = row["protein accession"]
        genome_accession = row["genome accession"]
        previous_genome = genome_by_protein.setdefault(protein_accession, genome_accession)
        if previous_genome != genome_accession:
            raise ValueError(f"Protein accession occurs in multiple genomes: {protein_accession}")
        sequence_accession = row["sequence accession"]
        if sequence_accession in merged_originals and sequence_accession != protein_accession:
            raise ValueError(
                f"Candidate accession conflicts with input protein accession {sequence_accession}"
            )

    ordered_rows = list(rows_by_identity.values())
    _write_fasta(input_path, merged_originals)
    _write_fasta(candidate_path, merged_candidates)

    def write_rows(temporary):
        with open(temporary, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=SEQUENCE_HEADERS, delimiter="\t")
            writer.writeheader()
            writer.writerows(ordered_rows)
    _atomic_write(row_path, write_rows)
