import csv
import os
import sys

from tangle.sequence import read_fasta_as_dict

from sieve.artifacts import (
    _atomic_write,
    input_fasta,
    merge_sequence_artifacts,
    read_sequence_rows,
    sequences_fasta,
)
from sieve.protein import LeaderSequenceCandidate, SEQUENCE_SOURCE_HMM_DETECTED


GENOMIC_LOCI_TSV = "genomic_loci.tsv"
LOCUS_HEADERS = [
    "protein accession", "genome accession", "contig accession", "sequence accession",
    "locus start 1b", "locus end 1b", "strand", "sequence length", "feature type",
    "feature index", "feature position 1b", "error",
]


def discover_candidates(proteins, rules):
    entries = []
    total = len(proteins)
    for index, protein in enumerate(proteins, start=1):
        row = {
            "protein accession": protein.protein_accession,
            "genome accession": protein.genome_accession,
        }
        print(
            f"[sequences {index}/{total}] {protein.protein_accession} "
            f"{protein.genome_accession}: discovering",
            file=sys.stderr,
        )
        candidates = rules.scoped_sequence_candidates_for_protein_row(protein, row)
        candidates = rules.bound_sequence_candidates(protein, candidates)
        print(
            f"[sequences {index}/{total}] {protein.protein_accession} "
            f"{protein.genome_accession}: done: candidates={len(candidates)}",
            file=sys.stderr,
        )
        entries.append({"row": row, "candidates": candidates})
    return entries


def _first_genome_proteins(proteins, artifacts_dir):
    existing_sequences = {}
    input_path = input_fasta(artifacts_dir)
    if os.path.exists(input_path):
        existing_sequences = read_fasta_as_dict(input_path)
    genome_by_accession = {}
    for row in read_sequence_rows(artifacts_dir):
        accession = row["protein accession"]
        genome = row["genome accession"]
        previous = genome_by_accession.setdefault(accession, genome)
        if previous != genome:
            raise ValueError(f"Protein accession occurs in multiple genomes: {accession}")

    accepted = []
    sequence_by_accession = dict(existing_sequences)
    for protein in proteins:
        accession = protein.protein_accession
        genome = protein.genome_accession
        sequence = protein.sequence()
        previous_sequence = sequence_by_accession.get(accession)
        if previous_sequence is not None and previous_sequence != sequence:
            raise ValueError(f"Conflicting input sequence for accession {accession}")
        previous_genome = genome_by_accession.get(accession)
        if previous_genome is not None and previous_genome != genome:
            print(
                f"Ignoring {accession} from genome {genome}: "
                f"accession already added from genome {previous_genome}",
                file=sys.stderr,
            )
            continue
        sequence_by_accession[accession] = sequence
        genome_by_accession[accession] = genome
        accepted.append(protein)
    return accepted


def build_artifacts(proteins, rules, artifacts_dir, write_loci=False):
    proteins = _first_genome_proteins(list(proteins), artifacts_dir)
    entries = discover_candidates(proteins, rules)
    originals = {}
    for protein in proteins:
        accession = protein.protein_accession
        sequence = protein.sequence()
        if accession in originals and originals[accession] != sequence:
            raise ValueError(f"Conflicting input sequence for accession {accession}")
        originals[accession] = sequence
    merge_sequence_artifacts(artifacts_dir, originals, entries)
    if write_loci:
        merge_locus_artifact(artifacts_dir, proteins, entries)
    return entries


def _locus_for_candidates(protein, candidates):
    if protein.sequence_source == SEQUENCE_SOURCE_HMM_DETECTED:
        if any(candidate.start_label.endswith("_anchor") for candidate in candidates):
            return protein.genomic_locus()
        starts = [candidate.start_aa_1b for candidate in candidates if candidate.start_label]
        earliest = min(starts, default=1)
        return protein._hmm_detected_locus(leader_prefix_len=max(0, 1 - earliest), start_trim_len=0)
    return protein.genomic_locus_with_leader()


def _locus_rows(protein, candidates):
    base = {
        "protein accession": protein.protein_accession,
        "genome accession": protein.genome_accession,
    }
    if not candidates:
        return []
    try:
        locus = _locus_for_candidates(protein, candidates)
        starts = [candidate.start_aa_1b for candidate in candidates if candidate.start_label]
        earliest = min(starts, default=1)
        base |= {
            "contig accession": locus.contig_accession,
            "locus start 1b": locus.start_1b,
            "locus end 1b": locus.end_1b,
            "strand": "+" if locus.strand == 1 else "-",
            "sequence length": len(locus.sequence()),
            "error": "",
        }
        rows = []
        for candidate in candidates:
            position = None
            if protein.sequence_source == SEQUENCE_SOURCE_HMM_DETECTED:
                if candidate.start_label:
                    position = (candidate.start_aa_1b - earliest) * 3 + 1
            elif not candidate.start_label:
                position = locus.start_codon_position_1b()
            elif candidate.start_aa_1b >= 1:
                try:
                    codon_start, _ = protein.protein_codon_interval_1b(candidate.start_aa_1b)
                    position = abs(codon_start - locus.start_1b) + 1
                except Exception:
                    pass
            if position is not None:
                rows.append(base | {
                    "sequence accession": candidate.accession,
                    "feature type": "start", "feature index": "", "feature position 1b": position,
                })
        feature_sets = [
            ("stop", [locus.stop_codon_position_1b()] if locus.stop_codon_position_1b() else []),
            ("dss", locus.dss_positions_1b()),
            ("ass", locus.ass_positions_1b()),
        ]
        for feature_type, positions in feature_sets:
            for index, position in enumerate(positions, start=1):
                rows.append(base | {
                    "sequence accession": "", "feature type": feature_type,
                    "feature index": index, "feature position 1b": position,
                })
        return rows or [base | {
            "sequence accession": candidates[0].accession,
            "feature type": "", "feature index": "", "feature position 1b": "",
        }]
    except Exception as error:
        return [base | {"error": str(error)}]


def merge_locus_artifact(artifacts_dir, proteins, entries):
    path = os.path.join(artifacts_dir, GENOMIC_LOCI_TSV)
    existing = []
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8", newline="") as f:
            existing = list(csv.DictReader(f, delimiter="\t"))
    candidate_sequences = read_fasta_as_dict(sequences_fasta(artifacts_dir))
    candidates_by_key = {}
    for row in read_sequence_rows(artifacts_dir):
        key = (row["protein accession"], row["genome accession"])
        if not row["sequence accession"]:
            candidates_by_key.setdefault(key, [])
            continue
        candidates_by_key.setdefault(key, []).append(LeaderSequenceCandidate(
            accession=row["sequence accession"],
            start_label=row["start label"],
            start_aa_1b=int(row["start aa 1b"]),
            sequence=candidate_sequences[row["sequence accession"]],
            end_aa_1b=int(row["end aa 1b"]) if row["end aa 1b"] else None,
        ))
    replaced = {(protein.protein_accession, protein.genome_accession) for protein in proteins}
    rows = [
        row for row in existing
        if (row["protein accession"], row["genome accession"]) not in replaced
    ]
    for protein in proteins:
        key = (protein.protein_accession, protein.genome_accession)
        rows.extend(_locus_rows(protein, candidates_by_key[key]))

    def write(temporary):
        with open(temporary, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=LOCUS_HEADERS, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    _atomic_write(path, write)
