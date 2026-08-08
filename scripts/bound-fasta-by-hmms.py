#!/usr/bin/env python3

import argparse
import os
import sys
import tempfile
from dataclasses import dataclass

from tangle import open_file_to_read

from sieve.hmmsearch import run_hmmsearch
from sieve.protein import hmm_align_sequences


@dataclass(frozen=True)
class DomainHit:
    sequence_accession: str
    model_name: str
    model_length: int
    hmm_start: int
    hmm_end: int
    sequence_start: int
    sequence_end: int


def positive_int(value):
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("position must be at least 1")
    return parsed


def accession_description(value):
    if any(character.isspace() for character in value):
        raise argparse.ArgumentTypeError("accession description cannot contain whitespace")
    return value


def read_unique_fasta(path):
    sequences = {}
    accession = None
    parts = []

    def store():
        if accession is None:
            return
        if accession in sequences:
            raise ValueError(f"Duplicate FASTA accession: {accession}")
        sequences[accession] = "".join(parts)

    with open_file_to_read(path) as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                store()
                header = line[1:].strip()
                if not header:
                    raise ValueError("FASTA header has no accession")
                accession = header.split(None, 1)[0]
                parts = []
            else:
                if accession is None:
                    raise ValueError("FASTA sequence occurs before its first header")
                parts.append(line)
    store()
    return sequences


def parse_domtblout(path):
    hits = []
    with open(path, "r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(maxsplit=22)
            if len(parts) < 22:
                raise ValueError(f"Expected at least 22 domtblout columns, got: {raw_line.rstrip()}")
            hits.append(DomainHit(
                sequence_accession=parts[0],
                model_name=parts[3],
                model_length=int(parts[5]),
                hmm_start=int(parts[15]),
                hmm_end=int(parts[16]),
                sequence_start=int(parts[17]),
                sequence_end=int(parts[18]),
            ))
    return hits


def search(hmm, fasta):
    with tempfile.TemporaryDirectory() as tmpd:
        domtblout = os.path.join(tmpd, "hmmsearch.domtblout")
        run_hmmsearch(hmm, fasta, domtblout, use_cut_ga=True)
        return parse_domtblout(domtblout)


def _hits_by_sequence(hits, known_accessions, label):
    by_sequence = {}
    for hit in hits:
        if hit.sequence_accession not in known_accessions:
            raise ValueError(
                f"{label} HMM result refers to unknown FASTA accession: {hit.sequence_accession}"
            )
        by_sequence.setdefault(hit.sequence_accession, []).append(hit)
    return by_sequence


def qualifying_hits(sequences, n_hits, c_hits, n_position=1, c_position=None, report=None):
    if n_position < 1 or (c_position is not None and c_position < 1):
        raise ValueError("HMM positions must be at least 1")
    report = report or (lambda message: print(message, file=sys.stderr))
    n_by_sequence = _hits_by_sequence(n_hits, sequences, "N-terminal")
    c_by_sequence = _hits_by_sequence(c_hits, sequences, "C-terminal")
    qualifying = {}

    for accession, sequence in sequences.items():
        sequence_n_hits = n_by_sequence.get(accession, [])
        sequence_c_hits = c_by_sequence.get(accession, [])
        if len(sequence_n_hits) != 1:
            report(f"Ignoring {accession}: expected one N-terminal HMM hit, found {len(sequence_n_hits)}")
            continue
        if len(sequence_c_hits) != 1:
            report(f"Ignoring {accession}: expected one C-terminal HMM hit, found {len(sequence_c_hits)}")
            continue

        n_hit = sequence_n_hits[0]
        c_hit = sequence_c_hits[0]
        required_c_position = c_hit.model_length if c_position is None else c_position
        if n_position > n_hit.model_length:
            report(
                f"Ignoring {accession}: N-terminal position {n_position} exceeds "
                f"model length {n_hit.model_length}"
            )
            continue
        if required_c_position > c_hit.model_length:
            report(
                f"Ignoring {accession}: C-terminal position {required_c_position} exceeds "
                f"model length {c_hit.model_length}"
            )
            continue
        if not n_hit.hmm_start <= n_position <= n_hit.hmm_end:
            report(
                f"Ignoring {accession}: N-terminal HMM hit at model positions "
                f"{n_hit.hmm_start}-{n_hit.hmm_end} does not span {n_position}"
            )
            continue
        if not c_hit.hmm_start <= required_c_position <= c_hit.hmm_end:
            report(
                f"Ignoring {accession}: C-terminal HMM hit at model positions "
                f"{c_hit.hmm_start}-{c_hit.hmm_end} does not span {required_c_position}"
            )
            continue
        qualifying[accession] = (n_position, required_c_position)
    return qualifying


def bounded_sequences(sequences, qualifying, n_alignments, c_alignments,
                      acc_desc=None, report=None):
    report = report or (lambda message: print(message, file=sys.stderr))
    bounded = {}
    for accession, (n_position, c_position) in qualifying.items():
        n_mapping = n_alignments[accession].aa_hmm_pos_1b(n_position)
        c_mapping = c_alignments[accession].aa_hmm_pos_1b(c_position)
        if n_mapping is None:
            report(
                f"Ignoring {accession}: N-terminal model position {n_position} "
                "aligns to a deletion"
            )
            continue
        if c_mapping is None:
            report(
                f"Ignoring {accession}: C-terminal model position {c_position} "
                "aligns to a deletion"
            )
            continue

        start = n_mapping[0]
        end = c_mapping[0]
        sequence = sequences[accession]
        if start > end:
            report(f"Ignoring {accession}: N-terminal boundary {start} is after C-terminal boundary {end}")
            continue
        if start < 1 or end > len(sequence):
            report(
                f"Ignoring {accession}: boundaries {start}-{end} are outside sequence length {len(sequence)}"
            )
            continue
        bounded_accession = f"{accession}_bounded_{start}-{end}"
        if acc_desc is not None:
            bounded_accession += f"|{acc_desc}"
        bounded[bounded_accession] = sequence[start - 1:end]
    return bounded


def align_and_bound(sequences, n_hits, c_hits, n_hmm, c_hmm,
                    n_position=1, c_position=None, acc_desc=None, report=None):
    report = report or (lambda message: print(message, file=sys.stderr))
    qualifying = qualifying_hits(
        sequences,
        n_hits,
        c_hits,
        n_position=n_position,
        c_position=c_position,
        report=report,
    )
    sequences_to_align = {
        accession: sequences[accession]
        for accession in qualifying
    }
    n_alignments = hmm_align_sequences(n_hmm, sequences_to_align)
    c_alignments = hmm_align_sequences(c_hmm, sequences_to_align)
    return bounded_sequences(
        sequences,
        qualifying,
        n_alignments,
        c_alignments,
        acc_desc=acc_desc,
        report=report,
    )


def write_fasta(sequences, stream):
    for accession, sequence in sequences.items():
        stream.write(f">{accession}\n{sequence}\n")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-hmm", required=True)
    parser.add_argument("--c-hmm", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--n-position", type=positive_int, default=1)
    parser.add_argument("--c-position", type=positive_int)
    parser.add_argument("--acc-desc", type=accession_description)
    args = parser.parse_args(argv)

    sequences = read_unique_fasta(args.fasta)
    n_hits = search(args.n_hmm, args.fasta)
    c_hits = search(args.c_hmm, args.fasta)
    bounded = align_and_bound(
        sequences,
        n_hits,
        c_hits,
        args.n_hmm,
        args.c_hmm,
        n_position=args.n_position,
        c_position=args.c_position,
        acc_desc=args.acc_desc,
    )
    write_fasta(bounded, sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
