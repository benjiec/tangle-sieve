#!/usr/bin/env python3

import argparse
import sys

from Bio import SeqIO

from tangle import open_file_to_read, open_file_to_write


def write_unique_sequences(input_fasta, output_fasta):
    seen_sequences = set()
    with (
        open_file_to_read(input_fasta) as input_f,
        open_file_to_write(output_fasta, "wt") as output_f,
    ):
        for record in SeqIO.parse(input_f, "fasta"):
            sequence = str(record.seq)
            if sequence in seen_sequences:
                continue
            seen_sequences.add(sequence)
            output_f.write(f">{record.id}\n{sequence}\n")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta")
    parser.add_argument("output_fasta")
    args = parser.parse_args(argv)

    write_unique_sequences(args.input_fasta, args.output_fasta)
    return 0


if __name__ == "__main__":
    sys.exit(main())
