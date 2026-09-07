#!/usr/bin/env python3

import argparse
import sys

from Bio import SeqIO

from tangle import open_file_to_read, open_file_to_write


def positive_int(value):
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("length must be at least 1")
    return parsed


def filter_by_length(input_fasta, minimum, maximum, output_fasta):
    if minimum < 1 or maximum < 1:
        raise ValueError("lengths must be at least 1")
    if minimum > maximum:
        raise ValueError("minimum length cannot exceed maximum length")

    with (
        open_file_to_read(input_fasta) as input_stream,
        open_file_to_write(output_fasta, "wt") as output_stream,
    ):
        for record in SeqIO.parse(input_stream, "fasta"):
            sequence = str(record.seq)
            if minimum <= len(sequence) <= maximum:
                output_stream.write(f">{record.description}\n{sequence}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Filter FASTA sequences by inclusive length limits."
    )
    parser.add_argument("input_fasta")
    parser.add_argument("minimum", type=positive_int)
    parser.add_argument("maximum", type=positive_int)
    parser.add_argument("output_fasta")
    args = parser.parse_args(argv)

    if args.minimum > args.maximum:
        parser.error("minimum length cannot exceed maximum length")

    filter_by_length(args.input_fasta, args.minimum, args.maximum, args.output_fasta)
    return 0


if __name__ == "__main__":
    sys.exit(main())
