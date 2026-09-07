#!/usr/bin/env python3

import argparse
import sys

from Bio import SeqIO

from tangle import open_file_to_read, open_file_to_write


def filter_accession_prefixes(big_fasta, filter_fasta, output_fasta):
    with open_file_to_read(filter_fasta) as filter_stream:
        filter_accessions = {
            record.id for record in SeqIO.parse(filter_stream, "fasta")
        }

    prefixes = {
        accession[:length]
        for accession in filter_accessions
        for length in range(1, len(accession) + 1)
    }

    with (
        open_file_to_read(big_fasta) as big_stream,
        open_file_to_write(output_fasta, "wt") as output_stream,
    ):
        for record in SeqIO.parse(big_stream, "fasta"):
            if record.id not in prefixes:
                continue
            output_stream.write(f">{record.description}\n{str(record.seq)}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Retain big-FASTA records whose accessions prefix accessions in a filter FASTA."
        )
    )
    parser.add_argument("big_fasta")
    parser.add_argument("filter_fasta")
    parser.add_argument("output_fasta")
    args = parser.parse_args(argv)

    filter_accession_prefixes(args.big_fasta, args.filter_fasta, args.output_fasta)
    return 0


if __name__ == "__main__":
    sys.exit(main())
