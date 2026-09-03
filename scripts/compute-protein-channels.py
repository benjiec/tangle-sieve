#!/usr/bin/env python3

import argparse
import csv
import os
import shutil
import sys
import tempfile

from sieve.channel_set import load_channel_set
from sieve.channels import iter_fasta_records


def write_channels_tsv(
    fasta_path,
    channel_set,
    backgrounds,
    output_path,
    on_record=None,
    on_invalid_sequence=None,
    on_duplicate_sequence_id=None,
):
    functions = channel_set.make_functions(backgrounds)
    short_names = [channel.short_name for channel in channel_set]
    header = ["filename", "sequence_id", "position", "residue", *short_names]
    filename = os.path.basename(os.fspath(fasta_path))
    output_directory = os.path.dirname(os.path.abspath(output_path))
    append = os.path.exists(output_path) and os.path.getsize(output_path) > 0
    if append:
        with open(output_path, encoding="utf-8", newline="") as existing_output:
            existing_header = next(
                csv.reader(existing_output, delimiter="\t"),
                None,
            )
        if existing_header != header:
            raise ValueError(
                "existing output TSV header does not match the configured ChannelSet"
            )
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            newline="",
            dir=output_directory,
            delete=False,
        ) as output:
            temporary_path = output.name
            writer = csv.writer(output, delimiter="\t", lineterminator="\n")
            if not append:
                writer.writerow(header)
            for sequence_id, sequence in iter_fasta_records(
                fasta_path,
                on_invalid_sequence=on_invalid_sequence,
                on_duplicate_sequence_id=on_duplicate_sequence_id,
            ):
                if on_record is not None:
                    on_record(sequence_id)
                values_by_channel = {
                    short_name: function(sequence)
                    for short_name, function in functions.items()
                }
                for short_name, values in values_by_channel.items():
                    if len(values) != len(sequence):
                        raise ValueError(
                            f"channel {short_name!r} produced {len(values)} values "
                            f"for a sequence of length {len(sequence)}"
                        )
                for index, residue in enumerate(sequence):
                    writer.writerow([
                        filename,
                        sequence_id,
                        index + 1,
                        residue,
                        *(values_by_channel[name][index] for name in short_names),
                    ])
        if append:
            with open(output_path, "rb") as existing_output:
                existing_output.seek(-1, os.SEEK_END)
                needs_newline = existing_output.read(1) not in (b"\n", b"\r")
            with open(output_path, "a", encoding="utf-8", newline="") as output:
                if needs_newline:
                    output.write("\n")
                with open(temporary_path, encoding="utf-8", newline="") as batch:
                    shutil.copyfileobj(batch, output)
            os.remove(temporary_path)
        else:
            os.replace(temporary_path, output_path)
        temporary_path = None
    finally:
        if temporary_path is not None and os.path.exists(temporary_path):
            os.remove(temporary_path)


def _print_record(sequence_id):
    print(f"Processing {sequence_id}", file=sys.stderr)


def _print_invalid_record(sequence_id, error):
    print(f"Skipping {sequence_id}: {error}", file=sys.stderr)


def _print_duplicate_record(sequence_id):
    print(f"Skipping duplicate accession {sequence_id}", file=sys.stderr)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta")
    parser.add_argument("-c", "--channel-set", required=True)
    parser.add_argument("-b", "--backgrounds", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args(argv)

    channel_set = load_channel_set(args.channel_set)
    backgrounds = channel_set.load_backgrounds(args.backgrounds)
    write_channels_tsv(
        args.fasta,
        channel_set,
        backgrounds,
        args.output,
        on_record=_print_record,
        on_invalid_sequence=_print_invalid_record,
        on_duplicate_sequence_id=_print_duplicate_record,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
