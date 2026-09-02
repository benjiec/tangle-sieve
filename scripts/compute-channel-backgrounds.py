#!/usr/bin/env python3

import argparse
import sys

from sieve.channel_set import load_channel_set


def _print_record(sequence_id):
    print(f"Processing {sequence_id}", file=sys.stderr)


def _print_invalid_record(sequence_id, error):
    print(f"Skipping {sequence_id}: {error}", file=sys.stderr)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta")
    parser.add_argument("-c", "--channel-set", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args(argv)

    channel_set = load_channel_set(args.channel_set)
    backgrounds = channel_set.compute_backgrounds(
        args.fasta,
        on_record=_print_record,
        on_invalid_sequence=_print_invalid_record,
    )
    channel_set.save_backgrounds(backgrounds, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
