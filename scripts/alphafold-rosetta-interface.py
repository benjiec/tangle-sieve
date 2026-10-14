#!/usr/bin/env python3
"""Analyze a selected AlphaFold model using Dockerized Rosetta InterfaceAnalyzer."""
import argparse
import sys
import zipfile

from sieve.alphafold_rosetta_interface import run


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('locator', help='fold locator, or - for locators/list-results TSV on stdin')
    parser.add_argument('--interface', required=True, help='two chain IDs separated by a colon, e.g. A:B')
    parser.add_argument('--model', type=int, help='required for a plain locator')
    parser.add_argument('--output', help='preserve files in a new or empty directory; otherwise use temporary files')
    parser.add_argument('--tsv', action='store_true', help='print headered TSV, one row per model')
    parser.add_argument('--repack', action='store_true',
                        help='repack side chains in both complex and separated partners (no backbone relaxation)')
    args = parser.parse_args(argv)
    try:
        run(args.locator, sys.stdin, args.output, args.interface, args.model, args.tsv, sys.stdout, repack=args.repack)
    except (OSError, ValueError, KeyError, RuntimeError, zipfile.BadZipFile) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
