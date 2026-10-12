#!/usr/bin/env python3
"""Report fixed-coordinate fold-integrity diagnostics with Dockerized PyRosetta."""
import argparse
import sys
import zipfile
from sieve.alphafold_rosetta_integrity import run


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('locator', help='fold locator or - for list-results/locator records on stdin')
    parser.add_argument('--model', type=int, help='required for a plain locator')
    parser.add_argument('--output', help='preserve files in a new or empty directory; otherwise use temporary files')
    parser.add_argument('--tsv', action='store_true', help='headered rows per model: complex plus each isolated chain')
    args = parser.parse_args(argv)
    try:
        run(args.locator, sys.stdin, args.output, args.model, args.tsv, sys.stdout)
    except (OSError, ValueError, KeyError, RuntimeError, zipfile.BadZipFile) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
