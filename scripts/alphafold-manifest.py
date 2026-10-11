#!/usr/bin/env python3
"""Recursively identify AlphaFold chains and checkpoint successful folds to TSV."""
import argparse
import sys
from sieve.alphafold_manifest import run


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='directory or ZIP archive to search recursively')
    parser.add_argument('candidates', help='master FASTA; full headers identify sequences')
    parser.add_argument('manifest', help='output/cache TSV; updated after each successful fold')
    args = parser.parse_args(argv)
    try:
        return run(args.input, args.candidates, args.manifest)
    except (OSError, ValueError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
