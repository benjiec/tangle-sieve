#!/usr/bin/env python3
"""Find folds by full FASTA headers; print headerless list-results TSV."""
import argparse
import csv
import sys

from sieve.alphafold_find import run
from sieve.alphafold_results import HIGHEST_CHOICES


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', help='manifest TSV')
    parser.add_argument('sequence_names', help='exact full FASTA headers, one per line, without >')
    parser.add_argument('--cross', action='store_true', help='require every input name to be listed')
    parser.add_argument('--max-chains', type=int, help='maximum total chain copies (inclusive)')
    parser.add_argument('--max-unique-inputs', type=int, help='maximum distinct input names (inclusive)')
    parser.add_argument('--highest', nargs='?', choices=HIGHEST_CHOICES, const='ipTM8_pTM2',
                        help='one best model per sequence ID; default metric 0.8*ipTM+0.2*pTM')
    args = parser.parse_args(argv)
    try:
        return run(args.manifest, args.sequence_names, cross=args.cross, max_chains=args.max_chains,
                   max_unique_inputs=args.max_unique_inputs, highest=args.highest)
    except (OSError, ValueError, csv.Error) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
