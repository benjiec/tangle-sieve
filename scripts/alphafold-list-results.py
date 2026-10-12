#!/usr/bin/env python3
"""Print headerless TSV: fold_locator, fold_id, model_id, ipTM, pTM."""
import argparse
import sys

from sieve.alphafold_results import HIGHEST_CHOICES, run


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', help='manifest TSV')
    parser.add_argument('sequence_id', help='exact sequence hash ID to select')
    parser.add_argument('--highest', nargs='?', choices=HIGHEST_CHOICES,
                        const='ipTM8_pTM2', default=None,
                        help='print only the highest-scoring row; without a value use 0.8*ipTM+0.2*pTM')
    args = parser.parse_args(argv)
    try:
        return run(args.manifest, args.sequence_id, highest=args.highest)
    except (OSError, ValueError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
