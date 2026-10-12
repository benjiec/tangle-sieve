#!/usr/bin/env python3
"""Print headerless TSV: fold_locator, fold_id, model_id, ipTM, pTM."""
import argparse
import sys

from sieve.alphafold_results import run


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', help='manifest TSV')
    parser.add_argument('sequence_id', help='exact sequence hash ID to select')
    args = parser.parse_args(argv)
    try:
        return run(args.manifest, args.sequence_id)
    except (OSError, ValueError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
