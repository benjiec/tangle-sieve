#!/usr/bin/env python3

import argparse
import sys

from sieve.alphafold_pdockq2 import collect_regions, score_zip, write_rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compute pairwise pDockQ2 for every model in an AlphaFold Server ZIP")
    parser.add_argument("zipfile")
    parser.add_argument("--region", action="append", default=[], metavar="CHAIN:START-END")
    parser.add_argument("--distance-cutoff", type=float, default=8.0)
    args = parser.parse_args(argv)
    try:
        write_rows(score_zip(args.zipfile, collect_regions(args.region), args.distance_cutoff), sys.stdout)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())
