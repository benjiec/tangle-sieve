#!/usr/bin/env python3

import argparse
import sys

from sieve.alphafold_pdockq2 import collect_regions, read_model_files, write_rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compute pairwise pDockQ2 for one AlphaFold Server model")
    parser.add_argument("cif")
    parser.add_argument("full_data_json")
    parser.add_argument("--region", action="append", default=[], metavar="CHAIN:START-END")
    parser.add_argument("--distance-cutoff", type=float, default=8.0)
    args = parser.parse_args(argv)
    try:
        rows = read_model_files(args.cif, args.full_data_json, collect_regions(args.region), args.distance_cutoff)
        write_rows(rows, sys.stdout)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())
