#!/usr/bin/env python3

import argparse
import csv
import glob
import os
import sys

from sieve.alphafold_pdockq2 import best_rows_by_pair, collect_regions, score_zip


COLUMNS = [
    "zipfile", "chain 1", "chain 2",
    "best full model", "best full pDockQ2 1 to 2", "best full pDockQ2 2 to 1",
    "best full pDockQ2 max",
    "best regional model", "best regional pDockQ2 1 to 2",
    "best regional pDockQ2 2 to 1", "best regional pDockQ2 max",
    "best regional region 1", "best regional region 2",
]


def main(argv=None):
    parser = argparse.ArgumentParser(description="Report the best pairwise pDockQ2 in each AlphaFold Server ZIP")
    parser.add_argument("directory")
    parser.add_argument("--region", action="append", default=[], metavar="CHAIN:START-END")
    parser.add_argument("--distance-cutoff", type=float, default=8.0)
    args = parser.parse_args(argv)
    try:
        regions = collect_regions(args.region)
        paths = sorted(glob.glob(os.path.join(args.directory, "*.zip")))
        if not paths:
            raise ValueError(f"directory contains no ZIP files: {args.directory}")
        writer = csv.DictWriter(sys.stdout, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for path in paths:
            rows = score_zip(path, regions, args.distance_cutoff)
            for (chain1, chain2), full, regional in best_rows_by_pair(rows):
                output_row = {
                    "zipfile": os.path.basename(path),
                    "chain 1": chain1,
                    "chain 2": chain2,
                    "best full model": full["model"],
                    "best full pDockQ2 1 to 2": f'{full["pDockQ2 1 to 2"]:.6f}',
                    "best full pDockQ2 2 to 1": f'{full["pDockQ2 2 to 1"]:.6f}',
                    "best full pDockQ2 max": f'{full["pDockQ2 max"]:.6f}',
                    "best regional model": "",
                    "best regional pDockQ2 1 to 2": "",
                    "best regional pDockQ2 2 to 1": "",
                    "best regional pDockQ2 max": "",
                    "best regional region 1": "",
                    "best regional region 2": "",
                }
                if regional is not None:
                    output_row.update({
                        "best regional model": regional["model"],
                        "best regional pDockQ2 1 to 2": f'{regional["pDockQ2 1 to 2"]:.6f}',
                        "best regional pDockQ2 2 to 1": f'{regional["pDockQ2 2 to 1"]:.6f}',
                        "best regional pDockQ2 max": f'{regional["pDockQ2 max"]:.6f}',
                        "best regional region 1": regional["region 1"],
                        "best regional region 2": regional["region 2"],
                    })
                writer.writerow(output_row)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())
