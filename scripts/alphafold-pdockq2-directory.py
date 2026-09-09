#!/usr/bin/env python3

import argparse
import csv
import glob
import os
import sys

from sieve.alphafold_pdockq2 import full_scores_by_model_and_pair, score_zip


def main(argv=None):
    parser = argparse.ArgumentParser(description="Report the best pairwise pDockQ2 in each AlphaFold Server ZIP")
    parser.add_argument("directory")
    parser.add_argument("--distance-cutoff", type=float, default=8.0)
    args = parser.parse_args(argv)
    try:
        paths = sorted(glob.glob(os.path.join(args.directory, "*.zip")))
        if not paths:
            raise ValueError(f"directory contains no ZIP files: {args.directory}")
        scores_by_zip = []
        all_pairs = set()
        for path in paths:
            scores = full_scores_by_model_and_pair(score_zip(path, cutoff=args.distance_cutoff))
            scores_by_zip.append((path, scores))
            all_pairs.update((chain1, chain2) for _model, chain1, chain2 in scores)
        pairs = sorted(all_pairs)
        columns = ["zipfile", "model", *(f"{chain1}_{chain2}_pDockQ2_max" for chain1, chain2 in pairs)]
        writer = csv.DictWriter(sys.stdout, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for path, scores in scores_by_zip:
            models = sorted({model for model, _chain1, _chain2 in scores})
            for model in models:
                output_row = {"zipfile": os.path.basename(path), "model": model}
                output_row.update({
                    f"{chain1}_{chain2}_pDockQ2_max": f"{scores[(model, chain1, chain2)]:.6f}"
                    for chain1, chain2 in pairs
                    if (model, chain1, chain2) in scores
                })
                writer.writerow(output_row)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())
