#!/usr/bin/env python3

import argparse
import csv
import glob
import os
import sys

from sieve.alphafold_pdockq2 import (
    collect_regions,
    score_column_name,
    score_zip,
    scores_by_model_and_column,
)


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
        scores_by_zip = []
        full_columns = set()
        regional_columns = set()
        for path in paths:
            rows = score_zip(path, regions=regions, cutoff=args.distance_cutoff)
            scores = scores_by_model_and_column(rows)
            scores_by_zip.append((path, scores))
            for row in rows:
                target = full_columns if row["scope"] == "full" else regional_columns
                target.add(score_column_name(row))
        score_columns = sorted(full_columns) + sorted(regional_columns)
        columns = ["zipfile", "model", *score_columns]
        writer = csv.DictWriter(sys.stdout, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for path, scores in scores_by_zip:
            models = sorted({model for model, _column in scores})
            for model in models:
                output_row = {"zipfile": os.path.basename(path), "model": model}
                output_row.update({
                    column: f"{scores[(model, column)]:.6f}"
                    for column in score_columns
                    if (model, column) in scores
                })
                writer.writerow(output_row)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())
