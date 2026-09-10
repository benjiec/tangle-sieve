#!/usr/bin/env python3

import argparse
import csv
import glob
import hashlib
import os
import sys
import tempfile

from sieve.alphafold_pdockq2 import (
    collect_regions,
    score_column_name,
    score_zip,
    scores_by_model_and_column,
)


FIXED_COLUMNS = ["zipfile", "zipfile_sha256", "model"]


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_cache(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None, []
    with open(path, encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        columns = reader.fieldnames
        if columns is None or any(column not in columns for column in FIXED_COLUMNS):
            raise ValueError(f"cache TSV has an invalid header: {path}")
        rows = list(reader)
    if any(None in row for row in rows):
        raise ValueError(f"cache TSV contains rows wider than its header: {path}")
    return columns, rows


def write_cache(path, columns, rows):
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            "w", encoding="utf-8", newline="", dir=directory, delete=False
        ) as stream:
            temporary_path = stream.name
            writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_path, path)
        temporary_path = None
    finally:
        if temporary_path is not None and os.path.exists(temporary_path):
            os.remove(temporary_path)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Report the best pairwise pDockQ2 in each AlphaFold Server ZIP")
    parser.add_argument("directory")
    parser.add_argument("--output", required=True, help="TSV cache and output path")
    parser.add_argument("--region", action="append", default=[], metavar="CHAIN:START-END")
    parser.add_argument("--distance-cutoff", type=float, default=8.0)
    args = parser.parse_args(argv)
    try:
        regions = collect_regions(args.region)
        paths = sorted(glob.glob(os.path.join(args.directory, "*.zip")))
        if not paths:
            raise ValueError(f"directory contains no ZIP files: {args.directory}")
        existing_columns, existing_rows = read_cache(args.output)
        cached_hashes = {row["zipfile_sha256"] for row in existing_rows}
        scores_by_zip = []
        full_columns = set()
        regional_columns = set()
        for path in paths:
            digest = sha256_file(path)
            if digest in cached_hashes:
                print(f"Skipping {os.path.basename(path)}: cached", file=sys.stderr)
                continue
            print(f"Processing {os.path.basename(path)}", file=sys.stderr)
            rows = score_zip(path, regions=regions, cutoff=args.distance_cutoff)
            scores = scores_by_model_and_column(rows)
            scores_by_zip.append((path, digest, scores))
            for row in rows:
                target = full_columns if row["scope"] == "full" else regional_columns
                target.add(score_column_name(row))
        new_score_columns = sorted(full_columns) + sorted(regional_columns)
        columns = FIXED_COLUMNS + new_score_columns
        if existing_columns is not None:
            if scores_by_zip and existing_columns != columns:
                raise ValueError(
                    "existing TSV header does not match columns required by this invocation"
                )
            columns = existing_columns
        score_columns = columns[len(FIXED_COLUMNS):]
        new_rows = []
        for path, digest, scores in scores_by_zip:
            models = sorted({model for model, _column in scores})
            for model in models:
                output_row = {
                    "zipfile": os.path.basename(path),
                    "zipfile_sha256": digest,
                    "model": model,
                }
                output_row.update({
                    column: f"{scores[(model, column)]:.6f}"
                    for column in score_columns
                    if (model, column) in scores
                })
                new_rows.append(output_row)
            print(f"Processed {os.path.basename(path)}: {len(models)} models", file=sys.stderr)
        if new_rows or existing_columns is None:
            write_cache(args.output, columns, existing_rows + new_rows)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())
