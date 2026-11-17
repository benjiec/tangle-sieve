#!/usr/bin/env python3

import argparse
from pathlib import Path
from tempfile import NamedTemporaryFile

from sieve.region_accession import parse_protein_accession


INPUT_HEADER = "Cluster Repr\tCluster Member"


def merge_cluster_tsvs(inputs, output):
    destination_exists = output.exists()
    for path in inputs:
        if path.resolve() == output.resolve() or (destination_exists and path.samefile(output)):
            raise ValueError(f"output refers to input file {path}")

    temporary = None
    try:
        with NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="\n", dir=output.parent, delete=False,
        ) as stream:
            temporary = Path(stream.name)
            stream.write(f"{INPUT_HEADER}\tProtein Accession\n")
            for path in inputs:
                with path.open(encoding="utf-8") as source:
                    if source.readline().rstrip("\n") != INPUT_HEADER:
                        raise ValueError(f"{path}:1: expected header {INPUT_HEADER!r}")
                    for line_number, line in enumerate(source, 2):
                        row = line.rstrip("\n")
                        fields = row.split("\t")
                        if len(fields) != 2:
                            raise ValueError(f"{path}:{line_number}: expected two tab-separated columns")
                        try:
                            protein = parse_protein_accession(fields[1])
                        except ValueError as error:
                            raise ValueError(f"{path}:{line_number}: {error}") from error
                        stream.write(f"{row}\t{protein}\n")
        temporary.replace(output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Merge cluster TSVs and recover protein accessions.")
    parser.add_argument("--output", required=True, type=Path, help="output TSV")
    parser.add_argument("inputs", nargs="+", type=Path, help="input cluster TSVs")
    args = parser.parse_args(argv)
    try:
        merge_cluster_tsvs(args.inputs, args.output)
    except (OSError, ValueError) as error:
        parser.exit(1, f"{parser.prog}: {error}\n")


if __name__ == "__main__":
    main()
