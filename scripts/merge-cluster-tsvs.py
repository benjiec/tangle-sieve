#!/usr/bin/env python3

import argparse
from pathlib import Path
from tempfile import NamedTemporaryFile

from sieve.region_accession import parse_protein_accession


INPUT_HEADER = "Cluster Repr\tCluster Member"


def fasta_records(path):
    identifier = None
    parts = []
    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(parts), header_line
                fields = line[1:].split(maxsplit=1)
                if not fields:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
                identifier = fields[0]
                header_line = line_number
                parts = []
            else:
                if identifier is None or len(line.split()) != 1:
                    raise ValueError(f"{path}:{line_number}: invalid FASTA sequence line")
                parts.append(line)
    if identifier is not None:
        yield identifier, "".join(parts), header_line


def cluster_sequences(path):
    representative = None
    members = {}

    def pairs():
        if representative not in members:
            raise ValueError(
                f"{path}:{cluster_line}: missing sequence for representative {representative!r}"
            )
        for member, (sequence, line_number) in members.items():
            yield (representative, member), (members[representative][0], sequence), line_number

    for identifier, sequence, line_number in fasta_records(path):
        if not sequence:
            if representative is not None:
                yield from pairs()
            representative = identifier
            cluster_line = line_number
            members = {}
            continue
        if representative is None:
            raise ValueError(f"{path}:{line_number}: member sequence before a cluster header")
        if identifier in members and members[identifier][0] != sequence:
            raise ValueError(
                f"{path}:{line_number}: conflicting sequences for {representative!r}, {identifier!r}"
            )
        members[identifier] = (sequence, line_number)
    if representative is not None:
        yield from pairs()


def merge_cluster_tsvs(inputs, output, cluster_fastas=()):
    destination_exists = output.exists()
    for path in (*inputs, *cluster_fastas):
        if path.resolve() == output.resolve() or (destination_exists and path.samefile(output)):
            raise ValueError(f"output refers to input file {path}")

    sequences_by_pair = {}
    for path in cluster_fastas:
        for pair, sequences, line_number in cluster_sequences(path):
            if sequences_by_pair.setdefault(pair, sequences) != sequences:
                raise ValueError(f"{path}:{line_number}: conflicting sequences for cluster pair {pair!r}")

    temporary = None
    try:
        with NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="\n", dir=output.parent, delete=False,
        ) as stream:
            temporary = Path(stream.name)
            header = f"{INPUT_HEADER}\tProtein Accession"
            if cluster_fastas:
                header += "\tSame as Cluster Repr\tContained"
            stream.write(header + "\n")
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
                        result = f"{row}\t{protein}"
                        if cluster_fastas:
                            pair = tuple(fields)
                            if pair not in sequences_by_pair:
                                raise ValueError(
                                    f"{path}:{line_number}: missing FASTA sequences for cluster pair {pair!r}"
                                )
                            representative, member = sequences_by_pair[pair]
                            result += "\ttrue" if representative == member else "\tfalse"
                            contained = representative in member or member in representative
                            result += "\ttrue" if contained else "\tfalse"
                        stream.write(result + "\n")
        temporary.replace(output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Merge cluster TSVs and recover protein accessions.")
    parser.add_argument("--output", required=True, type=Path, help="output TSV")
    parser.add_argument("inputs", nargs="+", type=Path, help="input cluster TSVs")
    parser.add_argument(
        "--cluster-fasta", nargs="+", type=Path, default=(), metavar="FASTA",
        help="MMseqs all-sequences cluster FASTAs (place after the input TSVs)",
    )
    args = parser.parse_args(argv)
    try:
        merge_cluster_tsvs(args.inputs, args.output, args.cluster_fasta)
    except (OSError, ValueError) as error:
        parser.exit(1, f"{parser.prog}: {error}\n")


if __name__ == "__main__":
    main()
