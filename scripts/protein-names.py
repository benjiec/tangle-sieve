#!/usr/bin/env python3
"""Write accession_id/protein_name TSV from a file of protein accessions.

Input contains one accession per line; blank lines are ignored. Names are FASTA
header descriptions with a trailing [organism] removed. Genome paths use the
current TANGLE_WORLD/TANGLE_AREA manifest. All lookups succeed before output.
"""

import argparse
import csv
import gzip
import re
import sys

from tangle.defaults import Defaults
from tangle.manifest import ManifestTable
from sieve.protein import CuratedProtein, _existing_path, _rows_from_table, _sql_string


def read_names(genome, accessions):
    path = Defaults.ncbi_genome_proteins(genome)
    if path is None:
        raise FileNotFoundError(f"Cannot find protein.faa for genome {genome}")
    path = _existing_path(path)
    names = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as source:
        for line in source:
            if not line.startswith(">"):
                continue
            parts = line[1:].strip().split(maxsplit=1)
            if not parts or parts[0] not in accessions:
                continue
            accession = parts[0]
            name = parts[1] if len(parts) == 2 else ""
            name = re.sub(r"\s*\[[^\[\]]*\]\s*$", "", name).strip()
            if not name:
                raise ValueError(f"Empty protein name for {accession} in {path}")
            if accession in names and names[accession] != name:
                raise ValueError(f"Conflicting protein names for {accession} in {path}")
            names[accession] = name
    missing = accessions - names.keys()
    if missing:
        raise ValueError(f"Cannot find proteins in {path}: {', '.join(sorted(missing))}")
    return names


def protein_names(accessions):
    if not accessions:
        return {}
    path = _existing_path(Defaults.area_sequence_manifest_tsv())
    genomes = {}
    for accession in dict.fromkeys(accessions):
        rows = _rows_from_table(ManifestTable, path, column_filters=[
            f"sequence_accession = {_sql_string(accession)}",
            "sequence_type = 'protein'",
        ])
        if not rows:
            raise ValueError(f"Cannot find protein {accession} in manifest")
        for row in rows:
            genome = row["sequence_database"]
            if not isinstance(genome, str) or not genome.strip():
                raise ValueError(f"Missing genome for {accession} in manifest")
            genomes.setdefault(genome, set()).add(accession)
    names = {}
    for genome in sorted(genomes):
        for accession, name in read_names(genome, genomes[genome]).items():
            if accession in names and names[accession] != name:
                raise ValueError(f"Conflicting protein names for {accession} across genomes")
            names[accession] = name
    return names


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("accessions_file", help="text file with one accession per line")
    args = parser.parse_args(argv)
    try:
        with open(args.accessions_file, encoding="utf-8") as source:
            accessions = [line.strip() for line in source if line.strip()]
        names = protein_names(accessions)
    except (OSError, ValueError) as error:
        print(f"{parser.prog}: {error}", file=sys.stderr)
        return 1
    finally:
        CuratedProtein.clear_cache()
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(("accession_id", "protein_name"))
    writer.writerows((accession, names[accession]) for accession in accessions)
    return 0


if __name__ == "__main__":
    sys.exit(main())
