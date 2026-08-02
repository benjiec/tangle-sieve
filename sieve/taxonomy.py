import csv
import os

from tangle.defaults import Defaults


TAXONOMY_FIELDS = {
    "domain",
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
}


def normalized(value):
    return str(value).strip().lower()


def genome_accession(row):
    for key in row:
        if normalized(key).replace("_", " ") in {"genome accession", "accession"}:
            return row[key]
    return ""


def read_taxonomy_rows():
    taxonomy_tsv = Defaults.area_genome_taxon_tsv()
    if not os.path.exists(taxonomy_tsv):
        return {}
    with open(taxonomy_tsv, "r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return {
        genome_accession(row): row
        for row in rows
        if genome_accession(row)
    }


def taxonomy_matches(row, taxon):
    expected = normalized(taxon)
    for key, value in row.items():
        if normalized(key) in TAXONOMY_FIELDS and normalized(value) == expected:
            return True
    return False
