"""FASTA identifiers shared by the region extractor and cluster TSV merger."""

import re


PFAM_PATTERN = r"PF[0-9]{5}(?:\.[1-9][0-9]*)?"
POSITIVE_INTEGER = r"[1-9][0-9]*"
REGION_ACCESSION = re.compile(
    rf"(?P<protein>\S+)_(?P<pfam>{PFAM_PATTERN})"
    rf"_(?P<minimum>{POSITIVE_INTEGER})_(?P<maximum>{POSITIVE_INTEGER})"
    rf"(?:_domain_(?P<domain>{POSITIVE_INTEGER}))?"
)


def format_region_accession(protein, pfam, minimum, maximum, domain=None):
    if not re.fullmatch(PFAM_PATTERN, pfam):
        raise ValueError(f"invalid Pfam accession {pfam!r}")
    accession = f"{protein}_{pfam}_{minimum}_{maximum}"
    if domain is not None:
        accession += f"_domain_{domain}"
    parse_protein_accession(accession)
    return accession


def parse_protein_accession(accession):
    match = REGION_ACCESSION.fullmatch(accession)
    if match is None:
        raise ValueError(f"invalid extracted region accession {accession!r}")
    minimum = int(match["minimum"])
    maximum = int(match["maximum"])
    if minimum > maximum or (match["domain"] is not None and int(match["domain"]) > maximum):
        raise ValueError(f"invalid match counts or domain index in {accession!r}")
    return match["protein"]
