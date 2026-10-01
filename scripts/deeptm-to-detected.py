#!/usr/bin/env python3
"""Convert DeepTMHMM three-line records to DetectedTable TSV."""

import argparse
import sys
from itertools import groupby

from tangle import open_file_to_read, unique_batch
from tangle.detected import DetectedTable


REGIONS = {"S": "signal", "O": "outside", "M": "membrane", "I": "inside"}


def parse_records(stream):
    """Read header, sequence, annotation triples; allow blank record separators."""
    lines = enumerate(stream, 1)
    for line_number, raw_header in lines:
        header = raw_header.strip()
        if not header:
            continue
        if not header.startswith(">"):
            raise ValueError(f"Line {line_number}: expected a > header")
        accession_text = header[1:].split("|", 1)[0].strip()
        if not accession_text:
            raise ValueError(f"Line {line_number}: header has no accession")
        accession = accession_text.split()[0]
        fields = []
        for field in ("sequence", "annotation"):
            entry = next(lines, None)
            if entry is None:
                raise ValueError(f"{accession}: missing {field} line")
            number, raw = entry
            value = raw.strip()
            if not value or value.startswith(">"):
                raise ValueError(f"Line {number}: {accession} has missing {field}")
            if any(char.isspace() for char in value):
                raise ValueError(f"Line {number}: whitespace within {field}")
            fields.append(value)
        sequence, annotation = fields
        if len(sequence) != len(annotation):
            raise ValueError(
                f"{accession}: sequence length {len(sequence)} differs from "
                f"annotation length {len(annotation)}"
            )
        unknown = set(annotation) - REGIONS.keys()
        if unknown:
            raise ValueError(f"{accession}: unknown annotation labels {sorted(unknown)}")
        yield accession, sequence, annotation


def detected_rows(records, batch):
    """Emit contiguous regions with 1-based, inclusive query coordinates."""
    for accession, _sequence, annotation in records:
        start = 1
        for label, region in groupby(annotation):
            length = sum(1 for _ in region)
            end = start + length - 1
            yield dict(
                detection_type="sequence",
                detection_method="other",
                batch=batch,
                query_accession=accession,
                query_database="_",
                target_accession=REGIONS[label],
                target_database="DeepTMHMM",
                query_start=start,
                query_end=end,
                target_start=1,
                target_end=length,
            )
            start = end + 1


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="DeepTMHMM header/sequence/annotation triples")
    parser.add_argument("output", help="DetectedTable TSV; coordinates are 1-based inclusive")
    args = parser.parse_args(argv)
    with open_file_to_read(args.input) as stream:
        # Validate all records before opening the output file.
        rows = list(detected_rows(parse_records(stream), unique_batch()))
    DetectedTable.write_tsv(args.output, rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
