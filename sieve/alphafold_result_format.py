"""Shared headerless TSV contract for AlphaFold result producers and consumers."""
import csv
import re
from dataclasses import dataclass, fields

from sieve.alphafold_confidence import model_confidence


@dataclass(frozen=True)
class ResultRecord:
    fold_locator: str
    fold_id: str
    model_id: int
    ipTM: object
    pTM: object


RESULT_COLUMNS = tuple(field.name for field in fields(ResultRecord))


def write_result(writer, record):
    writer.writerow(['NA' if getattr(record, name) is None else getattr(record, name)
                     for name in RESULT_COLUMNS])


def parse_result(values):
    if len(values) != len(RESULT_COLUMNS):
        raise ValueError(f'expected {len(RESULT_COLUMNS)} result columns')
    data = dict(zip(RESULT_COLUMNS, values))
    if not data['fold_locator'] or not re.fullmatch('[0-9a-f]{8}', data['fold_id']):
        raise ValueError('invalid fold locator or fold ID')
    if not re.fullmatch('[0-9]+', data['model_id']):
        raise ValueError('model_id must be a nonnegative integer')
    scores = model_confidence({name.lower(): None if data[name] == 'NA' else float(data[name])
                               for name in ('ipTM', 'pTM')})
    return ResultRecord(data['fold_locator'], data['fold_id'], int(data['model_id']),
                        scores['iptm'], scores['ptm'])


@dataclass(frozen=True)
class PatchInput:
    locator: str
    model_id: object = None
    record_number: object = None


def read_patch_inputs(stream):
    reader = csv.reader(stream, delimiter='\t', strict=True)
    number = 0
    try:
        for number, values in enumerate(reader, 1):
            if not values or len(values) == 1 and not values[0].strip():
                continue
            if len(values) == 1:
                yield PatchInput(values[0], record_number=number)
            else:
                try:
                    record = parse_result(values)
                except ValueError as exc:
                    raise ValueError(f'stdin record {number}: {exc}') from exc
                yield PatchInput(record.fold_locator, record.model_id, number)
    except csv.Error as exc:
        raise ValueError(f'stdin record near line {reader.line_num}: {exc}') from exc
