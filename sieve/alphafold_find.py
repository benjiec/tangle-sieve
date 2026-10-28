"""Find manifest folds by exact candidate headers and optional composition limits."""
import csv
import re
import sys
from collections import Counter

from sieve.alphafold_manifest import INPUT_COLUMNS
from sieve.alphafold_result_format import ResultRecord, write_result
from sieve.alphafold_results import HIGHEST_CHOICES, REQUIRED_COLUMNS, confidence_rows, selection_score


def read_names(path):
    # Only remove line terminators: spaces are part of the full FASTA header.
    with open(path, encoding='utf-8') as stream:
        return {name for line in stream if (name := line.rstrip('\r\n'))}


def selected_folds(manifest, names, *, cross=False, max_chains=None, max_unique_inputs=None):
    for label, value in [('max-chains', max_chains), ('max-unique-inputs', max_unique_inputs)]:
        if value is not None and (isinstance(value, bool) or not isinstance(value, int) or value < 1):
            raise ValueError(f'--{label} must be a positive integer')
    selected = []
    seen = {}
    with open(manifest, encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t', strict=True)
        columns = reader.fieldnames
        required = (*REQUIRED_COLUMNS, *INPUT_COLUMNS)
        if columns is None or len(set(columns)) != len(columns) or any(c not in columns for c in required):
            raise ValueError('manifest requires unique ID/locator columns and all five input description/copies pairs')
        for row in reader:
            line = reader.line_num
            if None in row or any(row.get(c) is None for c in columns) or any(not row[c] for c in REQUIRED_COLUMNS):
                raise ValueError(f'incomplete or malformed manifest row at line {line}')
            if any(not re.fullmatch('[0-9a-f]{8}', row[c]) for c in ('fold_id', 'sequence_id')):
                raise ValueError(f'invalid manifest ID at line {line}')
            composition = Counter()
            for i in range(1, 6):
                name, copies = row[f'input {i} description'], row[f'input {i} copies']
                if not name and not copies:
                    continue
                if not name or not re.fullmatch('[0-9]+', copies) or int(copies) < 1:
                    raise ValueError(f'invalid input {i} description/copies at line {line}')
                composition[name] += int(copies)
            if not composition:
                raise ValueError(f'fold has no inputs at line {line}')
            locator = row['fold_locator']
            signature = (row['fold_id'], row['sequence_id'], tuple(sorted(composition.items())))
            if locator in seen:
                if seen[locator] != signature:
                    raise ValueError(f'contradictory manifest rows for locator {locator} at line {line}')
                continue
            seen[locator] = signature
            inputs = set(composition)
            if not inputs.intersection(names) or (cross and not inputs.issubset(names)):
                continue
            if max_chains is not None and sum(composition.values()) > max_chains:
                continue
            if max_unique_inputs is not None and len(inputs) > max_unique_inputs:
                continue
            selected.append(row)
    return selected


def run(manifest, names_file, stdout=None, stderr=None, *, cross=False,
        max_chains=None, max_unique_inputs=None, highest=None):
    stdout = sys.stdout if stdout is None else stdout
    stderr = sys.stderr if stderr is None else stderr
    if highest is not None and highest not in HIGHEST_CHOICES:
        raise ValueError(f'invalid highest metric: {highest}')
    # Validate the entire manifest before opening any fold or writing results.
    folds = selected_folds(manifest, read_names(names_file), cross=cross,
                           max_chains=max_chains, max_unique_inputs=max_unique_inputs)
    failed = False
    def error(message):
        nonlocal failed
        failed = True
        print(f'ERROR: {message}', file=stderr)
    writer = csv.writer(stdout, delimiter='\t', lineterminator='\n')
    best = dict.fromkeys(fold['sequence_id'] for fold in folds)
    for fold in folds:
        for model, scores, _ in confidence_rows(fold['fold_locator'], error):
            record = ResultRecord(fold['fold_locator'], fold['fold_id'], model, scores['iptm'], scores['ptm'])
            if highest is None:
                write_result(writer, record)
            else:
                score = selection_score(scores, highest)
                previous = best[fold['sequence_id']]
                if score is not None and (previous is None or score > previous[0]):
                    best[fold['sequence_id']] = (score, record)
    if highest is not None:
        for winner in best.values():
            if winner is not None:
                write_result(writer, winner[1])
    return int(failed)
