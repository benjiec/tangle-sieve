"""List per-model confidence scores for matching manifest sequence IDs."""
import csv
import json
import re
import sys
import zipfile
from contextlib import contextmanager
from pathlib import Path

from sieve.alphafold_confidence import model_confidence
from sieve.alphafold_files import Fold, model_number, open_fold

REQUIRED_COLUMNS = ('fold_locator', 'fold_id', 'sequence_id')


@contextmanager
def open_result_fold(locator):
    if '::' in locator:
        with open_fold(locator) as fold:
            yield fold
    else:
        path = Path(locator)
        if not path.is_dir():
            raise ValueError(f'expected an extracted fold directory or ZIP fold locator: {locator}')
        names = tuple(p.name for p in path.iterdir() if p.is_file())
        models = sorted(n for n in names if model_number(n) is not None)
        if not models:
            raise ValueError(f'no model CIF files at fold locator: {locator}')
        yield Fold(str(path), '.', models,
                   lambda name: (path / name).read_text(encoding='utf-8'), names)


def selected_folds(manifest, sequence_id):
    if not re.fullmatch('[0-9a-f]{8}', sequence_id):
        raise ValueError('sequence ID must be eight lowercase hexadecimal characters')
    with open(manifest, encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        columns = reader.fieldnames
        if (columns is None or len(set(columns)) != len(columns) or
                any(c not in columns for c in REQUIRED_COLUMNS)):
            raise ValueError('manifest requires unique fold_locator, fold_id, and sequence_id columns')
        selected = []
        for row in reader:
            if None in row or any(not row.get(c) for c in REQUIRED_COLUMNS):
                raise ValueError(f'incomplete or malformed manifest row at line {reader.line_num}')
            if any(not re.fullmatch('[0-9a-f]{8}', row[c]) for c in ('fold_id', 'sequence_id')):
                raise ValueError(f'invalid manifest ID at line {reader.line_num}')
            if row['sequence_id'] == sequence_id:
                selected.append(row)
    return selected


def confidence_rows(locator, error):
    """Read only summaries; a malformed model does not suppress later models."""
    try:
        with open_result_fold(locator) as fold:
            models = sorted((model_number(name), name) for name in fold.models)
            if len({n for n, _ in models}) != len(models):
                raise ValueError('fold contains duplicate model numbers')
            for number, name in models:
                summary = name.removesuffix(f'model_{number}.cif') + f'summary_confidences_{number}.json'
                try:
                    count = fold.member_names.count(summary)
                    if count > 1:
                        raise ValueError(f'duplicate summary file: {summary}')
                    scores = model_confidence(json.loads(fold.read(summary)) if count else {})
                except (OSError, ValueError, KeyError, RuntimeError, zipfile.BadZipFile) as exc:
                    error(f'{locator} model {number}: {exc}')
                    continue
                yield number, scores
    except (OSError, ValueError, KeyError, RuntimeError, zipfile.BadZipFile) as exc:
        error(f'{locator}: {exc}')


def run(manifest, sequence_id, stdout=None, stderr=None):
    stdout = sys.stdout if stdout is None else stdout
    stderr = sys.stderr if stderr is None else stderr
    failed = False

    def error(message):
        nonlocal failed
        failed = True
        print(f'ERROR: {message}', file=stderr)

    folds = selected_folds(manifest, sequence_id)
    writer = csv.writer(stdout, delimiter='\t', lineterminator='\n')
    for fold in folds:
        locator = fold['fold_locator']
        for number, scores in confidence_rows(locator, error):
            writer.writerow([locator, fold['fold_id'], number,
                             *('NA' if scores[field] is None else str(scores[field])
                               for field in ('iptm', 'ptm'))])
    return int(failed)
