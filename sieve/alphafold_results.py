"""List per-model confidence scores for matching manifest sequence IDs."""
import csv
import json
import re
import sys
import zipfile
from decimal import Decimal

from sieve.alphafold_result_format import ResultRecord, write_result
from sieve.alphafold_confidence import model_confidence
from sieve.alphafold_files import model_number, read_seed, request_name

HIGHEST_CHOICES = ('ipTM', 'pTM', 'ipTM8_pTM2')

REQUIRED_COLUMNS = ('fold_locator', 'fold_id', 'sequence_id')


# Keep the public name for existing callers.
from sieve.alphafold_structure import open_structure_fold as open_result_fold, assembly_metadata


def selected_folds(manifest, sequence_id):
    if not re.fullmatch('[0-9a-f]{8}', sequence_id):
        raise ValueError('sequence ID must be eight lowercase hexadecimal characters')
    from sieve.alphafold_manifest import read_manifest_groups
    return [rows[0] for rows in read_manifest_groups(manifest).values()
            if rows[0]['sequence_id'] == sequence_id]


def confidence_rows(locator, error, *, include_seed=False):
    """Read only summaries; a malformed model does not suppress later models."""
    try:
        with open_result_fold(locator) as fold:
            models = sorted((model_number(name), name) for name in fold.models)
            if len({n for n, _ in models}) != len(models):
                raise ValueError('fold contains duplicate model numbers')
            seed = None
            metadata = assembly_metadata(fold)
            if metadata:
                for number, _ in models:
                    yield number, model_confidence({}), None
                return
            if include_seed:
                requests = {request_name(name) for _, name in models}
                if len(requests) != 1 or any(name in fold.member_names for name in requests):
                    seed = read_seed(fold.read, fold.models, fold.member_names)
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
                yield number, scores, seed
    except (OSError, ValueError, KeyError, RuntimeError, zipfile.BadZipFile) as exc:
        error(f'{locator}: {exc}')


def selection_score(scores, metric):
    fields = ('iptm', 'ptm') if metric == 'ipTM8_pTM2' else (metric.lower(),)
    if any(scores[field] is None for field in fields):
        return None
    values = [Decimal(str(scores[field])) for field in fields]
    if len(values) == 1:
        return values[0]
    return Decimal('0.8') * values[0] + Decimal('0.2') * values[1]


def run(manifest, sequence_id, stdout=None, stderr=None, *, highest=None, include_seed=False):
    stdout = sys.stdout if stdout is None else stdout
    stderr = sys.stderr if stderr is None else stderr
    failed = False

    def error(message):
        nonlocal failed
        failed = True
        print(f'ERROR: {message}', file=stderr)

    if highest is not None and highest not in HIGHEST_CHOICES:
        raise ValueError(f'invalid highest metric: {highest}')
    best_row = best_score = None
    folds = selected_folds(manifest, sequence_id)
    writer = csv.writer(stdout, delimiter='\t', lineterminator='\n')
    for fold in folds:
        locator = fold['fold_locator']
        for number, scores, seed in confidence_rows(locator, error, include_seed=include_seed):
            row = ResultRecord(locator, fold['fold_id'], number, scores['iptm'], scores['ptm'], seed)
            if highest is None:
                write_result(writer, row, include_seed=include_seed)
            else:
                score = selection_score(scores, highest)
                if score is not None and (best_score is None or score > best_score):
                    best_row, best_score = row, score
    if best_row is not None:
        write_result(writer, best_row, include_seed=include_seed)
    return int(failed)
