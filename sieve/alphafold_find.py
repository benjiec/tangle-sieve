"""Find manifest folds by exact candidate headers and optional composition limits."""
import csv
import sys

from tangle import open_file_to_read

from sieve.alphafold_manifest import read_manifest_groups
from sieve.alphafold_result_format import ResultRecord, write_result
from sieve.alphafold_results import HIGHEST_CHOICES, confidence_rows, selection_score


class MissingNames(ValueError):
    def __init__(self, names):
        self.names = sorted(names)
        super().__init__('sequence names not found in manifest: ' + ', '.join(repr(n) for n in self.names))


def read_names(path):
    # Only remove line terminators: spaces are part of the full FASTA header.
    with open_file_to_read(path) as stream:
        return {name for line in stream if (name := line.rstrip('\r\n'))}


def selected_folds(manifest, names, *, cross=False, max_chains=None, max_unique_inputs=None):
    for label, value in [('max-chains', max_chains), ('max-unique-inputs', max_unique_inputs)]:
        if value is not None and (isinstance(value, bool) or not isinstance(value, int) or value < 1):
            raise ValueError(f'--{label} must be a positive integer')
    groups = read_manifest_groups(manifest)
    available = {row['input description'] for rows in groups.values() for row in rows}
    missing = set(names) - available
    if missing:
        raise MissingNames(missing)
    selected = []
    for rows in groups.values():
        inputs = {row['input description'] for row in rows}
        if not inputs.intersection(names) or (cross and not inputs.issubset(names)):
            continue
        if max_chains is not None and len(rows) > max_chains:
            continue
        if max_unique_inputs is not None and len(inputs) > max_unique_inputs:
            continue
        selected.append(rows[0])
    return selected


def run(manifest, names_file, stdout=None, stderr=None, *, cross=False,
        max_chains=None, max_unique_inputs=None, highest=None):
    stdout = sys.stdout if stdout is None else stdout
    stderr = sys.stderr if stderr is None else stderr
    if highest is not None and highest not in HIGHEST_CHOICES:
        raise ValueError(f'invalid highest metric: {highest}')
    # Validate the entire manifest before opening any fold or writing results.
    try:
        folds = selected_folds(manifest, read_names(names_file), cross=cross,
                               max_chains=max_chains, max_unique_inputs=max_unique_inputs)
    except MissingNames as exc:
        for name in exc.names:
            print(f'ERROR: sequence name not found in manifest: {name!r}', file=stderr)
        return 1
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
