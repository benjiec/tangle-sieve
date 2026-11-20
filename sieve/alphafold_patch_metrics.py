"""Patch extraction shared by the human report and interface evaluation."""
from dataclasses import dataclass
from pathlib import Path
import subprocess
import zipfile

from sieve.alphafold_files import discover, fold_locator, model_number
from sieve.alphafold_pdockq2_patches import analyze
from sieve.alphafold_structure import open_structure_fold
from sieve.interface_description import DescriptionOptions
from sieve.patch_description_summary import PatchDescriptions


ERRORS = (OSError, ValueError, KeyError, TypeError, RuntimeError, zipfile.BadZipFile, subprocess.SubprocessError)


def jobs_for(request, error):
    path = Path(request.locator)
    if '::' in request.locator:
        source, member = request.locator.split('::', 1)
        locators = [f'{Path(source).resolve()}::{member}']
    elif path.suffix.lower() == '.cif':
        locators = [str(path.resolve())]
    elif path.is_dir() and any(model_number(p.name) is not None for p in path.iterdir() if p.is_file()):
        locators = [str(path.resolve())]
    else:
        locators = (fold_locator(f.source_path, f.fold_path) for f in discover(path, error))
    found = False
    for locator in locators:
        found = True
        try:
            with open_structure_fold(locator) as fold:
                models = [model_number(name) for name in fold.models]
            if len(models) != len(set(models)):
                raise ValueError('duplicate model numbers')
            if request.model_id is not None:
                if request.model_id not in models:
                    raise ValueError(f'model {request.model_id} not found')
                models = [request.model_id]
            for model in sorted(models, key=lambda number: -1 if number is None else number):
                yield locator, model
        except ERRORS as exc:
            error(f'{locator}: {exc}')
    if not found:
        error(f'{request.locator}: no folds found')


@dataclass
class PatchMetrics:
    rows: list
    descriptions: PatchDescriptions
    targets: dict
    model_scores: dict
    fold_ids: dict
    sequence_ids: dict

    def leaves(self):
        for row in self.rows:
            if row['scope'] == 'patch' and row['status'] == 'leaf':
                yield {
                    **row,
                    **self.descriptions.describe(row)['metrics'],
                    'sequence_id': self.sequence_ids[row['source']],
                    'fold_id': self.fold_ids[row['source']],
                    **self.model_scores.get((row['source'], row['model']), {}),
                }


def extract_patch_metrics(locator, *, options, chains=None, full_data=None, dssp=None,
                          description_options=None, selected_model=None, model_records=None,
                          model_seeds=None, model_scores=None):
    targets, inputs, fold_ids, sequence_ids = {}, {}, {}, {}
    model_scores = {} if model_scores is None else model_scores
    rows = analyze(
        locator, full_data, options, chains, summary_targets=targets, summary_inputs=inputs,
        model_scores=model_scores, fold_ids=fold_ids, sequence_ids=sequence_ids,
        selected_model=selected_model, model_records=model_records, model_seeds=model_seeds,
    )
    descriptions = PatchDescriptions(
        inputs, dssp, description_options or DescriptionOptions(cutoff=options.cutoff),
    )
    return PatchMetrics(rows, descriptions, targets, model_scores, fold_ids, sequence_ids)
