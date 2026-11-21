"""Evaluate AlphaFold patches using registered interface profiles."""
from itertools import combinations
from pathlib import Path
import sys

from sieve.alphafold_cache import write_cache
from sieve.alphafold_death_rules import read_domains
from sieve.alphafold_manifest import read_manifest_groups
from sieve.alphafold_patch_metrics import ERRORS, jobs_for
from sieve.alphafold_pdockq2 import parse_model
from sieve.alphafold_pdockq2_patches import PATCH_COLUMNS, input_models, patch_export_rows
from sieve.alphafold_secondary import find_dssp
from sieve.interface.registry import INTERFACES, resolve

SIDE_FIELDS = {
    'sequence name': 'sequence name {side}',
    'domain': 'domain {side}',
    'candidate face': 'candidate face {side}',
    'face status': 'face status {side}',
}
PROFILE_COLUMNS = tuple(dict.fromkeys(column for profile in INTERFACES for column in profile.columns))
EXTRA_COLUMNS = ('iptm', 'ptm', 'interface_family', *PROFILE_COLUMNS)
RESULT_COLUMNS = ('confidence', 'chemical feasibility', 'interface')
COLUMNS = (*PATCH_COLUMNS, *SIDE_FIELDS, *EXTRA_COLUMNS, *RESULT_COLUMNS)


def domain_overrides(values):
    assignments = {}
    for value in values:
        chain, separator, domain = value.partition(':')
        if not separator or not chain.strip() or not domain.strip() or ':' in domain:
            raise ValueError(f'expected --domain CHAIN:DOMAIN, got {value!r}')
        if chain in assignments:
            raise ValueError(f'duplicate domain assignment for chain {chain}')
        assignments[chain] = domain
    return assignments


def manifest_rows(groups, locator):
    path = Path(locator)
    if '::' not in locator and path.suffix.lower() == '.cif':
        return groups.get(locator, groups.get(str(path.parent.resolve()), ()))
    # A root ZIP can also be addressed by its explicit archive locator.
    return groups.get(locator, groups.get(locator.removesuffix('::.'), ()))


def evaluate_model(locator, number, request, *, groups, domain_table, overrides,
                   chains, full_data, dssp, error):
    seeds, model_scores = {}, {}
    record, = input_models(
        locator, full_data, selected_model=number, model_seeds=seeds, model_scores=model_scores,
    )
    residues, _ = parse_model(record[2], record[3])
    if chains and any(chain not in residues for chain in chains):
        raise ValueError(f'--chains requires two chains present in the model: {list(residues)}')
    if len(residues) < 2:
        raise ValueError('an interface requires at least two chains')
    if set(overrides) - residues.keys():
        unknown = ', '.join(sorted(set(overrides) - residues.keys()))
        raise ValueError(f'domain assignments reference chains absent from the model: {unknown}')
    group = manifest_rows(groups, locator)
    if group and {row['input chain'] for row in group} != residues.keys():
        raise ValueError('model/manifest chain mismatch')
    names = {row['input chain']: row['input description'] for row in group}
    names.update(request.sequence_names)
    pairs = [tuple(chains)] if chains else combinations(sorted(residues), 2)
    for pair in pairs:
        try:
            domains = {}
            for chain in pair:
                if chain in overrides:
                    domains[chain] = overrides[chain]
                elif names.get(chain) in domain_table:
                    domains[chain] = domain_table[names[chain]]
                else:
                    detail = f' for sequence {names[chain]!r}' if chain in names else ''
                    raise ValueError(f'missing domain assignment for chain {chain}{detail}')
            profile, oriented = resolve(pair, domains)
            metrics = profile.extract(
                locator, oriented, record=record, seeds=seeds, model_scores=model_scores, dssp=dssp,
            )
            if group:
                source = record[0]
                ids = (metrics.fold_ids[source], metrics.sequence_ids[source])
                if ids != (group[0]['fold_id'], group[0]['sequence_id']):
                    raise ValueError('model/manifest identity mismatch')
            current = []
            for patch in metrics.leaves():
                structure = metrics.descriptions.models[(patch['source'], patch['model'])]
                classification = profile.classify(patch, structure, domains, chains=oriented)
                patch.update(classification.settings)
                patch.update({
                    'interface_family': profile.name,
                    'confidence': profile.confidence(patch),
                    'chemical feasibility': profile.feasibility(patch),
                    'interface': classification.interface,
                })
                for side in (1, 2):
                    chain = patch[f'chain {side}']
                    details = classification.sides[chain]
                    patch.update({f'sequence name {side}': names.get(chain),
                                  f'domain {side}': domains[chain],
                                  f'candidate face {side}': details.get('candidate face'),
                                  f'face status {side}': details.get('face status')})
                current.append(patch)
            yield from current
        except ERRORS as exc:
            error(f'{locator} model {number} chains {"/".join(pair)}: {exc}')


def run(requests, output, *, manifest=None, domain_file=None, domains=(), chains=None,
        full_data=None, dssp=None, stderr=None):
    stderr = sys.stderr if stderr is None else stderr
    if chains is not None and (len(chains) != 2 or chains[0] == chains[1]):
        raise ValueError('--chains requires two distinct chains')
    overrides = domain_overrides(domains)
    explicit_mapping = chains is not None and all(chain in overrides for chain in chains)
    lookup_mapping = manifest is not None and domain_file is not None
    if not (explicit_mapping or lookup_mapping):
        raise ValueError(
            'domain mapping requires --chains with a --domain assignment for each selected chain, '
            'or both --manifest and --domain-file'
        )
    requests = list(requests)
    target = Path(output).resolve()
    protected = {Path(path).resolve() for path in (manifest, domain_file, full_data) if path}
    protected.update(Path(request.locator.split('::', 1)[0]).resolve() for request in requests)
    if target in protected or target.suffix.lower() != '.tsv':
        raise ValueError('output must be a separate .tsv file')
    if not requests:
        return 0
    groups = read_manifest_groups(manifest) if manifest else {}
    domain_table = read_domains(domain_file) if domain_file else {}
    executable = find_dssp(dssp)
    failed = False

    def error(message):
        nonlocal failed
        failed = True
        print(f'ERROR: {message}', file=stderr)

    rows = []
    seen = set()
    for request in requests:
        for locator, number in jobs_for(request, error):
            key = locator, number
            if key in seen:
                continue
            seen.add(key)
            try:
                rows.extend(evaluate_model(
                    locator, number, request, groups=groups, domain_table=domain_table,
                    overrides=overrides, chains=chains, full_data=full_data,
                    dssp=executable, error=error,
                ))
            except ERRORS as exc:
                error(f'{locator} model {number}: {exc}')
    if rows or not failed:
        write_cache(target, COLUMNS, patch_export_rows(
            rows, side_fields=SIDE_FIELDS, extra_columns=EXTRA_COLUMNS, result_columns=RESULT_COLUMNS,
        ))
    return int(failed)
