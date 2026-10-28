"""Spatial interface patches; algorithm specification: docs/pDockQ2.md."""
import csv
import json
import math
import re
import statistics
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import zipfile

from sieve.alphafold_confidence import model_confidence
from sieve.alphafold_files import open_fold, fold_locator, read_seed, read_job_request
from sieve.alphafold_manifest import chain_sequences, fold_id, sequence_id
from sieve.alphafold_composition import chain_compositions

from sieve.alphafold_pdockq2 import (
    OUTPUT_COLUMNS, _distance, model_number, parse_model, score_selections,
)

ALGORITHM_VERSION = '1.0'
EXTRA_COLUMNS = [
    'algorithm version', 'patch id', 'parent id', 'depth', 'status', 'rank',
    'graph residues', 'separator residues', 'discarded residues',
    'residue count 1', 'residue count 2', 'distance cutoff',
    'max bridge residues', 'min residues per chain', 'min contacts',
    'max separator checks',
    'contact distance min', 'contact distance max',
    'contact distance median', 'contact distance mean',
    *[f'nearest contact distance {stat} {side}'
      for side in (1, 2) for stat in ('min', 'max', 'median', 'mean')],
]


@dataclass(frozen=True)
class Options:
    cutoff: float = 8.0
    max_bridge_residues: int = 1
    min_residues_per_chain: int = 3
    min_contacts: int = 5
    max_separator_checks: int = 1000000

    def __post_init__(self):
        if not math.isfinite(self.cutoff) or self.cutoff <= 0:
            raise ValueError('distance cutoff must be finite and greater than zero')
        for name in ('max_bridge_residues', 'min_residues_per_chain', 'min_contacts', 'max_separator_checks'):
            value = getattr(self, name)
            minimum = 0 if name == 'max_bridge_residues' else 1
            if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
                raise ValueError(f'{name} must be an integer >= {minimum}')


def components(graph, nodes):
    remaining = set(nodes)
    result = []
    while remaining:
        seed = min(remaining)
        remaining.remove(seed)
        found, stack = {seed}, [seed]
        while stack:
            neighbors = graph[stack.pop()] & remaining
            remaining.difference_update(neighbors)
            found.update(neighbors)
            stack.extend(sorted(neighbors))
        result.append(frozenset(found))
    return result


def find_split(graph, nodes, qualifies, options, budget):
    """First qualifying separator by cardinality, then lexical residue order."""
    for size in range(1, min(options.max_bridge_residues, len(nodes) - 2) + 1):
        for separator in combinations(sorted(nodes), size):
            budget[0] += 1
            if budget[0] > options.max_separator_checks:
                raise ValueError('separator search limit exceeded; reduce --max-bridge-residues or increase --max-separator-checks')
            parts = components(graph, nodes - set(separator))
            if len(parts) < 2:
                continue
            accepted = [part for part in parts if qualifies(part)]
            if len(accepted) >= 2:
                discarded = set().union(*(part for part in parts if part not in accepted))
                return separator, accepted, discarded
    return (), [], set()


def ranges(numbers):
    result = []
    for number in sorted(set(numbers)):
        if result and number == result[-1][1] + 1:
            result[-1] = (result[-1][0], number)
        else:
            result.append((number, number))
    return result


def labels(nodes):
    return ','.join(f'{chain}:{number}' for chain, number in sorted(nodes))


def detect_patches(residues, pae, options=None, chains=None):
    options = options or Options()
    if chains is not None and (len(chains) != 2 or chains[0] == chains[1] or any(c not in residues for c in chains)):
        raise ValueError('--chains requires two distinct chains present in the model')
    pairs = [tuple(sorted(chains))] if chains else combinations(sorted(residues), 2)
    rows = []
    for chain1, chain2 in pairs:
        lookup = {(r.chain, r.number): r for c in (chain1, chain2) for r in residues[c]}
        contacts = [((a.chain, a.number), (b.chain, b.number))
                    for a in residues[chain1] for b in residues[chain2]
                    if _distance(a.ca, b.ca) <= options.cutoff]
        nodes = set(node for contact in contacts for node in contact)
        graph = {node: set() for node in nodes}
        for a, b in combinations(sorted(nodes), 2):
            if _distance(lookup[a].ca, lookup[b].ca) <= options.cutoff:
                graph[a].add(b)
                graph[b].add(a)

        def qualifies(part):
            edges = [(a, b) for a, b in contacts if a in part and b in part]
            return (len(edges) >= options.min_contacts
                    and len({a for a, _ in edges}) >= options.min_residues_per_chain
                    and len({b for _, b in edges}) >= options.min_residues_per_chain)

        def make_row(part, scope):
            selected = [[lookup[n] for n in sorted(part) if n[0] == c] for c in (chain1, chain2)]
            row = score_selections(chain1, chain2, *selected, pae,
                                   *[ranges(r.number for r in side) for side in selected], options.cutoff)
            contact_distances = [(a, b, _distance(lookup[a].ca, lookup[b].ca))
                                 for a, b in contacts if a in part and b in part]
            distances = [distance for _, _, distance in contact_distances]
            row.update({
                'contact distance min': min(distances) if distances else '',
                'contact distance max': max(distances) if distances else '',
                'contact distance mean': statistics.mean(distances) if distances else '',
                'contact distance median': statistics.median(distances) if distances else '',
            })
            row.update({'scope': scope, 'algorithm version': ALGORITHM_VERSION,
                        'graph residues': labels(part), 'separator residues': '', 'discarded residues': '',
                        'distance cutoff': options.cutoff, 'max bridge residues': options.max_bridge_residues,
                        'min residues per chain': options.min_residues_per_chain, 'min contacts': options.min_contacts,
                        'max separator checks': options.max_separator_checks})
            for side in (1, 2):
                text = row[f'interface residues {side}']
                row[f'residue count {side}'] = len(text.split(',')) if text else 0
                nearest = {}
                for a, b, distance in contact_distances:
                    node = a if side == 1 else b
                    nearest[node] = min(nearest.get(node, math.inf), distance)
                values = list(nearest.values())
                for stat, function in [('min', min), ('max', max),
                                       ('median', statistics.median), ('mean', statistics.mean)]:
                    row[f'nearest contact distance {stat} {side}'] = function(values) if values else ''
            return row

        reference = make_row(nodes, 'full')
        reference.update({'patch id': 'full', 'parent id': '', 'depth': 0, 'status': 'reference', 'rank': ''})
        rows.append(reference)
        patch_rows, budget = [], [0]
        stack = [(part, '', 0) for part in reversed(components(graph, nodes)) if qualifies(part)]
        while stack:
            part, parent, depth = stack.pop()
            patch_id = f'p{len(patch_rows) + 1}'
            row = make_row(part, 'patch')
            separator, children, discarded = find_split(graph, part, qualifies, options, budget)
            row.update({'patch id': patch_id, 'parent id': parent, 'depth': depth,
                        'status': 'split' if children else 'leaf',
                        'separator residues': labels(separator), 'discarded residues': labels(discarded)})
            patch_rows.append(row)
            stack.extend((child, patch_id, depth + 1) for child in reversed(children))
        ranked = sorted(patch_rows, key=lambda r: (-(r['pDockQ2 max'] if r['pDockQ2 max'] is not None else -1), -r['contact count'], r['graph residues']))
        for rank, row in enumerate(ranked, 1):
            row['rank'] = rank
        rows.extend(patch_rows)
    return rows


def input_models(path, full_data=None, *, model_scores=None, model_seeds=None, selected_model=None):
    def select_models(models):
        if selected_model is not None:
            models = [(number, name) for number, name in models if number == selected_model]
            if not models:
                raise ValueError(f'model {selected_model} not found in {path}')
        return models

    requests = {}

    def with_request(source, number, data):
        request = requests.get((source, number))
        if request is not None and 'sequences' in request:
            data = dict(data, _job_request=request)
        return data

    def record_scores(source, number, name, read, names):
        if name is not None:
            request_name = re.sub(r'summary_confidences_\d+\.json$', 'job_request.json', name)
            model_name = re.sub(r'summary_confidences_(\d+)\.json$', r'model_\1.cif', name)
            if request_name != name and request_name in names:
                requests[(source, number)] = read_job_request(read, [model_name], names)
        if model_seeds is not None:
            seed = None
            if name is not None:
                request = re.sub(r'summary_confidences_\d+\.json$', 'job_request.json', name)
                model = re.sub(r'summary_confidences_(\d+)\.json$', r'model_\1.cif', name)
                if request in names:
                    seed = read_seed(read, [model], names)
            model_seeds[(source, '' if number is None else number)] = seed
        if model_scores is None:
            return
        summary = {}
        if name is not None and name in names:
            if names.count(name) != 1:
                raise ValueError(f'expected exactly one summary confidence file: {name}')
            summary = json.loads(read(name))
        model_scores[(source, '' if number is None else number)] = model_confidence(summary)

    def summary_name(name, number):
        return str(name).removesuffix(f'model_{number}.cif') + f'summary_confidences_{number}.json'

    # Derived assemblies have coordinates and identity, but no joint PAE matrix.
    from sieve.alphafold_structure import open_structure_fold, assembly_metadata, COMPLEX_FILE
    candidate = Path(str(path))
    if '::' in str(path) or (candidate.is_dir() and (candidate/COMPLEX_FILE).is_file()) or candidate.name == COMPLEX_FILE:
        with open_structure_fold(path) as fold:
            meta = assembly_metadata(fold)
            if meta:
                if full_data:
                    raise ValueError('assembled inputs do not accept --full-data')
                select_models([(0, fold.models[0])])
                source = str(candidate.resolve()) if '::' not in str(path) else str(path)
                if model_scores is not None:
                    model_scores[(source,0)] = model_confidence({})
                if model_seeds is not None:
                    model_seeds[(source,0)] = None
                yield source, 0, fold.read(fold.models[0]), {'_assembly':meta}
                return
    if '::' in str(path):
        if full_data:
            raise ValueError('--full-data is only valid with CIF input')
        with open_fold(path) as fold:
            models = sorted((model_number(name), name) for name in fold.models)
            models = select_models(models)
            if len({number for number, _ in models}) != len(models):
                raise ValueError('fold must contain models with unique model numbers')
            for number, name in models:
                data = name.removesuffix(f'model_{number}.cif') + f'full_data_{number}.json'
                source = fold_locator(fold.source_path, fold.fold_path)
                record_scores(source, number, summary_name(name, number), fold.read, fold.member_names)
                yield (source, number,
                       fold.read(name), with_request(source, number, json.loads(fold.read(data))))
        return
    path = Path(path)
    if path.is_dir():
        if full_data:
            raise ValueError('--full-data is only valid with CIF input')
        models = sorted((model_number(p.name), p) for p in path.iterdir()
                        if p.is_file() and model_number(p.name) is not None)
        if models:
            models = select_models(models)
            if len({n for n, _ in models}) != len(models):
                raise ValueError('directory must contain models with unique model numbers')
            for number, cif in models:
                data = cif.with_name(cif.name.removesuffix(f'model_{number}.cif') + f'full_data_{number}.json')
                if not data.is_file():
                    raise ValueError(f'model {number} requires matching {data.name}')
                record_scores(str(path.resolve()), number, summary_name(cif.name, number),
                              lambda name: (path / name).read_text(), [p.name for p in path.iterdir()])
                yield str(path.resolve()), number, cif.read_text(), with_request(str(path.resolve()), number, json.loads(data.read_text()))
            return
        if selected_model is not None:
            raise ValueError(f'model {selected_model} not found in {path}')
        archives = sorted(p for p in path.glob('*.zip') if p.is_file() and not p.name.startswith('._'))
        if not archives:
            raise ValueError(f'directory contains no model CIF files or ZIP files: {path}')
        for archive in archives:
            yield from input_models(archive, model_scores=model_scores, model_seeds=model_seeds)
    elif path.suffix.lower() == '.cif':
        if not full_data:
            raise ValueError('CIF input requires --full-data')
        number = model_number(path.name)
        select_models([(number, path)])
        summary = path.with_name(summary_name(path.name, number)) if number is not None else None
        if summary is None:
            data_path = Path(full_data)
            match = re.search(r'full_data_(\d+)\.json$', data_path.name)
            if match:
                summary = data_path.with_name(data_path.name[:match.start()] + f'summary_confidences_{match[1]}.json')
        record_scores(str(path.resolve()), number, str(summary) if summary is not None else None,
                      lambda name: Path(name).read_text(),
                      [str(p) for p in summary.parent.iterdir()] if summary is not None else [])
        yield str(path.resolve()), number, path.read_text(), with_request(str(path.resolve()), number, json.loads(Path(full_data).read_text()))
    elif path.suffix.lower() == '.zip':
        if full_data:
            raise ValueError('--full-data is only valid with CIF input')
        with zipfile.ZipFile(path) as archive:
            models = sorted((model_number(n), n) for n in archive.namelist() if model_number(n) is not None)
            models = select_models(models)
            if not models or len({n for n, _ in models}) != len(models):
                raise ValueError('ZIP must contain models with unique model numbers')
            for number, name in models:
                data = name.removesuffix(f'model_{number}.cif') + f'full_data_{number}.json'
                if archive.namelist().count(data) != 1:
                    raise ValueError(f'model {number} requires exactly one matching {data}')
                record_scores(str(path.resolve()), number, summary_name(name, number),
                              archive.read, archive.namelist())
                yield str(path.resolve()), number, archive.read(name).decode(), with_request(str(path.resolve()), number, json.loads(archive.read(data)))
    else:
        raise ValueError('input must be a CIF, ZIP, or directory of extracted models or ZIPs')


def analyze(path, full_data=None, options=None, chains=None, summary_targets=None, summary_inputs=None, model_scores=None, fold_ids=None, sequence_ids=None, selected_model=None):
    rows = []
    seeds = {} if fold_ids is not None else None
    compositions = {}
    for source, model, cif, data in input_models(path, full_data, model_scores=model_scores, model_seeds=seeds, selected_model=selected_model):
        if fold_ids is not None or sequence_ids is not None:
            try:
                sequences = list(chain_compositions(cif, data.get('_job_request')).values())
                composition = sequence_id(sequences)
            except (ValueError, KeyError):
                if '_entity_poly_seq.' in cif or data.get('_job_request', {}).get('sequences'):
                    raise
                sequences, composition = [], None
            if source in compositions:
                previous = compositions[source]
                if previous is not None and composition is not None and previous != composition:
                    raise ValueError(f'fold sequence composition differs across models: {source}')
                if previous is None:
                    composition = None
            compositions[source] = composition
            if sequence_ids is not None:
                sequence_ids[source] = composition
            if fold_ids is not None:
                seed = seeds[(source, '' if model is None else model)]
                identifier = data['_assembly']['fold_id'] if '_assembly' in data else (fold_id(seed, sequences) if seed is not None and composition is not None else None)
                if source in fold_ids:
                    previous = fold_ids[source]
                    if previous is not None and identifier is not None and previous != identifier:
                        raise ValueError(f'fold seed differs across models: {source}')
                    if previous is None:
                        identifier = None
                fold_ids[source] = identifier
        residues, pae = parse_model(cif, data)
        if summary_inputs is not None:
            summary_inputs[(source, '' if model is None else model)] = cif
        for values in residues.values():
            for r in values:
                if not all(math.isfinite(v) for v in r.ca) or (r.plddt is not None and (not math.isfinite(r.plddt) or not 0 <= r.plddt <= 100)):
                    raise ValueError('coordinates must be finite and pLDDT must be in [0, 100]')
        if any(not math.isfinite(float(v)) or float(v) < 0 for row in (pae or []) for v in row):
            raise ValueError('PAE must be finite and nonnegative')
        patches = detect_patches(residues, pae, options, chains)
        if summary_targets is not None:
            summary_targets[(source, '' if model is None else model)] = (
                chains[1] if chains is not None else next(reversed(residues), None))
        for row in patches:
            row.update(source=source, model='' if model is None else model)
            if '_assembly' in data:
                row['confidence origin'] = 'inherited per-chain values; no assembled confidence'
            rows.append(row)
    return rows


def summary_residues(row, chain):
    return {(row[f'chain {side}'], int(number))
            for side in (1, 2) if row[f'chain {side}'] == chain
            for number in row[f'interface residues {side}'].split(',') if number}


def distinct_patches(candidates, targets, limit=5):
    """Greedy score-order selection; <=10% overlap on the summary target chain only."""
    selected, residue_sets = [], []
    for row in candidates:
        residues = summary_residues(row, targets[(row['source'], row['model'])])
        if not residues:
            continue
        if any(10 * len(residues & previous) > min(len(residues), len(previous))
               for previous in residue_sets):
            continue
        selected.append(row)
        residue_sets.append(residues)
        if len(selected) == limit:
            break
    return selected


def pae_confidence_text(row):
    if row['normalized PAE 1 to 2'] is None:
        return 'PAE-only confidence: NA (assembled model has no joint PAE)'
    if not row['contact count']:
        return 'PAE-only confidence: n/a (no contacts)'
    return (f"PAE-only confidence (0–1, higher better): "
            f"{row['chain 1']}→{row['chain 2']} {row['normalized PAE 1 to 2']:.6f} | "
            f"{row['chain 2']}→{row['chain 1']} {row['normalized PAE 2 to 1']:.6f}")


def write_summary(rows, targets, stream, descriptions=None, *, structure=False, model_scores=None, fold_ids=None, sequence_ids=None, chains=None):
    """Fold identity and confidence, then top five distinct target-chain patches."""
    for source in sorted({source for source, _ in targets}):
        last_chains = sorted({chain for (s, _), chain in targets.items() if s == source and chain is not None})
        identifier = (fold_ids or {}).get(source)
        print(f"\nFold ID: {identifier or 'n/a (seed or sequence metadata unavailable)'}", file=stream)
        sequence_identifier = (sequence_ids or {}).get(source)
        print(f"Sequence ID: {sequence_identifier or 'n/a (sequence metadata unavailable)'}", file=stream)
        print(source, file=stream)
        print('Model confidence (whole fold):', file=stream)
        for model in sorted((m for s, m in targets if s == source), key=lambda m: int(m) if m != '' else -1):
            scores = (model_scores or {}).get((source, model), {})
            def score_text(field):
                value = scores.get(field)
                return 'n/a (unavailable)' if value is None else f'{value:.6f}'
            label = model if model != '' else 'unnumbered'
            print(f"  Model {label} | ipTM {score_text('iptm')} | pTM {score_text('ptm')}", file=stream)
        if structure and descriptions is not None:
            descriptions.write_chains(source, stream)
        print(file=stream)
        if chains is not None:
            print(f'Top 5 distinct sites for chains {chains[0]}–{chains[1]}, across models '
                  f'(at most 10% residue overlap on chain {chains[1]}):', file=stream)
        else:
            print(f'Top 5 distinct sites on the last CIF chain ({", ".join(last_chains)}), across models (at most 10% last-chain residue overlap):', file=stream)
        candidates = [r for r in rows if r['source'] == source and r['scope'] == 'patch' and r['status'] == 'leaf'
                      and targets.get((source, r['model'])) in (r['chain 1'], r['chain 2'])
                      and (chains is None or {r['chain 1'], r['chain 2']} == set(chains))]
        if any(r['pDockQ2 max'] is None for r in candidates):
            print('  Geometric ranking by contact count; pDockQ2 unavailable for assembled coordinates.', file=stream)
        candidates.sort(key=lambda r: (-(r['pDockQ2 max'] if r['pDockQ2 max'] is not None else -1), -r['contact count'],
                                      str(r['model']), r['chain 1'], r['chain 2'], r['graph residues'], r['patch id']))
        if not candidates:
            subject = 'the selected chains' if chains is not None else 'the last chain'
            print(f'  No qualifying patches for {subject} under these settings.', file=stream)
        for rank, row in enumerate(distinct_patches(candidates, targets), 1):
            model = row['model'] if row['model'] != '' else 'unnumbered'
            patch_score = 'NA' if row['pDockQ2 max'] is None else f"{row['pDockQ2 max']:.6f}"
            print(f"  {rank}. Model {model} | patch {row['patch id']} ({row['status']}) | "
                  f"{row['chain 1']}–{row['chain 2']} | max pDockQ2 {patch_score}", file=stream)
            print(f"     Region 1 ({row['chain 1']}): {row['region 1']}", file=stream)
            print(f"     Region 2 ({row['chain 2']}): {row['region 2']}", file=stream)
            print(f'     {pae_confidence_text(row)}', file=stream)
            if descriptions is not None:
                descriptions.write_patch(row, stream)
            target = targets[(source, row['model'])]
            side = 1 if row['chain 1'] == target else 2
            print(f'     Nearest contact per {target} residue (Å; n={row[f"residue count {side}"]}): ' + ' | '.join(
                f"{label} {row[f'nearest contact distance {key} {side}']:.3f}"
                for key, label in [('min', 'min'), ('max', 'max'), ('median', 'median'), ('mean', 'average')]
            ), file=stream)
            print('     All interchain contact pairs (Å): ' + ' | '.join(
                f"{label} {row[f'contact distance {key}']:.3f}"
                for key, label in [('min', 'min'), ('max', 'max'), ('median', 'median'), ('mean', 'average')]
            ), file=stream)


# Export schema is independent of internal graph/hierarchy fields.
EXPORT_RENAMES = {
    **{f'residue count {i}': f'patch {i} residues' for i in (1,2)},
    **{f'contact distance {stat}': f'interchain CA distance {stat}'
       for stat in ('min','max','median','mean')},
    **{f'nearest contact distance {stat} {i}': f'nearest interchain CA distance {stat} {i}'
       for i in (1,2) for stat in ('min','max','median','mean')},
}
EXPORT_DROPPED = {'interface residues 1','interface residues 2','algorithm version',
                  'parent id','depth','graph residues','separator residues','discarded residues'}
METRIC_COLUMNS = ['interface structure 1','interface structure 2',
    'total buried surface area','total C/S contribution percent',
    'nonpolar contacts','salt bridge candidates','steric overlap candidates','phosphate/basic contacts',
    'salt bridge candidate heavy-atom distance min','salt bridge residues',
    'nominal patch charge difference','patch 1 hydropathy mean','patch 2 hydropathy mean',
    'candidate DD family face 1','candidate DD family face 2']
PATCH_COLUMNS = ['sequence_id','fold_id',
    *[EXPORT_RENAMES.get(c,c) for c in OUTPUT_COLUMNS+EXTRA_COLUMNS if c not in EXPORT_DROPPED],
    *METRIC_COLUMNS, 'dd enabled','dd helix presence threshold','confidence origin']


def write_patches(rows, stream):
    """Export only final leaf patches; internal reference/split rows remain available."""
    leaves = [r for r in rows if r['scope']=='patch' and r['status']=='leaf']
    groups = {}
    for row in leaves:
        groups.setdefault((row.get('source'),row.get('model'),row['chain 1'],row['chain 2']),[]).append(row)
    ranks = {}
    for group in groups.values():
        for rank,row in enumerate(sorted(group,key=lambda r:r['rank']),1):
            ranks[id(row)] = rank
    writer = csv.DictWriter(stream, fieldnames=PATCH_COLUMNS, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for row in leaves:
        exported = {EXPORT_RENAMES.get(k,k):v for k,v in row.items() if k not in EXPORT_DROPPED}
        exported['rank'] = ranks[id(row)]
        writer.writerow({k: 'NA' if exported.get(k) is None else f'{exported[k]:.6f}' if isinstance(exported[k],float)
                         else exported[k] for k in PATCH_COLUMNS})
