"""Spatial interface patches; algorithm specification: docs/alphafold-pdockq2-patches.md."""
import csv
import json
import math
import statistics
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import zipfile

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
        ranked = sorted(patch_rows, key=lambda r: (-r['pDockQ2 max'], -r['contact count'], r['graph residues']))
        for rank, row in enumerate(ranked, 1):
            row['rank'] = rank
        rows.extend(patch_rows)
    return rows


def input_models(path, full_data=None):
    path = Path(path)
    if path.is_dir():
        if full_data:
            raise ValueError('--full-data is only valid with CIF input')
        archives = sorted(path.glob('*.zip'))
        if not archives:
            raise ValueError(f'directory contains no ZIP files: {path}')
        for archive in archives:
            yield from input_models(archive)
    elif path.suffix.lower() == '.cif':
        if not full_data:
            raise ValueError('CIF input requires --full-data')
        yield str(path.resolve()), model_number(path.name), path.read_text(), json.loads(Path(full_data).read_text())
    elif path.suffix.lower() == '.zip':
        if full_data:
            raise ValueError('--full-data is only valid with CIF input')
        with zipfile.ZipFile(path) as archive:
            models = sorted((model_number(n), n) for n in archive.namelist() if model_number(n) is not None)
            if not models or len({n for n, _ in models}) != len(models):
                raise ValueError('ZIP must contain models with unique model numbers')
            for number, name in models:
                data = name.removesuffix(f'model_{number}.cif') + f'full_data_{number}.json'
                if archive.namelist().count(data) != 1:
                    raise ValueError(f'model {number} requires exactly one matching {data}')
                yield str(path.resolve()), number, archive.read(name).decode(), json.loads(archive.read(data))
    else:
        raise ValueError('input must be a CIF, ZIP, or directory of ZIPs')


def analyze(path, full_data=None, options=None, chains=None, summary_targets=None):
    rows = []
    for source, model, cif, data in input_models(path, full_data):
        residues, pae = parse_model(cif, data)
        if summary_targets is not None:
            summary_targets[(source, '' if model is None else model)] = next(reversed(residues), None)
        for values in residues.values():
            for r in values:
                if not all(math.isfinite(v) for v in (*r.ca, r.plddt)) or not 0 <= r.plddt <= 100:
                    raise ValueError('coordinates must be finite and pLDDT must be in [0, 100]')
        if any(not math.isfinite(float(v)) or float(v) < 0 for row in pae for v in row):
            raise ValueError('PAE must be finite and nonnegative')
        for row in detect_patches(residues, pae, options, chains):
            row.update(source=source, model='' if model is None else model)
            rows.append(row)
    return rows


def summary_residues(row):
    return {(row[f'chain {side}'], int(number))
            for side in (1, 2)
            for number in row[f'interface residues {side}'].split(',') if number}


def distinct_patches(candidates, limit=5):
    """Greedy score-order selection; <=10% overlap relative to either patch."""
    selected, residue_sets = [], []
    for row in candidates:
        residues = summary_residues(row)
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
    if not row['contact count']:
        return 'PAE-only confidence: n/a (no contacts)'
    return (f"PAE-only confidence (0–1, higher better): "
            f"{row['chain 1']}→{row['chain 2']} {row['normalized PAE 1 to 2']:.6f} | "
            f"{row['chain 2']}→{row['chain 1']} {row['normalized PAE 2 to 1']:.6f}")


def write_summary(rows, targets, stream):
    """Full pairwise interfaces, then top five distinct last-chain patches."""
    for source in sorted({source for source, _ in targets}):
        last_chains = sorted({chain for (s, _), chain in targets.items() if s == source and chain is not None})
        print(f'\n{source}', file=stream)
        print('Full-interface pDockQ2 (all analyzed chain pairs, before patch splitting):', file=stream)
        references = sorted((r for r in rows if r['source'] == source and r['scope'] == 'full'),
                            key=lambda r: (int(r['model']) if r['model'] != '' else -1,
                                           r['chain 1'], r['chain 2']))
        if not references:
            print('  No chain-pair results.', file=stream)
        for row in references:
            model = row['model'] if row['model'] != '' else 'unnumbered'
            a, b = row['chain 1'], row['chain 2']
            print(f"  Model {model} | {a}–{b} | "
                  f"{a}→{b} {row['pDockQ2 1 to 2']:.6f} | "
                  f"{b}→{a} {row['pDockQ2 2 to 1']:.6f} | "
                  f"max pDockQ2 {row['pDockQ2 max']:.6f}", file=stream)
            print(f"     Contacts: {row['contact count']} | residues {a}: {row['residue count 1']}, "
                  f"{b}: {row['residue count 2']} | {pae_confidence_text(row)}", file=stream)
        print(file=stream)
        print(f'Top 5 patches involving the last CIF chain ({", ".join(last_chains)}), across models (at most 10% residue overlap):', file=stream)
        candidates = [r for r in rows if r['source'] == source and r['scope'] == 'patch'
                      and targets.get((source, r['model'])) in (r['chain 1'], r['chain 2'])]
        candidates.sort(key=lambda r: (-r['pDockQ2 max'], -r['contact count'],
                                      str(r['model']), r['chain 1'], r['chain 2'], r['graph residues'], r['patch id']))
        if not candidates:
            print('  No qualifying patches for the last chain under these settings.', file=stream)
        for rank, row in enumerate(distinct_patches(candidates), 1):
            model = row['model'] if row['model'] != '' else 'unnumbered'
            print(f"  {rank}. Model {model} | patch {row['patch id']} ({row['status']}) | "
                  f"{row['chain 1']}–{row['chain 2']} | max pDockQ2 {row['pDockQ2 max']:.6f}", file=stream)
            print(f"     Region 1 ({row['chain 1']}): {row['region 1']}", file=stream)
            print(f"     Region 2 ({row['chain 2']}): {row['region 2']}", file=stream)
            print(f'     {pae_confidence_text(row)}', file=stream)
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


def write_patches(rows, stream):
    writer = csv.DictWriter(stream, fieldnames=OUTPUT_COLUMNS + EXTRA_COLUMNS, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for row in rows:
        writer.writerow({k: f'{v:.6f}' if isinstance(v, float) else v for k, v in row.items()})
