"""Validated sequence manifests with per-location, atomic TSV checkpoints."""
import csv
import hashlib
import io
import json
import re
import sys
import zipfile
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from sieve.alphafold_cache import read_cache, write_cache
from sieve.alphafold_files import discover, fold_locator, model_number
from sieve.alphafold_composition import chain_compositions, hash_payload

BASE_COLUMNS = ['source_path', 'fold_path', 'fold_locator', 'fold_id', 'sequence_id', 'timestamp', 'source_type']
INPUT_COLUMNS = ['input description', 'input chain']
COLUMNS = BASE_COLUMNS + INPUT_COLUMNS



class MissingSequences(ValueError):
    def __init__(self, chains):
        self.chains = chains
        super().__init__('; '.join(f'{sequence} not found for chain {chain}'
                                  for chain, sequence in chains))


def read_candidates(path, error):
    """Preserve headers verbatim; normalize sequence case and whitespace only."""
    sequences = {}
    headers = {}
    header = None
    parts = []
    start = 0
    count = 0

    def finish():
        nonlocal count
        if header is None:
            return
        count += 1
        if header in headers:
            error(f'{path}:{start}: duplicate FASTA header {header!r}; first seen at line {headers[header]}')
        else:
            headers[header] = start
        sequence = ''.join(parts).upper()
        label = f'{path}:{start} ({header})'
        if not header.strip():
            error(f'{label}: empty FASTA header')
        if not sequence:
            error(f'{label}: empty sequence')
            return
        if not re.fullmatch('[A-Z]+', sequence):
            error(f'{label}: invalid sequence characters: {sequence}')
        if sequence in sequences:
            error(f'{label}: duplicate sequence also present under {sequences[sequence]!r}: {sequence}')
        else:
            sequences[sequence] = header

    with open(path, encoding='utf-8') as stream:
        for number, line in enumerate(stream, 1):
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                finish()
                header, parts, start = line[1:], [], number
            elif line.strip():
                if header is None:
                    error(f'{path}:{number}: sequence before first FASTA header')
                else:
                    parts.append(''.join(line.split()))
    finish()
    if not count:
        error(f'{path}: no FASTA records')
    return sequences


def sequence_id(sequences):
    serialized = json.dumps(hash_payload(sequences), separators=(',', ':'), ensure_ascii=True)
    return hashlib.sha1(serialized.encode('utf-8')).hexdigest()[-8:]


def fold_id(seed, sequences):
    serialized = json.dumps((seed, hash_payload(sequences)), separators=(',', ':'), ensure_ascii=True)
    return hashlib.sha1(serialized.encode('utf-8')).hexdigest()[-8:]


def chain_sequences(text, request=None):
    """Return base sequences for FASTA matching; retain PTMs separately for identity."""
    return {chain: composition.sequence for chain, composition in chain_compositions(text, request).items()}


def identify(fold, candidates):
    numbers = [model_number(n) for n in fold.models]
    if len(numbers) != len(set(numbers)):
        raise ValueError('multiple model sets or duplicate model numbers in one fold directory')
    from sieve.alphafold_structure import request_for_fold
    request = request_for_fold(fold)
    compositions = None
    author_chains = None
    for name in sorted(fold.models):
        text = fold.read(name)
        current = chain_compositions(text, request)
        mapping = author_chain_ids(text, current)
        if author_chains is not None and mapping != author_chains:
            raise ValueError(f'chain IDs differ across models: {name}')
        author_chains = mapping
        if compositions is not None and current != compositions:
            raise ValueError(f'chain sequences differ across models (including PTMs): {name}')
        compositions = current
    chains = {c: value.sequence for c, value in (compositions or {}).items()}
    if not chains:
        raise ValueError('no chains found')
    missing = [(chain, sequence)
               for chain, sequence in chains.items() if sequence not in candidates]
    if missing:
        raise MissingSequences(missing)
    from sieve.alphafold_structure import identity_for_fold
    identity = identity_for_fold(fold, fold.read(fold.models[0]))
    return make_rows(fold.source_path, fold.fold_path, identity['fold_id'], identity['sequence_id'],
                    fold.timestamp, [(candidates[seq], author_chains[chain]) for chain, seq in
                                     sorted(chains.items(), key=lambda item: author_chains[item[0]])],
                    source_type=identity['source_type'])


def author_chain_ids(text, compositions):
    """Map CIF label chains to author chains used by the patch MMCIFParser."""
    data = MMCIF2Dict(io.StringIO(text))
    labels = data.get('_atom_site.label_asym_id', [])
    authors = data.get('_atom_site.auth_asym_id', [])
    if not labels or len(labels) != len(authors):
        raise ValueError('CIF lacks complete label/author chain mapping')
    mapping = {}
    reverse = {}
    for label, author in zip(labels, authors):
        if author in ('', '.', '?'):
            raise ValueError('CIF has missing author chain ID')
        if label in mapping and mapping[label] != author:
            raise ValueError('ambiguous CIF label/author chain mapping')
        if author in reverse and reverse[author] != label:
            raise ValueError('ambiguous CIF author chain ID shared by label chains')
        mapping[label] = author
        reverse[author] = label
    if set(compositions) - set(mapping):
        raise ValueError('CIF lacks author chain mapping for a sequence')
    return {chain: mapping[chain] for chain in compositions}


def make_rows(source, member, identifier, sequence_identifier, timestamp, inputs, source_type='alphafold'):
    base = dict(source_path=source, fold_path=member, fold_locator=fold_locator(source, member),
                fold_id=identifier, sequence_id=sequence_identifier, timestamp=timestamp, source_type=source_type)
    return [dict(base, **{'input description': description, 'input chain': chain})
            for description, chain in inputs]


def group_manifest_rows(rows):
    """Validate and group tall rows by locator, preserving first-seen fold order."""
    grouped = {}
    for row in rows:
        if None in row or any(value is None for value in row.values()):
            raise ValueError('incomplete or malformed manifest row')
        required = ('fold_locator', 'fold_id', 'sequence_id', *INPUT_COLUMNS)
        if any(not row.get(c, '').strip() for c in required):
            raise ValueError('incomplete manifest row: missing ID, locator, description, or chain')
        if any(not re.fullmatch('[0-9a-f]{8}', row[c]) for c in ('fold_id', 'sequence_id')):
            raise ValueError('invalid manifest ID')
        locator = row['fold_locator']
        group = grouped.setdefault(locator, [])
        if group and any(row.get(c) != group[0].get(c) for c in BASE_COLUMNS):
            raise ValueError(f'contradictory manifest metadata for locator {locator}')
        if any(r['input chain'] == row['input chain'] for r in group):
            raise ValueError(f'duplicate manifest chain {row["input chain"]} at {locator}')
        group.append(row)
    return grouped


def read_manifest_groups(path):
    with open(path, encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t', strict=True)
        columns = reader.fieldnames
        required = ('fold_locator', 'fold_id', 'sequence_id', *INPUT_COLUMNS)
        if (columns is None or len(set(columns)) != len(columns)
                or any(c not in columns for c in required)):
            raise ValueError('manifest uses an older or incomplete format; remove it and regenerate with input description and input chain')
        return group_manifest_rows(reader)


def validate_cache(rows):
    cached = {}
    for group in group_manifest_rows(rows).values():
        row = group[0]
        if any(c not in row for c in COLUMNS):
            raise ValueError('cache contains an incomplete row')
        if row['source_type'] not in ('alphafold', 'assembled'):
            raise ValueError('invalid manifest source_type')
        key = (row['source_path'], row['fold_path'])
        if row['fold_locator'] != fold_locator(*key) or key in cached:
            raise ValueError('inconsistent cached locator')
        cached[key] = group
    return cached


def load_manifest(path):
    columns, rows = read_cache(path, ['source_path', 'fold_path', 'fold_id'])
    if columns is None:
        return {}, False
    if len(set(columns)) != len(columns) or any(c not in columns for c in COLUMNS):
        raise ValueError('manifest uses an older or incomplete format; remove it and regenerate with input description and input chain')
    return validate_cache(rows), False


def checkpoint(output, cached):
    write_cache(output, COLUMNS, (row for group in cached.values() for row in group))


def run(input_path, fasta_path, output_path, stderr=None):
    stderr = stderr if stderr is not None else sys.stderr
    errors = []

    def error(message):
        errors.append(message)
        print(f'ERROR: {message}', file=stderr)

    output = Path(output_path).resolve()
    if output in (Path(input_path).resolve(), Path(fasta_path).resolve()):
        raise ValueError('manifest must not overwrite an input')
    if output.suffix.lower() != '.tsv':
        raise ValueError('manifest must have a .tsv extension')
    candidates = read_candidates(fasta_path, error)
    if errors:
        return 1
    cached, _ = load_manifest(output)
    processed = skipped = found = 0
    for fold in discover(input_path, error):
        found += 1
        key = (fold.source_path, fold.fold_path)
        if key in cached:
            if not cached[key][0]['timestamp']:
                for row in cached[key]:
                    row['timestamp'] = fold.timestamp
                checkpoint(output, cached)
            skipped += 1
            continue
        try:
            updated = identify(fold, candidates)
        except MissingSequences as exc:
            for chain, sequence in exc.chains:
                error(f'{sequence} not found in {fold.source_path}::{fold.fold_path} chain {chain}')
            continue
        except (OSError, ValueError, KeyError, TypeError, IndexError, UnicodeError, RuntimeError, zipfile.BadZipFile) as exc:
            error(f'{fold.source_path}::{fold.fold_path}: {exc}')
            continue
        cached[key] = updated
        checkpoint(output, cached)
        processed += 1
    if not found and not errors:
        error(f'{input_path}: no folds found')
    print(f'Folds: {found}; written: {processed}; cached: {skipped}; errors: {len(errors)}', file=stderr)
    return 1 if errors else 0
