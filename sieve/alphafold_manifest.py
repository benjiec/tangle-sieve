"""Validated sequence manifests with per-location, atomic TSV checkpoints."""
import hashlib
import io
import json
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.PDBData import protein_letters_3to1

from sieve.alphafold_cache import read_cache, write_cache
from sieve.alphafold_files import discover, fold_locator, model_number

COLUMNS = ['source_path', 'fold_path', 'fold_locator', 'fold_id', 'description', 'chain_ids',
           'copies']


class MissingSequences(ValueError):
    def __init__(self, chains):
        self.chains = chains
        super().__init__('; '.join(f'{sequence} not found for chain {chain}'
                                  for chain, sequence in chains))


def read_candidates(path, error):
    """Preserve headers verbatim; normalize sequence case and whitespace only."""
    sequences = {}
    header = None
    parts = []
    start = 0
    count = 0

    def finish():
        nonlocal count
        if header is None:
            return
        count += 1
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


def fold_id(sequences):
    serialized = json.dumps(sorted(sequences), separators=(',', ':'), ensure_ascii=True)
    return hashlib.sha1(serialized.encode('utf-8')).hexdigest()[-8:]


def chain_sequences(text):
    """Use declared entity sequences rather than only residues with coordinates."""
    cif = MMCIF2Dict(io.StringIO(text))

    def table(*columns):
        values = [cif.get(c, []) for c in columns]
        if not values[0] or len({len(v) for v in values}) != 1:
            raise ValueError('missing or inconsistent CIF fields: ' + ', '.join(columns))
        return zip(*values)

    entity_types = dict(table('_entity.id', '_entity.type'))
    residues = defaultdict(dict)
    for entity, number, monomer in table('_entity_poly_seq.entity_id', '_entity_poly_seq.num', '_entity_poly_seq.mon_id'):
        number = int(number)
        if number <= 0 or number in residues[entity]:
            raise ValueError(f'ambiguous entity {entity} residue {number}')
        if monomer not in protein_letters_3to1:
            raise ValueError(f'entity {entity}: unsupported residue {monomer} at {number}; cannot extract an exact protein sequence')
        residues[entity][number] = protein_letters_3to1[monomer]
    chains = {}
    for chain, entity in table('_struct_asym.id', '_struct_asym.entity_id'):
        if chain in chains:
            raise ValueError(f'duplicate CIF chain {chain}')
        if entity_types.get(entity) != 'polymer' or entity not in residues:
            raise ValueError(f'chain {chain}: no identifiable protein sequence (entity {entity})')
        sequence = residues[entity]
        if sorted(sequence) != list(range(1, len(sequence) + 1)):
            raise ValueError(f'chain {chain}: nonconsecutive entity sequence numbering')
        chains[chain] = ''.join(sequence[n] for n in sorted(sequence))
    return chains


def identify(fold, candidates):
    numbers = [model_number(n) for n in fold.models]
    if len(numbers) != len(set(numbers)):
        raise ValueError('multiple model sets or duplicate model numbers in one fold directory')
    chains = None
    for name in sorted(fold.models):
        current = chain_sequences(fold.read(name))
        if chains is not None and current != chains:
            raise ValueError(f'chain sequences differ across models: {name}')
        chains = current
    if not chains:
        raise ValueError('no chains found')
    missing = [(chain, sequence)
               for chain, sequence in chains.items() if sequence not in candidates]
    if missing:
        raise MissingSequences(missing)
    grouped = defaultdict(list)
    for chain, sequence in chains.items():
        grouped[sequence].append(chain)
    identifier = fold_id(chains.values())
    return [dict(source_path=fold.source_path, fold_path=fold.fold_path,
                 fold_locator=fold_locator(fold.source_path, fold.fold_path),
                 fold_id=identifier, description=candidates[sequence],
                 chain_ids=','.join(sorted(ids)), copies=str(len(ids)))
            for sequence, ids in sorted(grouped.items())]


def validate_cache(rows):
    groups = defaultdict(list)
    for row in rows:
        if any(row.get(c) is None for c in COLUMNS):
            raise ValueError('cache contains an incomplete row')
        groups[(row['source_path'], row['fold_path'])].append(row)
    for key, group in groups.items():
        chains = set()
        for row in group:
            ids = row['chain_ids'].split(',')
            if (int(row['copies']) != len(ids) or not all(ids) or
                    len(set(ids)) != len(ids) or chains.intersection(ids)):
                raise ValueError(f'invalid cached chains at {key}')
            chains.update(ids)
        if (len({r['fold_id'] for r in group}) != 1 or
                any(not re.fullmatch('[0-9a-f]{8}', r['fold_id']) for r in group)):
            raise ValueError(f'inconsistent cache at {key}')
    return groups


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
    columns, rows = read_cache(output, [c for c in COLUMNS if c != 'fold_locator'])
    for row in rows:
        row['fold_locator'] = fold_locator(row['source_path'], row['fold_path'])
    rows = [{column: row[column] for column in COLUMNS} for row in rows]
    cached = validate_cache(rows)
    if columns is not None and columns != COLUMNS:
        write_cache(output, COLUMNS, rows)
    processed = skipped = found = 0
    for fold in discover(input_path, error):
        found += 1
        key = (fold.source_path, fold.fold_path)
        if key in cached:
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
        write_cache(output, COLUMNS, [r for group in cached.values() for r in group])
        processed += 1
    if not found and not errors:
        error(f'{input_path}: no folds found')
    print(f'Folds: {found}; written: {processed}; cached: {skipped}; errors: {len(errors)}', file=stderr)
    return 1 if errors else 0
