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
from sieve.alphafold_files import discover, fold_locator, model_number, read_seed

BASE_COLUMNS = ['source_path', 'fold_path', 'fold_locator', 'fold_id', 'sequence_id', 'timestamp', 'source_type']
INPUT_COLUMNS = [f'input {i} {field}' for i in range(1, 6) for field in ('description', 'copies')]
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
    serialized = json.dumps(sorted(sequences), separators=(',', ':'), ensure_ascii=True)
    return hashlib.sha1(serialized.encode('utf-8')).hexdigest()[-8:]


def fold_id(seed, sequences):
    serialized = json.dumps((seed, sorted(sequences)), separators=(',', ':'), ensure_ascii=True)
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
    from sieve.alphafold_structure import identity_for_fold
    identity = identity_for_fold(fold, fold.read(fold.models[0]))
    return make_row(fold.source_path, fold.fold_path, identity['fold_id'], identity['sequence_id'],
                    fold.timestamp, [(candidates[seq], len(ids)) for seq, ids in sorted(grouped.items())],
                    source_type=identity['source_type'])


def make_row(source, member, identifier, sequence_identifier, timestamp, inputs, source_type='alphafold'):
    if len(inputs) > 5:
        raise ValueError(f'fold contains {len(inputs)} distinct sequences; maximum is 5')
    row = dict.fromkeys(COLUMNS, '')
    row.update(source_path=source, fold_path=member, fold_locator=fold_locator(source, member),
               fold_id=identifier, sequence_id=sequence_identifier, timestamp=timestamp, source_type=source_type)
    for index, (description, copies) in enumerate(inputs, 1):
        row[f'input {index} description'] = description
        row[f'input {index} copies'] = str(copies)
    return row


def validate_cache(rows):
    cached = {}
    for row in rows:
        if any(row.get(c) is None for c in COLUMNS):
            raise ValueError('cache contains an incomplete row')
        if row.get('source_type') not in ('alphafold', 'assembled'):
            raise ValueError('invalid manifest source_type')
        key = (row['source_path'], row['fold_path'])
        if key in cached:
            raise ValueError(f'duplicate cached fold: {key}')
        if any(not re.fullmatch('[0-9a-f]{8}', row[c]) for c in ('fold_id', 'sequence_id')):
            raise ValueError(f'invalid cached fold ID at {key}')
        ended = False
        count = 0
        for index in range(1, 6):
            description, copies = row[f'input {index} description'], row[f'input {index} copies']
            if description == '' and copies == '':
                ended = True
            elif ended or not description.strip() or not copies.isdigit() or int(copies) < 1:
                raise ValueError(f'invalid cached input {index} at {key}')
            else:
                count += 1
        if not count:
            raise ValueError(f'cached fold has no inputs: {key}')
        cached[key] = row
    return cached


def load_manifest(path):
    columns, rows = read_cache(path, ['source_path', 'fold_path', 'fold_id'])
    if columns is None:
        return {}, False
    if any(c not in columns for c in COLUMNS if c != 'source_type'):
        raise ValueError('manifest uses an older or incomplete format; remove it and regenerate with both fold_id and sequence_id')
    converted = [{c: row.get(c, 'alphafold' if c == 'source_type' else None) for c in COLUMNS} for row in rows]
    return validate_cache(converted), columns != COLUMNS


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
    cached, migrate = load_manifest(output)
    if migrate:
        write_cache(output, COLUMNS, cached.values())
    processed = skipped = found = 0
    for fold in discover(input_path, error):
        found += 1
        key = (fold.source_path, fold.fold_path)
        if key in cached:
            if not cached[key]['timestamp']:
                cached[key]['timestamp'] = fold.timestamp
                write_cache(output, COLUMNS, cached.values())
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
        write_cache(output, COLUMNS, cached.values())
        processed += 1
    if not found and not errors:
        error(f'{input_path}: no folds found')
    print(f'Folds: {found}; written: {processed}; cached: {skipped}; errors: {len(errors)}', file=stderr)
    return 1 if errors else 0
