"""Run containerized InterfaceAnalyzer on selected AlphaFold models."""
import csv
import io
import json
import math
import re
import subprocess
import tempfile
from pathlib import Path

from Bio.PDB import MMCIFParser

from sieve.alphafold_files import model_number, read_seed
from sieve.alphafold_manifest import chain_sequences, fold_id, sequence_id
from sieve.alphafold_result_format import PatchInput, read_patch_inputs
from sieve.alphafold_results import open_result_fold
from sieve.alphafold_structure import identity_for_fold

IMAGE = 'rosettacommons/rosetta:serial'
METRICS = ('dG_separated', 'dG_separated/dSASAx100', 'dSASA_int',
           'dSASA_hphobic', 'dSASA_polar', 'delta_unsatHbonds', 'hbonds_int',
           'sc_value', 'packstat', 'nres_int')
COLUMNS = ('fold_id', 'sequence_id', 'model_id', 'fold_locator', 'interface',
           'score_function', 'pack_input', 'pack_separated') + METRICS


def parse_interface(value):
    if not re.fullmatch(r'[A-Za-z0-9]:[A-Za-z0-9]', value) or value[0] == value[2]:
        raise ValueError('interface must name two distinct single-character chains, e.g. A:B')
    return value.split(':')


def inputs(source, stream, model):
    records = list(read_patch_inputs(stream)) if source == '-' else [PatchInput(source)]
    selected = []
    for record in records:
        if model is not None and record.model_id is not None and model != record.model_id:
            raise ValueError('--model conflicts with the model specified by stdin')
        number = record.model_id if record.model_id is not None else model
        if number is None or number < 0:
            raise ValueError('a plain locator requires a nonnegative --model')
        selected.append((record.locator, number))
    return selected


def read_scores(path):
    headers = None
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if not parts or parts[0] != 'SCORE:':
            continue
        if 'description' in parts:
            headers = parts[1:]
        else:
            if headers is None or len(parts) - 1 != len(headers):
                raise ValueError('malformed Rosetta scorefile')
            rows.append(dict(zip(headers, parts[1:])))
    if len(rows) != 1:
        raise ValueError('expected exactly one model in Rosetta scorefile')
    result = {}
    for key in METRICS:
        value = rows[0].get(key)
        if value is None or not math.isfinite(float(value)):
            raise ValueError(f'missing or invalid Rosetta metric: {key}')
        result[key] = value
    return result


def analyze(locator, number, interface, directory, image, repack=True):
    left, right = parse_interface(interface)
    with open_result_fold(locator) as fold:
        names = [n for n in fold.models if model_number(n) == number]
        if len(names) != 1:
            raise ValueError(f'expected exactly one model {number}: {locator}')
        content = fold.read(names[0])
        sequences = chain_sequences(content)
        identifier = identity_for_fold(fold, content)['fold_id']
    structure = MMCIFParser(QUIET=True).get_structure('model', io.StringIO(content))
    if len(structure) != 1:
        raise ValueError('CIF must contain exactly one coordinate model')
    chains = {c.id: len(list(c.get_residues())) for c in structure[0]}
    if left not in chains or right not in chains:
        raise ValueError(f'interface {interface} not present; coordinate chains: {list(chains)}')
    if any(len(c) != 1 or not c.isalnum() or not c.isascii() for c in chains):
        raise ValueError('Rosetta input requires single-character alphanumeric chain IDs')
    directory.mkdir()
    (directory / 'model.cif').write_text(content)
    row = dict(fold_id=identifier, sequence_id=sequence_id(sequences.values()),
               model_id=number, fold_locator=locator, interface=interface,
               score_function='ref2015', pack_input=repack, pack_separated=repack)
    command = ['docker', 'run', '--rm', '--network', 'none',
               '--mount', f'type=bind,source={directory},target=/work', '-w', '/work',
               image, 'InterfaceAnalyzer', '-s', 'model.cif',
               '-score:weights', 'ref2015', '-pack_input', str(repack).lower(),
               '-pack_separated', str(repack).lower(), '-compute_packstat', 'true',
               '-packstat:oversample', '100', '-compute_interface_sc', 'true',
               '-constant_seed', '-jran', '111111',
               '-out:file:score_only', 'interface.sc']
    # The executable ignores -interface for two-chain inputs and uses both chains.
    if len(chains) > 2:
        command += ['-interface', f'{left}_{right}']
    settings = dict(row, image=image, chain_residue_counts=chains, command=command,
                    status='running')
    settings_path = directory / 'settings.json'
    settings_path.write_text(json.dumps(settings, indent=2) + '\n')
    try:
        with (directory / 'run.log').open('w') as log:
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        settings['returncode'] = completed.returncode
        if completed.returncode:
            raise ValueError(f'InterfaceAnalyzer exited {completed.returncode}; see {directory / "run.log"}')
        row.update(read_scores(directory / 'interface.sc'))
        settings['status'] = 'complete'
    except Exception as exc:
        settings['status'] = 'failed'
        settings['error'] = str(exc)
        raise
    finally:
        settings_path.write_text(json.dumps(settings, indent=2) + '\n')
    return row


def run(source, stream, output, interface, model, tsv, stdout, repack=True):
    if output is not None:
        return _run(source, stream, output, interface, model, tsv, stdout, repack)
    with tempfile.TemporaryDirectory(prefix='alphafold-rosetta-') as temporary:
        try:
            return _run(source, stream, temporary, interface, model, tsv, stdout, repack)
        except (OSError, ValueError, KeyError, RuntimeError) as exc:
            # Preserve useful diagnostics before the temporary logs are deleted.
            details = []
            for log in sorted(Path(temporary).rglob('*.log')):
                tail = log.read_text(errors='replace')[-8000:].strip()
                if tail:
                    details.append(f'{log.relative_to(temporary)} (last 8000 characters):\n{tail}')
            message = f'{exc}\nTemporary working files removed.'
            if details:
                message += '\n' + '\n'.join(details)
            raise ValueError(message) from exc


def _run(source, stream, output, interface, model, tsv, stdout, repack=True):
    parse_interface(interface)
    records = inputs(source, stream, model)
    directory = Path(output).expanduser().resolve()
    if ',' in str(directory):
        raise ValueError('output directory cannot contain a comma (Docker mount syntax)')
    if directory.exists() and (not directory.is_dir() or any(directory.iterdir())):
        raise ValueError('output must be a new or empty directory; existing results are never overwritten')
    directory.mkdir(parents=True, exist_ok=True)
    if not records:
        if tsv:
            csv.writer(stdout, delimiter='\t', lineterminator='\n').writerow(COLUMNS)
        return
    inspected = subprocess.run(['docker', 'image', 'inspect', IMAGE],
                               capture_output=True, text=True)
    (directory / 'image-inspect.log').write_text(inspected.stderr)
    if inspected.returncode:
        raise ValueError(f'cannot inspect {IMAGE}; see {directory / "image-inspect.log"}')
    info = json.loads(inspected.stdout)[0]
    (directory / 'image.json').write_text(json.dumps(info, indent=2) + '\n')
    # Pin the local image ID so a moving tag cannot change between models.
    image = info['Id']
    rows = []
    for index, (locator, number) in enumerate(records, 1):
        rows.append(analyze(locator, number, interface,
                            directory / f'{index:04d}-model-{number}', image, repack))
    with (directory / 'summary.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, COLUMNS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    if tsv:
        writer = csv.DictWriter(stdout, COLUMNS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    else:
        for index, row in enumerate(rows):
            if index:
                print(file=stdout)
            for key in COLUMNS:
                print(f'{key}: {row[key]}', file=stdout)
