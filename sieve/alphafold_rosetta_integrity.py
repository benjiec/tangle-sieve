"""Host-side AlphaFold locator glue for the standalone PyRosetta integrity worker."""
import csv
import hashlib
import io
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio.PDB import MMCIFParser
from sieve.alphafold_files import model_number, read_seed
from sieve.alphafold_manifest import chain_sequences, fold_id, sequence_id
from sieve.alphafold_results import open_result_fold
from sieve.alphafold_structure import identity_for_fold
from sieve.alphafold_rosetta_interface import IMAGE, inputs
from sieve.rosetta_integrity_worker import MEASURES, PROTOCOL, check_finite

COLUMNS = ('fold_id', 'sequence_id', 'model_id', 'fold_locator', 'scope', 'chain',
           'score_function') + MEASURES + ('nonpolar_sasa_in_complex', 'nonpolar_shielded_by_partners')
WORKER = Path(__file__).with_name('rosetta_integrity_worker.py')


def prepare(locator, number, directory):
    with open_result_fold(locator) as fold:
        names = [n for n in fold.models if model_number(n) == number]
        if len(names) != 1:
            raise ValueError(f'expected exactly one model {number}: {locator}')
        content = fold.read(names[0])
        sequences = chain_sequences(content)
        metadata = dict(fold_id=identity_for_fold(fold, content)['fold_id'],
                        sequence_id=sequence_id(sequences.values()), model_id=number,
                        fold_locator=locator, score_function='ref2015')
    structure = MMCIFParser(QUIET=True).get_structure('model', io.StringIO(content))
    if len(structure) != 1:
        raise ValueError('CIF must contain exactly one coordinate model')
    expected = []
    for chain in structure[0]:
        if len(chain.id) != 1 or not chain.id.isascii() or not chain.id.isalnum():
            raise ValueError('Rosetta requires single-character alphanumeric chain IDs')
        for residue in chain:
            if residue.id[0] != ' ' or residue.is_disordered():
                raise ValueError('nonprotein residues or alternate conformations are unsupported')
            expected.append(dict(chain=chain.id, number=residue.id[1],
                                 icode=residue.id[2].strip(), name=residue.resname,
                                 atoms={a.name: [float(v) for v in a.coord] for a in residue
                                        if a.element not in ('H', 'D')}))
    if not expected:
        raise ValueError('model contains no residues')
    directory.mkdir()
    (directory / 'model.cif').write_text(content)
    (directory / 'expected.json').write_text(json.dumps(expected))
    shutil.copyfile(WORKER, directory / 'worker.py')
    return metadata


def validate_result(result):
    check_finite(result)
    rows = result.get('summaries')
    if not isinstance(rows, list) or len(rows) < 2:
        raise ValueError('integrity worker must report complex and chain summaries')
    if rows[0].get('scope') != 'complex' or rows[0].get('chain') != 'ALL':
        raise ValueError('missing complex integrity summary')
    chains = set()
    for index, row in enumerate(rows):
        for key in MEASURES + ('nonpolar_sasa_in_complex', 'nonpolar_shielded_by_partners'):
            if key not in row or isinstance(row[key], bool) or not isinstance(row[key], (int, float)):
                raise ValueError(f'missing or invalid integrity metric: {key}')
        if row['nres'] <= 0:
            raise ValueError('invalid residue count')
        if index:
            chain = row.get('chain')
            if row.get('scope') != 'chain' or not chain or chain in chains:
                raise ValueError('invalid or duplicate chain summary')
            chains.add(chain)
    if sum(row['nres'] for row in rows[1:]) != rows[0]['nres']:
        raise ValueError('complex and isolated residue counts disagree')
    return rows


def execute(locator, number, directory, image):
    metadata = prepare(locator, number, directory)
    command = ['docker', 'run', '--rm', '--network', 'none', '--mount',
               f'type=bind,source={directory},target=/work', '-w', '/work',
               image, 'python3', 'worker.py']
    settings = dict(metadata, protocol=PROTOCOL, image=image, command=command,
                    input_sha256=hashlib.sha256((directory/'model.cif').read_bytes()).hexdigest(),
                    worker_sha256=hashlib.sha256((directory/'worker.py').read_bytes()).hexdigest(),
                    status='running')
    path = directory / 'settings.json'
    path.write_text(json.dumps(settings, indent=2)+'\n')
    try:
        with (directory/'run.log').open('w') as log:
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        settings['returncode'] = completed.returncode
        if completed.returncode:
            raise ValueError(f'integrity worker exited {completed.returncode}; see {directory / "run.log"}')
        result = json.loads((directory/'integrity.json').read_text())
        rows = validate_result(result)
        settings.update(status='complete', pyrosetta_version=result['pyrosetta_version'])
        return [dict(metadata, **row) for row in rows]
    except Exception as exc:
        settings.update(status='failed', error=str(exc))
        raise
    finally:
        path.write_text(json.dumps(settings, indent=2)+'\n')


def run(source, stream, output, model, tsv, stdout):
    records = inputs(source, stream, model)
    if output is None:
        with tempfile.TemporaryDirectory(prefix='alphafold-integrity-') as name:
            try:
                return run_records(records, Path(name).resolve(), tsv, stdout)
            except (OSError, ValueError, KeyError, RuntimeError) as exc:
                logs = []
                for path in sorted(Path(name).rglob('*.log')):
                    tail = path.read_text(errors='replace')[-8000:].strip()
                    if tail:
                        logs.append(f'{path.relative_to(name)} (last 8000 characters):\n{tail}')
                raise ValueError(f'{exc}\nTemporary files removed.\n' + '\n'.join(logs)) from exc
    return run_records(records, Path(output).expanduser().resolve(), tsv, stdout)


def run_records(records, directory, tsv, stdout):
    if ',' in str(directory):
        raise ValueError('output directory cannot contain a comma (Docker mount syntax)')
    if directory.exists() and (not directory.is_dir() or any(directory.iterdir())):
        raise ValueError('output must be a new or empty directory; existing results are never overwritten')
    directory.mkdir(parents=True, exist_ok=True)
    rows = []
    if records:
        inspected = subprocess.run(['docker', 'image', 'inspect', IMAGE], capture_output=True, text=True)
        (directory/'image-inspect.log').write_text(inspected.stderr)
        if inspected.returncode:
            raise ValueError(f'cannot inspect {IMAGE}; see {directory / "image-inspect.log"}')
        info = json.loads(inspected.stdout)[0]
        (directory/'image.json').write_text(json.dumps(info, indent=2)+'\n')
        for index, (locator, number) in enumerate(records, 1):
            rows.extend(execute(locator, number, directory/f'{index:04d}-model-{number}', info['Id']))
    with (directory/'summary.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, COLUMNS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    if tsv:
        stdout.write((directory/'summary.tsv').read_text())
    else:
        for index, row in enumerate(rows):
            if index:
                print(file=stdout)
            for key in COLUMNS:
                print(f'{key}: {row[key]}', file=stdout)
