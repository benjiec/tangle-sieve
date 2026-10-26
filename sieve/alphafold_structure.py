"""Shared input/identity contract for AlphaFold folds and derived assemblies."""
import hashlib
import io
import json
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from sieve.alphafold_files import Fold, model_number, open_fold, read_seed, request_name

ASSEMBLY_FORMAT = 'sieve-assembly-v1'
COMPLEX_FILE = 'model0_complex.cif'


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True, allow_nan=False)


def digest(text):
    return hashlib.sha256(text.encode('utf-8')).hexdigest()


def assembly_id(identity):
    return hashlib.sha1(canonical([ASSEMBLY_FORMAT, identity]).encode()).hexdigest()[-8:]


@contextmanager
def open_structure_fold(locator):
    locator = str(locator)
    if '::' in locator:
        with open_fold(locator) as fold:
            yield fold
        return
    path = Path(locator).resolve()
    if path.is_file() and path.suffix.lower() == '.cif':
        directory = path.parent
        models = [path.name]
    elif path.is_dir():
        directory = path
        models = sorted(p.name for p in path.iterdir() if p.is_file() and model_number(p.name) is not None)
    else:
        raise ValueError(f'expected a CIF, extracted fold directory or ZIP fold locator: {locator}')
    if not models:
        raise ValueError(f'no model CIF files at fold locator: {locator}')
    names = tuple(p.name for p in directory.iterdir() if p.is_file())
    yield Fold(str(path), '.', models, lambda n: (directory/n).read_text(encoding='utf-8'), names)


def assembly_metadata(fold):
    """Validate assemblies before use; unrecognized provenance never means AF input."""
    complex_models = [n for n in fold.models if PurePosixPath(n).name == COMPLEX_FILE]
    if not complex_models:
        return None
    if len(complex_models) != 1 or len(fold.models) != 1:
        raise ValueError('assembled directory must contain exactly one model0_complex.cif')
    name = complex_models[0]
    metadata_name = str(PurePosixPath(name).with_name('provenance.json'))
    if fold.member_names.count(metadata_name) != 1:
        raise ValueError('assembled complex requires exactly one provenance.json')
    try:
        meta = json.loads(fold.read(metadata_name))
        if meta['format'] != ASSEMBLY_FORMAT or meta['source_type'] != 'assembled':
            raise ValueError('unsupported assembly provenance format')
        if meta['fold_id'] != assembly_id(meta['identity']):
            raise ValueError('assembly identity hash mismatch')
        text = fold.read(name)
        if meta['cif_sha256'] != digest(text):
            raise ValueError('assembled CIF content hash mismatch; regenerate after coordinate edits')
        from sieve.alphafold_manifest import chain_sequences, sequence_id
        sequences = chain_sequences(text)
        if sequence_id(sequences.values()) != meta['sequence_id']:
            raise ValueError('assembly sequence ID mismatch')
        mappings = meta['chains']
        if set(mappings) != set(sequences) or any(mappings[c]['sequence'] != s for c,s in sequences.items()):
            raise ValueError('assembly chain mapping disagrees with CIF sequences')
        if meta['identity']['output_sha256'] != meta['cif_sha256']:
            raise ValueError('assembly identity does not match output coordinates')
        return meta
    except (KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError(f'invalid assembly provenance: {exc}') from exc


def identity_for_fold(fold, content, *, require_seed=True):
    from sieve.alphafold_manifest import chain_sequences, sequence_id, fold_id
    meta = assembly_metadata(fold)
    if meta:
        return dict(fold_id=meta['fold_id'], sequence_id=meta['sequence_id'], source_type='assembled', seed=None)
    sequences = chain_sequences(content)
    seed = None
    if require_seed or any(request_name(n) != n and request_name(n) in fold.member_names for n in fold.models):
        seed = read_seed(fold.read, fold.models, fold.member_names)
    return dict(fold_id=fold_id(seed, sequences.values()) if seed is not None else None,
                sequence_id=sequence_id(sequences.values()), source_type='alphafold' if seed is not None else 'cif', seed=seed)


@dataclass
class StructureInput:
    locator: str
    number: int
    content: str
    metadata: object
    identity: dict
    sequences: dict
    atoms: dict
    monomers: dict
    residue_origins: dict
    plddt_known: dict


def load_structure(locator, number=0):
    """Read exact polymer metadata and normalize coordinate residue keys to sequence positions."""
    import math
    from sieve.alphafold_manifest import chain_sequences
    with open_structure_fold(locator) as fold:
        numbered = [(model_number(n) if model_number(n) is not None else 0, n) for n in fold.models]
        names = [n for i,n in numbered if i == number]
        if len(names) != 1:
            raise ValueError(f'expected exactly one model {number}: {locator}')
        content = fold.read(names[0])
        metadata = assembly_metadata(fold)
        identity = identity_for_fold(fold, content, require_seed=False)
        af = metadata is None and any(request_name(n) != n and request_name(n) in fold.member_names for n in fold.models)
    sequences = chain_sequences(content)
    data = MMCIF2Dict(io.StringIO(content))
    entities = dict(zip(data['_struct_asym.id'], data['_struct_asym.entity_id']))
    polymer = {}
    for e,n,r in zip(data['_entity_poly_seq.entity_id'],data['_entity_poly_seq.num'],data['_entity_poly_seq.mon_id']):
        polymer.setdefault(e,{})[int(n)] = r
    monomers = {c:polymer[e] for c,e in entities.items()}
    count = len(data.get('_atom_site.id', []))
    if not count:
        raise ValueError('CIF has no atom coordinates')
    def column(name, default=None):
        values = data.get('_atom_site.'+name)
        if values is None:
            if default is None:
                raise ValueError(f'CIF missing atom field {name}')
            return [default]*count
        if len(values) != count:
            raise ValueError(f'inconsistent atom field {name}')
        return values
    columns = {k:column(k) for k in ('label_asym_id','label_seq_id','label_atom_id','label_comp_id','type_symbol','Cartn_x','Cartn_y','Cartn_z')}
    for key,default in [('label_alt_id','.'),('pdbx_PDB_model_num','1'),('B_iso_or_equiv','0'),('occupancy','1'),('auth_seq_id','?'),('pdbx_PDB_ins_code','?')]:
        columns[key] = column(key,default)
    if len(set(columns['pdbx_PDB_model_num'])) != 1:
        raise ValueError('CIF must contain exactly one coordinate model')
    atoms = {c:[] for c in sequences};origins={c:{} for c in sequences};seen=set()
    for i in range(count):
        r = {k:v[i] for k,v in columns.items()}
        c=r['label_asym_id'];n=int(r['label_seq_id']);name=r['label_atom_id']
        if c not in sequences or n not in monomers[c] or monomers[c][n] != r['label_comp_id']:
            raise ValueError('atom coordinates disagree with declared polymer sequence')
        if r['label_alt_id'] not in ('.','?') or (c,n,name) in seen:
            raise ValueError('alternate or duplicate atom coordinates are unsupported')
        seen.add((c,n,name))
        coord=[float(r['Cartn_'+axis]) for axis in 'xyz']
        b=float(r['B_iso_or_equiv']);occupancy=float(r['occupancy'])
        if not all(math.isfinite(x) for x in [*coord,b,occupancy]):
            raise ValueError('nonfinite coordinates or atom properties')
        atoms[c].append(dict(number=n,name=name,resname=r['label_comp_id'],element=r['type_symbol'].upper(),coord=coord,bfactor=b,occupancy=occupancy))
        origins[c][str(n)] = dict(auth_seq_id=r['auth_seq_id'],insertion_code=r['pdbx_PDB_ins_code'])
    known={c: bool(metadata['chains'][c].get('plddt_known',False)) if metadata else af for c in sequences}
    for c, values in atoms.items():
        if known[c] and any(not 0 <= a['bfactor'] <= 100 for a in values):
            raise ValueError(f'chain {c}: source pLDDT must be in [0, 100]')
    return StructureInput(str(locator),number,content,metadata,identity,sequences,atoms,monomers,origins,known)
