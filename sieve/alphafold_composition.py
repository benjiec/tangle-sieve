"""Exact base sequences plus per-copy PTMs for versioned AlphaFold identities."""
import io
import re
from collections import Counter, defaultdict
from dataclasses import dataclass

from Bio.Data.PDBData import protein_letters_3to1, protein_letters_3to1_extended
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

STANDARD = protein_letters_3to1
AMINO_ACIDS = frozenset(STANDARD.values())


@dataclass(frozen=True, order=True)
class ChainComposition:
    sequence: str
    modifications: tuple = ()


def ccd_id(value):
    if not isinstance(value, str) or not value or value != value.strip():
        raise ValueError(f'invalid PTM CCD identifier: {value!r}')
    value = value.removeprefix('CCD_')
    if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_]*', value):
        raise ValueError(f'invalid PTM CCD identifier: {value!r}')
    return value


def modifications(values, sequence):
    if not isinstance(values, list):
        raise ValueError('protein modifications must be a list')
    result = {}
    for value in values:
        if not isinstance(value, dict):
            raise ValueError('PTM must specify ptmPosition and ptmType')
        position = value.get('ptmPosition')
        if isinstance(position, bool) or not isinstance(position, int) or not 1 <= position <= len(sequence):
            raise ValueError(f'PTM position out of range or not an integer: {position!r}')
        if position in result:
            raise ValueError(f'duplicate PTM position: {position}')
        result[position] = ccd_id(value.get('ptmType'))
    return tuple(sorted(result.items()))


def hash_payload(chains):
    """Legacy serialization is unchanged unless at least one copy is modified."""
    chains = [ChainComposition(c) if isinstance(c, str) else c for c in chains]
    if not all(isinstance(c, ChainComposition) for c in chains):
        raise ValueError('identity requires base sequences or chain compositions')
    if not any(c.modifications for c in chains):
        return sorted(c.sequence for c in chains)
    records = []
    for chain in chains:
        mods = modifications([dict(ptmPosition=n, ptmType=code) for n,code in chain.modifications], chain.sequence)
        records.append((chain.sequence, mods))
    return ['ptm-v1', sorted(records)]


def request_compositions(request):
    """Return explicit chain IDs (AF3) or a multiset of copies (AF Server)."""
    if request is None or 'sequences' not in request:
        return None
    if not isinstance(request['sequences'], list) or not request['sequences']:
        raise ValueError('job-request sequences must be a nonempty list')
    by_id = {}; copies = []; dialects = set()
    for entry in request['sequences']:
        if not isinstance(entry, dict) or len(entry) != 1:
            raise ValueError('ambiguous job-request sequence entry')
        kind = next(iter(entry))
        if kind not in ('protein', 'proteinChain') or not isinstance(entry[kind], dict):
            raise ValueError('job request contains unsupported nonprotein sequence')
        dialects.add(kind)
        protein = entry[kind]; sequence = protein.get('sequence')
        if not isinstance(sequence, str) or not sequence or any(c not in AMINO_ACIDS for c in sequence):
            raise ValueError('job request requires an unambiguous standard protein base sequence')
        composition = ChainComposition(sequence, modifications(protein.get('modifications', []), sequence))
        if kind == 'proteinChain':
            count = protein.get('count', 1)
            if isinstance(count, bool) or not isinstance(count, int) or count < 1:
                raise ValueError('proteinChain count must be a positive integer')
            copies.extend([composition]*count)
        else:
            ids = protein.get('id')
            ids = [ids] if isinstance(ids, str) else ids
            if not isinstance(ids, list) or not ids:
                raise ValueError('protein requires chain id or list of ids')
            for chain in ids:
                if not isinstance(chain, str) or not chain or chain in by_id:
                    raise ValueError('invalid or duplicate job-request chain ID')
                by_id[chain] = composition
    if len(dialects) != 1:
        raise ValueError('mixed protein/proteinChain job-request formats')
    return by_id if by_id else copies


def chain_compositions(text, request=None):
    """Use declared polymers, checking PTMs against the request when available.

    AF Server requests lack chain IDs. Match chemically compatible compositions
    as a multiset, never by incidental CIF chain order. Ambiguity is an error.
    """
    cif = MMCIF2Dict(io.StringIO(text))
    def table(*columns):
        values = [cif.get(c, []) for c in columns]
        if not values[0] or len({len(v) for v in values}) != 1:
            raise ValueError('missing or inconsistent CIF fields: ' + ', '.join(columns))
        return zip(*values)
    types = dict(table('_entity.id', '_entity.type'))
    residues = defaultdict(dict)
    for entity, number, monomer in table('_entity_poly_seq.entity_id', '_entity_poly_seq.num', '_entity_poly_seq.mon_id'):
        number = int(number)
        if number <= 0 or number in residues[entity]:
            raise ValueError(f'ambiguous entity {entity} residue {number}')
        residues[entity][number] = monomer
    polymers = {}
    for chain,entity in table('_struct_asym.id', '_struct_asym.entity_id'):
        if chain in polymers:
            raise ValueError(f'duplicate CIF chain {chain}')
        if types.get(entity) != 'polymer' or entity not in residues:
            raise ValueError(f'chain {chain}: no identifiable protein sequence (entity {entity})')
        values = residues[entity]
        if sorted(values) != list(range(1,len(values)+1)):
            raise ValueError(f'chain {chain}: nonconsecutive entity sequence numbering')
        polymers[chain] = tuple(values[i] for i in sorted(values))
    # Explicit CIF parents take precedence over the bundled standard CCD mapping.
    parents = {}
    if '_chem_comp.mon_nstd_parent_comp_id' in cif:
        for code,parent in table('_chem_comp.id', '_chem_comp.mon_nstd_parent_comp_id'):
            if parent not in ('.','?'):
                if parent not in STANDARD:
                    raise ValueError(f'ambiguous/unsupported parent residue for {code}: {parent}')
                if code in parents and parents[code] != STANDARD[parent]:
                    raise ValueError(f'conflicting parent residue for {code}')
                parents[code] = STANDARD[parent]
    def base(code):
        value = STANDARD.get(code) or parents.get(code) or protein_letters_3to1_extended.get(code)
        return value if value in AMINO_ACIDS else None
    requested = request_compositions(request)
    def matches(monomers, composition):
        if len(monomers) != len(composition.sequence):
            return False
        mods = dict(composition.modifications)
        for n,(monomer,letter) in enumerate(zip(monomers,composition.sequence),1):
            if n in mods:
                if monomer != mods[n] or (base(monomer) is not None and base(monomer) != letter):
                    return False
            elif monomer not in STANDARD or STANDARD[monomer] != letter:
                return False
        return True
    result = {}
    if requested is not None:
        if isinstance(requested, dict):
            if set(requested) != set(polymers):
                raise ValueError('job-request chain IDs disagree with CIF')
            for chain,monomers in polymers.items():
                if not matches(monomers, requested[chain]):
                    raise ValueError(f'chain {chain}: job-request sequence/PTMs disagree with CIF')
                result[chain] = requested[chain]
        else:
            counts = Counter(requested)
            for chain,monomers in polymers.items():
                candidates = [c for c in counts if matches(monomers,c)]
                if len(candidates) != 1:
                    raise ValueError(f'chain {chain}: ambiguous or inconsistent job-request sequence/PTM mapping')
                result[chain] = candidates[0]
            if Counter(result.values()) != counts:
                raise ValueError('job-request copy counts disagree with CIF')
    else:
        for chain,monomers in polymers.items():
            letters=[];mods=[]
            for n,code in enumerate(monomers,1):
                letter=base(code)
                if letter is None:
                    raise ValueError(f'chain {chain}: unsupported residue {code} at {n}; cannot extract an exact protein sequence without parent metadata')
                letters.append(letter)
                if code not in STANDARD:
                    mods.append((n,code))
            result[chain] = ChainComposition(''.join(letters),tuple(mods))
    # Polymer coordinates, if present, must agree with the declared monomers.
    atom_columns = ('_atom_site.label_asym_id','_atom_site.label_seq_id','_atom_site.label_comp_id')
    if any(k in cif for k in atom_columns):
        for chain,position,code in table(*atom_columns):
            try:
                n=int(position)
            except ValueError as exc:
                raise ValueError('unmapped nonprotein atom in protein fold') from exc
            if chain not in polymers or not 1 <= n <= len(polymers[chain]) or polymers[chain][n-1] != code:
                raise ValueError('atom residue/PTM disagrees with declared polymer sequence')
    return result
