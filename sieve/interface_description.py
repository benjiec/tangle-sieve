"""Coordinate-based protein-interface descriptors; see docs/interface-comparison.md."""
from collections import Counter
from dataclasses import asdict, dataclass
import copy
import io
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from Bio.PDB import DSSP, MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1

from sieve.alphafold_pdockq2 import parse_model
from sieve.alphafold_pdockq2_patches import Options, detect_patches, input_models
from sieve.channels import KYTE_DOOLITTLE_HYDROPATHY, NEUTRAL_PH_CHARGE

POSITIVE = {'LYS': {'NZ'}, 'ARG': {'NE', 'NH1', 'NH2'}}
NEGATIVE = {'ASP': {'OD1', 'OD2'}, 'GLU': {'OE1', 'OE2'}}
VDW = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80}
HYDROPHOBIC = frozenset('AVILMFWYPC')
BACKBONE = frozenset(('N', 'CA', 'C', 'O', 'OXT'))


@dataclass(frozen=True)
class DescriptionOptions:
    cutoff: float = 8.0
    flank: int = 5
    radius: float = 8.0
    sasa_points: int = 240
    salt_cutoff: float = 4.0
    nonpolar_cutoff: float = 4.5
    clash_overlap: float = 0.6

    def __post_init__(self):
        for field in ('cutoff', 'radius', 'salt_cutoff', 'nonpolar_cutoff', 'clash_overlap'):
            if not math.isfinite(getattr(self, field)) or getattr(self, field) <= 0:
                raise ValueError(f'{field} must be finite and positive')
        if (not isinstance(self.flank, int) or not isinstance(self.sasa_points, int)
                or self.flank < 0 or self.sasa_points < 20):
            raise ValueError('flank must be nonnegative and sasa_points >= 20')


def protein_residues(chain):
    return [r for r in chain if 'CA' in r]


def residue_name(residue):
    aa = seq1(residue.resname)
    if aa not in KYTE_DOOLITTLE_HYDROPATHY:
        raise ValueError(f'unsupported residue {residue.resname} at {residue.id}')
    return aa


def composition(residues):
    residues = list(residues)
    sequence = ''.join(residue_name(r) for r in residues)
    count = len(sequence)
    return {'residues': [r.id[1] for r in residues], 'sequence': sequence, 'count': count,
            'mean_hydropathy': sum(KYTE_DOOLITTLE_HYDROPATHY[a] for a in sequence) / count if count else None,
            'hydrophobic_fraction': sum(a in HYDROPHOBIC for a in sequence) / count if count else None,
            'positive_count': sum(a in 'KR' for a in sequence), 'negative_count': sum(a in 'DE' for a in sequence),
            'histidine_count': sequence.count('H'),
            'nominal_sidechain_charge': sum(NEUTRAL_PH_CHARGE[a] for a in sequence)}


def secondary_structure(model, cif, executable):
    if executable is None:
        return {}, [], {'status': 'unavailable', 'reason': 'DSSP explicitly disabled'}
    version = subprocess.run([executable, '--version'], capture_output=True, text=True, check=True, timeout=30)
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / 'model.cif'
        path.write_text(cif)
        assigned = DSSP(model, str(path), dssp=executable, file_type='MMCIF')
    codes = {(chain, rid[1]): assigned[(chain, rid)][2] for chain, rid in assigned.keys()}
    segments = []
    families = {'H': 'alpha', 'G': '3_10', 'I': 'pi', 'E': 'beta', 'B': 'beta_bridge', 'P': 'polyproline'}
    for chain in model:
        counters = Counter()
        current = None
        for r in protein_residues(chain):
            code = codes.get((chain.id, r.id[1]), '?')
            family = families.get(code, 'loop' if code != '?' else 'unassigned')
            if current is None or current['code'] != code or r.id[1] != current['end'] + 1:
                counters[family] += 1
                label = f'{family}_{counters[family]}' if family == '3_10' else f'{family}{counters[family]}'
                current = {'chain': chain.id, 'label': label, 'code': code,
                           'start': r.id[1], 'end': r.id[1], 'residues': [r.id[1]]}
                segments.append(current)
            else:
                current['end'] = r.id[1]
                current['residues'].append(r.id[1])
    expected = {(c.id, r.id[1]) for c in model for r in protein_residues(c)}
    return codes, segments, {'status': 'assigned', 'version': (version.stdout or version.stderr).strip(),
                              'unassigned_residues': [f'{c}:{r}' for c, r in sorted(expected - set(codes))]}


def sasa_burial(model, chains, points):
    """Rigid separation SASA, all residues retained; no mutation or relaxation."""
    pair = copy.deepcopy(model)
    for chain in list(pair):
        if chain.id not in chains:
            pair.detach_child(chain.id)
    sr = ShrakeRupley(n_points=points, probe_radius=1.4)
    sr.compute(pair, level='A')
    bound = {(c.id, r.id, a.name): float(a.sasa) for c in pair for r in c for a in r}
    result = {}
    for chain in pair:
        sr.compute(chain, level='A')
        for r in protein_residues(chain):
            values = [(a, float(a.sasa) - bound[(chain.id, r.id, a.name)]) for a in r]
            result[(chain.id, r.id[1])] = {
                'sasa_isolated': sum(float(a.sasa) for a in r),
                'sasa_complex': sum(bound[(chain.id, r.id, a.name)] for a in r),
                'buried_area': sum(delta for _, delta in values),
                'buried_nonpolar_area': sum(delta for a, delta in values if a.element in ('C', 'S'))}
    return result


def atom_interactions(model, chains, options):
    """Explicit geometric candidates, not bond energies or hydrogen-bond assignments."""
    result = []
    first = [(r, a) for r in protein_residues(model[chains[0]]) for a in r if a.element != 'H']
    second = [(r, a) for r in protein_residues(model[chains[1]]) for a in r if a.element != 'H']
    xyz = np.array([a.coord for _, a in second])
    if not len(xyz):
        return result
    for r1, a in first:
        distances = np.linalg.norm(xyz - a.coord, axis=1)
        for j in np.flatnonzero(distances <= max(options.salt_cutoff, options.nonpolar_cutoff, 4.0)):
            r2, b = second[j]
            d = float(distances[j])
            types = []
            salt = ((a.name in POSITIVE.get(r1.resname, ()) and b.name in NEGATIVE.get(r2.resname, ())) or
                    (b.name in POSITIVE.get(r2.resname, ()) and a.name in NEGATIVE.get(r1.resname, ())))
            if salt and d <= options.salt_cutoff:
                types.append('salt_bridge_candidate')
            if (residue_name(r1) in HYDROPHOBIC and residue_name(r2) in HYDROPHOBIC
                    and a.name not in BACKBONE and b.name not in BACKBONE
                    and a.element in ('C', 'S') and b.element in ('C', 'S') and d <= options.nonpolar_cutoff):
                types.append('nonpolar_contact')
            overlap = VDW.get(a.element, 0) + VDW.get(b.element, 0) - d
            if overlap > options.clash_overlap and not (a.name == 'SG' and b.name == 'SG' and 1.8 <= d <= 2.2):
                types.append('steric_overlap_candidate')
            for kind in types:
                result.append({'chain_1': chains[0], 'residue_1': r1.id[1], 'aa_1': residue_name(r1), 'atom_1': a.name,
                               'chain_2': chains[1], 'residue_2': r2.id[1], 'aa_2': residue_name(r2), 'atom_2': b.name,
                               'type': kind, 'distance': d, 'vdw_overlap': overlap if kind == 'steric_overlap_candidate' else None})
    return result


def neighborhood(chain, numbers, flank, radius):
    residues = protein_residues(chain)
    chosen = [r for r in residues if r.id[1] in numbers]
    index = {r.id[1]: i for i, r in enumerate(residues)}
    before, after = set(), set()
    for r in chosen:
        i = index[r.id[1]]
        before.update(range(max(0, i-flank), i))
        after.update(range(i+1, min(len(residues), i+flank+1)))
    selected_indices = {index[n] for n in numbers}
    # Union around discontinuous segments; no duplicated residues within each set.
    before -= selected_indices
    after -= selected_indices
    spatial = []
    selected_atoms = np.array([a.coord for r in chosen for a in r if a.element != 'H'])
    for r in residues:
        if r.id[1] in numbers or not len(selected_atoms):
            continue
        if any(np.any(np.linalg.norm(selected_atoms-a.coord, axis=1) <= radius) for a in r if a.element != 'H'):
            spatial.append(r)
    return {'patch': composition(chosen),
            'upstream': composition(residues[i] for i in sorted(before)),
            'downstream': composition(residues[i] for i in sorted(after)),
            'sequence_shell': composition(residues[i] for i in sorted(before | after)),
            'spatial_shell': composition(spatial)}


def describe_model(cif, data, chains, options, dssp):
    parsed, pae = parse_model(cif, data)
    if len(set(chains)) != 2 or any(c not in parsed for c in chains):
        raise ValueError('select two distinct protein chains present in each model')
    model = MMCIFParser(QUIET=True).get_structure('model', io.StringIO(cif))[0]
    for c in model:
        seen = set()
        for r in protein_residues(c):
            if r.id[2] != ' ' or r.id[1] in seen:
                raise ValueError('insertion codes or duplicate residue numbers are not supported')
            seen.add(r.id[1])
            residue_name(r)
            if not 0 <= float(r['CA'].bfactor) <= 100:
                raise ValueError('pLDDT must be finite and in [0,100]')
            if any(not np.isfinite(a.coord).all() for a in r):
                raise ValueError('nonfinite atom coordinates')
    if any(not math.isfinite(float(v)) or float(v) < 0 for row in pae for v in row):
        raise ValueError('PAE must be finite and nonnegative')
    codes, segments, dssp_status = secondary_structure(model, cif, dssp)
    burial = sasa_burial(model, chains, options.sasa_points)
    interactions = atom_interactions(model, chains, options)
    residue_rows = []
    for c in chains:
        for r in protein_residues(model[c]):
            residue_rows.append({'chain': c, 'number': r.id[1], 'aa': residue_name(r),
                                 'ca': [float(v) for v in r['CA'].coord], 'plddt': float(r['CA'].bfactor),
                                 'secondary_structure': codes.get((c, r.id[1])),
                                 'missing_backbone_atoms': [a for a in ('N', 'CA', 'C', 'O') if a not in r],
                                 **burial[(c, r.id[1])]})
    contacts = []
    for a in parsed[chains[0]]:
        for b in parsed[chains[1]]:
            distance = math.dist(a.ca, b.ca)
            if distance <= options.cutoff:
                contacts.append({'residue_1': a.number, 'residue_2': b.number, 'distance': distance,
                                 'pae_1_to_2': float(pae[a.token_index][b.token_index]),
                                 'pae_2_to_1': float(pae[b.token_index][a.token_index])})
    patches = []
    for row in detect_patches(parsed, pae, Options(cutoff=options.cutoff), chains):
        selected = {row[f'chain {i}']: {int(n) for n in row[f'interface residues {i}'].split(',') if n} for i in (1, 2)}
        features = [x for x in interactions if x['residue_1'] in selected[chains[0]] and x['residue_2'] in selected[chains[1]]]
        typed_pairs = {(x['residue_1'], x['residue_2'], x['type']) for x in features}
        patches.append({'score': row,
                        'description': {c: neighborhood(model[c], selected[c], options.flank, options.radius) for c in chains},
                        'secondary_elements': [s for s in segments if s['chain'] in selected and selected[s['chain']].intersection(s['residues'])],
                        'buried_area_by_chain': {c: sum(burial[(c,n)]['buried_area'] for n in selected[c]) for c in chains},
                        'typed_residue_pair_counts': dict(Counter(x[2] for x in typed_pairs)),
                        'interactions': features})
    return {'chains': list(chains), 'residues': residue_rows, 'contacts': contacts, 'patches': patches,
            'interactions': interactions, 'secondary_elements': segments, 'dssp': dssp_status,
            'buried_area_by_chain': {c: sum(v['buried_area'] for (ch,_), v in burial.items() if ch == c) for c in chains},
            'buried_nonpolar_area_by_chain': {c: sum(v['buried_nonpolar_area'] for (ch,_), v in burial.items() if ch == c) for c in chains},
            'parameters': asdict(options)}


def describe_input(path, chains, options, dssp):
    models = []
    for source, number, cif, data in input_models(path):
        result = describe_model(cif, data, chains, options, dssp)
        result.update(source=source, model=number)
        models.append(result)
    if not models:
        raise ValueError('no models')
    for chain in chains:
        signature = [(r['number'], r['aa']) for r in models[0]['residues'] if r['chain'] == chain]
        if any([(r['number'],r['aa']) for r in m['residues'] if r['chain'] == chain] != signature for m in models):
            raise ValueError('chain sequences or residue numbering differ across models')
    return models
