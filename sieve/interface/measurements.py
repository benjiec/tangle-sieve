"""Coordinate measurements shared across interface profiles."""
from collections import defaultdict
from dataclasses import dataclass
import math

import numpy as np

from sieve.alphafold_pdockq2 import _distance, ca_distance_metrics, region_confidence, score_selections
from sieve.interface_description import BACKBONE, hydrogen_bond_roles, protein_residues

TIGHT_CONTACT_CUTOFF = 4.5

def angle(first, second):
    scale = np.linalg.norm(first) * np.linalg.norm(second)
    if scale == 0:
        return None
    return float(np.degrees(np.arccos(np.clip(np.dot(first, second) / scale, -1, 1))))


def donor_hydrogens(residue, donor, previous):
    explicit = [a.coord for a in residue if a.element in ('H', 'D') and a - donor <= 1.25]
    if explicit:
        return explicit
    neighbors = [a.coord for a in residue if a.element != 'H' and a.element != 'D'
                 and a != donor and 0.5 < a - donor <= 1.75]
    if (donor.name == 'N' and previous is not None and 'C' in previous
            and previous['C'] - donor <= 1.75):
        neighbors.append(previous['C'].coord)
    if donor.element != 'N' or not neighbors or len(neighbors) > 2:
        return []
    vectors = np.array(neighbors) - donor.coord
    vectors /= np.linalg.norm(vectors, axis=1)[:, None]
    if len(vectors) == 2:
        direction = -vectors.sum(axis=0)
        norm = np.linalg.norm(direction)
        return [donor.coord + direction / norm] if norm > 1e-6 else []
    # Terminal amide/guanidinium nitrogens lie in the plane of their parent carbon.
    if donor.name not in ('ND2', 'NE2', 'NH1', 'NH2'):
        return []
    parent = next((a for a in residue if a.element == 'C' and a - donor <= 1.75), None)
    if parent is None:
        return []
    axis = vectors[0]
    for atom in residue:
        if atom == donor or atom == parent or atom.element in ('H', 'D') or atom - parent > 1.75:
            continue
        offset = atom.coord - parent.coord
        perpendicular = offset - np.dot(offset, axis) * axis
        norm = np.linalg.norm(perpendicular)
        if norm > 1e-6:
            return [donor.coord - 0.5 * axis + sign * np.sqrt(0.75) * perpendicular / norm
                    for sign in (-1, 1)]
    return []


@dataclass(frozen=True)
class ContactEvidence:
    backbone_bonds: frozenset
    hydrogen_bonds: frozenset
    salt_pairs: frozenset


def qualify_contacts(model, chains, features):
    """Measure contacts within the selected residues and retain atom-level evidence."""
    by_number = {c: {r.id[1]: r for r in protein_residues(model[c])} for c in chains}
    previous = {}
    for residues in by_number.values():
        ordered = list(residues.values())
        previous.update({right: left for left, right in zip(ordered, ordered[1:])})
    pairs = defaultdict(set)
    bonds, backbone, uncertain, clashes = set(), set(), set(), set()
    for item in features:
        first, second = item['residue_1'], item['residue_2']
        atoms = item['atom_1'], item['atom_2']
        kind = item['type']
        key = (first, second, *atoms)
        pairs[kind].add((first, second))
        if kind == 'steric_overlap_candidate':
            clashes.add(key)
        if kind != 'hydrogen_bond_candidate':
            continue
        a = by_number[chains[0]][first].child_dict.get(atoms[0])
        b = by_number[chains[1]][second].child_dict.get(atoms[1])
        if a is None or b is None:
            uncertain.add(key)
            continue
        for donor, acceptor in ((a, b), (b, a)):
            source = donor.parent
            if not (hydrogen_bond_roles(source, donor)[0]
                    and hydrogen_bond_roles(acceptor.parent, acceptor)[1]):
                continue
            hydrogens = donor_hydrogens(source, donor, previous.get(source))
            if not hydrogens:
                uncertain.add(key)
            for h in hydrogens:
                theta = angle(donor.coord - h, acceptor.coord - h)
                if theta is not None and theta >= 120:
                    bonds.add(key)
                    if all(atom in BACKBONE for atom in atoms):
                        backbone.add(key)
    metrics = {
        'nonpolar residue pairs': len(pairs['nonpolar_contact']),
        'salt bridge residue pairs': len(pairs['salt_bridge_candidate']),
        'steric clash residue pairs': len({(a, b) for a, b, _, _ in clashes - bonds}),
        'backbone hydrogen bonds': len(backbone),
        'side-chain hydrogen bonds': len(bonds - backbone),
        'unresolved hydrogen-bond candidates': len(uncertain - bonds),
        'hydrogen-bond atom pairs': ','.join(f'{chains[0]}:{a}/{c}-{chains[1]}:{b}/{d}'
                                             for a, b, c, d in sorted(bonds)),
    }
    return metrics, ContactEvidence(frozenset(backbone), frozenset(bonds),
                                    frozenset(pairs['salt_bridge_candidate']))


def heavy_contacts(first, second, cutoff):
    """Minimum heavy-atom distance for each residue pair within the cutoff."""
    atoms_a = [(r.id[1], a) for r in first for a in r if a.element not in ('H', 'D')]
    atoms_b = [(r.id[1], a) for r in second for a in r if a.element not in ('H', 'D')]
    if not atoms_a or not atoms_b:
        return {}
    xyz = np.array([a.coord for _, a in atoms_b])
    pairs = {}
    for number, atom in atoms_a:
        distances = np.linalg.norm(xyz - atom.coord, axis=1)
        for index in np.flatnonzero(distances <= cutoff):
            key = (number, atoms_b[index][0])
            pairs[key] = min(pairs.get(key, math.inf), float(distances[index]))
    return pairs


def atom_model_ca_metrics(model, chains, contacts):
    first, second = [{r.id[1]: tuple(map(float, r['CA'].coord))
                      for r in protein_residues(model[c])} for c in chains]
    return ca_distance_metrics({(a, b): _distance(first[a], second[b]) for a, b in contacts})


def selected_confidence(chains, selections, parsed, cutoff=8.0):
    """Score a structural residue selection using the model's original PAE indices."""
    residues, pae = parsed
    selected = [[r for r in residues[c] if r.number in ids] for c, ids in zip(chains, selections)]
    metrics = region_confidence(*selected)
    for side, (selection, ids) in enumerate(zip(selected, selections), 1):
        if len(selection) != len(ids):
            metrics[f'region mean pLDDT {side}'] = None
    if any(metrics[key] is None for key in ('region mean pLDDT 1', 'region mean pLDDT 2')):
        pae = None
    scores = score_selections(*chains, *selected, pae, None, None, cutoff)
    metrics['pDockQ2 distance cutoff'] = cutoff
    metrics['pDockQ2 contact count'] = scores['contact count']
    if scores['contact count']:
        for key in ('pDockQ2 1 to 2', 'pDockQ2 2 to 1', 'pDockQ2 max',
                    'normalized PAE 1 to 2', 'normalized PAE 2 to 1',
                    'interface pLDDT 1', 'interface pLDDT 2'):
            metrics[key] = scores[key]
    return metrics
