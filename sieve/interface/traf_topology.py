"""Infer the TRAF MATH strand topology from backbone coordinates."""
from dataclasses import dataclass
from itertools import combinations

import numpy as np


@dataclass(frozen=True)
class Groove:
    domain: tuple
    beta6: tuple
    beta7: tuple
    loop: tuple
    pocket: tuple


def strand_ladders(residues, donor_hydrogens):
    """Find antiparallel ladders of nonlocal backbone hydrogen bonds."""
    oxygen = np.array([r['O'].coord if 'O' in r else (np.inf,) * 3 for r in residues])
    pairs = set()
    for i, residue in enumerate(residues):
        if 'N' not in residue or residue.resname == 'PRO':
            continue
        n = residue['N']
        previous = residues[i - 1] if i else None
        hydrogens = donor_hydrogens(residue, n, previous)
        for j in np.flatnonzero(np.linalg.norm(oxygen - n.coord, axis=1) <= 3.5):
            if abs(i - j) < 8:
                continue
            if any(np.dot(n.coord - h, oxygen[j] - h) <=
                   -0.5 * np.linalg.norm(n.coord - h) * np.linalg.norm(oxygen[j] - h)
                   for h in hydrogens):
                pairs.add(tuple(sorted((i, int(j)))))

    # Allow a short interruption for a bulge. Nearby parallel contacts
    # cannot seed a ladder because they advance in the same direction.
    neighbors = {pair: set() for pair in pairs}
    for a, b in combinations(sorted(pairs), 2):
        di, dj = b[0] - a[0], a[1] - b[1]
        if 1 <= di <= 4 and 1 <= dj <= 4:
            neighbors[a].add(b)
            neighbors[b].add(a)
    ladders = []
    unseen = set(pairs)
    while unseen:
        queue = [min(unseen)]
        connected = set()
        while queue:
            pair = queue.pop()
            if pair in connected:
                continue
            connected.add(pair)
            queue.extend(neighbors[pair] - connected)
        unseen -= connected
        if len({i for i, _ in connected}) >= 2 and len({j for _, j in connected}) >= 2:
            ladders.append(tuple((min(positions), max(positions)) for positions in zip(*connected)))
    return sorted(ladders)


def near(first, second, margin=2):
    return max(first[0], second[0]) <= min(first[1], second[1]) + margin


def complete_backbone(residues):
    return bool(residues) and all(
        r.id[2] == ' ' and all(atom in r for atom in ('N', 'CA', 'C', 'O')) for r in residues
    ) and all(b.id[1] == a.id[1] + 1 for a, b in zip(residues, residues[1:]))


def find_groove(residues, donor_hydrogens, *, terminal):
    """Identify the two MATH sheets by their connectivity in sequence order."""
    ladders = strand_ladders(residues, donor_hydrogens) if residues else []
    candidates = set()
    for beta1, beta8 in ladders:
        if beta8[0] - beta1[1] < 75:
            continue
        if terminal and not terminal_cap(residues, beta8):
            continue
        for beta5, beta8_partner in ladders:
            if beta5[0] <= beta1[1] or beta5[1] >= beta8[0] or not near(beta8, beta8_partner):
                continue
            for beta5_partner, beta6 in ladders:
                if not near(beta5, beta5_partner, margin=3) or beta6[0] <= beta5[1]:
                    continue
                for beta4, beta7 in ladders:
                    if not beta1[1] < beta4[0] < beta5[0] < beta6[0] < beta7[0] < beta8[0]:
                        continue
                    if not any(
                        beta1[1] < beta2[0] < beta3[0] < beta4[0]
                        and near(beta3, beta3_partner) and near(beta4, beta4_partner)
                        for beta3, beta4_partner in ladders
                        for beta2, beta3_partner in ladders
                    ):
                        continue
                    candidates.add((beta1, beta5, beta6, beta7, beta8))
    if len(candidates) != 1:
        return None
    beta1, _, beta6, beta7, beta8 = next(iter(candidates))
    # Include β7 edge residues just upstream of the first observed ladder bond.
    beta7_start = max(beta6[1] + 1, beta7[0] - 2)
    start = beta1[0]
    domain = tuple(residues[start:] if terminal else residues)
    if terminal and not complete_backbone(domain):
        return None
    edge6 = tuple(residues[beta6[0]:beta6[1] + 1])
    edge7 = tuple(residues[beta7_start:beta7[1] + 1])
    loop = tuple(residues[beta6[1] + 1:beta7_start])
    # The TRAF2/3 pocket is formed at the β6 exit and β7 entry.
    pocket = (*edge6[-3:], *loop, *edge7[:9])
    return Groove(domain, edge6, edge7, loop, pocket)


def terminal_cap(residues, beta8):
    trailing = residues[beta8[1] + 1:]
    if len(trailing) > 10:
        return False
    terminal_strand = residues[beta8[0]:beta8[1] + 1]
    return all(min(r['CA'] - s['CA'] for s in terminal_strand) <= 22 for r in trailing)
