"""Evaluate TRAF MATH–SLiM contacts from model coordinates."""
from collections import defaultdict
from copy import deepcopy
from dataclasses import asdict, dataclass, replace
import io
from itertools import combinations
import json
import math

from Bio.PDB import MMCIFParser
import numpy as np

from sieve.alphafold_composition import chain_compositions
from sieve.alphafold_manifest import fold_id, sequence_id
from sieve.alphafold_patch_identity import attach_patch_identity, model_identity
from sieve.alphafold_pdockq2 import score_selections
from sieve.interface_description import (
    BACKBONE, HYDROPHOBIC, DescriptionOptions, atom_interactions, hydrogen_bond_roles,
    protein_residues, residue_label, sasa_burial,
)
from sieve.interface.traf_topology import complete_backbone, find_groove

TRAF_DOMAINS = frozenset(('MATH', 'TRAF', 'TRAF2', 'TRAF3', 'TRAF5', 'TRAF6'))
MOTIF_DOMAINS = frozenset(('SLiM', 'SLIM'))
COLUMNS = (
    'traf_paralog', 'traf_math_window', 'traf_math_source',
    'traf_motif_window', 'traf_motif_source', 'traf_motif_candidates',
    'traf_motif_register', 'traf_beta6', 'traf_beta7', 'traf_polar_loop',
    'traf_pocket', 'traf_pose_evidence',
    'traf_other_protomers', 'traf_status', 'traf_feasibility_rule',
    'traf_math_buried_area', 'traf_motif_buried_area', 'traf_tight_motif_residues',
    'traf_backbone_hbonds', 'traf_sidechain_hbonds', 'traf_hbond_atoms',
    'traf_unresolved_hbond_pairs',
    'traf_nonpolar_pairs', 'traf_clash_pairs', 'traf_salt_pairs',
    'traf_scoring_contacts', 'traf_motif_plddt', 'traf_thresholds',
)


@dataclass(frozen=True)
class Thresholds:
    mixed_area: float = 250
    packing_area: float = 400
    paired_area: float = 150


def parse_thresholds(values):
    thresholds = {}
    for item in values:
        key, sep, raw = item.partition('=')
        if not sep or key not in asdict(Thresholds()) or key in thresholds:
            raise ValueError(f'invalid or duplicate TRAF threshold {item!r}')
        value = float(raw)
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f'invalid TRAF threshold {item!r}')
        thresholds[key] = value
    return replace(Thresholds(), **thresholds)


def contact_pairs(first, second, cutoff):
    atoms_a = [(residue, atom) for residue in first for atom in residue if atom.element not in ('H', 'D')]
    atoms_b = [(residue, atom) for residue in second for atom in residue if atom.element not in ('H', 'D')]
    if not atoms_a or not atoms_b:
        return set()
    xyz = np.array([atom.coord for _, atom in atoms_b])
    return {(residue.id[1], atoms_b[j][0].id[1])
            for residue, atom in atoms_a
            for j in np.flatnonzero(np.linalg.norm(xyz - atom.coord, axis=1) <= cutoff)}


def angle(first, second):
    scale = np.linalg.norm(first) * np.linalg.norm(second)
    if scale == 0:
        return None
    return float(np.degrees(np.arccos(np.clip(np.dot(first, second) / scale, -1, 1))))


def residue_span(residues):
    return f'{residues[0].id[1]}-{residues[-1].id[1]}' if residues else None


def assess_pose(groove, motif, interactions):
    if groove is None or not (groove.loop and groove.pocket):
        return 'unresolved', None, 'groove landmarks unresolved'
    backbone, polar = interactions
    beta7 = {r.id[1]: r for r in groove.beta7}
    peptide = {r.id[1]: r for r in motif}
    ladder = [(a, b) for a, b, _, _ in backbone if a in beta7 and b in peptide]
    extended_pairing = set()
    for group in combinations(ladder, 3):
        ordered = sorted(set(group))
        first, last = ordered[0], ordered[-1]
        # An extended register spans four sequence steps on both edges.
        # A double bond at one residue pair alone adds no register span.
        if (last[0] - first[0] < 4 or not 4 <= first[1] - last[1] <= 5
                or any(a[1] < b[1] for a, b in zip(ordered, ordered[1:]) if a[0] != b[0])):
            continue
        beta_direction = beta7[last[0]]['CA'].coord - beta7[first[0]]['CA'].coord
        motif_direction = peptide[first[1]]['CA'].coord - peptide[last[1]]['CA'].coord
        strand_axis = groove.beta7[-1]['CA'].coord - groove.beta7[0]['CA'].coord
        if np.dot(beta_direction, motif_direction) < 0 and np.dot(strand_axis, motif_direction) < 0:
            extended_pairing.update(group)

    loop = {r.id[1] for r in groove.loop}
    polar_anchors = {b for a, b in polar if a in loop}
    # Side-chain contacts can remain plausible when the hydrogen orientation
    # cannot be resolved. Require a cluster, and keep it out of the H-bond count.
    for residue in motif:
        if len(polar_loop_contacts(residue, groove.loop)) >= 2:
            polar_anchors.add(residue.id[1])
    pocket_anchors = set()
    for r in motif:
        packing = contact_pairs_nonpolar(r, groove.pocket)
        if len(packing) >= 2:
            pocket_anchors.update((r.id[1], number) for number in packing)
    pocket_registers = {(start, stop) for start, _ in pocket_anchors for stop in polar_anchors
                        if 1 <= stop - start <= 4}
    evidence = []
    if extended_pairing:
        evidence.append('beta7=' + ','.join(f'{a}:{b}' for a, b in sorted(extended_pairing)))
    if pocket_registers:
        evidence.append('pocket/loop=' + ','.join(f'{a}:{b}' for a, b in sorted(pocket_registers)))
    if extended_pairing and pocket_registers:
        return 'unresolved', None, 'competing signatures ' + '; '.join(evidence)
    if extended_pairing:
        numbers = sorted({b for _, b in extended_pairing})
        if numbers[-1] - numbers[0] > 5:
            return 'unresolved', None, 'competing beta7 registers ' + evidence[0]
        return 'traf6', f'{numbers[0]}-{numbers[-1]}', evidence[0]
    if len(pocket_registers) == 1:
        start, stop = next(iter(pocket_registers))
        return 'traf2/3', f'{start}-{stop}', evidence[0]
    if pocket_registers:
        return 'unresolved', None, 'competing pocket registers ' + evidence[0]
    return 'atypical', None, 'no matching pose signature'


def contact_pairs_nonpolar(peptide, pocket):
    if residue_label(peptide) not in HYDROPHOBIC | {'S', 'T'}:
        return set()
    atoms = [a for a in peptide if a.name not in BACKBONE and a.element in ('C', 'S')]
    return {residue.id[1] for residue in pocket
            if any(a - b <= 4.5 for a in atoms for b in residue
                   if b.element in ('C', 'S'))}


def polar_loop_contacts(peptide, loop):
    return {
        residue.id[1] for residue in loop
        if any(a - b <= 3.8 and (donor_a and acceptor_b or donor_b and acceptor_a)
               for a in peptide if a.name not in BACKBONE and a.element in ('N', 'O')
               for donor_a, acceptor_a in (hydrogen_bond_roles(peptide, a),)
               for b in residue if b.name not in BACKBONE and b.element in ('N', 'O')
               for donor_b, acceptor_b in (hydrogen_bond_roles(residue, b),))
    }


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


def chemical_metrics(model, chains, math_residues, motif, tight, thresholds):
    pair = deepcopy(model)
    selected = dict(zip(chains, (math_residues, motif)))
    previous = {}
    for chain in list(pair):
        if chain.id not in selected:
            pair.detach_child(chain.id)
            continue
        residues = protein_residues(chain)
        previous[chain.id] = {right.id: left for left, right in zip(residues, residues[1:])}
        keep = {r.id for r in selected[chain.id]}
        for residue in list(chain):
            if residue.id not in keep:
                chain.detach_child(residue.id)
            else:
                for atom in list(residue):
                    if atom.element in ('H', 'D'):
                        # SASA uses heavy atoms. H-bond geometry reads the original model.
                        residue.detach_child(atom.id)

    features = atom_interactions(pair, chains, DescriptionOptions(cutoff=6))
    areas = sasa_burial(pair, chains, 240)
    math_area, motif_area = (max(0., sum(item['buried_area'] for (cid, _), item in areas.items()
                                       if cid == chain)) for chain in chains)
    pairs = defaultdict(set)
    bonds, backbone, uncertain, clashes = set(), set(), set(), []
    for item in features:
        first, second = item['residue_1'], item['residue_2']
        atoms = item['atom_1'], item['atom_2']
        kind = item['type']
        key = (first, second, *atoms)
        if kind == 'steric_overlap_candidate':
            clashes.append(key)
        elif kind == 'hydrogen_bond_candidate':
            left, right = (model[c][n] for c, n in zip(chains, (first, second)))
            a, b = left[atoms[0]], right[atoms[1]]
            for donor, acceptor in ((a, b), (b, a)):
                source = donor.parent
                if not (hydrogen_bond_roles(source, donor)[0]
                        and hydrogen_bond_roles(acceptor.parent, acceptor)[1]):
                    continue
                hydrogens = donor_hydrogens(source, donor, previous[source.parent.id].get(source.id))
                if not hydrogens:
                    uncertain.add(key)
                for h in hydrogens:
                    theta = angle(donor.coord - h, acceptor.coord - h)
                    if theta is not None and theta >= 120:
                        bonds.add(key)
                        if all(atom in BACKBONE for atom in atoms):
                            backbone.add(key)
        else:
            pairs[kind].add((first, second))
    clash_pairs = {(a, b) for a, b, c, d in clashes if (a, b, c, d) not in bonds}
    nonpolar = len(pairs['nonpolar_contact'])
    total_area = math_area + motif_area
    metrics = {
        'total buried surface area': total_area,
        'traf_math_buried_area': math_area,
        'traf_motif_buried_area': motif_area,
        'traf_tight_motif_residues': len(tight),
        'traf_backbone_hbonds': len(backbone),
        'traf_sidechain_hbonds': len(bonds - backbone),
        'traf_hbond_atoms': ','.join(f'{chains[0]}:{a}/{c}-{chains[1]}:{b}/{d}'
                                   for a, b, c, d in sorted(bonds)),
        'traf_unresolved_hbond_pairs': len(uncertain - bonds),
        'traf_nonpolar_pairs': nonpolar,
        'traf_clash_pairs': len(clash_pairs),
        'traf_salt_pairs': len(pairs['salt_bridge_candidate']),
    }
    span = min(len({bond[i] for bond in backbone}) for i in (0, 1))
    metrics['traf_feasibility_rule'] = feasibility_rule(metrics, span, thresholds)
    polar = {(a, b) for a, b, c, _ in bonds if c not in BACKBONE}
    polar.update(pairs['salt_bridge_candidate'])
    return metrics, (backbone, polar)


def feasibility_rule(metrics, span, limits):
    tight = metrics['traf_tight_motif_residues']
    clashes = metrics['traf_clash_pairs']
    nonpolar = metrics['traf_nonpolar_pairs']
    backbone = metrics['traf_backbone_hbonds']
    bonds = backbone + metrics['traf_sidechain_hbonds']
    area = metrics['total buried surface area']
    if clashes >= 3 or tight < 2:
        return 'clash_or_sparse_contact'
    if area >= limits.paired_area and backbone >= 3 and span >= 2 and clashes == 0:
        return 'small_but_well_paired'
    if area >= limits.mixed_area and bonds >= 2 and nonpolar >= 2:
        return 'hydrogen_bond_supported'
    if area >= limits.packing_area and tight >= 3 and nonpolar >= 5:
        return 'packing_supported_provisional'
    if metrics['traf_unresolved_hbond_pairs'] and area >= limits.mixed_area and nonpolar >= 2:
        return 'unresolved'
    return 'insufficient_support'


def confidence(row):
    return row.get('traf_confidence')


def feasibility(row):
    rule = row.get('traf_feasibility_rule')
    if rule in (None, 'unresolved'):
        return None
    return rule in ('small_but_well_paired', 'hydrogen_bond_supported', 'packing_supported_provisional')


@dataclass
class Evaluation:
    rows: list
    fold_ids: dict
    sequence_ids: dict
    descriptions: None = None

    def leaves(self):
        yield from self.rows


@dataclass
class Classification:
    interface: str
    sides: dict
    settings: dict


def classify_interface(row, structure, domains, *, chains):
    return Classification(row['traf_pose'], {chain: {} for chain in chains}, {})


def select_window(residues, window):
    start, end = window[0]
    selected = [r for r in residues if start <= r.id[1] <= end]
    return selected if len(selected) == end - start + 1 and complete_backbone(selected) else []


def candidate_windows(motif, contacts, tight):
    contacted = {b for _, b in contacts}
    packed = defaultdict(int)
    for _, number in tight:
        packed[number] += 1
    windows = []
    for size in (5, 6):
        for index in range(len(motif) - size + 1):
            window = motif[index:index + size]
            numbers = {r.id[1] for r in window}
            if contacted & numbers and complete_backbone(window):
                sites = packed.keys() & numbers
                density = sum(packed[n] for n in sites) / size
                windows.append((window, sites, density))
    if not windows and 1 <= len(motif) < 5 and complete_backbone(motif):
        return [(motif, set(packed), sum(packed.values()) / len(motif))]
    maximum = max((density for _, _, density in windows), default=0)
    return [item for item in windows if maximum - item[2] <= .25]


def interface_confidence(chains, math_residues, motif, parsed):
    residues, pae = parsed
    regions = [{r.id[1] for r in selection} for selection in (math_residues, motif)]
    selected = [[r for r in residues[c] if r.number in region]
                for c, region in zip(chains, regions)]
    values = [r.plddt for r in selected[1]]
    result = {}

    def valid(value):
        return value is not None and math.isfinite(value) and 0 <= value <= 100

    if values and len(values) == len(motif) and all(map(valid, values)):
        result['traf_motif_plddt'] = sum(values) / len(values)
    if pae is None or not all(valid(r.plddt) for selection in selected for r in selection):
        return result
    scores = score_selections(*chains, *selected, pae, None, None)
    result['traf_scoring_contacts'] = scores['contact count']
    if scores['contact count']:
        for key in ('pDockQ2 1 to 2', 'pDockQ2 2 to 1', 'pDockQ2 max',
                    'normalized PAE 1 to 2', 'normalized PAE 2 to 1',
                    'interface pLDDT 1', 'interface pLDDT 2'):
            result[key] = scores[key]
        result['traf_confidence'] = min(scores['pDockQ2 1 to 2'], scores['pDockQ2 2 to 1'])
    return result


@dataclass
class Candidate:
    window: list
    sites: set
    density: float
    metrics: dict
    interactions: tuple

    def ranking(self):
        m = self.metrics
        chemistry_known = m.get('traf_feasibility_rule') != 'unresolved' or 'traf_clash_pairs' in m
        return (self.density, chemistry_known, -m.get('traf_clash_pairs', 0),
                m.get('traf_backbone_hbonds', 0) + m.get('traf_sidechain_hbonds', 0),
                m.get('total buried surface area', 0) / len(self.window),
                m.get('traf_confidence', -1), m.get('traf_motif_plddt', -1))

    def comparable(self, other):
        first, second = self.ranking(), other.ranking()
        return (abs(first[0] - second[0]) <= .25 and first[1:4] == second[1:4]
                and abs(first[4] - second[4]) <= 10
                and (min(first[5], second[5]) < 0 or abs(first[5] - second[5]) <= .05)
                and (min(first[6], second[6]) < 0 or abs(first[6] - second[6]) <= 5))


def evaluate_sites(model, chains, domain, windows, thresholds, parsed):
    candidates = []
    for window, sites, density in windows:
        metrics = interface_confidence(chains, domain, window, parsed)
        contacts = contact_pairs(domain, window, 6)
        partners = {a for a, _ in contacts}
        if all(all(atom in r for atom in ('N', 'CA', 'C', 'O'))
               for r in domain if r.id[1] in partners):
            chemistry, interactions = chemical_metrics(model, chains, domain, window, sites, thresholds)
            metrics.update(chemistry)
        else:
            metrics['traf_feasibility_rule'] = 'unresolved'
            interactions = (set(), set())
        candidates.append(Candidate(window, sites, density, metrics, interactions))
    return sorted(candidates, key=Candidate.ranking, reverse=True)


def extract(locator, chains, *, record, seeds, model_scores, dssp, domains,
            motif_windows, math_windows, thresholds, parsed):
    source, number, cif, data = record
    model_id = '' if number is None else number
    residues, pae = parsed
    model = MMCIFParser(QUIET=True).get_structure('traf', io.StringIO(cif))[0]
    math_chain, peptide_chain = chains
    receptor = protein_residues(model[math_chain])
    peptide = protein_residues(model[peptide_chain])
    if any(not np.isfinite(a.coord).all() for r in (*receptor, *peptide) for a in r):
        raise ValueError('TRAF interface coordinates must be finite')
    explicit_math = math_windows.get(math_chain)
    explicit_motif = motif_windows.get(peptide_chain)
    if explicit_math:
        receptor = select_window(receptor, explicit_math)
    groove = find_groove(receptor, donor_hydrogens, terminal=not explicit_math)
    math_domain = receptor if explicit_math else list(groove.domain) if groove else []
    motif = select_window(peptide, explicit_motif) if explicit_motif else peptide

    seed = seeds.get((source, model_id))
    try:
        composition = list(chain_compositions(cif, data.get('_job_request')).values())
    except (ValueError, KeyError):
        if '_entity_poly_seq.' in cif or data.get('_job_request', {}).get('sequences'):
            raise
        composition = []
    sequence = sequence_id(composition) if composition else None
    fold = (data['_assembly']['fold_id'] if '_assembly' in data
            else fold_id(seed, composition) if seed is not None and composition else None)
    domain = domains[math_chain]
    row = {
        'source': source, 'model': model_id, 'fold_id': fold, 'sequence_id': sequence,
        'scope': 'patch', 'status': 'leaf', 'rank': 1, 'patch id': 'traf_slim',
        'chain 1': math_chain, 'chain 2': peptide_chain, 'distance cutoff': 6.0,
        'interface residues 1': '', 'interface residues 2': '',
        'traf_paralog': domain if domain.startswith('TRAF') and domain != 'TRAF' else None,
        'traf_math_window': residue_span(math_domain),
        'traf_math_source': 'supplied' if explicit_math else 'inferred',
        'traf_motif_source': 'supplied' if explicit_motif else 'inferred',
        'traf_motif_window': residue_span(motif) if explicit_motif else None,
        'traf_thresholds': json.dumps(asdict(thresholds), sort_keys=True),
        'traf_pose': 'unresolved',
        **model_scores.get((source, model_id), {}),
    }
    if groove:
        row.update(traf_beta6=residue_span(groove.beta6), traf_beta7=residue_span(groove.beta7),
                   traf_polar_loop=residue_span(groove.loop), traf_pocket=residue_span(groove.pocket))
    if not math_domain:
        row['traf_status'] = ('missing MATH window coordinates' if explicit_math
                              else 'complete terminal MATH topology unresolved')
    elif not motif:
        row['traf_status'] = 'missing motif coordinates'
    else:
        contacts = contact_pairs(math_domain, motif, 6.0)
        tight = contact_pairs(math_domain, motif, 4.5)
        if not contacts:
            row.update(traf_pose='no_contact', traf_status='no_contact',
                       traf_feasibility_rule='no_contact', **{'contact count': 0})
        else:
            windows = ([(motif, {b for _, b in tight}, len(tight) / len(motif))]
                       if explicit_motif else candidate_windows(motif, contacts, tight))
            if pae is not None and not all(math.isfinite(float(x)) and float(x) >= 0 for line in pae for x in line):
                parsed = residues, None
            candidates = evaluate_sites(model, chains, math_domain, windows, thresholds, parsed)
            if not candidates:
                row['traf_status'] = 'contacting motif lacks a contiguous backbone window'
            else:
                update_site(row, candidates, groove, contacts)
    attach_patch_identity(row, model_identity(cif, data, seed), number)
    return Evaluation([row], {source: fold}, {source: sequence})


def update_site(row, candidates, groove, contacts):
    best = candidates[0]
    motif_numbers = {r.id[1] for r in best.window}
    contacts = {(a, b) for a, b in contacts if b in motif_numbers}
    row['contact count'] = len(contacts)
    for side in (1, 2):
        numbers = sorted({pair[side - 1] for pair in contacts})
        row[f'interface residues {side}'] = ','.join(map(str, numbers))
        row[f'residue count {side}'] = len(numbers)
    row['region 1'] = row['traf_math_window']
    row['region 2'] = row['traf_motif_window'] = residue_span(best.window)
    row.update(best.metrics)
    if 'traf_clash_pairs' in best.metrics:
        pose, register, evidence = assess_pose(groove, best.window, best.interactions)
    else:
        pose, register, evidence = 'unresolved', None, 'incomplete contacting MATH backbone'
    competitors = []
    for candidate in candidates[1:]:
        if not best.comparable(candidate):
            continue
        alternate_pose, alternate_register, _ = assess_pose(groove, candidate.window, candidate.interactions)
        shared = len(best.sites & candidate.sites)
        same_site = shared >= max(1, min(len(best.sites), len(candidate.sites)) - 1)
        competing_register = (alternate_register is not None and register is not None
                              and (alternate_pose, alternate_register) != (pose, register))
        if not same_site or competing_register:
            competitors.append(residue_span(candidate.window))
    if competitors:
        row.update(traf_pose='unresolved', traf_status='competing inferred motif sites or registers',
                   traf_motif_candidates=','.join((residue_span(best.window), *competitors)))
    else:
        row.update(traf_pose=pose, traf_status='evaluated', traf_motif_register=register)
    row['traf_pose_evidence'] = evidence
