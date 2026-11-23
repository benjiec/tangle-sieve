"""Evaluate TRAF MATH–SLiM contacts from model coordinates."""
from collections import defaultdict
from dataclasses import asdict, dataclass, replace
from itertools import combinations
import json
import math

import numpy as np

from sieve.alphafold_composition import chain_compositions
from sieve.alphafold_manifest import fold_id, sequence_id
from sieve.alphafold_patch_identity import attach_patch_identity, model_identity
from sieve.interface_description import (
    BACKBONE, HYDROPHOBIC, hydrogen_bond_roles, protein_residues, residue_label,
)
from sieve.interface.measurements import (
    ContactEvidence, TIGHT_CONTACT_CUTOFF, atom_model_ca_metrics, donor_hydrogens,
    heavy_contacts, selected_confidence,
)
from sieve.interface.results import Classification, Evaluation
from sieve.patch_description_summary import PatchDescriptions
from sieve.interface.traf_topology import complete_backbone, find_groove

TRAF_PARALOGS = frozenset(('TRAF2', 'TRAF3', 'TRAF5', 'TRAF6'))
TRAF_DOMAINS = TRAF_PARALOGS | {'MATH', 'TRAF_MATH'}
MOTIF_DOMAINS = frozenset(('SLiM', 'SLIM'))
COLUMNS = (
    'traf_paralog', 'traf_math_window', 'traf_math_source',
    'traf_motif_window', 'traf_motif_source', 'traf_motif_candidates',
    'traf_motif_register', 'traf_beta6', 'traf_beta7', 'traf_polar_loop',
    'traf_pocket', 'traf_pose_evidence',
    'traf_other_protomers', 'traf_status', 'traf_feasibility_rule',
    'traf_thresholds',
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


def residue_span(residues):
    return f'{residues[0].id[1]}-{residues[-1].id[1]}' if residues else None


def assess_pose(groove, motif, interactions):
    if groove is None or not (groove.loop and groove.pocket):
        return 'unresolved', None, 'groove landmarks unresolved'
    backbone = interactions.backbone_bonds
    polar = {(a, b) for a, b, c, _ in interactions.hydrogen_bonds if c not in BACKBONE}
    polar.update(interactions.salt_pairs)
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


def feasibility_rule(metrics, span, limits):
    tight = metrics['tight-contact residues 2']
    clashes = metrics['steric clash residue pairs']
    nonpolar = metrics['nonpolar residue pairs']
    backbone = metrics['backbone hydrogen bonds']
    bonds = backbone + metrics['side-chain hydrogen bonds']
    area = metrics['total buried surface area']
    if clashes >= 3 or tight < 2:
        return 'clash_or_sparse_contact'
    if area >= limits.paired_area and backbone >= 3 and span >= 2 and clashes == 0:
        return 'small_but_well_paired'
    if area >= limits.mixed_area and bonds >= 2 and nonpolar >= 2:
        return 'hydrogen_bond_supported'
    if area >= limits.packing_area and tight >= 3 and nonpolar >= 5:
        return 'packing_supported_provisional'
    if metrics['unresolved hydrogen-bond candidates'] and area >= limits.mixed_area and nonpolar >= 2:
        return 'unresolved'
    return 'insufficient_support'


def confidence(row):
    return row.get('traf_confidence')


def feasibility(row):
    rule = row.get('traf_feasibility_rule')
    if rule in (None, 'unresolved'):
        return None
    return rule in ('small_but_well_paired', 'hydrogen_bond_supported', 'packing_supported_provisional')


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


@dataclass
class Candidate:
    window: list
    sites: set
    density: float
    metrics: dict
    evidence: ContactEvidence | None

    def ranking(self):
        m = self.metrics
        chemistry_known = self.evidence is not None
        confidence = m.get('traf_confidence')
        plddt = m.get('region mean pLDDT 2')
        return (self.density, chemistry_known,
                -m['steric clash residue pairs'] if chemistry_known else 0,
                m['backbone hydrogen bonds'] + m['side-chain hydrogen bonds'] if chemistry_known else 0,
                m['total buried surface area'] / len(self.window) if chemistry_known else 0,
                confidence if confidence is not None else -1,
                plddt if plddt is not None else -1)

    def comparable(self, other):
        first, second = self.ranking(), other.ranking()
        return (abs(first[0] - second[0]) <= .25 and first[1:4] == second[1:4]
                and abs(first[4] - second[4]) <= 10
                and (min(first[5], second[5]) < 0 or abs(first[5] - second[5]) <= .05)
                and (min(first[6], second[6]) < 0 or abs(first[6] - second[6]) <= 5))


def evaluate_sites(descriptions, row, domain, windows, thresholds, parsed, contacts, tight):
    chains = (row['chain 1'], row['chain 2'])
    model = descriptions.models[(row['source'], row['model'])][0]
    receptor = {r.id[1] for r in domain}
    candidates = []
    for window, sites, density in windows:
        motif = {r.id[1] for r in window}
        selected = {pair: distance for pair, distance in contacts.items() if pair[1] in motif}
        patch_row = with_contacts(row, selected)
        context = dict(zip(chains, (receptor, motif)))
        detail = descriptions.describe(
            patch_row, context=context, tight_contacts={pair for pair in tight if pair[1] in motif},
        )
        metrics = {**patch_row, **detail['metrics'],
                   **selected_confidence(chains, (receptor, motif), parsed),
                   **atom_model_ca_metrics(model, chains, selected)}
        scores = [metrics.get(key) for key in ('pDockQ2 1 to 2', 'pDockQ2 2 to 1')]
        if all(score is not None for score in scores):
            metrics['traf_confidence'] = min(scores)
        partners = {a for a, _ in selected}
        evidence = detail['evidence']
        if all(all(atom in r for atom in ('N', 'CA', 'C', 'O')) for r in domain if r.id[1] in partners):
            span = min(len({bond[i] for bond in evidence.backbone_bonds}) for i in (0, 1))
            metrics['traf_feasibility_rule'] = feasibility_rule(metrics, span, thresholds)
        else:
            metrics['traf_feasibility_rule'] = 'unresolved'
            evidence = None
        candidates.append(Candidate(window, sites, density, metrics, evidence))
    return sorted(candidates, key=Candidate.ranking, reverse=True)


def with_contacts(row, contacts):
    patch = {**row, 'contact count': len(contacts)}
    for side in (1, 2):
        numbers = sorted({pair[side - 1] for pair in contacts})
        patch[f'interface residues {side}'] = ','.join(map(str, numbers))
        patch[f'residue count {side}'] = len(numbers)
    return patch


def extract(locator, chains, *, record, seeds, model_scores, dssp, domains,
            motif_windows, math_windows, thresholds, parsed):
    source, number, cif, data = record
    model_id = '' if number is None else number
    residues, pae = parsed
    descriptions = PatchDescriptions({(source, model_id): cif}, dssp)
    model = descriptions.models[(source, model_id)][0]
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
        'contact basis': 'heavy atom',
        'interface residues 1': '', 'interface residues 2': '',
        'traf_paralog': domain if domain in TRAF_PARALOGS else None,
        'traf_math_window': residue_span(math_domain),
        'traf_math_source': 'supplied' if explicit_math else 'inferred',
        'traf_motif_source': 'supplied' if explicit_motif else 'inferred',
        'traf_motif_window': residue_span(motif) if explicit_motif else None,
        'region 1': residue_span(math_domain),
        'region 2': residue_span(motif) if explicit_motif else None,
        'traf_thresholds': json.dumps(asdict(thresholds), sort_keys=True),
        'traf_pose': 'unresolved',
        **model_scores.get((source, model_id), {}),
    }
    if groove:
        row.update(traf_beta6=residue_span(groove.beta6), traf_beta7=residue_span(groove.beta7),
                   traf_polar_loop=residue_span(groove.loop), traf_pocket=residue_span(groove.pocket))
    if '_assembly' in data:
        row['confidence origin'] = 'inherited per-chain values; no assembled confidence'
    if not math_domain:
        row['traf_status'] = ('missing MATH window coordinates' if explicit_math
                              else 'complete terminal MATH topology unresolved')
    elif not motif:
        row['traf_status'] = 'missing motif coordinates'
    else:
        contacts = heavy_contacts(math_domain, motif, 6.0)
        tight = {pair for pair, distance in contacts.items() if distance <= TIGHT_CONTACT_CUTOFF}
        if not contacts:
            row.update(with_contacts(row, contacts))
            row.update({'tight-contact cutoff': TIGHT_CONTACT_CUTOFF,
                        'tight-contact residues 1': 0, 'tight-contact residues 2': 0})
            row.update(traf_pose='no_contact', traf_status='no_contact', traf_feasibility_rule='no_contact')
        else:
            windows = ([(motif, {b for _, b in tight}, len(tight) / len(motif))]
                       if explicit_motif else candidate_windows(motif, contacts, tight))
            if pae is not None and not all(math.isfinite(float(x)) and float(x) >= 0 for line in pae for x in line):
                parsed = residues, None
            candidates = evaluate_sites(descriptions, row, math_domain, windows, thresholds, parsed, contacts, tight)
            if not candidates:
                row['traf_status'] = 'contacting motif lacks a contiguous backbone window'
            else:
                update_site(row, candidates, groove)
    attach_patch_identity(row, model_identity(cif, data, seed), number)
    return Evaluation([row], {source: fold}, {source: sequence}, descriptions.models)


def update_site(row, candidates, groove):
    best = candidates[0]
    row.update(best.metrics)
    row['region 1'] = row['traf_math_window']
    row['region 2'] = row['traf_motif_window'] = residue_span(best.window)
    if best.evidence is not None:
        pose, register, evidence = assess_pose(groove, best.window, best.evidence)
    else:
        pose, register, evidence = 'unresolved', None, 'incomplete contacting MATH backbone'
    competitors = []
    for candidate in candidates[1:]:
        if not best.comparable(candidate):
            continue
        alternate_pose, alternate_register, _ = (assess_pose(groove, candidate.window, candidate.evidence)
                                                 if candidate.evidence is not None else (None, None, None))
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
