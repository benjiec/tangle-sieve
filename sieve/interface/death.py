"""Death-family patch evaluation using the interactome face rules."""
from dataclasses import dataclass
from math import isfinite

from sieve.alphafold_death_rules import FACES, candidate_interface, classify
from sieve.alphafold_helix_normalization import normalize_helices
from sieve.alphafold_patch_metrics import extract_patch_metrics
from sieve.alphafold_pdockq2_patches import Options
from sieve.alphafold_structure import canonical, digest

# Transcribed from coral/notes/interactome/dd_face_rules.tsv.
# Domain, face, element, minimum contacting residues, requirement, notes.
RULE_ROWS = (
    ('DD', 'Ia', 'alpha1', 3, 'mandatory', 'core helix'),
    ('DD', 'Ia', 'alpha4', 3, 'mandatory', 'core helix'),
    ('DD', 'Ia', 'loop1-2', 0, 'optional', 'supporting loop'),
    ('DD', 'Ib', 'alpha2', 3, 'mandatory', 'core helix'),
    ('DD', 'Ib', 'alpha3', 3, 'mandatory', 'core helix'),
    ('DD', 'Ib', 'loop2-3', 0, 'optional', 'supporting loop'),
    ('DD', 'IIa', 'alpha4', 3, 'mandatory', 'core helix'),
    ('DD', 'IIa', 'loop4-5', 1, 'mandatory', 'anchoring loop'),
    ('DD', 'IIb', 'alpha6', 3, 'mandatory', 'core helix'),
    ('DD', 'IIb', 'loop1-2', 1, 'alternative_or', 'at least 1 residue from loop1-2 OR loop5-6'),
    ('DD', 'IIb', 'loop5-6', 1, 'alternative_or', 'at least 1 residue from loop1-2 OR loop5-6'),
    ('DD', 'IIIa', 'alpha3', 3, 'mandatory', 'core helix'),
    ('DD', 'IIIa', 'loop2-3', 1, 'mandatory', 'supporting loop'),
    ('DD', 'IIIb', 'loop1-2', 1, 'mandatory', 'core loop'),
    ('DD', 'IIIb', 'loop3-4', 1, 'mandatory', 'core loop'),
    ('DED', 'Ia', 'alpha1', 3, 'mandatory', 'core helix'),
    ('DED', 'Ia', 'alpha4', 3, 'mandatory', 'core helix'),
    ('DED', 'Ia', 'loop1-2', 0, 'optional', 'supporting loop'),
    ('DED', 'Ib', 'alpha2', 3, 'mandatory', 'core helix containing FL/RxDL motif'),
    ('DED', 'Ib', 'alpha3', 1, 'alternative_or', 'min 1 on alpha3 OR min 2 on loop2-3'),
    ('DED', 'Ib', 'loop2-3', 2, 'alternative_or', 'min 2 on loop2-3 OR min 1 on alpha3'),
    ('DED', 'Ib', 'alpha5', 0, 'optional', 'permitted peripheral contact if alpha3 <= 6 residues'),
    ('DED', 'IIa', 'alpha4', 2, 'mandatory', 'core helix'),
    ('DED', 'IIa', 'loop4-5', 1, 'mandatory', 'anchoring loop'),
    ('DED', 'IIb', 'alpha6', 3, 'mandatory', 'core helix'),
    ('DED', 'IIb', 'loop1-2', 1, 'alternative_or', 'at least 1 residue from loop1-2 OR loop5-6'),
    ('DED', 'IIb', 'loop5-6', 1, 'alternative_or', 'at least 1 residue from loop1-2 OR loop5-6'),
    ('DED', 'IIIa', 'alpha3', 1, 'alternative_or', 'min 1 on alpha3 OR min 2 on loop2-3'),
    ('DED', 'IIIa', 'loop2-3', 2, 'alternative_or', 'min 2 on loop2-3 OR min 1 on alpha3'),
    ('DED', 'IIIb', 'loop1-2', 1, 'mandatory', 'core loop'),
    ('DED', 'IIIb', 'loop3-4', 1, 'mandatory', 'core loop'),
    ('CARD', 'Ia', 'alpha1', 3, 'mandatory', 'core helix (predominantly alpha1b)'),
    ('CARD', 'Ia', 'alpha4', 3, 'mandatory', 'core helix'),
    ('CARD', 'Ia', 'loop1-2', 0, 'optional', 'supporting loop'),
    ('CARD', 'Ib', 'alpha2', 3, 'mandatory', 'core helix containing acidic/polar patch'),
    ('CARD', 'Ib', 'alpha3', 2, 'mandatory', 'core helix'),
    ('CARD', 'Ib', 'loop2-3', 0, 'optional', 'supporting loop'),
    ('CARD', 'IIa', 'alpha4', 2, 'mandatory', 'core helix'),
    ('CARD', 'IIa', 'loop4-5', 1, 'mandatory', 'anchoring loop'),
    ('CARD', 'IIb', 'alpha6', 2, 'alternative_or', 'min 2 on alpha6 OR min 2 on loop5-6'),
    ('CARD', 'IIb', 'loop5-6', 2, 'alternative_or', 'min 2 on alpha6 OR min 2 on loop5-6'),
    ('CARD', 'IIb', 'loop1-2', 1, 'alternative_or', 'at least 1 residue from loop1-2 OR loop5-6'),
    ('CARD', 'IIIa', 'alpha3', 2, 'mandatory', 'core helix'),
    ('CARD', 'IIIa', 'loop2-3', 1, 'mandatory', 'supporting loop'),
    ('CARD', 'IIIb', 'loop1-2', 1, 'mandatory', 'core loop (loop after alpha1b)'),
    ('CARD', 'IIIb', 'loop3-4', 1, 'mandatory', 'core loop'),
)

RULES = {}
for domain, face, element, minimum, requirement, notes in RULE_ROWS:
    RULES.setdefault(domain, {}).setdefault(face, []).append(
        dict(element=element, min_residues=minimum, requirement_type=requirement, notes=notes)
    )


def rules_id(rules):
    return digest(canonical([
        dict(domain=domain, face=face, **rule)
        for domain in sorted(rules)
        for face in FACES
        for rule in sorted(rules[domain].get(face, ()), key=lambda rule: rule['element'])
    ]))


RULES_ID = rules_id(RULES)


def extract(locator, chains, *, record, seeds, model_scores, dssp, distance_cutoff=10.0):
    return extract_patch_metrics(
        locator, options=Options(cutoff=distance_cutoff), chains=chains,
        model_records=[record], model_seeds=seeds, model_scores=model_scores, dssp=dssp,
    )


def available_number(value):
    return isinstance(value, (int, float)) and isfinite(value)


def confidence(patch):
    score = patch.get('pDockQ2 max')
    area = patch.get('total buried surface area')
    if not all(available_number(value) for value in (score, area)):
        return None
    if score < 0.23:
        return 'Reject'
    if score >= 0.70:
        return 'High'
    if (area < 1000 and score >= 0.40) or score >= 0.50:
        return 'Medium'
    return 'Low'


def feasibility(patch):
    area = patch.get('total buried surface area')
    carbon_sulfur = patch.get('total C/S contribution percent')
    contacts = patch.get('nonpolar contacts')
    salts = patch.get('salt bridge candidates')
    salt_distance = patch.get('salt bridge candidate heavy-atom distance min')
    charge = patch.get('nominal patch charge difference')
    overlaps = patch.get('steric overlap candidates')
    ca_distance = patch.get('contact distance min')
    if not all(available_number(value) for value in (
        area, carbon_sulfur, contacts, salts, salt_distance, charge, overlaps, ca_distance,
    )) or area <= 0:
        return None

    density = 100 * contacts / area
    hovering = (density <= 0.4 and ca_distance > 4.5 and salt_distance > 4) or overlaps >= 3
    short_salt = salts >= 1 and salt_distance <= 4
    lateral = (
        500 <= area < 800 and carbon_sulfur >= 40 and overlaps <= 1
        and ((density > 0.25 and short_salt) or density > 0.70)
    )
    cooperative = (
        800 <= area < 1400 and carbon_sulfur >= 40 and overlaps <= 1
        and ((density > 0.45 and short_salt) or density > 0.70)
    )
    zipper = (
        700 <= area < 1200 and carbon_sulfur >= 25 and salts >= 2 and salt_distance <= 3
        and abs(charge) >= 3 and density > 0.2 and overlaps <= 1
    )
    clamp = (
        area >= 1800 and carbon_sulfur >= 45 and density > 0.5 and overlaps <= 2
        and (density > 0.6 or (salts >= 1 and salt_distance <= 3.5))
    )
    hydrophobic = 1000 <= area < 1800 and carbon_sulfur >= 40 and salts >= 0 and overlaps <= 2
    return not hovering and (lateral or cooperative or zipper or clamp or hydrophobic)


@dataclass
class Classification:
    interface: str | None
    sides: dict
    settings: dict


def classify_interface(patch, structure, domains, *, chains, rules=RULES,
                       merge_kinked_helices=True, connector_length=1, angle=75):
    model, codes, segments, dssp_status = structure
    sides = {}
    for chain in chains:
        side = 1 if chain == patch['chain 1'] else 2
        if dssp_status['status'] != 'assigned':
            sides[chain] = {'candidate face': None, 'face status': str(dssp_status), 'evidence': []}
            continue
        selected = {int(value) for value in patch[f'interface residues {side}'].split(',') if value}
        chain_segments, report = segments, ''
        if merge_kinked_helices:
            chain_segments, report = normalize_helices(
                model, chain, codes, segments, connector_length=connector_length, angle=angle,
            )
        faces, evidence, status = classify(chain, selected, codes, chain_segments, rules[domains[chain]])
        if report:
            status += '; ' + report
        sides[chain] = {'candidate face': faces, 'face status': status, 'evidence': evidence}
    interface = candidate_interface(*(side['candidate face'] for side in sides.values()))
    settings = {
        'rules_id': RULES_ID if rules is RULES else rules_id(rules),
        'merge kinked helices': merge_kinked_helices,
        'kink connector length': connector_length,
        'kink angle': angle,
    }
    return Classification(interface, sides, settings)
