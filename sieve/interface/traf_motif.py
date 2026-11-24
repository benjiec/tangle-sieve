"""Sequence motifs evaluated at a selected TRAF binding site."""
from dataclasses import dataclass

from Bio.Data.PDBData import protein_letters_3to1

STANDARD_AMINO_ACIDS = frozenset(protein_letters_3to1.values())


@dataclass(frozen=True)
class MotifRule:
    name: str
    canonical: bool
    patterns: tuple[tuple[str, ...], ...]


@dataclass(frozen=True)
class MotifMatch:
    rule: MotifRule
    start: int
    end: int


# An empty selector is a wildcard for an observed standard amino acid.
RULES = {
    'traf6': (
        MotifRule('traf6_binding_strict', True, (('P', '', 'E', '', '', 'FYWDE'),)),
        MotifRule('traf6_binding_permissive', False, (('P', '', 'QE', '', '', 'FYWHDEILVKR'),)),
    ),
    'traf2/3': (
        MotifRule('traf2_binding_strict', True, (('P', '', 'QE', '', 'ST'),)),
        MotifRule('traf2_binding_major_consensus', True, (('PSAT', '', 'QE', 'E'),)),
        MotifRule('traf2_binding_variant', True, (('P', '', 'QEN', '', 'STEDAG'),)),
        MotifRule('traf2_binding_permissive', False, tuple(
            ('P', *(('',) * before), 'QEN', *(('',) * after), 'STEDAG')
            for before in (1, 2) for after in (1, 2)
        )),
    ),
}


def match_pattern(sequence, start, pattern, contacts):
    if contacts is not None and sum(start + offset in contacts for offset in range(len(pattern))) < 2:
        return False
    unknown = False
    for offset, allowed in enumerate(pattern):
        letter = sequence.get(start + offset)
        if letter not in STANDARD_AMINO_ACIDS:
            unknown = True
        elif allowed and letter not in allowed:
            return False
    return None if unknown else True


def match_motifs(sequence, start, contacts=None):
    """Return the strongest rule and shortest matching span per family."""
    matches = {}
    uncertain = set()
    for family, rules in RULES.items():
        possible = False
        for rule in rules:
            states = [match_pattern(sequence, start, pattern, contacts) for pattern in rule.patterns]
            lengths = [len(pattern) for pattern, state in zip(rule.patterns, states) if state is True]
            if lengths:
                matches[family] = MotifMatch(rule, start, start + min(lengths) - 1)
                if rule.canonical:
                    possible = False
                break
            possible |= None in states
        if possible:
            uncertain.add(family)
    return matches, uncertain


def motif_faces(sequence, contacts, *, anchor, window):
    if anchor is not None:
        return match_motifs(sequence, anchor)
    if not window:
        return {}, set()

    # Without a pose register, report a motif only when a single start in the
    # selected contact site has sequence support.
    first, last = window
    hits = []
    for start in range(first - 1, last + 1):
        matches, _ = match_motifs(sequence, start, contacts)
        if matches:
            hits.append(matches)
            if len(hits) > 1:
                return {}, set()
    return (hits[0] if hits else {}), set()
