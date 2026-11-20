"""Select one interface profile for a pair of annotated chains."""
from collections.abc import Callable
from dataclasses import dataclass

from sieve.interface import death


@dataclass(frozen=True)
class InterfaceProfile:
    name: str
    domains_a: frozenset[str]
    domains_b: frozenset[str]
    extract: Callable
    confidence: Callable
    feasibility: Callable
    classify: Callable
    columns: tuple[str, ...] = ()


DEATH_DOMAINS = frozenset(death.RULES)
INTERFACES = (
    InterfaceProfile('death', DEATH_DOMAINS, DEATH_DOMAINS, death.extract,
                     death.confidence, death.feasibility, death.classify_interface,
                     ('rules_id', 'merge kinked helices', 'kink connector length', 'kink angle')),
)


def resolve(chains, domains, profiles=INTERFACES):
    left, right = chains
    a, b = domains[left], domains[right]
    matches = []
    for profile in profiles:
        if a in profile.domains_a and b in profile.domains_b:
            matches.append((profile, chains))
        elif b in profile.domains_a and a in profile.domains_b:
            matches.append((profile, (right, left)))
    pair = f'{left}:{a}, {right}:{b}'
    if not matches:
        raise ValueError(f'unsupported interface domains ({pair})')
    if len(matches) > 1:
        names = ', '.join(profile.name for profile, _ in matches)
        raise ValueError(f'ambiguous interface domains ({pair}): {names}')
    return matches[0]
