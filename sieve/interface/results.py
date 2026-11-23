"""The result contracts used by interface profiles and their common consumer."""
from collections.abc import Iterable
from dataclasses import dataclass, field
from enum import Enum
from math import isfinite


class Confidence(str, Enum):
    HIGH = 'High'
    MEDIUM = 'Medium'
    LOW = 'Low'
    REJECT = 'Reject'


class InterfaceLabel(str, Enum):
    CANONICAL = 'Canonical'
    ATYPICAL = 'Atypical'
    VARIANT = 'Variant'


RESULT_COLUMNS = ('confidence', 'chemical feasibility', 'interface')


def available_number(value):
    return isinstance(value, (int, float)) and isfinite(value)


def confidence_level(score, *, medium_cutoff=0.50):
    if not available_number(score):
        return None
    if score >= 0.70:
        return Confidence.HIGH
    if score >= medium_cutoff:
        return Confidence.MEDIUM
    if score >= 0.23:
        return Confidence.LOW
    return Confidence.REJECT


def result_fields(confidence, feasibility, interface):
    """Validate the shared vocabulary and return the values for TSV export."""
    if feasibility is not None and not isinstance(feasibility, bool):
        raise ValueError('chemical feasibility must be a boolean or None')
    return {
        'confidence': Confidence(confidence).value if confidence is not None else None,
        'chemical feasibility': feasibility,
        'interface': InterfaceLabel(interface).value if interface is not None else '',
    }


@dataclass
class Classification:
    interface: InterfaceLabel | None
    sides: dict
    settings: dict


@dataclass
class Evaluation:
    rows: Iterable[dict]
    fold_ids: dict
    sequence_ids: dict
    structures: dict = field(default_factory=dict)
