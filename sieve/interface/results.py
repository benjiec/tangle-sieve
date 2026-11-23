"""The result contracts used by interface profiles and their common consumer."""
from collections.abc import Iterable
from dataclasses import dataclass, field


@dataclass
class Classification:
    interface: str | None
    sides: dict
    settings: dict


@dataclass
class Evaluation:
    rows: Iterable[dict]
    fold_ids: dict
    sequence_ids: dict
    structures: dict = field(default_factory=dict)
