"""Position-wise protein sequence channel definitions."""

import math
import os
import re
from abc import ABC, abstractmethod
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from numbers import Real
from types import MappingProxyType


CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

# Kyte and Doolittle, J Mol Biol. 1982;157:105-132.
KYTE_DOOLITTLE_HYDROPATHY = MappingProxyType({
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
})

# Eisenberg et al., Proc Natl Acad Sci USA. 1984;81:140-144.
EISENBERG_HYDROPHOBICITY = MappingProxyType({
    "A": 0.25,
    "R": -1.76,
    "N": -0.64,
    "D": -0.72,
    "C": 0.04,
    "Q": -0.69,
    "E": -0.62,
    "G": 0.16,
    "H": -0.40,
    "I": 0.73,
    "L": 0.53,
    "K": -1.10,
    "M": 0.26,
    "F": 0.61,
    "P": -0.07,
    "S": -0.26,
    "T": -0.18,
    "W": 0.37,
    "Y": 0.02,
    "V": 0.54,
})

# Campen et al., Protein Pept Lett. 2008;15:956-963.
TOP_IDP_DISORDER_PROPENSITY = MappingProxyType({
    "A": 0.060,
    "R": 0.180,
    "N": 0.007,
    "D": 0.192,
    "C": 0.020,
    "Q": 0.318,
    "E": 0.736,
    "G": 0.166,
    "H": 0.303,
    "I": -0.486,
    "L": -0.326,
    "K": 0.586,
    "M": -0.397,
    "F": -0.697,
    "P": 0.987,
    "S": 0.341,
    "T": 0.059,
    "W": -0.884,
    "Y": -0.510,
    "V": -0.121,
})

# Simplified side-chain charge convention used by localCIDER at neutral pH.
NEUTRAL_PH_CHARGE = MappingProxyType({
    amino_acid: 1.0 if amino_acid in "KR" else -1.0 if amino_acid in "DE" else 0.0
    for amino_acid in CANONICAL_AMINO_ACIDS
})


def _validate_radius(radius, *, allow_zero=True):
    if isinstance(radius, bool) or not isinstance(radius, int):
        raise TypeError("radius must be an integer")
    minimum = 0 if allow_zero else 1
    if radius < minimum:
        qualifier = "nonnegative" if allow_zero else "at least 1"
        raise ValueError(f"radius must be {qualifier}")
    return radius


def _validate_sequence(sequence):
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")
    sequence = sequence.upper()
    invalid = sorted(set(sequence) - CANONICAL_AMINO_ACIDS)
    if invalid:
        raise ValueError(
            "sequence contains noncanonical amino acids: " + ", ".join(invalid)
        )
    return sequence


def _validate_background_values(bg_mu, bg_sigma):
    if (
        isinstance(bg_mu, bool)
        or not isinstance(bg_mu, Real)
        or not math.isfinite(bg_mu)
    ):
        raise ValueError("bg_mu must be finite")
    if (
        isinstance(bg_sigma, bool)
        or not isinstance(bg_sigma, Real)
        or not math.isfinite(bg_sigma)
        or bg_sigma <= 0
    ):
        raise ValueError("bg_sigma must be finite and greater than zero")
    return float(bg_mu), float(bg_sigma)


def _validate_scale(scale):
    if not isinstance(scale, Mapping):
        raise TypeError("scale must be a mapping")
    normalized = {}
    for amino_acid, value in scale.items():
        if not isinstance(amino_acid, str) or len(amino_acid) != 1:
            raise ValueError("scale keys must be one-letter amino-acid codes")
        amino_acid = amino_acid.upper()
        if (
            isinstance(value, bool)
            or not isinstance(value, Real)
            or not math.isfinite(value)
        ):
            raise ValueError(f"scale value for {amino_acid!r} must be finite")
        normalized[amino_acid] = float(value)
    if set(normalized) != CANONICAL_AMINO_ACIDS:
        missing = sorted(CANONICAL_AMINO_ACIDS - set(normalized))
        extra = sorted(set(normalized) - CANONICAL_AMINO_ACIDS)
        details = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if extra:
            details.append("unexpected " + ", ".join(extra))
        raise ValueError(
            "scale must define exactly the 20 canonical amino acids ("
            + "; ".join(details)
            + ")"
        )
    return MappingProxyType(normalized)


def _validate_residues(residues):
    if isinstance(residues, str):
        residues = set(residues.upper())
    else:
        try:
            normalized = set()
            for residue in residues:
                if not isinstance(residue, str) or len(residue) != 1:
                    raise ValueError(
                        "residues must contain one-letter amino-acid codes"
                    )
                normalized.add(residue.upper())
            residues = normalized
        except TypeError as error:
            raise TypeError("residues must be a string or iterable") from error
    if not residues:
        raise ValueError("residues must not be empty")
    invalid = sorted(residues - CANONICAL_AMINO_ACIDS)
    if invalid:
        raise ValueError(
            "residues contains noncanonical amino acids: " + ", ".join(invalid)
        )
    return frozenset(residues)


def _validate_dipeptide(dipeptide):
    if not isinstance(dipeptide, str) or len(dipeptide) != 2:
        raise ValueError("dipeptide must contain exactly two amino acids")
    dipeptide = dipeptide.upper()
    invalid = sorted(set(dipeptide) - CANONICAL_AMINO_ACIDS)
    if invalid:
        raise ValueError(
            "dipeptide contains noncanonical amino acids: " + ", ".join(invalid)
        )
    return dipeptide


def _validate_angle(angle):
    if (
        isinstance(angle, bool)
        or not isinstance(angle, Real)
        or not math.isfinite(angle)
        or angle <= 0
        or angle > 360
    ):
        raise ValueError("angle must be finite and in the interval (0, 360]")
    return float(angle)


def _window_bounds(sequence_length, index, radius):
    return max(0, index - radius), min(sequence_length, index + radius + 1)


def _window_values(sequence, radius, score):
    values = []
    for index in range(len(sequence)):
        start, end = _window_bounds(len(sequence), index, radius)
        values.append(score(sequence, start, end))
    return values


def _mean_scale_values(sequence, start, end, scale):
    return sum(scale[amino_acid] for amino_acid in sequence[start:end]) / (end - start)


def iter_fasta_records(
    fasta_path,
    on_invalid_sequence=None,
    on_duplicate_sequence_id=None,
):
    """Yield ``(sequence_id, sequence)`` records from a protein FASTA file."""
    path = os.fspath(fasta_path)
    found_record = False
    sequence_id = None
    sequence_lines = []
    seen_sequence_ids = set()
    skip_record = False
    with open(path, encoding="utf-8") as fasta:
        for line_number, raw_line in enumerate(fasta, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sequence_id is not None:
                    if not sequence_lines:
                        raise ValueError(f"FASTA record {sequence_id!r} has no sequence")
                    if not skip_record:
                        sequence = "".join(sequence_lines)
                        try:
                            sequence = _validate_sequence(sequence)
                        except ValueError as error:
                            if on_invalid_sequence is None:
                                raise
                            on_invalid_sequence(sequence_id, error)
                        else:
                            yield sequence_id, sequence
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"FASTA header on line {line_number} is empty")
                sequence_id = header.split()[0]
                if sequence_id in seen_sequence_ids:
                    if on_duplicate_sequence_id is None:
                        raise ValueError(f"duplicate FASTA sequence ID: {sequence_id}")
                    on_duplicate_sequence_id(sequence_id)
                    skip_record = True
                else:
                    seen_sequence_ids.add(sequence_id)
                    skip_record = False
                found_record = True
                sequence_lines = []
            else:
                if sequence_id is None:
                    raise ValueError(
                        f"FASTA sequence data precedes the first header on line {line_number}"
                    )
                sequence_lines.append(line)
    if sequence_id is not None:
        if not sequence_lines:
            raise ValueError(f"FASTA record {sequence_id!r} has no sequence")
        if not skip_record:
            sequence = "".join(sequence_lines)
            try:
                sequence = _validate_sequence(sequence)
            except ValueError as error:
                if on_invalid_sequence is None:
                    raise
                on_invalid_sequence(sequence_id, error)
            else:
                yield sequence_id, sequence
    if not found_record:
        raise ValueError("FASTA file contains no records")


@dataclass(frozen=True)
class ChannelBackground:
    bg_mu: float
    bg_sigma: float

    def __post_init__(self):
        bg_mu, bg_sigma = _validate_background_values(self.bg_mu, self.bg_sigma)
        object.__setattr__(self, "bg_mu", bg_mu)
        object.__setattr__(self, "bg_sigma", bg_sigma)

    def as_dict(self):
        return {"bg_mu": self.bg_mu, "bg_sigma": self.bg_sigma}


class _BackgroundAccumulator:

    def __init__(self, channel):
        self.channel = channel
        self.count = 0
        self.mean = 0.0
        self.sum_squared_deviations = 0.0

    def add_sequence(self, sequence):
        for value in self.channel.raw_values(sequence):
            if not math.isfinite(value):
                raise ValueError(
                    f"channel {self.channel.short_name!r} produced a nonfinite value"
                )
            self.count += 1
            difference = value - self.mean
            self.mean += difference / self.count
            self.sum_squared_deviations += difference * (value - self.mean)

    def result(self):
        if self.count == 0:
            raise ValueError(
                f"channel {self.channel.short_name!r} has no usable positions"
            )
        sigma = math.sqrt(max(0.0, self.sum_squared_deviations / self.count))
        if not math.isfinite(sigma) or sigma <= 0:
            raise ValueError(
                f"channel {self.channel.short_name!r} background standard deviation "
                "is zero or nonfinite"
            )
        return ChannelBackground(self.mean, sigma)


class Channel(ABC):
    """Base class for one position-wise protein sequence channel."""

    kind = None
    requires_background = True

    def __init__(self, short_name):
        if not isinstance(short_name, str) or not re.fullmatch(
            r"[A-Za-z_][A-Za-z0-9_]*", short_name
        ):
            raise ValueError(
                "short_name must start with a letter or underscore and contain "
                "only letters, numbers, and underscores"
            )
        self.short_name = short_name

    def raw_values(self, sequence):
        """Return unnormalized values, one for each residue in ``sequence``."""
        return self._raw_values(_validate_sequence(sequence))

    @abstractmethod
    def _raw_values(self, sequence):
        raise NotImplementedError

    @abstractmethod
    def arguments(self):
        """Return a JSON-compatible description of channel arguments."""
        raise NotImplementedError

    def definition(self):
        return {
            "short_name": self.short_name,
            "type": self.kind,
            "arguments": self.arguments(),
        }

    def background_accumulator(self):
        if not self.requires_background:
            return None
        return _BackgroundAccumulator(self)

    def background_sequence_is_usable(self, sequence):
        return True

    def compute_background(self, fasta_path):
        """Compute residue-weighted population statistics from a FASTA file."""
        accumulator = self.background_accumulator()
        if accumulator is None:
            return None
        for _sequence_id, sequence in iter_fasta_records(fasta_path):
            if self.background_sequence_is_usable(sequence):
                accumulator.add_sequence(sequence)
        return accumulator.result()

    def make_function(self, background=None):
        """Return ``f(sequence)`` for this channel."""
        if not self.requires_background:
            if background is not None:
                raise ValueError(
                    f"channel {self.short_name!r} does not use a background"
                )
            return self.raw_values
        if not isinstance(background, ChannelBackground):
            raise TypeError(
                f"channel {self.short_name!r} requires a ChannelBackground"
            )

        def function(sequence):
            return [
                (value - background.bg_mu) / background.bg_sigma
                for value in self.raw_values(sequence)
            ]

        return function


class CompositionBiasChannel(Channel):
    kind = "composition_bias"

    def __init__(self, short_name, radius, residues):
        super().__init__(short_name)
        self.radius = _validate_radius(radius)
        self.residues = _validate_residues(residues)

    def _raw_values(self, sequence):
        return _window_values(
            sequence,
            self.radius,
            lambda seq, start, end: sum(
                amino_acid in self.residues for amino_acid in seq[start:end]
            ) / (end - start),
        )

    def arguments(self):
        return {"radius": self.radius, "residues": "".join(sorted(self.residues))}


class ShortMotifChannel(Channel):
    kind = "short_motif"
    requires_background = False

    def __init__(self, short_name, pattern, flags=0):
        super().__init__(short_name)
        if isinstance(flags, bool) or not isinstance(flags, int):
            raise TypeError("flags must be an integer")
        try:
            self._compiled = re.compile(pattern, flags)
        except (TypeError, re.error) as error:
            raise ValueError(f"invalid motif pattern: {error}") from error
        self.pattern = self._compiled.pattern
        self.flags = flags
        empty_match = self._compiled.search("")
        if empty_match is not None and empty_match.start() == empty_match.end():
            raise ValueError("motif pattern must not produce zero-length matches")

    def _raw_values(self, sequence):
        result = [0] * len(sequence)
        for start in range(len(sequence) + 1):
            match = self._compiled.match(sequence, start)
            if match is None:
                continue
            match_start, match_end = match.span()
            if match_start == match_end:
                raise ValueError("motif pattern produced a zero-length match")
            for index in range(match_start, match_end):
                result[index] = 1
        return result

    def arguments(self):
        return {"pattern": self.pattern, "flags": self.flags}


class NetChargeChannel(Channel):
    kind = "net_charge"

    def __init__(self, short_name, radius):
        super().__init__(short_name)
        self.radius = _validate_radius(radius)

    def _raw_values(self, sequence):
        return _window_values(
            sequence,
            self.radius,
            lambda seq, start, end: _mean_scale_values(
                seq, start, end, NEUTRAL_PH_CHARGE
            ),
        )

    def arguments(self):
        return {"radius": self.radius}


class HydropathyChannel(Channel):
    kind = "hydropathy"

    def __init__(self, short_name, radius, scale=KYTE_DOOLITTLE_HYDROPATHY):
        super().__init__(short_name)
        self.radius = _validate_radius(radius)
        self.scale = _validate_scale(scale)

    def _raw_values(self, sequence):
        return _window_values(
            sequence,
            self.radius,
            lambda seq, start, end: _mean_scale_values(
                seq, start, end, self.scale
            ),
        )

    def arguments(self):
        return {"radius": self.radius, "scale": dict(sorted(self.scale.items()))}


def _hydrophobic_moment(sequence, start, end, scale, angle_radians):
    x = 0.0
    y = 0.0
    for offset, amino_acid in enumerate(sequence[start:end]):
        hydrophobicity = scale[amino_acid]
        phase = offset * angle_radians
        x += hydrophobicity * math.cos(phase)
        y += hydrophobicity * math.sin(phase)
    return math.hypot(x, y) / (end - start)


class HydrophobicMomentChannel(Channel):
    kind = "hydrophobic_moment"

    def __init__(
        self,
        short_name,
        radius,
        angle=100.0,
        scale=EISENBERG_HYDROPHOBICITY,
    ):
        super().__init__(short_name)
        self.radius = _validate_radius(radius)
        self.angle = _validate_angle(angle)
        self.scale = _validate_scale(scale)

    def _raw_values(self, sequence):
        angle_radians = math.radians(self.angle)
        return _window_values(
            sequence,
            self.radius,
            lambda seq, start, end: _hydrophobic_moment(
                seq, start, end, self.scale, angle_radians
            ),
        )

    def arguments(self):
        return {
            "radius": self.radius,
            "angle": self.angle,
            "scale": dict(sorted(self.scale.items())),
        }


class DisorderPropensityChannel(Channel):
    kind = "disorder_propensity"

    def __init__(
        self,
        short_name,
        radius,
        scale=TOP_IDP_DISORDER_PROPENSITY,
    ):
        super().__init__(short_name)
        self.radius = _validate_radius(radius)
        self.scale = _validate_scale(scale)

    def _raw_values(self, sequence):
        return _window_values(
            sequence,
            self.radius,
            lambda seq, start, end: _mean_scale_values(
                seq, start, end, self.scale
            ),
        )

    def arguments(self):
        return {"radius": self.radius, "scale": dict(sorted(self.scale.items()))}


def _sequence_entropy(sequence, start, end):
    length = end - start
    counts = Counter(sequence[start:end])
    return -sum(
        (count / length) * math.log2(count / length)
        for count in counts.values()
    )


class SequenceEntropyChannel(Channel):
    kind = "sequence_entropy"

    def __init__(self, short_name, radius):
        super().__init__(short_name)
        self.radius = _validate_radius(radius)

    def _raw_values(self, sequence):
        return _window_values(sequence, self.radius, _sequence_entropy)

    def arguments(self):
        return {"radius": self.radius}


def _dipeptide_frequency(sequence, start, end, dipeptide):
    pair_count = end - start - 1
    return sum(
        sequence[index:index + 2] == dipeptide
        for index in range(start, end - 1)
    ) / pair_count


class DipeptideFrequencyChannel(Channel):
    kind = "dipeptide_frequency"

    def __init__(self, short_name, radius, dipeptide):
        super().__init__(short_name)
        self.radius = _validate_radius(radius, allow_zero=False)
        self.dipeptide = _validate_dipeptide(dipeptide)

    def _raw_values(self, sequence):
        if len(sequence) < 2:
            raise ValueError(
                "dipeptide frequency requires a sequence of at least two residues"
            )
        return _window_values(
            sequence,
            self.radius,
            lambda seq, start, end: _dipeptide_frequency(
                seq, start, end, self.dipeptide
            ),
        )

    def arguments(self):
        return {"radius": self.radius, "dipeptide": self.dipeptide}

    def background_sequence_is_usable(self, sequence):
        return len(sequence) >= 2
