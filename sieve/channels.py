"""Position-wise protein sequence channels and proteome backgrounds."""

import math
import os
import re
from collections import Counter
from collections.abc import Callable, Iterable, Mapping
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

NEUTRAL_PH_CHARGE = MappingProxyType({
    aa: 1.0 if aa in "KR" else -1.0 if aa in "DE" else 0.0
    for aa in CANONICAL_AMINO_ACIDS
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


def _validate_background(bg_mu, bg_sigma):
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
    return normalized


def _validate_residues(residues):
    if isinstance(residues, str):
        residues = set(residues.upper())
    elif isinstance(residues, Iterable):
        normalized = set()
        for residue in residues:
            if not isinstance(residue, str) or len(residue) != 1:
                raise ValueError("residues must contain one-letter amino-acid codes")
            normalized.add(residue.upper())
        residues = normalized
    else:
        raise TypeError("residues must be a string or iterable")
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


def _window_bounds(sequence_length, index, radius):
    return max(0, index - radius), min(sequence_length, index + radius + 1)


def _window_values(sequence, radius, score):
    values = []
    for index in range(len(sequence)):
        start, end = _window_bounds(len(sequence), index, radius)
        values.append(score(sequence, start, end))
    return values


def _zscore_channel(raw_values, bg_mu, bg_sigma):
    bg_mu, bg_sigma = _validate_background(bg_mu, bg_sigma)

    def channel(sequence):
        sequence = _validate_sequence(sequence)
        return [(value - bg_mu) / bg_sigma for value in raw_values(sequence)]

    return channel


def _mean_scale_values(sequence, start, end, scale):
    return sum(scale[amino_acid] for amino_acid in sequence[start:end]) / (end - start)


def _composition_values(sequence, radius, residues):
    return _window_values(
        sequence,
        radius,
        lambda seq, start, end: sum(
            amino_acid in residues for amino_acid in seq[start:end]
        ) / (end - start),
    )


def _scale_values(sequence, radius, scale):
    return _window_values(
        sequence,
        radius,
        lambda seq, start, end: _mean_scale_values(seq, start, end, scale),
    )


def mk_composition_bias(radius, residues, bg_mu, bg_sigma):
    """Return a channel for the local fraction of a specified residue set."""
    radius = _validate_radius(radius)
    residues = _validate_residues(residues)

    def raw_values(sequence):
        return _composition_values(sequence, radius, residues)

    return _zscore_channel(raw_values, bg_mu, bg_sigma)


def mk_short_motif(pattern):
    """Return a binary channel marking residues covered by regex matches."""
    try:
        compiled = re.compile(pattern)
    except (TypeError, re.error) as error:
        raise ValueError(f"invalid motif pattern: {error}") from error
    empty_match = compiled.search("")
    if empty_match is not None and empty_match.start() == empty_match.end():
        raise ValueError("motif pattern must not produce zero-length matches")

    def channel(sequence):
        sequence = _validate_sequence(sequence)
        result = [0] * len(sequence)
        for start in range(len(sequence) + 1):
            match = compiled.match(sequence, start)
            if match is None:
                continue
            match_start, match_end = match.span()
            if match_start == match_end:
                raise ValueError("motif pattern produced a zero-length match")
            for index in range(match_start, match_end):
                result[index] = 1
        return result

    return channel


def mk_net_charge(radius, bg_mu, bg_sigma):
    """Return a channel for local net charge per residue at neutral pH."""
    radius = _validate_radius(radius)

    def raw_values(sequence):
        return _scale_values(sequence, radius, NEUTRAL_PH_CHARGE)

    return _zscore_channel(raw_values, bg_mu, bg_sigma)


def mk_hydropathy(
    radius,
    bg_mu,
    bg_sigma,
    scale=KYTE_DOOLITTLE_HYDROPATHY,
):
    """Return a channel for mean local hydropathy."""
    radius = _validate_radius(radius)
    scale = _validate_scale(scale)

    def raw_values(sequence):
        return _scale_values(sequence, radius, scale)

    return _zscore_channel(raw_values, bg_mu, bg_sigma)


def _validate_angle(angle):
    if (
        isinstance(angle, bool)
        or not isinstance(angle, Real)
        or not math.isfinite(angle)
        or angle <= 0
        or angle > 360
    ):
        raise ValueError("angle must be finite and in the interval (0, 360]")
    return math.radians(float(angle))


def _hydrophobic_moment(sequence, start, end, scale, angle_radians):
    x = 0.0
    y = 0.0
    for offset, amino_acid in enumerate(sequence[start:end]):
        hydrophobicity = scale[amino_acid]
        phase = offset * angle_radians
        x += hydrophobicity * math.cos(phase)
        y += hydrophobicity * math.sin(phase)
    return math.hypot(x, y) / (end - start)


def _hydrophobic_moment_values(sequence, radius, angle_radians, scale):
    return _window_values(
        sequence,
        radius,
        lambda seq, start, end: _hydrophobic_moment(
            seq, start, end, scale, angle_radians
        ),
    )


def mk_hydrophobic_moment(
    radius,
    bg_mu,
    bg_sigma,
    angle=100.0,
    scale=EISENBERG_HYDROPHOBICITY,
):
    """Return an Eisenberg-style local hydrophobic-moment channel."""
    radius = _validate_radius(radius)
    angle_radians = _validate_angle(angle)
    scale = _validate_scale(scale)

    def raw_values(sequence):
        return _hydrophobic_moment_values(sequence, radius, angle_radians, scale)

    return _zscore_channel(raw_values, bg_mu, bg_sigma)


def mk_disorder_propensity(
    radius,
    bg_mu,
    bg_sigma,
    scale=TOP_IDP_DISORDER_PROPENSITY,
):
    """Return a channel for mean local TOP-IDP disorder propensity."""
    radius = _validate_radius(radius)
    scale = _validate_scale(scale)

    def raw_values(sequence):
        return _scale_values(sequence, radius, scale)

    return _zscore_channel(raw_values, bg_mu, bg_sigma)


def _sequence_entropy(sequence, start, end):
    length = end - start
    counts = Counter(sequence[start:end])
    return -sum(
        (count / length) * math.log2(count / length)
        for count in counts.values()
    )


def mk_sequence_entropy(radius, bg_mu, bg_sigma):
    """Return a channel for local Shannon entropy in bits."""
    radius = _validate_radius(radius)

    def raw_values(sequence):
        return _window_values(sequence, radius, _sequence_entropy)

    return _zscore_channel(raw_values, bg_mu, bg_sigma)


def _dipeptide_frequency(sequence, start, end, dipeptide):
    pair_count = end - start - 1
    return sum(
        sequence[index:index + 2] == dipeptide
        for index in range(start, end - 1)
    ) / pair_count


def _dipeptide_values(sequence, radius, dipeptide, *, allow_too_short=False):
    if len(sequence) < 2:
        if allow_too_short:
            return []
        raise ValueError("dipeptide frequency requires a sequence of at least two residues")
    return _window_values(
        sequence,
        radius,
        lambda seq, start, end: _dipeptide_frequency(seq, start, end, dipeptide),
    )


def mk_dipeptide_frequency(radius, dipeptide, bg_mu, bg_sigma):
    """Return a channel for the local frequency of one ordered dipeptide."""
    radius = _validate_radius(radius, allow_zero=False)
    dipeptide = _validate_dipeptide(dipeptide)

    def raw_values(sequence):
        return _dipeptide_values(sequence, radius, dipeptide)

    return _zscore_channel(raw_values, bg_mu, bg_sigma)


def _iter_fasta_sequences(fasta_path):
    path = os.fspath(fasta_path)
    found_record = False
    header = None
    sequence_lines = []
    with open(path, encoding="utf-8") as fasta:
        for line_number, raw_line in enumerate(fasta, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    if not sequence_lines:
                        raise ValueError(f"FASTA record {header!r} has no sequence")
                    yield _validate_sequence("".join(sequence_lines))
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"FASTA header on line {line_number} is empty")
                found_record = True
                sequence_lines = []
            else:
                if header is None:
                    raise ValueError(
                        f"FASTA sequence data precedes the first header on line {line_number}"
                    )
                sequence_lines.append(line)
    if header is not None:
        if not sequence_lines:
            raise ValueError(f"FASTA record {header!r} has no sequence")
        yield _validate_sequence("".join(sequence_lines))
    if not found_record:
        raise ValueError("FASTA file contains no records")


def _compute_background(fasta_path, raw_values: Callable[[str], Iterable[float]]):
    count = 0
    mean = 0.0
    sum_squared_deviations = 0.0
    for sequence in _iter_fasta_sequences(fasta_path):
        for value in raw_values(sequence):
            if not math.isfinite(value):
                raise ValueError("channel produced a nonfinite background value")
            count += 1
            difference = value - mean
            mean += difference / count
            sum_squared_deviations += difference * (value - mean)
    if count == 0:
        raise ValueError("FASTA file contains no usable positions")
    sigma = math.sqrt(max(0.0, sum_squared_deviations / count))
    if not math.isfinite(sigma) or sigma <= 0:
        raise ValueError("background standard deviation is zero or nonfinite")
    return mean, sigma


def compute_composition_background(fasta_path, radius, residues):
    radius = _validate_radius(radius)
    residues = _validate_residues(residues)

    return _compute_background(
        fasta_path,
        lambda sequence: _composition_values(sequence, radius, residues),
    )


def compute_net_charge_background(fasta_path, radius):
    radius = _validate_radius(radius)
    return _compute_background(
        fasta_path,
        lambda sequence: _scale_values(sequence, radius, NEUTRAL_PH_CHARGE),
    )


def compute_hydropathy_background(
    fasta_path,
    radius,
    scale=KYTE_DOOLITTLE_HYDROPATHY,
):
    radius = _validate_radius(radius)
    scale = _validate_scale(scale)
    return _compute_background(
        fasta_path,
        lambda sequence: _scale_values(sequence, radius, scale),
    )


def compute_hydrophobic_moment_background(
    fasta_path,
    radius,
    angle=100.0,
    scale=EISENBERG_HYDROPHOBICITY,
):
    radius = _validate_radius(radius)
    angle_radians = _validate_angle(angle)
    scale = _validate_scale(scale)
    return _compute_background(
        fasta_path,
        lambda sequence: _hydrophobic_moment_values(
            sequence, radius, angle_radians, scale
        ),
    )


def compute_disorder_background(
    fasta_path,
    radius,
    scale=TOP_IDP_DISORDER_PROPENSITY,
):
    radius = _validate_radius(radius)
    scale = _validate_scale(scale)
    return _compute_background(
        fasta_path,
        lambda sequence: _scale_values(sequence, radius, scale),
    )


def compute_entropy_background(fasta_path, radius):
    radius = _validate_radius(radius)
    return _compute_background(
        fasta_path,
        lambda sequence: _window_values(sequence, radius, _sequence_entropy),
    )


def compute_dipeptide_background(fasta_path, radius, dipeptide):
    radius = _validate_radius(radius, allow_zero=False)
    dipeptide = _validate_dipeptide(dipeptide)
    return _compute_background(
        fasta_path,
        lambda sequence: _dipeptide_values(
            sequence, radius, dipeptide, allow_too_short=True
        ),
    )
