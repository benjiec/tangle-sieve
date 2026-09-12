"""Cached genomic background windows for Gimme FPR calibration."""
import bisect
import hashlib
import os
import random
import tempfile

import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser


def valid_start_intervals(sequence, length):
    """Half-open ranges of starts whose windows contain at most 10% Ns."""
    if length < 1:
        raise ValueError("Background length must be positive")
    if len(sequence) < length:
        return []
    bases = np.frombuffer(sequence.upper().encode("ascii"), dtype=np.uint8)
    prefix = np.empty(len(bases) + 1, dtype=np.int64)
    prefix[0] = 0
    np.cumsum(bases == ord("N"), out=prefix[1:])
    valid = prefix[length:] - prefix[:-length] <= length // 10
    edges = np.flatnonzero(valid[1:] != valid[:-1]) + 1
    boundaries = [0, *edges.tolist(), len(valid)]
    return [(start, end) for start, end in zip(boundaries, boundaries[1:]) if valid[start]]


def genomic_background(source, length, nseq=10000):
    """Uniformly sample valid genomic starts with replacement, using seed zero.

    Cache identity includes source metadata, window length, sample count and
    sampler version (which fixes the N limit and seed). No genomepy dependency.
    """
    if type(length) is not int or length < 1 or type(nseq) is not int or nseq < 1:
        raise ValueError("Background length and sample count must be positive integers")
    source = os.path.realpath(source)
    stat = os.stat(source)
    identity = f"v1:{source}:{stat.st_size}:{stat.st_mtime_ns}:{length}:{nseq}"
    digest = hashlib.sha256(identity.encode()).hexdigest()
    root = os.environ.get("XDG_CACHE_HOME") or os.path.expanduser("~/.cache")
    directory = os.path.join(root, "sieve", "gimme-backgrounds", digest)
    destination = os.path.join(directory, "background.fna")
    if os.path.isfile(destination):
        return destination
    intervals, cumulative = [], []
    total = 0
    with open(source) as fasta:
        for index, (_title, sequence) in enumerate(SimpleFastaParser(fasta)):
            for start, end in valid_start_intervals(sequence, length):
                total += end - start
                intervals.append((index, start))
                cumulative.append(total)
    if not total:
        raise ValueError(f"No {length}-base background windows with at most 10% Ns in {source}")
    rng = random.Random(0)
    selected = {}
    for sample in range(nseq):
        position = rng.randrange(total)
        interval = bisect.bisect_right(cumulative, position)
        index, start = intervals[interval]
        offset = position - (cumulative[interval - 1] if interval else 0)
        selected.setdefault(index, []).append((sample, start + offset))
    os.makedirs(directory, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=directory, delete=False) as output:
            temporary = output.name
            with open(source) as fasta:
                for index, (title, sequence) in enumerate(SimpleFastaParser(fasta)):
                    for sample, start in selected.get(index, []):
                        output.write(f">background_{sample} {title.split()[0]}:{start}-{start + length}\n")
                        output.write(sequence[start:start + length].upper() + "\n")
        os.replace(temporary, destination)
    finally:
        if temporary is not None and os.path.exists(temporary):
            os.unlink(temporary)
    return destination
