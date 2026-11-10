"""Optional geometric normalization of DSSP H segments for face classification."""
from copy import deepcopy
import math
import numpy as np


def validate_kink_settings(connector_length, angle):
    if isinstance(connector_length, bool) or not isinstance(connector_length, int) or connector_length < 1:
        raise ValueError('kink connector length must be a positive integer')
    if not math.isfinite(angle) or not 0 < angle <= 180:
        raise ValueError('kink angle must be finite and in (0, 180]')


def local_axis(points):
    """Directed local helical axis from four consecutive C-alpha positions."""
    p = np.asarray(points, dtype=float)
    if p.shape != (4, 3) or not np.isfinite(p).all():
        raise ValueError('requires four finite Cα coordinates')
    d = np.diff(p, axis=0)
    axis = np.cross(d[1] - d[0], d[2] - d[1])
    norm = np.linalg.norm(axis)
    if norm < 1e-10:
        raise ValueError('degenerate local helix axis')
    axis /= norm
    direction = np.dot(axis, p[-1] - p[0])
    if abs(direction) < 1e-10:
        raise ValueError('indeterminate N→C axis direction')
    return axis if direction > 0 else -axis


def normalize_helices(model, chain, codes, segments, *, connector_length=1, angle=70):
    """Evaluate original neighbors once; merge connected qualifying H segments.

    Returns copied segments for this chain and human-readable decisions. Raw DSSP
    codes and the caller's segments are never modified.
    """
    validate_kink_settings(connector_length, angle)
    original = sorted((s for s in segments if s['chain'] == chain), key=lambda s:s['start'])
    helices = [s for s in original if s['code'] == 'H']
    accepted = set(); reports = []
    def axis(segment, last):
        positions = list(segment['residues'])
        if len(positions) < 4:
            raise ValueError('fewer than four helix residues')
        positions = positions[-4:] if last else positions[:4]
        try:
            points = [model[chain][(' ', n, ' ')]['CA'].coord for n in positions]
        except (KeyError, TypeError):
            raise ValueError('missing terminal helix Cα coordinates') from None
        return local_axis(points)
    for i, (left, right) in enumerate(zip(helices, helices[1:])):
        gap = list(range(left['end'] + 1, right['start']))
        if not 1 <= len(gap) <= connector_length:
            continue
        if any(codes.get((chain,n)) not in ('-', 'T', 'S') for n in gap):
            continue
        label = f"{left['start']}–{left['end']} + {right['start']}–{right['end']}"
        try:
            degrees = float(np.degrees(np.arccos(np.clip(np.dot(axis(left,True),axis(right,False)),-1,1))))
        except ValueError as exc:
            reports.append(f'not merged {label}: {exc}')
            continue
        if degrees < angle:
            accepted.add(i)
            reports.append(f'merged {label}: {degrees:.3f}°')
    merged = []; removed = set(); i = 0
    while i < len(helices):
        first = i
        while i in accepted:
            i += 1
        segment = deepcopy(helices[first])
        if i != first:
            segment['end'] = helices[i]['end']
            segment['residues'] = list(range(segment['start'],segment['end']+1))
            removed.update(range(segment['start'],segment['end']+1))
        merged.append(segment)
        i += 1
    result = merged + [deepcopy(s) for s in original if s['code'] != 'H' and not set(s['residues']) & removed]
    result.sort(key=lambda s:s['start'])
    counts = {}
    for s in result:
        prefix = 'alpha' if s['code'] == 'H' else ('loop' if s['code'] in ('-', 'T', 'S') else None)
        if prefix:
            counts[prefix] = counts.get(prefix,0)+1
            s['label'] = f'{prefix}{counts[prefix]}'
    return result, '; '.join(reports)
