"""Provisional DD face candidates from six sequential DSSP H segments.

Structural elements follow https://pmc.ncbi.nlm.nih.gov/articles/PMC2908332/.
Coverage is a screening heuristic, not a canonical topology assignment.
"""
import math

FACES = {
    'Ia': ('H1', 'H4'), 'Ib': ('H2', 'H3'),
    'IIa': ('H4', 'L45'), 'IIb': ('L56', 'H6'),
    'IIIa': ('H3',), 'IIIb': ('L12', 'L34'),
}


def validate_threshold(threshold):
    if not math.isfinite(threshold) or not 0 < threshold <= 1:
        raise ValueError('--dd-helix-presence-threshold must be greater than 0 and at most 1')


def candidate_faces(chain, selected, codes, segments, threshold=.5):
    validate_threshold(threshold)
    helices = sorted((s for s in segments if s['chain'] == chain and s['code'] == 'H'),
                     key=lambda s: s['start'])
    if len(helices) != 6:
        return None
    if any(codes.get((chain, n), '?') == '?' for n in range(helices[0]['start'], helices[-1]['end']+1)):
        return None
    elements = {f'H{i}': set(s['residues']) for i,s in enumerate(helices,1)}
    for i in range(1,6):
        interval = set(range(helices[i-1]['end']+1, helices[i]['start']))
        # Do not call another structured segment (e.g. G or E) a canonical loop.
        elements[f'L{i}{i+1}'] = interval if all(codes.get((chain,n)) in ('-', 'T', 'S') for n in interval) else set()
    ranked = []
    for order,(face,required) in enumerate(FACES.items()):
        if any(not elements[e] for e in required):
            continue
        coverage = [len(selected & elements[e])/len(elements[e]) for e in required]
        if min(coverage) >= threshold:
            ranked.append((-min(coverage), -sum(coverage)/len(coverage), order, face))
    return '; '.join(row[-1] for row in sorted(ranked)) or 'none'
