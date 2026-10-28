"""Provisional DD face candidates from six sequential DSSP H segments.

Structural elements follow https://pmc.ncbi.nlm.nih.gov/articles/PMC2908332/.
Coverage is a screening heuristic, not a canonical topology assignment.
"""

FACES = {
    'Ia': ('H1', 'H4'), 'Ib': ('H2', 'H3'),
    'IIa': ('H4', 'L45'), 'IIb': ('L56', 'H6'),
    'IIIa': ('H3',), 'IIIb': ('L12', 'L34'),
}


def validate_threshold(threshold, name='helix'):
    if isinstance(threshold, bool) or not isinstance(threshold, int) or threshold < 1:
        raise ValueError(f'--dd-{name}-presence-threshold must be a positive integer residue count')


def candidate_faces(chain, selected, codes, segments, threshold=3, loop_threshold=1):
    validate_threshold(threshold)
    validate_threshold(loop_threshold, 'loop')
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
        counts = [len(selected & elements[e]) for e in required]
        passed = [count >= (threshold if e.startswith('H') else loop_threshold)
                  for e, count in zip(required, counts)]
        if all(counts) and any(passed):
            partial = not all(passed)
            label = f'{face} (partial)' if partial else face
            ranked.append((partial, -min(coverage), -sum(coverage)/len(coverage), order, label))
    return '; '.join(row[-1] for row in sorted(ranked)) or 'none'
