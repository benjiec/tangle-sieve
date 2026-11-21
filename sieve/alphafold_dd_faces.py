"""Classify pairs of death-family faces."""


def candidate_interface(face_1, face_2):
    """Classify complementary heuristic faces; None means unavailable."""
    if face_1 in (None, 'NA') or face_2 in (None, 'NA'):
        return None
    if 'none' in (face_1, face_2):
        return 'None'
    sides = [dict((label.removesuffix(' (partial)'), label.endswith(' (partial)'))
                  for label in faces.split('; ')) for faces in (face_1, face_2)]
    variant = False
    for a, b in (('Ia', 'Ib'), ('IIa', 'IIb'), ('IIIa', 'IIIb')):
        for first, second in ((a, b), (b, a)):
            if first in sides[0] and second in sides[1]:
                if not sides[0][first] and not sides[1][second]:
                    return 'Canonical'
                variant = True
    return 'Variant' if variant else 'Atypical'
