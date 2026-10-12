"""Validation shared by AlphaFold confidence readers."""
import math


def model_confidence(summary):
    if not isinstance(summary, dict):
        raise ValueError('summary confidence JSON must be an object')
    scores = {}
    for field in ('iptm', 'ptm'):
        value = summary.get(field)
        if value is not None and (isinstance(value, bool) or not isinstance(value, (int, float))
                                  or not math.isfinite(value) or not 0 <= value <= 1):
            raise ValueError(f'summary confidence {field} must be a finite number in [0, 1] or null')
        scores[field] = value
    return scores

