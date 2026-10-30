"""Stable patch membership identity, independent of detection settings and paths."""
from sieve.alphafold_composition import chain_compositions, hash_payload
from sieve.alphafold_structure import canonical, digest, ASSEMBLY_FORMAT


def model_identity(cif, data, seed):
    if '_assembly' in data:
        return [ASSEMBLY_FORMAT, data['_assembly']['identity']]
    if seed is not None:
        try:
            compositions = chain_compositions(cif, data.get('_job_request'))
        except (ValueError, KeyError):
            if '_entity_poly_seq.' in cif or data.get('_job_request', {}).get('sequences'):
                raise
        else:
            return ['alphafold', seed, hash_payload(compositions.values())]
    # Bare CIFs have no seeded fold identity; retain exact content identity.
    return ['cif', seed, digest(cif)]


def patch_identifier(identity, model, sides):
    members = sorted((chain, sorted(set(positions))) for chain, positions in sides)
    return digest(canonical(['sieve-patch-v1', identity, model, members]))[:32]


def attach_patch_identity(row, identity, model):
    row['patch_id'] = patch_identifier(identity, model,
        [(row[f'chain {i}'], [int(n) for n in row[f'interface residues {i}'].split(',') if n]) for i in (1,2)])
