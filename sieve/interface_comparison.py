"""Transparent sequence-mapped interface comparisons; no fitted composite score."""
from collections import Counter
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices


def chain_rows(model, chain):
    return [r for r in model['residues'] if r['chain'] == chain]


def residue_mapping(reference, query, reference_chain, query_chain, gap_open=10):
    ref, qry = chain_rows(reference, reference_chain), chain_rows(query, query_chain)
    aligner = PairwiseAligner(mode='global')
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score = -gap_open
    aligner.extend_gap_score = -0.5
    alignments = iter(aligner.align(''.join(r['aa'] for r in ref), ''.join(r['aa'] for r in qry)))
    alignment = next(alignments)
    tied = next(alignments, None) is not None
    rows = []
    for rb, qb in zip(*alignment.aligned):
        for i, j in zip(range(*rb), range(*qb)):
            rows.append({'reference_chain': reference_chain, 'query_chain': query_chain,
                         'reference_residue': ref[i]['number'], 'query_residue': qry[j]['number'],
                         'reference_aa': ref[i]['aa'], 'query_aa': qry[j]['aa']})
    return {'rows': rows, 'gap_open': gap_open, 'gap_extend': 0.5,
            'multiple_optimal_alignments': tied, 'alignment': str(alignment),
            'aligned_pairs': len(rows), 'reference_length': len(ref), 'query_length': len(qry),
            'reference_coverage': len(rows)/len(ref), 'query_coverage': len(rows)/len(qry),
            'identity': sum(r['reference_aa'] == r['query_aa'] for r in rows)/len(rows) if rows else None}


def overlap_scores(reference, query):
    """Both sets must be expressed in reference numbering and restricted to mapped pairs."""
    shared = len(reference & query)
    return {'reference_contacts': len(reference), 'query_contacts': len(query), 'shared_contacts': shared,
            'precision': shared/len(query) if query else None,
            'recall': shared/len(reference) if reference else None,
            'f1': 2*shared/(len(reference)+len(query)) if reference and query else None}


def fit_coordinates(fixed, moving):
    fixed, moving = np.asarray(fixed, float), np.asarray(moving, float)
    if len(fixed) < 3 or np.linalg.matrix_rank(fixed-fixed.mean(axis=0)) < 2 or np.linalg.matrix_rank(moving-moving.mean(axis=0)) < 2:
        return None
    f, m = fixed.mean(axis=0), moving.mean(axis=0)
    u, _, vt = np.linalg.svd((moving-m).T @ (fixed-f))
    correction = np.eye(3)
    correction[-1,-1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    translation = f-m@rotation
    rmsd = float(np.sqrt(np.mean(np.sum((moving@rotation+translation-fixed)**2, axis=1))))
    return rotation, translation, rmsd


def compare_models(reference, query, mappings):
    maps = [{r['query_residue']: r['reference_residue'] for r in mapping['rows']} for mapping in mappings]
    eligible = [set(m.values()) for m in maps]
    rc = {(c['residue_1'],c['residue_2']) for c in reference['contacts']}
    qc = {(c['residue_1'],c['residue_2']) for c in query['contacts']}
    rs = {(a,b) for a,b in rc if a in eligible[0] and b in eligible[1]}
    qs = {(maps[0][a],maps[1][b]) for a,b in qc if a in maps[0] and b in maps[1]}
    result = {'reference_model': reference['model'], 'query_model': query['model'],
              **overlap_scores(rs,qs), 'reference_unmapped_contacts': len(rc)-len(rs),
              'query_unmapped_contacts': len(qc)-len(qs)}
    coordinates = []
    for mapping in mappings:
        ref = {r['number']:r['ca'] for r in chain_rows(reference,mapping['rows'][0]['reference_chain'])} if mapping['rows'] else {}
        qry = {r['number']:r['ca'] for r in chain_rows(query,mapping['rows'][0]['query_chain'])} if mapping['rows'] else {}
        coordinates.append((np.array([ref[r['reference_residue']] for r in mapping['rows']]),
                            np.array([qry[r['query_residue']] for r in mapping['rows']])))
    fit = fit_coordinates(*coordinates[0]) if len(coordinates[0][0]) else None
    result['anchor_rmsd'] = fit[2] if fit else None
    result['partner_rmsd_after_anchor_fit'] = float(np.sqrt(np.mean(np.sum((coordinates[1][1]@fit[0]+fit[1]-coordinates[1][0])**2,axis=1)))) if fit and len(coordinates[1][0]) else None
    partner_fit = fit_coordinates(*coordinates[1]) if len(coordinates[1][0]) else None
    result['partner_independent_rmsd'] = partner_fit[2] if partner_fit else None
    result['typed_contacts'] = {}
    for kind in ('salt_bridge_candidate','nonpolar_contact','steric_overlap_candidate'):
        rset = {(x['residue_1'],x['residue_2']) for x in reference['interactions'] if x['type']==kind and x['residue_1'] in eligible[0] and x['residue_2'] in eligible[1]}
        qset = {(maps[0][x['residue_1']],maps[1][x['residue_2']]) for x in query['interactions'] if x['type']==kind and x['residue_1'] in maps[0] and x['residue_2'] in maps[1]}
        result['typed_contacts'][kind] = overlap_scores(rset,qset)
    return result


def contact_recurrence(models):
    counts = Counter(pair for m in models for pair in {(x['residue_1'],x['residue_2']) for x in m['contacts']})
    return [{'residue_1':a,'residue_2':b,'models':count,'fraction':count/len(models)} for (a,b),count in sorted(counts.items())]


def compare_ensembles(reference, query):
    variants = []
    for gap in (10,8,12):
        mappings = [residue_mapping(reference[0],query[0],rc,qc,gap) for rc,qc in zip(reference[0]['chains'],query[0]['chains'])]
        variants.append({'gap_open':gap, 'mappings': mappings,
                         'comparisons':[compare_models(r,q,mappings) for r in reference for q in query]})
    return {'method':'sequence-mapped contact overlap, version 1', 'primary':variants[0],
            'mapping_sensitivity':variants[1:], 'reference_recurrence':contact_recurrence(reference),
            'query_recurrence':contact_recurrence(query)}
