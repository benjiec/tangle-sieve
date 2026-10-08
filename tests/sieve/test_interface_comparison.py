import unittest
import numpy as np
from sieve.interface_comparison import fit_coordinates, overlap_scores, residue_mapping, compare_models, contact_recurrence


def model(chains=('A','B'), offset=0):
    coords = [[0,0,0],[2,0,0],[0,2,0],[0,0,3]]
    return {'model':0,'chains':list(chains),'residues':[
        {'chain':c,'number':i+1+offset,'aa':aa,'ca':[x,y,z+ci*5]}
        for ci,c in enumerate(chains) for i,(aa,(x,y,z)) in enumerate(zip('AKDW',coords))],
        'contacts':[{'residue_1':1+offset,'residue_2':2+offset}], 'interactions':[]}


class ComparisonTests(unittest.TestCase):
    def test_set_metrics_and_empty(self):
        s=overlap_scores({(1,2),(2,3)},{(1,2)})
        self.assertEqual((s['precision'],s['recall'],s['f1']),(1,.5,2/3))
        self.assertIsNone(overlap_scores(set(),set())['f1'])
        self.assertIsNone(overlap_scores({(1,2)},set())['f1'])
        self.assertEqual(overlap_scores({(1,2)},{(4,5)})['f1'],0)

    def test_mapping_renumbered_chains_and_rigid_transform(self):
        ref,qry=model(),model(('X','Y'),100)
        rot=np.array([[0,-1,0],[1,0,0],[0,0,1]])
        for r in qry['residues']:
            r['ca']=(np.array(r['ca'])@rot+[20,30,-8]).tolist()
        mappings=[residue_mapping(ref,qry,a,b) for a,b in zip(ref['chains'],qry['chains'])]
        result=compare_models(ref,qry,mappings)
        self.assertEqual(result['f1'],1)
        self.assertAlmostEqual(result['anchor_rmsd'],0)
        self.assertAlmostEqual(result['partner_rmsd_after_anchor_fit'],0)
        qry['residues'][4]['ca'][0]+=2
        result=compare_models(ref,qry,mappings)
        self.assertAlmostEqual(result['partner_rmsd_after_anchor_fit'],1)

    def test_unmapped_excluded_not_false_mismatch(self):
        r,q=model(),model()
        mappings=[residue_mapping(r,q,c,c) for c in r['chains']]
        q['contacts'].append({'residue_1':999,'residue_2':2})
        result=compare_models(r,q,mappings)
        self.assertEqual(result['f1'],1)
        self.assertEqual(result['query_unmapped_contacts'],1)

    def test_degenerate_fit_and_reflection(self):
        self.assertIsNone(fit_coordinates([[0,0,0],[1,0,0]],[[0,0,0],[1,0,0]]))
        self.assertIsNone(fit_coordinates([[0,0,0],[1,0,0],[2,0,0]],[[0,0,0],[1,0,0],[2,0,0]]))
        fixed=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
        fitted=fit_coordinates(fixed,fixed*[-1,1,1])
        self.assertGreater(fitted[2],0)
        self.assertAlmostEqual(np.linalg.det(fitted[0]),1)

    def test_recurrence_counts_models_not_duplicate_rows(self):
        m=model()
        m['contacts']*=2
        other=model(); other['contacts']=[]
        self.assertEqual(contact_recurrence([m,other])[0]['fraction'],.5)

if __name__=='__main__': unittest.main()
