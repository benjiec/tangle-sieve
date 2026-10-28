import unittest
from sieve.alphafold_dd_faces import candidate_faces, validate_threshold


def structure(count=6):
    segments=[];codes={}
    for i in range(count):
        start=1+i*6
        residues=list(range(start,start+4))
        segments.append(dict(chain='A',code='H',label=f'alpha{i+1}',start=start,end=start+3,residues=residues))
        codes.update({('A',n):'H' for n in residues})
        if i<count-1:codes.update({('A',n):'-' for n in (start+4,start+5)})
    return codes,segments


class DDFaceTests(unittest.TestCase):
    def test_exact_six_and_missing_assignments(self):
        for n in (0,5,7):
            c,s=structure(n);self.assertIsNone(candidate_faces('A',set(range(1,50)),c,s))
        c,s=structure();c.pop(('A',5))
        self.assertIsNone(candidate_faces('A',set(range(1,50)),c,s))

    def test_threshold_boundary_ranking_and_ties(self):
        c,s=structure()
        selected={1,2,3,19,20,21,*range(13,17)}
        self.assertEqual(candidate_faces('A',selected,c,s),'IIIa; Ia')
        self.assertEqual(candidate_faces('A',selected,c,s,4),'IIIa')
        self.assertEqual(candidate_faces('A',set(range(1,35)),c,s,1),'Ia; Ib; IIa; IIb; IIIa; IIIb')
        self.assertEqual(candidate_faces('A',set(),c,s),'none')
        self.assertIsNone(candidate_faces('B',selected,c,s))

    def test_loop_faces_and_no_structured_loop_substitution(self):
        c,s=structure()
        self.assertEqual(candidate_faces('A',{5,17},c,s),'IIIb')
        c[('A',5)]='G'
        self.assertEqual(candidate_faces('A',{5,17},c,s),'none')
        c,s=structure()
        self.assertEqual(candidate_faces('A',{19,20,21,23},c,s),'IIa')
        self.assertEqual(candidate_faces('A',{29,31,32,33},c,s),'IIb')

    def test_invalid_threshold(self):
        for x in (0,-1,.5,3.0,True,float('nan'),float('inf')):
            with self.assertRaises(ValueError):validate_threshold(x)
            with self.assertRaises(ValueError):validate_threshold(x,'loop')

    def test_independent_loop_count_and_short_elements(self):
        c,s=structure()
        self.assertEqual(candidate_faces('A',{19,20,21,23},c,s), 'IIa')
        self.assertEqual(candidate_faces('A',{19,20,21,23},c,s,3,2), 'none')
        self.assertEqual(candidate_faces('A',{19,20,21,23,24},c,s,3,2), 'IIa')
        self.assertEqual(candidate_faces('A',{5,6,17,18},c,s,3,3), 'none')
        self.assertEqual(candidate_faces('A',set(range(1,35)),c,s,5,1), 'IIIb')
