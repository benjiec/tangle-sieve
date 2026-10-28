import unittest
from sieve.alphafold_dd_faces import candidate_faces, candidate_interface, validate_threshold


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
    def test_interface_complements_both_orientations(self):
        for a,b in (('Ia','Ib'),('IIa','IIb'),('IIIa','IIIb')):
            for first,second in ((a,b),(b,a)):
                for p,q in ((False,False),(True,False),(False,True),(True,True)):
                    with self.subTest(first=first,second=second,p=p,q=q):
                        self.assertEqual(candidate_interface(
                            first+(' (partial)' if p else ''),
                            second+(' (partial)' if q else '')),
                            'Variant' if p or q else 'Canonical')

    def test_interface_multiple_candidates_and_precedence(self):
        for a,b,expected in (
                ('Ia (partial); IIa','Ib; IIb','Canonical'),
                ('IIIa; Ia','IIb; Ib (partial)','Variant'),
                ('Ia (partial)','IIb','Atypical'),
                ('Ia','Ia','Atypical'),
                ('Ia; IIIa','IIa; IIb','Atypical'),
                ('none','Ia','None'), ('none','none','None'),
                (None,'Ia',None), ('NA','Ia',None),
                (None,'none',None), ('NA','none',None), (None,None,None)):
            for first,second in ((a,b),(b,a)):
                with self.subTest(first=first,second=second):
                    self.assertEqual(candidate_interface(first,second),expected)

    def test_exact_six_and_missing_assignments(self):
        for n in (0,5,7):
            c,s=structure(n);self.assertIsNone(candidate_faces('A',set(range(1,50)),c,s))
        c,s=structure();c.pop(('A',5))
        self.assertIsNone(candidate_faces('A',set(range(1,50)),c,s))

    def test_threshold_boundary_ranking_and_ties(self):
        c,s=structure()
        selected={1,2,3,19,20,21,*range(13,17)}
        self.assertEqual(candidate_faces('A',selected,c,s),'IIIa; Ia; Ib (partial); IIa (partial)')
        self.assertEqual(candidate_faces('A',selected,c,s,4),'IIIa; Ia (partial); Ib (partial); IIa (partial)')
        self.assertEqual(candidate_faces('A',set(range(1,35)),c,s,1),'Ia; Ib; IIa; IIb; IIIa; IIIb')
        self.assertEqual(candidate_faces('A',set(),c,s),'none')
        self.assertIsNone(candidate_faces('B',selected,c,s))

    def test_loop_faces_and_no_structured_loop_substitution(self):
        c,s=structure()
        self.assertEqual(candidate_faces('A',{5,17},c,s),'IIIb')
        c[('A',5)]='G'
        self.assertEqual(candidate_faces('A',{5,17},c,s),'IIIb')
        c[('A',6)]='G'
        self.assertEqual(candidate_faces('A',{5,17},c,s),'IIIb')
        self.assertEqual(candidate_faces('A',{5,17,18},c,s,3,2),'IIIb (partial)')
        self.assertEqual(candidate_faces('A',{17,18},c,s),'IIIb (partial)')
        for code in ('E','B','I','P'):
            c[('A',6)]=code
            self.assertEqual(candidate_faces('A',{5,17},c,s),'IIIb (partial)')
        c,s=structure()
        self.assertEqual(candidate_faces('A',{19,20,21,23},c,s),'IIa; Ia (partial)')
        self.assertEqual(candidate_faces('A',{29,31,32,33},c,s),'IIb')

    def test_reported_chain_with_310_connectors(self):
        codes={('B',n):'-' for n in range(1,108)}
        segments=[]
        for start,end,code in ((4,14,'H'),(15,17,'G'),(18,24,'H'),
                               (29,38,'H'),(42,57,'H'),(58,60,'G'),
                               (63,72,'H'),(76,105,'H')):
            residues=list(range(start,end+1))
            segments.append(dict(chain='B',code=code,start=start,end=end,residues=residues))
            codes.update({('B',n):code for n in residues})
        selected={11,12,13,14,15,16,40,41,42,43}
        self.assertEqual(candidate_faces('B',selected,codes,segments),'IIIb; Ia (partial); IIa (partial)')

    def test_invalid_threshold(self):
        for x in (0,-1,.5,3.0,True,float('nan'),float('inf')):
            with self.assertRaises(ValueError):validate_threshold(x)
            with self.assertRaises(ValueError):validate_threshold(x,'loop')

    def test_independent_loop_count_and_short_elements(self):
        c,s=structure()
        self.assertEqual(candidate_faces('A',{19,20,21,23},c,s), 'IIa; Ia (partial)')
        self.assertEqual(candidate_faces('A',{19,20,21,23},c,s,3,2), 'IIa (partial); Ia (partial)')
        self.assertEqual(candidate_faces('A',{19,20,21,23,24},c,s,3,2), 'IIa; Ia (partial)')
        self.assertEqual(candidate_faces('A',{5,6,17,18},c,s,3,3), 'IIIb (partial)')
        self.assertEqual(candidate_faces('A',set(range(1,35)),c,s,5,1), 'IIIb; Ia (partial); Ib (partial); IIa (partial); IIb (partial); IIIa (partial)')

    def test_partial_helix_and_full_before_partial(self):
        c,s=structure()
        self.assertEqual(candidate_faces('A',{1,2,3,19},c,s), 'Ia (partial); IIa (partial)')
        # Full IIb has lower fractional coverage than partial Ia, but ranks first.
        selected={1,2,3,4,19,20,29,31,32,33}
        self.assertEqual(candidate_faces('A',selected,c,s), 'IIb; Ia (partial); IIa (partial)')

    def test_partial_requires_any_contact(self):
        c,s=structure()
        for selected,expected in (({1}, 'Ia (partial)'),
                                  ({1,2,3}, 'Ia (partial)'),
                                  ({1,2,19,20}, 'Ia (partial); IIa (partial)'),
                                  (set(), 'none'), ({11,12}, 'none'),
                                  ({99}, 'none')):
            self.assertEqual(candidate_faces('A',selected,c,s), expected)

    def test_partial_loops_and_single_element(self):
        c,s=structure()
        self.assertEqual(candidate_faces('A',{5,6,17},c,s,3,2), 'IIIb (partial)')
        self.assertEqual(candidate_faces('A',{5,17},c,s,3,2), 'IIIb (partial)')
        self.assertEqual(candidate_faces('A',{13,14},c,s), 'IIIa (partial); Ib (partial)')
        self.assertEqual(candidate_faces('A',{13,14,15},c,s), 'IIIa; Ib (partial)')

    def test_partial_ties_and_threshold_one(self):
        c,s=structure()
        selected={1,2,3,19,7,8,9,13}
        self.assertEqual(candidate_faces('A',selected,c,s), 'Ia (partial); Ib (partial); IIIa (partial); IIa (partial)')
        self.assertEqual(candidate_faces('A',{1,19},c,s,1), 'Ia; IIa (partial)')
