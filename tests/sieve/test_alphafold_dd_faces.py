import unittest
from sieve.alphafold_dd_faces import candidate_interface


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
