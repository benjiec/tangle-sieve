import copy
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from sieve.alphafold_helix_normalization import normalize_helices, local_axis, validate_kink_settings


def fixture(count=2, connector='T'):
    model={'A':{}}; codes={}; segments=[]
    for i in range(count):
        start=1+i*5
        residues=list(range(start,start+4))
        segments.append(dict(chain='A',code='H',start=start,end=start+3,residues=residues,label=f'alpha{i+1}'))
        for j,n in enumerate(residues):
            theta=j*np.deg2rad(100)
            model['A'][(' ',n,' ')]={'CA':SimpleNamespace(coord=np.array([np.cos(theta),np.sin(theta),j*1.5]))}
            codes['A',n]='H'
        if i<count-1:
            n=start+4;codes['A',n]=connector
            segments.append(dict(chain='A',code=connector,start=n,end=n,residues=[n],label='loop1'))
    return model,codes,segments


class NormalizationTests(unittest.TestCase):
    def test_merge_copy_chained(self):
        for count in (2,3):
            m,c,s=fixture(count); before=copy.deepcopy(s)
            out,report=normalize_helices(m,'A',c,s)
            self.assertEqual(len(out),1)
            self.assertEqual(out[0]['residues'],list(range(1,5*count)))
            self.assertEqual(report.count('merged'),count-1)
            self.assertEqual(s,before)
            self.assertEqual(c['A',5],'T')

    def test_connector_codes(self):
        for code in ('-','T','S','G','I','E','B','P','?'):
            m,c,s=fixture(connector=code)
            out,_=normalize_helices(m,'A',c,s)
            self.assertEqual(sum(x['code']=='H' for x in out),1 if code in '-TS' else 2)
        del c['A',5]
        self.assertEqual(len(normalize_helices(m,'A',c,s)[0]),3)

    def test_connector_length(self):
        m,c,s=fixture()
        for n in range(9,5,-1):
            m['A'][(' ',n+1,' ')]=m['A'].pop((' ',n,' '))
        s[2].update(start=7,end=10,residues=[7,8,9,10])
        c['A',6]='T'
        self.assertEqual(len(normalize_helices(m,'A',c,s)[0]),3)
        self.assertEqual(len(normalize_helices(m,'A',c,s,connector_length=2)[0]),1)

    def test_angle_strict(self):
        m,c,s=fixture()
        for angle,merged in ((69,True),(70,False),(130,False)):
            v=np.array([np.sin(np.deg2rad(angle)),0,np.cos(np.deg2rad(angle))])
            measured=float(np.degrees(np.arccos(v[2])))
            with patch('sieve.alphafold_helix_normalization.local_axis',side_effect=[np.array([0,0,1]),v]):
                out,_=normalize_helices(m,'A',c,s,angle=70 if angle!=70 else measured)
            self.assertEqual(len(out)==1,merged)

    def test_geometry_failure(self):
        for kind in ('missing','nonfinite','degenerate','short'):
            m,c,s=fixture()
            if kind=='missing': del m['A'][(' ',1,' ')]
            elif kind=='nonfinite': m['A'][(' ',1,' ')]['CA'].coord[:]=np.nan
            elif kind=='short': s[0]['residues']=[1,2,3]
            else:
                for n in range(1,5):m['A'][(' ',n,' ')]['CA'].coord[:]=n
            out,report=normalize_helices(m,'A',c,s)
            self.assertIn('not merged',report)
            self.assertEqual(len(out),3)

    def test_rigid_transform(self):
        m,c,s=fixture(); p=np.array([m['A'][(' ',n,' ')]['CA'].coord for n in range(1,5)])
        rotation=np.array([[0,-1,0],[1,0,0],[0,0,1]])
        np.testing.assert_allclose(local_axis(p@rotation.T+17),rotation@local_axis(p),atol=1e-12)

    def test_validation(self):
        for length in (0,-1,1.5,True):
            with self.assertRaises(ValueError):validate_kink_settings(length,70)
        for angle in (0,-1,181,float('inf'),float('nan')):
            with self.assertRaises(ValueError):validate_kink_settings(1,angle)

    def test_renumber_and_eligibility(self):
        from sieve.alphafold_death_rules import elements_for_chain
        m,c,s=fixture(7)
        for n in (10,15,20,25,30):c['A',n]='G'
        s.append(dict(chain='B',code='H',start=1,end=4,residues=[1,2,3,4],label='alpha1'))
        out,_=normalize_helices(m,'A',c,s)
        self.assertTrue(all(x['chain']=='A' for x in out))
        elements,status=elements_for_chain('A',c,out)
        self.assertEqual(status,'eligible')
        self.assertEqual(elements['alpha1'],set(range(1,10)))
        self.assertEqual(elements['alpha6'],set(range(31,35)))
        m,c,s=fixture(6)
        out,_=normalize_helices(m,'A',c,s)
        self.assertIsNone(elements_for_chain('A',c,out)[0])
