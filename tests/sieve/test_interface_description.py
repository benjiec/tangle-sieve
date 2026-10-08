import unittest
import numpy as np
from Bio.PDB import Model, Chain, Residue, Atom
from sieve.interface_description import DescriptionOptions, composition, neighborhood, atom_interactions, sasa_burial, secondary_structure, describe_model
from tests.sieve.test_alphafold_pdockq2_patches import CIF, DATA


def make_model():
    model=Model.Model(0)
    for c,name,x in [('A','LYS',0),('B','ASP',4)]:
        chain=Chain.Chain(c); model.add(chain)
        r=Residue.Residue((' ',1,' '),name,' '); chain.add(r)
        for atom,coord,element in [('CA',[x,0,0],'C'),('NZ' if c=='A' else 'OD1',[x,2,0],'N' if c=='A' else 'O')]:
            r.add(Atom.Atom(atom,np.array(coord,float),90,1,' ',atom,1,element=element))
    return model


class DescriptorTests(unittest.TestCase):
    def test_charge_and_empty(self):
        m=make_model()
        self.assertEqual(composition(m['A'])['nominal_sidechain_charge'],1)
        self.assertEqual(composition(m['B'])['nominal_sidechain_charge'],-1)
        self.assertIsNone(composition([])['mean_hydropathy'])

    def test_salt_boundary_and_atom_detail(self):
        m=make_model()
        interactions=atom_interactions(m,('A','B'),DescriptionOptions())
        salt=[x for x in interactions if x['type']=='salt_bridge_candidate']
        self.assertEqual(len(salt),1)
        self.assertEqual((salt[0]['atom_1'],salt[0]['atom_2']),('NZ','OD1'))
        m['B'][1]['OD1'].coord[0]+=.01
        self.assertFalse(any(x['type']=='salt_bridge_candidate' for x in atom_interactions(m,('A','B'),DescriptionOptions())))

    def test_discontinuous_shell_deduplicated(self):
        m=make_model(); chain=m['A']
        for n in range(2,6):
            r=Residue.Residue((' ',n,' '),'ALA',' '); chain.add(r)
            r.add(Atom.Atom('CA',np.array([n,0,0],float),90,1,' ','CA',n,element='C'))
        d=neighborhood(chain,{2,4},1,8)
        self.assertEqual(d['sequence_shell']['residues'],[1,3,5])
        self.assertEqual(d['spatial_shell']['residues'],[1,3,5])
        self.assertEqual(neighborhood(chain,set(),1,8)['spatial_shell']['count'],0)

    def test_sasa_separation(self):
        m=make_model(); b=sasa_burial(m,('A','B'),120)
        self.assertGreater(b[('A',1)]['buried_area'],0)
        for a in m['B'].get_atoms(): a.coord+=100
        b=sasa_burial(m,('A','B'),120)
        self.assertAlmostEqual(b[('A',1)]['buried_area'],0)

    def test_invalid_options_and_explicit_no_dssp(self):
        for kwargs in ({'cutoff':0},{'radius':float('nan')},{'flank':-1},{'sasa_points':1},{'flank':1.5}):
            with self.assertRaises(ValueError): DescriptionOptions(**kwargs)
        self.assertEqual(secondary_structure(make_model(),'',None)[2]['status'],'unavailable')

    def test_empty_contacts_and_reversed_chain_order(self):
        result=describe_model(CIF,DATA,('B','A'),DescriptionOptions(cutoff=1),None)
        self.assertEqual(result['contacts'],[])
        self.assertEqual(result['patches'][0]['description']['B']['patch']['count'],0)
        result=describe_model(CIF,DATA,('B','A'),DescriptionOptions(),None)
        self.assertEqual(result['contacts'][0]['pae_1_to_2'],10)
        self.assertEqual(result['contacts'][0]['pae_2_to_1'],2)

    def test_invalid_chain_and_pae(self):
        with self.assertRaises(ValueError): describe_model(CIF,DATA,('A','Z'),DescriptionOptions(),None)
        for value in (-1,float('nan'),float('inf')):
            bad=dict(DATA,pae=[[0,value],[10,0]])
            with self.assertRaises(ValueError): describe_model(CIF,bad,('A','B'),DescriptionOptions(),None)

if __name__=='__main__': unittest.main()
