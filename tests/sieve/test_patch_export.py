import io
import unittest
from unittest.mock import patch

from sieve.alphafold_pdockq2_patches import PATCH_COLUMNS, Options, detect_patches, patch_export_rows
from sieve.alphafold_pdockq2 import parse_model
from sieve.patch_description_summary import PatchDescriptions
from tests.sieve.test_alphafold_pdockq2_patches import CIF, DATA


class ExportTests(unittest.TestCase):
    def test_metrics_strict_patch_scope_dedup_and_zero_area(self):
        descriptions=PatchDescriptions({('fold',0):CIF},None)
        residues,pae=parse_model(CIF,DATA)
        row=detect_patches(residues,pae,Options(min_residues_per_chain=1,min_contacts=1))[-1]
        row.update(source='fold',model=0)
        base=dict(type='salt_bridge_candidate',residue_1=1,residue_2=1,
                  aa_1='K',aa_2='D',atom_1='NZ',atom_2='OD1',distance=3.5)
        interactions=[base,dict(base,atom_2='OD2',distance=3),dict(base,residue_2=9,distance=2),
                      dict(base,type='nonpolar_contact'),dict(base,type='nonpolar_contact'),
                      dict(base,type='steric_overlap_candidate'),dict(base,type='phosphate_basic_contact')]
        hydrogen_bond = dict(base, type='hydrogen_bond_candidate')
        interactions.extend([
            hydrogen_bond, hydrogen_bond, dict(hydrogen_bond, atom_2='OD2'),
            dict(hydrogen_bond, residue_1=9), dict(hydrogen_bond, residue_2=9),
        ])
        burial={('A',1):dict(buried_area=20,buried_nonpolar_area=10),
                ('B',1):dict(buried_area=30,buried_nonpolar_area=5)}
        descriptions.pairs[('fold',0,'A','B')]=(interactions,burial)
        metrics=descriptions.describe(row)['metrics']
        self.assertEqual(metrics['total buried surface area'],50)
        self.assertEqual(metrics['total C/S contribution percent'],30)
        self.assertEqual(metrics['nonpolar contacts'],1)
        self.assertEqual(metrics['salt bridge candidates'],1)
        self.assertEqual(metrics['hydrogen-bond candidates'], 2)
        row.update(metrics)
        self.assertEqual([side['hydrogen-bond candidates'] for side in patch_export_rows([row])], [2, 2])
        report = io.StringIO()
        descriptions.write_patch(row, report)
        self.assertIn('Hydrogen-bond candidates (unique atom pairs within patch, ≤ 3.5 Å): 2',
                      report.getvalue())
        self.assertEqual(metrics['steric overlap candidates'],1)
        self.assertEqual(metrics['phosphate/basic contacts'],1)
        self.assertEqual(metrics['salt bridge candidate heavy-atom distance min'],3)
        self.assertIn('OD2',metrics['salt bridge residues']);self.assertNotIn('D9',metrics['salt bridge residues'])
        self.assertEqual(metrics['nominal patch charge difference'],0)
        self.assertEqual(metrics['patch 1 hydropathy mean'],1.8)
        self.assertNotIn('candidate DD interface',metrics)
        descriptions.details.clear()
        descriptions.pairs[('fold',0,'A','B')]=([],{(c,1):dict(buried_area=0,buried_nonpolar_area=0) for c in 'AB'})
        metrics=descriptions.describe(row)['metrics']
        self.assertIsNone(metrics['total C/S contribution percent'])
        self.assertIsNone(metrics['salt bridge candidate heavy-atom distance min'])
        self.assertIsNone(metrics['salt bridge residues'])

    def test_salt_proximity_scope_boundary_and_dedup(self):
        import copy
        import numpy as np
        from sieve.patch_description_summary import salt_patch_residues
        from sieve.interface_description import DescriptionOptions
        descriptions=PatchDescriptions({('fold',0):CIF},None)
        model=descriptions.models[('fold',0)][0]
        for chain in model:
            original=chain[1]
            for atom in original:atom.coord=np.array([0.,0.,0.])
            for n,z in ((2,3.5),(3,7.0),(4,3.5001)):
                residue=copy.deepcopy(original);residue.id=(' ',n,' ')
                for atom in residue:atom.coord=np.array([0.,0.,z])
                chain.add(residue)
            self.assertEqual(salt_patch_residues(chain,{1},3.5),{1,2})
            self.assertEqual(salt_patch_residues(chain,set(),3.5),set())
            self.assertEqual(salt_patch_residues(chain,{1},3.4),{1})
            # Nearby hydrogen/deuterium must not bring residue 3 into the shell.
            from Bio.PDB.Atom import Atom
            for element in ('H','D'):
                chain[3].add(Atom(element,np.array([0.,0.,0.]),0,1,' ',element,99,element=element))
            self.assertEqual(salt_patch_residues(chain,{1},3.5),{1,2})
        residues,pae=parse_model(CIF,DATA)
        row=detect_patches(residues,pae,Options(min_residues_per_chain=1,min_contacts=1))[-1]
        row.update(source='fold',model=0)
        base=dict(type='salt_bridge_candidate',residue_1=1,residue_2=2,
                  aa_1='K',aa_2='D',atom_1='NZ',atom_2='OD1',distance=3.)
        interactions=[base,dict(base,atom_2='OD2',distance=2.7),dict(base,residue_1=2),
                      dict(base,residue_2=3,distance=2.),dict(base,type='nonpolar_contact')]
        burial={(c,1):dict(buried_area=0,buried_nonpolar_area=0) for c in 'AB'}
        descriptions.pairs[('fold',0,'A','B')]=(interactions,burial)
        metrics=descriptions.describe(row)['metrics']
        self.assertEqual(metrics['salt bridge candidates'],2)
        self.assertEqual(metrics['salt bridge candidate heavy-atom distance min'],2.7)
        self.assertEqual(metrics['nonpolar contacts'],0)
        self.assertNotIn('D3',metrics['salt bridge residues'])
        output=io.StringIO();descriptions.write_patch(row,output)
        self.assertIn('salt_bridge_candidate 2',output.getvalue())
        self.assertEqual(output.getvalue().count('Salt candidate ('),2)
        self.assertIn('within 3.5 Å of same-chain patch; counted',output.getvalue())
        for invalid in (0,-1,float('nan'),float('inf')):
            with self.assertRaises(ValueError):DescriptionOptions(salt_patch_radius=invalid)

    def test_export_only_leaves_and_removes_internal_fields(self):
        residues,pae=parse_model(CIF,DATA)
        rows=detect_patches(residues,pae,Options(min_residues_per_chain=1,min_contacts=1))
        rows.append(dict(rows[-1],status='split'))
        parsed = list(patch_export_rows(rows))
        self.assertEqual(len(parsed),2)
        row=parsed[0]
        for field in ('parent id','depth','graph residues','algorithm version','separator residues','discarded residues','interface residues 1','residue count 1'):
            self.assertNotIn(field,row)
        self.assertEqual(row['patch residues'], 1)
        self.assertIn('interchain CA distance min',row)
        self.assertIn('nearest interchain CA distance min',row)

    def test_tall_header_without_leaves(self):
        self.assertEqual(list(patch_export_rows([])), [])
        self.assertEqual(tuple(PATCH_COLUMNS), (
            'fold_id', 'sequence_id', 'source', 'model', 'patch_id', 'patch label',
            'sequence_index', 'chain', 'region', 'patch residues', 'interface pLDDT',
            'normalized PAE to partner', 'pDockQ2 to partner',
            'nearest interchain CA distance min', 'nearest interchain CA distance max',
            'nearest interchain CA distance median', 'nearest interchain CA distance mean',
            'interface structure', 'patch hydropathy mean', 'scope', 'status', 'rank',
            'contact count', 'pDockQ2 max', 'distance cutoff', 'max bridge residues',
            'min residues per chain', 'min contacts', 'max separator checks',
            'interchain CA distance min', 'interchain CA distance max',
            'interchain CA distance median', 'interchain CA distance mean',
            'total buried surface area', 'total C/S contribution percent', 'nonpolar contacts',
            'salt bridge candidates', 'steric overlap candidates', 'phosphate/basic contacts',
            'salt bridge candidate heavy-atom distance min', 'salt bridge residues',
            'nominal patch charge difference', 'hydrogen-bond candidates', 'confidence origin',
        ))

    def test_tall_pairs_preserve_direction_rank_and_missing_fold_ids(self):
        leaves = [
            {'source': source, 'patch_id': patch_id, 'rank': rank, 'scope': 'patch',
             'status': 'leaf', 'chain 1': 'B', 'chain 2': 'A',
             'normalized PAE 1 to 2': 0.25, 'normalized PAE 2 to 1': 0.75,
             'pDockQ2 1 to 2': 0.2, 'pDockQ2 2 to 1': 0.8, 'contact count': 4}
            for source, patch_id, rank in [('first', 'p2', 5), ('first', 'p1', 3), ('second', 'p1', 8)]
        ]
        rows = list(patch_export_rows(leaves))
        self.assertEqual(len(rows), 6)
        for i, (leaf, rank) in enumerate(zip(leaves, (2, 1, 1))):
            pair = rows[2*i:2*i+2]
            self.assertEqual([row['sequence_index'] for row in pair], [1, 2])
            self.assertEqual([row['chain'] for row in pair], ['B', 'A'])
            self.assertEqual([row['normalized PAE to partner'] for row in pair], ['0.250000', '0.750000'])
            self.assertEqual([row['pDockQ2 to partner'] for row in pair], ['0.200000', '0.800000'])
            for row in pair:
                self.assertEqual((row['source'], row['patch_id'], row['rank'], row['fold_id'], row['contact count']),
                                 (leaf['source'], leaf['patch_id'], rank, 'NA', 4))

    def test_empty_export_after_narrower_cutoff(self):
        residues, pae = parse_model(CIF, DATA)
        for cutoff, expected in ((8.0, 2), (0.1, 0)):
            rows = detect_patches(residues, pae, Options(cutoff=cutoff, min_residues_per_chain=1, min_contacts=1))
            self.assertEqual(len(list(patch_export_rows(rows))), expected)
