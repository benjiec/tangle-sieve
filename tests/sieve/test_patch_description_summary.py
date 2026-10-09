import io
import unittest
from unittest.mock import patch
from sieve.patch_description_summary import PatchDescriptions, residue_ranges
from sieve.alphafold_pdockq2 import parse_model
from sieve.alphafold_pdockq2_patches import detect_patches, Options
from tests.sieve.test_alphafold_pdockq2_patches import CIF, DATA


class DescriptionSummaryTests(unittest.TestCase):
    def test_identical_sequence_different_structures_not_collapsed(self):
        segments=[{'chain':'A','label':'alpha1','code':'H','start':1,'end':1,'residues':[1]},
                  {'chain':'B','label':'beta1','code':'E','start':1,'end':1,'residues':[1]}]
        with patch('sieve.patch_description_summary.secondary_structure', return_value=({('A',1):'H',('B',1):'E'},segments,{})):
            descriptions=PatchDescriptions({('one',0):CIF,('two',0):CIF},'fake')
        output=io.StringIO(); descriptions.write_chains('one',output)
        text=output.getvalue()
        self.assertEqual(text.count('Sequence group'),1)
        self.assertIn('chains A, B',text)
        self.assertIn('alpha1[H](1–1)',text)
        self.assertIn('beta1[E](1–1)',text)
        self.assertEqual(text.count('model 0 chain A'),1)

    def test_patch_uses_actual_selected_residues_and_local_confidence(self):
        descriptions=PatchDescriptions({('one',0):CIF},None)
        residues,pae=parse_model(CIF,DATA)
        row=detect_patches(residues,pae,Options(min_residues_per_chain=1,min_contacts=1))[1]
        row.update(source='one',model=0)
        output=io.StringIO(); descriptions.write_patch(row,output)
        text=output.getvalue()
        self.assertIn('A 90.00 | B 80.00',text)
        self.assertIn('unassigned:1',text)
        self.assertIn('spatial_shell: n=0',text)
        self.assertIn('salt_bridge_candidate 0',text)
        self.assertNotIn('alignment',text.lower())
        self.assertEqual(residue_ranges({1,2,4}),'1–2,4')

    def test_salt_scope_and_atom_deduplication(self):
        descriptions=PatchDescriptions({('one',0):CIF},None)
        residues,pae=parse_model(CIF,DATA)
        row=detect_patches(residues,pae)[0]; row.update(source='one',model=0)
        base={'type':'salt_bridge_candidate','residue_1':1,'residue_2':1,
              'aa_1':'K','aa_2':'D','atom_1':'NZ','atom_2':'OD1','distance':3.5}
        interactions=[base,dict(base,atom_2='OD2',distance=3.0),dict(base,residue_2=9),dict(base,residue_1=8,residue_2=9)]
        burial={(c,1):{'buried_area':10,'buried_nonpolar_area':5} for c in ('A','B')}
        descriptions.pairs[('one',0,'A','B')]=(interactions,burial)
        output=io.StringIO(); descriptions.write_patch(row,output)
        text=output.getvalue()
        self.assertIn('salt_bridge_candidate 1',text)
        self.assertEqual(text.count('Salt candidate ('),2)
        self.assertIn('within patch',text)
        self.assertIn('one endpoint in patch',text)
        self.assertIn('/OD2 3.00 Å',text)

if __name__=='__main__': unittest.main()

class ExecutableDiscoveryTests(unittest.TestCase):
    def test_precedence_and_missing(self):
        from sieve.patch_description_summary import find_dssp
        with patch.dict('os.environ', {'SIEVE_DSSP':'/env/mkdssp'}, clear=True), patch('sieve.patch_description_summary.shutil.which', side_effect=lambda value: value) as which:
            self.assertEqual(find_dssp('/explicit/mkdssp'),'/explicit/mkdssp')
            which.assert_called_once_with('/explicit/mkdssp')
            self.assertEqual(find_dssp(),'/env/mkdssp')
        with patch.dict('os.environ', {}, clear=True), patch('sieve.patch_description_summary.shutil.which', return_value=None) as which:
            self.assertIsNone(find_dssp())
            which.assert_called_once_with('mkdssp')

    def test_invalid_and_empty_configuration_do_not_fall_back(self):
        from sieve.patch_description_summary import find_dssp
        for value in ('', '/missing/mkdssp'):
            with patch.dict('os.environ', {'SIEVE_DSSP':value}, clear=True), patch('sieve.patch_description_summary.shutil.which', return_value=None):
                with self.assertRaisesRegex(ValueError,'SIEVE_DSSP'):
                    find_dssp()
                with self.assertRaisesRegex(ValueError,'--dssp'):
                    find_dssp(value)

    def test_executable_path_with_spaces(self):
        import tempfile
        from pathlib import Path
        from sieve.patch_description_summary import find_dssp
        with tempfile.TemporaryDirectory() as directory:
            executable=Path(directory)/'my dssp'
            executable.write_text('#!/bin/sh\nexit 0\n')
            executable.chmod(0o755)
            self.assertEqual(find_dssp(str(executable)),str(executable))
            executable.chmod(0o644)
            with self.assertRaises(ValueError):
                find_dssp(str(executable))
