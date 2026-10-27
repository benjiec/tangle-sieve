import io
import unittest
from unittest.mock import patch

import numpy as np
from Bio.PDB import Atom, Chain, Model, Residue

from sieve.interface_description import composition, atom_interactions, DescriptionOptions, neighborhood
from sieve.patch_description_summary import PatchDescriptions, charge_text
from tests.sieve.test_alphafold_ptm import cif


def model(ptm='SEP', basic='LYS'):
    result = Model.Model(0)
    for c, name, atoms in [
        ('A', ptm, [('CA', (0, 5, 0), 'C'), ('P', (-1.5, 0, 0), 'P'),
                    ('O1P', (0, 0, 0), 'O'), ('O2P', (-1.5, 1.5, 0), 'O'),
                    ('O3P', (-1.5, -1.5, 0), 'O')]),
        ('B', basic, [('CA', (4, 5, 0), 'C'),
                      ('NZ' if basic == 'LYS' else 'NH1', (4, 0, 0), 'N')])]:
        chain = Chain.Chain(c); result.add(chain)
        residue = Residue.Residue((' ', 1, ' '), name, ' '); chain.add(residue)
        for n, (a, coord, element) in enumerate(atoms):
            residue.add(Atom.Atom(a, np.array(coord, float), 90, 1, ' ', a, n, element=element))
    return result


class PhosphorylationTests(unittest.TestCase):
    def test_charge_scenarios_and_no_parent_hydropathy(self):
        for ptm in ('SEP', 'TPO', 'PTR'):
            with self.subTest(ptm=ptm):
                m = model(ptm)
                stats = composition([m['A'][1], m['B'][1]])
                self.assertEqual(stats['sequence'], f'[{ptm}]K')
                self.assertEqual(stats['count'], 2)
                self.assertIsNone(stats['mean_hydropathy'])
                self.assertIsNone(stats['hydrophobic_fraction'])
                self.assertIsNone(stats['nominal_sidechain_charge'])
                self.assertEqual(stats['charge_phosphate_minus1'], 0)
                self.assertEqual(stats['charge_phosphate_minus2'], -1)
                self.assertEqual(stats['positive_count'], 1)
                self.assertIn('phosphate −2', charge_text(stats))
                self.assertEqual(composition(m['B'])['nominal_sidechain_charge'], 1)
        unknown = composition(model('ZZZ')['A'])
        self.assertIsNone(unknown['charge_phosphate_minus1'])
        self.assertIn('unsupported PTM', charge_text(unknown))

    def test_phosphate_basic_boundary_and_reverse(self):
        for ptm in ('SEP', 'TPO', 'PTR'):
            for basic in ('LYS', 'ARG'):
                with self.subTest(ptm=ptm, basic=basic):
                    m = model(ptm, basic)
                    rows = atom_interactions(m, ('A', 'B'), DescriptionOptions())
                    contacts = [r for r in rows if r['type'] == 'phosphate_basic_contact']
                    self.assertEqual(len(contacts), 1)
                    self.assertEqual(contacts[0]['atom_1'], 'O1P')
                    self.assertEqual(contacts[0]['distance'], 4)
                    reverse = atom_interactions(m, ('B', 'A'), DescriptionOptions())
                    self.assertEqual(len([r for r in reverse if r['type'] == 'phosphate_basic_contact']), 1)
                    m['B'][1]['NZ' if basic == 'LYS' else 'NH1'].coord[0] = 4.001
                    self.assertFalse(any(r['type'] == 'phosphate_basic_contact'
                                         for r in atom_interactions(m, ('A', 'B'), DescriptionOptions())))
        m = model('ZZZ')
        self.assertFalse(any(r['type'] == 'phosphate_basic_contact'
                             for r in atom_interactions(m, ('A', 'B'), DescriptionOptions())))

    def test_phosphorus_clash_and_unsupported_element(self):
        m = model()
        m['B'][1]['CA'].coord = np.array((0, 0, 0), float)
        clashes = [r for r in atom_interactions(m, ('A', 'B'), DescriptionOptions())
                   if r['type'] == 'steric_overlap_candidate' and r['atom_1'] == 'P' and r['atom_2'] == 'CA']
        self.assertEqual(len(clashes), 1)
        self.assertAlmostEqual(clashes[0]['vdw_overlap'], 2.0)
        m['B'][1]['CA'].element = 'XX'
        with self.assertRaisesRegex(ValueError, 'unsupported element'):
            atom_interactions(m, ('A', 'B'), DescriptionOptions())

    def test_unmodified_patch_with_phosphorylated_shell(self):
        m = model(); chain = m['A']
        r = Residue.Residue((' ', 2, ' '), 'LYS', ' '); chain.add(r)
        r.add(Atom.Atom('CA', np.array((1, 5, 0), float), 90, 1, ' ', 'CA', 20, element='C'))
        values = neighborhood(chain, {2}, 1, 8)
        self.assertEqual(values['patch']['nominal_sidechain_charge'], 1)
        self.assertIsNotNone(values['patch']['mean_hydropathy'])
        self.assertIsNone(values['sequence_shell']['mean_hydropathy'])
        self.assertEqual(values['sequence_shell']['phosphate_count'], 1)

    def test_dssp_not_disabled_for_phosphorylation_and_coverage(self):
        text = cif({'A': ('ASST', ((2, 'SEP'), (3, 'SEP')))})
        with patch('sieve.patch_description_summary.secondary_structure', return_value=(
                {('A', 2): 'H'}, [], {})) as dssp:
            descriptions = PatchDescriptions({('fold', 0): text}, 'mkdssp-test')
            self.assertEqual(dssp.call_args.args[2], 'mkdssp-test')
        out = io.StringIO(); descriptions.write_chains('fold', out)
        self.assertIn('1/2 assigned', out.getvalue())
        self.assertIn('| 4 residues', out.getvalue())

    def test_summary_deduplicates_phosphate_contacts_and_keeps_normal_chemistry(self):
        from Bio.PDB import MMCIFIO
        m = model()
        m['A'][1]['O2P'].coord = np.array((0, 1, 0), float)
        m['B'][1]['NZ'].coord = np.array((3, 0, 0), float)
        writer = MMCIFIO(); writer.set_structure(m)
        text = io.StringIO(); writer.save(text)
        descriptions = PatchDescriptions({('fold', 0): text.getvalue()}, None)
        row = {'source': 'fold', 'model': 0, 'chain 1': 'A', 'chain 2': 'B',
               'interface residues 1': '1', 'interface residues 2': '1',
               'interface pLDDT 1': 90, 'interface pLDDT 2': 90}
        out = io.StringIO(); descriptions.write_patch(row, out)
        result = out.getvalue()
        self.assertIn('Phosphate/basic contacts (unique residue pairs within patch): 1', result)
        self.assertEqual(result.count('Phosphate/basic candidate ('), 1)
        self.assertIn('A:[SEP]1/O1P–B:K1/NZ 3.00 Å', result)
        self.assertIn('nominal charge +1 | K/R 1', result)
        self.assertIn('mean hydropathy n/a', result)
        self.assertIn('difference (A − B): -2 (phosphate −1); -3 (phosphate −2)', result)
        self.assertNotIn('WARNING: incomplete', result)
        m['A'][1].detach_child('O3P'); writer.set_structure(m)
        text = io.StringIO(); writer.save(text)
        descriptions = PatchDescriptions({('fold', 0): text.getvalue()}, None)
        out = io.StringIO(); descriptions.write_patch(row, out)
        self.assertIn('missing O3P', out.getvalue())
