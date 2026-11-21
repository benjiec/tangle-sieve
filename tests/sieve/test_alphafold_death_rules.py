from pathlib import Path
import tempfile
import unittest

from sieve.alphafold_death_rules import read_domains, read_rules, classify, candidate_interface
from sieve.alphafold_patch_identity import patch_identifier
from tests.sieve.test_alphafold_dd_faces import structure

ROOT=Path(__file__).resolve().parents[2]
RULES=ROOT/'data/death-face-rules.tsv'
SEQUENCE_HEADER = 'Sequence Name\tProtein Accession\tDomain Type\tDomain Accession\n'


class RuleTests(unittest.TestCase):
    def test_domains_extra_columns_exact_names_duplicate_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'domains.tsv'
            path.write_text(
                'sequence name\tProtein Accession\tdomain type\tDomain Accession\n'
                ' name \tp.1\tDD\tp.1_PF00531_1_1\nB\tq.1\tDED\n'
            )
            self.assertEqual(read_domains(path), {' name ': 'DD', 'B': 'DED'})
            for text, message in (
                (SEQUENCE_HEADER + 'A\tp.1\tDD\nA\tp.2\tDD\n', 'duplicate sequence assignment'),
                (SEQUENCE_HEADER + 'A\tp.1\t\n', 'incomplete domain assignment'),
                (SEQUENCE_HEADER + 'A\tp.1\n', 'incomplete domain assignment'),
                (SEQUENCE_HEADER + ' \tp.1\tDD\n', 'incomplete domain assignment'),
                (SEQUENCE_HEADER + 'A\tp.1\tDD\tp.1_PF00531_1_1\textra\n', 'incomplete domain assignment'),
                ('Name\tDomain Type\nA\tDD\n', 'requires unique Sequence Name and Domain Type'),
                ('Sequence Name\tDomain\nA\tDD\n', 'requires unique Sequence Name and Domain Type'),
                ('Sequence Name\tDomain Type\tdomain type\nA\tDD\tDED\n',
                 'requires unique Sequence Name and Domain Type'),
            ):
                with self.subTest(text=text):
                    path.write_text(text)
                    with self.assertRaisesRegex(ValueError, message):
                        read_domains(path)

    def test_mandatory_or_partial_optional_and_mixed_pairing(self):
        rules,_=read_rules(RULES);c,s=structure()
        # DED Ib: H2 + either H3 or L23; optional H5 alone does nothing.
        self.assertIn('Ib', classify('A',{7,8,9,13},c,s,rules['DED'])[0].split('; '))
        self.assertIn('Ib', classify('A',{7,8,9,11,12},c,s,rules['DED'])[0].split('; '))
        self.assertIn('Ib (partial)',classify('A',{7},c,s,rules['DED'])[0])
        self.assertEqual(classify('A',{25},c,s,{'Ib':rules['DED']['Ib']})[0],'none')
        self.assertEqual(classify('A',{13},c,s,{'IIIa':rules['DED']['IIIa']})[0],'IIIa')
        self.assertEqual(classify('A',{13},c,s,{'IIIa':rules['DD']['IIIa']})[0],'IIIa (partial)')
        self.assertEqual(candidate_interface('Ia','Ib'),'Canonical')
        c[('A',5)]='G'
        self.assertEqual(classify('A',{5,17},c,s,{'IIIb':rules['DD']['IIIb']})[0],'IIIb')
        c.pop(('A',6))
        faces,_,status=classify('A',{5,17},c,s,rules['DD'])
        self.assertIsNone(faces);self.assertIn('missing',status)
        self.assertIsNone(classify('B',set(),c,s,rules['DD'])[0])

    def test_rules_validation_and_hash_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'rules.tsv';lines=RULES.read_text().splitlines()
            original=read_rules(RULES)[1]
            path.write_text('\n'.join([lines[0],*reversed(lines[1:])])+'\n')
            self.assertEqual(read_rules(path)[1],original)
            for row in ('DD\tIa\talpha7\t3\tmandatory\tx',
                        'DD\tIa\talpha1\t0\tmandatory\tx',
                        'DD\tIa\talpha1\t1.5\tmandatory\tx',
                        'DD\tIa\talpha1\t3\tunknown\tx',
                        'DD\tIa\talpha1\t0\toptional\tx'):
                path.write_text(lines[0]+'\n'+row+'\n')
                with self.assertRaises(ValueError):read_rules(path)
            path.write_text('\n'.join([*lines,lines[1]])+'\n')
            with self.assertRaises(ValueError):read_rules(path)

    def test_patch_hash_membership_and_identity(self):
        sides=[('A',[3,2,1,1]),('B',[5,6])]
        value=patch_identifier(['fold',42],0,sides)
        self.assertEqual(len(value),32)
        self.assertEqual(value,patch_identifier(['fold',42],0,[('B',[6,5]),('A',[1,2,3])]))
        for identity,model,members in ((['fold',43],0,sides),(['fold',42],1,sides),
                                      (['fold',42],0,[('C',[1,2,3]),('B',[5,6])]),
                                      (['fold',42],0,[('A',[1,2]),('B',[5,6])])):
            self.assertNotEqual(value,patch_identifier(identity,model,members))
