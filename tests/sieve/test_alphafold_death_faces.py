import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from sieve.alphafold_death_rules import read_domains, read_rules, classify, candidate_interface
from sieve.alphafold_patch_identity import patch_identifier
from sieve.alphafold_death_faces import run
from sieve.alphafold_manifest import run as manifest
from sieve.alphafold_result_format import PatchInput
from sieve.alphafold_pdockq2_patches import analyze, Options
from tests.sieve.test_alphafold_dd_faces import structure
from tests.sieve.test_alphafold_overlap import fixture
from tests.scripts.helpers import load_script

ROOT=Path(__file__).resolve().parents[2]
RULES=ROOT/'data/death-face-rules.tsv'


class RuleTests(unittest.TestCase):
    def test_domains_extra_columns_exact_names_duplicate_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'domains.tsv'
            path.write_text('sequence name\tdomain\tOther\n name \tDD\tignored\n')
            self.assertEqual(read_domains(path),{' name ':'DD'})
            for text in ('Sequence Name\tDomain\nA\tDD\nA\tDD\n',
                         'Sequence Name\tDomain\nA\t\n', 'Name\tDomain\nA\tDD\n'):
                path.write_text(text)
                with self.assertRaises(ValueError):read_domains(path)

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


class RunnerTests(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory();self.addCleanup(temp.cleanup)
        self.root=Path(temp.name);self.fold=self.root/'fold';self.fold.mkdir()
        text=fixture([('A','AAAA',[0,0,0]),('B','CCCC',[6,0,0])])
        for n in (0,1):
            (self.fold/f'fold_model_{n}.cif').write_text(text)
            (self.fold/f'fold_full_data_{n}.json').write_text(json.dumps(dict(
                token_chain_ids=['A']*4+['B']*4,token_res_ids=[1,2,3,4]*2,pae=[[1]*8 for _ in range(8)])))
        (self.fold/'fold_job_request.json').write_text('{"modelSeeds":[42]}')
        fasta=self.root/'seq.faa';fasta.write_text('>first\nAAAA\n>second\nCCCC\n')
        self.manifest=self.root/'manifest.tsv'
        self.assertEqual(manifest(self.fold,fasta,self.manifest,io.StringIO()),0)
        self.domains=self.root/'domains.tsv';self.domains.write_text('Sequence Name\tDomain\nfirst\tDD\nsecond\tDED\n')
        self.output=self.root/'faces.tsv';self.err=io.StringIO()
        self.options=Options(min_residues_per_chain=1,min_contacts=1)

    def execute(self, requests=None):
        return run(requests or [PatchInput(str(self.fold))],self.manifest,self.domains,RULES,self.output,
                   options=self.options,stderr=self.err)

    def test_shared_hash_model_selection_and_removed_patch_columns(self):
        with patch('sieve.alphafold_death_faces.find_dssp',return_value='fake'), patch(
                'sieve.alphafold_death_faces.load_secondary',return_value=(None,{},[],{'status':'assigned'})):
            self.assertEqual(self.execute([PatchInput(str(self.fold),1)]),0)
        with self.output.open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
        expected=analyze(str(self.fold),options=self.options,selected_model=1)
        self.assertEqual([r['patch_id'] for r in rows],[r['patch_id'] for r in expected if r['scope']=='patch' and r['status']=='leaf'])
        wider=analyze(str(self.fold),options=Options(cutoff=10,min_residues_per_chain=1,min_contacts=1),selected_model=1)
        self.assertEqual(rows[0]['patch_id'],next(r['patch_id'] for r in wider if r['scope']=='patch' and r['status']=='leaf'))
        self.assertNotIn('face evidence 1',rows[0])
        self.assertNotIn('face evidence 2',rows[0])
        self.assertEqual(rows[0]['candidate interface'],'NA')
        self.assertEqual(rows[0]['model_id'],'1')
        self.assertEqual(rows[0]['domain 2'],'DED')
        patch_script=load_script(ROOT/'scripts/alphafold-pdockq2-patches.py')
        out=self.root/'patches.tsv'
        with patch('sys.stdout',io.StringIO()):
            self.assertEqual(patch_script.main([str(self.fold),'--no-dssp','--silent','--output',str(out),
                '--min-residues-per-chain','1','--min-contacts','1']),0)
        with out.open() as stream:patch_rows=list(csv.DictReader(stream,delimiter='\t'))
        self.assertNotIn('candidate DD interface',patch_rows[0])
        self.assertEqual(rows[0]['patch_id'],next(r['patch_id'] for r in patch_rows if r['model']=='1'))

    def test_mapping_errors_preserve_output_before_dssp(self):
        self.output.write_text('existing')
        for domains in ('Sequence Name\tDomain\nfirst\tDD\n',
                        'Sequence Name\tDomain\nfirst\tXX\nsecond\tDED\n'):
            self.domains.write_text(domains)
            with patch('sieve.alphafold_death_faces.find_dssp',return_value='fake'), patch('sieve.alphafold_death_faces.load_secondary') as dssp:
                with self.assertRaises(ValueError):self.execute()
                dssp.assert_not_called()
            self.assertEqual(self.output.read_text(),'existing')

    def test_model_load_and_dssp_errors_continue(self):
        (self.fold/'fold_full_data_0.json').write_text('broken')
        with patch('sieve.alphafold_death_faces.find_dssp',return_value='fake'), patch(
                'sieve.alphafold_death_faces.load_secondary',return_value=(None,{},[],{'status':'assigned'})):
            self.assertEqual(self.execute(),1)
        with self.output.open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
        self.assertEqual({r['model_id'] for r in rows},{'1'})
        self.assertIn('ERROR:',self.err.getvalue())
        with patch('sieve.alphafold_death_faces.find_dssp',return_value='fake'), patch(
                'sieve.alphafold_death_faces.load_secondary',side_effect=RuntimeError('DSSP failed')):
            self.assertEqual(self.execute(),1)
        with self.output.open() as stream:self.assertEqual(list(csv.DictReader(stream,delimiter='\t')),[])

    def test_pipeline_and_atomic_failure(self):
        script=load_script(ROOT/'scripts/alphafold-death-faces.py')
        with self.manifest.open() as stream:identity=next(csv.DictReader(stream,delimiter='\t'))
        stdin=io.StringIO(f"{self.fold}\t{identity['fold_id']}\t1\t0.5\t0.6\n")
        with patch('sys.stdin',stdin), patch('sys.stdout',io.StringIO()) as stdout, patch(
                'sieve.alphafold_death_faces.find_dssp',return_value='fake'), patch(
                'sieve.alphafold_death_faces.load_secondary',return_value=(None,{},[],{'status':'assigned'})):
            self.assertEqual(script.main(['-','--manifest',str(self.manifest),'--domains',str(self.domains),
                '--rules',str(RULES),'--output',str(self.output),'--min-residues-per-chain','1','--min-contacts','1']),0)
            self.assertEqual(stdout.getvalue(),'')
            before=self.output.read_bytes()
            with patch('sieve.alphafold_cache.os.replace',side_effect=OSError('disk failure')):
                with self.assertRaises(OSError):self.execute()
            self.assertEqual(self.output.read_bytes(),before)

    def test_manifest_identity_mismatch_is_fatal(self):
        text=self.manifest.read_text()
        with self.manifest.open() as stream:row=next(csv.DictReader(stream,delimiter='\t'))
        self.manifest.write_text(text.replace(row['fold_id'],'00000000'))
        self.output.write_text('existing')
        with patch('sieve.alphafold_death_faces.find_dssp',return_value='fake'), patch('sieve.alphafold_death_faces.load_secondary') as dssp:
            with self.assertRaisesRegex(ValueError,'identity mismatch'):self.execute()
            dssp.assert_not_called()
        self.assertEqual(self.output.read_text(),'existing')
