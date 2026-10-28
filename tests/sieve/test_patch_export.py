import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from sieve.alphafold_manifest import sequence_id, fold_id
from sieve.alphafold_pdockq2_patches import write_patches, detect_patches, Options
from sieve.alphafold_pdockq2 import parse_model
from sieve.patch_description_summary import PatchDescriptions
from tests.scripts.test_alphafold_pdockq2_patches import SCRIPT
from tests.sieve.test_alphafold_pdockq2_patches import CIF, DATA
from tests.sieve.test_alphafold_overlap import fixture


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
        burial={('A',1):dict(buried_area=20,buried_nonpolar_area=10),
                ('B',1):dict(buried_area=30,buried_nonpolar_area=5)}
        descriptions.pairs[('fold',0,'A','B')]=(interactions,burial)
        metrics=descriptions.describe(row)['metrics']
        self.assertEqual(metrics['total buried surface area'],50)
        self.assertEqual(metrics['total C/S contribution percent'],30)
        self.assertEqual(metrics['nonpolar contacts'],1)
        self.assertEqual(metrics['salt bridge candidates'],1)
        self.assertEqual(metrics['steric overlap candidates'],1)
        self.assertEqual(metrics['phosphate/basic contacts'],1)
        self.assertEqual(metrics['salt bridge candidate heavy-atom distance min'],3)
        self.assertIn('OD2',metrics['salt bridge residues']);self.assertNotIn('D9',metrics['salt bridge residues'])
        self.assertEqual(metrics['nominal patch charge difference'],0)
        self.assertEqual(metrics['patch 1 hydropathy mean'],1.8)
        self.assertIsNone(metrics['candidate DD family face 1'])
        self.assertIsNone(metrics['candidate DD interface'])
        descriptions.dd=True
        descriptions.details.clear()
        with patch('sieve.patch_description_summary.candidate_faces', side_effect=['Ia','Ib (partial)']):
            enriched=dict(row,**descriptions.describe(row)['metrics'])
        out=io.StringIO();write_patches([enriched],out)
        exported=list(csv.DictReader(io.StringIO(out.getvalue()),delimiter='\t'))
        self.assertEqual(exported[0]['candidate DD interface'],'Variant')
        descriptions.dd=False
        descriptions.details.clear()
        descriptions.pairs[('fold',0,'A','B')]=([],{(c,1):dict(buried_area=0,buried_nonpolar_area=0) for c in 'AB'})
        metrics=descriptions.describe(row)['metrics']
        self.assertIsNone(metrics['total C/S contribution percent'])
        self.assertIsNone(metrics['salt bridge candidate heavy-atom distance min'])
        self.assertIsNone(metrics['salt bridge residues'])

    def test_export_only_leaves_and_removes_internal_fields(self):
        residues,pae=parse_model(CIF,DATA)
        rows=detect_patches(residues,pae,Options(min_residues_per_chain=1,min_contacts=1))
        rows.append(dict(rows[-1],status='split'))
        out=io.StringIO();write_patches(rows,out)
        parsed=list(csv.DictReader(io.StringIO(out.getvalue()),delimiter='\t'))
        self.assertEqual(len(parsed),1)
        row=parsed[0]
        for field in ('parent id','depth','graph residues','algorithm version','separator residues','discarded residues','interface residues 1','residue count 1'):
            self.assertNotIn(field,row)
        self.assertEqual(row['patch 1 residues'],'1')
        self.assertIn('interchain CA distance min',row)
        self.assertIn('nearest interchain CA distance min 2',row)

    def test_cli_silent_identity_recompute_and_empty(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);fold=root/'fold';fold.mkdir();output=root/'patches.tsv'
            (fold/'fold_model_0.cif').write_text(fixture([('A','AAAA',[0,0,0]),('B','CCCC',[6,0,0])]))
            (fold/'fold_job_request.json').write_text('{"modelSeeds":[42]}')
            data=dict(token_chain_ids=['A']*4+['B']*4,token_res_ids=list(range(1,5))*2,pae=[[1]*8 for _ in range(8)])
            (fold/'fold_full_data_0.json').write_text(json.dumps(data))
            args=[str(fold),'--silent','--no-dssp','--output',str(output),'--min-residues-per-chain','1','--min-contacts','1','--dd']
            with contextlib.redirect_stdout(io.StringIO()) as stdout:
                self.assertEqual(SCRIPT.main(args),0)
            self.assertEqual(stdout.getvalue(),'')
            with output.open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
            self.assertTrue(rows)
            self.assertEqual(rows[0]['sequence_id'],sequence_id(['AAAA','CCCC']))
            self.assertEqual(rows[0]['fold_id'],fold_id(42,['AAAA','CCCC']))
            self.assertEqual(rows[0]['candidate DD family face 1'],'NA')
            self.assertEqual(rows[0]['candidate DD interface'],'NA')
            self.assertEqual(rows[0]['dd enabled'],'True')
            self.assertEqual(rows[0]['dd helix presence threshold'],'3')
            self.assertEqual(rows[0]['dd loop presence threshold'],'1')
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(SCRIPT.main(args+['--dd-helix-presence-threshold','4',
                                                  '--dd-loop-presence-threshold','2']),0)
            with output.open() as stream:alias_rows=list(csv.DictReader(stream,delimiter='\t'))
            self.assertEqual(alias_rows[0]['dd helix presence threshold'],'4')
            self.assertEqual(alias_rows[0]['dd loop presence threshold'],'2')
            with contextlib.redirect_stdout(io.StringIO()) as stdout:
                self.assertEqual(SCRIPT.main(args+['--distance-cutoff','.1']),0)
            with output.open() as stream:self.assertEqual(list(csv.DictReader(stream,delimiter='\t')),[])
            self.assertEqual(stdout.getvalue(),'')
            with contextlib.redirect_stdout(io.StringIO()) as stdout:
                self.assertEqual(SCRIPT.main([str(fold),'--silent','--no-dssp']),0)
            self.assertEqual(stdout.getvalue(),'')
            for value in ('0','1.1','nan'):
                with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):
                    SCRIPT.main(args+['--dd-helix-presence-threshold',value])
                with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):
                    SCRIPT.main(args+['--dd-loop-presence-threshold',value])
