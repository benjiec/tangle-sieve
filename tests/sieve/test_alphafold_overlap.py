import copy
import csv
import io
import json
import shutil
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
from Bio.PDB import MMCIFIO

from sieve import alphafold_overlap as app
from sieve.alphafold_structure import load_structure, open_structure_fold, assembly_metadata
from sieve.alphafold_manifest import chain_sequences, run as manifest, load_manifest
from sieve.alphafold_results import run as results
from sieve.alphafold_pdockq2_patches import analyze, write_summary, Options as PatchOptions
from sieve.alphafold_rosetta_integrity import prepare
from tests.scripts.helpers import load_script


POINTS=np.array([[0,0,0],[3,0,0],[3,3,0],[3,3,3]],dtype=float)
NAMES={'A':'ALA','C':'CYS','D':'ASP','E':'GLU','G':'GLY','K':'LYS'}


def fixture(chains, transform=None, missing=(), auth_offset=0):
    """Independent minimal mmCIF fixture with full sequence and observed coordinates."""
    data={'data_':'test'}
    def add(category,record):
        for k,v in record.items(): data.setdefault(category+'.'+k,[]).append(str(v))
    serial=0
    for entity,(c,seq,offset) in enumerate(chains,1):
        add('_entity',dict(id=entity,type='polymer'));add('_struct_asym',dict(id=c,entity_id=entity))
        for n,aa in enumerate(seq,1):
            add('_entity_poly_seq',dict(entity_id=entity,num=n,mon_id=NAMES[aa]))
            if (c,n) in missing: continue
            xyz=POINTS[n-1]+np.array(offset)
            if transform is not None: xyz=xyz@transform[0]+transform[1]
            for name,shift,element in [('CA',[0,0,0],'C'),('N',[-1,0,0],'N'),('O',[0,1,0],'O')]:
                serial+=1;p=xyz+shift
                add('_atom_site',dict(group_PDB='ATOM',id=serial,label_asym_id=c,auth_asym_id=c,label_entity_id=entity,
                    label_seq_id=n,auth_seq_id=n+auth_offset,label_comp_id=NAMES[aa],label_atom_id=name,type_symbol=element,
                    label_alt_id='.',pdbx_PDB_ins_code='?',Cartn_x=p[0],Cartn_y=p[1],Cartn_z=p[2],
                    B_iso_or_equiv=90,occupancy=1,pdbx_PDB_model_num=1))
    writer=MMCIFIO();writer.set_dict(data);out=io.StringIO();writer.save(out);return out.getvalue()


class OverlapTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.a=self.make('a',[('A','AAAA',[-7,0,0]),('B','CCCC',[0,0,0])])
        self.b=self.make('b',[('X','CCCC',[0,0,0]),('Y','DDDD',[7,0,0])])

    def make(self,name,chains,**kwargs):
        path=self.root/name;path.mkdir()
        (path/'fold_model_0.cif').write_text(fixture(chains,**kwargs))
        (path/'fold_job_request.json').write_text('{"modelSeeds":[42]}')
        return path

    def assemble(self,name='out',**kwargs):
        out=self.root/name
        meta=app.run(str(self.a),str(self.b),out,stdout=io.StringIO(),**kwargs)
        return out,meta

    def test_fit_proper_rotation_and_degeneracy(self):
        rng=np.random.default_rng(5);p=rng.normal(size=(12,3));rot=np.array([[0,-1,0],[1,0,0],[0,0,1]])
        q=p@rot+[4,6,9];r,t=app.fit(p,q)
        np.testing.assert_allclose(q@r+t,p,atol=1e-12)
        r,t=app.fit(p,p*[-1,1,1]);self.assertAlmostEqual(np.linalg.det(r),1)
        self.assertGreater(np.linalg.norm((p*[-1,1,1])@r+t-p),.1)
        for p in [POINTS[:2],np.array([[0,0,0],[1,0,0],[2,0,0]])]:
            with self.assertRaises(ValueError):app.fit(p,p)

    def test_output_roundtrip_identity_and_rerun(self):
        out,meta=self.assemble()
        model=load_structure(out)
        self.assertEqual(model.sequences,{'A':'AAAA','B':'CCCC','Y':'DDDD'})
        self.assertEqual(meta['fold_id'],model.identity['fold_id'])
        self.assertEqual(meta['sequence_id'],model.identity['sequence_id'])
        with self.assertRaisesRegex(ValueError,'never overwritten'): self.assemble()
        out2,meta2=self.assemble('rerun',options=app.Options(clash_overlap=.5))
        self.assertEqual(meta['fold_id'],meta2['fold_id'])
        self.assertEqual((out/'model0_complex.cif').read_bytes(),(out2/'model0_complex.cif').read_bytes())
        shutil.move(out2,self.root/'moved')
        self.assertEqual(load_structure(self.root/'moved').identity['fold_id'],meta['fold_id'])
        for name in ['comparisons','cross_partner_clashes','shared_interface_changes','alignment_sensitivity','shared_chain_deviations']:
            self.assertTrue((out/(name+'.tsv')).read_text().splitlines()[0])
        with self.assertRaises(ValueError):app.run(self.a,self.b,self.a,stdout=io.StringIO())

    def test_unique_shared_errors(self):
        a,b=load_structure(self.a),load_structure(self.b)
        for seqs in [{'X':'EEEE','Y':'DDDD'}, {'X':'CCCC','Y':'AAAA'}, {'X':'CCCC','Y':'CCCC'}]:
            bad=copy.deepcopy(b);bad.sequences=seqs
            with self.assertRaises(ValueError):app.unique_shared(a,bad)
        a.sequences['C']='CCCC'
        with self.assertRaisesRegex(ValueError,'duplicated'):app.unique_shared(a,b)

    def test_multichain_and_recursive(self):
        self.a=self.make('multi',[('A','AAAA',[-7,0,0]),('B','CCCC',[0,0,0]),('Z','EEEE',[0,-7,0])])
        out,meta=self.assemble()
        self.assertEqual(len(load_structure(out).sequences),4)
        with (out/'cross_partner_metrics.tsv').open() as stream:
            self.assertEqual(len(list(csv.DictReader(stream,delimiter='\t'))),2)
        third=self.make('third',[('B','DDDD',[7,0,0]),('Z','GGGG',[14,0,0])])
        final=self.root/'recursive'
        app.run(out,third,final,stdout=io.StringIO())
        reloaded=load_structure(final)
        self.assertEqual(len(reloaded.sequences),5)
        self.assertEqual(reloaded.metadata['inputs'][0]['history']['fold_id'],meta['fold_id'])
        previous=load_structure(out)
        for c in previous.atoms:
            np.testing.assert_allclose(app.coordinates(reloaded.atoms[c]),app.coordinates(previous.atoms[c]),atol=1e-6)
        # Earlier paths are not required after assembly.
        shutil.rmtree(out)
        self.assertEqual(len(load_structure(final).sequences),5)

    def test_missing_coordinates_author_numbering_standalone(self):
        self.a=self.make('missing',[('A','AAAA',[-7,0,0]),('B','CCCC',[0,0,0])],missing=[('B',4)],auth_offset=100)
        standalone=self.a/'plain.cif';standalone.write_text((self.a/'fold_model_0.cif').read_text())
        # No matching request name => unknown confidence in a standalone non-AF CIF.
        self.assertFalse(load_structure(standalone).plddt_known['B'])
        out,meta=self.assemble()
        model=load_structure(out)
        self.assertEqual(len(model.sequences['B']),4)
        self.assertEqual(len([a for a in model.atoms['B'] if a['name']=='CA']),3)
        self.assertEqual(meta['chains']['B']['residue_origins']['1']['auth_seq_id'],'101')
        self.assertIn('coverage', ';'.join(meta['assessment']['warnings']))

    def test_corruption_and_missing_provenance(self):
        out,meta=self.assemble()
        cif=out/'model0_complex.cif';original=cif.read_text();cif.write_text(original+'\n')
        with self.assertRaisesRegex(ValueError,'hash mismatch'):load_structure(out)
        cif.write_text(original)
        data=out/'provenance.json';data.unlink()
        with self.assertRaisesRegex(ValueError,'provenance'):load_structure(out)
        data.write_text(json.dumps(dict(meta,format='unknown')))
        with self.assertRaisesRegex(ValueError,'unsupported'):load_structure(out)

    def test_warning_classification_and_sensitivity(self):
        a,b=load_structure(self.a),load_structure(self.b)
        for atom in b.atoms['Y']:atom['coord'][0]-=14
        result=app.compare(a,b,app.Options())
        self.assertGreater(result['summary']['cross_partner_clashes'],0)
        self.assertGreater(result['summary']['cross_partner_backbone_clashes'],0)
        self.assertIn('WARNING',app.quick_summary(result,app.Options()))
        for atom in b.atoms['X']:
            if atom['number']==4:atom['coord'][2]+=20
        result=app.compare(a,b,app.Options())
        self.assertTrue(any('RMSD' in w for w in result['warnings']))
        for option in [dict(clash_overlap=0),dict(rmsd_warning=float('nan')),dict(coverage_warning=2)]:
            with self.assertRaises(ValueError):app.Options(**option)

    def test_manifest_results_patch_rosetta_and_zip(self):
        out,meta=self.assemble()
        fasta=self.root/'candidates.faa';fasta.write_text('>a\nAAAA\n>c\nCCCC\n>d\nDDDD\n')
        manifest_path=self.root/'manifest.tsv'
        self.assertEqual(manifest(out,fasta,manifest_path,io.StringIO()),0)
        cached,_=load_manifest(manifest_path)
        self.assertEqual(next(iter(cached.values()))[0]['source_type'],'assembled')
        text=io.StringIO()
        self.assertEqual(results(manifest_path,meta['sequence_id'],text,io.StringIO(),include_seed=True),0)
        self.assertTrue(text.getvalue().endswith('\t0\tNA\tNA\tNA\n'))
        text=io.StringIO();results(manifest_path,meta['sequence_id'],text,highest='ipTM')
        self.assertEqual(text.getvalue(),'')
        archive=self.root/'derived.zip'
        with zipfile.ZipFile(archive,'w') as z:
            for path in out.iterdir():z.write(path,'complex/'+path.name)
        for locator in [str(out),f'{archive}::complex']:
            targets={};ids={}
            rows=analyze(locator,options=PatchOptions(min_residues_per_chain=1,min_contacts=1),summary_targets=targets,fold_ids=ids)
            self.assertTrue(rows)
            self.assertTrue(all(r['pDockQ2 max'] is None for r in rows))
            self.assertEqual(next(iter(ids.values())),meta['fold_id'])
            stream=io.StringIO();write_summary(rows,targets,stream)
            self.assertIn('Geometric ranking',stream.getvalue())
            with tempfile.TemporaryDirectory() as tmp:
                metadata=prepare(locator,0,Path(tmp)/'prepared')
                self.assertEqual(metadata['fold_id'],meta['fold_id'])
        archive_manifest=self.root/'zip.tsv'
        self.assertEqual(manifest(archive,fasta,archive_manifest,io.StringIO()),0)
        script=load_script(Path(__file__).resolve().parents[2]/'scripts/alphafold-pdockq2-patches.py')
        import contextlib
        with contextlib.redirect_stdout(io.StringIO()) as printed:
            self.assertEqual(script.main([str(out),'--no-dssp','--min-residues-per-chain','1','--min-contacts','1']),0)
        self.assertIn('Inherited source pLDDT',printed.getvalue())
        self.assertIn('pDockQ2 NA',printed.getvalue())

    def test_incomplete_manifest_rejected(self):
        out,meta=self.assemble()
        fasta=self.root/'candidates.faa';fasta.write_text('>a\nAAAA\n>c\nCCCC\n>d\nDDDD\n')
        path=self.root/'cache.tsv'
        self.assertEqual(manifest(self.a,fasta,path,io.StringIO()),0)
        with path.open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
        old_id=rows[0]['fold_id'];columns=[c for c in rows[0] if c!='source_type']
        with path.open('w') as stream:
            w=csv.DictWriter(stream,columns,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
        original=path.read_bytes()
        with self.assertRaisesRegex(ValueError,'remove it and regenerate'):
            manifest(out,fasta,path,io.StringIO())
        self.assertEqual(path.read_bytes(),original)

    def test_all_models_error_report_and_atomic_failure(self):
        (self.a/'fold_model_1.cif').write_text('broken')
        # A malformed nonselected model is recorded without losing the selected assembly.
        broken,_=self.assemble('broken-extra',all_models=True)
        self.assertIn('ERROR', (broken/'comparisons.tsv').read_text())
        (self.a/'fold_model_1.cif').write_text(fixture([('A','AAAA',[-7,0,0]),('B','EEEE',[0,0,0])]))
        out,meta=self.assemble(all_models=True)
        with (out/'comparisons.tsv').open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
        self.assertEqual(len(rows),2);self.assertIn('ERROR',rows[1]['status'])
        with patch('sieve.alphafold_overlap.os.replace',side_effect=OSError('disk error')):
            with self.assertRaises(OSError):self.assemble('failed')
        self.assertFalse((self.root/'failed').exists())
        self.assertFalse(list(self.root.glob('.overlap-*')))

    def test_clash_threshold_boundaries_and_unknown_elements(self):
        atom=dict(number=1,name='CA',resname='ALA',element='C',coord=[0.,0.,0.],bfactor=90,occupancy=1)
        for distance,count in [(3.01,0),(2.99,1),(20,0)]:
            summary,details,_=app.pair_geometry([atom],[dict(atom,coord=[distance,0,0])],app.Options())
            self.assertEqual(summary['clash_atom_pairs'],count)
            self.assertEqual(summary['backbone_backbone_clashes'],count)
        with self.assertRaisesRegex(ValueError,'elements'):
            app.pair_geometry([dict(atom,element='XX')],[atom],app.Options())

    def test_nested_zip_and_pipeline(self):
        import contextlib
        out,meta=self.assemble()
        inner=io.BytesIO()
        with zipfile.ZipFile(inner,'w') as z:
            for path in out.iterdir():z.write(path,'complex/'+path.name)
        archive=self.root/'outer.zip'
        with zipfile.ZipFile(archive,'w') as z:z.writestr('inner.zip',inner.getvalue())
        locator=f'{archive}::inner.zip!/complex'
        self.assertEqual(load_structure(locator).identity['fold_id'],meta['fold_id'])
        script=load_script(Path(__file__).resolve().parents[2]/'scripts/alphafold-pdockq2-patches.py')
        record=f"{locator}\t{meta['fold_id']}\t0\tNA\tNA\tNA\n"
        with patch('sys.stdin',io.StringIO(record)), contextlib.redirect_stdout(io.StringIO()) as printed:
            self.assertEqual(script.main(['-','--no-dssp','--min-residues-per-chain','1','--min-contacts','1']),0)
        self.assertIn(meta['fold_id'],printed.getvalue())
        self.assertIn('pDockQ2 NA',printed.getvalue())

    def test_cli_errors_and_missing_primary_model(self):
        import contextlib
        script=load_script(Path(__file__).resolve().parents[2]/'scripts/alphafold-overlap.py')
        with contextlib.redirect_stderr(io.StringIO()) as error:
            code=script.main([str(self.a),str(self.b),'--output',str(self.root/'bad'),'--model-a','4'])
        self.assertEqual(code,1)
        self.assertIn('ERROR — cannot assemble',error.getvalue())
        self.assertFalse((self.root/'bad').exists())

    def test_specific_assessment_messages(self):
        self.assertEqual(app.assessment_status(), 'No detected geometric conflict')
        self.assertEqual(app.assessment_status(replacement=True),
                         'WARNING — Partners do not clash; merged shared-chain interface requires refinement')
        for field, phrase in [('alignment','alignment requires review'),
                              ('cross_clashes','cross-partner clashes detected'),
                              ('alternative_clashes','alternative alignments'),
                              ('inherited','input assembly warnings')]:
            for replacement in (False,True):
                with self.subTest(field=field,replacement=replacement):
                    status=app.assessment_status(**{field:True},replacement=replacement)
                    self.assertIn(phrase,status)
                    self.assertNotIn('Partners do not clash',status)
                    self.assertEqual('interface requires refinement' in status,replacement)

    def test_combination_summary_identifies_selected_pair_and_table(self):
        output=io.StringIO()
        app.run(self.a,self.b,self.root/'summary',all_models=True,stdout=output)
        text=output.getvalue()
        self.assertIn('Selected models: A=0, B=0',text)
        self.assertIn('Model combinations: 1 assessed;',text)
        self.assertIn('combinations with warnings; 0 errors.',text)
        self.assertIn(f"All combinations: {(self.root/'summary').resolve() / 'comparisons.tsv'}",text)

    def test_input_identities_in_stdout_and_readme(self):
        # Selected nonzero models, with --all-models, must identify the saved pair.
        shutil.copyfile(self.a/'fold_model_0.cif', self.a/'fold_model_2.cif')
        shutil.copyfile(self.b/'fold_model_0.cif', self.b/'fold_model_3.cif')
        stream=io.StringIO()
        out=self.root/'identities'
        app.run(self.a,self.b,out,model_a=2,model_b=3,all_models=True,stdout=stream)
        for label,path,number in [('A',self.a,2),('B',self.b,3)]:
            identity=load_structure(path,number).identity
            expected=f"Input {label}: Fold ID {identity['fold_id']} | Sequence ID {identity['sequence_id']} | Model {number}"
            self.assertIn(expected,stream.getvalue())
            self.assertIn(expected,(out/'README.md').read_text())
        self.assertLess(stream.getvalue().index('Locator:'),stream.getvalue().index('Input A:'))

    def test_recursive_and_unknown_input_identities(self):
        out,meta=self.assemble()
        third=self.make('identity-third',[('A','DDDD',[7,0,0]),('B','EEEE',[14,0,0])])
        plain=self.root/'standalone.cif'
        plain.write_text((third/'fold_model_0.cif').read_text())
        stream=io.StringIO();final=self.root/'identity-recursive'
        app.run(out,plain,final,stdout=stream)
        expected_a=f"Input A: Fold ID {meta['fold_id']} | Sequence ID {meta['sequence_id']} | Model 0"
        expected_b=f"Input B: Fold ID NA | Sequence ID {load_structure(plain).identity['sequence_id']} | Model 0"
        for expected in [expected_a,expected_b]:
            self.assertIn(expected,stream.getvalue())
            self.assertIn(expected,(final/'README.md').read_text())

class ExplicitSequenceTests(unittest.TestCase):
    setUp = OverlapTests.setUp
    make = OverlapTests.make
    assemble = OverlapTests.assemble

    def test_first_copy_uses_cif_order_and_preserves_other_copies(self):
        left=self.make('duplicates-left',[('Z','CCCC',[0,0,0]),('A','CCCC',[20,0,0]),('D','DDDD',[-7,0,0])])
        right=self.make('duplicates-right',[('Y','CCCC',[0,0,0]),('B','CCCC',[40,0,0]),('D','DDDD',[7,0,0])])
        a,b=load_structure(left),load_structure(right)
        with self.assertRaises(ValueError):app.unique_shared(a,b)
        self.assertEqual(app.unique_shared(a,b,'CCCC'),('Z','Y'))
        out=self.root/'explicit';stream=io.StringIO()
        meta=app.run(left,right,out,sequence='CCCC',stdout=stream)
        self.assertEqual(meta['identity']['shared_chains'],['Z','Y'])
        self.assertEqual(meta['alignment_selection']['requested_sequence'],'CCCC')
        self.assertIn('Alignment chains: input A=Z, input B=Y',stream.getvalue())
        self.assertIn('Alignment chains: input A=Z, input B=Y',(out/'README.md').read_text())
        values=list(load_structure(out).sequences.values())
        self.assertEqual(values.count('CCCC'),3)
        self.assertEqual(values.count('DDDD'),2)

    def test_missing_invalid_sequence_and_no_output(self):
        for sequence,label in [('AAAA','B'),('DDDD','A'),('GGGG','A')]:
            with self.subTest(sequence=sequence):
                out=self.root/'missing'
                with self.assertRaisesRegex(ValueError,f'not found in input {label}'):
                    app.run(self.a,self.b,out,sequence=sequence,stdout=io.StringIO())
                self.assertFalse(out.exists())
        for sequence in ('','cccc','CC CC','CCCC*','X'):
            with self.subTest(sequence=sequence),self.assertRaises(ValueError):
                app.unique_shared(load_structure(self.a),load_structure(self.b),sequence)

    def test_explicit_same_selection_preserves_assembly_ids(self):
        _,automatic=self.assemble('auto')
        _,explicit=self.assemble('specified',sequence='CCCC')
        self.assertEqual(automatic['fold_id'],explicit['fold_id'])
        self.assertEqual(automatic['sequence_id'],explicit['sequence_id'])
        self.assertEqual(automatic['alignment_selection']['mode'],'automatic')

    def test_all_models_selects_independently_and_records_missing_sequence(self):
        (self.a/'fold_model_1.cif').write_text(fixture([('Z','CCCC',[0,0,0]),('A','AAAA',[-7,0,0])]))
        (self.b/'fold_model_1.cif').write_text(fixture([('X','GGGG',[0,0,0]),('Y','DDDD',[7,0,0])]))
        out,_=self.assemble('all-explicit',sequence='CCCC',all_models=True)
        with (out/'comparisons.tsv').open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
        self.assertEqual(len(rows),4)
        self.assertEqual(sum(r['status'].startswith('ERROR') for r in rows),2)
        self.assertEqual(next(r for r in rows if r['model_a']=='1' and r['model_b']=='0')['shared_chain_a'],'Z')
        self.assertTrue(all('not found in input B' in r['warnings'] for r in rows if r['model_b']=='1'))

    def test_cli_sequence_option(self):
        script=load_script(Path(__file__).resolve().parents[2]/'scripts/alphafold-overlap.py')
        with patch('sys.stdout',new_callable=io.StringIO):
            self.assertEqual(script.main([str(self.a),str(self.b),'--sequence','CCCC',
                                          '--output',str(self.root/'cli-explicit')]),0)

class SubsequenceTests(unittest.TestCase):
    setUp = OverlapTests.setUp
    make = OverlapTests.make

    def test_different_offsets_and_retained_flanks(self):
        a=self.make('sub-a',[('Z','ACCC',[0,0,0]),('A','DDDD',[-10,0,0])])
        b=self.make('sub-b',[('Y','CCCA',[0,0,0]),('B','GGGG',[10,0,0])])
        first,second=load_structure(a),load_structure(b)
        result=app.compare(first,second,app.Options(),sequence='CCC')
        summary=result['summary']
        self.assertEqual((summary['alignment_start_a'],summary['alignment_end_a']),(2,4))
        self.assertEqual((summary['alignment_start_b'],summary['alignment_end_b']),(1,3))
        self.assertEqual(summary['matched_CA'],3)
        self.assertEqual(summary['coverage'],1)
        self.assertAlmostEqual(summary['CA_RMSD_A'],0,places=6)
        self.assertEqual([r['source_position_a'] for r in result['deviations']],[2,3,4])
        self.assertEqual([r['source_position_b'] for r in result['deviations']],[1,2,3])
        self.assertTrue(any('flanks are discarded' in w for w in result['warnings']))
        out=self.root/'sub-output';stream=io.StringIO()
        meta=app.run(a,b,out,sequence='CCC',stdout=stream)
        self.assertEqual(meta['alignment_selection']['range_a'],[2,4])
        self.assertEqual(meta['alignment_selection']['range_b'],[1,3])
        self.assertEqual(meta['identity']['alignment'],'matched-subsequence-position-CA')
        self.assertEqual(list(load_structure(out).sequences.values()),['ACCC','DDDD','GGGG'])
        self.assertIn('Alignment ranges: A=2-4, B=1-3',stream.getvalue())

    def test_first_chain_and_first_occurrence(self):
        a=self.make('repeat-a',[('Z','CCCC',[0,0,0]),('A','CCCA',[-10,0,0])])
        b=self.make('repeat-b',[('Y','CCCC',[0,0,0]),('B','CCCA',[10,0,0])])
        result=app.compare(load_structure(a),load_structure(b),app.Options(),sequence='CCC')
        s=result['summary']
        self.assertEqual((s['shared_chain_a'],s['shared_chain_b']),('Z','Y'))
        self.assertEqual((s['alignment_start_a'],s['alignment_start_b']),(1,1))
        self.assertEqual(s['matched_CA'],3)

    def test_short_or_missing_coordinates_no_fallback(self):
        a,b=load_structure(self.a),load_structure(self.b)
        with self.assertRaisesRegex(ValueError,'at least three'):
            app.compare(a,b,app.Options(),sequence='CC')
        a.atoms['B']=[x for x in a.atoms['B'] if x['number']!=2]
        with self.assertRaisesRegex(ValueError,'at least three'):
            app.compare(a,b,app.Options(),sequence='CCC')

    def test_replacement_contacts_use_correspondence_not_raw_numbers(self):
        a=self.make('mapped-a',[('Z','ACCC',[0,0,0]),('D','DDDD',[-100,0,0])])
        b=self.make('mapped-b',[('Y','CCCA',[0,0,0]),('G','GGGG',[0,0,0])])
        first,second=load_structure(a),load_structure(b)
        lookup={(x['number'],x['name']):x for x in second.atoms['Y']}
        for atom in first.atoms['Z']:
            if atom['number']>1:
                atom['coord']=list(lookup[(atom['number']-1,atom['name'])]['coord'])
            else:
                atom['coord']=[v+100 for v in atom['coord']]
        for atom in second.atoms['Y']:
            if atom['number']==4:atom['coord']=[v+100 for v in atom['coord']]
        result=app.compare(first,second,app.Options(clash_overlap=.5),sequence='CCC')
        change=result['changes'][0]
        self.assertGreater(change['original_clashes'],0)
        self.assertEqual(change['introduced_clashes'],0)
        self.assertEqual(change['worsened_clashes'],0)
        self.assertEqual(change['resolved_clashes'],0)
        self.assertEqual(change['contacts_lost'],0)
        self.assertEqual(change['contacts_gained'],0)
