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
        self.assertEqual(next(iter(cached.values()))['source_type'],'assembled')
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

    def test_legacy_manifest_upgrade(self):
        out,meta=self.assemble()
        fasta=self.root/'candidates.faa';fasta.write_text('>a\nAAAA\n>c\nCCCC\n>d\nDDDD\n')
        path=self.root/'cache.tsv'
        self.assertEqual(manifest(self.a,fasta,path,io.StringIO()),0)
        with path.open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
        old_id=rows[0]['fold_id'];columns=[c for c in rows[0] if c!='source_type']
        with path.open('w') as stream:
            w=csv.DictWriter(stream,columns,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
        self.assertEqual(manifest(out,fasta,path,io.StringIO()),0)
        cached,_=load_manifest(path)
        self.assertEqual(len(cached),2)
        self.assertEqual(cached[(str(self.a.resolve()),'.')]['fold_id'],old_id)
        self.assertEqual(cached[(str(self.a.resolve()),'.')]['source_type'],'alphafold')

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
