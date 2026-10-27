import copy
import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

from Bio.PDB import MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from sieve.alphafold_composition import ChainComposition, chain_compositions, hash_payload, STANDARD
from sieve.alphafold_files import open_fold
from sieve.alphafold_manifest import sequence_id, fold_id, chain_sequences, run
from sieve.alphafold_structure import identity_for_fold, open_structure_fold, load_structure
from sieve.alphafold_pdockq2_patches import analyze, Options
from sieve.alphafold_results import run as list_results
from sieve.alphafold_rosetta_integrity import prepare
from sieve.alphafold_overlap import compare, Options as OverlapOptions
from tests.sieve.test_alphafold_overlap import fixture

THREE = {v:k for k,v in STANDARD.items()}


def cif(specs, parents=None):
    """Four-residue chains with a real atom-to-polymer map, including HETATM PTMs."""
    data = MMCIF2Dict(io.StringIO(fixture([(c,'AAAA',[i*6,0,0]) for i,c in enumerate(specs)])))
    entities = dict(zip(data['_struct_asym.entity_id'],data['_struct_asym.id']))
    def code(chain, position):
        seq,mods = specs[chain]
        return dict(mods).get(position,THREE[seq[position-1]])
    for i,(entity,position) in enumerate(zip(data['_entity_poly_seq.entity_id'],data['_entity_poly_seq.num'])):
        data['_entity_poly_seq.mon_id'][i] = code(entities[entity],int(position))
    for i,(chain,position) in enumerate(zip(data['_atom_site.label_asym_id'],data['_atom_site.label_seq_id'])):
        data['_atom_site.label_comp_id'][i] = code(chain,int(position))
        if int(position) in dict(specs[chain][1]):data['_atom_site.group_PDB'][i] = 'HETATM'
    if parents:
        data['_chem_comp.id'] = list(parents)
        data['_chem_comp.mon_nstd_parent_comp_id'] = list(parents.values())
    writer=MMCIFIO();writer.set_dict(data);out=io.StringIO();writer.save(out);return out.getvalue()


def server_request(entries):
    return dict(modelSeeds=['42'], sequences=[dict(proteinChain=dict(sequence=s,count=n,
        modifications=[dict(ptmPosition=p,ptmType='CCD_'+c) for p,c in mods])) for s,mods,n in entries])


class PTMTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name).resolve()
        self.specs={'A':('ASST',((2,'SEP'),(3,'SEP'))),'B':('ASST',((2,'SEP'),(3,'SEP')))}
        self.request=server_request([('ASST',((2,'SEP'),(3,'SEP')),2)])

    def folder(self, name='fold', specs=None, request=None):
        path=self.root/name;path.mkdir()
        (path/'fold_job_request.json').write_text(json.dumps([self.request if request is None else request]))
        for n in range(2):(path/f'fold_model_{n}.cif').write_text(cif(self.specs if specs is None else specs))
        return path

    def test_legacy_ids_byte_for_byte_and_empty_ptms(self):
        for seqs in [[],['ASST'],['ASST','CCCC','ASST']]:
            payload=json.dumps(sorted(seqs),separators=(',',':'),ensure_ascii=True)
            self.assertEqual(sequence_id(seqs),hashlib.sha1(payload.encode()).hexdigest()[-8:])
            seeded=json.dumps((42,sorted(seqs)),separators=(',',':'),ensure_ascii=True)
            self.assertEqual(fold_id(42,seqs),hashlib.sha1(seeded.encode()).hexdigest()[-8:])
            self.assertEqual(sequence_id([ChainComposition(s) for s in seqs]),sequence_id(seqs))

    def test_ptm_ids_order_multiplicity_position_type_and_assignment(self):
        a=ChainComposition('ASST',((2,'SEP'),(3,'SEP')))
        b=ChainComposition('ASST',((3,'CCD_SEP'),(2,'CCD_SEP')))
        self.assertEqual(sequence_id([a,a,'CCCC']),sequence_id(['CCCC',b,b]))
        self.assertNotEqual(sequence_id([a]),sequence_id(['ASST']))
        variants=[ChainComposition('ASST',((2,'SEP'),)),ChainComposition('ASST',((3,'SEP'),)),
                  ChainComposition('ASST',((4,'TPO'),)),a]
        self.assertEqual(len({sequence_id([c]) for c in variants}),4)
        self.assertNotEqual(sequence_id([a,a]),sequence_id([a,'ASST']))
        self.assertNotEqual(sequence_id([ChainComposition('ASST',((2,'SEP'),)),'AASS']),
                            sequence_id(['ASST',ChainComposition('AASS',((3,'SEP'),))]))
        self.assertEqual(sequence_id([a]),sequence_id([b]))
        self.assertNotEqual(fold_id(42,[a]),fold_id(42,['ASST']))
        self.assertNotEqual(fold_id(42,[a]),fold_id(43,[a]))

    def test_server_copies_and_explicit_ids(self):
        text=cif(self.specs)
        compositions=chain_compositions(text,self.request)
        self.assertEqual(chain_sequences(text,self.request),{'A':'ASST','B':'ASST'})
        self.assertEqual(compositions['A'].modifications,((2,'SEP'),(3,'SEP')))
        request=dict(sequences=[dict(protein=dict(id=['B','A'],sequence='ASST',modifications=[
            dict(ptmType='SEP',ptmPosition=3),dict(ptmType='SEP',ptmPosition=2)]))])
        self.assertEqual(chain_compositions(text,request),compositions)
        # Reordered server sequence entries and chains have no effect on IDs.
        specs={'Z':('CCCC',()),'Q':('ASST',((3,'SEP'),))}
        req=server_request([('ASST',((3,'SEP'),),1),('CCCC',(),1)])
        expected=sequence_id(chain_compositions(cif(specs),req).values())
        req['sequences'].reverse()
        self.assertEqual(sequence_id(chain_compositions(cif(specs),req).values()),expected)

    def test_bad_modifications_and_mismatch_rejected(self):
        bad=[None,[dict(ptmPosition=0,ptmType='SEP')],[dict(ptmPosition=5,ptmType='SEP')],
             [dict(ptmPosition=True,ptmType='SEP')],[dict(ptmPosition=2.0,ptmType='SEP')],
             [dict(ptmPosition=2,ptmType='')],[dict(ptmPosition=2,ptmType=' SEP')],
             [dict(ptmPosition=2,ptmType='SEP')]*2]
        for mods in bad:
            request=copy.deepcopy(self.request);request['sequences'][0]['proteinChain']['modifications']=mods
            with self.subTest(mods=mods),self.assertRaises(ValueError):chain_compositions(cif(self.specs),request)
        for field,value in [('count',1),('count',True),('sequence','ATST'),('modifications',[])]:
            request=copy.deepcopy(self.request);request['sequences'][0]['proteinChain'][field]=value
            with self.subTest(field=field),self.assertRaises(ValueError):chain_compositions(cif(self.specs),request)
        text=cif(self.specs)
        data=MMCIF2Dict(io.StringIO(text));data['_atom_site.label_comp_id'][3]='SER'
        writer=MMCIFIO();writer.set_dict(data);out=io.StringIO();writer.save(out)
        with self.assertRaisesRegex(ValueError,'atom residue/PTM'):chain_compositions(out.getvalue(),self.request)

    def test_unknown_parent_and_ambiguous_server_assignment(self):
        specs={'A':('ASST',((2,'MYMOD'),))}
        request=server_request([('ASST',((2,'MYMOD'),),1)])
        with self.assertRaisesRegex(ValueError,'unsupported residue'):chain_compositions(cif(specs))
        self.assertEqual(chain_compositions(cif(specs),request)['A'].sequence,'ASST')
        self.assertEqual(chain_compositions(cif(specs,{'MYMOD':'SER'}))['A'].sequence,'ASST')
        ambiguous=server_request([('ASST',((2,'MYMOD'),),1),('ATST',((2,'MYMOD'),),1)])
        with self.assertRaisesRegex(ValueError,'ambiguous'):chain_compositions(cif(specs),ambiguous)
        with self.assertRaisesRegex(ValueError,'parent'):chain_compositions(cif(specs,{'MYMOD':'SER,THR'}))

    def test_manifest_identity_fasta_and_cache(self):
        folder=self.folder();fasta=self.root/'candidates.faa';fasta.write_text('>base sequence header\nASST\n')
        output=self.root/'manifest.tsv';errors=io.StringIO()
        self.assertEqual(run(folder,fasta,output,errors),0,errors.getvalue())
        with output.open() as stream:row=next(csv.DictReader(stream,delimiter='\t'))
        self.assertEqual(row['input 1 description'],'base sequence header');self.assertEqual(row['input 1 copies'],'2')
        expected=sequence_id(chain_compositions(cif(self.specs),self.request).values())
        self.assertEqual(row['sequence_id'],expected)
        with open_structure_fold(folder) as fold:
            identity=identity_for_fold(fold,fold.read(fold.models[0]))
        self.assertEqual(row['fold_id'],identity['fold_id'])
        saved=output.read_bytes()
        self.assertEqual(run(folder,fasta,output,io.StringIO()),0);self.assertEqual(output.read_bytes(),saved)
        missing=self.root/'missing.faa';missing.write_text('>other\nCCCC\n')
        errors=io.StringIO()
        self.assertEqual(run(folder,missing,self.root/'bad.tsv',errors),1)
        self.assertIn('ASST not found',errors.getvalue())
        listing=io.StringIO();self.assertEqual(list_results(output,expected,listing),0)
        self.assertEqual(len(listing.getvalue().splitlines()),2)

    def test_all_models_agree_and_nested_zip(self):
        folder=self.folder();fasta=self.root/'candidates.faa';fasta.write_text('>base\nASST\n')
        inner=io.BytesIO()
        with zipfile.ZipFile(inner,'w') as z:
            for p in folder.iterdir():z.writestr('fold/'+p.name,p.read_bytes())
        archive=self.root/'outer.zip'
        with zipfile.ZipFile(archive,'w') as z:z.writestr('inner.zip',inner.getvalue())
        out=self.root/'zip.tsv';self.assertEqual(run(archive,fasta,out,io.StringIO()),0)
        with open_fold(f'{archive}::inner.zip!/fold') as fold:
            identity=identity_for_fold(fold,fold.read(fold.models[0]))
        with out.open() as stream:self.assertEqual(next(csv.DictReader(stream,delimiter='\t'))['fold_id'],identity['fold_id'])
        different=dict(self.specs);different['B']=('ASST',((2,'SEP'),))
        (folder/'fold_model_1.cif').write_text(cif(different))
        errors=io.StringIO();self.assertEqual(run(folder,fasta,self.root/'bad.tsv',errors),1)
        self.assertIn('mapping',errors.getvalue())
        self.assertFalse((self.root/'bad.tsv').exists())

    def test_patch_ids_and_explicit_analysis_limits(self):
        folder=self.folder()
        for n in range(2):
            (folder/f'fold_full_data_{n}.json').write_text(json.dumps(dict(token_chain_ids=['A']*4+['B']*4,
                token_res_ids=list(range(1,5))*2,pae=[[1]*8 for _ in range(8)])))
        ids={};seqs={};analyze(folder,options=Options(min_residues_per_chain=1,min_contacts=1),fold_ids=ids,sequence_ids=seqs)
        source=load_structure(folder)
        self.assertEqual(next(iter(ids.values())),source.identity['fold_id'])
        self.assertEqual(next(iter(seqs.values())),source.identity['sequence_id'])
        with self.assertRaisesRegex(ValueError,'PTM'):prepare(str(folder),0,self.root/'rosetta')
        with self.assertRaisesRegex(ValueError,'PTM'):compare(source,source,OverlapOptions())
        from sieve.patch_description_summary import PatchDescriptions
        descriptions=PatchDescriptions({('fold',0):cif(self.specs)},None)
        self.assertIn(('fold',0),descriptions.modified)

    def test_multitoken_pae_and_geometry(self):
        from sieve.alphafold_pdockq2 import parse_model, score_selections
        from sieve.patch_description_summary import PatchDescriptions
        text=cif(self.specs)
        data=MMCIF2Dict(io.StringIO(text))
        keys=list(zip(data['_atom_site.label_asym_id'],map(int,data['_atom_site.label_seq_id'])))
        tokens=[]
        for c in self.specs:
            for n in range(1,5):
                tokens.extend([(c,n)]*(keys.count((c,n)) if n in (2,3) else 1))
        size=len(tokens)
        full=dict(token_chain_ids=[c for c,n in tokens],token_res_ids=[n for c,n in tokens],
                  pae=[[float((i*3+j)%27) for j in range(size)] for i in range(size)])
        residues,pae=parse_model(text,full)
        self.assertEqual([len(residues[c]) for c in residues],[4,4])
        first=residues['A'][0:3];second=residues['B'][1:2]
        row=score_selections('A','B',first,second,pae,None,None,100)
        self.assertEqual(row['contact count'],3)
        expected=sum(sum(1/(1+(pae[i][j]/10)**2) for i in a.token_indices for j in b.token_indices)
                     /len(a.token_indices)/len(b.token_indices) for a in first for b in second)/len(first)
        self.assertAlmostEqual(row['normalized PAE 1 to 2'],expected)
        expected_reverse=sum(sum(1/(1+(pae[j][i]/10)**2) for i in a.token_indices for j in b.token_indices)
                     /len(a.token_indices)/len(b.token_indices) for a in first for b in second)/len(first)
        self.assertAlmostEqual(row['normalized PAE 2 to 1'],expected_reverse)
        self.assertEqual(row['interface pLDDT 1'],sum(r.plddt for r in first)/len(first))
        perm=list(reversed(range(size)))
        reordered=dict(token_chain_ids=[full['token_chain_ids'][i] for i in perm],
                       token_res_ids=[full['token_res_ids'][i] for i in perm],
                       pae=[[pae[i][j] for j in perm] for i in perm])
        r,p=parse_model(text,reordered)
        again=score_selections('A','B',r['A'][0:3],r['B'][1:2],p,None,None,100)
        self.assertAlmostEqual(row['pDockQ2 max'],again['pDockQ2 max'])
        descriptions=PatchDescriptions({('fold',0):text},None)
        row.update(source='fold',model=0)
        out=io.StringIO();descriptions.write_patch(row,out)
        self.assertIn('[SEP]',out.getvalue());self.assertIn('PTM atoms retained',out.getvalue())
        self.assertIn('chemistry',out.getvalue());self.assertNotIn('nominal charge +',out.getvalue())
        out=io.StringIO();descriptions.write_chains('fold',out)
        self.assertIn('secondary structure unavailable',out.getvalue())
        self.assertIn('| 4 residues',out.getvalue())
        self.assertIn('4/4 (100.0%)',out.getvalue())
        bad=copy.deepcopy(full);bad['token_res_ids'][1]=1
        with self.assertRaisesRegex(ValueError,'duplicate token mapping'):parse_model(text,bad)
        bad=copy.deepcopy(full);bad['token_res_ids'][2]=3
        with self.assertRaisesRegex(ValueError,'PTM token count'):parse_model(text,bad)
        bad=copy.deepcopy(full);bad['pae'].pop()
        with self.assertRaisesRegex(ValueError,'dimensions'):parse_model(text,bad)
