import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from sieve.alphafold_find import run, selected_folds
from sieve.alphafold_manifest import INPUT_COLUMNS
from sieve.alphafold_result_format import read_patch_inputs
from tests.scripts.helpers import load_script


class FindTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(); self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.manifest = self.root/'manifest.tsv'
        self.names = self.root/'names.txt'; self.names.write_text('B\n')
        self.out, self.err = io.StringIO(), io.StringIO()

    def row(self, name, inputs, sequence='aaaaaaaa', fold='11111111', scores=None):
        folder = self.root/name; folder.mkdir()
        for model, values in (scores or {0: {'iptm': .5, 'ptm': .6}}).items():
            (folder/f'fold_model_{model}.cif').write_text('coordinates must not be read')
            (folder/f'fold_summary_confidences_{model}.json').write_text(json.dumps(values))
        row = dict.fromkeys(INPUT_COLUMNS, '')
        row.update(fold_locator=str(folder), fold_id=fold, sequence_id=sequence)
        row['_inputs'] = inputs
        return row

    def write(self, rows):
        with self.manifest.open('w', newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=['fold_locator','fold_id','sequence_id',*INPUT_COLUMNS],delimiter='\t')
            writer.writeheader()
            for row in rows:
                chain = 0
                for description,copies in row.get('_inputs', [(row.get('input description',''),1)]):
                    for _ in range(copies):
                        current={k:v for k,v in row.items() if k != '_inputs'}
                        current['input description']=row.get('input description') or description
                        current['input chain']=row.get('input chain') or chr(65+chain)
                        writer.writerow(current);chain += 1

    def execute(self, **kwargs):
        code=run(self.manifest,self.names,self.out,self.err,**kwargs)
        return code,list(csv.reader(io.StringIO(self.out.getvalue()),delimiter='\t'))

    def test_filters_self_cross_duplicates_and_limits(self):
        self.names.write_text('A full header\nB\n')
        a=self.row('a',[('A full header',2)])
        b=self.row('b',[('A full header',1),('A full header',1)])
        c=self.row('c',[('A full header',3),('B',3)])
        d=self.row('d',[('B',1),('C',1)])
        self.write([a,b,c,d])
        cases=[({},[a,b,c,d]),({'cross':True},[a,b,c]),
               ({'max_chains':2},[a,b,d]),({'max_unique_inputs':1},[a,b]),
               ({'cross':True,'max_chains':6,'max_unique_inputs':2},[a,b,c]),
               ({'cross':True,'max_chains':2,'max_unique_inputs':1},[a,b])]
        for options,expected in cases:
            with self.subTest(options=options):
                selected=selected_folds(self.manifest,{'A full header','B'},**options)
                self.assertEqual([r['fold_locator'] for r in selected],[r['fold_locator'] for r in expected])
        code,rows=self.execute();self.assertEqual(code,0);self.assertEqual(len(rows),4)
        self.assertEqual(len(list(read_patch_inputs(io.StringIO(self.out.getvalue())))),4)

    def test_exact_names_and_empty_list(self):
        rows=[self.row('spaced',[(' A full header ',1)]),self.row('exact',[('A full header',1)])]
        self.write(rows);self.names.write_text(' A full header \n\n A full header \n')
        _,out=self.execute();self.assertEqual([r[0] for r in out],[rows[0]['fold_locator']])
        self.names.write_text('');self.out=io.StringIO()
        self.assertEqual(self.execute(),(0,[]))

    def test_highest_per_sequence_and_filtered_folds(self):
        a=self.row('a',[('B',1)],scores={2:{'iptm':.8,'ptm':.1},10:{'iptm':.8,'ptm':.1}})
        b=self.row('b',[('B',1)],fold='22222222',scores={0:{'iptm':.6,'ptm':1}})
        c=self.row('c',[('B',1)],sequence='bbbbbbbb',scores={0:{'iptm':.4,'ptm':.4}})
        d=self.row('excluded',[('C',1)],scores={0:{'iptm':1,'ptm':1}})
        self.write([a,b,c,d])
        for metric,expected in [('ipTM',a),('pTM',b),('ipTM8_pTM2',b)]:
            with self.subTest(metric=metric):
                self.out=io.StringIO();code,rows=self.execute(highest=metric)
                self.assertEqual(code,0);self.assertEqual([r[0] for r in rows],[expected['fold_locator'],c['fold_locator']])
                if metric=='ipTM':self.assertEqual(rows[0][2],'2')

    def test_missing_scores_errors_and_later_valid_fold(self):
        a=self.row('a',[('B',1)],scores={0:{}})
        missing=dict(a,fold_locator=str(self.root/'absent'))
        b=self.row('b',[('B',1)],sequence='bbbbbbbb',scores={0:{'iptm':0,'ptm':0}})
        self.write([a,missing,b])
        code,rows=self.execute();self.assertEqual(code,1);self.assertEqual(rows[0][3:],['NA','NA'])
        self.out=io.StringIO();code,rows=self.execute(highest='ipTM8_pTM2')
        self.assertEqual(code,1);self.assertEqual(len(rows),1);self.assertEqual(rows[0][0],b['fold_locator'])
        self.assertIn('ERROR:',self.err.getvalue())

    def test_validation_before_fold_reads_and_conflicting_duplicates(self):
        a=self.row('a',[('B',1)])
        bads=[dict(a,fold_id='bad'),dict(a,sequence_id='bad'),
              dict(a,**{'input chain':'A'})]
        for bad in bads:
            self.write([a,bad])
            with patch('sieve.alphafold_find.confidence_rows') as read:
                with self.assertRaises(ValueError):self.execute()
                read.assert_not_called()
        self.write([a,dict(a,sequence_id='bbbbbbbb')])
        with self.assertRaisesRegex(ValueError,'contradictory'):self.execute()
        for value in (0,-1,1.5,True):
            with self.assertRaises(ValueError):selected_folds(self.manifest,{'B'},max_chains=value)
            with self.assertRaises(ValueError):selected_folds(self.manifest,{'B'},max_unique_inputs=value)

    def test_cli_defaults_and_optional_highest(self):
        a=self.row('a',[('B',1)],scores={0:{'iptm':.5,'ptm':.5},1:{'iptm':.9,'ptm':.9}})
        self.write([a])
        script=load_script(Path(__file__).resolve().parents[2]/'scripts/alphafold-find.py')
        for flags,count in [([],2),(['--highest'],1),(['--cross','--max-chains','1','--max-unique-inputs','1'],2)]:
            with patch('sys.stdout',new_callable=io.StringIO) as output:
                self.assertEqual(script.main([str(self.manifest),str(self.names),*flags]),0)
                self.assertEqual(len(output.getvalue().splitlines()),count)

    def test_multiple_chain_rows_and_highest_tie(self):
        a=self.row('a',[('B',2)])
        b=self.row('b',[('B',2)])
        self.write([a,b])
        code,rows=self.execute(highest='ipTM')
        self.assertEqual(code,0);self.assertEqual(len(rows),1)
        self.assertEqual(rows[0][0],a['fold_locator'])

    def test_malformed_headers_and_short_rows(self):
        for text in ('fold_locator\tfold_id\tsequence_id\n',
                     'fold_locator\tfold_locator\n', ''):
            self.manifest.write_text(text)
            with self.assertRaises(ValueError):self.execute()
        a=self.row('a',[('B',1)]);self.write([a])
        with self.manifest.open('a') as stream:stream.write('incomplete\n')
        with self.assertRaisesRegex(ValueError,'incomplete'):self.execute()

    def test_names_from_stdin_preserves_headers_and_stream(self):
        a=self.row('stdin-match',[(' A full header ',1)])
        b=self.row('stdin-other',[('A full header',1)])
        self.write([a,b])
        script=load_script(Path(__file__).resolve().parents[2]/'scripts/alphafold-find.py')
        stdin=io.StringIO(' A full header \r\n\n A full header \n')
        with patch('sys.stdin',stdin), patch('sys.stdout',new_callable=io.StringIO) as output:
            self.assertEqual(script.main([str(self.manifest),'-']),0)
        self.assertFalse(stdin.closed)
        rows=list(csv.reader(io.StringIO(output.getvalue()),delimiter='\t'))
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0][0],a['fold_locator'])
        with patch('sys.stdin',io.StringIO()), patch('sys.stdout',new_callable=io.StringIO) as output:
            self.assertEqual(script.main([str(self.manifest),'-']),0)
        self.assertEqual(output.getvalue(),'')

    def test_gzip_names_using_shared_reader(self):
        import gzip
        a=self.row('gzip-match',[('B',1)]);self.write([a])
        compressed=self.root/'names.gz'
        with gzip.open(compressed,'wt',encoding='utf-8') as stream:stream.write('B\n')
        self.assertEqual(run(self.manifest,compressed,self.out,self.err),0)
        self.assertIn(a['fold_locator'],self.out.getvalue())

    def test_missing_names_report_all_before_model_reads(self):
        self.write([self.row('valid',[('B',1)])])
        self.names.write_text('B\nmissing Z\nmissing A\nmissing Z\n')
        with patch('sieve.alphafold_find.confidence_rows') as read:
            code,rows=self.execute(cross=True,max_chains=1,highest='ipTM')
            read.assert_not_called()
        self.assertEqual((code,rows),(1,[]))
        self.assertEqual(self.err.getvalue().splitlines(),[
            "ERROR: sequence name not found in manifest: 'missing A'",
            "ERROR: sequence name not found in manifest: 'missing Z'"])

    def test_present_names_excluded_by_filters_are_not_errors(self):
        self.write([self.row('large',[('B',3),('C',1)])])
        for options in ({'cross':True},{'max_chains':2},{'max_unique_inputs':1}):
            self.assertEqual(self.execute(**options),(0,[]))
        self.assertEqual(self.err.getvalue(),'')

    def test_missing_names_exact_match_empty_manifest_and_stdin(self):
        self.write([self.row('exact',[('B',1)])])
        self.names.write_text('b\n B \n')
        self.assertEqual(self.execute(),(1,[]))
        self.write([])
        self.names.write_text('B\n')
        self.assertEqual(self.execute(),(1,[]))
        script=load_script(Path(__file__).resolve().parents[2]/'scripts/alphafold-find.py')
        with patch('sys.stdin',io.StringIO('missing\n')), patch('sys.stdout',io.StringIO()) as out, patch('sys.stderr',io.StringIO()) as err:
            self.assertEqual(script.main([str(self.manifest),'-']),1)
            self.assertEqual(out.getvalue(),'')
            self.assertIn('missing',err.getvalue())
