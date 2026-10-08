import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import zipfile
from tests.scripts.helpers import load_script
from tests.sieve.test_alphafold_pdockq2_patches import CIF, DATA

SCRIPT=load_script(str(Path(__file__).resolve().parents[2]/'scripts/alphafold-interface-compare.py'))


class CommandTests(unittest.TestCase):
    def test_end_to_end_small_zip_and_existing_output_protection(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); path=root/'model.zip'; out=root/'report'
            with zipfile.ZipFile(path,'w') as z:
                z.writestr('fold_model_0.cif',CIF)
                z.writestr('fold_full_data_0.json',json.dumps(DATA))
            args=[str(path),str(path),'--reference-chains','A','B','--query-chains','A','B','--no-dssp','--output-dir',str(out)]
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(SCRIPT.main(args),0)
            comp=json.loads((out/'comparison.json').read_text())
            self.assertEqual(comp['primary']['comparisons'][0]['f1'],1)
            self.assertIsNone(comp['primary']['comparisons'][0]['anchor_rmsd'])
            self.assertEqual(len(comp['provenance']['inputs'][0]['sha256']),64)
            ref=json.loads((out/'reference.json').read_text())
            self.assertEqual(ref[0]['dssp']['status'],'unavailable')
            self.assertTrue(ref[0]['residues'][0]['missing_backbone_atoms'])
            before=(out/'report.md').read_text()
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                SCRIPT.main(args)
            self.assertEqual((out/'report.md').read_text(),before)

    def test_help(self):
        with contextlib.redirect_stdout(io.StringIO()) as out, self.assertRaises(SystemExit):
            SCRIPT.main(['--help'])
        self.assertIn('interface-comparison.md',out.getvalue())

if __name__=='__main__': unittest.main()
