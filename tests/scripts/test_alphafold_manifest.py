import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from tests.scripts.helpers import load_script
from tests.sieve.test_alphafold_manifest import cif


class ManifestCLITests(unittest.TestCase):
    def test_cli_success_cache_and_errors(self):
        script = load_script(Path(__file__).resolve().parents[2] / 'scripts' / 'alphafold-manifest.py')
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'fold'
            source.mkdir()
            (source / 'fold_model_0.cif').write_text(cif())
            fasta = root / 'candidate.faa'
            fasta.write_text('>first full header\nAC\n>second\nDE\n')
            out = root / 'manifest.tsv'
            args = list(map(str, (source, fasta, out)))
            log = io.StringIO()
            with contextlib.redirect_stderr(log):
                self.assertEqual(script.main(args), 0)
                self.assertEqual(script.main(args), 0)
                fasta.write_text('>first\nAC\n')
                self.assertEqual(script.main(args), 0)
                out.write_text(out.read_text().splitlines()[0] + '\n')
                self.assertEqual(script.main(args), 1)
                self.assertEqual(script.main([str(source), str(root/'missing'), str(out)]), 1)
            self.assertIn('cached: 1', log.getvalue())
            self.assertIn('ERROR: DE not found in ', log.getvalue())
