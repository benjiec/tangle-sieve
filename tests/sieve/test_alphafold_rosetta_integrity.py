import io
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from sieve import alphafold_rosetta_integrity as app
from sieve import rosetta_integrity_worker as worker
from tests.scripts.helpers import load_script


def result():
    base = dict.fromkeys(worker.MEASURES, 1.0)
    base.update(nonpolar_sasa_in_complex=1.0, nonpolar_shielded_by_partners=0.0)
    return dict(pyrosetta_version='test', summaries=[
        dict(base, scope='complex', chain='ALL', nres=2),
        dict(base, scope='chain', chain='A', nres=1),
        dict(base, scope='chain', chain='B', nres=1)])


class IntegrityTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def test_classification_and_finite(self):
        for element, charge, expected in [('C', 0, True), ('S', -.4, True),
                                          ('C', .5, False), ('N', 0, False), ('H', 0, False)]:
            self.assertEqual(worker.nonpolar(element, charge), expected)
        for bad in [float('nan'), float('inf'), -float('inf')]:
            with self.assertRaises(ValueError):
                worker.check_finite({'nested': [bad]})

    def test_validate_results(self):
        self.assertEqual(len(app.validate_result(result())), 3)
        for modify in [lambda r: r['summaries'][1].pop('packstat'),
                       lambda r: r['summaries'][1].update(fa_rep=float('nan')),
                       lambda r: r['summaries'][1].update(nres=0),
                       lambda r: r['summaries'][1].update(scope='complex'),
                       lambda r: r['summaries'][2].update(chain='A'),
                       lambda r: r['summaries'][0].update(nres=10)]:
            data = result()
            modify(data)
            with self.assertRaises(ValueError):
                app.validate_result(data)

    def fake_execute(self, locator, number, directory, image):
        self.assertEqual(image, 'sha256:pinned')
        directory.mkdir()
        return [dict(row, fold_id='12345678', sequence_id='abcdefgh', model_id=number,
                     fold_locator=locator, score_function='ref2015') for row in result()['summaries']]

    def test_pipeline_tsv_and_overwrite(self):
        output = io.StringIO()
        inspected = SimpleNamespace(returncode=0, stdout=json.dumps([{'Id': 'sha256:pinned'}]), stderr='')
        with patch.object(app.subprocess, 'run', return_value=inspected), patch.object(app, 'execute', self.fake_execute):
            app.run('-', io.StringIO('loc\t12345678\t3\tNA\tNA\nloc\t12345678\t4\tNA\tNA\n'),
                    self.root, None, True, output)
        self.assertEqual(len(output.getvalue().splitlines()), 7)
        self.assertEqual(output.getvalue(), (self.root/'summary.tsv').read_text())
        with self.assertRaisesRegex(ValueError, 'never overwritten'):
            app.run('loc', None, self.root, 0, False, io.StringIO())

    def test_temporary_success_failure_interrupt(self):
        for exception in [None, ValueError('worker failed'), KeyboardInterrupt()]:
            seen = []
            def records(items, directory, tsv, stdout):
                seen.append(directory)
                (directory/'run.log').write_text('specific diagnostic')
                if exception is not None:
                    raise exception
                print('summary', file=stdout)
            out = io.StringIO()
            with patch.object(app, 'run_records', records):
                if isinstance(exception, ValueError):
                    with self.assertRaisesRegex(ValueError, 'specific diagnostic'):
                        app.run('loc', None, None, 0, False, out)
                elif exception:
                    with self.assertRaises(KeyboardInterrupt):
                        app.run('loc', None, None, 0, False, out)
                else:
                    app.run('loc', None, None, 0, False, out)
                    self.assertEqual(out.getvalue(), 'summary\n')
            self.assertFalse(seen[0].exists())

    def test_empty_stdin_no_docker(self):
        with patch.object(app.subprocess, 'run') as run:
            out = io.StringIO()
            app.run('-', io.StringIO(''), None, None, True, out)
        run.assert_not_called()
        self.assertEqual(len(out.getvalue().splitlines()), 1)

    def test_cli_and_missing_model(self):
        script = load_script('scripts/alphafold-rosetta-integrity.py')
        with patch.object(script, 'run') as run:
            self.assertEqual(script.main(['-', '--tsv']), 0)
        self.assertIsNone(run.call_args.args[2])
        with self.assertRaisesRegex(ValueError, 'requires'):
            app.run('loc', None, None, None, False, io.StringIO())
        with self.assertRaisesRegex(ValueError, 'conflicts'):
            app.run('-', io.StringIO('loc\t12345678\t1\tNA\tNA\n'), None, 0, False, io.StringIO())

    def test_execution_artifacts(self):
        def prepare(locator, number, directory):
            directory.mkdir()
            (directory/'model.cif').write_text('CIF')
            (directory/'worker.py').write_text('worker')
            return dict(fold_id='12345678', sequence_id='12345678', model_id=0,
                        fold_locator=locator, score_function='ref2015')
        directory = self.root/'success'
        def execute(command, **kwargs):
            (directory/'integrity.json').write_text(json.dumps(result()))
            return SimpleNamespace(returncode=0)
        with patch.object(app, 'prepare', prepare), patch.object(app.subprocess, 'run', execute):
            self.assertEqual(len(app.execute('loc', 0, directory, 'image')), 3)
        settings = json.loads((directory/'settings.json').read_text())
        self.assertEqual(settings['status'], 'complete')
        self.assertEqual(len(settings['input_sha256']), 64)
        directory = self.root/'failure'
        with patch.object(app, 'prepare', prepare), patch.object(app.subprocess, 'run', return_value=SimpleNamespace(returncode=9)):
            with self.assertRaisesRegex(ValueError, 'exited 9'):
                app.execute('loc', 0, directory, 'image')
        self.assertEqual(json.loads((directory/'settings.json').read_text())['status'], 'failed')
