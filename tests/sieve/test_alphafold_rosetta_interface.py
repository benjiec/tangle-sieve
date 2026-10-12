import io
import json
from types import SimpleNamespace

import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch
from sieve import alphafold_rosetta_interface as r


def test_bad_interface():
    for value in ['A_A', 'A:A', 'AA:B', 'A:', 'A:B:C', ' :B']:
        with unittest.TestCase().assertRaises(ValueError):
            r.parse_interface(value)


def test_input_contract():
    assert r.inputs('-', io.StringIO('/fold\tfc0f972f\t0\tNA\t0.7\n'), None) == [('/fold', 0)]
    assert r.inputs('/fold', None, 2) == [('/fold', 2)]
    assert r.inputs('-', io.StringIO('\n'), None) == []
    for text, model in [('/fold\n', None), ('/fold\n', -1),
                        ('/fold\tfc0f972f\t0\tNA\tNA\n', 1), ('x\ty\n', 0)]:
        with unittest.TestCase().assertRaises(ValueError):
            r.inputs('-', io.StringIO(text), model)


def scorefile(path, values=None):
    values = ['1'] * len(r.METRICS) if values is None else values
    path.write_text('SCORE: ' + ' '.join(r.METRICS) + ' description\nSCORE: ' +
                    ' '.join(values) + ' model_0001\n')


def test_scores(tmp_path):
    path = tmp_path / 'score.sc'
    scorefile(path)
    assert r.read_scores(path)['dG_separated'] == '1'
    for bad in ['nan', 'inf', 'bad']:
        scorefile(path, [bad] + ['1'] * (len(r.METRICS) - 1))
        with unittest.TestCase().assertRaises(ValueError):
            r.read_scores(path)
    for content in ['', 'SCORE: 1 2\n', 'SCORE: description\nSCORE: model\n']:
        path.write_text(content)
        with unittest.TestCase().assertRaises(ValueError):
            r.read_scores(path)
    scorefile(path)
    path.write_text(path.read_text() + path.read_text().splitlines()[1] + '\n')
    with unittest.TestCase().assertRaises(ValueError):
        r.read_scores(path)


def test_multiple_outputs(tmp_path, monkeypatch, tsv):
    def fake_run(command, **kwargs):
        return SimpleNamespace(returncode=0, stderr='', stdout=json.dumps([{'Id': 'sha256:test'}]))
    monkeypatch.setattr(r.subprocess, 'run', fake_run)
    def fake_analyze(locator, number, interface, directory, image):
        assert image == 'sha256:test'
        directory.mkdir()
        row = dict.fromkeys(r.COLUMNS, '1')
        row.update(model_id=number, fold_locator=locator)
        return row
    monkeypatch.setattr(r, 'analyze', fake_analyze)
    out = io.StringIO()
    r.run('-', io.StringIO('/fold\tfc0f972f\t0\tNA\tNA\n/fold\tfc0f972f\t1\tNA\tNA\n'),
          tmp_path, 'A:B', None, tsv, out)
    assert len(list(tmp_path.glob('*-model-*'))) == 2
    if tsv:
        assert out.getvalue() == (tmp_path / 'summary.tsv').read_text()
        assert len(out.getvalue().splitlines()) == 3
    else:
        assert 'model_id: 0\n' in out.getvalue()
        assert 'model_id: 1\n' in out.getvalue()
    with unittest.TestCase().assertRaisesRegex(ValueError, 'never overwritten'):
        r.run('/fold', None, tmp_path, 'A:B', 0, tsv, out)


def test_empty_and_docker_failure(tmp_path, monkeypatch):
    out = io.StringIO()
    r.run('-', io.StringIO(''), tmp_path, 'A:B', None, True, out)
    assert out.getvalue().startswith('fold_id\t')
    monkeypatch.setattr(r.subprocess, 'run', lambda *a, **k: SimpleNamespace(returncode=1, stderr='unavailable'))
    with unittest.TestCase().assertRaisesRegex(ValueError, 'cannot inspect'):
        r.run('/fold', None, tmp_path, 'A:B', 0, False, out)
    assert (tmp_path / 'image-inspect.log').read_text() == 'unavailable'

class MonkeyPatch:
    def __init__(self, case):
        self.case = case
    def setattr(self, obj, name, value):
        patcher = patch.object(obj, name, value)
        patcher.start()
        self.case.addCleanup(patcher.stop)

class InterfaceTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.monkey = MonkeyPatch(self)
    def test_interface(self):
        test_bad_interface()
    def test_inputs(self):
        test_input_contract()
    def test_scores(self):
        test_scores(self.root)
    def test_text(self):
        test_multiple_outputs(self.root, self.monkey, False)
    def test_tsv(self):
        test_multiple_outputs(self.root, self.monkey, True)
    def test_failures(self):
        test_empty_and_docker_failure(self.root, self.monkey)

class AnalysisTests(unittest.TestCase):
    def test_execution_and_failure_artifacts(self):
        from contextlib import contextmanager
        @contextmanager
        def opened(locator):
            yield SimpleNamespace(models=['x_model_0.cif'], member_names=(), read=lambda n: 'cif')
        chains = [SimpleNamespace(id=c, get_residues=lambda: iter([1, 2])) for c in 'AB']
        with tempfile.TemporaryDirectory() as temp, \
             patch.object(r, 'open_result_fold', opened), \
             patch.object(r, 'chain_sequences', return_value={'A': 'AA', 'B': 'GG'}), \
             patch.object(r, 'read_seed', return_value=1), \
             patch.object(r, 'MMCIFParser') as parser:
            parser.return_value.get_structure.return_value = [chains]
            root = Path(temp)
            def execute(command, **kwargs):
                self.assertNotIn('-interface', command)
                scorefile(root / 'success' / 'interface.sc')
                return SimpleNamespace(returncode=0)
            with patch.object(r.subprocess, 'run', execute):
                row = r.analyze('locator', 0, 'A:B', root / 'success', 'image')
            self.assertEqual(row['dG_separated'], '1')
            settings = json.loads((root / 'success/settings.json').read_text())
            self.assertEqual(settings['status'], 'complete')
            self.assertEqual(settings['chain_residue_counts'], {'A': 2, 'B': 2})
            with patch.object(r.subprocess, 'run', return_value=SimpleNamespace(returncode=2)):
                with self.assertRaisesRegex(ValueError, 'exited 2'):
                    r.analyze('locator', 0, 'A:B', root / 'failure', 'image')
            self.assertEqual(json.loads((root / 'failure/settings.json').read_text())['status'], 'failed')
            with self.assertRaisesRegex(ValueError, 'not present'):
                r.analyze('locator', 0, 'A:C', root / 'badchain', 'image')
            with self.assertRaisesRegex(ValueError, 'exactly one model'):
                r.analyze('locator', 9, 'A:B', root / 'missing', 'image')

class TemporaryOutputTests(unittest.TestCase):
    def test_success_cleanup_and_summary(self):
        seen = []
        def execute(source, stream, directory, interface, model, tsv, stdout):
            root = Path(directory)
            seen.append(root)
            (root / 'result').write_text('test')
            print('summary', file=stdout)
        out = io.StringIO()
        with patch.object(r, '_run', execute):
            r.run('locator', None, None, 'A:B', 0, False, out)
        self.assertEqual(out.getvalue(), 'summary\n')
        self.assertFalse(seen[0].exists())

    def test_failure_cleanup_and_diagnostics(self):
        seen = []
        def execute(*args):
            root = Path(args[2])
            seen.append(root)
            (root / 'image-inspect.log').write_text('Docker daemon unavailable')
            raise ValueError('cannot inspect image')
        with patch.object(r, '_run', execute):
            with self.assertRaisesRegex(ValueError, 'Docker daemon unavailable'):
                r.run('locator', None, None, 'A:B', 0, False, io.StringIO())
        self.assertFalse(seen[0].exists())

    def test_interrupt_cleanup(self):
        seen = []
        def execute(*args):
            seen.append(Path(args[2]))
            raise KeyboardInterrupt
        with patch.object(r, '_run', execute):
            with self.assertRaises(KeyboardInterrupt):
                r.run('locator', None, None, 'A:B', 0, False, io.StringIO())
        self.assertFalse(seen[0].exists())

    def test_cli_output_optional(self):
        from tests.scripts.helpers import load_script
        script = load_script('alphafold-rosetta-interface.py')
        with patch.object(script, 'run') as run:
            self.assertEqual(script.main(['locator', '--interface', 'A:B', '--model', '0']), 0)
        self.assertIsNone(run.call_args.args[2])
