"""Opt-in real Docker pipeline check; run with the sieve Python environment.

Usage: python tests/scripts/check_rosetta_integrity_pipeline.py MANIFEST SEQUENCE_ID REPORT_DIR
"""
import csv
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time


def main():
    manifest, sequence, report = sys.argv[1:]
    report = Path(report).resolve()
    report.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[2]
    producer_cmd = [sys.executable, str(root/'scripts/alphafold-list-results.py'), manifest, sequence, '--highest']
    consumer_cmd = [sys.executable, str(root/'scripts/alphafold-rosetta-integrity.py'), '-', '--tsv']
    env = dict(os.environ, PYTHONPATH=str(root))
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix='pipeline-work-', dir=report) as temporary:
        env['TMPDIR'] = temporary
        with (report/'producer.stderr').open('w') as errors:
            producer = subprocess.Popen(producer_cmd, stdout=subprocess.PIPE, stderr=errors, env=env)
            try:
                consumer = subprocess.run(consumer_cmd, stdin=producer.stdout, capture_output=True,
                                          text=True, env=env, timeout=600)
            finally:
                producer.stdout.close()
                producer.wait(timeout=30)
        (report/'summary.tsv').write_text(consumer.stdout)
        (report/'consumer.stderr').write_text(consumer.stderr)
        assert producer.returncode == 0, 'producer failed'
        assert consumer.returncode == 0, consumer.stderr
        rows = list(csv.DictReader(io.StringIO(consumer.stdout), delimiter='\t'))
        assert len(rows) >= 2 and rows[0]['scope'] == 'complex', 'missing summaries'
        assert all(row['sequence_id'] == sequence for row in rows), 'wrong sequence'
        assert len({(row['fold_id'], row['model_id'], row['fold_locator']) for row in rows}) == 1
        assert int(rows[0]['nres']) == sum(int(row['nres']) for row in rows[1:])
        assert not list(Path(temporary).iterdir()), 'temporary working files survived success'
        # Real CLI parsing and cleanup errors, without launching model computations.
        cases = [('empty', '', [], 0),
                 ('malformed', 'not\ta\tresult\n', [], 1),
                 ('missing-model', '/missing/fold\n', [], 1),
                 ('missing-fold', '/missing/fold\n', ['--model', '0'], 1)]
        for name, text, options, code in cases:
            result = subprocess.run(consumer_cmd+options, input=text, capture_output=True,
                                    text=True, env=env, timeout=60)
            (report/f'{name}.stderr').write_text(result.stderr)
            assert result.returncode == code, (name, result.stderr)
            if code:
                assert not result.stdout, name
                assert result.stderr.startswith('ERROR:'), name
            else:
                assert len(result.stdout.splitlines()) == 1, name
            assert not list(Path(temporary).iterdir()), f'temporary files survived {name}'
    evidence = dict(status='passed', producer=producer_cmd, consumer=consumer_cmd,
                    elapsed_seconds=round(time.monotonic()-started, 2),
                    fold_id=rows[0]['fold_id'], sequence_id=sequence,
                    model_id=rows[0]['model_id'], summaries=len(rows),
                    checks=['actual subprocess pipe', 'TSV schema/identity', 'residue counts',
                            'temporary cleanup on success and failure', 'empty stdin',
                            'malformed input', 'missing model', 'missing fold'])
    (report/'report.json').write_text(json.dumps(evidence, indent=2)+'\n')
    print(json.dumps(evidence, indent=2))


if __name__ == '__main__':
    main()
