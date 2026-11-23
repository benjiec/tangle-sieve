import csv
from dataclasses import replace
from functools import partial
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from sieve.alphafold_interface import COLUMNS, run
from sieve.alphafold_manifest import fold_id, run as build_manifest, sequence_id
from sieve.alphafold_pdockq2_patches import PATCH_COLUMNS, Options, analyze, input_models, patch_export_rows
from sieve.alphafold_patch_metrics import extract_patch_metrics
from sieve.alphafold_result_format import PatchInput, ResultRecord, write_result
from sieve.alphafold_secondary import load_secondary
from sieve.interface import death
from sieve.interface.registry import INTERFACES, resolve
from tests.scripts.helpers import load_script
from tests.sieve.test_alphafold_dd_faces import structure
from tests.sieve.test_alphafold_overlap import fixture
from tests.sieve.test_alphafold_results import zip_bytes

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = load_script(ROOT / 'scripts/alphafold-interface.py')


class DeathTests(unittest.TestCase):
    def test_confidence_boundaries(self):
        for score, area, expected in (
            (0.229, 700, 'Reject'), (0.23, 700, 'Low'), (0.399, 700, 'Low'),
            (0.40, 999, 'Medium'), (0.40, 1000, 'Low'), (0.499, 1000, 'Low'),
            (0.50, 1000, 'Medium'), (0.699, 1000, 'Medium'), (0.70, 1000, 'High'),
            (None, 700, None), (0.70, None, None),
        ):
            with self.subTest(score=score, area=area):
                self.assertEqual(death.confidence({
                    'pDockQ2 max': score, 'total buried surface area': area,
                }), expected)

    def test_feasibility_predicates_and_missing_metrics(self):
        zipper = {
            'total buried surface area': 700, 'total C/S contribution percent': 25,
            'nonpolar contacts': 2, 'salt bridge candidates': 2,
            'salt bridge candidate heavy-atom distance min': 3,
            'nominal patch charge difference': -3, 'steric overlap candidates': 1,
            'contact distance min': 4.5,
        }
        self.assertTrue(death.feasibility(zipper))
        for key in zipper:
            if key == 'salt bridge candidate heavy-atom distance min':
                continue
            with self.subTest(missing=key):
                self.assertIsNone(death.feasibility({**zipper, key: None}))
        for key, value in (
            ('total buried surface area', 699), ('total buried surface area', 1200),
            ('total C/S contribution percent', 24.9), ('nonpolar contacts', 1),
            ('salt bridge candidates', 1), ('salt bridge candidate heavy-atom distance min', 3.01),
            ('nominal patch charge difference', -2.99), ('steric overlap candidates', 2),
        ):
            with self.subTest(key=key, value=value):
                self.assertFalse(death.feasibility({**zipper, key: value}))
        self.assertTrue(death.feasibility({**zipper, 'nominal patch charge difference': 3}))
        self.assertIsNone(death.feasibility({**zipper, 'total buried surface area': 0}))
        self.assertIsNone(death.feasibility({**zipper, 'nominal patch charge difference': '+2 / +3'}))

        for area, cs, contacts, salts, distance, overlaps, ca, expected in (
            (499, 40, 2, 1, 4, 1, 4.5, False),
            (500, 40, 2, 1, 4, 1, 4.5, True),  # lateral
            (600, 40, 5, 0, 5, 1, 4.5, True),  # lateral, density branch
            (800, 40, 4, 1, 4, 1, 4.5, True),  # cooperative
            (800, 40, 3, 1, 4, 1, 4.5, False),
            (800, 40, 7, 0, 5, 1, 4.5, True),  # cooperative, density branch
            (1000, 40, 4, 0, 4.01, 2, 4.5, True),  # hydrophobic
            (1000, 40, 4, 0, 4.01, 2, 4.51, False),  # hovering
            (1000, 40, 4, 0, 4, 2, 4.51, True),
            (1000, 40, 5, 0, 4, 3, 4.5, False),  # overlaps
            (1799, 40, 5, 0, 4, 2, 4.5, True),
            (1800, 40, 5, 0, 4, 2, 4.5, False),
            (1800, 45, 10, 1, 3.5, 2, 4.5, True),  # clamp
            (1800, 45, 9, 1, 3.5, 2, 4.5, False),
            (1800, 45, 11, 0, 5, 2, 4.5, True),  # clamp, density branch
        ):
            with self.subTest(area=area, cs=cs, contacts=contacts, distance=distance, ca=ca):
                self.assertIs(death.feasibility({
                    **zipper, 'total buried surface area': area, 'total C/S contribution percent': cs,
                    'nonpolar contacts': contacts, 'salt bridge candidates': salts,
                    'salt bridge candidate heavy-atom distance min': distance,
                    'steric overlap candidates': overlaps, 'contact distance min': ca,
                }), expected)

    def test_missing_salt_distance_uses_fallback(self):
        patch_row = {
            'total buried surface area': 598.121699, 'total C/S contribution percent': 41.306076,
            'nonpolar contacts': 1, 'salt bridge candidates': 0,
            'salt bridge candidate heavy-atom distance min': None,
            'nominal patch charge difference': 2, 'steric overlap candidates': 0,
            'contact distance min': 5.185086,
        }
        self.assertIs(death.feasibility(patch_row), False)
        self.assertIs(death.feasibility({**patch_row, 'nonpolar contacts': 5}), True)

    def test_classification_uses_rule_evidence_and_normalization_defaults(self):
        codes, segments = structure()
        both_codes = {**codes, **{('B', n): code for (_, n), code in codes.items()}}
        both_segments = [*segments, *[dict(segment, chain='B') for segment in segments]]
        patch_row = {
            'chain 1': 'A', 'chain 2': 'B',
            'interface residues 1': '1,2,3,19,20,21', 'interface residues 2': '7,8,9,13,14,15',
        }
        with patch('sieve.interface.death.normalize_helices', return_value=(both_segments, '')) as normalize:
            result = death.classify_interface(
                patch_row, (None, both_codes, both_segments, {'status': 'assigned'}),
                {'A': 'DD', 'B': 'DD'}, chains=('B', 'A'),
            )
        self.assertEqual(result.interface, 'Canonical')
        self.assertTrue(result.sides['A']['candidate face'].startswith('Ia'))
        self.assertTrue(result.sides['B']['candidate face'].startswith('Ib'))
        self.assertTrue(result.sides['A']['evidence'])
        self.assertEqual([call.kwargs for call in normalize.call_args_list],
                         [{'connector_length': 1, 'angle': 75}] * 2)

    def test_registry_orientation_ambiguity_and_unsupported_domains(self):
        death_profile = next(profile for profile in INTERFACES if profile.name == 'death')
        asymmetric = replace(death_profile, name='test', domains_a=frozenset({'TRAF'}),
                             domains_b=frozenset({'SLiM'}),
                             extract=partial(death.extract, distance_cutoff=6.0))
        self.assertEqual(resolve(('A', 'B'), {'A': 'SLiM', 'B': 'TRAF'}, (death_profile, asymmetric)),
                         (asymmetric, ('B', 'A')))
        domains = {'A': 'CARD', 'B': 'DED'}
        self.assertEqual(resolve(('A', 'B'), domains), (death_profile, ('A', 'B')))
        with self.assertRaisesRegex(ValueError, 'ambiguous'):
            resolve(('A', 'B'), domains, (death_profile, replace(death_profile, name='other')))
        with self.assertRaisesRegex(ValueError, 'unsupported'):
            resolve(('A', 'B'), {'A': 'CARD', 'B': 'SLiM'})


class PipelineTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name).resolve()
        self.archive = self.root / 'predictions.zip'
        self.locator = f'{self.archive}::inner.zip!/.'
        self.manifest = self.root / 'manifest.tsv'
        self.domains = self.root / 'domains.tsv'
        self.output = self.root / 'evaluated.tsv'
        self.errors = io.StringIO()
        cif = fixture([('A', 'AAAA', [0, 0, 0]), ('B', 'CCCC', [6, 0, 0]),
                       ('C', 'GGGG', [70, 0, 0])])
        data = json.dumps(dict(token_chain_ids=['A']*4 + ['B']*4 + ['C']*4,
                               token_res_ids=[1, 2, 3, 4]*3, pae=[[1]*12 for _ in range(12)]))
        members = {'fold_job_request.json': '{"modelSeeds":[42]}'}
        for number in (0, 1):
            members[f'fold_model_{number}.cif'] = cif
            members[f'fold_full_data_{number}.json'] = data if number == 1 else 'invalid JSON'
            members[f'fold_summary_confidences_{number}.json'] = '{"iptm":0.6,"ptm":0.7}'
        self.archive.write_bytes(zip_bytes({'inner.zip': zip_bytes(members)}))
        fasta = self.root / 'sequences.faa'
        fasta.write_text('>first sequence\nAAAA\n>second sequence\nCCCC\n>third sequence\nGGGG\n')
        self.assertEqual(build_manifest(self.archive, fasta, self.manifest, io.StringIO()), 0)
        self.domains.write_text('Sequence Name\tDomain Type\nfirst sequence\tDD\nsecond sequence\tDED\n')
        self.dssp = patch('sieve.alphafold_interface.find_dssp', return_value=None)
        self.dssp.start()
        self.addCleanup(self.dssp.stop)

    def read(self):
        with self.output.open() as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            self.assertEqual(reader.fieldnames[:len(COLUMNS)], list(COLUMNS))
            self.assertNotIn('traf_pose_evidence', reader.fieldnames)
            return list(reader)

    def test_piped_model_and_domain_lookup_keep_legacy_columns(self):
        stdin = io.StringIO()
        write_result(csv.writer(stdin, delimiter='\t'), ResultRecord(self.locator, '12345678', 1, None, None))
        stdin.seek(0)
        with patch('sys.stdin', stdin):
            self.assertEqual(SCRIPT.main([
                '-', '--chains', 'B', 'A', '--manifest', str(self.manifest),
                '--domain-file', str(self.domains), '--output', str(self.output),
            ]), 0)
        rows = self.read()
        self.assertEqual([row['chain'] for row in rows], ['A', 'B'])
        self.assertEqual([row['sequence name'] for row in rows], ['first sequence', 'second sequence'])
        self.assertEqual([row['domain'] for row in rows], ['DD', 'DED'])
        self.assertEqual({row['model'] for row in rows}, {'1'})
        self.assertEqual({row['sequence_id'] for row in rows}, {sequence_id(['AAAA', 'CCCC', 'GGGG'])})
        self.assertEqual({row['fold_id'] for row in rows}, {fold_id(42, ['AAAA', 'CCCC', 'GGGG'])})
        self.assertEqual({row['iptm'] for row in rows}, {'0.600000'})
        self.assertEqual({row['interface_family'] for row in rows}, {'death'})
        self.assertEqual({row['interface'] for row in rows}, {'NULL'})
        self.assertEqual({row['salt bridge candidate heavy-atom distance min'] for row in rows}, {''})
        self.assertEqual({row['chemical feasibility'] for row in rows}, {'false'})
        self.assertEqual({row['kink angle'] for row in rows}, {'75'})
        self.assertEqual({row['distance cutoff'] for row in rows}, {'10.000000'})

        args = dict(options=Options(cutoff=10.0), chains=('A', 'B'), selected_model=1)
        expected = extract_patch_metrics(self.locator, **args)
        self.assertEqual(
            [{key: row[key] for key in PATCH_COLUMNS} for row in rows],
            [{key: str(value) for key, value in row.items()} for row in patch_export_rows(expected.leaves(), missing_value='')],
        )
        self.assertEqual(rows[0]['patch_id'], next(
            row['patch_id'] for row in analyze(self.locator, **args) if row['status'] == 'leaf'
        ))

    def test_explicit_and_supplied_names_with_batch_errors(self):
        request = PatchInput(self.locator, 1, sequence_names={'B': 'second sequence'})
        self.assertEqual(run([PatchInput(self.locator, 0), request], self.output,
                             manifest=self.manifest, domain_file=self.domains,
                             domains=['A:CARD'], chains=('A', 'B'),
                             stderr=self.errors), 1)
        self.assertIn('model 0', self.errors.getvalue())
        self.assertEqual([row['domain'] for row in self.read()], ['CARD', 'DED'])

    def test_incomplete_mapping_fails_before_opening_inputs(self):
        self.output.write_text('original')
        for flags in (
            {}, {'chains': ('A', 'B')}, {'domains': ['A:DD', 'B:DED']},
            {'chains': ('A', 'B'), 'domains': ['A:DD', 'C:CARD']},
            {'manifest': self.manifest}, {'domain_file': self.domains},
            {'chains': ('A', 'B'), 'domains': ['A:DD'], 'domain_file': self.domains},
        ):
            with self.subTest(flags=flags), patch('sieve.alphafold_interface.jobs_for') as jobs:
                with self.assertRaisesRegex(ValueError, 'domain mapping requires'):
                    run([PatchInput(self.locator, 1)], self.output, **flags)
                jobs.assert_not_called()
        self.assertEqual(self.output.read_text(), 'original')

    def test_lookup_error_identifies_fold_and_unmapped_chain(self):
        self.output.write_text('original')
        self.assertEqual(run([PatchInput(self.locator, 1)], self.output,
                             manifest=self.manifest, domain_file=self.domains,
                             chains=('A', 'C'), stderr=self.errors), 1)
        self.assertIn(self.locator, self.errors.getvalue())
        self.assertIn("missing domain assignment for chain C for sequence 'third sequence'",
                      self.errors.getvalue())
        self.assertEqual(self.output.read_text(), 'original')

    def test_registry_profile_controls_cutoff(self):
        death_profile = next(profile for profile in INTERFACES if profile.name == 'death')
        custom = replace(death_profile, name='narrow', extract=partial(death.extract, distance_cutoff=0.1))
        with patch('sieve.alphafold_interface.resolve', return_value=(custom, ('A', 'B'))):
            self.assertEqual(run([PatchInput(self.locator, 1)], self.output,
                                 domains=['A:DD', 'B:DED'], chains=('A', 'B')), 0)
        self.assertEqual(self.read(), [])

    def test_manifest_identity_mismatch_preserves_output(self):
        self.output.write_text('original')
        text = self.manifest.read_text()
        expected_id = fold_id(42, ['AAAA', 'CCCC', 'GGGG'])
        self.manifest.write_text(text.replace(expected_id, '00000000'))
        self.assertEqual(run([PatchInput(self.locator, 1)], self.output,
                             manifest=self.manifest, domain_file=self.domains,
                             chains=('A', 'B'), stderr=self.errors), 1)
        self.assertIn('model/manifest identity mismatch', self.errors.getvalue())
        self.assertEqual(self.output.read_text(), 'original')

    def test_tall_rows_keep_side_specific_face_results_and_patch_identity(self):
        request = PatchInput(self.locator, 1)
        flags = dict(domains=['A:DD', 'B:DED'], chains=('B', 'A'))
        self.assertEqual(run([request], self.output, **flags), 0)
        before = self.read()
        with patch('sieve.patch_description_summary.load_secondary') as load, patch(
                'sieve.interface.death.normalize_helices', return_value=([], 'merged example')), patch(
                'sieve.interface.death.classify',
                side_effect=lambda chain, *_: ('Ia' if chain == 'A' else 'Ib', [], 'eligible')):
            record = next(input_models(self.locator, selected_model=1))
            model, _, _, _ = load_secondary(record[2], None)
            load.return_value = (model, {}, [], {'status': 'assigned'})
            self.assertEqual(run([request], self.output, **flags), 0)
        after = self.read()
        self.assertEqual([row['patch_id'] for row in after], [row['patch_id'] for row in before])
        self.assertEqual([(row['chain'], row['candidate face']) for row in after],
                         [('A', 'Ia'), ('B', 'Ib')])
        self.assertEqual({row['face status'] for row in after}, {'eligible; merged example'})
        self.assertEqual({row['interface'] for row in after}, {'Canonical'})

    def test_unmatched_pair_and_failed_atomic_write_preserve_output(self):
        self.output.write_text('original')
        request = PatchInput(self.locator, 1)
        self.assertEqual(run([request], self.output, chains=('A', 'B'),
                             domains=['A:CARD', 'B:unknown'], stderr=self.errors), 1)
        self.assertIn('unsupported interface domains', self.errors.getvalue())
        self.assertEqual(self.output.read_text(), 'original')
        with patch('sieve.alphafold_cache.os.replace', side_effect=OSError('disk failure')):
            with self.assertRaisesRegex(OSError, 'disk failure'):
                run([request], self.output, chains=('A', 'B'), domains=['A:CARD', 'B:DED'])
        self.assertEqual(self.output.read_text(), 'original')
