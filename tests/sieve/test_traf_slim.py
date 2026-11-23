import csv
import gzip
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from Bio.PDB import MMCIFIO, MMCIFParser
from Bio.PDB.Chain import Chain
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from sieve.alphafold_interface import COLUMNS, run
from sieve.alphafold_manifest import fold_id, sequence_id
from sieve.alphafold_result_format import PatchInput
from sieve.interface import traf_slim as traf
from sieve.interface.measurements import donor_hydrogens, heavy_contacts
from sieve.patch_description_summary import PatchDescriptions
from sieve.interface.registry import resolve
from sieve.interface.traf_topology import find_groove
from tests.scripts.helpers import load_script
from tests.sieve.test_alphafold_results import zip_bytes


FIXTURES = Path(__file__).resolve().parents[1] / 'fixtures' / 'traf_slim'
CASES = (
    ('1LB5', 'A', 'B', 603, 608, 'traf6'),
    ('1LB6', 'A', 'B', 604, 609, 'traf6'),
    ('1CZZ', 'A', 'D', 250, 254, 'traf2/3'),
    ('1CA9', 'B', 'H', 422, 426, 'traf2/3'),
    ('1L0A', 'A', 'B', 3, 7, 'traf2/3'),
    ('4GHU', 'A', 'B', 143, 147, 'traf2/3'),
)


def load_model(code):
    with gzip.open(FIXTURES / f'{code}.cif.gz', 'rt') as stream:
        return MMCIFParser(QUIET=True).get_structure(code, stream)[0]


def append_copy(chain, residues, first, offset):
    for n, original in enumerate(residues, first):
        clone = original.copy()
        clone.id = (' ', n, ' ')
        for atom in clone:
            atom.coord += offset
        chain.add(clone)


class TrafTests(TestCase):
    def test_confidence_threshold_boundaries(self):
        for score, label in (
            (None, None), (float('nan'), None), (float('inf'), None),
            (0, 'Reject'), (0.229999, 'Reject'), (0.23, 'Low'),
            (0.499999, 'Low'), (0.50, 'Medium'), (0.699999, 'Medium'), (0.70, 'High'),
        ):
            with self.subTest(score=score):
                self.assertEqual(traf.confidence({'pDockQ2 max': score}), label)

    def test_cli_confidence_uses_higher_direction(self):
        model = load_model('1LB5')
        for atom in model.get_atoms():
            atom.set_bfactor(95)
        rows = self.run_cli(model, '--motif-window', 'B:603-608', reverse_pae=30)
        forward, reverse = (float(row['pDockQ2 to partner']) for row in rows)
        self.assertGreaterEqual(forward, 0.70)
        self.assertLess(reverse, 0.23)
        self.assertEqual([row['confidence'] for row in rows], ['High', 'High'])
        self.assertEqual(float(rows[0]['pDockQ2 max']), forward)

    def test_pose_labels_map_to_shared_classification(self):
        for pose, label, face in (
            ('traf6', 'Canonical', 'traf6'),
            ('traf2/3', 'Canonical', 'traf2/3'),
            ('atypical', 'Atypical', None),
            ('unresolved', None, None),
            ('no_contact', None, None),
        ):
            with self.subTest(pose=pose):
                result = traf.classify_interface({'traf_pose': pose}, None, {}, chains=('T', 'P'))
                self.assertEqual(result.interface, label)
                self.assertEqual(result.sides['T'].get('candidate face'), face)
                self.assertIsNone(result.sides['P'].get('candidate face'))

    def test_registry_orients_traf5(self):
        profile, chains = resolve(('P', 'T'), {'P': 'SLiM', 'T': 'TRAF5'})
        self.assertEqual((profile.name, chains), ('traf_slim', ('T', 'P')))
        self.assertEqual(resolve(('T', 'P'), {'T': 'TRAF_MATH', 'P': 'SLIM'})[0], profile)
        with self.assertRaisesRegex(ValueError, 'unsupported interface domains'):
            resolve(('T', 'P'), {'T': 'TRAF', 'P': 'SLIM'})

    def test_feasibility_thresholds_and_no_upper_limit(self):
        limits = traf.Thresholds()
        base = {'tight-contact residues 2': 3, 'steric clash residue pairs': 0,
                'nonpolar residue pairs': 0, 'backbone hydrogen bonds': 0,
                'side-chain hydrogen bonds': 0, 'unresolved hydrogen-bond candidates': 0}
        for area, span, changes, expected in (
            (149.99, 2, {'backbone hydrogen bonds': 3}, 'insufficient_support'),
            (150, 2, {'backbone hydrogen bonds': 3}, 'small_but_well_paired'),
            (150, 1, {'backbone hydrogen bonds': 3}, 'insufficient_support'),
            (250, 0, {'side-chain hydrogen bonds': 2, 'nonpolar residue pairs': 2}, 'hydrogen_bond_supported'),
            (400, 0, {'nonpolar residue pairs': 5}, 'packing_supported_provisional'),
            (1800, 0, {'nonpolar residue pairs': 5}, 'packing_supported_provisional'),
            (1800, 0, {'nonpolar residue pairs': 5, 'steric clash residue pairs': 3}, 'clash_or_sparse_contact'),
            (400, 0, {'nonpolar residue pairs': 5, 'tight-contact residues 2': 1}, 'clash_or_sparse_contact'),
            (250, 0, {'nonpolar residue pairs': 2, 'unresolved hydrogen-bond candidates': 1}, 'unresolved'),
        ):
            with self.subTest(area=area, changes=changes):
                metrics = {**base, **changes, 'total buried surface area': area}
                self.assertEqual(traf.feasibility_rule(metrics, span, limits), expected)

    def test_threshold_validation(self):
        self.assertEqual(traf.parse_thresholds(['paired_area=200']).paired_area, 200)
        for args in (['paired_area=NaN'], ['missing=10'], ['mixed_area=1', 'mixed_area=2']):
            with self.subTest(args=args), self.assertRaises(ValueError):
                traf.parse_thresholds(args)

    def test_renumbered_complexes_recover_pose_from_structure(self):
        for code, a, b, start, stop, label in CASES:
            with self.subTest(code=code):
                model = load_model(code)
                for chain in model:
                    original = list(chain)
                    for r in original:
                        chain.detach_child(r.id)
                        r.id = (' ', r.id[1] + 1000, ' ')
                        chain.add(r)
                domain = find_groove(traf.protein_residues(model[a]), donor_hydrogens, terminal=True)
                self.assertIsNotNone(domain)
                motif = [model[b][n] for n in range(start + 1000, stop + 1001)]
                candidates, row = self.evaluate(model, (a, b), domain, motif, supplied=True)
                self.assertTrue(traf.feasibility(candidates[0].metrics))
                traf.update_site(row, candidates, domain)
                self.assertEqual(row['traf_pose'], label)

    def test_topology_terminal_extension_truncation_and_explicit_window(self):
        chain = load_model('1LB5')['A']
        native = traf.protein_residues(chain)
        for length in (0, 8, 20):
            with self.subTest(upstream=length):
                coil = [r.copy() for r in native[:length]]
                for n, r in enumerate(coil, native[0].id[1] - length):
                    r.id = (' ', n, ' ')
                    for a in r:
                        a.coord += 100
                groove = find_groove(coil + native, donor_hydrogens, terminal=True)
                self.assertEqual(traf.residue_span(groove.beta7), '463-479')
                self.assertEqual(traf.residue_span(groove.domain), '351-501')
        append_copy(chain, native[:12], 502, 100)
        self.assertIsNone(find_groove(traf.protein_residues(chain), donor_hydrogens, terminal=True))
        self.assertIsNotNone(find_groove(native, donor_hydrogens, terminal=False))
        self.assertIsNone(find_groove(native[:-16], donor_hydrogens, terminal=True))
        self.assertIsNone(find_groove(native[30:], donor_hydrogens, terminal=True))

    def test_cli_infers_motif_on_long_tail_and_accepts_explicit_window(self):
        model = load_model('1LB5')
        core = [model['B'][n] for n in range(603, 609)]
        append_copy(model['B'], core * 3, 609, 100)
        inferred, supplied = self.run_cli(model), self.run_cli(
            model, '--math-window', 'A:351-501', '--motif-window', 'B:603-608')
        self.assertEqual([row['chain'] for row in inferred], ['A', 'B'])
        self.assertEqual(inferred[0]['traf_motif_source'], 'inferred')
        self.assertEqual(supplied[0]['traf_motif_source'], 'supplied')
        self.assertEqual(inferred[0]['traf_motif_register'], supplied[0]['traf_motif_register'])
        self.assertEqual([row['interface'] for row in inferred + supplied], ['Canonical'] * 4)
        self.assertEqual([row['candidate face'] for row in inferred + supplied], ['traf6', ''] * 2)
        forward, reverse = (float(row['pDockQ2 to partner']) for row in inferred)
        self.assertEqual(float(inferred[0]['pDockQ2 max']), max(forward, reverse))
        self.assertEqual({row['contact basis'] for row in inferred}, {'heavy atom'})
        self.assertEqual({row['backbone hydrogen bonds'] for row in supplied}, {'5'})
        self.assertEqual(supplied[0]['hydrogen-bond atom pairs'].count('A:472/'), 2)
        self.assertEqual({row['region mean pLDDT'] == '' for row in supplied}, {False})
        self.assertNotIn('face status', inferred[0])
        self.assertNotIn('traf_nonpolar_pairs', inferred[0])
        self.assertTrue(all(row['patch hydropathy mean'] and row['interface structure']
                            for row in inferred))
        first, second = ([model[c][n] for n in range(lo, hi + 1)]
                         for c, lo, hi in (('A', 351, 501), ('B', 603, 608)))
        contact_distances = {(a.id[1], b.id[1]): float(a['CA'] - b['CA'])
                             for a in first for b in second
                             if any(x - y <= 6 for x in a for y in b
                                    if x.element not in ('H', 'D') and y.element not in ('H', 'D'))}
        self.assertEqual(int(supplied[0]['contact count']), len(contact_distances))
        self.assertAlmostEqual(float(supplied[0]['interchain CA distance max']),
                               max(contact_distances.values()), places=5)
        self.assertAlmostEqual(sum(float(row['buried surface area']) for row in supplied),
                               float(supplied[0]['total buried surface area']), places=5)

    def test_mixed_manifest_export_keeps_one_header_and_two_rows_per_patch(self):
        model = load_model('1LB5')
        names = {'A': 'receptor', 'B': 'ligand', 'C': 'first death', 'D': 'second death'}
        domains = {'A': 'TRAF_MATH', 'B': 'SLIM', 'C': 'DD', 'D': 'DED'}
        for chain_id, offset in (('C', 100), ('D', 104)):
            chain = Chain(chain_id)
            model.add(chain)
            append_copy(chain, list(model['A'])[:4], 1, offset)
        with TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            archive, manifest, domain_file, output = (
                root / name for name in ('mixed.zip', 'manifest.tsv', 'domains.tsv', 'out.tsv')
            )
            writer = MMCIFIO()
            writer.set_structure(model)
            buffer = io.StringIO()
            writer.save(buffer)
            cif = MMCIF2Dict(io.StringIO(buffer.getvalue()))
            ids = {c: str(i) for i, c in enumerate(names, 1)}
            cif['_entity.id'] = list(ids.values())
            cif['_entity.type'] = ['polymer'] * len(ids)
            cif['_struct_asym.id'] = list(ids)
            cif['_struct_asym.entity_id'] = list(ids.values())
            cif['_atom_site.label_entity_id'] = [ids[c] for c in cif['_atom_site.label_asym_id']]
            polymers = [(c.id, i, r.resname) for c in model for i, r in enumerate(c, 1)]
            cif['_entity_poly_seq.entity_id'] = [ids[c] for c, _, _ in polymers]
            cif['_entity_poly_seq.num'] = [str(i) for _, i, _ in polymers]
            cif['_entity_poly_seq.mon_id'] = [residue for _, _, residue in polymers]
            writer.set_dict(cif)
            buffer = io.StringIO()
            writer.save(buffer)
            keys = [(c.id, r.id[1]) for c in model for r in c]
            archive.write_bytes(zip_bytes({
                'fold_model_0.cif': buffer.getvalue(),
                'fold_full_data_0.json': json.dumps(dict(
                    token_chain_ids=[c for c, _ in keys], token_res_ids=[n for _, n in keys],
                    pae=[[1] * len(keys) for _ in keys],
                )),
                'fold_job_request.json': '{"modelSeeds":[42]}',
            }))
            sequences = [''.join(traf.residue_label(r) for r in c) for c in model]
            with manifest.open('w') as stream:
                rows = csv.writer(stream, delimiter='\t')
                rows.writerow(['fold_locator', 'fold_id', 'sequence_id', 'input description', 'input chain'])
                for c, name in names.items():
                    rows.writerow([archive, fold_id(42, sequences), sequence_id(sequences), name, c])
            domain_file.write_text('Sequence Name\tDomain Type\n' + ''.join(
                f'{name}\t{domains[c]}\n' for c, name in names.items()
            ))
            errors = io.StringIO()
            self.assertEqual(run([PatchInput(str(archive))], output, manifest=manifest,
                                 domain_file=domain_file, stderr=errors), 0, errors.getvalue())
            with output.open() as stream:
                reader = csv.DictReader(stream, delimiter='\t')
                rows = list(reader)
                self.assertEqual(reader.fieldnames[:len(COLUMNS)], list(COLUMNS))
            self.assertEqual(len(rows), 4)
            self.assertEqual(len({row['fold_id'] for row in rows}), 1)
            self.assertEqual(len({row['patch_id'] for row in rows}), 2)
            self.assertEqual(reader.fieldnames.count('candidate face'), 1)
            for family in ('death', 'traf_slim'):
                pair = [r for r in rows if r['interface_family'] == family]
                self.assertEqual([r['sequence_index'] for r in pair], ['1', '2'])
                self.assertTrue(all(r['total buried surface area'] for r in pair))
                self.assertTrue(all(r['kink angle'] == ('75' if family == 'death' else '') for r in pair))
                self.assertTrue(all(bool(r['traf_pose_evidence']) == (family == 'traf_slim') for r in pair))
                self.assertTrue(all(r['confidence'] in ('High', 'Medium', 'Low', 'Reject') for r in pair))
                self.assertTrue(all(r['chemical feasibility'] in ('true', 'false', 'NULL') for r in pair))
                if family == 'traf_slim':
                    self.assertEqual([r['interface'] for r in pair], ['Canonical', 'Canonical'])
                    self.assertEqual([r['candidate face'] for r in pair], ['traf6', ''])
            self.assertNotIn('\tNA\t', output.read_text())

    def test_inferred_windows_recover_fixture_pose_without_labels(self):
        for code, a, b, _, _, label in CASES:
            with self.subTest(code=code):
                model = load_model(code)
                groove = find_groove(traf.protein_residues(model[a]), donor_hydrogens, terminal=True)
                motif = traf.protein_residues(model[b])
                candidates, row = self.evaluate(model, (a, b), groove, motif)
                traf.update_site(row, candidates, groove)
                self.assertEqual(row['traf_pose'], label)

    def test_ambiguous_distinct_sites_report_both_windows(self):
        model = load_model('1LB5')
        peptide = model['B']
        original = list(peptide)
        append_copy(peptide, original, 609, 100)
        append_copy(peptide, original, 617, 0)
        result = self.run_cli(model)[0]
        self.assertEqual(result['interface'], '')
        self.assertIn('competing inferred motif sites', result['traf_status'])
        self.assertIn(',', result['traf_motif_candidates'])

    def test_reversed_and_shifted_poses_lack_traf6_signature(self):
        for transform in ('reverse', 'shift'):
            with self.subTest(transform=transform):
                model = load_model('1LB5')
                motif = [model['B'][n] for n in range(603, 609)]
                coordinates = [r.copy() for r in motif]
                shift = coordinates[2]['CA'].coord - coordinates[0]['CA'].coord
                for residue, other in zip(motif, coordinates[::-1]):
                    for atom in residue:
                        if transform == 'shift':
                            atom.coord += shift
                        elif atom.name in other:
                            atom.coord = other[atom.name].coord
                row = self.run_cli(model, '--motif-window', 'B:603-608')[0]
                self.assertNotEqual(row['candidate face'], 'traf6')

    def test_no_contact_and_missing_coverage(self):
        model = load_model('1LB5')
        for atom in model['B'].get_atoms():
            atom.coord += 100
        row = self.run_cli(model)[0]
        self.assertEqual((row['interface'], row['confidence']), ('', 'NULL'))
        self.assertEqual(row['traf_status'], 'no_contact')
        self.assertEqual((row['contact count'], row['tight-contact residues']), ('0', '0'))
        row = self.run_cli(model, '--motif-window', 'B:603-609')[0]
        self.assertEqual((row['interface'], row['traf_status']), ('', 'missing motif coordinates'))

    def test_shared_descriptors_agree_for_identical_full_chain_selections(self):
        model = load_model('1LB5')
        descriptions = PatchDescriptions({}, None)
        descriptions.models[('fixture', 0)] = model, {}, [], {'status': 'unavailable'}
        chains = ('A', 'B')
        selections = {c: {r.id[1] for r in model[c]} for c in chains}
        row = {'source': 'fixture', 'model': 0, 'chain 1': 'A', 'chain 2': 'B',
               **{f'interface residues {side}': ','.join(map(str, sorted(selections[c])))
                  for side, c in enumerate(chains, 1)}}
        full = descriptions.describe(row)
        cropped = descriptions.describe(row, context=selections)
        self.assertEqual(full['metrics'], cropped['metrics'])
        self.assertEqual(full['evidence'], cropped['evidence'])

    def test_explicit_partial_math_window_retains_confidence_and_chemistry(self):
        row = self.run_cli(load_model('1LB5'), '--math-window', 'A:405-479',
                           '--motif-window', 'B:603-608')[0]
        self.assertEqual(row['interface'], '')
        self.assertEqual(row['traf_pose_evidence'], 'groove landmarks unresolved')
        self.assertNotEqual(row['confidence'], 'NULL')
        self.assertNotEqual(row['chemical feasibility'], 'NULL')

    @staticmethod
    def evaluate(model, chains, groove, motif, supplied=False):
        descriptions = PatchDescriptions({}, None)
        descriptions.models[('fixture', 0)] = model, {}, [], {'status': 'unavailable'}
        row = {'source': 'fixture', 'model': 0, 'chain 1': chains[0], 'chain 2': chains[1],
               'traf_math_window': traf.residue_span(groove.domain)}
        contacts = heavy_contacts(groove.domain, motif, 6.0)
        tight = {pair for pair, distance in contacts.items() if distance <= 4.5}
        windows = ([(motif, {b for _, b in tight}, len(tight) / len(motif))] if supplied
                   else traf.candidate_windows(motif, contacts, tight))
        candidates = traf.evaluate_sites(descriptions, row, groove.domain, windows, traf.Thresholds(),
                                         ({c: [] for c in chains}, None), contacts, tight)
        return candidates, row

    def run_cli(self, model, *extra, reverse_pae=4):
        script = load_script(Path(__file__).resolve().parents[2] / 'scripts/alphafold-interface.py')
        with TemporaryDirectory() as temp:
            root = Path(temp)
            cif, data, output = (root / name for name in ('model.cif', 'confidence.json', 'results.tsv'))
            writer = MMCIFIO()
            writer.set_structure(model)
            writer.save(str(cif))
            keys = [(c.id, r.id[1]) for c in model for r in traf.protein_residues(c)]
            data.write_text(json.dumps(dict(token_chain_ids=[c for c, _ in keys],
                                            token_res_ids=[n for _, n in keys],
                                            pae=[[1 if c == d or c == 'A' else reverse_pae for d, _ in keys]
                                                 for c, _ in keys])))
            self.assertEqual(script.main([str(cif), '--full-data', str(data), '--output', str(output),
                                          '--chains', 'B', 'A', '--domain', 'A:TRAF_MATH',
                                          '--domain', 'B:SLiM', *extra]), 0)
            with output.open() as stream:
                return list(csv.DictReader(stream, delimiter='\t'))
