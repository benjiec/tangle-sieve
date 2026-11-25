import csv
import gzip
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

from Bio.PDB import MMCIFIO, MMCIFParser
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB.Chain import Chain
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from sieve.alphafold_interface import COLUMNS, run
from sieve.alphafold_manifest import fold_id, sequence_id
from sieve.alphafold_pdockq2 import parse_model
from sieve.alphafold_result_format import PatchInput
from sieve.alphafold_secondary import load_secondary
from sieve.interface import traf_slim as traf
from sieve.interface.measurements import donor_hydrogens, heavy_contacts, qualify_contacts
from sieve.interface.traf_motif import match_motifs, motif_faces
from sieve.interface_description import DescriptionOptions, atom_interactions
from sieve.patch_description_summary import PatchDescriptions
from sieve.interface.registry import resolve
from sieve.interface.traf_topology import find_dssp_groove, find_groove
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
DSSP_STRANDS = {
    '1LB5': [(351, 357), (373, 376), (380, 381), (388, 395), (406, 414), (429, 434),
             (446, 451), (466, 466), (469, 478), (488, 489), (492, 500)],
    '1LB6': [(351, 357), (373, 376), (380, 381), (388, 395), (406, 414), (429, 434),
             (446, 451), (466, 466), (469, 478), (489, 489), (492, 500)],
    '1CZZ': [(353, 359), (376, 377), (381, 382), (389, 395), (406, 414), (430, 434),
             (443, 447), (463, 463), (467, 474), (486, 486), (489, 496)],
    '1CA9': [(353, 359), (376, 377), (381, 382), (389, 395), (406, 414), (430, 434),
             (443, 447), (448, 448), (463, 463), (467, 474), (486, 486), (489, 496)],
    '1L0A': [(353, 359), (375, 377), (381, 382), (389, 395), (406, 414), (430, 434),
             (444, 448), (464, 464), (468, 475), (486, 486), (489, 496)],
    '4GHU': [(416, 421), (439, 440), (444, 445), (452, 458), (469, 477), (493, 497),
             (507, 511), (527, 527), (531, 538), (549, 549), (552, 552), (553, 559)],
}
RANK_MOTIF_FACES = 'traf6_binding_strict; traf2_binding_variant'


def fixture_dssp(residues, code, chain):
    codes = {(chain, r.id[1]): '-' for r in residues}
    for start, stop in DSSP_STRANDS[code]:
        for number in range(start, stop + 1):
            codes[chain, number] = 'B' if start == stop else 'E'
    return codes


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

    def test_motif_rules_and_precedence(self):
        for sequence, family, name in (
            ('PQEINF', 'traf6', 'strict'),
            ('PQQINH', 'traf6', 'permissive'),
            ('PAQAT', 'traf2/3', 'strict'),
            ('PFEES', 'traf2/3', 'strict'),
            ('PANAA', 'traf2/3', 'variant'),
            ('PAAQFLE', 'traf2/3', 'permissive'),
            ('PAEE', 'traf2/3', 'major_consensus'),
            ('SAQE', 'traf2/3', 'major_consensus'),
            ('AAQE', 'traf2/3', 'major_consensus'),
            ('TAQE', 'traf2/3', 'major_consensus'),
        ):
            with self.subTest(sequence=sequence):
                matches, _ = match_motifs(dict(enumerate(sequence, 1)), 1)
                prefix = 'traf2' if family == 'traf2/3' else family
                match = matches[family]
                self.assertEqual(match.rule.name, f'{prefix}_binding_{name}')
                self.assertEqual((match.start, match.end), (1, len(sequence)))
                self.assertEqual(match.sequence, sequence)

    def test_pose_and_local_motif_jointly_determine_interface(self):
        template = load_model('1LB5')['B'][603]
        names = {aa: name for name, aa in protein_letters_3to1.items()}
        for pose, sequence, label, motif_face in (
            ('traf6', 'PQEINF', 'Canonical', 'traf6_binding_strict'),
            ('traf2/3', 'SAQE', 'Canonical', 'traf2_binding_major_consensus'),
            ('traf2/3', 'PANAA', 'Canonical', 'traf2_binding_variant'),
            ('traf6', 'PQQINH', 'Variant', 'traf6_binding_permissive'),
            ('traf2/3', 'PAAQFLE', 'Variant', 'traf2_binding_permissive'),
            ('traf6', 'PAEEAF', 'Canonical', 'traf6_binding_strict; traf2_binding_major_consensus'),
            ('traf6', 'PAQATG', 'Atypical', 'traf2_binding_strict'),
            ('traf6', 'AAAAAAAPQEINF', 'Atypical', None),
            ('traf6', 'PQEIN', None, None),
            ('atypical', 'PQEINF', 'Atypical', 'traf6_binding_strict'),
            ('unresolved', 'PQEINF', None, 'traf6_binding_strict'),
            ('no_contact', 'PQEINF', None, None),
        ):
            with self.subTest(pose=pose, sequence=sequence):
                peptide = Chain('P')
                append_copy(peptide, [template] * len(sequence), 1, 0)
                for residue, aa in zip(peptide, sequence):
                    residue.resname = names[aa]
                row = {
                    'traf_pose': pose, 'traf_motif_bounds': (1, 6),
                    'traf_motif_anchor': 1 if pose in ('traf6', 'traf2/3') else None,
                    'interface residues 2': '' if pose == 'no_contact' else '2,3,4,5,6',
                }
                result = traf.classify_interface(row, ({'P': peptide},), {}, chains=('T', 'P'))
                self.assertEqual(result.interface, label)
                self.assertEqual(result.sides['T'].get('candidate face'),
                                 pose if pose in ('traf6', 'traf2/3') else None)
                self.assertEqual(result.sides['P'].get('candidate face'), motif_face)
                if motif_face is None:
                    self.assertNotIn('traf_motif_matches', result.settings)
                elif sequence == 'PAEEAF':
                    self.assertEqual(result.settings['traf_motif_matches'],
                                     'traf6_binding_strict:P:1-6:PAEEAF; '
                                     'traf2_binding_major_consensus:P:1-4:PAEE')

    def test_unanchored_motif_requires_a_unique_contacting_register(self):
        sequence = dict(enumerate('PQEINFAPQEINF', 1))
        matches, _ = motif_faces(sequence, {3, 5}, anchor=None, window=(1, 13))
        self.assertEqual(matches['traf6'].rule.name, 'traf6_binding_strict')
        matches, _ = motif_faces(sequence, {3, 5, 10, 12}, anchor=None, window=(1, 13))
        self.assertEqual(matches, {})

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

    def test_dssp_backward_scan_recovers_fixture_poses(self):
        for code, a, b, start, stop, label in CASES:
            with self.subTest(code=code):
                model = load_model(code)
                receptor = traf.protein_residues(model[a])
                groove = find_dssp_groove(receptor, fixture_dssp(receptor, code, a), a, terminal=True)
                self.assertIsNotNone(groove)
                motif = [model[b][n] for n in range(start, stop + 1)]
                candidates, row = self.evaluate(model, (a, b), groove, motif, supplied=True)
                self.assertEqual(row['traf_pose'], label)
                result = traf.classify_interface(row, (model,), {}, chains=(a, b))
                self.assertEqual(result.interface, 'Canonical')
                if code == '1CA9':
                    self.assertEqual(result.sides[b]['candidate face'], 'traf2_binding_major_consensus')

    def test_dssp_scan_can_pass_a_terminal_helix_and_split_strands(self):
        model = load_model('1LB5')
        native = list(model['A'])
        append_copy(model['A'], native[:21], 502, 100)

        def reported_dssp(cif, executable):
            parsed, _, _, _ = load_secondary(cif, None)
            codes = fixture_dssp(list(parsed['A']), '1LB5', 'A')
            codes.update({('A', n): 'H' for n in range(502, 523)})
            return parsed, codes, [], {'status': 'assigned'}

        with patch('sieve.patch_description_summary.load_secondary', side_effect=reported_dssp), patch(
                'sieve.interface.traf_slim.find_groove', side_effect=AssertionError('used coordinate fallback')):
            row = self.run_cli(model, '--motif-window', 'B:603-608')[0]
        self.assertEqual((row['traf_math_window'], row['traf_beta7']), ('351-522', '464-478'))
        self.assertEqual((row['interface'], row['candidate face']), ('Canonical', 'traf6'))

    def test_dssp_scan_rejects_two_backward_assignments(self):
        residues = list(load_model('1LB5')['A'])
        spans = [(i, i + 1) for i in (145, 125, 105, 85, 65, 45, 25, 10, 0)]
        neighbors = {i: set() for i in range(len(spans))}
        for a, b in ((0, 3), (0, 7), (0, 8), (1, 4), (2, 3), (4, 5), (5, 6)):
            neighbors[a].add(b)
            neighbors[b].add(a)
        with patch('sieve.interface.traf_topology.dssp_strands', return_value=spans), patch(
                'sieve.interface.traf_topology.strand_neighbors', return_value=neighbors):
            self.assertIsNone(find_dssp_groove(residues, {}, 'A', terminal=True))

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
        self.assertEqual([row['candidate face'] for row in inferred + supplied],
                         ['traf6', RANK_MOTIF_FACES] * 2)
        self.assertEqual({row['traf_motif_matches'] for row in inferred + supplied},
                         {'traf6_binding_strict:B:603-608:PTEDEY; '
                          'traf2_binding_variant:B:603-607:PTEDE'})
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

    def test_cli_rescores_traf6_chemistry_without_changing_geometry(self):
        model = load_model('1LB5')
        automatic = self.run_cli(model)
        geometry = automatic[1]['region']
        original = self.run_cli(model, '--motif-window', f'B:{geometry}')
        complete = self.run_cli(model, '--motif-window', 'B:603-608')
        self.assertEqual(geometry, '604-608')
        for auto, core, full in zip(automatic, original, complete):
            self.assertEqual(auto['traf_scoring_sequence'], 'TEDEY')
            self.assertEqual(full['traf_scoring_sequence'], 'PTEDEY')
            self.assertEqual((auto['traf_chemistry_window'], core['traf_chemistry_window']),
                             ('603-608', '604-608'))
            for key in ('region', 'contact count', 'pDockQ2 max', 'interface', 'candidate face',
                        'interface structure', 'patch hydropathy mean'):
                self.assertEqual(auto[key], core[key], key)
            for key in ('total buried surface area', 'chemical feasibility', 'nonpolar residue pairs',
                        'steric clash residue pairs', 'tight-contact residues', 'traf_feasibility_rule'):
                self.assertEqual(auto[key], full[key], key)
            self.assertGreater(float(auto['total buried surface area']), float(core['total buried surface area']))

    def test_cli_uses_short_traf2_consensus_for_chemistry(self):
        model = load_model('1CA9')
        rows = self.run_cli(model, chains=('B', 'H'))
        self.assertEqual(rows[1]['candidate face'], 'traf2_binding_major_consensus')
        self.assertEqual(rows[1]['interface'], 'Canonical')
        self.assertEqual(rows[1]['traf_chemistry_window'], '424-427')
        self.assertNotEqual(rows[1]['region'], '424-427')

    def test_cli_permissive_motif_rescores_and_distant_motif_does_not(self):
        for name, residue, label, bounds in (
            ('permissive', (608, 'HIS'), 'Variant', '603-608'),
            ('distant', (603, 'ALA'), 'Atypical', '604-608'),
        ):
            with self.subTest(name=name):
                model = load_model('1LB5')
                chain = model['B']
                append_copy(chain, [chain[n] for n in range(603, 609)], 609, 100)
                number, amino_acid = residue
                chain[number].resname = amino_acid
                rows = self.run_cli(model)
                self.assertEqual(rows[1]['region'], '604-608')
                self.assertEqual(rows[1]['interface'], label)
                self.assertEqual(rows[1]['traf_chemistry_window'], bounds)
                if name == 'distant':
                    self.assertEqual(rows[1]['traf_scoring_sequence'], 'TEDEY')
                    self.assertEqual(rows[1]['traf_motif_matches'], '')

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
                    self.assertEqual([r['candidate face'] for r in pair], ['traf6', RANK_MOTIF_FACES])
            self.assertNotIn('\tNA\t', output.read_text())

    def test_inferred_windows_recover_fixture_pose_without_labels(self):
        for code, a, b, _, _, label in CASES:
            with self.subTest(code=code):
                model = load_model(code)
                groove = find_groove(traf.protein_residues(model[a]), donor_hydrogens, terminal=True)
                motif = traf.protein_residues(model[b])
                candidates, row = self.evaluate(model, (a, b), groove, motif)
                if code == '1L0A':
                    self.assertEqual(row['traf_pose'], 'unresolved')
                    self.assertEqual(row['traf_pose_evidence'],
                                     'competing pocket registers pocket/loop=1:5,3:5')
                else:
                    self.assertEqual(row['traf_pose'], label)

    def test_full_patch_recovers_packing_anchor_outside_scoring_window(self):
        model = load_model('1CZZ')

        def clipped_core(peptide, contacts, tight):
            window = [r for r in peptide if 251 <= r.id[1] <= 255]
            pairs = {(a, b) for a, b in tight if 251 <= b <= 255}
            return [(window, {b for _, b in pairs}, len(pairs) / len(window))]

        with patch('sieve.interface.traf_slim.candidate_windows', side_effect=clipped_core):
            _, full = self.run_cli(model, chains=('A', 'D'))
        _, explicit = self.run_cli(model, '--motif-window', 'D:251-255', chains=('A', 'D'))
        self.assertEqual(full['region'], explicit['region'])
        self.assertEqual(full['pDockQ2 max'], explicit['pDockQ2 max'])
        self.assertEqual(full['traf_pose_window'], '250-256')
        self.assertEqual(full['traf_pose_evidence'], 'pocket/loop=250:252')
        self.assertEqual(full['interface'], 'Canonical')
        self.assertEqual(explicit['traf_pose_window'], '251-255')
        self.assertEqual(explicit['interface'], 'Atypical')

    def test_pose_cannot_combine_anchors_across_contact_runs(self):
        model = load_model('1CZZ')
        chains = ('A', 'D')
        groove = find_groove(traf.protein_residues(model['A']), donor_hydrogens, terminal=True)
        motif = [model['D'][250], model['D'][252]]
        contacts = heavy_contacts(groove.domain, motif, 6)
        selections = {chain: {pair[i] for pair in contacts} for i, chain in enumerate(chains)}
        options = DescriptionOptions()
        _, evidence = qualify_contacts(model, chains, atom_interactions(model, chains, options, selections=selections))
        self.assertEqual(traf.assess_pose(groove, motif, evidence)[0], 'traf2/3')
        row = {'chain 1': 'A', 'chain 2': 'D', 'traf_pose': 'unresolved'}
        traf.assign_patch_pose(row, model, groove, motif, contacts, options)
        self.assertEqual(row['traf_pose_window'], '250,252')
        self.assertEqual(row['traf_pose'], 'atypical')

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
        self.assertIn('competing pose registers', result['traf_pose_evidence'])
        self.assertEqual(result['traf_motif_matches'], '')

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
        self.assertEqual(row['traf_scoring_sequence'], '')
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

    def test_unresolved_topology_keeps_metrics_on_full_chain(self):
        def unmatched_dssp(cif, executable):
            model, _, _, _ = load_secondary(cif, None)
            codes = {('A', r.id[1]): 'H' for r in model['A']}
            return model, codes, [], {'status': 'assigned'}

        with patch('sieve.patch_description_summary.load_secondary', side_effect=unmatched_dssp), patch(
                'sieve.interface.traf_slim.find_groove', side_effect=AssertionError('used coordinate fallback')):
            model = load_model('1LB5')
            rows = self.run_cli(model, '--motif-window', 'B:603-608')
            for row in rows:
                self.assertEqual((row['interface'], row['traf_math_window']), ('', ''))
                self.assertIn('full-chain scoring fallback', row['traf_status'])
                self.assertIn('groove landmarks unresolved', row['traf_pose_evidence'])
                self.assertTrue(all(row[field] for field in (
                    'normalized PAE to partner', 'pDockQ2 to partner', 'interface pLDDT',
                    'total buried surface area',
                )))
                self.assertEqual(row['chemical feasibility'], 'true')
                self.assertNotEqual(row['confidence'], 'NULL')
            self.assertEqual(rows[0]['region'], '347-501')
            self.assertEqual([row['candidate face'] for row in rows], ['', RANK_MOTIF_FACES])

            row = self.run_cli(
                model, '--motif-window', 'B:603-608',
                '--traf-threshold', 'paired_area=100000',
                '--traf-threshold', 'mixed_area=100000',
                '--traf-threshold', 'packing_area=100000',
            )[0]
            self.assertEqual(row['chemical feasibility'], 'false')
            self.assertEqual(row['traf_feasibility_rule'], 'insufficient_support')

            incomplete = load_model('1LB5')
            incomplete['A'][472].detach_child('N')
            row = self.run_cli(incomplete, '--motif-window', 'B:603-608')[0]
            self.assertIn('full-chain scoring fallback', row['traf_status'])
            self.assertEqual(row['chemical feasibility'], 'NULL')
            self.assertEqual(row['traf_feasibility_rule'], 'unresolved')
            self.assertNotEqual(row['confidence'], 'NULL')

            with patch('sieve.alphafold_interface.parse_model',
                       side_effect=lambda cif, data: (parse_model(cif, data)[0], None)):
                row = self.run_cli(model)[0]
            self.assertEqual((row['normalized PAE to partner'], row['pDockQ2 to partner']), ('', ''))
            self.assertTrue(row['interface pLDDT'])
            self.assertEqual(row['chemical feasibility'], 'true')
            for atom in model['B'].get_atoms():
                atom.coord += 100
            row = self.run_cli(model)[0]
            self.assertIn('full-chain scoring fallback; no_contact', row['traf_status'])
            self.assertEqual((row['contact count'], row['pDockQ2 to partner']), ('0', ''))

    @staticmethod
    def evaluate(model, chains, groove, motif, supplied=False):
        descriptions = PatchDescriptions({}, None)
        descriptions.models[('fixture', 0)] = model, {}, [], {'status': 'unavailable'}
        row = {'source': 'fixture', 'model': 0, 'chain 1': chains[0], 'chain 2': chains[1],
               'traf_math_window': traf.residue_span(groove.domain), 'traf_pose': 'unresolved'}
        contacts = heavy_contacts(groove.domain, motif, 6.0)
        tight = {pair for pair, distance in contacts.items() if distance <= 4.5}
        windows = ([(motif, {b for _, b in tight}, len(tight) / len(motif))] if supplied
                   else traf.candidate_windows(motif, contacts, tight))
        candidates = traf.evaluate_sites(descriptions, row, groove.domain, windows, traf.Thresholds(),
                                         ({c: [] for c in chains}, None), contacts, tight)
        traf.update_site(row, candidates)
        traf.assign_patch_pose(row, model, groove, motif, contacts, descriptions.options)
        return candidates, row

    def run_cli(self, model, *extra, reverse_pae=4, chains=('A', 'B')):
        script = load_script(Path(__file__).resolve().parents[2] / 'scripts/alphafold-interface.py')
        with TemporaryDirectory() as temp:
            root = Path(temp)
            cif, data, output = (root / name for name in ('model.cif', 'confidence.json', 'results.tsv'))
            writer = MMCIFIO()
            writer.set_structure(model)
            writer.save(str(cif))
            keys = [(c.id, r.id[1]) for c in model for r in traf.protein_residues(c)]
            receptor, peptide = chains
            data.write_text(json.dumps(dict(token_chain_ids=[c for c, _ in keys],
                                            token_res_ids=[n for _, n in keys],
                                            pae=[[1 if c == d or c == receptor else reverse_pae for d, _ in keys]
                                                 for c, _ in keys])))
            self.assertEqual(script.main([str(cif), '--full-data', str(data), '--output', str(output),
                                          '--chains', peptide, receptor, '--domain', f'{receptor}:TRAF_MATH',
                                          '--domain', f'{peptide}:SLiM', *extra]), 0)
            with output.open() as stream:
                return list(csv.DictReader(stream, delimiter='\t'))
