import csv
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from Bio.PDB import MMCIFIO, MMCIFParser

from sieve.interface import traf_slim as traf
from sieve.interface.registry import resolve
from sieve.interface.traf_topology import find_groove
from tests.scripts.helpers import load_script


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
    def test_registry_orients_traf5(self):
        profile, chains = resolve(('P', 'T'), {'P': 'SLiM', 'T': 'TRAF5'})
        self.assertEqual((profile.name, chains), ('traf_slim', ('T', 'P')))

    def test_feasibility_thresholds_and_no_upper_limit(self):
        limits = traf.Thresholds()
        base = dict(traf_tight_motif_residues=3, traf_clash_pairs=0, traf_nonpolar_pairs=0,
                    traf_backbone_hbonds=0, traf_sidechain_hbonds=0,
                    traf_unresolved_hbond_pairs=0)
        for area, span, changes, expected in (
            (149.99, 2, {'traf_backbone_hbonds': 3}, 'insufficient_support'),
            (150, 2, {'traf_backbone_hbonds': 3}, 'small_but_well_paired'),
            (150, 1, {'traf_backbone_hbonds': 3}, 'insufficient_support'),
            (250, 0, {'traf_sidechain_hbonds': 2, 'traf_nonpolar_pairs': 2}, 'hydrogen_bond_supported'),
            (400, 0, {'traf_nonpolar_pairs': 5}, 'packing_supported_provisional'),
            (1800, 0, {'traf_nonpolar_pairs': 5}, 'packing_supported_provisional'),
            (1800, 0, {'traf_nonpolar_pairs': 5, 'traf_clash_pairs': 3}, 'clash_or_sparse_contact'),
            (400, 0, {'traf_nonpolar_pairs': 5, 'traf_tight_motif_residues': 1}, 'clash_or_sparse_contact'),
            (250, 0, {'traf_nonpolar_pairs': 2, 'traf_unresolved_hbond_pairs': 1}, 'unresolved'),
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
                domain = find_groove(traf.protein_residues(model[a]), traf.donor_hydrogens, terminal=True)
                self.assertIsNotNone(domain)
                motif = [model[b][n] for n in range(start + 1000, stop + 1001)]
                tight = {y for _, y in traf.contact_pairs(domain.domain, motif, 4.5)}
                metrics, interactions = traf.chemical_metrics(model, (a, b), domain.domain,
                                                               motif, tight, traf.Thresholds())
                self.assertTrue(traf.feasibility(metrics))
                self.assertEqual(traf.assess_pose(domain, motif, interactions)[0], label)

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
                groove = find_groove(coil + native, traf.donor_hydrogens, terminal=True)
                self.assertEqual(traf.residue_span(groove.beta7), '463-479')
                self.assertEqual(traf.residue_span(groove.domain), '351-501')
        append_copy(chain, native[:12], 502, 100)
        self.assertIsNone(find_groove(traf.protein_residues(chain), traf.donor_hydrogens, terminal=True))
        self.assertIsNotNone(find_groove(native, traf.donor_hydrogens, terminal=False))
        self.assertIsNone(find_groove(native[:-16], traf.donor_hydrogens, terminal=True))
        self.assertIsNone(find_groove(native[30:], traf.donor_hydrogens, terminal=True))

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
        self.assertEqual([row['interface'] for row in inferred + supplied], ['traf6'] * 4)
        forward, reverse = (float(row['pDockQ2 to partner']) for row in inferred)
        self.assertEqual(float(inferred[0]['confidence']), min(forward, reverse))

    def test_inferred_windows_recover_fixture_pose_without_labels(self):
        for code, a, b, _, _, label in CASES:
            with self.subTest(code=code):
                model = load_model(code)
                groove = find_groove(traf.protein_residues(model[a]), traf.donor_hydrogens, terminal=True)
                motif = traf.protein_residues(model[b])
                contacts = traf.contact_pairs(groove.domain, motif, 6.0)
                tight = traf.contact_pairs(groove.domain, motif, 4.5)
                candidates = traf.evaluate_sites(model, (a, b), groove.domain,
                                                  traf.candidate_windows(motif, contacts, tight),
                                                  traf.Thresholds(), ({a: [], b: []}, None))
                row = {'traf_math_window': traf.residue_span(groove.domain)}
                traf.update_site(row, candidates, groove, contacts)
                self.assertEqual(row['traf_pose'], label)

    def test_ambiguous_distinct_sites_report_both_windows(self):
        model = load_model('1LB5')
        peptide = model['B']
        original = list(peptide)
        append_copy(peptide, original, 609, 100)
        append_copy(peptide, original, 617, 0)
        result = self.run_cli(model)[0]
        self.assertEqual(result['interface'], 'unresolved')
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
                self.assertNotEqual(row['interface'], 'traf6')

    def test_no_contact_and_missing_coverage(self):
        model = load_model('1LB5')
        for atom in model['B'].get_atoms():
            atom.coord += 100
        row = self.run_cli(model)[0]
        self.assertEqual((row['interface'], row['confidence']), ('no_contact', 'NULL'))
        row = self.run_cli(model, '--motif-window', 'B:603-609')[0]
        self.assertEqual((row['interface'], row['traf_status']), ('unresolved', 'missing motif coordinates'))

    def test_explicit_partial_math_window_retains_confidence_and_chemistry(self):
        row = self.run_cli(load_model('1LB5'), '--math-window', 'A:405-479',
                           '--motif-window', 'B:603-608')[0]
        self.assertEqual(row['interface'], 'unresolved')
        self.assertEqual(row['traf_pose_evidence'], 'groove landmarks unresolved')
        self.assertNotEqual(row['confidence'], 'NULL')
        self.assertNotEqual(row['chemical feasibility'], 'NULL')

    def run_cli(self, model, *extra):
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
                                            pae=[[1 if c == d or c == 'A' else 4 for d, _ in keys]
                                                 for c, _ in keys])))
            self.assertEqual(script.main([str(cif), '--full-data', str(data), '--output', str(output),
                                          '--chains', 'B', 'A', '--domain', 'A:TRAF6',
                                          '--domain', 'B:SLiM', *extra]), 0)
            with output.open() as stream:
                return list(csv.DictReader(stream, delimiter='\t'))
