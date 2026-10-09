import io
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

from sieve.alphafold_pdockq2 import Residue, score_pair
from sieve.alphafold_pdockq2_patches import (
    Options, analyze, components, detect_patches, find_split, input_models, write_patches,
)


def graph_from_edges(edges):
    graph = {}
    for a, b in edges:
        graph.setdefault(a, set()).add(b)
        graph.setdefault(b, set()).add(a)
    return graph


class SplitTests(unittest.TestCase):
    def test_transitive_components(self):
        graph = graph_from_edges([(1, 2), (2, 3), (4, 5)])
        self.assertEqual(components(graph, set(graph)), [frozenset({1, 2, 3}), frozenset({4, 5})])

    def test_articulation_and_disabled_splitting(self):
        graph = graph_from_edges([(0, 1), (1, 2), (2, 3), (3, 4)])
        eligible = lambda part: len(part) >= 2
        self.assertEqual(find_split(graph, set(graph), eligible, Options(), [0])[0], (2,))
        self.assertFalse(find_split(graph, set(graph), eligible, Options(max_bridge_residues=0), [0])[1])

    def test_two_vertex_separator(self):
        graph = graph_from_edges([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)])
        eligible = lambda part: len(part) >= 2
        self.assertFalse(find_split(graph, set(graph), eligible, Options(), [0])[1])
        self.assertEqual(find_split(graph, set(graph), eligible, Options(max_bridge_residues=2), [0])[0], (0, 3))

    def test_dangling_fragment_cannot_be_second_patch(self):
        graph = graph_from_edges([(0, 1), (1, 2), (2, 0), (2, 3)])
        self.assertFalse(find_split(graph, set(graph), lambda part: len(part) >= 2, Options(), [0])[1])

    def test_search_budget_fails_explicitly(self):
        graph = graph_from_edges([(0, 1), (1, 2), (2, 0)])
        with self.assertRaisesRegex(ValueError, 'limit exceeded'):
            find_split(graph, set(graph), lambda part: True, Options(max_separator_checks=1), [0])


class PatchTests(unittest.TestCase):
    def fixture(self):
        # Three paired sites linked by same-chain A residues, each of which
        # contacts a B residue in its neighboring site. Separators are A:2/4.
        entries = [('A', 1, 0), ('A', 2, 7), ('A', 3, 14), ('A', 4, 21), ('A', 5, 28),
                   ('B', 60, 0), ('B', 136, 14), ('B', 200, 28)]
        residues = {'A': [], 'B': []}
        for i, (chain, number, x) in enumerate(entries):
            residues[chain].append(Residue(chain, number, i, 90 if chain == 'A' else 80,
                                            (x, 0 if chain == 'A' else 3, 0)))
        pae = [[2 if i < j else 10 for j in range(len(entries))] for i in range(len(entries))]
        return residues, pae

    def test_recursive_split_and_parent_scores(self):
        residues, pae = self.fixture()
        options = Options(min_residues_per_chain=1, min_contacts=1)
        rows = detect_patches(residues, pae, options)
        patches = rows[1:]
        self.assertEqual(len(patches), 5)
        self.assertEqual(max(r['depth'] for r in patches), 2)
        self.assertEqual(sum(r['status'] == 'leaf' for r in patches), 3)
        self.assertEqual(patches[0]['separator residues'], 'A:2')
        direct = score_pair('A', 'B', residues, pae, {})
        self.assertEqual(rows[0]['pDockQ2 max'], direct['pDockQ2 max'])
        for row in patches:
            region = {c: [(int(n), int(n)) for n in row[f'interface residues {i}'].split(',')]
                      for i, c in enumerate(('A', 'B'), 1)}
            expected = score_pair('A', 'B', residues, pae, region)
            self.assertEqual(row['pDockQ2 max'], expected['pDockQ2 max'])
        self.assertEqual(sorted(r['rank'] for r in patches), list(range(1, 6)))

    def test_discontinuous_sequence_and_order_invariance(self):
        residues = {'A': [Residue('A', 1, 0, 90, (0, 0, 0))],
                    'B': [Residue('B', 60, 1, 90, (0, 0, 7)), Residue('B', 136, 2, 90, (0, 0, 8))]}
        pae = [[2]*3 for _ in range(3)]
        options = Options(max_bridge_residues=0, min_residues_per_chain=1, min_contacts=1)
        rows = detect_patches(residues, pae, options)
        self.assertEqual(rows[1]['interface residues 2'], '60,136')
        self.assertEqual(rows[1]['contact count'], 2)  # Inclusive 8 A boundary.
        self.assertEqual([rows[1][f'contact distance {key}'] for key in ('min', 'max', 'median', 'mean')],
                         [7.0, 8.0, 7.5, 7.5])
        reverse = {k: list(reversed(v)) for k, v in reversed(list(residues.items()))}
        self.assertEqual(rows, detect_patches(reverse, pae, options, ['B', 'A']))

    def test_adjacent_sequence_can_be_separate_and_no_contacts(self):
        residues = {'A': [Residue('A', 1, 0, 90, (0, 0, 0)), Residue('A', 2, 1, 90, (100, 0, 0))],
                    'B': [Residue('B', 1, 2, 90, (0, 0, 7)), Residue('B', 2, 3, 90, (100, 0, 7))]}
        pae = [[2]*4 for _ in range(4)]
        options = Options(min_residues_per_chain=1, min_contacts=1)
        self.assertEqual(len(detect_patches(residues, pae, options)), 3)
        rows = detect_patches(residues, pae, Options(cutoff=6))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['pDockQ2 max'], 0)
        self.assertTrue(all(rows[0][f'contact distance {key}'] == ''
                            for key in ('min', 'max', 'median', 'mean')))

    def test_nearest_contact_weights_each_residue_once(self):
        residues = {'A': [Residue('A', 1, 0, 90, (0, 0, 0)),
                          Residue('A', 2, 1, 90, (0, 0, 2))],
                    'B': [Residue('B', 1, 2, 90, (0, 0, 3)),
                          Residue('B', 2, 3, 90, (0, 0, 9))]}
        rows = detect_patches(residues, [[2]*4 for _ in range(4)],
                              Options(max_bridge_residues=0, min_residues_per_chain=1, min_contacts=1))
        # Included pairs have distances 3,1,7. B minima are 1,7;
        # A minima are 3,1. A1-B2=9 is excluded.
        row = rows[1]
        self.assertEqual(row['contact count'], 3)
        self.assertAlmostEqual(row['contact distance mean'], 11/3)
        self.assertEqual(row['nearest contact distance mean 2'], 4)
        self.assertEqual(row['nearest contact distance median 2'], 4)
        self.assertEqual(row['nearest contact distance max 2'], 7)
        self.assertEqual(row['nearest contact distance mean 1'], 2)
        empty = detect_patches(residues, [[2]*4 for _ in range(4)], Options(cutoff=0.5))[0]
        self.assertEqual(empty['nearest contact distance mean 2'], '')

    def test_size_filter_and_invalid_options(self):
        residues, pae = self.fixture()
        self.assertEqual(len(detect_patches(residues, pae, Options(min_contacts=100))), 1)
        for kwargs in ({'cutoff': float('nan')}, {'cutoff': float('inf')}, {'cutoff': 0},
                       {'max_bridge_residues': -1}, {'min_contacts': 0}, {'min_residues_per_chain': 0}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                Options(**kwargs)
        with self.assertRaises(ValueError):
            detect_patches(residues, pae, chains=['A', 'X'])

    def test_tsv_schema(self):
        residues, pae = self.fixture()
        output = io.StringIO()
        write_patches(detect_patches(residues, pae), output)
        self.assertIn('algorithm version', output.getvalue())
        self.assertIn('separator residues', output.getvalue())


CIF = '''data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 1 ? 0 0 0 1 90 1 A 1
ATOM 2 C CA . ALA B 2 1 ? 0 0 7 1 80 1 B 1
#
'''
DATA = {'token_chain_ids': ['A', 'B'], 'token_res_ids': [1, 1], 'pae': [[0, 2], [10, 0]]}


class InputTests(unittest.TestCase):
    def test_macos_metadata_is_not_a_model(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'x_model_0.cif').write_text(CIF)
            (root/'x_full_data_0.json').write_text(json.dumps(DATA))
            (root/'._x_model_0.cif').write_bytes(b'not a CIF')
            (root/'._x_full_data_0.json').write_bytes(b'not JSON')
            self.assertEqual([m[1] for m in input_models(root)], [0])
            archive = root/'x.zip'
            with zipfile.ZipFile(archive, 'w') as stream:
                stream.writestr('x_model_0.cif', CIF)
                stream.writestr('x_full_data_0.json', json.dumps(DATA))
                stream.writestr('._x_model_0.cif', 'metadata')
                stream.writestr('__MACOSX/._x_model_0.cif', 'metadata')
                stream.writestr('__MACOSX/x_model_0.cif', 'metadata')
            self.assertEqual([m[1] for m in input_models(archive)], [0])
            (root/'x_model_0.cif').unlink()
            (root/'._x.zip').write_bytes(b'not a ZIP')
            self.assertEqual([m[1] for m in input_models(root)], [0])

    def test_model_zip_directory_and_missing_data(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cif, data, archive = root/'x_model_0.cif', root/'x_full_data_0.json', root/'x.zip'
            cif.write_text(CIF)
            data.write_text(json.dumps(DATA))
            with zipfile.ZipFile(archive, 'w') as stream:
                for n in range(2):
                    stream.writestr(f'x_model_{n}.cif', CIF)
                    stream.writestr(f'x_full_data_{n}.json', json.dumps(DATA))
            options = Options(min_residues_per_chain=1, min_contacts=1)
            single = analyze(cif, data, options)
            self.assertEqual(len(single), 2)
            self.assertEqual(single[1]['contact distance median'], 7.0)
            self.assertEqual(single[1]['contact distance mean'], 7.0)
            zipped = analyze(archive, options=options)
            self.assertEqual(len(zipped), 4)
            extracted = analyze(root, options=options)
            self.assertEqual(len(extracted), 2)
            self.assertEqual(extracted[0]['source'], str(root.resolve()))
            self.assertEqual(extracted[1]['pDockQ2 max'], single[1]['pDockQ2 max'])
            self.assertEqual(single[1]['pDockQ2 max'], zipped[1]['pDockQ2 max'])
            with self.assertRaises(ValueError):
                list(input_models(cif))
            with self.assertRaises(ValueError):
                list(input_models(archive, data))
            with zipfile.ZipFile(root/'bad.zip', 'w') as stream:
                stream.writestr('x_model_0.cif', CIF)
            with self.assertRaises(ValueError):
                analyze(root/'bad.zip')
            bad = dict(DATA, pae=[[0, float('nan')], [2, 0]])
            data.write_text(json.dumps(bad))
            with self.assertRaisesRegex(ValueError, 'PAE must'):
                analyze(cif, data)

    def test_extracted_directory_order_and_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(ValueError, 'no model CIF'):
                list(input_models(root))
            for n in (10, 2):
                (root/f'x_model_{n}.cif').write_text(CIF)
                (root/f'x_full_data_{n}.json').write_text(json.dumps(DATA))
            self.assertEqual([m[1] for m in input_models(root)], [2, 10])
            with self.assertRaisesRegex(ValueError, 'only valid with CIF'):
                list(input_models(root, 'data.json'))
            (root/'x_full_data_2.json').unlink()
            with self.assertRaisesRegex(ValueError, 'requires matching'):
                list(input_models(root))
            (root/'y_model_2.cif').write_text(CIF)
            with self.assertRaisesRegex(ValueError, 'unique model numbers'):
                list(input_models(root))

    def test_zip_only_directory_remains_supported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = root/'x.zip'
            with zipfile.ZipFile(archive, 'w') as stream:
                stream.writestr('x_model_0.cif', CIF)
                stream.writestr('x_full_data_0.json', json.dumps(DATA))
            self.assertEqual(list(input_models(root)), list(input_models(archive)))


if __name__ == '__main__':
    unittest.main()
