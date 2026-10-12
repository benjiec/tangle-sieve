"""Standalone PyRosetta worker, executed inside the pinned Docker image.

No PyRosetta import at module scope, so schema and validation remain testable
on the host. Measurements do not relax or repack coordinates.
"""
import csv
import json
import math
from pathlib import Path

ENERGY_TERMS = ('fa_atr', 'fa_rep', 'fa_sol', 'fa_elec', 'lk_ball_wtd',
                'hbond_sr_bb', 'hbond_lr_bb', 'hbond_bb_sc', 'hbond_sc',
                'rama_prepro', 'p_aa_pp', 'omega', 'fa_dun')
PROTOCOL = dict(score_function='ref2015', decompose_bb_hb_into_pair_energies=True, repacking=False, relaxation=False,
                probe_radius=1.4, nonpolar_definition='C/S with abs(partial charge) <= 0.4, including attached H surface',
                unsat_calculator='BuriedUnsatisfiedPolarsCalculator',
                unsat_burial_cutoff=0.1, unsat_generous=True, unsat_legacy=False,
                unsat_vsasa=False, packstat_oversample=100, random_seed=111111)
MEASURES = ('nres', 'sasa', 'nonpolar_sasa', 'buried_unsat_heavy', 'packstat',
            'total_score', 'score_per_residue') + ENERGY_TERMS


def nonpolar(element, charge):
    return element in ('C', 'S') and abs(charge) <= 0.4


def write_table(path, rows):
    if not rows:
        raise ValueError('empty diagnostic table')
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def check_finite(value):
    if isinstance(value, dict):
        for item in value.values():
            check_finite(item)
    elif isinstance(value, list):
        for item in value:
            check_finite(item)
    elif isinstance(value, float) and not math.isfinite(value):
        raise ValueError('nonfinite integrity result')


def unsatisfied_atoms(pose, r):
    """Use Rosetta's calculator, selecting atoms to expose its unbound atom map."""
    calc = r.protocols.simple_pose_metric_calculators.BuriedUnsatisfiedPolarsCalculator(
        'integrity_sasa', 'integrity_hbonds', 0.1, True, False, False)
    total = int(calc.get('all_heavy_unsats', pose))
    counts = json.loads(calc.get('residue_bur_unsat_polars', pose))
    selection = r.core.id.AtomID_Map_bool_t()
    selection.resize(pose.size())
    for i in range(1, pose.size() + 1):
        selection.resize(i, pose.residue(i).natoms(), False)
    result = set()
    for i, count in enumerate(counts, 1):
        if not count:
            continue
        residue = pose.residue(i)
        for a in range(1, residue.nheavyatoms() + 1):
            atom_type = residue.atom_type(a)
            if not (atom_type.is_donor() or atom_type.is_acceptor()):
                continue
            atom_id = r.core.id.AtomID(a, i)
            selection.set(atom_id, True)
            calc.set_special_region(selection)
            calc.notify_energy_change()
            if int(calc.get('all_heavy_unsats', pose)):
                result.add((i, a))
            selection.set(atom_id, False)
    if len(result) != total:
        raise ValueError('atom-level buried-unsatisfied count disagrees with Rosetta total')
    return result


def measure(pose, score, r, state):
    print(f'Measuring {state}: {pose.size()} residues', flush=True)
    total = score(pose)
    atom_sasa = r.core.id.AtomID_Map_double_t()
    residue_sasa = r.utility.vector1_double()
    sasa = r.core.scoring.calc_per_atom_sasa(pose, atom_sasa, residue_sasa, 1.4)
    unsats = unsatisfied_atoms(pose, r)
    weights = score.weights()
    terms = {name: getattr(r.core.scoring, name) for name in ENERGY_TERMS}
    rows, atoms = [], []
    for i in range(1, pose.size() + 1):
        residue = pose.residue(i)
        identity = dict(state=state, chain=pose.pdb_info().chain(i),
                        residue_number=pose.pdb_info().number(i),
                        insertion_code=pose.pdb_info().icode(i).strip(), residue=residue.name3())
        nonpolar_area = 0.0
        hydrogen_sasa = {}
        for h in range(residue.nheavyatoms()+1, residue.natoms()+1):
            parent = residue.atom_base(h)
            hydrogen_sasa[parent] = hydrogen_sasa.get(parent, 0.0) + atom_sasa[r.core.id.AtomID(h, i)]
        for a in range(1, residue.nheavyatoms() + 1):
            area = atom_sasa[r.core.id.AtomID(a, i)]
            group_area = area + hydrogen_sasa.get(a, 0.0)
            is_nonpolar = nonpolar(residue.atom_type(a).element(), residue.atomic_charge(a))
            if is_nonpolar:
                nonpolar_area += group_area
            atoms.append(dict(identity, atom=residue.atom_name(a).strip(), sasa=area, group_sasa=group_area,
                              nonpolar=is_nonpolar, buried_unsat=(i, a) in unsats))
        energies = pose.energies().residue_total_energies(i)
        rows.append(dict(identity, sasa=residue_sasa[i], nonpolar_sasa=nonpolar_area,
                         buried_unsat_heavy=sum((i, a) in unsats for a in range(1, residue.nheavyatoms()+1)),
                         total_score=pose.energies().residue_total_energy(i),
                         **{name: energies[term] * weights[term] for name, term in terms.items()}))
    energy_map = pose.energies().total_energies()
    summary = dict(nres=pose.size(), sasa=sasa,
                   nonpolar_sasa=sum(row['nonpolar_sasa'] for row in rows),
                   buried_unsat_heavy=len(unsats),
                   packstat=r.core.scoring.packstat.compute_packing_score(pose, 100),
                   total_score=total, score_per_residue=total / pose.size(),
                   **{name: energy_map[term] * weights[term] for name, term in terms.items()})
    return summary, rows, atoms


def residue_key(row):
    return row['chain'], row['residue_number'], row['insertion_code']


def main():
    import pyrosetta
    from pyrosetta import rosetta as r
    pyrosetta.init('-mute all -constant_seed -jran 111111 -pose_metrics:sasa_calculator_probe_radius 1.4')
    pose = pyrosetta.pose_from_file('model.cif')
    expected = json.loads(Path('expected.json').read_text())
    actual = []
    for i in range(1, pose.size()+1):
        res = pose.residue(i)
        if not res.is_protein():
            raise ValueError('integrity analysis currently supports protein-only models')
        # Compare heavy-atom names and coordinates: reject silent residue loss or rebuilding.
        actual.append(dict(chain=pose.pdb_info().chain(i), number=pose.pdb_info().number(i),
                           icode=pose.pdb_info().icode(i).strip(), name=res.name3(),
                           atoms={res.atom_name(a).strip(): [res.xyz(a).x, res.xyz(a).y, res.xyz(a).z]
                                  for a in range(1, res.nheavyatoms()+1)
                                  if not res.atom_type(a).is_virtual()}))
    if len(actual) != len(expected):
        raise ValueError('Rosetta changed the input residue count')
    for observed, original in zip(actual, expected):
        if any(observed[k] != original[k] for k in ('chain', 'number', 'icode', 'name')):
            raise ValueError('Rosetta changed residue identities or numbering')
        # Terminal OXT may be added by Rosetta; all original heavy atoms must be retained.
        for atom, xyz in original['atoms'].items():
            if atom not in observed['atoms'] or max(abs(a-b) for a,b in zip(xyz, observed['atoms'][atom])) > 0.002:
                raise ValueError('Rosetta changed an input heavy atom or its coordinates')
        extras = set(observed['atoms']) - set(original['atoms'])
        if extras - {'OXT'}:
            raise ValueError(f'Rosetta rebuilt missing heavy atoms: {extras}')
    score = pyrosetta.create_score_function('ref2015')
    options = score.energy_method_options()
    options.hbond_options().decompose_bb_hb_into_pair_energies(True)
    score.set_energy_method_options(options)
    complex_metrics, residues, atoms = measure(pose, score, r, 'complex')
    complex_rows = {residue_key(row): row for row in residues}
    summaries = [dict(scope='complex', chain='ALL', **complex_metrics,
                      nonpolar_sasa_in_complex=complex_metrics['nonpolar_sasa'],
                      nonpolar_shielded_by_partners=0.0)]
    exposure = []
    for index in range(1, pose.num_chains()+1):
        chain_pose = pose.split_by_chain(index)
        chain_id = chain_pose.pdb_info().chain(1)
        metrics, chain_rows, chain_atoms = measure(chain_pose, score, r, 'isolated')
        bound = sum(complex_rows[residue_key(row)]['nonpolar_sasa'] for row in chain_rows)
        summaries.append(dict(scope='chain', chain=chain_id, **metrics,
                              nonpolar_sasa_in_complex=bound,
                              nonpolar_shielded_by_partners=metrics['nonpolar_sasa']-bound))
        for row in chain_rows:
            other = complex_rows[residue_key(row)]
            exposure.append(dict(chain=chain_id, residue_number=row['residue_number'],
                                 insertion_code=row['insertion_code'], residue=row['residue'],
                                 isolated_nonpolar_sasa=row['nonpolar_sasa'],
                                 complex_nonpolar_sasa=other['nonpolar_sasa'],
                                 partner_shielded_nonpolar_sasa=row['nonpolar_sasa']-other['nonpolar_sasa'],
                                 isolated_unsat_heavy=row['buried_unsat_heavy'],
                                 complex_unsat_heavy=other['buried_unsat_heavy']))
        residues.extend(chain_rows)
        atoms.extend(chain_atoms)
    # Complex shielding compares all isolated chains to the assembled complex.
    summaries[0]['nonpolar_shielded_by_partners'] = sum(s['nonpolar_shielded_by_partners'] for s in summaries[1:])
    result = dict(protocol=PROTOCOL, pyrosetta_version=pyrosetta.version().splitlines()[-1], summaries=summaries)
    check_finite([result, residues, atoms, exposure])
    write_table('residues.tsv', residues)
    write_table('atoms.tsv', atoms)
    write_table('exposure.tsv', exposure)
    Path('integrity.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
