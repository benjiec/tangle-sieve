"""Run explicitly inside the Rosetta image; no host PyRosetta dependency."""
import math
import sys
import tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import pyrosetta
from pyrosetta import rosetta as r
from sieve.rosetta_integrity_worker import measure, ENERGY_TERMS

pyrosetta.init('-mute all -constant_seed -jran 111111 -pose_metrics:sasa_calculator_probe_radius 1.4')
pose = pyrosetta.pose_from_sequence('AAAAAAAA', 'fa_standard')
for i in range(1, pose.size()+1):
    pose.set_phi(i, -60)
    pose.set_psi(i, -45)
    pose.set_omega(i, 180)
# Reload coordinates to normalize terminal torsion caches after synthetic construction.
with tempfile.NamedTemporaryFile(suffix='.pdb') as pdb:
    pose.dump_pdb(pdb.name)
    pose = pyrosetta.pose_from_file(pdb.name)
score = pyrosetta.create_score_function('ref2015')
options = score.energy_method_options()
options.hbond_options().decompose_bb_hb_into_pair_energies(True)
score.set_energy_method_options(options)
single, rows, atoms = measure(pose, score, r, 'single')
assert math.isclose(sum(row['sasa'] for row in rows), single['sasa'], abs_tol=1e-6)
assert math.isclose(sum(row['total_score'] for row in rows), single['total_score'], abs_tol=1e-6)
assert sum(row['buried_unsat'] for row in atoms) == single['buried_unsat_heavy']
assert math.isclose(sum(a['group_sasa'] for a in atoms if a['nonpolar']), single['nonpolar_sasa'], abs_tol=1e-6)
for term in ENERGY_TERMS:
    assert math.isclose(sum(row[term] for row in rows), single[term], abs_tol=1e-6), term
far = pose.clone()
for i in range(1, far.size()+1):
    for a in range(1, far.residue(i).natoms()+1):
        atom_id = r.core.id.AtomID(a, i)
        far.set_xyz(atom_id, far.xyz(atom_id) + r.numeric.xyzVector_double_t(1000, 0, 0))
assembly = pose.clone()
assembly.append_pose_by_jump(far, 1)
combined, _, _ = measure(assembly, score, r, 'separated_control')
assert math.isclose(combined['sasa'], 2*single['sasa'], abs_tol=1e-5)
assert math.isclose(combined['nonpolar_sasa'], 2*single['nonpolar_sasa'], abs_tol=1e-5)
assert combined['buried_unsat_heavy'] == 2*single['buried_unsat_heavy']
assert math.isclose(combined['total_score'], 2*single['total_score'], abs_tol=1e-5), {k: (single[k],combined[k],combined[k]-2*single[k]) for k in ('total_score',)+ENERGY_TERMS}
print('PASS: residue/atom reconciliation, weighted energy sums, separated-chain exposure and energy controls')
