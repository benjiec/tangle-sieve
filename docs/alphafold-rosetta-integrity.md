# AlphaFold fold-integrity diagnostics

```bash
sieve-py ~/git/tangle/sieve/scripts/alphafold-list-results.py \
  ~/git/tangle/coral/notes/af_manifest.tsv b911edee --highest | \
  sieve-py ~/git/tangle/sieve/scripts/alphafold-rosetta-integrity.py -
```

The installed `rosettacommons/rosetta:serial` image supplies PyRosetta. No host
PyRosetta installation is needed. Plain extracted-directory or ZIP/nested-ZIP
locators require `--model NUMBER`. Stdin accepts the same locator/result-row
contract as the interface wrapper. Multiple rows select multiple models.

Default output is one metric per line, grouped by model and scope. `--tsv`
prints a header and one row for the complex (`scope=complex, chain=ALL`), then
one row per isolated chain (`scope=chain`). Each row includes both hash IDs,
model number, and locator. Chain metrics describe the isolated conformation;
`nonpolar_sasa_in_complex` separately reports its exposure within the assembly.
No interface argument is required.

`--output DIRECTORY` optionally preserves artifacts in a new or empty directory.
Without it, temporary working files are removed on success or failure; error
log tails are copied into the stderr message before removal. Docker output never
pollutes summary stdout. Processing stops on the first error with nonzero status;
the combined summary is printed only after every model succeeds. Empty stdin
runs no Docker jobs (TSV mode prints the header).

## Measurements and definitions

- `sasa`: solvent-accessible surface area in Å², with a 1.4 Å water probe,
  using Rosetta `calc_per_atom_sasa` with explicit hydrogen atoms.
- `nonpolar_sasa`: exposed area of C/S atoms with absolute Rosetta partial charge
  <=0.4, **including their attached hydrogen surface**. This is an explicit
  operational definition, not a count of hydrophobic amino acids. It need not
  exactly match InterfaceAnalyzer's hydrophobic SASA implementation.
- `nonpolar_sasa_in_complex`: remaining exposed nonpolar area in the complex.
- `nonpolar_shielded_by_partners`: isolated minus complex nonpolar exposure.
  For the complex row, this is summed over chains. This measures protection by
  other chains, not the amount of hydrophobic core inside an individual chain.
- `buried_unsat_heavy`: absolute number of buried unsatisfied heavy polar atoms
  from Rosetta `BuriedUnsatisfiedPolarsCalculator`. Explicit settings: burial
  cutoff 0.1 Å², generous hydrogen-bond counting, legacy counting off, VSASA off,
  surface residues not skipped (calculator default). This is **not** the
  interface wrapper's delta count. Hydrogen-only unsatisfied sites are excluded.
  Atom selections are checked against the calculator's whole-pose total.
- `packstat`: Rosetta whole-pose packing statistic, oversampling 100. Chain rows
  use isolated chains; the complex row uses the assembly. This is a packing
  diagnostic, not an explicit cavity-volume measurement.
- `total_score` and `score_per_residue`: REF2015 score and score divided by
  residue count. These are not folding free energies or cross-sequence rankings.
- Energy components are **weighted** Rosetta energy units: `fa_atr`, `fa_rep`,
  `fa_sol`, `fa_elec`, `lk_ball_wtd`, four hydrogen-bond terms, `rama_prepro`,
  `p_aa_pp`, `omega`, and `fa_dun`. These selected components do not sum to total
  score because additional REF2015 terms exist. `fa_rep` is a repulsion
  diagnostic, not a literal count of clashes. Conformational terms are not
  binary torsion/rotamer outlier calls.

Backbone hydrogen-bond energies are decomposed into pair energies so residue
energy contributions sum to the whole-pose score.

The protocol does not relax, minimize, or repack. Rosetta adds hydrogens and may
add terminal OXT atoms or virtual atoms. All original physical heavy-atom
coordinates and residue identities are checked after import; unsupported
residue loss, real heavy-atom rebuilding, and altered input coordinates fail.
Only protein-only inputs with single-character alphanumeric chain IDs and no
alternate conformations are supported. Chain isolation preserves the modeled
conformation; it does not search for an unbound fold. Seed is fixed at 111111;
image metadata, worker hash, input hash, protocol, and command are recorded.

## Saved artifacts

The root contains `summary.tsv`, `image.json`, and `image-inspect.log`.
Each numbered model directory contains:

- `model.cif`, `expected.json`, and a copy of `worker.py` for reproducibility.
- `settings.json`, `run.log`, and `integrity.json`.
- `residues.tsv`: every residue in complex and isolated states, with exposure,
  unsatisfied counts, total score, and selected energy terms.
- `atoms.tsv`: heavy-atom SASA, heavy-atom plus attached-H `group_sasa`, nonpolar
  classification, and buried-unsatisfied flags in each state. Numbering and
  insertion codes follow input coordinate identities. Virtual atoms, if present
  in Rosetta's heavy-atom list, have no physical surface contribution.
- `exposure.tsv`: aligned per-residue isolated/complex nonpolar exposure,
  partner shielding, and absolute unsatisfied heavy-polar counts in both states.

To inspect a long helix, examine its residues in `exposure.tsv`, identify
unsatisfied atoms using `atoms.tsv`, and inspect energy terms in `residues.tsv`.
Large exposed areas and unfavorable components are investigation targets, not
validated pass/fail cutoffs. Exposure alone does not quantify internal-core
burial or prove a native topology. Compare alternatives of the same sequence
under the same preparation protocol; absolute totals depend on protein size.

## Validation

Host tests cover input selection, result validation, output formats, persistence,
cleanup, and failure handling. The optional Docker calculation control is:

```bash
docker run --rm --network none --mount type=bind,source="$PWD",target=/repo,readonly \
  rosettacommons/rosetta:serial \
  python3 /repo/tests/sieve/rosetta_integrity_container_checks.py
```

Run from the `sieve` repository. It checks atom/residue sums and verifies that
well-separated duplicate chains have additive SASA and energy. Real-model
integration has been checked on `fc0f972f`, model 0. These validate calculations
and loading, not biological accuracy of the model.

Sources: [REF2015](https://doi.org/10.1021/acs.jctc.7b00125),
[Rosetta buried-polar calculator](https://github.com/RosettaCommons/rosetta/blob/main/source/src/protocols/simple_pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.cc),
[packing statistic](https://docs.rosettacommons.org/docs/latest/application_documentation/analysis/interface-analyzer).

## Derived assemblies

An `alphafold-overlap.py` output directory is a supported locator. Select `--model 0`
or pipe a matching list-results row (without `--highest`, since assemblies have no
ipTM/pTM). The loader verifies the CIF/provenance hashes and uses the derived
fold ID. Coordinate analysis and preparation settings are unchanged. See
[assembly documentation](alphafold-overlap.md).
