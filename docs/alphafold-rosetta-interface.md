# Rosetta interface analysis

Analyze a selected AlphaFold model with the installed
`rosettacommons/rosetta:serial` Docker image:

```bash
sieve-py ~/git/tangle/sieve/scripts/alphafold-list-results.py \
  ~/git/tangle/coral/notes/af_manifest.tsv b911edee --highest | \
  sieve-py ~/git/tangle/sieve/scripts/alphafold-rosetta-interface.py - \
    --interface A:B --output results/b911edee
```

For a plain locator, specify `--model 0`. Both extracted fold directories and
ZIP/nested-ZIP locators use the existing fold loader. Stdin accepts the shared
five-column list-results format or one locator per record. Multiple records
produce one analysis each. Conflicting `--model` selections are errors.

`--output DIRECTORY` is optional. When provided, it must be new or empty.
When omitted, working files are created in a temporary directory and removed
after success or failure; only the stdout summary is retained. On failure,
log tails are included in the stderr error message before cleanup. Existing results are
never overwritten. Each numbered model directory contains `model.cif`,
`interface.sc`, `run.log`, and `settings.json`. The output root contains image
metadata (including available registry digests), the Docker inspection log, and
`summary.tsv`. The installed image ID is pinned for all container executions.
Failed executions retain their logs and settings; processing stops at the first
error with nonzero status. Stdout and the combined summary are emitted only when
all requested analyses succeed. Empty stdin runs no containers; with `--tsv` it
prints only the header.

The default stdout summary prints identification, protocol settings, and one
metric per line. `--tsv` prints headers and one row per model, matching
`summary.tsv`. Missing, malformed, or nonfinite required metrics are errors.

| Metric | Meaning |
| --- | --- |
| dG_separated | Rosetta energy difference between complex and separated partners |
| dG_separated/dSASAx100 | Energy difference normalized by buried area, times 100 |
| dSASA_int | Total buried SASA across both sides, Å² |
| dSASA_hphobic | Hydrophobic buried SASA, Å² |
| dSASA_polar | Polar buried SASA, Å² |
| delta_unsatHbonds | Change in buried unsatisfied polar count |
| hbonds_int | Interchain hydrogen-bond count |
| sc_value | Shape complementarity |
| packstat | Packing statistic |
| nres_int | Interface residue count |

Protocol: `ref2015`, input and separated repacking on by default, shape complementarity
and packstat enabled, packstat oversampling 100, fixed random seed 111111.
The interface syntax is two distinct single-character alphanumeric chain IDs,
e.g. `A:B`. Chains are checked against CIF coordinate chain IDs. Two-chain
inputs use the executable's automatic pair selection; larger inputs pass
Rosetta's `A_B` interface option. General chain-group syntax is not supported.

Both `pack_input` and `pack_separated` are enabled by default. Use `--no-repack`
for the unrepacked baseline; explicit `--repack` remains supported. Repacking selects
side-chain conformations in the complex and separated partners; it does not
perform backbone relaxation or minimization. Both settings appear in the text/TSV
summary and saved settings. Compare protocols using separate output directories:

```bash
afl b911edee --highest | sieve-py ~/git/tangle/sieve/scripts/alphafold-rosetta-interface.py - \
  --interface A:B --repack --output results/b911edee-repacked
```

Repacking may change scores in either direction; it does not guarantee favorable
association. Keep preparation consistent across candidate/control comparisons.

With `--no-repack`, scores describe a fixed-coordinate baseline. Neither protocol
provides measured binding free
energies or a feasibility verdict. Preparation sensitivity requires separate
validation. The initial integration test covered a two-chain AF CIF; larger
complexes require additional validation.

Reference: [Rosetta InterfaceAnalyzer documentation](https://docs.rosettacommons.org/docs/latest/application_documentation/analysis/interface-analyzer).
