# Death-family patch faces

`alphafold-death-faces.py` independently opens models and uses the same patch
detection library as `alphafold-pdockq2-patches.py`. It does not read patch TSVs.
No intermediate cache is used. Domain assignments and rules are explicit inputs.

```sh
sieve-py ~/git/tangle/sieve/scripts/alphafold-list-results.py manifest.tsv b911edee --highest |
  sieve-py ~/git/tangle/sieve/scripts/alphafold-death-faces.py - \
    --manifest manifest.tsv --domains sequence_domains.tsv \
    --rules ~/git/tangle/sieve/data/death-face-rules.tsv \
    --output faces.tsv --distance-cutoff 10
```

The positional input also accepts a locator, extracted fold directory, ZIP, or
CIF (`--full-data` required). `--chains A B` restricts the pair. Match
`--distance-cutoff`, `--max-bridge-residues`, `--min-residues-per-chain`,
`--min-contacts`, and `--max-separator-checks` between scripts. Defaults match the
patch script. DSSP is required: `--dssp`, then `SIEVE_DSSP`, then `mkdssp` on PATH.
Normal output is the required `--output` TSV, with no stdout summary.

## Assignments and rules

Domain TSV requires `Sequence Name` and `Domain` headers, case-insensitively.
Extra columns are ignored. Values match exactly and case-sensitively: sequence
names are full FASTA headers, including spaces; `DD` and `dd` are different.
Duplicate names, blank required values, and malformed rows are errors.
Unused sequence assignments are permitted. Selected chains must have assignments
whose domains exist in the rule file.

Rule TSV requires `domain`, `face`, `element`, `min_residues`,
`requirement_type`, and `notes`. The example `data/death-face-rules.tsv` contains
the proposed DD and DED heuristics, not calibrated structural rules.
Elements are `alpha1` through `alpha6` and `loop1-2` through `loop5-6`.
Faces are Ia, Ib, IIa, IIb, IIIa, IIIb. A domain may specify a subset of faces.
Duplicate domain/face/element rows are errors. Active minimum counts must be
positive integers; optional counts may be zero.

- `mandatory`: all such rows must pass.
- `alternative_or`: one OR group per domain/face; at least one must pass.
- `optional`: descriptive only; does not affect calling or ranking.
- `notes`: descriptive only. Motifs and conditional statements are not executed.

A full call meets all mandatory requirements and its OR group, if present.
Otherwise any contact with an active element gives `(partial)`. Optional-only
contacts give no face. Full calls rank first; within each group rank by minimum,
then mean fractional element coverage, using all mandatory elements and the
best-covered alternative. Ties follow Ia, Ib, IIa, IIb, IIIa, IIIb. Rank order
is not a probability and does not affect interface classification.

Exactly six DSSP H segments are numbered sequentially. G is permitted in the
intervening connectors with -, T, S. Other connector structures contribute zero
coverage. Missing DSSP assignments in the six-helix span or a different number
of H segments produce `NA` and a status reason. These assumptions do not validate
canonical topology or relative orientation; compare against known structures
before interpreting labels biologically.

Each side uses its own domain rules. Complementary pairs Ia/Ib, IIa/IIb, IIIa/IIIb
work in either orientation, including DD/DED: any full pair yields `Canonical`,
otherwise any complementary partial pair yields `Variant`. Other candidates
produce `Atypical`; either side `none` produces `None`; ineligibility produces
`NA` with precedence over `None`.

## Output and joins

One row per leaf patch:

- `fold_locator`, `fold_id`, `sequence_id`, `model_id`, `patch_id`, `patch label`.
- `chain 1/2`, `sequence name 1/2`, `domain 1/2`.
- `candidate face 1/2`, `candidate interface`.
- `status 1/2`, `rules_id`, and all five patch-detection settings.

Join each side to the tall manifest by locator plus chain ID. Fold and sequence
IDs are checked against the model. CIF author chain IDs are used. A direct CIF
is matched to its parent fold directory in the manifest.

Join to patch metrics by locator/model/chain pair and `patch_id`, checking
settings separately. The patch script names locator/model columns `source` and
`model`. Both scripts use lexically ordered chain pairs.

`patch_id` is the first 32 hexadecimal characters of SHA-256 over canonical JSON
of a version tag, full fold identity payload, model number, and sorted chain IDs
with sorted unique contacting residue sets. AF identity includes seed and the
full sequence/PTM payload, not the short fold ID. Assemblies use the full
provenance identity. Bare CIFs without full seeded identity use a content digest
(and seed if present). Paths, traversal labels, rules, and detection settings
are excluded. The same membership can have different contact pairs/scores under
different cutoffs. This identifies membership, not biological site equivalence
across models. `patch label` retains the readable traversal label (`p1`, etc.).
`rules_id` is SHA-256 of parsed rules sorted by domain/face/element, including notes.

## Errors and migration

Mapping/rule errors stop before annotation or output replacement. Selected model
loading, patch computation, and DSSP failures are reported to stderr; other
models continue, and exit status is nonzero. The output contains successful
models and is replaced atomically, including a header-only file if none succeed.
Model data are retained during preflight to avoid reopening them during annotation.

The patch script no longer accepts `--dd`, `--dd-helix-presence-threshold`, or
`--dd-loop-presence-threshold`, and no longer emits DD face/interface/settings
columns. Use this script and the rule file instead. Both outputs now have
`patch_id`; the former patch TSV `patch id` traversal column is `patch label`.
