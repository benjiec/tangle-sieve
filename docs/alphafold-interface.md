# AlphaFold interface evaluation

`alphafold-interface.py` writes a tall TSV of leaf patches with confidence,
chemical feasibility, and interface classification. It shares patch detection and
descriptors with the report-only `alphafold-pdockq2-patches.py`.

Supply domains explicitly to evaluate a selected pair:

```sh
alphafold-interface.py predictions.zip --chains A C \
  --domain A:DD --domain C:CARD --output interfaces.tsv
```

For batches with different chain annotations, pipe result rows and resolve names
through the manifest:

```sh
alphafold-find.py manifest.tsv names.txt |
  alphafold-interface.py - --manifest manifest.tsv \
    --domain-file sequences.tsv --output interfaces.tsv
```

The input also accepts a fold locator, extracted fold directory, or CIF
(with `--full-data`). Stdin accepts one locator per line or result rows with a
selected model. Without `--chains`, the command evaluates all chain pairs.
Explicit chain assignments take precedence over domain-file lookup. The domain
TSV uses case-insensitive `Sequence Name` and `Domain Type` headers and exact
values. Additional columns are allowed. The registry reports an error for a
missing assignment or a domain pair with zero or multiple matching profiles.

## Death-family profile

The death profile accepts DD, DED, and CARD on either side. It detects patches at
10 Å, and normalizes helices for classification with a one-residue maximum
connector and an angle strictly below 75°. The interface field comes from the
face classifier. The confidence and chemical feasibility calculations use
unrounded metrics. The electrostatic predicate uses the magnitude of the signed
charge difference. Each unavailable result is written as `NULL`.

DSSP discovery uses `--dssp`, then `SIEVE_DSSP`, then `mkdssp` on PATH. When DSSP
assignments are unavailable, the interface result is `NULL` and the face status
records the reason.

The Python rule table in `sieve.interface.death` contains the DD, DED, and
CARD coverage heuristics transcribed from `coral/notes/interactome/dd_face_rules.tsv`.
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

## Output

The first columns retain the [patch TSV layout](pDockQ2.md). There are two rows
per leaf patch, in chain order. The rows also include the family, chain domain,
face and status, model confidence scores, rule and normalization settings,
and the three evaluation results. Pair-level results repeat on both rows.

Output is replaced atomically. Model and pair errors are reported on stderr and
return a nonzero exit status. Successful rows are still written if available.

## Kink normalization

Only H segments separated by 1 through the configured maximum residues, all
assigned -, T, or S, qualify. G, I, E, B, P, and missing assignments do not.
The local N→C axes use the last four Cα positions of the first helix and the
first four of the second. For successive coordinate differences d0,d1,d2,
the axis is the normalized cross product (d1−d0) × (d2−d1), oriented N→C.
Missing/nonfinite coordinates, short helices, and degenerate axes prevent a merge.
No pLDDT threshold is used.

Original adjacent pairs are evaluated once. Connected qualifying pairs become
one H element, including connector residues; elements are then renumbered.
Exactly six resulting H elements are still required: a six-to-five merge yields
NA. Accepted original ranges/angles and geometry failures appear in `status 1/2`.
The TSV records `merge kinked helices`, `kink connector length`, and `kink angle`.

Raw DSSP summaries, patch membership, patch IDs, and patch metrics are unchanged.
These thresholds are uncalibrated heuristics. Review normalized assignments
against known domain structures before treating face labels as biological evidence.
