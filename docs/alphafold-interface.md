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
selected model. Domain mapping requires either `--chains` with a `--domain`
assignment for each selected chain or both `--manifest` and `--domain-file`.
With manifest lookup and no `--chains`, the command evaluates chain pairs
that match a registered interface profile.
Explicit chain assignments take precedence over domain-file lookup. The domain
TSV uses case-insensitive `Sequence Name` and `Domain Type` headers and exact
values. Additional columns are allowed. Missing or unknown domain assignments
and ambiguous profile matches are errors. An explicit `--chains` pair must
match a profile.

## TSV layout and shared measurements

Each patch has two rows, identified by `patch_id` and `sequence_index`.
Every export starts with the same identity, chain, shared measurement, and
call columns. A death-only or TRAF-only run appends that profile's columns.
A mixed run appends both blocks in registry order. Unavailable measurements
and inapplicable profile fields are empty TSV fields. The shared `interface`
column reports `Canonical`, `Atypical`, or `Variant` when classified, and is
empty otherwise. The shared `candidate face` column gives the profile's
face or pose on the corresponding chain row. All profiles use `High`,
`Medium`, `Low`, or `Reject` for `confidence`, and `true` or `false` for
`chemical feasibility`. Both fields use `NULL` for unavailable calls.

Shared measurements include per-partner burial, tight-contact residue counts
at 4.5 Å, region mean pLDDT, and Cα distance summaries. `contact basis`
identifies the residue pairs used for the Cα summaries. TRAF selects pairs
by heavy-atom distance and the death profile selects them by Cα distance.
`pDockQ2 contact count` and `pDockQ2 distance cutoff` describe the
independent scoring selection. Interface structure comes from DSSP, with
unassigned positions marked in the summary.

Hydrogen-bond candidates are counted by atom pair. The shared qualified
bond columns use a donor–hydrogen–acceptor angle of at least 120° and
report unresolved donor geometry separately. The new nonpolar, salt, and
clash residue-pair columns use only the measured selection. Clash pairs
exclude geometry-qualified hydrogen-bond atom overlaps. The legacy contact
columns retain their existing counting rules, including the death profile's
expanded salt-candidate neighborhood. In the death profile, burial is
measured in the full chain-pair context and summed over patch residues. In
TRAF, burial, tight-contact counts, and chemistry use the selected receptor
region and `traf_chemistry_window`. Contact counts and Cα distances,
structure, hydropathy, and confidence summarize the original geometry window
reported in `region`.

## Death-family profile

The death profile accepts DD, DED, and CARD on either side. It detects patches at
10 Å, and normalizes helices for classification with a one-residue maximum
connector and an angle strictly below 75°. The interface field comes from the
face classifier. The confidence and chemical feasibility calculations use
unrounded metrics. The electrostatic predicate uses the magnitude of the signed
charge difference. Feasibility uses 100 Å when the minimum salt-bridge candidate
distance is unavailable. The raw distance column is empty.

DSSP discovery uses `--dssp`, then `SIEVE_DSSP`, then `mkdssp` on PATH. When DSSP
assignments are unavailable, the interface result is empty and the face status
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

## TRAF–SLiM interfaces

Annotate the receptor chain as `TRAF_MATH`, `TRAF2`, `TRAF3`, `TRAF5`,
`TRAF6`, or `MATH`, and the peptide chain as `SLiM`. The paralog
annotation is independent of pose classification.

```sh
python scripts/alphafold-interface.py fold.zip --output traf.tsv \
  --chains A B --domain A:TRAF_MATH --domain B:SLiM
```

Without `--math-window`, use an input construct with a complete MATH domain
at its C terminus. The evaluator scans DSSP β segments backward from the
last strand and checks their order and antiparallel Cα contacts to assign the
MATH landmarks. It can join short interruptions in a strand and scan past
terminal helices and loops. If DSSP has no assignments for the receptor,
the backbone hydrogen-bond detector is used instead. Use
`--math-window A:50-201` for known boundaries, tags, or extensions.

If automatic MATH assignment is unresolved, the evaluator searches the full
annotated receptor chain and reports the measured contacts, burial, pLDDT,
PAE, and pDockQ2 whenever the inputs support them. `region` records the
scoring span and `traf_status` identifies the full-chain fallback. Confidence
and chemical feasibility describe the selected patch and use the usual
thresholds when the inputs support them. The inferred MATH window and pose
classification are empty.

Use `--motif-window B:24-32` to supply an inclusive peptide interval for
scoring, pose determination, and chemistry. Otherwise the evaluator finds
peptide contacts within 6 Å and selects a compact five- or six-residue
`traf_motif_window` for scoring. Pose determination uses all peptide residues
in those 6 Å heavy-atom contacts and records their ranges in `traf_pose_window`.
Separated contacting runs are assessed independently. Competing poses or
registers leave the MATH face unresolved.

Pose classes use backbone and side-chain geometry from the input model.
`traf6` requires an ordered antiparallel β7 backbone ladder with at least
three hydrogen bonds spanning at least four sequence steps on both partners.
`traf2/3` requires a packing anchor followed by a β6–β7 loop polar anchor
within a five-residue span. A polar anchor can be supported by a qualified
hydrogen bond, a salt bridge, or two compatible side-chain contacts within
3.8 Å. The last case supplies pose evidence only. The MATH row puts
`traf6` or `traf2/3` in `candidate face` when the pose is assigned.
The SLiM row puts the strongest matching rule per family in `candidate face`.
The sequence check starts one residue before the TRAF6 β7 register or
at the TRAF2/3 packing anchor. It can read adjacent sequence outside the
scored window. For an unresolved pose with a selected contact site, a unique
motif supported by contacts can still appear on the SLiM row.

| SLiM rule | Pattern | Interface label with a matching MATH pose |
| --- | --- | --- |
| `traf6_binding_strict` | `P.E..[FYWDE]` | `Canonical` |
| `traf6_binding_permissive` | `P.[QE]..[FYWHDEILVKR]` | `Variant` |
| `traf2_binding_strict` | `P.[QE].[ST]` | `Canonical` |
| `traf2_binding_major_consensus` | `[PSAT].[QE]E` | `Canonical` |
| `traf2_binding_variant` | `P.[QEN].[STEDAG]` | `Canonical` |
| `traf2_binding_permissive` | `P.{1,2}[QEN].{1,2}[STEDAG]` | `Variant` |

Rules are checked in the order shown within each family. An atypical pose
or a named pose without a same-family motif exports `Atypical` when there
is enough evidence for the call. Unresolved poses, no-contact results, and
insufficient sequence evidence for a named pose leave `interface` empty.
See `traf_status` and `traf_pose_evidence` for the structural diagnostics.

The automatic selector evaluates feasibility on the geometry window first.
After pose selection, a match in the same family sets `traf_chemistry_window`
to the exact motif span and recalculates feasibility and its input metrics.
For a rule with multiple matching lengths, the shortest matching span is
used. An explicit `--motif-window` is used as supplied for chemistry.

Both directional pDockQ2 scores use the existing 8 Å Cα scoring contacts.
TRAF confidence uses `pDockQ2 max` and assigns `High` at 0.70 or higher,
`Medium` at 0.50 or higher, `Low` at 0.23 or higher, and `Reject` below 0.23.
Chemical feasibility uses summed buried area and geometric contacts.
The provisional area floors are 250 Å² for mixed support, 400 Å² for
packing-only support, and 150 Å² for backbone pairing, with no upper area
bound. Override these floors using `--traf-threshold NAME=VALUE` with
`mixed_area`, `packing_area`, or `paired_area`.

With manifest-driven evaluation of all annotated pairs,
`traf_other_protomers` lists other evaluated MATH chains contacting the same
motif. Reference structures are stored only under `tests/fixtures/traf_slim`.
