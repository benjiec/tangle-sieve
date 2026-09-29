# Spatial interface patches with regional pDockQ2

Algorithm version: **1.0**. Implementation: `sieve/alphafold_pdockq2_patches.py`.
Command: `scripts/alphafold-pdockq2-patches.py`.

This algorithm locates and describes candidate contact sites in a predicted
protein complex. It detects spatially connected interface residues, optionally
splits narrow connections recursively, and scores each resulting residue set
using the existing regional pDockQ2 implementation. It does not search sequence
windows or optimize the physical contact cutoff. Results are observed geometric
patches in a model, not experimentally established binding sites.

## Reusable Methods description

For each chain pair, we identified interchain residue pairs with Cα–Cα distance
≤8 Å. Residues participating in these contacts formed the vertices of an
undirected graph; edges joined any two vertices within 8 Å, including same-chain
pairs. Connected components defined initial interface patches. Components were
recursively split by removing qualifying vertex separators of at most one
residue. A split required at least two resulting components, each retaining at
least three contacting residues per chain and five interchain contacts.
Separators were searched in increasing cardinality and then lexicographic
(chain ID, residue number) order; the first qualifying separator was selected.
Separator residues were retained in the parent report and excluded from child
patches. Each parent and qualifying descendant was scored using mean directional
PAE weights multiplied by mean interface pLDDT, followed by the pDockQ2 sigmoid.
Both directional scores, their maximum, residue identities, and contact counts
were reported. These are regional adaptations of pDockQ2, without a new fit or
validation of the original calibration for automatically selected patches.

This paragraph describes the defaults. Replace the cutoff, separator limit, and
minimum-size values when using other settings. Cite the pDockQ2 publication below
for the scoring equation; cite this versioned specification for patch detection.

## Invocation and parameters

```bash
alphafold-pdockq2-patches.py model.cif --full-data full_data.json --output patches.tsv
alphafold-pdockq2-patches.py predictions.zip --output patches.tsv
alphafold-pdockq2-patches.py ./folds --chains A B --output patches.tsv
```

Run in the installed sieve Python environment, as with the other sieve scripts.
For a repository checkout, use its Python environment with `PYTHONPATH=sieve`
from the workspace root. A directory means all immediate `*.zip` children,
processed in sorted order; subdirectories are not searched.

| Argument | Default | Meaning |
|---|---:|---|
| positional input | required | CIF, ZIP, or directory of ZIPs |
| `--full-data` | none | Matching full-data JSON; required only for CIF input |
| `--output` | optional | TSV path; recomputed and overwritten, no cache; omit for stdout summary only |
| `--chains A B` | all pairs | Analyze one distinct pair; output chain order is lexical |
| `--distance-cutoff` | 8.0 Å | Inclusive cutoff for interchain contacts and all graph edges |
| `--max-bridge-residues` | 1 | Maximum vertex-separator cardinality; 0 disables splitting |
| `--min-residues-per-chain` | 3 | Minimum distinct contacting residues on each chain |
| `--min-contacts` | 5 | Minimum number of interchain residue-pair contacts |
| `--max-separator-checks` | 1,000,000 | Search budget per model and chain pair, shared across roots/descendants |

The minimum sizes are practical user-adjustable filters, not calibrated
biological thresholds. No patch-radius argument or distance sweep is used.
All analysis completes before the output is opened, so input-validation and
search-budget failures do not replace an existing output. Filesystem failures
while writing can still leave a partial file.

## Input data and identity mapping

The command reuses `parse_model` from `alphafold_pdockq2.py`:

1. Parse the CIF using Biopython's `MMCIFParser` with its default author chain
   and residue numbering. Use only the first structural model (`structure[0]`).
2. For each residue containing a Cα atom, store the chain ID, integer residue
   number, Cα Cartesian coordinates, and Cα `B_iso_or_equiv` value. The latter is
   interpreted as AlphaFold pLDDT on the 0–100 scale, **not** a crystallographic
   temperature factor and not an atom-averaged confidence.
3. Read JSON `token_chain_ids` and `token_res_ids`. Their zipped order maps
   `(chain ID, residue number)` to a matrix index. Match each CIF residue to
   that index; never renumber a region or slice PAE by local residue number.
4. Read the square JSON `pae` matrix. For residue tokens i and j, use `pae[i][j]`
   for the direction labeled “1 to 2” and `pae[j][i]` for “2 to 1”. These are
   explicit matrix-index conventions, not a reinterpretation of PAE axes.

Duplicate JSON mappings, matrix dimension mismatches, missing CIF mappings, and
JSON tokens without matching Cα residues are rejected by the shared parser.
This command additionally rejects nonfinite coordinates, pLDDT outside [0,100],
and nonfinite or negative PAE. It is intended for ordinary AlphaFold protein
chains with unique integer residue IDs, not structures needing insertion-code
or alternate-conformer analysis. Biopython's selected atom is used for
alternate locations. Nonprotein tokens without a matching Cα are unsupported.

**Not used:** other atom coordinates, side-chain chemistry, hydrogen bonds,
sequence annotations, atom-level confidence arrays in JSON, ipTM, pTM, ranking
scores, MSAs, templates, or summary-confidence files. Only the fields listed
above contribute to this analysis. A low score is not a binding-energy estimate.

ZIPs must contain unique `*_model_N.cif` files and exactly one matching
`*_full_data_N.json` per model, with the same member path prefix. Every model
is analyzed independently. Patch IDs are local to a source/model/chain pair;
matching patch IDs across models does **not** imply matching contact sites.
Automatic cross-model patch matching is not performed in version 1.0.

## Graph construction

For chain pair A,B and cutoff d:

```
C = {(a,b): a in A, b in B, EuclideanDistance(Cα(a), Cα(b)) <= d}
V = all endpoints of contacts in C
E = {(u,v): u,v in V, u != v, EuclideanDistance(Cα(u), Cα(v)) <= d}
```

Edges E include both same-chain and interchain edges. Only C contributes to
pDockQ2. Same-chain edges provide spatial connectivity, never additional scored
contacts. No noninterface residues are added to connect patches through the
protein interior. Sequence adjacency is irrelevant: residues B:60 and B:136
can be in the same component if joined spatially, directly or transitively.

Connected components of (V,E) are initial patches. No maximum diameter is
imposed: a transitive component can extend much farther than d.

## Recursive bridge splitting

For a component P:

1. Enumerate subsets S of P by size 1,2,...,N, with residue identities sorted
   lexicographically by chain ID and numerically by residue number.
2. Remove S and its incident edges; compute components of the induced graph
   on P minus S.
3. Within each component, count **remaining interchain contacts** and their
   distinct endpoints on each chain. Graph nodes that lost all opposite-chain
   neighbors do not count toward the size filter or pLDDT average.
4. Accept the first S yielding at least two components meeting both size
   filters. Retain all qualifying components as children. Record nonqualifying
   fragments as discarded residues on the parent row.
5. Remove S from every child; do not duplicate or assign separator residues to
   either side. They remain represented and scored in the parent.
6. Repeat independently for every child. Stop when no acceptable separator of
   size ≤N exists. Traversal uses an explicit stack, not Python recursion.

Size-filtered initial components are omitted; the full-interface reference row
still includes their contacts. After a split, graph residues that no longer
have an interchain contact remain in that child's connectivity graph. Their
presence is transparent in `graph residues`; they are absent from its actual
interface-residue lists and confidence averages. Connectivity is always the
induced subgraph of the original interface graph, not a graph rebuilt from
newly eligible endpoints after each split.

N=1 finds qualifying articulation points. N>1 finds small vertex separators;
it is not a bridge-length criterion. Separator identity is chosen geometrically
and deterministically, without consulting confidence scores. Only one split
is chosen at each node; alternative split histories are not enumerated.
A larger N permits more splits but cannot bypass minimum-size constraints.

For N>1 the search can be combinatorial. With v vertices, up to
sum(k=1..N) binomial(v,k) subsets may be checked at a node, each requiring a
connectivity traversal and size assessment. The default N=1 implementation
also explicitly checks each residue rather than using a specialized linear-time
articulation algorithm. The search budget raises an error before output is
written; it never silently truncates a result or labels it complete.

### Worked topology example

Suppose two spatial contact groups are linked only through A:42. Each group
has ≥3 contacting residues on each chain and ≥5 contacts. With N=1:

- Parent `p1` reports both groups and A:42, with status `split`.
- `separator residues` is `A:42`.
- Children `p2` and `p3` exclude A:42 and are scored independently.
- If either group fails the filters after removal, this separator is rejected.

If A:42 and B:90 provide two independent connecting routes, removing just one
may leave the graph connected; N=2 can test removing both. These examples define
graph behavior, not a claim that either separator is a biological boundary.

## Regional pDockQ2 equations

For a reported residue set P, take all original contacts C with both endpoints
in P. Let I_A and I_B be the distinct endpoints of these remaining contacts.

```
w_AB = mean over (a,b) in C_P of 1 / (1 + (PAE[a,b] / 10)^2)
w_BA = mean over (a,b) in C_P of 1 / (1 + (PAE[b,a] / 10)^2)
p_A  = mean Cα pLDDT over unique residues in I_A
p_B  = mean Cα pLDDT over unique residues in I_B
x_AB = w_AB * p_A
x_BA = w_BA * p_B
q(x) = L / (1 + exp(-k * (x-x0))) + b
L=1.31034849; x0=84.7326239; k=0.0747157696; b=0.00501886443
```

Report q(x_AB), q(x_BA), and their maximum. PAE is averaged over contact
**pairs**; pLDDT is averaged over unique interface **residues**. Distances select
pairs but do not continuously weight them. There is no hard PAE threshold,
patch-size bonus, score clipping, or refit. Following the existing scorer, zero
PAE weight or zero mean pLDDT returns score zero rather than the sigmoid floor;
an empty interface therefore returns zero. PAE=0 itself gives weight 1, not 0.

A full-interface reference is always reported for every selected chain pair,
even when no component qualifies. It includes all interchain contacts before
splitting. It is not assigned a patch rank.

## TSV column specification

| Columns | Meaning |
|---|---|
| `source`, `model` | Absolute input CIF/ZIP path and model number (blank if CIF filename has no model number) |
| `scope` | `full` reference or `patch` |
| `chain 1`, `chain 2` | Canonically ordered chain IDs |
| `region 1`, `region 2` | Compact ranges of selected graph nodes on each chain; gaps remain gaps |
| `contact count` | Remaining interchain residue-pair count |
| `interface residues 1`, `interface residues 2` | Exact sorted residue numbers contributing contacts |
| `interface pLDDT 1`, `interface pLDDT 2` | Means over those unique contacting residues |
| `normalized PAE 1 to 2`, `normalized PAE 2 to 1` | Means of PAE weights, not mean PAE in Å |
| `pDockQ2 1 to 2`, `pDockQ2 2 to 1`, `pDockQ2 max` | Directional scores and their maximum within this row |
| `algorithm version` | Specification/algorithm version |
| `patch id`, `parent id`, `depth` | Preorder IDs p1,p2,...; root has blank parent and depth 0; reference ID is `full` |
| `status` | `reference`, `split`, or `leaf` |
| `rank` | Across all qualifying parent/child patches in this model/chain pair; blank for reference |
| `graph residues` | Chain-qualified selected graph vertices, including nodes no longer in interchain contacts |
| `separator residues` | Nodes removed to produce this row's children; blank for leaves |
| `discarded residues` | Nonqualifying fragments removed during this row's split, excluding the separator |
| `residue count 1`, `residue count 2` | Counts of actual contacting residues |
| `distance cutoff`, `max bridge residues`, `min residues per chain`, `min contacts`, `max separator checks` | Invocation parameters repeated for reproducibility |

Ranking uses decreasing maximum pDockQ2, then decreasing contact count, then
lexical `graph residues`. TSV rows retain tree order, not rank order. Scores
are ranked before formatting to six decimals, so apparently tied printed
scores can have different ranks. Parent and child rows overlap deliberately;
never sum their contact counts as independent sites.

## Validation and limits

Automated tests cover transitivity, sequence-discontinuous patches, separate
patches with consecutive residue numbers, inclusive cutoff boundaries,
recursive splitting, one- and two-vertex separators, dangling branches,
minimum-size rejection, deterministic order, directionality/index preservation
through equality with explicit regional scoring, search-budget failure,
CIF/ZIP/directory input, invalid PAE/options, and output schema.

The algorithm's search space is connected components and one deterministic
separator hierarchy, not all residue subsets or all possible binding-site
boundaries. A broad interface lacking a small separator will remain one patch
even if a small subregion has higher pDockQ2. Small qualifying patches can rank
highly. Record counts and inspect parent/child results rather than treating the
maximum as a calibrated probability of interaction. Multiple predicted models
are not independent experimental evidence. Compare residue sets across models
for reproducibility and known structures or contact-directed experiments for
biological accuracy.

## References

- Zhu et al., *Evaluation of AlphaFold-Multimer prediction on multi-chain
  protein complexes*, Bioinformatics (2023),
  [doi:10.1093/bioinformatics/btad424](https://doi.org/10.1093/bioinformatics/btad424).
  Source for the pDockQ2 concept and fitted sigmoid, not this graph algorithm.
- Local scoring implementation: `sieve/alphafold_pdockq2.py`, functions
  `parse_model`, `score_selections`, and `_pdockq2`.
- Local detection implementation: `sieve/alphafold_pdockq2_patches.py`.
- Tests: `tests/sieve/test_alphafold_pdockq2_patches.py`.

## Human-readable standard-output summary

After analysis (and writing the TSV if requested), the command prints up to five highest-scoring
patch rows per input source, ranked across its models and chain pairs involving
the **last chain encountered in the first CIF structural model**. This is CIF
chain order, not lexical chain-ID order. For example, with CIF chains A,B,C,
A–C and B–C are eligible; A–B is not. If chain order varies by model, each model's
own last chain determines eligibility. `--chains` still restricts the analysis;
if it excludes the last chain there are no eligible summary rows.

All qualifying parent, intermediate, and leaf patches are eligible. The `full`
reference is excluded from this summary only; the TSV remains unchanged. The
summary greedily filters overlapping or repeated sites across models. It
shows model, patch ID/status, chain pair, both selected regions, and maximum
directional pDockQ2. Fewer than five matches are all printed; no matches produce
an explicit message. Ranking uses descending score, descending contact count,
then model identifier as text, chain IDs, graph residue text, and patch ID for
stable ties. A directory gets a separate summary for every ZIP.

### Residue-overlap filtering

Summary candidates are visited in score order. Retain a candidate only when
its overlap with every previously retained patch is at most 10%, stopping at
five retained candidates or exhaustion. Overlap is
`|R1 intersection R2| / min(|R1|, |R2|)`, where R contains the actual contacting
(chain ID, residue number) pairs from both chains. Thus overlap cannot exceed
10% of either patch, exact 10% is allowed, and containment counts as 100%.
Graph-only residues are excluded. This filtering compares models within the
same input source, assuming consistent chain and residue identities across
models. Separate ZIPs are summarized independently. The TSV retains all rows.
This is greedy ranking, not optimization of the combined score of five patches;
a high-scoring parent can suppress its descendants. Fewer than five distinct
patches is a valid result.

## Contact-distance statistics

Each TSV row includes `contact distance min`, `contact distance max`,
`contact distance median`, and `contact distance mean`, all in Å. The stdout
summary prints the same statistics to three decimal places (TSV: six).
These describe exactly the **interchain Cα–Cα residue pairs included in that
row's pDockQ2 score**, with each pair counted once and equally weighted.
Same-chain connectivity edges and pairs exceeding `--distance-cutoff` are
excluded. Statistics are recalculated independently for each parent and child.
Median is the middle sorted distance, or the average of the two middle
distances for an even number of contacts. Mean is the arithmetic average
(labeled “average” in stdout). For one contact all four statistics equal its distance. Rows without contacts have blank statistics, not artificial zero distances.

These describe a distribution truncated at the chosen cutoff. Lowering the
cutoff can change graph connectivity, splits, qualifying patch sizes, and
scores; it does not guarantee a higher score or a particular motif size.

### Nearest contact per residue

The summary additionally reports **Nearest contact per last-chain residue**,
with the residue count, minimum, maximum, median, and arithmetic average in Å.
Each actual contacting residue on that chain contributes exactly one value:
its shortest Cα–Cα distance to a residue on the partner chain **within this
reported patch**. A residue with five partners contributes only its minimum,
not five observations. Noncontacting graph nodes contribute no observation.
The partner chain is the other chain of the current pair, not all other chains
combined. Child patches use their own remaining partners after bridge removal.

The TSV reports both directions using columns
`nearest contact distance min 1`, `nearest contact distance max 1`,
`nearest contact distance median 1`, `nearest contact distance mean 1`,
and the corresponding columns ending in `2`. Side 1/2 means `chain 1`/`chain 2`,
not CIF order. Stdout selects whichever side is the last CIF chain. The existing
`contact distance ...` columns are unchanged and stdout labels them **All
interchain contact pairs**. Both use blanks for zero observations; for one
residue all four nearest-contact statistics equal its minimum distance.

For example, if one B residue has contact distances 1 and 3 Å and another has
one contact at 7 Å, all-pair distances are [1,3,7] (mean 3.667 Å), while nearest
contacts per B residue are [1,7] (mean and median 4 Å).

For cutoff selection, nearest-contact statistics describe **residue coverage**,
not docking accuracy. They do not alter graph construction, pDockQ2, ranking,
or overlap filtering. To compare coverage across cutoffs, hold a starting
residue set and its partner set fixed (for example, the patch from an 8 Å run)
and count how many per-residue minima are ≤ each candidate cutoff. This command
still reports each run's detected patches; it does not freeze a baseline across
runs or run a cutoff sweep. Recomputed patches can lose residues and partners,
so their medians are not directly comparable as fixed-set coverage measurements.

### Full-interface results before the patch summary

For each input source, stdout first reports every analyzed chain pair in every
model, before any patch splitting. Each entry includes both directional pDockQ2
scores and their maximum on the first line, and contact-pair and contacting-
residue counts on the second. Distance statistics are omitted from this full-
interface stdout section; they remain in the TSV and top-five patch summary.
Models are ordered numerically, then chain pairs lexically. These are pairwise full-interface scores, not a single
aggregate score for the entire multimer and not chain-versus-all-other-chains
scores. By default all chain pairs are analyzed; `--chains` still restricts
which pairs are computed. The subsequent top-five summary retains its
last-chain and ≤10% overlap filters. TSV contents and scoring are unchanged.

### PAE-only confidence in stdout

Both full-interface and patch summaries now show the two directional PAE-only
confidence values, using the existing TSV `normalized PAE 1 to 2` and
`normalized PAE 2 to 1` columns without changing their calculations. Each is
`mean(1 / (1 + (PAE / 10 Å)^2))` over the row's interchain contact pairs. The
transformation is applied to each pair before averaging; it is not a transform
of the mean PAE. It includes neither pLDDT nor the fitted pDockQ2 sigmoid.
Values are dimensionless, with higher values indicating lower predicted error:
PAE 0, 5, 10, and 20 Å map respectively to 1, 0.8, 0.5, and 0.2 for one pair.
This is an inverse-like transformation, not the reciprocal 1/PAE, and not a
calibrated probability of binding. No-contact rows display `n/a` in stdout
(the existing TSV convention of zero normalized PAE is retained).

Full-interface output retains two lines per model/pair, appending PAE-only
confidence to the contact-count line. Patch output uses a separate confidence
line. Ranking and non-overlap selection still use regular pDockQ2, not PAE-only
confidence. The arrow labels retain the matrix-index conventions defined above.
