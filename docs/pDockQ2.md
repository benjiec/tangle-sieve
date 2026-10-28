# pDockQ2: scoring, spatial patches, and interpretation

This document describes the **regional/pairwise pDockQ2 implementation in sieve**,
including the quantities it measures and how we use it. It uses the published
pDockQ2 transformation and fitted coefficients without recalibration. Applying
that transformation to selected residue regions and automatically detected
patches is our adaptation, not a separately validated scoring method.

Implementation references:

- [Input parsing, scoring, and coefficients](../sieve/alphafold_pdockq2.py):
  `parse_model`, `score_selections`, `_pdockq2`, `PDOCKQ2_PARAMETERS`.
- [Spatial patch detection and reporting](../sieve/alphafold_pdockq2_patches.py).
- [Patch-search command](../scripts/alphafold-pdockq2-patches.py).
- [Scoring and graph tests](../tests/sieve/test_alphafold_pdockq2_patches.py).

Patch algorithm version: **1.0**. This document is the consolidated scoring,
command-line, patch-detection, and output specification.

## What goes into the calculation

| Input | Source | Role |
|---|---|---|
| Cα coordinates | First structural model in CIF | Select physical residue-pair contacts |
| Residue pLDDT | CIF Cα B-factor field | Estimate local structural confidence at contacting residues |
| Directional PAE | Full-data JSON `pae` matrix | Estimate uncertainty in relative residue positioning |
| Chain/residue-to-token mapping | JSON `token_chain_ids`, `token_res_ids` | Look up the correct original PAE indices |
| Residue selection and distance cutoff | User region or detected patch; cutoff defaults to 8 Å | Define the interface being assessed |

The CIF B-factor is interpreted as AlphaFold pLDDT, not a crystallographic
thermal-motion parameter. pLDDT is used on its original 0–100 scale. Coordinates
and PAE are in Å. Cropping a region preserves original PAE token indices.

PAE and pLDDT describe different aspects of confidence: pLDDT describes local
structural accuracy; PAE describes positional uncertainty for one residue when
aligned on another. Neither is a physical inter-residue distance. See the
[AlphaFold confidence documentation](https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/).

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
| `--silent` | off | Suppress stdout; retain TSV output and stderr diagnostics |
| `--dd` | off | Enable provisional six-helix DD face candidates |
| `--dd-helix-presence-threshold` | 3 | Minimum contacting residue count per required helix; positive integer |
| `--dd-loop-presence-threshold` | 1 | Minimum contacting residue count per required loop; positive integer |
| `--chains A B` | all pairs | Analyze one distinct pair; summarize sites on the second argument; output chain order is lexical |
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

## 1. Define contacts using physical distance

For selected residue sets A and B, let

\[
C=\{(a,b): a\in A,\ b\in B,\ \|\mathbf r_{C\alpha,a}-\mathbf r_{C\alpha,b}\|\le d_c\}.
\]

The default cutoff is **d_c = 8 Å**, inclusive. Let I_A and I_B be the unique
residues on each side participating in at least one pair in C.

Physical distance controls **inclusion**, not continuous weighting. At an 8 Å
cutoff, contacts at 5 Å and 7 Å have the same distance-based inclusion status;
a pair at 9 Å is excluded. Each included pair contributes once to PAE averaging.
There is no additional contact-count bonus or reward for shorter distances.

Changing the cutoff can change the selected contact pairs, interface residues,
PAE averages, pLDDT averages, and patch topology. Increasing the cutoff therefore
does not guarantee a larger score; decreasing it does not guarantee improvement
or a specific motif size.

## 2. Transform and average PAE

Define the dimensionless weight

\[
w(p)=\frac{1}{1+(p/10\text{ Å})^2}.
\]

The two directional PAE-only confidence values are

\[
W_{AB}=\frac{1}{|C|}\sum_{(a,b)\in C}w(PAE[a,b]),\qquad
W_{BA}=\frac{1}{|C|}\sum_{(a,b)\in C}w(PAE[b,a]).
\]

The transform is applied **before averaging**; this is not w(mean PAE).
The divisor 10 Å is the PAE scaling parameter, separate from the physical
contact cutoff of 8 Å.

| PAE for one pair | PAE weight |
|---:|---:|
| 0 Å | 1.00 |
| 5 Å | 0.80 |
| 10 Å | 0.50 |
| 20 Å | 0.20 |
| 30 Å | 0.10 |

Higher W means lower predicted error. W is not a binding probability, and 0.5
must not be described as “50% probability of binding.” It is an inverse-like
transformation, not 1/PAE. There is no hard PAE exclusion threshold in this
implementation. All physically selected contacts contribute.

The output labels A→B and B→A identify matrix lookup directions: A→B uses
`pae[token_A][token_B]`. For AlphaFold Server, matrix element [i,j] describes
error in token j when aligned using token i's frame. These arrows are not
causal interaction directions. See the
[AlphaFold Server output definition](https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/alphafold-server-your-gateway-to-alphafold-3/interpreting-results-from-alphafold-server/).

The TSV fields `normalized PAE 1 to 2` and `normalized PAE 2 to 1` contain W,
not PAE in Å. Stdout calls them **PAE-only confidence**.

## 3. Average pLDDT over unique interface residues

\[
P_A=\frac{1}{|I_A|}\sum_{a\in I_A}pLDDT(a),\qquad
P_B=\frac{1}{|I_B|}\sum_{b\in I_B}pLDDT(b).
\]

Each contacting residue contributes once, even if it has several partners.
Thus PAE is averaged over **pairs**, whereas pLDDT is averaged over **unique
residues**. Low-pLDDT residues outside the selected interface do not directly
lower its score. Low pLDDT within the interface does lower the corresponding
score, holding the PAE term constant.

Our directional input values are

\[
X_{AB}=W_{AB}P_A,\qquad X_{BA}=W_{BA}P_B.
\]

This explicitly describes our pairing of matrix direction and chain pLDDT;
the two directional inputs should be retained when comparing results rather
than interpreting only their maximum.

## 4. Apply the fitted sigmoid

\[
q(X)=\frac{L}{1+\exp[-k(X-X_0)]}+b,
\]

with the constants used in the code:

| Parameter | Value | Meaning |
|---|---:|---|
| L | 1.31034849 | Vertical amplitude |
| X_0 | 84.7326239 | Horizontal midpoint |
| k | 0.0747157696 | Steepness |
| b | 0.00501886443 | Lower offset |

Report q(X_AB), q(X_BA), and their maximum. The maximum is over the **two
directions for that row**, not automatically over models. The code does not
clip the fitted result to [0,1].

Special handling inherited from our scorer: if mean pLDDT or mean PAE weight
is zero, return zero instead of applying the sigmoid. An interface with no
contacts has zero score. The TSV uses zero PAE weight for no contacts, but
stdout displays PAE-only confidence as `n/a`, distinguishing missing contacts
from measured low confidence. PAE equal to zero gives weight one, not zero.

### Worked example

If W = 0.80 and interface pLDDT = 85, X = 68 and q(X) is approximately 0.297.
A PAE weight near 0.8 therefore need not give a pDockQ2 near 0.8.

Conversely, q = 0.272623 implies X approximately 66.529. This could arise from
W = 0.805586 and pLDDT approximately 82.585 in the winning direction. Inferring
pLDDT this way is conditional on knowing which direction supplied the maximum;
it does not establish both chain-specific pLDDT values.

This is not the standard logistic function centered at zero. Its midpoint is
84.733, where the output is L/2+b, approximately 0.660. A sharp band in a PAE
heatmap indicates contrast with neighboring regions, not necessarily the low
absolute error and high interface pLDDT needed for a large sigmoid output.

## Origin of the fit and scope of this adaptation

Zhu et al. fitted the pDockQ2 sigmoid to actual interface DockQ_i values using
SciPy. Their benchmark comprised 1,928 complexes with two to six chains,
including homo- and heteromers, modeled using AlphaFold-Multimer v2.2.0. The
analysis used its top-ranked predictions. DockQ_i assessed one chain against
all remaining chains together using experimental reference structures. The
fitted coefficients, including X_0=84.733, are empirical parameters—not physical
constants or per-model quantities. See Methods §§2.1–2.6 in the
[original publication](https://doi.org/10.1093/bioinformatics/btad424).

Our implementation applies the published form and coefficients to chain pairs,
explicit regions, or graph-detected patches, using the contact and averaging
rules specified above. We have not established reference-code equivalence for
all implementation conventions, nor recalibrated the score for selected
six-residue motifs, long disordered tails, or AlphaFold Server predictions.
The file named `scripts/pDockQ2.py` is a separate legacy pDockQ-style script:
it uses interface size and pLDDT without PAE. It should not be confused with
the implementation described here.

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

## TSV column specification

Only final **leaf patches** are exported and considered for the top-five summary.
Full-interface references and intermediate split parents remain internal to the
algorithm. A single unsplit component can coincide with the full interface, but
multiple or split patches differ. No qualifying leaves produces a header-only TSV.

| Columns | Meaning |
|---|---|
| `sequence_id`, `fold_id` | Shared composition and seed-aware identity; `NA` if unavailable |
| `source`, `model` | Input locator/path and model number |
| `scope`, `status` | Always `patch`, `leaf` in the export |
| `chain 1`, `chain 2` | Canonically ordered chain IDs |
| `region 1`, `region 2` | Selected graph-node ranges; gaps remain gaps |
| `interface structure 1`, `interface structure 2` | DSSP elements and contacting residue ranges, as in stdout; missing assignments labeled `unassigned` |
| `patch 1 residues`, `patch 2 residues` | Unique contacting residue counts |
| `contact count` | Interchain Cα contact pairs |
| `interface pLDDT 1/2`, `normalized PAE 1 to 2/2 to 1`, `pDockQ2 1 to 2/2 to 1`, `pDockQ2 max` | Existing confidence statistics |
| `patch id`, `rank` | Existing preorder patch ID; rank among exported leaves within model/chain pair |
| `interchain CA distance min/max/median/mean` | Interchain Cα–Cα contact-pair distances, Å |
| `nearest interchain CA distance min/max/median/mean 1/2` | Per-residue nearest Cα-contact distributions, Å |
| `total buried surface area` | Sum of selected-residue burial on both chains, Å²; not divided by two |
| `total C/S contribution percent` | 100 × summed C/S burial / summed total burial; `NA` if denominator is zero |
| `nonpolar contacts`, `salt bridge candidates`, `steric overlap candidates` | Unique residue-pair counts, with both endpoints in the patch |
| `phosphate/basic contacts` | Separate unique residue-pair count for phosphate/basic geometric candidates |
| `salt bridge candidate heavy-atom distance min` | Minimum distance among within-patch ordinary salt candidates, Å; `NA` if absent |
| `salt bridge residues` | Semicolon-separated chain/residue/atom pairs, using the nearest atom pair per residue pair, with distance; `NA` if absent |
| `nominal patch charge difference` | Side 1 minus side 2; numeric for ordinary residues, labeled −1/−2 phosphate scenarios when needed, `NA` for unsupported chemistry |
| `patch 1 hydropathy mean`, `patch 2 hydropathy mean` | Means over unique contacting residues; `NA` for sets containing PTMs |
| `candidate DD family face 1/2` | Ranked heuristic labels with `--dd`; `NA` when disabled/ineligible; `none` when eligible but no face qualifies |
| `dd enabled`, `dd helix presence threshold`, `dd loop presence threshold` | DD settings for this run |
| `confidence origin` | Inherited confidence note for assembled models; otherwise `NA` |
| `distance cutoff`, `max bridge residues`, `min residues per chain`, `min contacts`, `max separator checks` | Patch-detection parameters |

Algorithm version, parent ID, depth, graph/separator/discarded residue lists,
and raw interface-residue columns are omitted from the export. Descriptor
calculations are shared between stdout and TSV. The stdout salt-candidate list
also includes nearby candidates with only one endpoint in the patch, explicitly
labeled; those contextual candidates do not contribute to TSV counts or distances.
Phosphate/basic contacts are not counted as ordinary Asp/Glu salt candidates.

Leaf rank uses decreasing pDockQ2, then contact count and deterministic residue
ordering. TSV rows retain traversal order, not rank order. Values are formatted
to six decimal places after ranking. Changing settings always recomputes results
and replaces `--output`; there is no cache.

### Silent output

`--silent` suppresses stdout, including model and structure summaries. `--output`
still writes the enriched TSV. With no output path, no file is created. Errors
and library warnings remain on stderr.

```sh
sieve-py ~/git/tangle/sieve/scripts/alphafold-pdockq2-patches.py locator \
  --silent --output patches.tsv --dd --dd-helix-presence-threshold 3 --dd-loop-presence-threshold 1
```

### Candidate DD faces

`--dd` enables a **screening heuristic**, not a validated DD topology classifier.
Exactly six DSSP `H` segments are required on a chain, provisionally numbered
H1–H6 in sequence order. Five/seven helices are not repaired or merged. Missing
assignments inside that span make the chain ineligible. G/I helices are not
counted as H segments. A required loop is the complete intervening sequence
between its two H segments and must contain only DSSP `-`, `T`, or `S` codes;
empty loops or structured insertions cannot satisfy a loop requirement.

| Candidate | Required elements |
|---|---|
| Ia | H1, H4 |
| Ib | H2, H3 |
| IIa | H4, H4–H5 loop |
| IIb | H5–H6 loop, H6 |
| IIIa | H3 |
| IIIb | H1–H2 loop, H3–H4 loop |

Each required helix must contain at least `--dd-helix-presence-threshold`
contacting patch residues (default 3). Each required loop must contain at least
`--dd-loop-presence-threshold` contacting residues (default 1). Both are positive
integer counts, not fractions. A full match requires every required element to
meet its threshold. A label followed by `(partial)` requires at least one
contacting residue in every required element, at least one element meeting its
threshold, and at least one falling below it. If any element has zero contacts,
or none meets its threshold, the face is omitted. Elements shorter than their
threshold can contribute to partial matches but cannot pass that threshold.
The single-element IIIa face cannot be partial under this rule.

Full matches precede partial matches. Within each group, labels are ordered by minimum required-element
fractional coverage (contacting residues divided by element length),
then mean coverage; remaining ties use Ia, Ib, IIa, IIb, IIIa, IIIb order.
This ordering is not a likelihood or calibrated probability.
Partial labels expose threshold edge cases; classifications still depend on
the contact distance cutoff.

Element identities follow the [PIDDosome structural study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2908332/)
and [DD assembly description](https://pmc.ncbi.nlm.nih.gov/articles/PMC3320699/).
Using whole helices/loops, sequential H labels, and these residue-count cutoffs are explicit
implementation approximations. Six helices do not establish canonical topology
or the orientation of the contacted helix surface. Chains containing multiple
domains or fragmented/fused helices may fail or be misclassified. Validate labels
against known structural face assignments and check stability under alternative helix and loop count thresholds
before using them as evidence. No reference-based classification is performed.

## Human-readable standard-output summary

After analysis (and writing the TSV if requested), the command prints fold identity and model confidence followed by up to five highest-scoring
patch rows per input source, ranked across its models and chain pairs involving
the **last chain encountered in the first CIF structural model**. This is CIF
chain order, not lexical chain-ID order. For example, with CIF chains A,B,C,
A–C and B–C are eligible; A–B is not. If chain order varies by model, each model's
own last chain determines eligibility. `--chains A B` restricts analysis and summary to A–B and uses B for overlap
filtering, even if B is not the last CIF chain. Reversing the arguments to
`--chains B A` uses A for overlap filtering; output pair order remains lexical.

All qualifying parent, intermediate, and leaf patches are eligible. The `full`
reference is excluded from this summary only; the TSV remains unchanged. The
summary greedily filters overlapping or repeated sites across models. It
shows model, patch ID/status, chain pair, both selected regions, and maximum
directional pDockQ2. Fewer than five matches are all printed; no matches produce
an explicit message. Ranking uses descending score, descending contact count,
then model identifier as text, chain IDs, graph residue text, and patch ID for
stable ties. A directory gets a separate summary for every ZIP.

### Target-chain residue-overlap filtering

Summary candidates are visited in score order. Retain a candidate only when
its overlap with every previously retained patch is at most 10%, stopping at
five retained candidates or exhaustion. Overlap is
`|R1 intersection R2| / min(|R1|, |R2|)`, where R contains the actual contacting
(chain ID, residue number) pairs from the **summary target chain only** (the second `--chains` argument, or the last CIF chain by default). Thus overlap cannot exceed
10% of either patch, exact 10% is allowed, and containment counts as 100%.
Graph-only residues are excluded. This filtering compares models within the
same input source, assuming consistent chain and residue identities across
models. Separate ZIPs are summarized independently. The TSV retains all rows.
This is greedy ranking, not optimization of the combined score of five patches;
a high-scoring parent can suppress its descendants. Shared partner-chain
residues do not suppress distinct target-chain sites: two B sites contacting
the same A surface remain eligible. This identifies sites on the target
chain, not distinct partner-specific interactions; overlapping target sites
can suppress each other even when their partner chains differ. Fewer than five distinct
patches is a valid result.

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
selected target-chain and ≤10% target-chain residue-overlap filters. TSV contents and scoring are unchanged.

### PAE-only confidence in stdout

Both sections display the directional PAE-only weights defined in the scoring
formula above. Full-interface output appends them to the contact-count line;
patch output uses a separate line. No-contact entries print `n/a` (the TSV
retains zero). Ranking still uses regular pDockQ2 rather than PAE-only confidence.

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

## Contact-distance statistics

Each TSV row includes `interchain CA distance min`, `interchain CA distance max`,
`interchain CA distance median`, and `interchain CA distance mean`, all in Å. The stdout
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

The summary additionally reports **Nearest contact per target-chain residue**,
with the residue count, minimum, maximum, median, and arithmetic average in Å.
Each actual contacting residue on that chain contributes exactly one value:
its shortest Cα–Cα distance to a residue on the partner chain **within this
reported patch**. A residue with five partners contributes only its minimum,
not five observations. Noncontacting graph nodes contribute no observation.
The partner chain is the other chain of the current pair, not all other chains
combined. Child patches use their own remaining partners after bridge removal.

The TSV reports both directions using columns
`nearest interchain CA distance min 1`, `nearest interchain CA distance max 1`,
`nearest interchain CA distance median 1`, `nearest interchain CA distance mean 1`,
and the corresponding columns ending in `2`. Side 1/2 means `chain 1`/`chain 2`,
not CIF order. Stdout selects whichever side is the last CIF chain. The existing
`interchain CA distance ...` columns and stdout labels them **All
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

## Intended utility and interpretation

Our chosen use is **ranking and prioritizing candidate interfaces for further
structural assessment**, without recalibration or claims of binding probability.
Use the same modeling settings, contact cutoff, region/patch rules, and score
selection when comparing domain versions or homologous systems. Report PAE-only
confidence and interface pLDDT alongside pDockQ2 to expose why a score is low.

A larger regional score means a larger combined confidence signal under these
rules. It does not by itself establish stronger affinity, preference between
competing binding sites, a chemically favorable interaction, or experimental
binding. Nor do pDockQ2 predictions automatically inherit the interpretation
of thresholds on actual reference-based DockQ. All ranking/filtering thresholds
chosen for these regional analyses should be disclosed as analysis choices.

A proposed follow-up assessment is to examine:

1. Contact location relative to the candidate motif and domain surface.
2. Similarity of the predicted pose and contact map to an experimentally
   established homologous interface.
3. Chemically plausible side-chain packing, hydrogen-bond geometry, salt
   bridges, hydrophobic contacts, and possible steric clashes.
4. Recurrence of the site and pose across predicted models and sensitivity to
   analysis settings.

These checks support or challenge a structural hypothesis; predicted coordinates
do not prove chemical bonds or biological binding. Model recurrence is a
reproducibility check, not independent experimental replication. Targeted
mutations and binding measurements can test functional relevance.

For a known reference organism, generate a comparable prediction to obtain PAE
and pLDDT, and use its experimental structure for geometric comparison.
Experimental PDB coordinates alone do not provide AlphaFold confidence values.
Similar scores between organisms are supporting context, not calibration or
proof of a conserved interaction. Incorrect-site or other relevant controls
help test whether the approach distinguishes the proposed site from alternatives.

### Ranking versus calibration

With positive L and k, the sigmoid is monotonic: changing its coefficients alone
preserves rankings for nonempty, positive-input scores. Recalibration could
change the mapping to expected experimental interface accuracy, but cannot
repair a reversed ordering of the combined input X. Testing PAE and pLDDT as
separate predictors would be a different model requiring validation.

The main unresolved limitation is whether these regional scores reliably rank
the specific class of interfaces being studied. Compare predictions with known
structures, including incorrect predictions, and evaluate held-out protein
families to test that. Binding preference requires comparative binding evidence,
not just structural-confidence calibration.

## Suggested reporting language

> Predicted interfaces were prioritized using a regional adaptation of pDockQ2.
> Interface residues were clustered as connected components of a graph with
> Cα–Cα edges ≤8 Å, including same-chain edges. Components were recursively
> split using vertex separators of at most one residue, requiring at least
> two children with three contacting residues per chain and five interchain
> contacts each. The first qualifying separator in size and residue-ID order
> was used; separator residues remained in the parent and were excluded from
> children. Parent and descendant patches were retained.
> Interchain contacts were defined by Cα–Cα distances ≤8 Å. Directional mean
> PAE weights were multiplied by mean pLDDT over the corresponding unique
> interface residues and transformed with the published pDockQ2 sigmoid.
> Both directional scores, their maximum, PAE-only confidence, and contact
> geometry were evaluated. The coefficients were not recalibrated for regional
> patches, and the scores were used for ranking and assessment rather than
> interpreted as binding probabilities. Candidate interfaces were further
> examined for agreement with experimentally established interfaces and
> chemically plausible interactions.

Adapt this paragraph to the actual cutoff, bridge/size filters, selection method, and analyses
performed; do not claim validation steps that were not carried out.

## Additional interface descriptors

The [interface descriptor definitions](interface-comparison.md) document DSSP,
hydropathy/charge context, geometric interaction candidates, and buried surface
area used by the patch stdout summary. These calculations do not change pDockQ2.

### Structural descriptors in the patch stdout summary

`alphafold-pdockq2-patches.py` now includes structural descriptors by default.
The usual CIF plus `--full-data`, ZIP, and directory inputs still work. No second
reference structure, alignment, output archive, or report directory is needed.
The optional `--output` TSV exports enriched leaf-patch rows; see the column specification above.

Before the full-interface scores, the summary groups chains by exact modeled
amino-acid sequence within each input. Each group lists its chain IDs and length,
then N-to-C secondary-structure strings and DSSP-code composition percentages.
Identical sequence **does not imply identical structure**: distinct assignments
or residue numbering are displayed separately with their model and chain IDs.
Identical strings/compositions share a line. Counts use all modeled Cα residues;
missing DSSP assignments are explicitly `?`. Strand labels describe runs, not
whole sheets. Secondary-structure labels are local to each model and chain.

Each of the existing top-five distinct sites now also reports:

- Mean local pLDDT over unique contacting residues on each chain.
- Contacting-residue sequence (concatenated in residue order) and each residue
  range's DSSP element/code; use the region list to recognize sequence gaps.
- Patch, upstream, downstream, combined sequence-shell, and spatial-shell mean
  Kyte–Doolittle hydropathy, nominal charge, K/R and D/E counts, and histidine count.
- Patch buried area per chain and its carbon/sulfur contribution, in Å².
- Unique residue-pair counts for nonpolar contacts, candidate salt bridges, and
  candidate steric overlaps. Atom pairs are deduplicated within each type.
- Interchain salt-bridge candidates touching the patch, with the shortest charged
  atom pair per residue pair. Labels distinguish both endpoints inside the patch
  from one endpoint outside the Cα-defined patch. This is not a search for
  intrachain salt bridges or arbitrary nearby charged residues.
- The difference between the two patches' nominal net charges (chain 1 minus
  chain 2). The separate signed charges are the primary quantities; their
  difference is not an electrostatic-complementarity or binding score.

Descriptor equations, atom criteria, conventions, and sources are in
[interface-comparison.md](interface-comparison.md). The patch graph, ranking,
last-chain overlap filter, confidence calculations, and distance statistics
are unchanged. The descriptor uses the exact selected patch row, including
custom bridge/minimum-size settings; it does not rerun patch detection with
other defaults. Pairwise atom interactions and surface burial are cached for
reported patches sharing a model/chain pair.

Additional options:

```text
--dssp PATH          mkdssp executable (otherwise SIEVE_DSSP, then PATH)
--no-dssp            Explicitly omit secondary structure; retain other descriptors
--sequence-flank 5   Upstream/downstream residues around each contacting residue
--spatial-radius 8   Same-chain heavy-atom neighborhood radius in Å
--sasa-points 240    Surface sample points per atom
```

Set `export SIEVE_DSSP="/path/to/dssp/bin/mkdssp"` to configure the executable.
`--dssp` overrides this variable; `--no-dssp` bypasses discovery. An invalid or
empty explicit setting is an error. There is no repository-local fallback.
See [executable configuration](../README.md#external-executable-configuration).

If DSSP cannot be located, secondary structure is explicitly reported as
unassigned; other descriptors remain available. If a located DSSP fails, the run
fails rather than inventing assignments. A functioning installation requires
its runtime dictionaries, not only the executable. Side chains are not rebuilt;
missing atoms can reduce detected contacts. Hydropathy, nominal charge, and
geometric contact candidates do not establish binding energy or confirmed bonds.

Example, without writing a TSV:

```sh
sieve-py scripts/alphafold-pdockq2-patches.py \
  data/fold_hs_myd88_dd_1x_hs_irak4_dd_1x.zip
```

### Extracted AlphaFold directory input

Pass an extracted fold directory directly:

```sh
sieve-py scripts/alphafold-pdockq2-patches.py /path/to/extracted_fold
```

The script reads immediate `*model_N.cif` files and matches each to
`*full_data_N.json` with the same prefix. Models are ordered numerically and
reported together under the directory's path. Duplicate model numbers or missing
confidence files cause an error. Subdirectories are not searched.

When extracted models are present, ZIP files in that directory are ignored to
avoid analyzing the same fold twice. If no extracted models are present, the
existing directory-of-ZIPs behavior is retained. `--full-data` is only needed
and accepted when providing a single CIF file.

## Manifest fold locators

The patch script accepts the `fold_locator` column from `alphafold-manifest.py`
as its positional input. For example:

```sh
sieve-py ~/git/tangle/sieve/scripts/alphafold-pdockq2-patches.py \
  '/path/batch.zip::inner.zip!/folder/fold'
```

The locator selects exactly one directory within an archive, including nested
ZIPs, and pairs its model CIFs with their matching full-data JSON files. Other
folds and unrelated nested archives are not analyzed. The output `source` retains
the complete locator so folds within the same archive remain distinguishable.
Plain CIF, ZIP, and directory inputs retain their existing behavior. See the
[manifest documentation](alphafold-manifest.md) for locator syntax and migration.

## Model confidence and optional chain structure summary

Every stdout summary prints whole-fold ipTM and pTM once for each analyzed model,
including models without qualifying interfaces or patches. Values come from the
matching `*_summary_confidences_N.json` beside each model, including inside ZIPs
and nested fold locators. Missing files, fields, or null values print
`n/a (unavailable)`; invalid score values or malformed JSON produce an error.
These scores are reported as supplied, not calculated from PAE. For a renamed
standalone CIF, a numbered `*_full_data_N.json` path can identify its adjacent
summary file. The patch TSV includes both IDs; see the column specification above.

The per-chain secondary structure/composition summary prints only with
`--structure`. Patch-level descriptions still print by default. Existing
`--dssp` and `--no-dssp` options continue to control secondary structure assignment.

## Fold header and printed patch results

Each fold's stdout block starts with `Fold ID: HASH` and `Sequence ID: HASH`.
Both use the manifest's shared calculation: the last eight hex characters of
SHA-1 over compact JSON. `sequence_id` hashes the sorted chain sequences with
copies; `fold_id` hashes `(seed, sorted_sequences_with_copies)`. The seed comes
from the matching job-request JSON and is normalized to an integer. Models with
different sequence-composition hashes or seeds are rejected.

For standalone inputs with no job-request file, the fold ID is reported as
unavailable; the sequence ID still prints when sequence metadata is available.
Missing sequence metadata makes both unavailable. A present job request with
invalid or multiple seeds produces an error. The manifest requires a valid seed
and skips folds without one.

The full-interface pDockQ2 section is no longer printed. Per-model ipTM/pTM,
optional `--structure` summaries, and the top distinct patch results remain.
Full-interface rows remain internal and are not exported; only leaf patches are written.

## Read locators or result rows from stdin

Pass `-` as the positional input to read headerless TSV from stdin. Each record
is either a single locator/path (existing direct-input behavior), or a result
record produced by `alphafold-list-results.py` (five columns, or six with
`--seed`). Result rows select exactly their
`model_id` within their `fold_locator`, before model data is read. The shared
`alphafold_result_format` module defines the column order and reader/writer.

```sh
sieve-py ~/git/tangle/sieve/scripts/alphafold-list-results.py manifest.tsv b911edee --highest | \
  sieve-py ~/git/tangle/sieve/scripts/alphafold-pdockq2-patches.py -
```

Multiple records are analyzed separately in input order. Blank records are
ignored; empty input is a successful no-op. Existing options apply to every
record. `--full-data` remains limited to direct CIF inputs, including one-field
CIF records. Errors identify their stdin record and abort before stdout or the
optional TSV is written. Model IDs must exist; a missing model never falls back
to all models. Piped scores/IDs do not replace source metadata.

## Assembled complexes without PAE

Directories and ZIP locators produced by `alphafold-overlap.py` support geometric
patch analysis without a joint PAE matrix. PAE-dependent scores and pDockQ2 print
`NA`; ranking uses contact count and deterministic ties. Known source pLDDT is
labeled inherited. The TSV always includes `confidence origin`, populated for derived rows and `NA` otherwise. Pairwise contacts, burial, chemistry and
secondary structure retain their existing definitions; this does not add combined
multichain patch scoring. See [assembly documentation](alphafold-overlap.md).


### Phosphorylation in patch descriptions

The patch summary supports SEP/TPO/PTR geometry, phosphate/basic contact
candidates, and explicit −1/−2 phosphate charge scenarios. DSSP uses the original
CIF; `--structure` reports assignment coverage for modified residues. Hydropathy
remains unavailable for residue sets containing PTMs. See
[PTM tokens and phosphorylation descriptors](alphafold-manifest.md#ptm-tokens-in-patch-analysis)
for aggregation rules, chemical assumptions, and validation limits.
