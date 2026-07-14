# Rules

Sieve rules are Python objects, usually collected in a `Rules(...)` wrapper and
referenced from scripts with `-r module.path.rule_name`.

For example:

```python
from sieve.rules import Leader, Rules

rule = Rules(
    Leader().upstreamOfPfam("PF00081").betweenAA(-45, -15).is_mTP()
)
```

This document currently covers `Leader` rules. Other rule types can be added
later.


## Leader Rules

`Leader()` rules discover possible N-terminal leader starts, classify those
candidate sequences, and pass proteins whose candidate leader classification
matches the requested prediction.

Supported predictions are:

```python
Leader().is_mTP()
Leader().is_SP()
Leader().is_noTP()
```

`is_mTP()` and `is_SP()` can use either TargetP or DeepLoc. TargetP is the
default. DeepLoc is selected per rule:

```python
Leader().is_mTP(deeploc=True)
Leader().is_SP(deeploc=True)
```

By default, leader starts are considered near the beginning of the protein:

```python
Leader().is_mTP()
```

You can restrict the search window with `betweenAA(start, end)`:

```python
Leader().betweenAA(-30, 3).is_mTP()
```

The default window is equivalent to `betweenAA(-30, 3)`.


## Leader Start Coordinates

Leader starts are methionines (`M`) in a candidate sequence context.

For rules without a Pfam anchor, coordinates are relative to the annotated or
detected protein start:

* `1` means the first amino acid of the protein.
* `2` means the second amino acid of the protein.
* `-1` means one amino acid upstream of the protein start.
* `-30` means thirty amino acids upstream of the protein start.
* `0` is not considered a valid candidate coordinate.

For rules with a Pfam anchor, coordinates are relative to the inferred start of
the Pfam hit:

```python
Leader().upstreamOfPfam("PF00081").betweenAA(-45, -15).is_mTP()
```

This searches for methionine starts from 45 to 15 amino acids upstream of the
inferred `PF00081` anchor. Pfam accessions are prefix-matched, so `PF00081`
matches versioned hits such as `PF00081.28`.

For Pfam-anchored rules, `betweenAA(...)` is an enforced constraint on every
candidate sequence. If Sieve falls back to the original protein sequence
because no suitable methionine start was found, the original sequence is still
kept only when its start is within the requested window relative to the Pfam
anchor.


## Pfam-Anchored Leader Discovery

Use `upstreamOfPfam(accession)` when the leader should be discovered relative
to a domain instead of the protein's annotated start.

```python
Leader().upstreamOfPfam("PF00081").betweenAA(-45, -15).is_mTP()
```

If multiple matching Pfam hits are detected, Sieve uses the earliest inferred
anchor. The anchor is inferred from the query coordinates and, when available,
the target/model start coordinate.

For curated NCBI proteins, the search context follows the protein sequence and
the available genomic annotation. This means that when a domain falls in a later
CDS/exon, leader discovery moves upstream through the spliced protein sequence
rather than counting intronic sequence.

For HMM-detected proteins, Sieve also uses detected protein fragments as a
spliced protein context. For Pfam-anchored HMM-detected proteins, Sieve may also
evaluate an additional raw upstream genomic context at the anchor. Candidate
accessions from this context include an `_anchor` suffix.

If no matching Pfam hit is found, there is no Pfam-relative coordinate to
evaluate, so the rule produces no candidates for that protein.


## Original Sequence Candidates

Sieve sometimes includes the original protein sequence, without any discovered
leader extension or alternative start, as a candidate.

For FASTA-only proteins, the original sequence is always included alongside any
discovered leader candidates for protein-start leader discovery. For
Pfam-anchored leader discovery, the original sequence is included only when its
start is within the requested Pfam-relative `betweenAA(...)` window.

For curated NCBI proteins, the original sequence is always included alongside
any discovered leader candidates for protein-start leader discovery. For
Pfam-anchored leader discovery, the original sequence is included only when its
start is within the requested Pfam-relative `betweenAA(...)` window.

For curated HMM-detected proteins, the original sequence is a fallback
candidate. If leader discovery finds methionine candidates, those candidates are
used without automatically including the original sequence. If no methionine
candidate can be found for a Pfam-anchored rule, Sieve classifies the original
protein sequence only when its start is within the requested Pfam-relative
`betweenAA(...)` window. If the Pfam anchor is missing, the protein is rejected.

Rules that do not use leader discovery use the original protein sequence.


## TargetP Classification

By default, leader candidates are classified with TargetP. Sieve writes one
candidate row for each scoped leader candidate and records the TargetP
probabilities in rule result columns:

```text
Leader.call('noTP')
Leader.call('mTP')
Leader.call('SP')
```

The winning call is the prediction with the highest available probability. A
candidate passes a `Leader().is_mTP()`, `Leader().is_SP()`, or
`Leader().is_noTP()` rule when that winning call matches the requested
prediction.

When a protein has multiple leader candidates, the protein-level leader rule
passes if any scoped candidate has the requested prediction. Candidate-level rows
then show which specific sequence accessions matched the requested leader
classification.


## DeepLoc Classification

DeepLoc classification is enabled by passing `deeploc=True` to `is_mTP()` or
`is_SP()`:

```python
Leader().upstreamOfPfam("PF00081").betweenAA(-45, -15).is_mTP(deeploc=True)
Leader().is_SP(deeploc=True)
```

When a rule uses DeepLoc, the filtering script must be given a DeepLoc result
CSV with `--deeploc-csv`. Sieve does not run DeepLoc itself.

The CSV must include these columns:

```text
Protein_ID
Localizations
Signals
```

`Protein_ID` must match the candidate sequence accession. That includes
discovered leader accessions such as `p1_with_leader_2_M`, not just the
original protein accession.

For signal rules, Sieve reads the `Signals` column:

* `Mitochondrial transit peptide` is treated as `Leader.call('mTP') == 100`.
* `Signal peptide` is treated as `Leader.call('SP') == 100`.
* A blank or unrecognized signal gives both `mTP` and `SP` a score of `0`.

Other numeric DeepLoc columns are copied into rule results as percentage-valued
leader calls. For example, a CSV column named `Endoplasmic reticulum` with a
value of `0.2588` becomes:

```text
Leader.call('Endoplasmic reticulum') = 26
```

Non-numeric values in these additional columns are ignored.

DeepLoc localization can be matched directly with `localize_at(...)`:

```python
Leader().localize_at("Endoplasmic reticulum")
```

`localize_at(...)` uses an exact string match against the `Localizations`
column. It still evaluates scoped leader candidates in the same way as other
`Leader()` rules.

DeepLoc score columns can also be referenced with `LeaderCall(...)` in
positive-rule expressions:

```python
LeaderCall("Endoplasmic reticulum").ge(10)
```


## Candidate Accessions

Discovered leader candidate accessions are derived from the original protein
accession and the candidate start coordinate.

Examples:

```text
p1_with_leader_1_M
p1_with_leader_3_M
p1_with_leader_u2_M
p1_with_leader_u15_PF00081_M
p1_with_leader_u15_PF00081_anchor_M
```

Upstream coordinates use a `u` prefix in accessions, so coordinate `-15`
appears as `u15`.


## Examples

Classify mitochondrial targeting peptides near the protein start:

```python
from sieve.rules import Leader, Rules

rule = Rules(Leader().is_mTP())
```

Classify signal peptides among methionines from amino acid 1 through 10:

```python
from sieve.rules import Leader, Rules

rule = Rules(Leader().betweenAA(1, 10).is_SP())
```

Classify mitochondrial targeting peptides upstream of a Pfam domain:

```python
from sieve.rules import Leader, Rules

rule = Rules(
    Leader().upstreamOfPfam("PF00081").betweenAA(-45, -15).is_mTP()
)
```

Classify mitochondrial targeting peptides with DeepLoc results supplied on the
command line:

```python
from sieve.rules import Leader, Rules

rule = Rules(
    Leader().upstreamOfPfam("PF00081").betweenAA(-45, -15).is_mTP(deeploc=True)
)
```

Match a DeepLoc localization:

```python
from sieve.rules import Leader, Rules

rule = Rules(Leader().localize_at("Endoplasmic reticulum"))
```

Combine leader rules with standard Python operators supported by Sieve rules:

```python
from sieve.rules import Leader, Rules

rule = Rules(
    Leader().is_mTP() | Leader().is_SP()
)
```
