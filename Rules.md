# Rules

Sieve rules are Python objects that describe which proteins should pass a
filter. In normal script usage, you define a rule in a Python module and pass it
with `-r module.path.rule_name`.

```python
from sieve.rules import Leader, Pfam, Rules

mnsod_rule = Rules(
    Pfam.matches("PF00081")
    & Leader().upstreamOfPfam("PF00081").betweenAA(-45, -15).is_mTP()
)
```

The object referenced by `-r` may be either a single rule or an already wrapped
`Rules(...)` instance. Sieve wraps raw rule objects automatically.


## Rule Results

Rule results are written to `rule-results.tsv`. The standard columns are:

```text
protein accession
sequence accession
genome accession
contig accession
pass all
```

Each atomic rule also gets its own column. Some rules add annotation columns,
such as `Leader.call('mTP')`.

`pass all` is the final result for the whole rule expression on a candidate
sequence row. For rules that discover alternative leader starts, a single input
protein can produce multiple output rows, one per scoped sequence candidate.

Common result values are:

```text
true
false
maybe
error
yes
too_far
missing_<value>
missing_<value>_and_<value>
```

`true`, `yes`, and `too_far` are considered passing values. `false`, `maybe`,
`error`, and `missing_*` values are not passing values.

For `&`, any `false` makes the combined rule false. If there is no false value,
`error` and then `maybe` are preserved. For `|`, any passing value makes the
combined rule true; otherwise `maybe` and then `error` are preserved when
present.


## Combining Rules

Rules compose with Python operators:

```python
from sieve.rules import KO, Leader, Pfam, Rules

mnsod_rule = Rules(
    (Pfam.matches("PF00081") | Pfam.matches("PF02777"))
    & KO.matches("K04564")
    & Leader().is_mTP()
)
```

Use parentheses whenever mixing `&` and `|`; normal Python precedence applies.


## Pfam Rules

Use `Pfam.matches(accession)` to require a detected Pfam hit.

```python
from sieve.rules import Pfam, Rules

rule = Rules(Pfam.matches("PF00081"))
```

Pfam accessions are prefix-matched before the first dot. This means
`PF00081` matches versioned hits such as `PF00081.28`. The rule returns `true`
when any detected Pfam row for the protein matches, otherwise `false`.

During artifact evaluation, Pfam rows are regenerated from `input.faa` by
`check-artifacts-by-rules.py`. The required `--pfam-hmm` search uses HMMER
gathering thresholds (`--cut_ga`). Curated Pfam tables are used only while
building curated artifacts when leader discovery is anchored to a Pfam hit.


## KO Rules

Use `KO.matches(accession)` to require a detected KO assignment.

```python
from sieve.rules import KO, Rules

rule = Rules(KO.matches("K04564"))
```

KO accessions are matched exactly. The rule returns `true` when any detected KO
row for the protein matches, otherwise `false`.

During artifact evaluation, KO rows are regenerated from `input.faa` by
`check-artifacts-by-rules.py`. KO searches require `--ko-hmm` and
`--ko-thresholds`; they use the model-specific thresholds instead of
`--cut_ga`.


## HMM Alignment Rules

`HMMAlignment(profile)` rules align each protein to an HMM profile and evaluate
positions in profile coordinates. The profile path is passed to the rule; the
rule label uses the basename of the profile path.

```python
from sieve.rules import HMMAlignment, Rules

rule = Rules(HMMAlignment("/models/profile.hmm").covers(1, 100))
```

The alignment is cached per protein and profile while evaluating a rule set.


## HMMAlignment().is_at(...)

Use `is_at(expected, hmm_pos)` to require an exact amino acid string beginning
at a 1-based HMM position.

```python
rule = Rules(HMMAlignment("profile.hmm").is_at("H", 27))
rule = Rules(HMMAlignment("profile.hmm").is_at("HAD", 27))
```

The rule returns `true` only if every requested HMM position is covered by a
non-gap residue and the observed amino acids match `expected`. If any requested
position is absent or mismatched, the rule returns `false`.


## HMMAlignment().matches_regex(...)

Use `matches_regex(pattern, hmm_pos)` to match a regular expression against the
aligned amino acid string beginning at a 1-based HMM position.

```python
rule = Rules(HMMAlignment("profile.hmm").matches_regex("EFN[AGST]G", 5))
```

The regex is evaluated with Python regular expression syntax. Alignment gaps
inside the HMM-covered region are represented as `-`. Matching stops when the
rule reaches a position outside the HMM alignment. If there is no aligned text
to test, the rule returns `false`.


## HMMAlignment().covers(...)

Use `covers(start, end)` to require every HMM position in an inclusive 1-based
range to be covered by a protein residue.

```python
rule = Rules(HMMAlignment("profile.hmm").covers(1, 45))
```

The rule returns `true` when every position from `start` through `end` has a
non-gap residue, otherwise `false`.

Use `between(protein_start, protein_end)` after `covers(...)` to additionally
require the protein residues mapped to that HMM range to lie within an inclusive
1-based protein-coordinate window:

```python
rule = Rules(HMMAlignment("profile.hmm").covers(2, 50).between(1, 60))
```


## HMMAlignment().spans(...)

Use `spans(start, end)` when the alignment must extend across an inclusive HMM
profile region but internal deletions are allowed. The rule returns `true` when
the first covered HMM position is at or before `start` and the last covered HMM
position is at or after `end`.

```python
rule = Rules(HMMAlignment("profile.hmm").spans(2, 50))
```

Unlike `covers(2, 50)`, `spans(2, 50)` can pass when an internal HMM position,
such as position 37, maps to a gap.

`between(...)` can also constrain the protein coordinates mapped within the
requested HMM span:

```python
rule = Rules(HMMAlignment("profile.hmm").spans(2, 50).between(1, 60))
```

When a rule set performs leader discovery, HMM alignment rules are evaluated
separately for every discovered candidate sequence. Candidates using the same
profile are aligned together in one `hmmalign` batch. Without leader discovery,
the HMM rule is evaluated against the original protein sequence.


## Leader Rules

`Leader()` rules discover possible N-terminal leader starts, classify those
candidate sequences, and pass proteins whose candidate leader classification
matches the requested prediction.

A `Leader()` expression can also be used without a classification method to
perform candidate discovery only:

```python
Leader().upstreamOfPfam("PF00081").betweenAA(-45, -15)
```

This emits every methionine-started candidate in the inclusive coordinate
window without running TargetP or requiring DeepLoc. Each discovered candidate
passes the bare `Leader()` rule and continues through the other rules in the
rule set. If the Pfam anchor is missing or the window contains no methionine
start, the protein produces no candidate rows. The original input sequence is
not used as a fallback for a bare `Leader()` rule.

Supported TargetP predictions are:

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


## TF Motif Rules

`TFMotifs.has_within(distance, motif_a, motif_b, min_score_threshold=8)` scans
genomic locus sequence with `gimme scan` and evaluates motif hits in genomic
locus coordinates.

```python
from sieve.rules import TFMotifs, Rules

rule = Rules(
    TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP")
)
```

Motif names are prefix-matched, so `GM.5.0.Rel` matches a hit named
`GM.5.0.Rel.0001`. Hits below `min_score_threshold` are ignored. Strand is
recorded by Gimme but does not prevent two hits from pairing.

The rule searches for at least one qualifying hit for each motif. If both are
present and the nearest edges of any motif pair are within `distance`, the rule
returns `yes`. If both are present but no pair is close enough, the rule returns
`too_far`, which is treated as passing. Missing motifs return non-passing
`missing_*` values.

Possible TF motif outcomes include:

```text
yes
too_far
missing_GM.5.0.Rel
missing_GM.5.0.bZIP
missing_GM.5.0.Rel_and_GM.5.0.bZIP
false
error
```

`false` is used when the requested scope has no intervals to scan, such as
asking for introns in a single-exon gene.


## TF Motif Scopes

Without an explicit scope, TF motif rules consider all CDS/exon intervals and
all introns in the locus:

```python
TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP")
```

Use `.in_intron()` to search any intron, or `.in_intron(n)` to search a
specific 1-based intron:

```python
TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron()
TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_intron(2)
```

Use `.in_exon()` to search any CDS/exon interval, or `.in_exon(n)` to search a
specific 1-based CDS/exon:

```python
TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_exon()
TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").in_exon(2)
```

Use `.between(start, end)` to search a locus-relative genomic interval:

```python
TFMotifs.has_within(20, "GM.5.0.Rel", "GM.5.0.bZIP").between(81, 90)
```

For `.between(...)`, start and end may be supplied in either order. For intron
and exon scopes, numbering follows the gene's CDS order. On reverse-strand
genes, intron 1 and exon 1 are still counted from the gene direction rather than
from the leftmost genomic coordinate.

A TF motif rule can only be scoped once. For example,
`.in_intron(2).between(1, 10)` is invalid.


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

Require Pfam, KO, and a leader call:

```python
from sieve.rules import KO, Leader, Pfam, Rules

rule = Rules(
    Pfam.matches("PF00081")
    & KO.matches("K04564")
    & Leader().upstreamOfPfam("PF00081").betweenAA(-45, -15).is_mTP()
)
```

Require HMM coverage and a conserved motif at profile positions:

```python
from sieve.rules import HMMAlignment, Rules

profile = HMMAlignment("models/mnsod.hmm")

rule = Rules(
    profile.covers(1, 180)
    & profile.is_at("H", 27)
    & profile.matches_regex("D.E", 150)
)
```

Find proteins with a TF motif pair in the second intron:

```python
from sieve.rules import Rules, TFMotifs

rule = Rules(
    TFMotifs.has_within(
        20,
        "GM.5.0.Rel",
        "GM.5.0.bZIP",
        min_score_threshold=8,
    ).in_intron(2)
)
```

Combine leader rules with `|`:

```python
from sieve.rules import Leader, Rules

rule = Rules(
    Leader().is_mTP() | Leader().is_SP()
)
```
