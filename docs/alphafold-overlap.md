# AlphaFold overlap and derived complexes

```sh
sieve-py ~/git/tangle/sieve/scripts/alphafold-overlap.py locatorA locatorB --output complex
sieve-py ~/git/tangle/sieve/scripts/alphafold-overlap.py complex locatorC --output next_complex
```

Inputs are AF ZIP locators (including nested ZIPs), extracted AF fold directories,
assembled directories, or standalone CIFs with full polymer sequence metadata.
Exactly one sequence must be shared, occurring exactly once in each input.
Both inputs must have a nonshared partner. All chains are preserved except the
second copy of the shared chain. Input A's shared coordinates are retained;
all input B coordinates undergo the same proper rigid transformation.
No repacking, relaxation, sequence mutation, or fabricated confidence is performed.

Default model is 0. `--model-a N` and `--model-b N` choose alternatives.
`--all-models` adds every combination to comparisons.tsv, but saves only the
selected pair as output model 0. Other combination errors are reported in the
table and stdout; an error in the primary pair prevents output. A successful
run, including a run with structural warnings, exits 0; input/primary alignment
errors exit 1. Output must be new or empty; publication is atomic.

## Quick assessment and thresholds

The first stdout line after identifiers and the beginning of the generated README
show a specific warning assessment or `No detected geometric conflict`.
When shared-chain replacement conflicts are the only warning, the assessment is
`WARNING — Partners do not clash; merged shared-chain interface requires refinement`.
Alignment problems, cross-partner clashes, alternative-alignment clashes, and
inherited input warnings receive separate messages; multiple issues are listed
together. Stdout identifies the selected models. With `--all-models`, the final
count means **combinations with warnings**, not individual warning messages, and
points to `comparisons.tsv` for every combination's details. Errors print `ERROR — cannot assemble`. Warnings do not suppress output:
the coordinates remain useful for inspection. Neither label establishes binding,
affinity, thermodynamic stability, or complete stereochemical validity.

All thresholds below are configurable screening choices, **not validated binary
acceptance criteria for protein binding**:

| Option | Default | Meaning |
| --- | --- | --- |
| `--rmsd-warning` | 2 Å | Warn if whole shared-chain Cα RMSD exceeds this |
| `--site-rmsd-warning` | 1.5 Å | Warn on a site's Cα RMSD under the whole-chain fit |
| `--coverage-warning` | 0.9 | Warn if shared or observed binding-site Cα coverage falls below this fraction |
| `--clash-overlap` | 0.4 Å | Report atom pairs with radius sum minus distance above this |
| `--site-distance` | 8 Å | A shared Cα belongs to a site when near a partner heavy atom |
| `--contact-distance` | 4.5 Å | Heavy-atom distance defining a residue-pair contact for before/after comparisons |

Also warn when >10% of matched Cα deviations exceed 3 Å, when replacing the shared
chain introduces/worsens overlaps, when alternative fits reveal clashes, or when
an input assembly already carries warnings. The 3 Å/10% rule is an additional
heuristic. High terminal deviations are shown even when they do not trigger it.
Fewer than three matched Cα atoms or collinear/degenerate coordinates are errors.
No unobserved coordinates are generated; sequence coverage uses full declared sequences.

**Sources and limitations.** MolProbity uses approximately 0.4 Å overlap to identify
serious nonbonded conflicts, with explicit hydrogen and interaction-aware treatment
([Chen et al., 2010](https://journals.iucr.org/d/issues/2010/01/00/dz5180/dz5180.pdf)).
This tool uses only elemental heavy-atom radii (C 1.70, N 1.55, O 1.52, S 1.80 Å),
without hydrogen-bond or covalent/disulfide exceptions. Its output is geometric
**overlap candidates, not MolProbity clashscore**; favorable donor–acceptor contacts
can be flagged and must be inspected. Atom-pair counts are not independent clashes
or energetic penalties. Backbone counts are reported separately. Nonprotein
entities, unsupported elements, alternate conformations, and ambiguous atom/sequence
mapping produce errors rather than silent exclusion.

CombFold provides precedent for assembling predicted subcomplexes and screening
clashes, but uses different confidence-aware filters and alignment policies
([Shor and Schneidman-Duhovny, 2024](https://www.nature.com/articles/s41592-024-02174-0)).
Our numerical RMSD and coverage cutoffs are heuristics, not thresholds validated by
that paper. No single RMSD threshold guarantees compatible local interfaces.
Compare the local deviations and sensitivity fits, then test whether restrained
joint refinement preserves both interfaces. Experimental simultaneous binding is
still needed to validate biological occurrence.

## File reference

Every TSV has a header, including empty clash tables. `NA` means unavailable.
Chain IDs in analytical TSVs refer to the original inputs unless stated otherwise;
residue numbers are **declared sequence positions**, not original author numbering.

| File | Contents and interpretation |
| --- | --- |
| `model0_complex.cif` | Coordinates of the selected assembly, complete entity/sequence/chain tables; safe for recursive loading |
| `provenance.json` | Versioned format, IDs, CIF checksum, inputs and embedded ancestry, chain/residue mappings, known confidence origins, rigid transform, settings, assessment |
| `comparisons.tsv` | Global shared-chain Cα RMSD, coverage, largest deviation, cross-partner overlap totals and status for each combination; selected combination is also summarized in stdout |
| `shared_chain_deviations.tsv` | Selected pair: deviation at each matched shared sequence position and known source pLDDT; locate local differences and terminal outliers |
| `shared_chain_sites.tsv` | Per-partner shared-site RMSD **under the whole-chain fit**, matched/observed site counts and coverage |
| `alignment_sensitivity.tsv` | Independent fits on whole chain, terminal-trimmed chain, each site, and source pLDDT≥80 subset when known; fitted RMSD and resulting partner clashes. These fits do not alter the saved assembly |
| `cross_partner_metrics.tsv` | One row per nonshared input-A/input-B chain pair: closest distance, greatest signed overlap, contact count and overlap counts |
| `cross_partner_clashes.tsv` | Detailed new cross-input atom overlaps; negative max overlap in summary means separation, not a clash |
| `shared_interface_changes.tsv` | Each input-B partner against its original shared chain versus the retained input-A chain: overlaps introduced/worsened/resolved, contacts lost/gained |
| `shared_interface_clashes.tsv` | Atoms in assembled shared-chain interfaces; distinguishes introduced, worsened and preexisting overlap candidates |
| `shared_interface_contact_changes.tsv` | Exact residue pairs lost/gained at the heavy-atom contact cutoff |
| `README.md` | Run assessment, output-chain mapping, and this reference |

**Where to assess fold similarity:** start with shared-chain RMSD and coverage in
stdout/comparisons.tsv, then shared_chain_sites.tsv and shared_chain_deviations.tsv.
An independently fitted site RMSD in alignment_sensitivity.tsv can be low even
when that site has moved relative to the rest of the chain; it must not replace
the site's deviation under the global fit. Missing-coordinate regions cannot be
assessed. Structural similarity is not evidence that either prediction is correct.

## Identity, provenance, and recursive use

The directory path is the locator. The sequence ID uses the unchanged sorted,
copy-counted sequence hash. The fold ID is the last eight SHA1 hex characters of
canonical JSON containing the `sieve-assembly-v1` domain marker and assembly identity.
Identity includes ordered input CIF SHA256 hashes, selected model numbers, alignment
policy, shared-chain choices and output CIF hash. Paths, timestamps, and diagnostic
thresholds do not change identity. Input order can change coordinates and identity.
Full hashes are preserved because eight-character IDs are labels, not collision-proof
content verification. Edited CIFs fail checksum validation and must be regenerated.

Output chain IDs are unique single ASCII alphanumeric characters for Rosetta
compatibility (maximum 62). Sequence positions become consistent label/auth residue
numbers; original author numbering is retained in provenance. Full sequences include
residues without coordinates. Standalone CIFs need `_entity`, `_entity_poly_seq`,
`_struct_asym` and unambiguous atom-to-sequence mapping. Ordinary B-factors are never
assumed to be pLDDT; known AF source values are inherited and labeled accordingly.

Recursive loading needs only the output directory; input paths are historical
references. Standalone exported CIFs retain their full sequence metadata, but the
reserved model0_complex.cif requires its provenance to verify derived identity.
Assembly-order sensitivity and input warnings remain relevant at every step.

## Existing scripts

```sh
sieve-py ~/git/tangle/sieve/scripts/alphafold-manifest.py complex candidates.faa manifest.tsv
sieve-py ~/git/tangle/sieve/scripts/alphafold-list-results.py manifest.tsv SEQUENCE_ID
sieve-py ~/git/tangle/sieve/scripts/alphafold-pdockq2-patches.py complex --chains A C
sieve-py ~/git/tangle/sieve/scripts/alphafold-rosetta-integrity.py complex --model 0
```

Manifest `source_type` distinguishes `alphafold` and `assembled`. Older current-format
manifests gain that column without losing rows or IDs. The existing maximum of five
distinct sequences remains; larger complexes work directly but cannot be registered.
Use a new immutable output directory for each assembly; manifest cache identity
remains location-based. Directory and ZIP locators both support derived complexes.

Result listing emits `NA` for ipTM/pTM/seed. `--highest` excludes rows without its
required confidence values; it will not select an assembled structure. Five- and
six-column pipes remain supported.

Patch analysis on assemblies is geometric: PAE/pDockQ2 are unavailable (`NA`), and
patches are ranked by contact count with deterministic ties. Known source pLDDT is
explicitly inherited, never joint confidence. Pairwise burial still has the existing
pairwise meaning; this change does not implement joint multichain patch scoring.
Rosetta integrity/interface use the derived identity but otherwise keep their existing
coordinate-analysis settings and limitations.
