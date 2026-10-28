# AlphaFold sequence manifest

```sh
sieve-py ~/git/tangle/sieve/scripts/alphafold-manifest.py INPUT_DIRECTORY candidates.faa manifest.tsv
```

The input can also be a ZIP. The script recursively searches directories and
ZIPs (including ZIPs inside ZIPs), grouping `*_model_N.cif` files by directory.
Template CIFs and macOS metadata are ignored. Multiple locations remain separate,
even when their sequence composition is identical. Symlink children are skipped.
Nested ZIPs are spooled to temporary storage without extracting member paths.

The entire FASTA is validated before scanning or changing the manifest. Errors
include duplicate full FASTA headers (exact text after `>`), duplicate sequences, empty records/headers, sequence text before a header,
and nonalphabetic sequence characters. All detected errors are reported before
stopping. Sequence whitespace is removed and letters are uppercased; headers are
preserved in full after `>` (excluding the line ending), without accession parsing.
No substring matching, alignment, or terminal-stop removal is performed.

Sequences are read from CIF entity sequence tables and mapped to `_struct_asym`
chain IDs. All models in a fold must agree on base sequences and PTMs. Modified residues
are matched to the unmodified candidate FASTA sequence using the job request
and CIF metadata. Ambiguous mappings and chains without protein sequences
cause errors rather than silently disappearing. Five models are typical,
but a different nonzero number is accepted. Duplicate model numbers are rejected.
The manifest does not require confidence JSON or MSA files.

Each unmatched chain is reported on its own line, sequence first:
`ERROR: SEQUENCE not found in SOURCE_PATH::FOLD_PATH chain CHAIN_ID`. A fold
with any error contributes no rows; other folds continue. Exit status is 1 if any
errors occurred, otherwise 0. Errors and a final count go to stderr.

## TSV columns

One row per fold location, with 16 columns:

- `source_path`: absolute extracted fold directory, or outermost filesystem ZIP.
- `fold_path`: directory within that source, `.` at its root. Nested archive
  boundaries use `!/`, for example `inner.zip!/batch/fold`.
- `fold_locator`: value to pass directly as the positional input to
  `alphafold-pdockq2-patches.py` (quote it in the shell). Extracted folds use their
  absolute directory path; ZIP folds use `SOURCE.zip::FOLD_PATH`.
- `fold_id`: last eight lowercase hex characters of SHA-1 over UTF-8 compact JSON
  of `(seed, sorted_sequences_with_copies)`, e.g. `[42,["AC","AC","DE"]]`.
- `sequence_id`: for unmodified folds, the previous sequence-only fold hash: last eight lowercase hex
  characters of SHA-1 over compact JSON of the sorted sequences with copies,
  e.g. `["AC","AC","DE"]`. PTM-containing folds use the conditional payload
  below for both IDs. Both IDs are reference labels, not cache keys.
- `timestamp`: latest modification time among the fold's model CIF files.
- `input 1 description`, `input 1 copies`, through `input 5 description`,
  `input 5 copies`: all ten columns are always present. Each input represents a
  distinct sequence, with its full FASTA header and chain copy count in one model.
  Inputs are ordered by sequence lexically, matching the previous row ordering.
  Unused pairs are blank. More than five distinct sequences is a fold error;
  the fold is skipped without truncation. Individual chain IDs are no longer stored.

Standard TSV quoting preserves tabs and quotes within FASTA headers.

## Checkpoints and reruns

The two location fields identify a cache entry. After each successful fold, its
row is committed using a temporary file and atomic replacement.
Completed folds are skipped solely by the two paths. Changes to FASTA sequences,
headers, or model files do not refresh cached rows. Delete the row for a fold's
two paths to force reprocessing, or use a fresh TSV. FASTA validation still runs
first on every invocation; invalid FASTA leaves the TSV untouched.

Existing entries outside the current input, or no longer discoverable, are retained;
the script does not prune historical locations. Do not run multiple writers against
one manifest simultaneously. Discovery still lists directories and opens archives,
including nested archives, but cached model contents are not read or hashed.
Output-write failures stop the run.

Manifests from before `sequence_id` was added must be removed and regenerated;
there is no ID migration. Existing files with an old schema are rejected unchanged.

## Seed validation

Each fold needs its matching `*_job_request.json`, containing exactly one job and
exactly one `modelSeeds` value. Nonnegative integer seeds and digit strings are
accepted and normalized to integers, so `"42"` and `42` produce the same fold ID.
Missing, invalid, or multiple seeds cause a fold error; processing continues with
other folds. Duplicate job-request archive members are rejected. The manifest's
sequence ID preserves the previous hash definition exactly for unmodified folds.

## PTM-aware identity

When any chain has a PTM, the sequence hash input is compact JSON of
`["ptm-v1", sorted_chain_records]`. Each record contains the unmodified sequence
and its sorted list of `[position, CCD_code]` pairs. Positions are one-based.
Every chain copy contributes a record, including unmodified copies. Chain names
and the order of chains or modification entries do not affect the IDs. AF Server
codes such as `CCD_SEP` are normalized to `SEP`. The fold hash uses
`[seed, sequence_hash_input]`. Folds without PTMs retain exactly their previous
sequence and fold IDs, including when the modification list is empty.

AF Server `proteinChain` requests and AF3 `protein` requests are checked against
CIF monomers, chain copies, and available atom residue metadata. Without request
sequence metadata, a standard parent must be identifiable from CIF parent
metadata or Biopython's residue mapping. Unknown or ambiguous parents, invalid
positions, duplicate modification positions, and inconsistent models are errors.

FASTA matching and input description/copy columns still use the base sequence;
modified and unmodified copies of that sequence share the same description.
Their different modification states nevertheless contribute separately to IDs.
Cached rows are not refreshed: delete any previously cached PTM row to reprocess
it. No manifest columns are added.

Identity support does not establish support for PTM chemistry. Rosetta
preparation and overlap assembly currently reject PTM-containing structures
explicitly rather than discard modified residues. Patch analysis supports PTM
confidence and geometry as described below.

### PTM tokens in patch analysis

A modified residue can have several confidence tokens. The patch parser retains
all of them, accepting multiple tokens only for nonstandard residues whose token
count matches their observed heavy-atom count. Duplicate tokens for standard
residues and inconsistent mappings remain errors.

Contacts and pLDDT retain the existing Cα-based definitions. For each contacting
residue pair, the script averages `1 / (1 + (PAE / 10)^2)` over all token pairs,
separately in each direction. It then averages those values across residue pairs,
so an atom-tokenized residue does not receive extra weight. This is a documented
extension of the scoring implementation, not a PTM-specific pDockQ2 calibration.
Single-token residues retain the previous calculation.

For PTM-containing models, the summary retains contact geometry, Cα pLDDT,
confidence, and SASA burial with all modified atoms present. Modified residues
appear as bracketed CCD codes, such as `[SEP]`. DSSP runs on the original CIF;
`--structure` reports how many modified residues received assignments. Unassigned
residues remain `?`, and DSSP warnings remain visible.

### Phosphorylation descriptors

SEP, TPO, and PTR are recognized as phosphorylated residues. For each patch and
sequence/spatial shell:

- Ordinary-residue chemistry remains available. Hydropathy and hydrophobic
  fraction are `n/a` only for sets containing modified residues; modified residues
  are never assigned the parent amino acid's hydropathy.
- Nominal charge is reported under two explicit scenarios: every phosphate is
  −1, or every phosphate is −2, added to the existing ordinary-residue charges.
  These are sensitivity scenarios, not inferred protonation states or a pH model.
  Unknown PTMs make charge unavailable for the affected set. K/R, D/E, and H
  counts refer to ordinary residues; phosphate charges are additional.
- `Phosphate/basic contacts` count unique residue pairs whose observed terminal
  phosphate O1P/O2P/O3P atoms approach Lys NZ or Arg NE/NH1/NH2 within 4 Å.
  The nearest atom pair is printed when either endpoint belongs to the patch,
  with `within patch` or `one endpoint in patch` scope. These are distance-only
  candidates, not assigned hydrogen bonds or interaction energies. Ordinary
  Asp/Glu salt-bridge candidates remain separately reported.
- Steric-overlap screening includes phosphorus with a 1.80 Å radius. Unsupported
  elements cause an explicit error instead of receiving a zero radius. Missing
  phosphate atoms produce a warning; contacts use observed atoms only.
- Nonpolar-contact classification still covers ordinary residues only. The
  carbon/sulfur burial calculation retains all PTM atoms.

Chemical names and atom definitions follow the wwPDB CCD entries for
[SEP](https://files.rcsb.org/ligands/view/SEP.cif),
[TPO](https://www.rcsb.org/ligand/TPO), and
[PTR](https://files.rcsb.org/ligands/view/PTR.cif).
The −1/−2 alternatives reflect the monoanionic/dianionic phosphate states studied
in [phosphorylated model peptides](https://pubmed.ncbi.nlm.nih.gov/7529751/);
this script does not estimate their populations in a protein.
Secondary structure uses [DSSP](https://github.com/PDB-REDO/dssp), without residue
renaming or coordinate changes. Coverage confirms assignments were returned,
not that every assignment is experimentally correct.

Validation on the supplied SEP archive with mkdssp 4.6.1 returned assignments
for all six SEP residues in each of five models. That executable also emitted
`Unknown compound: SEP`; the script preserves the warning. TPO/PTR chemistry
rules have synthetic tests; their DSSP behavior still needs representative
real structures for validation.

## Timestamp meaning

Filesystem models use modification time in UTC. ZIP models prefer the Unix mtime
in the extended timestamp field when available for every model, also in UTC.
UTC values use ISO 8601 with `Z`, for example `2026-10-07T17:47:34Z`. Otherwise,
the latest ZIP date/time is used without a timezone suffix because the ZIP DOS
field does not specify a timezone. Nested ZIPs use their model members' timestamps,
not the containing archives' timestamps. Templates and macOS sidecars are ignored.

These are file modification times, not verified prediction-completion times.
Copying, repacking, or extracting files can affect them. Model contents need not
be read to obtain timestamps.

## Load one fold into the patch script

Copy a `fold_locator` value from the TSV and pass it as one quoted argument:

```sh
sieve-py ~/git/tangle/sieve/scripts/alphafold-pdockq2-patches.py \
  '/path/batch.zip::folder/fold'
```

Nested example: `/path/batch.zip::inner.zip!/folder/fold`; a fold at the ZIP root
uses `/path/batch.zip::.`. Only the selected fold is loaded. Each model requires
its matching `*_full_data_N.json`; the manifest itself does not verify those JSON
files. Cache identity remains the original two paths.

The separators `::` and `!/` are reserved for locator syntax; source/member names
containing these delimiters are not supported by this representation.

## List model confidence results

```sh
sieve-py ~/git/tangle/sieve/scripts/alphafold-list-results.py manifest.tsv b911edee
```

Selects rows with the exact `sequence_id` and prints headerless TSV to stdout:
`fold_locator`, `fold_id`, `model_id`, `ipTM`, `pTM`. Folds retain manifest order;
models are sorted numerically. IDs and locators are copied from the manifest,
without recomputing hashes. Distinct fold locations remain separate.

Only model filenames and matching `*_summary_confidences_N.json` files are needed;
CIF contents and full-data JSON are not read; job requests are read only with `--seed`. Missing summary files,
missing fields, and null values print `NA`. No matches produce no output and a
successful exit status. Malformed JSON, invalid score values, duplicate model
numbers/summary members, and inaccessible folds are reported to stderr. Other
models/folds continue; any such error yields exit status 1. Invalid manifest
headers/rows or sequence IDs fail before any results are printed. The input
manifest is not modified.

### Include the seed

Add `--seed` to append a sixth, headerless `seed` column. The seed is read once
per fold from `modelSeeds` in its matching `*_job_request.json`. Missing job
requests print `NA`; malformed, duplicate, or multiple-seed requests report an
error, skip that fold, and yield exit status 1 while other folds continue.
`--highest` preserves the selected row's seed. Without `--seed`, output remains
five columns. Shared pipeline readers accept either format.

```sh
sieve-py ~/git/tangle/sieve/scripts/alphafold-list-results.py manifest.tsv b911edee --highest --seed
```

### Select the highest-scoring model

Without `--highest`, all matching model rows are printed as above. With it, exactly
one eligible row is selected across all matching fold locations and models:

| Option | Selection score |
| --- | --- |
| `--highest` | `0.8*ipTM + 0.2*pTM` |
| `--highest ipTM` | ipTM |
| `--highest pTM` | pTM |
| `--highest ipTM8_pTM2` | `0.8*ipTM + 0.2*pTM` |

```sh
# All matching models
sieve-py ~/git/tangle/sieve/scripts/alphafold-list-results.py manifest.tsv b911edee

# Highest combined score (bare flag defaults to ipTM8_pTM2)
sieve-py ~/git/tangle/sieve/scripts/alphafold-list-results.py manifest.tsv b911edee --highest

# Highest ipTM
sieve-py ~/git/tangle/sieve/scripts/alphafold-list-results.py manifest.tsv b911edee --highest ipTM

# Highest pTM
sieve-py ~/git/tangle/sieve/scripts/alphafold-list-results.py manifest.tsv b911edee --highest pTM

# Explicit combined score
sieve-py ~/git/tangle/sieve/scripts/alphafold-list-results.py manifest.tsv b911edee --highest ipTM8_pTM2
```

Place a bare `--highest` after both positional arguments so the next positional
argument is not interpreted as its optional value. Metric names are case-sensitive.

A row missing a required score is ineligible: ipTM selection requires only ipTM,
pTM selection requires only pTM, and combined selection requires both. Zero is a
valid score. If no rows qualify, stdout is empty and the exit status is zero unless
an input error occurred. Ties retain the first row in manifest order, then numeric
model order. Combined scores use decimal arithmetic on the loaded score values.

Output remains headerless, with five columns (six with `--seed`); the combined score
is not added. Errors still go to stderr and yield exit status 1. If some inputs
fail, the best row among successfully read eligible models is still printed.

### Pipe selected results into patch analysis

```sh
sieve-py ~/git/tangle/sieve/scripts/alphafold-list-results.py manifest.tsv b911edee --highest | \
  sieve-py ~/git/tangle/sieve/scripts/alphafold-pdockq2-patches.py -
```

A result row selects only its `model_id` at its `fold_locator`. Other models'
CIFs and confidence files are not loaded. Omitting `--highest` sends all result
rows, each analyzed separately in input order. The five-column headerless TSV
schema, serialization, and parsing are shared by both scripts.

The patch script also accepts one-field locator records on stdin; these analyze
all models exactly as a direct locator argument would. Blank records are ignored.
Empty stdin produces no output and leaves any existing `--output` file untouched.
TSV quoting supports tabs, quotes, and newlines within paths.

Piped IDs and scores are validated as fields but do not override the fold's
metadata: displayed IDs and confidence scores are read/computed from fold files.
Malformed records and missing requested models cause errors, with stdin record
context. The patch script remains fail-fast and buffers output until all requested
analyses succeed. An optional `--output` TSV contains results from all records.

For shell pipelines where an upstream failure must affect the pipeline exit
status, enable `set -o pipefail`; the consumer cannot infer upstream errors from
an empty or partial valid stream.

## Assembled complex directories

Directories written by `alphafold-overlap.py` can be registered through the same
manifest command. Their `source_type` is `assembled`; AF rows are `alphafold`.
Current-format manifests without this column are upgraded while preserving cached
rows and IDs. Assembly IDs come from validated provenance, not a synthetic seed.
The sequence-ID algorithm and five-distinct-sequence limit are unchanged.

Result listing supports these directories and ZIP locators, printing `NA` for
ipTM/pTM/seed. `--highest` excludes these unscored rows. See
[overlap and derived complexes](alphafold-overlap.md) for creation, recursive use,
provenance, and artifact verification.

## Find folds by candidate names

```sh
sieve-py ~/git/tangle/sieve/scripts/alphafold-find.py manifest.tsv sequence_names.txt
sieve-py ~/git/tangle/sieve/scripts/alphafold-find.py manifest.tsv sequence_names.txt \
  --cross --max-chains 6 --max-unique-inputs 2 --highest
```

The names file contains one full FASTA header per line, without the leading `>`.
Matching is exact and case-sensitive; spaces are preserved, and only line endings
are removed. Empty lines and repeated names are ignored. An empty names file or
no matches produces no output. Names absent from the manifest do not cause errors.

By default, a fold qualifies if any input description is listed. `--cross`
requires every input description to be listed; it does not require every listed
name to occur in each fold. Self-interactions qualify. Optional positive integer
limits are inclusive and can be combined with each other and `--cross`:

- `--max-chains N`: sum of all input copy counts must be at most N.
- `--max-unique-inputs N`: number of distinct full descriptions must be at most N.

Neither limit has a default restriction. A×3 + B×3 has six chains and two unique
inputs. A×2 has two chains and one unique input, whether A appears in one column
with two copies or two columns with one copy each.

Output is the shared headerless five-column TSV: `fold_locator`, `fold_id`,
`model_id`, `ipTM`, `pTM`, suitable for piping into the patch and Rosetta scripts.
Each locator is processed once even if multiple names or repeated manifest rows
match. Different locators remain separate, including those with the same fold ID.
Contradictory IDs or aggregated input compositions for one locator are errors.

Without `--highest`, every model in each qualifying fold is output, in manifest
order then numeric model order; missing confidence values are `NA`.
`--highest [ipTM|pTM|ipTM8_pTM2]` outputs one best model per sequence ID across
qualifying folds. Omitting the value selects `0.8*ipTM + 0.2*pTM`. Models missing
values required by the selected metric are ineligible; a sequence ID with no
eligible model produces no row. Ties keep the first manifest fold, then the lowest
model number. Sequence-ID groups are output in their first qualifying manifest
occurrence order.

The entire manifest is validated before folds are opened or results written.
All five description/copies column pairs are required. Unused pairs must be blank;
used pairs must contain a description and a positive integer copy count. Errors
reading individual folds/models go to stderr; other results continue and the
exit status is nonzero. Normal output has no headers, summaries, or titles.
