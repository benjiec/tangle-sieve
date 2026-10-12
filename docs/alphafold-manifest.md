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
chain IDs. All models in a fold must agree. Only standard protein residues are
currently supported; ambiguous/modified residues and chains without protein
sequences cause errors rather than silently disappearing. Five models are typical,
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
- `sequence_id`: previous sequence-only fold hash: last eight lowercase hex
  characters of SHA-1 over compact JSON of the sorted sequences with copies,
  e.g. `["AC","AC","DE"]`. Both IDs are reference labels, not cache keys.
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
sequence-only ID preserves the previous hash definition exactly.

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
CIF contents, full-data JSON, and job requests are not read. Missing summary files,
missing fields, and null values print `NA`. No matches produce no output and a
successful exit status. Malformed JSON, invalid score values, duplicate model
numbers/summary members, and inaccessible folds are reported to stderr. Other
models/folds continue; any such error yields exit status 1. Invalid manifest
headers/rows or sequence IDs fail before any results are printed. The input
manifest is not modified.

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

Output remains the original five columns without a header; the combined score
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
