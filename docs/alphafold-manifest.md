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
include duplicate sequences, empty records/headers, sequence text before a header,
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

One row per distinct sequence within a fold location:

- `source_path`: absolute extracted fold directory, or outermost filesystem ZIP.
- `fold_path`: directory within that source, `.` at its root. Nested archive
  boundaries use `!/`, for example `inner.zip!/batch/fold`.
- `fold_locator`: value to pass directly as the positional input to
  `alphafold-pdockq2-patches.py` (quote it in the shell). Extracted folds use their
  absolute directory path; ZIP folds use `SOURCE.zip::FOLD_PATH`.
- `fold_id`: last eight lowercase hex characters of SHA-1 over UTF-8 compact JSON
  of the lexically sorted sequence list, including one entry per chain copy.
  Example hash input: `["AC","AC","DE"]`. This is a reference label, not a
  unique identifier or cache key.
- `description`: full FASTA header.
- `chain_ids`: comma-separated sorted CIF chain IDs sharing the sequence.
- `copies`: number of those chains in one model (not summed across models).
Standard TSV quoting preserves tabs and quotes within FASTA headers.

## Checkpoints and reruns

The two location fields identify a cache entry. After each successful fold, all
its rows are committed together using a temporary file and atomic replacement.
Completed folds are skipped solely by the two paths. Changes to FASTA sequences,
headers, or model files do not refresh cached rows. Delete all rows for a fold's
two paths to force reprocessing, or use a fresh TSV. FASTA validation still runs
first on every invocation; invalid FASTA leaves the TSV untouched.

Existing entries outside the current input, or no longer discoverable, are retained;
the script does not prune historical locations. Do not run multiple writers against
one manifest simultaneously. Discovery still lists directories and opens archives,
including nested archives, but cached model contents are not read or hashed.
Output-write failures stop the run.

Older manifests containing `sequence`, `sequence_sha256`, or `source_fingerprint`
are automatically rewritten to the seven columns above after FASTA and cache
validation, preserving cached locations.

## Load one fold into the patch script

Copy a `fold_locator` value from the TSV and pass it as one quoted argument:

```sh
sieve-py ~/git/tangle/sieve/scripts/alphafold-pdockq2-patches.py \
  '/path/batch.zip::folder/fold'
```

Nested example: `/path/batch.zip::inner.zip!/folder/fold`; a fold at the ZIP root
uses `/path/batch.zip::.`. Only the selected fold is loaded. Each model requires
its matching `*_full_data_N.json`; the manifest itself does not verify those JSON
files. Existing manifests gain the locator column on rerun without reprocessing
cached folds. Cache identity remains the original two paths.

The separators `::` and `!/` are reserved for locator syntax; source/member names
containing these delimiters are not supported by this representation.
