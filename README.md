# Sieve

Untangle protein sequences curated using the Tangle schema.


## Setup

Requires Python 3.13 (3.14 won't work because pyarrow issues).

Create a virtual env and run the following to install dependencies

```
python -m pip install --only-binary=:all: pyarrow==20.0.0
pip install -e .
```

You will need the following two external tools

  * TargetP: e.g. you can use a docker image, by default it is assumed to be "local-targetp:2.0"
  * Gimme: installed via pip

Aliase `sieve-py` to the python executable in the virtualenv. E.g.

```
alias sieve-py='venv-sieve/bin/python3'
```


## Workflows

Define rules for proteins using abstractions in `sieve/rules.py`.

In a Tangle setup (see `tangle/README.md` for environment variables and
directory structure), use the following to select HMM-detected and NCBI
curated proteins matching a specific KO (K04564 is an example), then build an
artifact directory containing the original and leader-candidate sequences:

```
sieve-py sieve/scripts/ko-find-matches.py K04564 --taxon cnidaria | \
  PYTHONPATH=coral sieve-py sieve/scripts/build-protein-artifacts.py \
    -r coral.rules.mt_MnSOD_cnidaria.rule \
    --artifacts-dir res-mt_MnSOD_cnidaria
```

`build-protein-artifacts.py` reads tab-separated protein and genome accessions
from standard input, so it can be composed with `ko-find-matches.py` or
`pfam-find-matches.py`. It uses curated Pfam detections when a leader rule is
anchored to a Pfam domain.

An artifact directory contains:

* `input.faa`: original protein sequences.
* `sequences.faa`: discovered candidate sequences.
* `sequences.tsv`: candidate-to-original associations and start metadata.
* `genomic_loci.tsv`: genomic locus metadata for curated inputs.

`sequences.tsv` also records an inclusive `end aa 1b` when candidate sequences
are C-terminally bounded.

Builders are incremental. Re-adding identical sequences and metadata is
idempotent, and newly discovered candidates are added. An accession that has a
different sequence is rejected. If a curated protein accession occurs in more
than one genome with the same sequence, the first genome is retained and later
genomes are ignored; their candidates and locus rows are not added.

If using DeepLoc to classify leaders, then submit the FASTA in
`res-mt_MnSOD_cnidaria/sequences.faa` to DeepLoc and download results to
`deeploc.csv`.

`ko-find-matches.py --taxon` filters KO matches by an exact,
case-insensitive taxonomy value at any supported rank before candidate
discovery and rule evaluation. Genomes without a matching taxonomy row are
excluded.

Use `--match-starts-before POSITION` to require the KO hit to start at or
before a 1-based protein position (`query_start <= POSITION`). Use
`--match-ends-before POSITION` to require the hit to end at or before a
position (`query_end <= POSITION`). When both options are supplied, the same
hit must satisfy both constraints. A protein with multiple hits is included if
any one hit satisfies all supplied constraints. For example:

```
sieve-py sieve/scripts/ko-find-matches.py K04564 \
  --taxon cnidaria \
  --match-starts-before 50 \
  --match-ends-before 200
```

Pass `-o OUTPUT.faa` to write the full matched protein sequences as FASTA
instead of printing the protein and genome accessions as tab-separated rows.

Use `--max-evalue-rank VALUE` to retain only detection rows whose
`custom_metric_name` is `evalue-rank` and whose numeric `custom_metric_value`
is at most `VALUE`. The default is `1`, so rows without an `evalue-rank` metric
are excluded unless a different maximum is specified. Position constraints are
applied to these same rows. For example, this requires one hit to have e-value
rank at most 3 and start at or before protein position 50:

```
sieve-py sieve/scripts/ko-find-matches.py K04564 \
  --max-evalue-rank 3 \
  --match-starts-before 50
```

Evaluate the completed artifact directory with fresh Pfam, KO, and HMM
alignment results:

```
PYTHONPATH=coral sieve-py sieve/scripts/check-artifacts-by-rules.py \
  -r coral.rules.mt_MnSOD.rule \
  --artifacts-dir res-mt_MnSOD_cnidaria \
  --pfam-hmm assets/pfam.hmm \
  --ko-hmm assets/k04564.hmm \
  --ko-thresholds assets/k04564_threshold.tsv \
  --deeploc-csv deeploc.csv
```

Rule evaluation never discovers leaders. It uses only the candidates recorded
by the build step. Pfam searches use HMMER gathering thresholds (`--cut_ga`).
KO searches disable gathering thresholds and require the model-specific
threshold file supplied with `--ko-thresholds`.

Dump sequences (original protein sequence plus discovered leaders, if leader
discovery was performed) that passed all the rules to a FASTA file

```
sieve-py sieve/scripts/rule-results-to-fasta.py \
  --artifacts-dir res-mt_MnSOD_cnidaria --output mt_MnSOD_cnidaria.faa --taxon cnidaria
```

The --taxon option can filter results by taxonomy (this can be the literal
value - not just a word - matching any of domain, phylum, family, class, order,
genus)

Construct MSA and HMM profiles

```
docker run --platform linux/amd64 \
  --rm -v .:/app pegi3s/muscle \
  -in /app/mt_MnSOD_cnidaria.faa -out /app/mt_MnSOD_cnidaria.msa

hmmbuild mt_MnSOD_cnidaria.hmm mt_MnSOD_cnidaria.msa
```

Use the following to search every candidate in an artifact directory and
compute a threshold for a new HMM profile. Candidate rows whose `pass all`
value is `true` are labeled positive; all other rows are labeled negative.
Threshold discovery deliberately searches without `--cut_ga`.
The output appends the HMMER `--domtblout` used for scoring as comment-prefixed
lines under `# hmmsearch domtblout`, so individual domain scores remain
available without making the threshold TSV invalid.

```
PYTHONPATH=coral sieve-py sieve/scripts/hmmsearch-threshold.py \
  --hmm mt_MnSOD_cnidaria.hmm \
  --artifacts-dir mt_MnSOD_cnidaria \
  --output mt_MnSOD_cnidaria_thresholds.tsv
```

Bound FASTA sequences between one N-terminal and one C-terminal profile-HMM
domain hit. Both searches use `--cut_ga`, then `hmmalign` maps the requested HMM
positions to exact sequence coordinates. A sequence is omitted with a report on
stderr unless it has exactly one hit from each model and both hits span their
requested model positions. Output accessions record the inclusive source
coordinates as `X_bounded_A-B`:

```
PYTHONPATH=coral sieve-py sieve/scripts/bound-fasta-by-hmms.py \
  --n-hmm n-terminal.hmm \
  --c-hmm c-terminal.hmm \
  --fasta proteins.faa > bounded.faa
```

Use `--n-position 2` or `--c-position 197` to require alternative HMM
alignment endpoints. Use `--acc-desc DESC` to append `|DESC` to every output
accession; for example, `X_bounded_A-B|DESC`. Descriptions cannot contain
whitespace.

Build the same artifact format from a protein FASTA. Supply `--pfam-hmm` when
Pfam detections are needed for anchored leader discovery:

```
PYTHONPATH=../coral sieve-py ../sieve/scripts/build-fasta-artifacts.py \
  --fasta transcript_proteins.faa \
  -r coral.rules.mt_MnSOD_cnidaria.rule_fasta \
  --artifacts-dir tmp \
  --pfam-hmm assets/pfam.hmm
```

If the construction rule contains `KO.matches("K04564", bound_cterm=True)`,
also supply `--ko-hmm`. The builder searches without `--cut_ga` or a threshold
file and truncates candidates at the greatest matching inclusive `query_end`:

```
PYTHONPATH=../coral sieve-py ../sieve/scripts/build-fasta-artifacts.py \
  --fasta transcript_proteins.faa \
  -r coral.rules.mt_MnSOD_cnidaria.rule_fasta \
  --artifacts-dir tmp \
  --ko-hmm assets/k04564.hmm
```

The FASTA and curated builders can add records to the same artifact directory.
Candidate accessions must remain globally unique; repeated protein accessions
follow the first-genome policy described above.

To use DeepLoc, submit `tmp/sequences.faa` after all build steps and download
the results to a file such as `deeploc.csv`. Then evaluate the artifacts:

```
PYTHONPATH=../coral sieve-py ../sieve/scripts/check-artifacts-by-rules.py \
  -r coral.rules.mt_MnSOD_cnidaria.rule_fasta \
  --artifacts-dir tmp \
  --pfam-hmm assets/pfam.hmm \
  --ko-hmm assets/k04564.hmm \
  --ko-thresholds assets/k04564_threshold.tsv \
  --deeploc-csv deeploc.csv
```

The Pfam and KO HMM files are also registered for `HMMAlignment(...)` rules
that refer to their basenames. Additional alignment profiles can be supplied
through repeatable, non-recursive `--hmm-dir` options.


## Calling External Tools

Some of the above scripts call the following tools.

Use targetP for leader classification,

```
docker run --rm --platform linux/amd64 -v .:/data local-targetp:2.0 -fasta query.faa -org non-pl -format short -stdout
```

Use gimmemotifs to detect TF motifs

```
gimme scan locus.fna -b -c 0.85
```
