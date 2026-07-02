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

## Calling External Tools

You can call these separately. E.g. targetP

```
docker run --rm --platform linux/amd64 -v .:/data local-targetp:2.0 -fasta query.faa -org non-pl -format short -stdout
```

or gimmemotifs

```
gimme scan locus.fna -b -c 0.85
```
