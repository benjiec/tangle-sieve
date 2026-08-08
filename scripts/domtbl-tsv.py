import re
import argparse
ap = argparse.ArgumentParser()
ap.add_argument("domtblout")
args = ap.parse_args()

columns = (
  "target name",
  "accession",
  "tlen",
  "query name",
  "accession",
  "qlen",
  "seq e-value",
  "seq score",
  "seq bias",
  "dom #",
  "dom total",
  "dom c-Evalue",
  "dom i-Evalue",
  "dom score",
  "dom bias",
  "hmm from",
  "hmm to",
  "ali from",
  "ali to",
  "env from",
  "env to",
  "acc",
  "description"
)

print("\t".join(columns))

with open(args.domtblout) as f:
    for line in f:
        if line[0] == "#":
            line = line[1:]
        line = line.strip()
        if line.startswith("--"):
            continue
        if line.startswith("target"):
            continue
        if line.startswith("Program"):
            break
        line = re.sub(r"\s+", '\t', line)
        print(line)
