import re
import argparse
ap = argparse.ArgumentParser()
ap.add_argument("domtblout")
args = ap.parse_args()

columns = (
  "target name",
  "target accession",
  "query name",
  "query accession",
  "seq e-value",
  "seq score",
  "seq bias",
  "best dom e-value",
  "best dom score",
  "best dom bias",
  "doms exp",
  "doms reg",
  "doms clu",
  "doms ov",
  "doms env",
  "doms dom",
  "doms rep",
  "doms inc",
  "description of target"
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
