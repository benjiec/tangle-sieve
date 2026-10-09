import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("input_fasta")
parser.add_argument("output_a3m")
parser.add_argument("--left", default=None, type=int)
parser.add_argument("--right", default=None, type=int)
args = parser.parse_args()

records = list(SeqIO.parse(args.input_fasta, "fasta"))
query_seq = records[0].seq

# Find columns where the query sequence has non-gap characters
keep_indices = [i for i, char in enumerate(query_seq) if char != "-"]

with open(args.output_a3m, "w") as out:
    for rec in records:
        # Keep only the columns corresponding to the query backbone
        a3m_seq = "".join([rec.seq[i] for i in keep_indices])
        if args.right:
            a3m_seq = a3m_seq[:args.right]
        if args.left:
            a3m_seq = a3m_seq[args.left-1:]
        out.write(f">{rec.id}\n{a3m_seq}\n")
