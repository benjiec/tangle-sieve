import os


RULE_RESULTS_TSV = "rule-results.tsv"
SEQUENCES_FASTA = "sequences.faa"


def rule_results_tsv(artifacts_dir):
    return os.path.join(artifacts_dir, RULE_RESULTS_TSV)


def sequences_fasta(artifacts_dir):
    return os.path.join(artifacts_dir, SEQUENCES_FASTA)
