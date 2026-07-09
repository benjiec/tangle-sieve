from tangle import open_file_to_write


def write_rule_rows(rows, output_tsv, rules):
    rules.write_rows(output_tsv, rows)


def write_rule_fasta(rows, fasta_output, candidates_for_row):
    with open_file_to_write(fasta_output, "wt") as f:
        for row in rows:
            for candidate in candidates_for_row(row):
                f.write(f">{candidate.accession}\n{candidate.sequence}\n")
