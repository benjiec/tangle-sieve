import csv
import os
import tempfile
import unittest

from tangle.sequence import read_fasta_as_dict, write_fasta_from_dict

from sieve.protein import CuratedProtein
from tests.fixtures import DefaultsFixture
from tests.scripts.helpers import load_script


class TestRuleResultsToFastaScript(unittest.TestCase):

    def setUp(self):
        CuratedProtein.clear_cache()
        self.fx = DefaultsFixture(self)
        self.repo = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    def tearDown(self):
        CuratedProtein.clear_cache()
        self.fx.cleanup()

    def write_artifacts(self, artifacts):
        os.makedirs(artifacts, exist_ok=True)
        rule_name = "Leader().upstreamOfPfam('PF00081').betweenAA(-45, 0).is_mTP()"
        with open(os.path.join(artifacts, "rule-results.tsv"), "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=[
                "protein accession",
                "sequence accession",
                "genome accession",
                "contig accession",
                "pass all",
                rule_name,
                "OtherRule",
            ], delimiter="\t")
            writer.writeheader()
            writer.writerows([
                {
                    "protein accession": "p1",
                    "sequence accession": "seq1",
                    "genome accession": "g1",
                    "contig accession": "c1",
                    "pass all": "true",
                    rule_name: "true",
                    "OtherRule": "yes",
                },
                {
                    "protein accession": "p2",
                    "sequence accession": "seq2",
                    "genome accession": "g2",
                    "contig accession": "c2",
                    "pass all": "false",
                    rule_name: "true",
                    "OtherRule": "yes",
                },
                {
                    "protein accession": "p3",
                    "sequence accession": "seq3",
                    "genome accession": "g3",
                    "contig accession": "c3",
                    "pass all": "true",
                    rule_name: "false",
                    "OtherRule": "yes",
                },
            ])
        write_fasta_from_dict({
            "seq1": "MSEQONE",
            "seq2": "MSEQTWO",
            "seq3": "MSEQTHREE",
        }, os.path.join(artifacts, "sequences.faa"))
        return rule_name

    def write_taxonomy(self):
        self.fx.write_taxonomy_rows([
            {
                "Genome Accession": "g1",
                "Genome Name": "g1 genome",
                "TaxID": "101",
                "Organism": "Example one",
                "Domain": "Eukaryota",
                "Kingdom": "Metazoa",
                "Phylum": "Cnidaria",
                "Class": "Hydrozoa",
                "Order": "Anthoathecata",
                "Family": "Hydractiniidae",
                "Genus": "Hydractinia",
                "Species": "Hydractinia symbiolongicarpus",
            },
            {
                "Genome Accession": "g3",
                "Genome Name": "g3 genome",
                "TaxID": "103",
                "Organism": "Example three",
                "Domain": "Eukaryota",
                "Kingdom": "Metazoa",
                "Phylum": "Arthropoda",
                "Class": "Insecta",
                "Order": "Diptera",
                "Family": "Drosophilidae",
                "Genus": "Drosophila",
                "Species": "Drosophila melanogaster",
            },
        ])

    def test_defaults_to_pass_all_true(self):
        script = load_script(os.path.join(self.repo, "scripts", "rule-results-to-fasta.py"))
        with tempfile.TemporaryDirectory() as tmpd:
            artifacts = os.path.join(tmpd, "artifacts")
            self.write_artifacts(artifacts)
            output = os.path.join(tmpd, "passed.faa")

            script.main(["--artifacts-dir", artifacts, "--output", output])

            self.assertEqual(read_fasta_as_dict(output), {
                "seq1": "MSEQONE",
                "seq3": "MSEQTHREE",
            })

    def test_explicit_pass_all_rule_overrides_default(self):
        script = load_script(os.path.join(self.repo, "scripts", "rule-results-to-fasta.py"))
        with tempfile.TemporaryDirectory() as tmpd:
            artifacts = os.path.join(tmpd, "artifacts")
            self.write_artifacts(artifacts)
            output = os.path.join(tmpd, "failed.faa")

            script.main(["--artifacts-dir", artifacts, "--output", output, "--rule", "pass all=false"])

            self.assertEqual(read_fasta_as_dict(output), {
                "seq2": "MSEQTWO",
            })

    def test_ands_rule_filters_with_punctuated_rule_names(self):
        script = load_script(os.path.join(self.repo, "scripts", "rule-results-to-fasta.py"))
        with tempfile.TemporaryDirectory() as tmpd:
            artifacts = os.path.join(tmpd, "artifacts")
            rule_name = self.write_artifacts(artifacts)
            output = os.path.join(tmpd, "filtered.faa")

            script.main([
                "--artifacts-dir", artifacts,
                "--output", output,
                "--rule", f"{rule_name}=true",
                "--rule", "OtherRule=yes",
            ])

            self.assertEqual(read_fasta_as_dict(output), {
                "seq1": "MSEQONE",
            })

    def test_filters_by_any_taxonomy_rank(self):
        script = load_script(os.path.join(self.repo, "scripts", "rule-results-to-fasta.py"))
        self.write_taxonomy()
        with tempfile.TemporaryDirectory() as tmpd:
            artifacts = os.path.join(tmpd, "artifacts")
            self.write_artifacts(artifacts)
            output = os.path.join(tmpd, "cnidaria.faa")

            script.main(["--artifacts-dir", artifacts, "--output", output, "--taxon", "cnidaria"])

            self.assertEqual(read_fasta_as_dict(output), {
                "seq1": "MSEQONE",
            })

