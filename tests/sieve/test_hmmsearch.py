import os
import tempfile
import unittest
from subprocess import CompletedProcess
from unittest.mock import patch

from sieve.hmmsearch import (
    detected_rows_from_hits,
    parse_domtblout,
    read_ko_thresholds,
    run_hmmsearch,
)


def domtblout_line(sequence, model_name, model_accession, full_score, domain_score, ali_from=7):
    return " ".join([
        sequence,
        "-",
        "100",
        model_name,
        model_accession,
        "200",
        "1e-20",
        str(full_score),
        "0.0",
        "1",
        "1",
        "1e-10",
        "1e-10",
        str(domain_score),
        "0.0",
        "3",
        "50",
        str(ali_from),
        str(ali_from + 47),
        str(ali_from - 1),
        str(ali_from + 48),
        "0.98",
        "description with spaces",
    ])


class TestHmmsearch(unittest.TestCase):

    def write_lines(self, path, lines):
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    def test_parse_domtblout_keeps_full_and_domain_scores(self):
        with tempfile.TemporaryDirectory() as tmpd:
            path = os.path.join(tmpd, "hits.domtblout")
            self.write_lines(path, [
                "# comment",
                domtblout_line("p1", "SOD_Fe_N", "PF00081.28", 120.0, 88.5),
            ])

            hits = parse_domtblout(path)

            self.assertEqual(len(hits), 1)
            self.assertEqual(hits[0].sequence_accession, "p1")
            self.assertEqual(hits[0].model_accession, "PF00081.28")
            self.assertEqual(hits[0].full_bitscore, 120.0)
            self.assertEqual(hits[0].domain_bitscore, 88.5)
            self.assertEqual(hits[0].sequence_start, 7)
            self.assertEqual(hits[0].hmm_start, 3)

    def test_parse_domtblout_uses_model_name_when_accession_is_missing(self):
        with tempfile.TemporaryDirectory() as tmpd:
            path = os.path.join(tmpd, "hits.domtblout")
            self.write_lines(path, [domtblout_line("p1", "K04564", "-", 120.0, 88.5)])

            hits = parse_domtblout(path)

            self.assertEqual(hits[0].model_accession, "K04564")

    def test_read_ko_thresholds_uses_first_three_columns(self):
        with tempfile.TemporaryDirectory() as tmpd:
            path = os.path.join(tmpd, "thresholds.tsv")
            self.write_lines(path, [
                "model\tthreshold\tscore_type\tdefinition",
                "K00001\t345.37\tdomain\talcohol dehydrogenase",
                "K00002\t453.33\tfull\tother definition with spaces",
            ])

            thresholds = read_ko_thresholds(path)

            self.assertEqual(thresholds["K00001"].threshold, 345.37)
            self.assertEqual(thresholds["K00001"].score_type, "domain")
            self.assertEqual(thresholds["K00002"].score_type, "full")

    def test_detected_rows_filter_ko_using_domain_or_full_score_type(self):
        with tempfile.TemporaryDirectory() as tmpd:
            domtblout = os.path.join(tmpd, "hits.domtblout")
            thresholds = os.path.join(tmpd, "thresholds.tsv")
            self.write_lines(domtblout, [
                domtblout_line("p1", "K00001", "-", 500.0, 40.0),
                domtblout_line("p2", "K00002", "-", 500.0, 40.0),
                domtblout_line("p3", "K00003", "-", 40.0, 500.0),
            ])
            self.write_lines(thresholds, [
                "model threshold score_type definition",
                "K00001 100 domain domain threshold",
                "K00002 100 full full threshold",
                "K00003 100 full full threshold",
            ])

            rows = detected_rows_from_hits(
                parse_domtblout(domtblout),
                "KO",
                threshold_by_model=read_ko_thresholds(thresholds),
            )

            self.assertEqual(
                [(row["query_accession"], row["target_accession"], row["bitscore"]) for row in rows],
                [("p2", "K00002", 500.0)],
            )

    def test_detected_rows_error_for_missing_ko_threshold(self):
        with tempfile.TemporaryDirectory() as tmpd:
            domtblout = os.path.join(tmpd, "hits.domtblout")
            self.write_lines(domtblout, [domtblout_line("p1", "K99999", "-", 500.0, 400.0)])

            with self.assertRaisesRegex(ValueError, "Missing KO threshold for K99999"):
                detected_rows_from_hits(parse_domtblout(domtblout), "KO", threshold_by_model={})

    def test_run_hmmsearch_uses_cut_ga_only_when_requested(self):
        with patch("subprocess.run", return_value=CompletedProcess(["hmmsearch"], 0)) as run:
            run_hmmsearch("pfam.hmm", "proteins.faa", "pfam.domtblout", use_cut_ga=True)
            run_hmmsearch("ko.hmm", "proteins.faa", "ko.domtblout", use_cut_ga=False)

        self.assertEqual(
            run.call_args_list[0].args[0],
            ["hmmsearch", "--cut_ga", "--domtblout", "pfam.domtblout", "pfam.hmm", "proteins.faa"],
        )
        self.assertEqual(
            run.call_args_list[1].args[0],
            ["hmmsearch", "--domtblout", "ko.domtblout", "ko.hmm", "proteins.faa"],
        )

    def test_run_hmmsearch_uses_cut_ga_by_default(self):
        with patch("subprocess.run") as run:
            run_hmmsearch("pfam.hmm", "proteins.faa", "pfam.domtblout")
        self.assertEqual(
            run.call_args.args[0],
            ["hmmsearch", "--cut_ga", "--domtblout", "pfam.domtblout", "pfam.hmm", "proteins.faa"],
        )

    def test_run_hmmsearch_passes_cpu_count(self):
        with patch("subprocess.run") as run:
            run_hmmsearch("pfam.hmm", "proteins.faa", "pfam.domtblout", cpus=8)
        self.assertEqual(
            run.call_args.args[0],
            [
                "hmmsearch", "--cut_ga", "--cpu", "8", "--domtblout",
                "pfam.domtblout", "pfam.hmm", "proteins.faa",
            ],
        )

    def test_run_hmmsearch_rejects_invalid_cpu_count(self):
        for cpus in (0, -1, 1.5, "2"):
            with self.subTest(cpus=cpus), self.assertRaisesRegex(ValueError, "cpus"):
                run_hmmsearch("pfam.hmm", "proteins.faa", "pfam.domtblout", cpus=cpus)
