import csv
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass

from sieve.protein import CuratedProtein
from tangle.sequence import write_fasta_from_dict


RULE_TRUE = "true"
RULE_FALSE = "false"
RULE_MAYBE = "maybe"
RULE_ERROR = "error"


def _rule_bool(value):
    return RULE_TRUE if value else RULE_FALSE


def _merge_and(values):
    if RULE_FALSE in values:
        return RULE_FALSE
    if RULE_ERROR in values:
        return RULE_ERROR
    if RULE_MAYBE in values:
        return RULE_MAYBE
    return RULE_TRUE


def _merge_or(values):
    if RULE_TRUE in values:
        return RULE_TRUE
    if RULE_MAYBE in values:
        return RULE_MAYBE
    if RULE_ERROR in values:
        return RULE_ERROR
    return RULE_FALSE


def _result_counts(values):
    counts = {RULE_TRUE: 0, RULE_FALSE: 0, RULE_MAYBE: 0, RULE_ERROR: 0}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    return counts


def _format_result_counts(counts):
    return " ".join([
        f"{key}={counts[key]}"
        for key in (RULE_TRUE, RULE_FALSE, RULE_MAYBE, RULE_ERROR)
        if counts.get(key, 0)
    ]) or "no results"


def _target_prefix(value):
    return str(value).split(".", 1)[0]


def _motif_matches(feature, motif):
    return feature == motif or feature.startswith(motif + ".")


def _edge_distance(a_start, a_end, b_start, b_end):
    if a_start <= b_end and b_start <= a_end:
        return 0
    if a_end < b_start:
        return b_start - a_end
    return a_start - b_end


def _safe_filename(value):
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
    return value[:120] or "rule"


def _leader_sequence_id(context):
    return f"{context.protein.protein_accession}_with_leader"


def _rule_artifacts_dir(artifacts_dir, rule):
    if artifacts_dir is None:
        return None
    return os.path.join(artifacts_dir, _safe_filename(rule.label))


def _run_command(cmd, artifacts_dir=None):
    if artifacts_dir is not None:
        os.makedirs(artifacts_dir, exist_ok=True)
        with open(os.path.join(artifacts_dir, "command.txt"), "w", encoding="utf-8") as f:
            f.write(" ".join(cmd) + "\n")
    completed = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if artifacts_dir is not None:
        with open(os.path.join(artifacts_dir, "stdout.txt"), "w", encoding="utf-8") as f:
            f.write(completed.stdout)
        with open(os.path.join(artifacts_dir, "stderr.txt"), "w", encoding="utf-8") as f:
            f.write(completed.stderr)
    return completed


class Rule(object):

    def __and__(self, other):
        return AndRule(self, _as_rule(other))

    def __or__(self, other):
        return OrRule(self, _as_rule(other))

    def evaluate(self, context):
        raise NotImplementedError

    def evaluate_many(self, contexts, artifacts_dir=None):
        results = {}
        for context in contexts:
            try:
                results[context.key] = self.evaluate(context)
            except Exception as e:
                print(f"{self.label} failed for {context.key}: {e}", file=sys.stderr)
                results[context.key] = RULE_ERROR
        return results

    def annotation_columns(self):
        return []

    def annotations_many(self, contexts, rule_results):
        return {}

    def atomic_rules(self):
        return [self]

    def resolve(self, context, atomic_results):
        return atomic_results[self.label][context.key]


def _as_rule(value):
    if not isinstance(value, Rule):
        raise TypeError(f"Expected Rule, got {type(value)}")
    return value


class CompositeRule(Rule):

    def __init__(self, *rules):
        self.rules = [_as_rule(rule) for rule in rules]

    def atomic_rules(self):
        rules = []
        seen = set()
        for rule in self.rules:
            for atomic in rule.atomic_rules():
                if atomic.label not in seen:
                    rules.append(atomic)
                    seen.add(atomic.label)
        return rules


class AndRule(CompositeRule):

    @property
    def label(self):
        return " & ".join([rule.label for rule in self.rules])

    def evaluate(self, context):
        return _merge_and([rule.evaluate(context) for rule in self.rules])

    def resolve(self, context, atomic_results):
        return _merge_and([rule.resolve(context, atomic_results) for rule in self.rules])


class OrRule(CompositeRule):

    @property
    def label(self):
        return " | ".join([rule.label for rule in self.rules])

    def evaluate(self, context):
        return _merge_or([rule.evaluate(context) for rule in self.rules])

    def resolve(self, context, atomic_results):
        return _merge_or([rule.resolve(context, atomic_results) for rule in self.rules])


class RuleContext(object):

    def __init__(self, protein):
        self.protein = protein
        self.key = (protein.protein_accession, protein.genome_accession)
        self._hmm_alignments = {}

    def hmm_alignment(self, profile):
        if profile not in self._hmm_alignments:
            self._hmm_alignments[profile] = self.protein.hmm_align(profile)
        return self._hmm_alignments[profile]


class Rules(object):

    def __init__(self, rule):
        self.rule = _as_rule(rule)

    def atomic_rules(self):
        return self.rule.atomic_rules()

    def headers(self):
        atomic_rules = self.atomic_rules()
        annotation_columns = []
        seen = set()
        for rule in atomic_rules:
            for column in rule.annotation_columns():
                if column not in seen:
                    annotation_columns.append(column)
                    seen.add(column)
        return (
            ["protein accession", "genome accession", "pass all"]
            + [rule.label for rule in atomic_rules]
            + annotation_columns
        )

    def write_rows(self, output_tsv, rows):
        with open(output_tsv, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=self.headers(), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def check(self, protein_keys, output_tsv, artifacts_dir=None, trace=True):
        contexts = [
            RuleContext(CuratedProtein(protein_accession, genome_accession))
            for protein_accession, genome_accession in protein_keys
        ]
        atomic_rules = self.atomic_rules()
        if artifacts_dir is not None:
            artifacts_dir = os.path.abspath(artifacts_dir)
            os.makedirs(artifacts_dir, exist_ok=True)
        atomic_results = {}
        total_rules = len(atomic_rules)
        total_proteins = len(contexts)
        for i, rule in enumerate(atomic_rules, start=1):
            if trace:
                print(f"[rules {i}/{total_rules}] {rule.label}: {total_proteins} proteins", file=sys.stderr)
            rule_results = rule.evaluate_many(contexts, _rule_artifacts_dir(artifacts_dir, rule))
            atomic_results[rule.label] = rule_results
            if trace:
                counts = _result_counts(rule_results.values())
                print(f"[rules {i}/{total_rules}] done: {_format_result_counts(counts)}", file=sys.stderr)

        atomic_annotations = {}
        for rule in atomic_rules:
            annotations = rule.annotations_many(contexts, atomic_results[rule.label])
            for context_key, context_annotations in annotations.items():
                atomic_annotations.setdefault(context_key, {}).update(context_annotations)

        rows = []
        for context in contexts:
            row = {
                "protein accession": context.protein.protein_accession,
                "genome accession": context.protein.genome_accession,
                "pass all": self.rule.resolve(context, atomic_results),
            }
            for rule in atomic_rules:
                row[rule.label] = atomic_results[rule.label][context.key]
            row.update(atomic_annotations.get(context.key, {}))
            rows.append(row)

        self.write_rows(output_tsv, rows)
        return rows


class DetectedTargetRule(Rule):

    def __init__(self, label, row_getter, accession, prefix_match=False):
        self.label = label
        self.row_getter = row_getter
        self.accession = accession
        self.prefix_match = prefix_match

    def evaluate(self, context):
        for row in self.row_getter(context.protein):
            target = row["target_accession"]
            if self.prefix_match:
                target = _target_prefix(target)
            if target == self.accession:
                return RULE_TRUE
        return RULE_FALSE


class Pfam(object):

    @staticmethod
    def matches(accession):
        return DetectedTargetRule(
            label=f"Pfam.matches('{accession}')",
            row_getter=lambda protein: protein.detected_pfam(),
            accession=accession,
            prefix_match=True,
        )


class KO(object):

    @staticmethod
    def matches(accession):
        return DetectedTargetRule(
            label=f"KO.matches('{accession}')",
            row_getter=lambda protein: protein.detected_ko(),
            accession=accession,
        )


class HMMAlignment(object):

    def __init__(self, profile):
        self.profile = profile

    def is_at(self, expected, hmm_pos):
        return HMMPositionRule(self.profile, expected, hmm_pos)

    def covers(self, start, end):
        return HMMCoverageRule(self.profile, start, end)


class HMMPositionRule(Rule):

    def __init__(self, profile, expected, hmm_pos):
        self.profile = profile
        self.expected = expected
        self.hmm_pos = hmm_pos
        self.label = f"HMMAlignment('{os.path.basename(profile)}').is_at('{expected}', {hmm_pos})"

    def evaluate(self, context):
        alignment = context.hmm_alignment(self.profile)
        for offset, expected_aa in enumerate(self.expected):
            aa = alignment.aa_at_hmm_pos_1b(self.hmm_pos + offset)
            if aa is None or aa[1] != expected_aa:
                return RULE_FALSE
        return RULE_TRUE


class HMMCoverageRule(Rule):

    def __init__(self, profile, start, end):
        self.profile = profile
        self.start = start
        self.end = end
        self.label = f"HMMAlignment('{os.path.basename(profile)}').covers({start}, {end})"

    def evaluate(self, context):
        alignment = context.hmm_alignment(self.profile)
        for hmm_pos in range(self.start, self.end + 1):
            if alignment.aa_at_hmm_pos_1b(hmm_pos) is None:
                return RULE_FALSE
        return RULE_TRUE


class Leader(object):

    @staticmethod
    def is_mTP():
        return LeaderRule("mTP")

    @staticmethod
    def is_SP():
        return LeaderRule("SP")

    @staticmethod
    def is_noTP():
        return LeaderRule("noTP")


class LeaderRule(Rule):

    def __init__(self, prediction):
        self.prediction = prediction
        self.label = f"Leader.is_{prediction}()"
        self._predictions_by_key = {}

    def annotation_columns(self):
        return ["Leader.call"]

    def evaluate_many(self, contexts, artifacts_dir=None):
        sequence_ids = {_leader_sequence_id(context): context for context in contexts}
        results = {context.key: RULE_ERROR for context in contexts}
        self._predictions_by_key = {}
        try:
            with tempfile.TemporaryDirectory() as tmpd:
                working_dir = artifacts_dir if artifacts_dir is not None else tmpd
                if artifacts_dir is not None:
                    os.makedirs(artifacts_dir, exist_ok=True)
                fasta_path = os.path.join(working_dir, "query.faa")
                fasta = {
                    sequence_id: context.protein.sequence_with_leader()
                    for sequence_id, context in sequence_ids.items()
                }
                write_fasta_from_dict(fasta, fasta_path)
                cmd = [
                    "docker", "run", "--rm", "--platform", "linux/amd64",
                    "-v", f"{working_dir}:/data",
                    "local-targetp:2.0",
                    "-fasta", "/data/query.faa",
                    "-org", "non-pl",
                    "-format", "short",
                    "-stdout",
                ]
                completed = _run_command(cmd, artifacts_dir)
                if completed.returncode != 0:
                    raise subprocess.CalledProcessError(
                        completed.returncode,
                        cmd,
                        output=completed.stdout,
                        stderr=completed.stderr,
                    )
                predictions = _parse_targetp_output(completed.stdout)
        except Exception as e:
            print(f"{self.label} batch failed: {e}", file=sys.stderr)
            return results

        for sequence_id, context in sequence_ids.items():
            prediction = predictions.get(sequence_id)
            if prediction is None:
                print(f"{self.label} missing TargetP row for {context.key}", file=sys.stderr)
                results[context.key] = RULE_ERROR
            else:
                self._predictions_by_key[context.key] = prediction
                results[context.key] = _rule_bool(prediction == self.prediction)
        return results

    def annotations_many(self, contexts, rule_results):
        return {
            context.key: {"Leader.call": self._predictions_by_key[context.key]}
            for context in contexts
            if context.key in self._predictions_by_key
        }


def _parse_targetp_output(text):
    predictions = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            predictions[parts[0]] = parts[1]
    return predictions


class TFMotifs(object):

    @staticmethod
    def has_within(distance, motif_a, motif_b, min_score_threshold=8):
        return TFMotifWithinRule(distance, motif_a, motif_b, min_score_threshold)


class TFMotifWithinRule(Rule):

    def __init__(self, distance, motif_a, motif_b, min_score_threshold, scope=None):
        self.distance = distance
        self.motif_a = motif_a
        self.motif_b = motif_b
        self.min_score_threshold = min_score_threshold
        self.scope = scope
        self.label = self._label()

    def _label(self):
        base = (
            f"TFMotifs.has_within({self.distance}, '{self.motif_a}', '{self.motif_b}', "
            f"min_score_threshold={self.min_score_threshold})"
        )
        if self.scope is None:
            return base
        scope_type = self.scope[0]
        if scope_type == "intron":
            intron_number = self.scope[1]
            if intron_number is None:
                return f"{base}.in_intron()"
            return f"{base}.in_intron({intron_number})"
        if scope_type == "exon":
            exon_number = self.scope[1]
            if exon_number is None:
                return f"{base}.in_exon()"
            return f"{base}.in_exon({exon_number})"
        if scope_type == "between":
            return f"{base}.between({self.scope[1]}, {self.scope[2]})"
        raise ValueError(f"Unsupported TFMotif scope: {self.scope}")

    def in_intron(self, intron_number=None):
        self._require_unscoped()
        return TFMotifWithinRule(
            self.distance,
            self.motif_a,
            self.motif_b,
            self.min_score_threshold,
            scope=("intron", intron_number),
        )

    def in_exon(self, exon_number=None):
        self._require_unscoped()
        return TFMotifWithinRule(
            self.distance,
            self.motif_a,
            self.motif_b,
            self.min_score_threshold,
            scope=("exon", exon_number),
        )

    def between(self, start, end):
        self._require_unscoped()
        return TFMotifWithinRule(
            self.distance,
            self.motif_a,
            self.motif_b,
            self.min_score_threshold,
            scope=("between", start, end),
        )

    def _require_unscoped(self):
        if self.scope is not None:
            raise ValueError("TFMotifs rules can only be scoped once")

    def evaluate_many(self, contexts, artifacts_dir=None):
        sequence_ids = {f"seq{i}": context for i, context in enumerate(contexts)}
        results = {context.key: RULE_ERROR for context in contexts}
        loci = {}
        try:
            for sequence_id, context in sequence_ids.items():
                loci[sequence_id] = context.protein.genomic_locus_with_leader()
            with tempfile.TemporaryDirectory() as tmpd:
                working_dir = artifacts_dir if artifacts_dir is not None else tmpd
                if artifacts_dir is not None:
                    os.makedirs(artifacts_dir, exist_ok=True)
                fasta_path = os.path.join(working_dir, "locus.fna")
                fasta = {
                    sequence_id: locus.sequence()
                    for sequence_id, locus in loci.items()
                }
                write_fasta_from_dict(fasta, fasta_path)
                cmd = ["gimme", "scan", fasta_path, "-b", "-c", "0.85"]
                completed = _run_command(cmd, artifacts_dir)
                if completed.returncode != 0:
                    raise subprocess.CalledProcessError(
                        completed.returncode,
                        cmd,
                        output=completed.stdout,
                        stderr=completed.stderr,
                    )
                hits_by_sequence = _parse_gimme_scan_output(completed.stdout)
        except Exception as e:
            print(f"{self.label} batch failed: {e}", file=sys.stderr)
            return results

        for sequence_id, context in sequence_ids.items():
            try:
                locus = loci[sequence_id]
                hits = hits_by_sequence.get(sequence_id, [])
                results[context.key] = self._evaluate_locus(locus, hits)
            except Exception as e:
                print(f"{self.label} failed for {context.key}: {e}", file=sys.stderr)
                results[context.key] = RULE_ERROR
        return results

    def _evaluate_locus(self, locus, hits):
        intervals = _scope_intervals(locus, self.scope)
        if not intervals:
            return RULE_FALSE
        motif_a_hits = []
        motif_b_hits = []
        for hit in hits:
            hit_start, hit_end = hit.normalized_interval()
            if hit.score < self.min_score_threshold:
                continue
            if not _hit_in_any_interval(hit_start, hit_end, intervals):
                continue
            if _motif_matches(hit.feature, self.motif_a):
                motif_a_hits.append((hit_start, hit_end))
            if _motif_matches(hit.feature, self.motif_b):
                motif_b_hits.append((hit_start, hit_end))
        for a_start, a_end in motif_a_hits:
            for b_start, b_end in motif_b_hits:
                if _edge_distance(a_start, a_end, b_start, b_end) <= self.distance:
                    return RULE_TRUE
        return RULE_MAYBE


@dataclass
class GimmeHit:
    sequence: str
    start: int
    end: int
    feature: str
    score: float
    strand: str

    def normalized_interval(self):
        return (min(self.start, self.end), max(self.start, self.end))


def _intron_interval(locus, intron_number):
    if intron_number < 1:
        raise ValueError("intron_number must be 1 or greater")
    if len(locus.cds_intervals_1b) <= intron_number:
        return None
    left_cds = locus.cds_intervals_1b[intron_number - 1]
    right_cds = locus.cds_intervals_1b[intron_number]
    start = left_cds[1] + 1
    end = right_cds[0] - 1
    if start > end:
        return None
    return (start, end)


def _all_intron_intervals(locus):
    return [
        interval
        for interval in (_intron_interval(locus, i) for i in range(1, len(locus.cds_intervals_1b)))
        if interval is not None
    ]


def _exon_interval(locus, exon_number):
    if exon_number < 1:
        raise ValueError("exon_number must be 1 or greater")
    index = exon_number - 1
    if index >= len(locus.cds_intervals_1b):
        return None
    return locus.cds_intervals_1b[index]


def _scope_intervals(locus, scope):
    if scope is None:
        return list(locus.cds_intervals_1b) + _all_intron_intervals(locus)
    scope_type = scope[0]
    if scope_type == "intron":
        intron_number = scope[1]
        if intron_number is None:
            return _all_intron_intervals(locus)
        interval = _intron_interval(locus, intron_number)
        return [] if interval is None else [interval]
    if scope_type == "exon":
        exon_number = scope[1]
        if exon_number is None:
            return list(locus.cds_intervals_1b)
        interval = _exon_interval(locus, exon_number)
        return [] if interval is None else [interval]
    if scope_type == "between":
        return [(min(scope[1], scope[2]), max(scope[1], scope[2]))]
    raise ValueError(f"Unsupported TFMotif scope: {scope}")


def _hit_in_any_interval(hit_start, hit_end, intervals):
    return any(hit_start >= start and hit_end <= end for start, end in intervals)


def _parse_gimme_scan_output(text):
    hits_by_sequence = {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if parts[0] == "sequence":
            continue
        if len(parts) < 6:
            continue
        hit = GimmeHit(
            sequence=parts[0],
            start=int(parts[1]),
            end=int(parts[2]),
            feature=parts[3],
            score=float(parts[4]),
            strand=parts[5],
        )
        hits_by_sequence.setdefault(hit.sequence, []).append(hit)
    return hits_by_sequence
