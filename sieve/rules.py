import csv
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass

from sieve.protein import CuratedProtein, LeaderSequenceCandidate
from tangle.sequence import write_fasta_from_dict


RULE_TRUE = "true"
RULE_FALSE = "false"
RULE_MAYBE = "maybe"
RULE_ERROR = "error"
RULE_YES = "yes"
RULE_TOO_FAR = "too_far"

RULE_PASSING = {RULE_TRUE, RULE_YES, RULE_TOO_FAR}
RULE_STANDARD_ORDER = [RULE_TRUE, RULE_FALSE, RULE_MAYBE, RULE_ERROR, RULE_YES, RULE_TOO_FAR]
TARGETP_PROBABILITY_THRESHOLD = 0.7
TARGETP_PROBABILITY_ORDER = ["noTP", "mTP", "SP"]
TARGETP_PROBABILITY_COLUMNS = {
    prediction: f"Leader.call('{prediction}')"
    for prediction in TARGETP_PROBABILITY_ORDER
}
DEEPLOC_PROTEIN_ID_COLUMN = "Protein_ID"
DEEPLOC_LOCALIZATIONS_COLUMN = "Localizations"
DEEPLOC_SIGNALS_COLUMN = "Signals"
DEEPLOC_MTP_SIGNAL = "Mitochondrial transit peptide"
DEEPLOC_SP_SIGNAL = "Signal peptide"
DEEPLOC_SIGNAL_PREDICTIONS = ["mTP", "SP"]
LEADER_LOCALIZATION_COLUMN = "Leader.localization"


def _rule_bool(value):
    return RULE_TRUE if value else RULE_FALSE


def _rule_passes(value):
    return value in RULE_PASSING


def _merge_and(values):
    if RULE_FALSE in values:
        return RULE_FALSE
    if RULE_ERROR in values:
        return RULE_ERROR
    if RULE_MAYBE in values:
        return RULE_MAYBE
    if any(not _rule_passes(value) for value in values):
        return RULE_FALSE
    return RULE_TRUE


def _merge_or(values):
    if any(_rule_passes(value) for value in values):
        return RULE_TRUE
    if RULE_MAYBE in values:
        return RULE_MAYBE
    if RULE_ERROR in values:
        return RULE_ERROR
    return RULE_FALSE


def _normalize_pass_all(value):
    if _rule_passes(value):
        return RULE_TRUE
    if value in (RULE_FALSE, RULE_MAYBE, RULE_ERROR):
        return value
    return RULE_FALSE


def _result_counts(values):
    counts = {RULE_TRUE: 0, RULE_FALSE: 0, RULE_MAYBE: 0, RULE_ERROR: 0}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    return counts


def _format_result_counts(counts):
    ordered_keys = [key for key in RULE_STANDARD_ORDER if counts.get(key, 0)]
    ordered_keys.extend(sorted(key for key in counts if key not in RULE_STANDARD_ORDER))
    return " ".join([f"{key}={counts[key]}" for key in ordered_keys]) or "no results"


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


def _safe_result_value(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "value"


def _locus_sequence_id(context):
    return f"{context.protein.protein_accession}_locus"


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

    def evaluate_many(self, contexts, artifacts_dir=None, **kwargs):
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

    def sequence_result(self, row, candidate):
        return row[self.label]

    def atomic_rules(self):
        return [self]

    def uses_leader_candidates(self):
        return False

    def uses_deeploc(self):
        return False

    def resolve(self, context, atomic_results):
        return atomic_results[self.label][context.key]

    def resolve_row(self, row):
        return row[self.label]

    def filter_sequence_candidates(self, protein, candidates, row):
        return candidates

    def scope_sequence_candidates(self, protein, candidates, row):
        return candidates


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

    def uses_leader_candidates(self):
        return any(rule.uses_leader_candidates() for rule in self.rules)

    def uses_deeploc(self):
        return any(rule.uses_deeploc() for rule in self.rules)


class AndRule(CompositeRule):

    @property
    def label(self):
        return " & ".join([rule.label for rule in self.rules])

    def evaluate(self, context):
        return _merge_and([rule.evaluate(context) for rule in self.rules])

    def resolve(self, context, atomic_results):
        return _merge_and([rule.resolve(context, atomic_results) for rule in self.rules])

    def resolve_row(self, row):
        return _merge_and([rule.resolve_row(row) for rule in self.rules])

    def filter_sequence_candidates(self, protein, candidates, row):
        for rule in self.rules:
            candidates = rule.filter_sequence_candidates(protein, candidates, row)
        return candidates

    def scope_sequence_candidates(self, protein, candidates, row):
        for rule in self.rules:
            candidates = rule.scope_sequence_candidates(protein, candidates, row)
        return candidates


class OrRule(CompositeRule):

    @property
    def label(self):
        return " | ".join([rule.label for rule in self.rules])

    def evaluate(self, context):
        return _merge_or([rule.evaluate(context) for rule in self.rules])

    def resolve(self, context, atomic_results):
        return _merge_or([rule.resolve(context, atomic_results) for rule in self.rules])

    def resolve_row(self, row):
        return _merge_or([rule.resolve_row(row) for rule in self.rules])

    def filter_sequence_candidates(self, protein, candidates, row):
        filtered = {}
        for rule in self.rules:
            if _rule_passes(rule.resolve_row(row)):
                for candidate in rule.filter_sequence_candidates(protein, candidates, row):
                    filtered[candidate.accession] = candidate
        return list(filtered.values())

    def scope_sequence_candidates(self, protein, candidates, row):
        scoped = {}
        rules = [rule for rule in self.rules if rule.uses_leader_candidates()]
        if not rules:
            rules = self.rules
        for rule in rules:
            for candidate in rule.scope_sequence_candidates(protein, candidates, row):
                scoped[candidate.accession] = candidate
        return list(scoped.values())


class RuleContext(object):

    def __init__(self, protein, hmm_profiles=None):
        self.protein = protein
        self.key = (protein.protein_accession, protein.genome_accession)
        self._hmm_alignments = {}
        self._hmm_profiles = _hmm_profiles_by_basename(hmm_profiles or [])

    def _resolve_hmm_profile(self, profile):
        if os.path.basename(profile) != profile:
            return profile
        return self._hmm_profiles.get(profile, profile)

    def hmm_alignment(self, profile):
        resolved_profile = self._resolve_hmm_profile(profile)
        if resolved_profile not in self._hmm_alignments:
            self._hmm_alignments[resolved_profile] = self.protein.hmm_align(resolved_profile)
        return self._hmm_alignments[resolved_profile]


def _hmm_profiles_by_basename(profiles):
    by_basename = {}
    for profile in profiles:
        profile = os.path.abspath(profile)
        basename = os.path.basename(profile)
        if basename in by_basename and by_basename[basename] != profile:
            raise ValueError(
                f"Ambiguous HMM profile basename {basename!r}: "
                f"{by_basename[basename]!r} and {profile!r}"
            )
        by_basename[basename] = profile
    return by_basename


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
            ["protein accession", "sequence accession", "genome accession", "contig accession", "pass all"]
            + [rule.label for rule in atomic_rules]
            + annotation_columns
        )

    def write_rows(self, output_tsv, rows):
        with open(output_tsv, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=self.headers(), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def sequence_candidates_for_row(self, row):
        protein = CuratedProtein(row["protein accession"], row["genome accession"])
        return self.sequence_candidates_for_protein_row(protein, row)

    def sequence_candidates_for_protein_row(self, protein, row):
        candidates = self._row_sequence_candidates(protein, row)
        return self.rule.filter_sequence_candidates(protein, candidates, row)

    def scoped_sequence_candidates_for_row(self, row):
        protein = CuratedProtein(row["protein accession"], row["genome accession"])
        return self.scoped_sequence_candidates_for_protein_row(protein, row)

    def scoped_sequence_candidates_for_protein_row(self, protein, row):
        return self._row_sequence_candidates(protein, row)

    def uses_leader_candidates(self):
        return self.rule.uses_leader_candidates()

    def uses_deeploc(self):
        return self.rule.uses_deeploc()

    def _scoped_sequence_candidates(self, protein, row):
        if not self.uses_leader_candidates():
            return _original_sequence_candidate(protein)
        candidates = protein.sequences_with_leader()
        return _dedupe_sequence_candidates(self.rule.scope_sequence_candidates(protein, candidates, row))

    def _row_sequence_candidates(self, protein, row):
        candidates = self._scoped_sequence_candidates(protein, row)
        sequence_accession = row.get("sequence accession")
        if not sequence_accession:
            return candidates
        return [
            candidate for candidate in candidates
            if candidate.accession == sequence_accession
        ]

    def check(self, protein_keys, output_tsv, artifacts_dir=None, trace=True, deeploc_csv=None,
              hmm_profiles=None):
        contexts = [
            RuleContext(
                CuratedProtein(protein_accession, genome_accession),
                hmm_profiles=hmm_profiles,
            )
            for protein_accession, genome_accession in protein_keys
        ]
        return self.check_proteins(
            contexts,
            output_tsv,
            artifacts_dir=artifacts_dir,
            trace=trace,
            deeploc_csv=deeploc_csv,
        )

    def check_proteins(self, proteins, output_tsv, artifacts_dir=None, trace=True, deeploc_csv=None,
                       hmm_profiles=None):
        if self.uses_deeploc() and deeploc_csv is None:
            raise ValueError("--deeploc-csv is required for DeepLoc-backed Leader rules")
        contexts = [
            protein if isinstance(protein, RuleContext) else RuleContext(
                protein,
                hmm_profiles=hmm_profiles,
            )
            for protein in proteins
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
            rule_results = rule.evaluate_many(
                contexts,
                _rule_artifacts_dir(artifacts_dir, rule),
                deeploc_csv=deeploc_csv,
            )
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
                "contig accession": self._contig_accession(context),
            }
            for rule in atomic_rules:
                row[rule.label] = atomic_results[rule.label][context.key]
            row.update(atomic_annotations.get(context.key, {}))
            for candidate in self._scoped_sequence_candidates(context.protein, row):
                candidate_row = row.copy()
                candidate_row["sequence accession"] = candidate.accession
                _apply_sequence_annotations(candidate_row, candidate)
                for rule in atomic_rules:
                    candidate_row[rule.label] = rule.sequence_result(candidate_row, candidate)
                candidate_row["pass all"] = _normalize_pass_all(self.rule.resolve_row(candidate_row))
                rows.append(candidate_row)

        self.write_rows(output_tsv, rows)
        return rows

    def _contig_accession(self, context):
        try:
            return context.protein.genomic_locus_with_leader().contig_accession
        except Exception:
            return ""


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

    def matches_regex(self, pattern, hmm_pos):
        return HMMRegexRule(self.profile, pattern, hmm_pos)

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


class HMMRegexRule(Rule):

    def __init__(self, profile, pattern, hmm_pos):
        self.profile = profile
        self.pattern = pattern
        self.hmm_pos = hmm_pos
        self.label = f"HMMAlignment('{os.path.basename(profile)}').matches_regex('{pattern}', {hmm_pos})"

    def evaluate(self, context):
        alignment = context.hmm_alignment(self.profile)
        aligned = []
        hmm_pos = self.hmm_pos
        while True:
            aa = alignment.aa_at_hmm_pos_1b(hmm_pos)
            if aa is None:
                if hmm_pos not in alignment.hmm_positions_1b():
                    break
                aligned.append("-")
            else:
                aligned.append(aa[1])
            hmm_pos += 1
        if not aligned:
            return RULE_FALSE
        return _rule_bool(re.match(self.pattern, "".join(aligned)) is not None)


class Leader(object):

    def __init__(self, window_start=-30, window_end=3, pfam_accession=None):
        self.window_start = window_start
        self.window_end = window_end
        self.pfam_accession = pfam_accession

    def upstreamOfPfam(self, accession):
        return Leader(self.window_start, self.window_end, accession)

    def betweenAA(self, start, end):
        return Leader(start, end, self.pfam_accession)

    def is_mTP(self, deeploc=False):
        source = "deeploc" if deeploc else "targetp"
        return LeaderRule("mTP", self.window_start, self.window_end, self.pfam_accession, source=source)

    def is_SP(self, deeploc=False):
        source = "deeploc" if deeploc else "targetp"
        return LeaderRule("SP", self.window_start, self.window_end, self.pfam_accession, source=source)

    def is_noTP(self):
        return LeaderRule("noTP", self.window_start, self.window_end, self.pfam_accession)

    def localize_at(self, localization):
        return LeaderRule(
            localization,
            self.window_start,
            self.window_end,
            self.pfam_accession,
            source="deeploc",
            call_type="localization",
        )


class LeaderRule(Rule):

    def __init__(
        self,
        prediction,
        window_start=-30,
        window_end=3,
        pfam_accession=None,
        source="targetp",
        call_type="signal",
    ):
        self.prediction = prediction
        self.window_start = window_start
        self.window_end = window_end
        self.pfam_accession = pfam_accession
        self.source = source
        self.call_type = call_type
        self.label = self._label()
        self._predictions_by_key = {}
        self._annotation_columns = []

    def _label(self):
        base = "Leader()"
        if self.pfam_accession is not None:
            base += f".upstreamOfPfam('{self.pfam_accession}')"
        base += f".betweenAA({self.window_start}, {self.window_end})"
        if self.call_type == "localization":
            return f"{base}.localize_at('{self.prediction}')"
        if self.source == "deeploc":
            return f"{base}.is_{self.prediction}(deeploc=True)"
        return f"{base}.is_{self.prediction}()"

    def annotation_columns(self):
        if self.source == "targetp":
            return [TARGETP_PROBABILITY_COLUMNS[prediction] for prediction in TARGETP_PROBABILITY_ORDER]
        return self._annotation_columns

    def uses_leader_candidates(self):
        return True

    def uses_deeploc(self):
        return self.source == "deeploc"

    def sequence_result(self, row, candidate):
        if row[self.label] == RULE_ERROR:
            return RULE_ERROR
        if self._leader_call_matches(row):
            return RULE_TRUE
        return RULE_FALSE

    def evaluate_many(self, contexts, artifacts_dir=None, deeploc_csv=None, **kwargs):
        sequence_contexts = {}
        sequence_candidates = {}
        candidate_counts_by_key = {}
        for context in contexts:
            candidates = self._scoped_candidates(context.protein)
            candidate_counts_by_key[context.key] = len(candidates)
            for candidate in candidates:
                sequence_contexts[candidate.accession] = context
                sequence_candidates[candidate.accession] = candidate
        results = {
            context.key: RULE_ERROR if candidate_counts_by_key[context.key] else RULE_FALSE
            for context in contexts
        }
        self._predictions_by_key = {}
        if not sequence_candidates:
            return results
        if self.source == "deeploc":
            return self._evaluate_many_deeploc(
                contexts,
                sequence_contexts,
                sequence_candidates,
                candidate_counts_by_key,
                results,
                deeploc_csv,
            )
        try:
            with tempfile.TemporaryDirectory() as tmpd:
                working_dir = artifacts_dir if artifacts_dir is not None else tmpd
                if artifacts_dir is not None:
                    os.makedirs(artifacts_dir, exist_ok=True)
                fasta_path = os.path.join(working_dir, "query.faa")
                fasta = {
                    sequence_id: candidate.sequence
                    for sequence_id, candidate in sequence_candidates.items()
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

        calls_by_key = {context.key: [] for context in contexts}
        for sequence_id, context in sequence_contexts.items():
            call = predictions.get(sequence_id)
            if call is None:
                print(f"{self.label} missing TargetP row for {context.key}", file=sys.stderr)
                results[context.key] = RULE_ERROR
            else:
                candidate = sequence_candidates[sequence_id]
                calls_by_key[context.key].append((candidate.start_label, call))

        for context in contexts:
            if len(calls_by_key[context.key]) != candidate_counts_by_key[context.key]:
                continue
            self._predictions_by_key[context.key] = calls_by_key[context.key]
            results[context.key] = _rule_bool(
                any(call.prediction == self.prediction for _start, call in calls_by_key[context.key])
            )
        return results

    def _evaluate_many_deeploc(
        self,
        contexts,
        sequence_contexts,
        sequence_candidates,
        candidate_counts_by_key,
        results,
        deeploc_csv,
    ):
        try:
            predictions = _parse_deeploc_csv(deeploc_csv)
        except Exception as e:
            print(f"{self.label} DeepLoc parse failed: {e}", file=sys.stderr)
            return results

        self._annotation_columns = _deeploc_annotation_columns(predictions.values())
        calls_by_key = {context.key: [] for context in contexts}
        for sequence_id, context in sequence_contexts.items():
            call = predictions.get(sequence_id)
            if call is None:
                print(f"{self.label} missing DeepLoc row for {context.key}: {sequence_id}", file=sys.stderr)
                results[context.key] = RULE_ERROR
            else:
                candidate = sequence_candidates[sequence_id]
                calls_by_key[context.key].append((candidate.start_label, call))

        for context in contexts:
            if len(calls_by_key[context.key]) != candidate_counts_by_key[context.key]:
                continue
            self._predictions_by_key[context.key] = calls_by_key[context.key]
            results[context.key] = _rule_bool(
                any(self._call_matches_prediction(call) for _start, call in calls_by_key[context.key])
            )
        return results

    def annotations_many(self, contexts, rule_results):
        return {
            context.key: {"_Leader.calls_by_start": _format_leader_calls(self._predictions_by_key[context.key])}
            for context in contexts
            if context.key in self._predictions_by_key
        }

    def filter_sequence_candidates(self, protein, candidates, row):
        if not _rule_passes(row[self.label]):
            return []
        if "sequence accession" in row:
            if not self._leader_call_matches(row):
                return []
            return candidates
        matching_starts = {
            start
            for start, call in _parse_leader_calls(row)
            if self._columns_match_prediction(call)
        }
        if self.pfam_accession is not None:
            candidates = self._scoped_candidates(protein)
        return [
            candidate for candidate in candidates
            if candidate.start_label in matching_starts
        ]

    def _leader_call_matches(self, row):
        if self.call_type == "localization":
            return row.get(LEADER_LOCALIZATION_COLUMN, "") == self.prediction
        return _leader_call_prediction_from_columns(row) == self.prediction

    def _columns_match_prediction(self, columns):
        if self.call_type == "localization":
            return columns.get(LEADER_LOCALIZATION_COLUMN, "") == self.prediction
        return _leader_call_prediction_from_columns(columns) == self.prediction

    def _call_matches_prediction(self, call):
        if self.call_type == "localization":
            return call.localization == self.prediction
        return call.prediction == self.prediction

    def scope_sequence_candidates(self, protein, candidates, row):
        return self._scoped_candidates(protein)

    def _scoped_candidates(self, protein):
        anchor = 1
        if self.pfam_accession is not None:
            anchor = _earliest_pfam_anchor(protein, self.pfam_accession)
            if anchor is None:
                return []
            candidates = protein.leader_sequence_candidates_at_anchor(
                anchor,
                _target_prefix(self.pfam_accession),
                self.window_start,
                self.window_end,
            )
            return _dedupe_sequence_candidates([
                candidate for candidate in candidates
                if candidate.start_label or _original_candidate_in_anchor_window(
                    candidate,
                    anchor,
                    self.window_start,
                    self.window_end,
                )
            ])

        candidates = [
            candidate for candidate in protein.sequences_with_leader()
            if candidate.start_label
        ]
        original_candidates = [
            candidate for candidate in protein.sequences_with_leader()
            if not candidate.start_label
        ]
        filtered = [
            candidate for candidate in candidates
            if _coord_in_window(
                _candidate_coord_relative_to_anchor(candidate.start_aa_1b, anchor),
                self.window_start,
                self.window_end,
            )
        ]
        if filtered:
            return _dedupe_sequence_candidates(original_candidates + filtered)
        if original_candidates:
            return original_candidates
        return _original_sequence_candidate(protein)


@dataclass
class TargetPCall:
    prediction: str
    probabilities: dict

    def probability(self, prediction):
        return self.probabilities.get(prediction)


@dataclass
class DeepLocCall:
    prediction: str
    localization: str
    probabilities: dict

    def probability(self, prediction):
        return self.probabilities.get(prediction)


def _format_leader_calls(calls):
    return {
        start: _format_leader_call_columns(call)
        for start, call in calls
    }


def _probability_percent(probability):
    return int(round(probability * 100))


def _format_leader_call_columns(call):
    if isinstance(call, DeepLocCall):
        return _format_deeploc_columns(call)
    return _format_targetp_columns(call)


def _format_targetp_columns(call):
    return {
        TARGETP_PROBABILITY_COLUMNS[prediction]: _format_probability_percent(call.probability(prediction))
        for prediction in TARGETP_PROBABILITY_ORDER
    }


def _format_deeploc_columns(call):
    columns = {
        LEADER_LOCALIZATION_COLUMN: call.localization,
    }
    for prediction in DEEPLOC_SIGNAL_PREDICTIONS:
        columns[TARGETP_PROBABILITY_COLUMNS[prediction]] = _format_probability_percent(call.probability(prediction))
    for prediction in sorted(call.probabilities):
        if prediction in DEEPLOC_SIGNAL_PREDICTIONS:
            continue
        columns[_leader_call_column(prediction)] = _format_probability_percent(call.probability(prediction))
    return columns


def _format_probability_percent(probability):
    if probability is None:
        return ""
    return str(_probability_percent(probability))


def _leader_call_prediction_from_columns(row):
    scored = []
    for prediction in TARGETP_PROBABILITY_ORDER:
        probability = _parse_probability(row.get(TARGETP_PROBABILITY_COLUMNS[prediction], ""))
        if probability is not None:
            scored.append((probability, prediction))
    if not scored:
        return ""
    if max(probability for probability, _prediction in scored) <= 0:
        return ""
    return max(scored)[1]


def _apply_sequence_annotations(row, candidate):
    calls = row.pop("_Leader.calls_by_start", None)
    if calls is None:
        return
    row.update(calls.get(candidate.start_label, _empty_leader_call_columns()))


def _empty_leader_call_columns():
    return {
        column: ""
        for column in TARGETP_PROBABILITY_COLUMNS.values()
    }


def _targetp_call(prediction, no_tp=None, sp=None, mtp=None):
    return TargetPCall(
        prediction=prediction,
        probabilities={
            key: value
            for key, value in [
                ("noTP", no_tp),
                ("SP", sp),
                ("mTP", mtp),
            ]
            if value is not None
        },
    )


def _parse_probability(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _parse_deeploc_csv(path):
    predictions = {}
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        for column in (DEEPLOC_PROTEIN_ID_COLUMN, DEEPLOC_LOCALIZATIONS_COLUMN, DEEPLOC_SIGNALS_COLUMN):
            if column not in fieldnames:
                raise ValueError(f"DeepLoc CSV is missing required column: {column}")
        score_columns = [
            column for column in fieldnames
            if column not in {
                DEEPLOC_PROTEIN_ID_COLUMN,
                DEEPLOC_LOCALIZATIONS_COLUMN,
                DEEPLOC_SIGNALS_COLUMN,
            }
        ]
        for row in reader:
            sequence_id = row.get(DEEPLOC_PROTEIN_ID_COLUMN, "")
            if not sequence_id:
                continue
            if sequence_id in predictions:
                raise ValueError(f"Duplicate DeepLoc Protein_ID: {sequence_id}")
            signal = row.get(DEEPLOC_SIGNALS_COLUMN, "")
            prediction = ""
            probabilities = {
                "mTP": 0.0,
                "SP": 0.0,
            }
            if signal == DEEPLOC_MTP_SIGNAL:
                prediction = "mTP"
                probabilities["mTP"] = 1.0
            elif signal == DEEPLOC_SP_SIGNAL:
                prediction = "SP"
                probabilities["SP"] = 1.0
            for column in score_columns:
                probability = _parse_probability(row.get(column, ""))
                if probability is not None:
                    probabilities[column] = probability
            predictions[sequence_id] = DeepLocCall(
                prediction=prediction,
                localization=row.get(DEEPLOC_LOCALIZATIONS_COLUMN, ""),
                probabilities=probabilities,
            )
    return predictions


def _deeploc_annotation_columns(calls):
    columns = [
        TARGETP_PROBABILITY_COLUMNS[prediction]
        for prediction in DEEPLOC_SIGNAL_PREDICTIONS
    ]
    columns.append(LEADER_LOCALIZATION_COLUMN)
    seen = set(columns)
    for call in calls:
        for prediction in sorted(call.probabilities):
            if prediction in DEEPLOC_SIGNAL_PREDICTIONS:
                continue
            column = _leader_call_column(prediction)
            if column not in seen:
                columns.append(column)
                seen.add(column)
    return columns


def _leader_call_column(prediction):
    return f"Leader.call('{prediction}')"


def _candidate_coord_relative_to_anchor(candidate_start_aa_1b, anchor_aa_1b):
    if candidate_start_aa_1b >= anchor_aa_1b:
        return candidate_start_aa_1b - anchor_aa_1b + 1
    return candidate_start_aa_1b - anchor_aa_1b


def _coord_in_window(coord, start, end):
    low = min(start, end)
    high = max(start, end)
    return coord != 0 and low <= coord <= high


def _original_candidate_in_anchor_window(candidate, anchor_aa_1b, window_start, window_end):
    return _coord_in_window(
        _candidate_coord_relative_to_anchor(candidate.start_aa_1b, anchor_aa_1b),
        window_start,
        window_end,
    )


def _original_sequence_candidate(protein):
    for candidate in protein.sequences_with_leader():
        if not candidate.start_label:
            return [candidate]
    return [
        LeaderSequenceCandidate(
            accession=protein.protein_accession,
            start_label="",
            start_aa_1b=1,
            sequence=protein.sequence(),
        )
    ]


def _dedupe_sequence_candidates(candidates):
    by_sequence = {}
    for candidate in candidates:
        current = by_sequence.get(candidate.sequence)
        if current is None or len(candidate.accession) < len(current.accession):
            by_sequence[candidate.sequence] = candidate
    return [
        candidate for candidate in candidates
        if by_sequence[candidate.sequence].accession == candidate.accession
    ]


def _earliest_pfam_anchor(protein, accession):
    anchors = []
    for row in protein.detected_pfam():
        target = row["target_accession"]
        if not _motif_matches(target, accession):
            continue
        query_start = min(row["query_start"], row["query_end"])
        target_start = row.get("target_start")
        if target_start is None:
            target_start = 1
        anchors.append(query_start - (target_start - 1))
    if not anchors:
        return None
    return min(anchors)


def _parse_leader_calls(calls_by_start):
    calls = []
    for start, columns in calls_by_start.items():
        calls.append((start, columns))
    return calls


def _parse_targetp_output(text):
    predictions = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            predictions[parts[0]] = _targetp_call(
                parts[1],
                _parse_probability(parts[2]) if len(parts) > 2 else None,
                _parse_probability(parts[3]) if len(parts) > 3 else None,
                _parse_probability(parts[4]) if len(parts) > 4 else None,
            )
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

    def evaluate_many(self, contexts, artifacts_dir=None, **kwargs):
        sequence_ids = {_locus_sequence_id(context): context for context in contexts}
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
        if not motif_a_hits and not motif_b_hits:
            return f"missing_{_safe_result_value(self.motif_a)}_and_{_safe_result_value(self.motif_b)}"
        if not motif_a_hits:
            return f"missing_{_safe_result_value(self.motif_a)}"
        if not motif_b_hits:
            return f"missing_{_safe_result_value(self.motif_b)}"
        for a_start, a_end in motif_a_hits:
            for b_start, b_end in motif_b_hits:
                if _edge_distance(a_start, a_end, b_start, b_end) <= self.distance:
                    return RULE_YES
        return RULE_TOO_FAR


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
