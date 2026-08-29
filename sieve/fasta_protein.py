import os
import subprocess
import tempfile
from io import StringIO

from Bio import AlignIO

from sieve.protein import (
    LEADER_CANDIDATE_MAX_PROTEIN_START_1B,
    LeaderSequenceCandidate,
    ProteinHMMAlignment,
    accession_with_suffix,
    _leader_relative_label,
    _leader_start_label,
    _protein_start_at_anchor,
)


class FastaProtein(object):

    def __init__(self, protein_accession, sequence, pfam_rows=None, ko_rows=None):
        self.protein_accession = protein_accession
        self.genome_accession = ""
        self._sequence = sequence
        self._pfam_rows = list(pfam_rows or [])
        self._ko_rows = list(ko_rows or [])
        self._leader_sequence_candidates_cache = None
        self._leader_sequence_candidates_at_anchor_cache = {}

    def sequence(self):
        return self._sequence

    def _candidate_accession(self, suffix):
        return accession_with_suffix(self.protein_accession, suffix)

    def detected_pfam(self):
        return self._pfam_rows

    def detected_ko(self):
        return self._ko_rows

    def _original_sequence_candidate(self):
        return [
            LeaderSequenceCandidate(
                accession=self.protein_accession,
                start_label="",
                start_aa_1b=1,
                sequence=self.sequence(),
                protein_start_aa_1b=1,
            )
        ]

    def leader_sequence_candidates(self):
        if self._leader_sequence_candidates_cache is not None:
            return self._leader_sequence_candidates_cache

        sequence = self.sequence()
        last_candidate_index = min(len(sequence), LEADER_CANDIDATE_MAX_PROTEIN_START_1B)
        candidates = self._original_sequence_candidate()
        for index, aa in enumerate(sequence[:last_candidate_index]):
            if aa != "M":
                continue
            start_aa_1b = index + 1
            start_label = _leader_start_label(start_aa_1b)
            candidates.append(LeaderSequenceCandidate(
                accession=self._candidate_accession(f"with_leader_{start_label}_M"),
                start_label=start_label,
                start_aa_1b=start_aa_1b,
                sequence=sequence[index:],
                protein_start_aa_1b=start_aa_1b,
            ))
        self._leader_sequence_candidates_cache = candidates
        return self._leader_sequence_candidates_cache

    def sequences_with_leader(self):
        return self.leader_sequence_candidates()

    def leader_sequence_candidates_at_anchor(self, anchor_aa_1b, anchor_label, window_start, window_end):
        cache_key = (anchor_aa_1b, anchor_label, window_start, window_end)
        if cache_key in self._leader_sequence_candidates_at_anchor_cache:
            return self._leader_sequence_candidates_at_anchor_cache[cache_key]

        low = min(window_start, window_end)
        high = max(window_start, window_end)
        candidates = self._original_sequence_candidate()
        seen_start_labels = set()
        sequence = self.sequence()
        anchor_index = anchor_aa_1b - 1
        if anchor_index < 0 or anchor_index >= len(sequence):
            self._leader_sequence_candidates_at_anchor_cache[cache_key] = candidates
            return self._leader_sequence_candidates_at_anchor_cache[cache_key]

        for index, aa in enumerate(sequence):
            if aa != "M":
                continue
            if index < anchor_index:
                relative_start = index - anchor_index
            else:
                relative_start = index - anchor_index + 1
            if relative_start < low or relative_start > high:
                continue
            start_label = f"{_leader_relative_label(relative_start)}_{anchor_label}"
            if start_label in seen_start_labels:
                continue
            seen_start_labels.add(start_label)
            candidates.append(LeaderSequenceCandidate(
                accession=self._candidate_accession(f"with_leader_{start_label}_M"),
                start_label=start_label,
                start_aa_1b=relative_start,
                sequence=sequence[index:],
                protein_start_aa_1b=_protein_start_at_anchor(anchor_aa_1b, relative_start),
            ))
        self._leader_sequence_candidates_at_anchor_cache[cache_key] = candidates
        return self._leader_sequence_candidates_at_anchor_cache[cache_key]

    def genomic_locus(self):
        raise ValueError(f"Genomic locus is not available for FASTA-only protein {self.protein_accession}")

    def genomic_locus_with_leader(self):
        raise ValueError(f"Genomic locus is not available for FASTA-only protein {self.protein_accession}")

    def hmm_align(self, hmm_profile_fn):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".faa", mode="wt") as fasta:
            fasta.write(f">{self.protein_accession}\n{self.sequence()}\n")
            fasta_path = fasta.name
        try:
            cmd = ["hmmalign", hmm_profile_fn, fasta_path]
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            alignment = AlignIO.read(StringIO(result.stdout), "stockholm")
            return ProteinHMMAlignment(alignment, self.protein_accession)
        finally:
            os.remove(fasta_path)
