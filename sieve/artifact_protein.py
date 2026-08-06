from sieve.protein import CuratedProtein


class ArtifactProtein(object):

    def __init__(self, protein_accession, genome_accession, sequence, pfam_rows=None, ko_rows=None):
        self.protein_accession = protein_accession
        self.genome_accession = genome_accession
        self._sequence = sequence
        self._pfam_rows = list(pfam_rows or [])
        self._ko_rows = list(ko_rows or [])
        self._curated = CuratedProtein(protein_accession, genome_accession) if genome_accession else None

    def sequence(self):
        return self._sequence

    def detected_pfam(self):
        return self._pfam_rows

    def detected_ko(self):
        return self._ko_rows

    def sequences_with_leader(self):
        raise RuntimeError("Leader discovery is not available during artifact evaluation")

    def leader_sequence_candidates_at_anchor(self, *args, **kwargs):
        raise RuntimeError("Leader discovery is not available during artifact evaluation")

    def genomic_locus(self):
        if self._curated is None:
            raise ValueError(f"Genomic locus is not available for FASTA-only protein {self.protein_accession}")
        return self._curated.genomic_locus()

    def genomic_locus_with_leader(self):
        if self._curated is None:
            raise ValueError(f"Genomic locus is not available for FASTA-only protein {self.protein_accession}")
        return self._curated.genomic_locus_with_leader()
