import os
import subprocess
import tempfile
from dataclasses import dataclass
from io import StringIO
from typing import Optional

import duckdb
from Bio import AlignIO
from Bio.Seq import Seq

from tangle.defaults import Defaults
from tangle.detected import DetectedTable
from tangle.manifest import ManifestTable
from tangle.models import CSVSource, Schema
from tangle import unique_batch
from tangle.sequence import read_fasta_as_dict, extract_subsequence_strand_sensitive


SEQUENCE_SOURCE_HMM_DETECTED = "hmm-detected"
SEQUENCE_SOURCE_NCBI = "ncbi"

SEQUENCE_TYPE_PROTEIN = "protein"
GFF_TYPE_CDS = "CDS"
GFF_TYPE_GENE = "gene"
GFF_TYPE_MRNA = "mRNA"
GFF_TYPE_TRANSCRIPT = "transcript"
GFF_TYPE_START_CODON = "start_codon"
GFF_TYPE_STOP_CODON = "stop_codon"

_DUCKDB_TABLE_CACHE = {}
_FASTA_FILE_CACHE = {}


def _existing_path(path):
    if path is None:
        raise FileNotFoundError("Path is not available")
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return str(path)


def _sql_string(value):
    return "'" + str(value).replace("'", "''") + "'"


def _file_cache_key(table, path):
    if path is None or not os.path.exists(path):
        return None
    path = str(path)
    stat = os.stat(path)
    return (table.name, path, stat.st_mtime_ns, stat.st_size)


def _cached_table_name(table, path):
    key = _file_cache_key(table, path)
    if key is None:
        return None
    if key not in _DUCKDB_TABLE_CACHE:
        schema = Schema("__curated_protein__" + unique_batch())
        schema.add_table(CSVSource(table, str(path)))
        schema.duckdb_load()
        _DUCKDB_TABLE_CACHE[key] = (schema.name, table.name)
    schema_name, table_name = _DUCKDB_TABLE_CACHE[key]
    return f"{schema_name}.{table_name}"


def _file_cache_key_for_path(path):
    if path is None or not os.path.exists(path):
        return None
    path = str(path)
    stat = os.stat(path)
    return (path, stat.st_mtime_ns, stat.st_size)


def _cached_fasta(path):
    path = _existing_path(path)
    key = _file_cache_key_for_path(path)
    if key not in _FASTA_FILE_CACHE:
        _FASTA_FILE_CACHE[key] = read_fasta_as_dict(path)
    return _FASTA_FILE_CACHE[key]


def _rows_from_table(table, path, column_filters=None):
    table_name = _cached_table_name(table, path)
    if table_name is None:
        return []
    if column_filters:
        cond = f" WHERE {' AND '.join(column_filters)}"
    else:
        cond = ""
    query = f"SELECT * FROM {table_name}{cond}"
    return duckdb.execute(query).fetchdf().to_dict("records")


def _relative_pos(locus_start, locus_end, genomic_pos):
    if locus_start <= locus_end:
        return genomic_pos - locus_start + 1
    return locus_start - genomic_pos + 1


def _interval_relative_to_locus(locus_start, locus_end, feature_start, feature_end):
    rel_start = _relative_pos(locus_start, locus_end, feature_start)
    rel_end = _relative_pos(locus_start, locus_end, feature_end)
    return (min(rel_start, rel_end), max(rel_start, rel_end))


def _parse_gff_attributes(raw):
    attrs = {}
    for part in raw.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
            value = value.strip('"')
        else:
            key, value = part, ""
        attrs[key] = value
    return attrs


def _parse_gff_line(line):
    if not line.strip() or line.startswith("#"):
        return None
    parts = line.rstrip("\n").split("\t")
    if len(parts) != 9:
        return None
    return dict(
        seqid=parts[0],
        type=parts[2],
        start=int(parts[3]),
        end=int(parts[4]),
        strand=parts[6],
        attrs=_parse_gff_attributes(parts[8]),
    )


def _gff_attr_values(attrs, key):
    value = attrs.get(key)
    if value is None:
        return []
    return [v for v in value.split(",") if v]


@dataclass
class GenomicLocus:
    genome_accession: str
    contig_accession: str
    start_1b: int
    end_1b: int
    strand: int
    _sequence: str
    cds_intervals_1b: list
    start_codon_interval_1b: Optional[tuple] = None
    stop_codon_interval_1b: Optional[tuple] = None

    def sequence(self):
        return self._sequence

    def start_codon_position_1b(self):
        if self.start_codon_interval_1b is None:
            return None
        return self.start_codon_interval_1b[0]

    def stop_codon_position_1b(self):
        if self.stop_codon_interval_1b is None:
            return None
        return self.stop_codon_interval_1b[0]

    def dss_positions_1b(self):
        if len(self.cds_intervals_1b) < 2:
            return []
        return [right for _left, right in self.cds_intervals_1b[:-1]]

    def ass_positions_1b(self):
        if len(self.cds_intervals_1b) < 2:
            return []
        return [left for left, _right in self.cds_intervals_1b[1:]]


class CuratedProtein(object):

    @staticmethod
    def clear_cache():
        for schema_name, _table_name in _DUCKDB_TABLE_CACHE.values():
            duckdb.execute(f"DROP SCHEMA IF EXISTS {schema_name} CASCADE")
        _DUCKDB_TABLE_CACHE.clear()
        _FASTA_FILE_CACHE.clear()

    def __init__(self, protein_accession, genome_accession):
        self.protein_accession = protein_accession
        self.genome_accession = genome_accession
        self._manifest_entry = None
        self._sequence = None
        self._detected_rows = None
        self._genomic_fasta = None
        self._locus = None
        self._leader_prefix_cache = None
        self._start_trim_len_cache = None

    @property
    def manifest_entry(self):
        if self._manifest_entry is None:
            matches = _rows_from_table(ManifestTable, Defaults.area_sequence_manifest_tsv(), column_filters=[
                f"sequence_accession = {_sql_string(self.protein_accession)}",
                f"sequence_database = {_sql_string(self.genome_accession)}",
                f"sequence_type = {_sql_string(SEQUENCE_TYPE_PROTEIN)}",
            ])
            if not matches:
                raise ValueError(f"Cannot find protein {self.protein_accession} from {self.genome_accession} in manifest")
            if len(matches) > 1:
                raise ValueError(f"Multiple manifest rows found for protein {self.protein_accession} from {self.genome_accession}")
            self._manifest_entry = matches[0]
        return self._manifest_entry

    @property
    def sequence_source(self):
        return self.manifest_entry["sequence_source"]

    def _protein_fasta_path(self):
        if self.sequence_source == SEQUENCE_SOURCE_HMM_DETECTED:
            return Defaults.area_detected_proteins(self.genome_accession)
        if self.sequence_source == SEQUENCE_SOURCE_NCBI:
            return Defaults.ncbi_genome_proteins(self.genome_accession)
        raise ValueError(f"Unsupported sequence source: {self.sequence_source}")

    def sequence(self):
        if self._sequence is None:
            fasta = _cached_fasta(self._protein_fasta_path())
            if self.protein_accession not in fasta:
                raise ValueError(f"Cannot find protein sequence {self.protein_accession}")
            self._sequence = fasta[self.protein_accession]
        return self._sequence

    def _genomic_sequences(self):
        if self._genomic_fasta is None:
            self._genomic_fasta = _cached_fasta(Defaults.ncbi_genome_fna(self.genome_accession))
        return self._genomic_fasta

    def _detected_protein_rows(self):
        if self._detected_rows is None:
            self._detected_rows = _rows_from_table(
                DetectedTable,
                Defaults.area_detected_proteins_tsv_path(self.genome_accession),
                column_filters=[
                    f"target_accession = {_sql_string(self.protein_accession)}",
                    f"target_database = {_sql_string(self.genome_accession)}",
                ]
            )
        return self._detected_rows

    def _detected_rows_in_protein_order(self):
        rows = self._detected_protein_rows()
        if not rows:
            raise ValueError(f"Cannot find detected fragments for {self.protein_accession}")

        def sort_key(row):
            if row.get("target_start") is None:
                return min(row["query_start"], row["query_end"])
            return min(row["target_start"], row["target_end"])

        return sorted(rows, key=sort_key)

    def _hmm_detected_locus(self, leader_prefix_len=0, start_trim_len=0):
        rows = self._detected_rows_in_protein_order()
        contig_accessions = {row["query_accession"] for row in rows}
        if len(contig_accessions) != 1:
            raise ValueError(f"Detected fragments for {self.protein_accession} span multiple contigs")

        contig_accession = rows[0]["query_accession"]
        first = rows[0]
        last = rows[-1]
        start = first["query_start"]
        end = last["query_end"]
        if leader_prefix_len:
            if start <= end:
                start -= leader_prefix_len * 3
            else:
                start += leader_prefix_len * 3
        if start_trim_len:
            if start <= end:
                start += start_trim_len * 3
            else:
                start -= start_trim_len * 3
        if start < 1 or end < 1:
            raise ValueError(f"Leader extension for {self.protein_accession} extends before the contig")

        contig_seq = self._genomic_sequences().get(contig_accession)
        if contig_seq is None:
            raise ValueError(f"Cannot find contig sequence {contig_accession}")
        sequence = extract_subsequence_strand_sensitive(contig_seq, start, end)

        cds_intervals = []
        locus_left = min(start, end)
        locus_right = max(start, end)
        for row in rows:
            feature_left = min(row["query_start"], row["query_end"])
            feature_right = max(row["query_start"], row["query_end"])
            clipped_left = max(locus_left, feature_left)
            clipped_right = min(locus_right, feature_right)
            if clipped_left <= clipped_right:
                cds_intervals.append(_interval_relative_to_locus(start, end, clipped_left, clipped_right))
        cds_intervals = sorted(cds_intervals)

        start_codon = (1, 3) if self.sequence()[start_trim_len:].startswith("M") or leader_prefix_len else None

        return GenomicLocus(
            genome_accession=self.genome_accession,
            contig_accession=contig_accession,
            start_1b=start,
            end_1b=end,
            strand=1 if start <= end else -1,
            _sequence=sequence,
            cds_intervals_1b=cds_intervals,
            start_codon_interval_1b=start_codon,
        )

    def _iter_gff_rows(self):
        gff_path = _existing_path(Defaults.ncbi_genome_gff(self.genome_accession))
        with open(gff_path, "rt", encoding="utf-8") as f:
            for line in f:
                row = _parse_gff_line(line)
                if row is not None:
                    yield row

    def _gff_gene_blocks(self):
        block = []
        for row in self._iter_gff_rows():
            if row["type"] == GFF_TYPE_GENE and block:
                yield block
                block = []
            block.append(row)
        if block:
            yield block

    def _gff_row_matches_protein(self, row):
        attrs = row["attrs"]
        candidates = []
        for key in ("protein_id", "protein_id_", "Name", "ID", "Derives_from"):
            candidates.extend(_gff_attr_values(attrs, key))
        return self.protein_accession in candidates

    def _gff_rows_for_protein_from_block(self):
        for block in self._gff_gene_blocks():
            gene_ids = {
                row["attrs"].get("ID") for row in block
                if row["type"] == GFF_TYPE_GENE and row["attrs"].get("ID")
            }
            if not gene_ids:
                continue
            transcript_rows = [
                row for row in block
                if row["type"] in (GFF_TYPE_MRNA, GFF_TYPE_TRANSCRIPT)
            ]
            if not transcript_rows:
                continue
            if not any(set(_gff_attr_values(row["attrs"], "Parent")) & gene_ids for row in transcript_rows):
                continue
            if any(row["type"] == GFF_TYPE_CDS and self._gff_row_matches_protein(row) for row in block):
                return block
        return None

    def _gff_rows_for_protein_by_two_pass_scan(self):
        cds_rows = [
            row for row in self._iter_gff_rows()
            if row["type"] == GFF_TYPE_CDS and self._gff_row_matches_protein(row)
        ]
        if not cds_rows:
            return None
        parent_ids = set()
        for row in cds_rows:
            for parent_id in _gff_attr_values(row["attrs"], "Parent"):
                parent_ids.add(parent_id)

        context_rows = []
        for row in self._iter_gff_rows():
            if row in cds_rows:
                context_rows.append(row)
            elif row["type"] in (GFF_TYPE_MRNA, GFF_TYPE_TRANSCRIPT) and row["attrs"].get("ID") in parent_ids:
                context_rows.append(row)
            elif row["type"] in (GFF_TYPE_START_CODON, GFF_TYPE_STOP_CODON) and set(_gff_attr_values(row["attrs"], "Parent")) & parent_ids:
                context_rows.append(row)
        return context_rows

    def _gff_rows_for_protein(self):
        rows = self._gff_rows_for_protein_from_block()
        if rows is not None:
            return rows
        rows = self._gff_rows_for_protein_by_two_pass_scan()
        if rows is not None:
            return rows
        raise ValueError(f"Cannot find CDS rows for protein {self.protein_accession}")

    def _ncbi_locus(self):
        gff_rows = self._gff_rows_for_protein()
        cds_rows = [row for row in gff_rows if row["type"] == GFF_TYPE_CDS and self._gff_row_matches_protein(row)]

        parent_ids = set()
        for row in cds_rows:
            for parent_id in _gff_attr_values(row["attrs"], "Parent"):
                parent_ids.add(parent_id)

        transcript_rows = [
            row for row in gff_rows
            if row["type"] in (GFF_TYPE_MRNA, GFF_TYPE_TRANSCRIPT)
            and row["attrs"].get("ID") in parent_ids
        ]
        if transcript_rows:
            transcript = transcript_rows[0]
            contig_accession = transcript["seqid"]
            start = transcript["start"]
            end = transcript["end"]
            strand = transcript["strand"]
        else:
            contig_accession = cds_rows[0]["seqid"]
            start = min(row["start"] for row in cds_rows)
            end = max(row["end"] for row in cds_rows)
            strand = cds_rows[0]["strand"]

        if strand == "-":
            locus_start, locus_end = end, start
        else:
            locus_start, locus_end = start, end

        contig_seq = self._genomic_sequences().get(contig_accession)
        if contig_seq is None:
            raise ValueError(f"Cannot find contig sequence {contig_accession}")
        sequence = extract_subsequence_strand_sensitive(contig_seq, locus_start, locus_end)

        cds_intervals = [
            _interval_relative_to_locus(locus_start, locus_end, row["start"], row["end"])
            for row in cds_rows
        ]
        cds_intervals = sorted(cds_intervals)

        feature_rows = [
            row for row in gff_rows
            if row["type"] in (GFF_TYPE_START_CODON, GFF_TYPE_STOP_CODON)
            and set(_gff_attr_values(row["attrs"], "Parent")) & parent_ids
        ]
        start_codon = None
        stop_codon = None
        for row in feature_rows:
            interval = _interval_relative_to_locus(locus_start, locus_end, row["start"], row["end"])
            if row["type"] == GFF_TYPE_START_CODON:
                start_codon = interval
            elif row["type"] == GFF_TYPE_STOP_CODON:
                stop_codon = interval

        return GenomicLocus(
            genome_accession=self.genome_accession,
            contig_accession=contig_accession,
            start_1b=locus_start,
            end_1b=locus_end,
            strand=1 if strand != "-" else -1,
            _sequence=sequence,
            cds_intervals_1b=cds_intervals,
            start_codon_interval_1b=start_codon,
            stop_codon_interval_1b=stop_codon,
        )

    def genomic_locus(self):
        if self._locus is None:
            if self.sequence_source == SEQUENCE_SOURCE_HMM_DETECTED:
                self._locus = self._hmm_detected_locus()
            elif self.sequence_source == SEQUENCE_SOURCE_NCBI:
                self._locus = self._ncbi_locus()
            else:
                raise ValueError(f"Unsupported sequence source: {self.sequence_source}")
        return self._locus

    def _leader_prefix(self):
        if self._leader_prefix_cache is not None:
            return self._leader_prefix_cache

        if self.sequence_source != SEQUENCE_SOURCE_HMM_DETECTED or self.sequence().startswith("M") or self._start_trim_len():
            self._leader_prefix_cache = ""
            return self._leader_prefix_cache

        locus = self.genomic_locus()
        contig_seq = self._genomic_sequences()[locus.contig_accession]
        if locus.strand == 1:
            upstream_end = locus.start_1b - 1
            frame_start = upstream_end % 3
            upstream_dna = contig_seq[frame_start:upstream_end]
        else:
            upstream_start = locus.start_1b + 1
            right_len = len(contig_seq) - upstream_start + 1
            usable_len = right_len - (right_len % 3)
            upstream_dna = str(Seq(contig_seq[upstream_start - 1:upstream_start - 1 + usable_len]).reverse_complement())

        usable_len = (len(upstream_dna) // 3) * 3
        upstream_dna = upstream_dna[len(upstream_dna) - usable_len:]
        upstream_aa = str(Seq(upstream_dna).translate(table="Standard", to_stop=False)) if upstream_dna else ""

        tail = upstream_aa.split("*")[-1]
        last_m = tail.rfind("M")
        if last_m >= 0:
            self._leader_prefix_cache = tail[last_m:]
        else:
            self._leader_prefix_cache = tail
        return self._leader_prefix_cache

    def _start_trim_len(self):
        if self._start_trim_len_cache is not None:
            return self._start_trim_len_cache

        self._start_trim_len_cache = 0
        if self.sequence_source == SEQUENCE_SOURCE_HMM_DETECTED:
            sequence = self.sequence()
            if not sequence.startswith("M"):
                for trim_len in (1, 2):
                    if len(sequence) > trim_len and sequence[trim_len] == "M":
                        self._start_trim_len_cache = trim_len
                        break
        return self._start_trim_len_cache

    def sequence_with_leader(self):
        return self._leader_prefix() + self.sequence()[self._start_trim_len():]

    def genomic_locus_with_leader(self):
        if self.sequence_source != SEQUENCE_SOURCE_HMM_DETECTED:
            return self.genomic_locus()
        return self._hmm_detected_locus(
            leader_prefix_len=len(self._leader_prefix()),
            start_trim_len=self._start_trim_len(),
        )

    def _detected_rows_for_query(self, path):
        return _rows_from_table(DetectedTable, path, column_filters=[
            f"query_accession = {_sql_string(self.protein_accession)}",
            f"query_database = {_sql_string(self.genome_accession)}",
        ])

    def detected_pfam(self):
        return self._detected_rows_for_query(Defaults.area_protein_pfam_tsv())

    def detected_ko(self):
        return self._detected_rows_for_query(Defaults.area_protein_ko_assigned_tsv())

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


class ProteinHMMAlignment(object):

    def __init__(self, alignment, protein_accession):
        self.alignment = alignment
        self.protein_accession = protein_accession
        self._by_hmm_pos = self._build_map()

    def _build_map(self):
        records = {record.id: record for record in self.alignment}
        if self.protein_accession not in records:
            raise ValueError(f"Cannot find {self.protein_accession} in HMM alignment")
        record = records[self.protein_accession]
        rf = self.alignment.column_annotations.get("reference_annotation")
        if rf is None:
            rf = self.alignment.column_annotations.get("RF")
        if rf is None:
            rf = "x" * self.alignment.get_alignment_length()

        by_hmm_pos = {}
        hmm_pos = 0
        aa_pos = 0
        for i, rf_char in enumerate(rf):
            aa = record.seq[i]
            if aa not in ".-":
                aa_pos += 1
            if rf_char not in ".-":
                hmm_pos += 1
                if aa not in ".-":
                    by_hmm_pos[hmm_pos] = (aa_pos, str(aa))
                else:
                    by_hmm_pos[hmm_pos] = None
        return by_hmm_pos

    def aa_hmm_pos_1b(self, hmm_pos_1b):
        return self._by_hmm_pos.get(hmm_pos_1b)

    def aa_at_hmm_pos_1b(self, hmm_pos_1b):
        return self.aa_hmm_pos_1b(hmm_pos_1b)

    def hmm_positions_1b(self):
        return self._by_hmm_pos.keys()
