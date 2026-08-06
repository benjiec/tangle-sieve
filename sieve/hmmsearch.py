import subprocess
from dataclasses import dataclass


@dataclass
class DomtbloutHit:
    sequence_accession: str
    model_accession: str
    model_name: str
    full_evalue: float
    full_bitscore: float
    domain_evalue: float
    domain_bitscore: float
    hmm_start: int
    hmm_end: int
    sequence_start: int
    sequence_end: int

    def score(self, score_type):
        if score_type == "full":
            return self.full_bitscore
        if score_type == "domain":
            return self.domain_bitscore
        raise ValueError(f"Unsupported score_type: {score_type}")


@dataclass
class KOThreshold:
    threshold: float
    score_type: str


def run_hmmsearch(hmm_file, fasta_file, domtblout_path, use_cut_ga=True):
    cmd = ["hmmsearch"]
    if use_cut_ga:
        cmd.append("--cut_ga")
    cmd.extend(["--domtblout", domtblout_path, hmm_file, fasta_file])
    return subprocess.run(cmd, check=True, capture_output=True, text=True)


def parse_domtblout(path):
    hits = []
    with open(path, "r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(maxsplit=22)
            if len(parts) < 22:
                raise ValueError(f"Expected at least 22 domtblout columns, got: {raw_line.rstrip()}")
            model_accession = parts[4]
            if model_accession == "-":
                model_accession = parts[3]
            hits.append(DomtbloutHit(
                sequence_accession=parts[0],
                model_accession=model_accession,
                model_name=parts[3],
                full_evalue=float(parts[6]),
                full_bitscore=float(parts[7]),
                domain_evalue=float(parts[12]),
                domain_bitscore=float(parts[13]),
                hmm_start=int(parts[15]),
                hmm_end=int(parts[16]),
                sequence_start=int(parts[17]),
                sequence_end=int(parts[18]),
            ))
    return hits


def read_ko_thresholds(path):
    thresholds = {}
    with open(path, "r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(maxsplit=3)
            if parts[0] == "model":
                continue
            if len(parts) < 3:
                raise ValueError(f"Expected model, threshold, and score_type columns: {raw_line.rstrip()}")
            score_type = parts[2]
            if score_type not in {"domain", "full"}:
                raise ValueError(f"Unsupported KO threshold score_type for {parts[0]}: {score_type}")
            thresholds[parts[0]] = KOThreshold(float(parts[1]), score_type)
    return thresholds


def detected_rows_from_hits(hits, database, query_database="", threshold_by_model=None):
    rows = []
    for hit in hits:
        score_type = "domain"
        if threshold_by_model is not None:
            if hit.model_accession not in threshold_by_model:
                raise ValueError(f"Missing KO threshold for {hit.model_accession}")
            threshold = threshold_by_model[hit.model_accession]
            score_type = threshold.score_type
            if hit.score(score_type) < threshold.threshold:
                continue
        rows.append({
            "detection_type": "model",
            "detection_method": "hmmsearch",
            "batch": "",
            "query_accession": hit.sequence_accession,
            "query_database": query_database,
            "query_type": "protein",
            "target_accession": hit.model_accession,
            "target_database": database,
            "target_type": "protein",
            "query_start": hit.sequence_start,
            "query_end": hit.sequence_end,
            "target_start": hit.hmm_start,
            "target_end": hit.hmm_end,
            "evalue": hit.domain_evalue if score_type == "domain" else hit.full_evalue,
            "bitscore": hit.score(score_type),
        })
    return rows
