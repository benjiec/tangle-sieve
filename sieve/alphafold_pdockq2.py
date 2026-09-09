import csv
import io
import json
import math
import os
import re
import zipfile
from dataclasses import dataclass
from itertools import combinations

from Bio.PDB import MMCIFParser


PDOCKQ2_PARAMETERS = (1.31034849, 84.7326239, 0.0747157696, 0.00501886443)
OUTPUT_COLUMNS = [
    "source", "model", "scope", "chain 1", "chain 2", "region 1", "region 2",
    "contact count", "interface residues 1", "interface residues 2",
    "interface pLDDT 1", "interface pLDDT 2", "normalized PAE 1 to 2",
    "normalized PAE 2 to 1", "pDockQ2 1 to 2", "pDockQ2 2 to 1",
    "pDockQ2 max",
]


@dataclass(frozen=True)
class Residue:
    chain: str
    number: int
    token_index: int
    plddt: float
    ca: tuple[float, float, float]


def parse_region(value):
    match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", value)
    if not match:
        raise ValueError(f"invalid region {value!r}; expected CHAIN:START-END")
    chain, start, end = match.group(1), int(match.group(2)), int(match.group(3))
    if start < 1 or end < start:
        raise ValueError(f"invalid region {value!r}; require 1 <= START <= END")
    return chain, start, end


def collect_regions(values):
    regions = {}
    for value in values:
        chain, start, end = parse_region(value)
        if chain in regions:
            raise ValueError(f"more than one region specified for chain {chain}")
        regions[chain] = [(start, end)]
    return regions


def region_label(ranges):
    return "all" if ranges is None else ",".join(f"{start}-{end}" for start, end in ranges)


def _selected(number, ranges):
    return ranges is None or any(start <= number <= end for start, end in ranges)


def _distance(first, second):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(first, second)))


def _mean(values):
    return sum(values) / len(values) if values else 0.0


def _pdockq2(normalized_pae, interface_plddt):
    if normalized_pae == 0.0 or interface_plddt == 0.0:
        return 0.0
    length, midpoint, slope, intercept = PDOCKQ2_PARAMETERS
    x = normalized_pae * interface_plddt
    return length / (1.0 + math.exp(-slope * (x - midpoint))) + intercept


def parse_model(cif_text, full_data):
    chain_ids = full_data.get("token_chain_ids")
    residue_ids = full_data.get("token_res_ids")
    pae = full_data.get("pae")
    if not isinstance(chain_ids, list) or not isinstance(residue_ids, list):
        raise ValueError("full-data JSON lacks token_chain_ids or token_res_ids")
    if len(chain_ids) != len(residue_ids):
        raise ValueError("token_chain_ids and token_res_ids have different lengths")
    token_by_residue = {}
    for index, key in enumerate(zip(chain_ids, residue_ids)):
        key = (str(key[0]), int(key[1]))
        if key in token_by_residue:
            raise ValueError(f"duplicate token mapping for chain {key[0]} residue {key[1]}")
        token_by_residue[key] = index
    if not isinstance(pae, list) or len(pae) != len(chain_ids) or any(
        not isinstance(row, list) or len(row) != len(chain_ids) for row in pae
    ):
        raise ValueError("PAE matrix dimensions do not match token metadata")

    structure = MMCIFParser(QUIET=True).get_structure("alphafold", io.StringIO(cif_text))
    residues = {}
    observed_tokens = set()
    for chain in structure[0]:
        chain_residues = []
        for residue in chain:
            if "CA" not in residue:
                continue
            number = int(residue.id[1])
            key = (chain.id, number)
            if key not in token_by_residue:
                raise ValueError(f"CIF residue {chain.id}:{number} has no token mapping")
            atom = residue["CA"]
            observed_tokens.add(token_by_residue[key])
            chain_residues.append(Residue(
                chain.id,
                number,
                token_by_residue[key],
                float(atom.bfactor),
                tuple(float(value) for value in atom.coord),
            ))
        residues[chain.id] = chain_residues
    missing = sorted(set(chain_ids) - set(residues))
    if missing:
        raise ValueError(f"CIF lacks chains present in token metadata: {', '.join(missing)}")
    missing_tokens = sorted(set(range(len(chain_ids))) - observed_tokens)
    if missing_tokens:
        index = missing_tokens[0]
        raise ValueError(
            f"token mapping {chain_ids[index]}:{residue_ids[index]} has no CIF C-alpha atom"
        )
    return residues, pae


def validate_regions(regions, residues):
    for chain, ranges in regions.items():
        if chain not in residues:
            raise ValueError(f"region refers to unknown chain {chain}")
        known = {residue.number for residue in residues[chain]}
        for start, end in ranges:
            absent = next((number for number in range(start, end + 1) if number not in known), None)
            if absent is not None:
                raise ValueError(f"region {chain}:{start}-{end} includes absent residue {absent}")


def score_pair(chain1, chain2, residues, pae, regions, cutoff=8.0):
    selected1 = [r for r in residues[chain1] if _selected(r.number, regions.get(chain1))]
    selected2 = [r for r in residues[chain2] if _selected(r.number, regions.get(chain2))]
    contacts = [(r1, r2) for r1 in selected1 for r2 in selected2 if _distance(r1.ca, r2.ca) <= cutoff]
    interface1 = {r1.number: r1 for r1, _ in contacts}
    interface2 = {r2.number: r2 for _, r2 in contacts}
    plddt1 = _mean([residue.plddt for residue in interface1.values()])
    plddt2 = _mean([residue.plddt for residue in interface2.values()])
    norm12 = _mean([1.0 / (1.0 + (float(pae[r1.token_index][r2.token_index]) / 10.0) ** 2) for r1, r2 in contacts])
    norm21 = _mean([1.0 / (1.0 + (float(pae[r2.token_index][r1.token_index]) / 10.0) ** 2) for r1, r2 in contacts])
    score12 = _pdockq2(norm12, plddt1)
    score21 = _pdockq2(norm21, plddt2)
    return {
        "chain 1": chain1,
        "chain 2": chain2,
        "region 1": region_label(regions.get(chain1)),
        "region 2": region_label(regions.get(chain2)),
        "contact count": len(contacts),
        "interface residues 1": ",".join(str(number) for number in sorted(interface1)),
        "interface residues 2": ",".join(str(number) for number in sorted(interface2)),
        "interface pLDDT 1": plddt1,
        "interface pLDDT 2": plddt2,
        "normalized PAE 1 to 2": norm12,
        "normalized PAE 2 to 1": norm21,
        "pDockQ2 1 to 2": score12,
        "pDockQ2 2 to 1": score21,
        "pDockQ2 max": max(score12, score21),
    }


def score_model(cif_text, full_data, regions=None, cutoff=8.0):
    if cutoff <= 0:
        raise ValueError("distance cutoff must be greater than zero")
    regions = regions or {}
    residues, pae = parse_model(cif_text, full_data)
    validate_regions(regions, residues)
    return score_pairs(residues, pae, regions, cutoff)


def score_pairs(residues, pae, regions=None, cutoff=8.0):
    regions = regions or {}
    rows = []
    for chain1, chain2 in combinations(residues, 2):
        full_row = score_pair(chain1, chain2, residues, pae, {}, cutoff)
        full_row["scope"] = "full"
        rows.append(full_row)
        if chain1 in regions or chain2 in regions:
            regional_row = score_pair(chain1, chain2, residues, pae, regions, cutoff)
            regional_row["scope"] = "regional"
            rows.append(regional_row)
    return rows


def model_number(name):
    match = re.search(r"_model_(\d+)\.cif$", name)
    return int(match.group(1)) if match else None


def score_zip(path, regions=None, cutoff=8.0):
    rows = []
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        models = {model_number(name): name for name in names if model_number(name) is not None}
        if not models:
            raise ValueError(f"archive contains no *_model_N.cif files: {path}")
        for number, cif_name in sorted(models.items()):
            prefix = cif_name.removesuffix(f"model_{number}.cif")
            data_name = f"{prefix}full_data_{number}.json"
            if data_name not in names:
                raise ValueError(f"model {number} lacks matching {data_name}")
            cif_text = archive.read(cif_name).decode("utf-8")
            full_data = json.loads(archive.read(data_name))
            for row in score_model(cif_text, full_data, regions, cutoff):
                row["source"] = os.path.basename(path)
                row["model"] = number
                rows.append(row)
    return rows


def read_model_files(cif_path, full_data_path, regions=None, cutoff=8.0):
    with open(cif_path, encoding="utf-8") as stream:
        cif_text = stream.read()
    with open(full_data_path, encoding="utf-8") as stream:
        full_data = json.load(stream)
    rows = score_model(cif_text, full_data, regions, cutoff)
    number = model_number(os.path.basename(cif_path))
    for row in rows:
        row["source"] = os.path.basename(cif_path)
        row["model"] = "" if number is None else number
    return rows


def write_rows(rows, output):
    writer = csv.DictWriter(output, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        formatted = dict(row)
        for column in OUTPUT_COLUMNS:
            if isinstance(formatted.get(column), float):
                formatted[column] = f"{formatted[column]:.6f}"
        writer.writerow(formatted)


def best_row(rows, scope=None):
    candidates = rows if scope is None else [row for row in rows if row["scope"] == scope]
    return max(candidates, key=lambda row: row["pDockQ2 max"]) if candidates else None


def best_rows_by_pair(rows):
    grouped = {}
    for row in rows:
        key = (row["chain 1"], row["chain 2"])
        grouped.setdefault(key, []).append(row)
    return [
        (key, best_row(grouped[key], "full"), best_row(grouped[key], "regional"))
        for key in sorted(grouped)
    ]


def best_full_scores_by_pair(rows):
    return {
        key: full["pDockQ2 max"]
        for key, full, _regional in best_rows_by_pair(rows)
        if full is not None
    }


def full_scores_by_model_and_pair(rows):
    scores = {}
    for row in rows:
        if row["scope"] != "full":
            continue
        key = (row["model"], row["chain 1"], row["chain 2"])
        if key in scores:
            raise ValueError(
                f"duplicate full score for model {row['model']} pair "
                f"{row['chain 1']}-{row['chain 2']}"
            )
        scores[key] = row["pDockQ2 max"]
    return scores


def score_column_name(row):
    if row["scope"] == "full":
        return f"{row['chain 1']}_{row['chain 2']}_pDockQ2_max"
    return (
        f"{row['chain 1']}_{row['region 1']}_"
        f"{row['chain 2']}_{row['region 2']}_pDockQ2_max"
    )


def scores_by_model_and_column(rows):
    scores = {}
    for row in rows:
        key = (row["model"], score_column_name(row))
        if key in scores:
            raise ValueError(f"duplicate score for model {key[0]} column {key[1]}")
        scores[key] = row["pDockQ2 max"]
    return scores
