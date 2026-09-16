"""Shared content hashing and atomic TSV cache I/O."""

import csv
import hashlib
import os
import tempfile


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_cache(path, required_columns=("zipfile", "zipfile_sha256", "model")):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None, []
    with open(path, encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        columns = reader.fieldnames
        if columns is None or any(column not in columns for column in required_columns):
            raise ValueError(f"cache TSV has an invalid header: {path}")
        rows = list(reader)
    if any(None in row for row in rows):
        raise ValueError(f"cache TSV contains rows wider than its header: {path}")
    return columns, rows


def write_cache(path, columns, rows):
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            "w", encoding="utf-8", newline="", dir=directory, delete=False
        ) as stream:
            temporary_path = stream.name
            writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_path, path)
        temporary_path = None
    finally:
        if temporary_path is not None and os.path.exists(temporary_path):
            os.remove(temporary_path)


