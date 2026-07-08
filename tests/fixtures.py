import os
import tempfile
from pathlib import Path

from tangle.detected import DetectedTable
from tangle.sequence import write_fasta_from_dict


class DefaultsFixture(object):

    def __init__(self, test_case):
        self.test_case = test_case
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.old_world = os.environ.get("TANGLE_WORLD")
        self.old_area = os.environ.get("TANGLE_AREA")
        os.environ["TANGLE_WORLD"] = str(self.root)
        os.environ["TANGLE_AREA"] = "area1"

        self.area_genomics = self.root / "areas" / "area1" / "genomics"
        self.area_metadata = self.root / "areas" / "area1" / "metadata"
        self.ncbi_data = self.root / "ncbi" / "ncbi_dataset" / "data"
        self.area_genomics.mkdir(parents=True)
        self.area_metadata.mkdir(parents=True)
        self.ncbi_data.mkdir(parents=True)

    def cleanup(self):
        if self.old_world is None:
            os.environ.pop("TANGLE_WORLD", None)
        else:
            os.environ["TANGLE_WORLD"] = self.old_world
        if self.old_area is None:
            os.environ.pop("TANGLE_AREA", None)
        else:
            os.environ["TANGLE_AREA"] = self.old_area
        self.tmp.cleanup()

    def genome_dir(self, genome_accession):
        d = self.ncbi_data / genome_accession
        d.mkdir(exist_ok=True)
        return d

    def area_genome_dir(self, genome_accession):
        d = self.area_genomics / genome_accession
        d.mkdir(exist_ok=True)
        return d

    def write_manifest(self, rows):
        from tangle.manifest import ManifestTable

        ManifestTable.write_tsv(str(self.area_genomics / "sequences.tsv"), rows)

    def write_genomic_fasta(self, genome_accession, sequences):
        write_fasta_from_dict(sequences, str(self.genome_dir(genome_accession) / "genomic.fna"))

    def write_ncbi_proteins(self, genome_accession, sequences):
        write_fasta_from_dict(sequences, str(self.genome_dir(genome_accession) / "protein.faa"))

    def write_detected_proteins(self, genome_accession, sequences):
        write_fasta_from_dict(sequences, str(self.area_genome_dir(genome_accession) / "proteins.faa"))

    def write_detected_rows(self, genome_accession, rows):
        DetectedTable.write_tsv(str(self.area_genome_dir(genome_accession) / "proteins.tsv"), rows)

    def write_taxonomy_rows(self, rows):
        import csv

        fieldnames = [
            "Genome Accession",
            "Genome Name",
            "TaxID",
            "Organism",
            "Domain",
            "Kingdom",
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Genus",
            "Species",
        ]
        with open(self.area_metadata / "genomes.tsv", "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def write_gff(self, genome_accession, text):
        with open(self.genome_dir(genome_accession) / "genomic.gff", "w") as f:
            f.write(text)
