#!/usr/bin/env python3
"""
01_prepare_reference.py — build the combined mapping reference.

Downloads the Chlamydomonas reinhardtii v5.5 nuclear genome (GCF_000002595.2)
from the NCBI FTP site, fetches the chloroplast (NC_005353.1) and mitochondrial
(NC_001638.1) genomes via NCBI Entrez, concatenates everything into one FASTA,
then builds the minimap2 index and BLASTN database used by the downstream steps.

Requires: biopython (pip); minimap2 and makeblastdb on PATH
          (tested with minimap2 2.31 and BLASTN 2.16.0+).

Usage:
    python 01_prepare_reference.py --outdir reference --email you@example.com
"""

import argparse
import gzip
import shutil
import subprocess
import urllib.request
from pathlib import Path

from Bio import Entrez, SeqIO

GENOME_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/595/"
    "GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5/"
    "GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.fna.gz"
)
ORGANELLE_ACCESSIONS = ["NC_005353.1", "NC_001638.1"]  # chloroplast, mitochondrion


def download_genome(dest: Path) -> Path:
    """Download and decompress the v5.5 nuclear genome FASTA."""
    gz = dest / "Cr_v5.5_genomic.fna.gz"
    fasta = dest / "Cr_v5.5_genomic.fna"
    if not fasta.exists():
        print(f"Downloading {GENOME_URL}")
        urllib.request.urlretrieve(GENOME_URL, gz)
        with gzip.open(gz, "rb") as fi, open(fasta, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        gz.unlink()
    print(f"Nuclear genome: {fasta}")
    return fasta


def fetch_organelles(dest: Path, email: str) -> Path:
    """Fetch chloroplast and mitochondrial genomes via Entrez."""
    out = dest / "organelles.fna"
    if not out.exists():
        Entrez.email = email
        recs = []
        for acc in ORGANELLE_ACCESSIONS:
            with Entrez.efetch(db="nucleotide", id=acc, rettype="fasta",
                               retmode="text") as h:
                rec = SeqIO.read(h, "fasta")
            recs.append(rec)
            print(f"Fetched {acc}: {len(rec.seq):,} bp")
        SeqIO.write(recs, out, "fasta")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="reference", help="output directory")
    ap.add_argument("--email", required=True, help="email for NCBI Entrez")
    args = ap.parse_args()

    dest = Path(args.outdir)
    dest.mkdir(parents=True, exist_ok=True)

    nuclear = download_genome(dest)
    organelles = fetch_organelles(dest, args.email)

    combined = dest / "combined_reference.fna"
    with open(combined, "wb") as fo:
        for f in (nuclear, organelles):
            with open(f, "rb") as fi:
                shutil.copyfileobj(fi, fo)
    n_seqs = sum(1 for line in open(combined) if line.startswith(">"))
    print(f"Combined reference: {combined} ({n_seqs} sequences)")

    print("Building minimap2 index ...")
    subprocess.run(["minimap2", "-d", str(dest / "combined_reference.mmi"),
                    str(combined)], check=True)
    print("Building BLASTN database ...")
    subprocess.run(["makeblastdb", "-in", str(combined), "-dbtype", "nucl",
                    "-out", str(dest / "combined_reference_db")], check=True)
    print("Done.")


if __name__ == "__main__":
    main()
