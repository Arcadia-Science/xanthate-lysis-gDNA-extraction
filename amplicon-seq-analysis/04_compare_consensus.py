#!/usr/bin/env python3
"""
04_compare_consensus.py — verify vendor consensus sequences against the reference.

For each vendor consensus FASTA (one per successfully assembled sample):
  1. Locus confirmation: BLASTN against the combined reference database
     (from 01_prepare_reference.py); reports the top hit. This catches
     off-target amplification (e.g. a reaction whose product maps to the
     wrong chromosome).
  2. Within-reaction identity: pairwise alignment of each consensus to the
     reaction's anchor sample (the CRC_cc124 control when present) with MAFFT.
     Long-read amplicon consensuses come back in arbitrary orientation, so
     each comparison is scored against both the sequence and its reverse
     complement and the better-scoring orientation is kept.

Expects FASTAs named <order>_<num>_<reaction>_<prep>_<sample>.fasta
(e.g. AI-8632_1_A_HF-7_cc124_A3.fasta). Samples that failed de novo assembly
have no consensus file and are simply absent.

Requires: biopython, pandas (pip); mafft and blastn on PATH
          (tested with MAFFT v7.526, BLASTN 2.16.0+).

Usage:
    python 04_compare_consensus.py --fasta-dir fastas \
        --blast-db reference/combined_reference_db --out results
"""

import argparse
import glob
import os
import subprocess
import tempfile

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# Expected full-length amplicon sizes (bp) per reaction and genotype.
# Vendor consensuses are primer-trimmed, so observed lengths run ~20-45 bp shorter.
# Validated panel: A (chr9), D (rbcL), E (cox1). C (chr13 candidate) is listed
# because the sequenced run contains those samples; it amplified off-target.
EXPECTED_BP = {
    "A": "1937 (CC-124) / ~149 (CC-1373)",   # chr9 deletion marker
    "C": "2167 (CC-124) / ~204 (CC-1373)",   # chr13 candidate (failed: off-target)
    "D": "389 (both)",                        # rbcL (chloroplast)
    "E": "455 (both)",                        # cox1 (mitochondrion)
}


def parse_sample(fasta_path):
    """<order>_<num>_<reaction>_<prep>_<sample>.fasta -> (num, reaction, sample)."""
    base = os.path.basename(fasta_path)
    base = base.rsplit(".", 1)[0]
    parts = base.split("_")
    return int(parts[1]), parts[2], "_".join(parts[4:])


def load_consensus(fasta_dir):
    records = {}
    for f in sorted(glob.glob(os.path.join(fasta_dir, "*.fasta"))):
        num, rxn, sample = parse_sample(f)
        records[(rxn, sample)] = dict(num=num, file=os.path.basename(f),
                                      seq=str(SeqIO.read(f, "fasta").seq))
    return records


def blastn_top_hit(seq, db):
    """Top BLASTN hit (by bitscore) of a consensus against the combined reference."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tf:
        tf.write(f">q\n{seq}\n")
        path = tf.name
    out = subprocess.run(
        ["blastn", "-query", path, "-db", db, "-max_target_seqs", "3",
         "-outfmt", "6 sseqid pident length bitscore"],
        capture_output=True, text=True, check=True).stdout
    os.unlink(path)
    hits = [l.split("\t") for l in out.splitlines() if l.strip()]
    if not hits:
        return dict(top_hit="no_hit", pident=None, aln_len=0)
    top = max(hits, key=lambda h: float(h[3]))
    return dict(top_hit=top[0], pident=float(top[1]), aln_len=int(top[2]))


def pairwise_identity(seq_a, seq_b):
    """Identity (%) with strand auto-detection: try seq_b and revcomp(seq_b), keep best.

    Identity = matched columns / columns where both sequences have a base,
    from a MAFFT --auto pairwise alignment.
    """
    best = None
    for strand, s in (("fwd", seq_b), ("revcomp", str(Seq(seq_b).reverse_complement()))):
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tf:
            tf.write(f">A\n{seq_a}\n>B\n{s}\n")
            path = tf.name
        out = subprocess.run(["mafft", "--quiet", "--auto", path],
                             capture_output=True, text=True, check=True).stdout
        os.unlink(path)
        seqs, cur = {}, None
        for line in out.splitlines():
            if line.startswith(">"):
                cur = line[1:].strip()
                seqs[cur] = ""
            else:
                seqs[cur] += line.strip()
        a, b = seqs["A"].upper(), seqs["B"].upper()
        match = sum(1 for x, y in zip(a, b) if x == y and x != "-")
        both = sum(1 for x, y in zip(a, b) if x != "-" and y != "-")
        pct = 100 * match / both if both else 0.0
        if best is None or pct > best["pct_identity"]:
            best = dict(strand=strand, pct_identity=round(pct, 2),
                        mismatches=both - match)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta-dir", required=True, help="directory of vendor consensus FASTAs")
    ap.add_argument("--blast-db", required=True, help="BLASTN db prefix from 01_prepare_reference.py")
    ap.add_argument("--out", default="results", help="output directory for CSV table")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    records = load_consensus(args.fasta_dir)
    rows = []
    for rxn in sorted({k[0] for k in records}):
        sub = {k: v for k, v in records.items() if k[0] == rxn}
        anchor_name = "CRC_cc124" if (rxn, "CRC_cc124") in sub else sorted(sub)[0][1]
        anchor = sub[(rxn, anchor_name) if (rxn, anchor_name) in sub
                     else (rxn, sorted(sub)[0][1])]
        for (r, sample), rec in sorted(sub.items(), key=lambda kv: kv[1]["num"]):
            hit = blastn_top_hit(rec["seq"], args.blast_db)
            if sample == anchor_name:
                ident = dict(strand="anchor", pct_identity=100.0, mismatches=0)
            else:
                ident = pairwise_identity(anchor["seq"], rec["seq"])
            rows.append(dict(num=rec["num"], rxn=r, sample=sample,
                             length=len(rec["seq"]),
                             expected_bp=EXPECTED_BP.get(r, ""),
                             top_hit=hit["top_hit"], hit_pident=hit["pident"],
                             identity_vs_anchor=ident["pct_identity"],
                             mismatches=ident["mismatches"], strand=ident["strand"]))
            pid = f"{hit['pident']}%" if hit["pident"] is not None else "no hit"
            print(f"{r} {sample:14s} {len(rec['seq']):5d} bp  hit={hit['top_hit']:14s} "
                  f"({pid} id)  vs anchor: {ident['pct_identity']}% "
                  f"({ident['mismatches']} mm, {ident['strand']})")

    df = pd.DataFrame(rows).sort_values(["rxn", "num"]).reset_index(drop=True)
    out_csv = os.path.join(args.out, "consensus_comparison.csv")
    df.to_csv(out_csv, index=False)
    print(f"\nTable written to {out_csv}")


if __name__ == "__main__":
    main()
