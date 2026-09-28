#!/usr/bin/env python3
"""
design/04_organellar_qc_primers.py — organellar presence/absence QC primer design.

Documents and reproduces the organellar marker analysis:

  1. Coverage-gap screen (normalized: >= 20% of median vs <= 10% of median,
     100 bp windows) on the chloroplast (NC_005353.1) and mitochondrial
     (NC_001638.1) alignments. Result for this strain pair: no usable large
     deletions in single-copy regions (candidates fall in the chloroplast
     inverted repeats, which are repetitive and unsuitable for primers).
  2. CIGAR-level indel scan of the CC-1373 alignments. Result: no indels
     >= 30 bp; the largest are 21 bp (chloroplast LSC) and 6 bp
     (mitochondrion) — both below gel resolution. Organellar strain typing by
     gel shift is therefore not feasible; the markers are presence/absence QC.
  3. Variant calling in the rbcL and cox1 gene regions (+/- 500 bp) per strain
     (bcftools mpileup | call, QUAL >= 20, DP >= 5) to place primers in
     SNP-free windows.
  4. Primer3 design (same settings as the nuclear panel) and blastn-short
     uniqueness verification. rbcL is on the minus strand, so its template is
     reverse-complemented before Primer3 and primer positions are mapped back.

Note on reproduction: with these Primer3 settings the cox1 pair (RAF066/067)
is recovered exactly; the rbcL pair depends on the called variant set and
primer3-py version, so an equivalent SNP-free pair may be returned instead of
the historical RAF064/065. Final panel sequences are fixed and verified in
05_make_amplicon_genbank.py.

Requires: primer3-py, biopython, numpy, pandas (pip); samtools, bcftools,
          blastn on PATH.

Usage:
    python design/04_organellar_qc_primers.py \
        --reference reference/combined_reference.fna \
        --ref-bam workdir/bam/CC124.bam --query-bam workdir/bam/CC1373.bam \
        --blast-db reference/combined_reference_db --out organellar_primers.tsv
"""

import argparse
import subprocess
import os
import tempfile

import numpy as np
import pandas as pd
import primer3
from Bio import SeqIO
from Bio.Seq import Seq

# Gene regions (1-based) used for variant calling and primer placement.
LOCI = {
    "rbcL": dict(acc="NC_005353.1", start=122489, end=123917, strand=-1),
    "cox1": dict(acc="NC_001638.1", start=5041, end=6559, strand=+1),
}
FLANK = 500
MIN_INDEL_BP = 30  # CIGAR scan threshold; gel-detectable would need ~50 bp

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module
_p3 = import_module("03_design_typing_primers")
PRIMER3_GLOBAL_ARGS = _p3.PRIMER3_GLOBAL_ARGS
call_flank_variants = _p3.call_flank_variants
blastn_short_uniqueness = _p3.blastn_short_uniqueness


def depth_array(bam, acc):
    out = subprocess.run(["samtools", "depth", "-a", "-r", acc, bam],
                         capture_output=True, text=True, check=True).stdout
    depths = {}
    for line in out.splitlines():
        _, pos, d = line.split("\t")
        depths[int(pos)] = int(d)
    n = max(depths) if depths else 0
    return np.array([depths.get(p, 0) for p in range(1, n + 1)])


def coverage_gap_normalized(ref_depths, query_depths, window=100, min_gap_bp=100):
    """Windows where ref >= 20% of its median and query <= 10% of its median."""
    ref_t = np.median(ref_depths[ref_depths > 0]) * 0.20
    query_t = np.median(query_depths[query_depths > 0]) * 0.10
    gaps, in_gap, gap_start = [], False, 0
    for i in range(0, len(ref_depths) - window, window):
        is_gap = (ref_depths[i:i + window].mean() >= ref_t and
                  query_depths[i:i + window].mean() <= query_t)
        if is_gap and not in_gap:
            gap_start, in_gap = i, True
        elif not is_gap and in_gap:
            in_gap = False
            if i - gap_start >= min_gap_bp:
                gaps.append((gap_start, i))
    return gaps


def cigar_indels(bam, acc, min_len=MIN_INDEL_BP):
    """CIGAR D/I events >= min_len in primary alignments."""
    import pysam
    events = []
    bf = pysam.AlignmentFile(bam, "rb")
    for r in bf.fetch(acc):
        if r.is_secondary or r.is_supplementary or r.is_unmapped:
            continue
        pos = r.reference_start
        for op, ln in r.cigartuples:
            if op == 2 and ln >= min_len:  # deletion from reference
                events.append((acc, pos, "DEL", ln))
            elif op == 1 and ln >= min_len:  # insertion
                events.append((acc, pos, "INS", ln))
            if op in (0, 2, 3, 7, 8):
                pos += ln
    bf.close()
    return events


def design_gene_primers(ref_seqs, locus, variants_by_strain, product_range=(300, 550)):
    """Primer3 design in SNP-free windows of a gene region; handles minus strand."""
    info = LOCI[locus]
    acc, s, e, strand = info["acc"], info["start"], info["end"], info["strand"]
    seq = str(ref_seqs[acc].seq[s - 1:e])
    if strand == -1:
        seq = str(Seq(seq).reverse_complement())
    var = set()
    for positions in variants_by_strain.values():
        var |= positions
    # map variant positions to template coordinates
    tvar = set()
    for p in var:
        if s <= p <= e:
            tvar.add((e - p) if strand == -1 else (p - s))  # 0-based on template
    ok = []
    w = 0
    for p in sorted(tvar):
        if p - w >= 25:
            ok.append((w, p - w))
        w = max(w, p + 1)
    if len(seq) - w >= 25:
        ok.append((w, len(seq) - w))
    if not ok:
        return None
    res = primer3.design_primers(
        seq_args={"SEQUENCE_ID": f"Cr_{locus}", "SEQUENCE_TEMPLATE": seq,
                  "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [ok[0][0], ok[0][1],
                                                          ok[-1][0], ok[-1][1]]},
        global_args={**PRIMER3_GLOBAL_ARGS,
                     "PRIMER_PRODUCT_SIZE_RANGE": [list(product_range)]})
    if "PRIMER_LEFT_0_SEQUENCE" not in res:
        return None
    return dict(F=res["PRIMER_LEFT_0_SEQUENCE"], R=res["PRIMER_RIGHT_0_SEQUENCE"],
                F_tm=round(res["PRIMER_LEFT_0_TM"], 1),
                R_tm=round(res["PRIMER_RIGHT_0_TM"], 1),
                product=res["PRIMER_PAIR_0_PRODUCT_SIZE"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--ref-bam", required=True)
    ap.add_argument("--query-bam", required=True)
    ap.add_argument("--blast-db", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ref_seqs = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))

    print("=== 1. Organellar coverage-gap screen (normalized thresholds) ===")
    for acc in ("NC_005353.1", "NC_001638.1"):
        rd = depth_array(args.ref_bam, acc)
        qd = depth_array(args.query_bam, acc)
        gaps = coverage_gap_normalized(rd, qd)
        print(f"  {acc}: {len(gaps)} candidate gap(s) >= 100 bp "
              f"(expect only repetitive IR regions in chloroplast)")

    print("=== 2. CIGAR indel scan of CC-1373 alignments (>= 30 bp) ===")
    for acc in ("NC_005353.1", "NC_001638.1"):
        ev = cigar_indels(args.query_bam, acc)
        print(f"  {acc}: {len(ev)} indel(s) >= {MIN_INDEL_BP} bp")

    print("=== 3. Variant calling in rbcL / cox1 regions (+/- 500 bp) ===")
    variants = {}
    for locus, info in LOCI.items():
        for strain, bam in (("CC124", args.ref_bam), ("CC1373", args.query_bam)):
            variants[(locus, strain)] = call_flank_variants(
                bam, args.reference, info["acc"],
                max(1, info["start"] - FLANK), info["end"] + FLANK)
            print(f"  {locus} {strain}: {len(variants[(locus, strain)])} variant(s)")

    print("=== 4. Primer3 design + blastn-short verification ===")
    rows = []
    for locus in LOCI:
        vbs = {"CC124": variants[(locus, "CC124")], "CC1373": variants[(locus, "CC1373")]}
        res = design_gene_primers(ref_seqs, locus, vbs)
        if not res:
            print(f"  {locus}: no primer pair found")
            continue
        uniq = blastn_short_uniqueness({"F": res["F"], "R": res["R"]}, args.blast_db)
        rows.append(dict(locus=locus, primer_F=res["F"], primer_R=res["R"],
                         F_tm=res["F_tm"], R_tm=res["R_tm"], product_bp=res["product"],
                         F_perfect_hits=uniq.get("F", (0, 0))[0],
                         R_perfect_hits=uniq.get("R", (0, 0))[0]))
        print(f"  {locus}: F {res['F']} / R {res['R']} -> {res['product']} bp, "
              f"perfect hits F={uniq.get('F')}, R={uniq.get('R')}")

    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    print(f"Written to {args.out}")


if __name__ == "__main__":
    main()
