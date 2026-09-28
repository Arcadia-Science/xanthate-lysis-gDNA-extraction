#!/usr/bin/env python3
"""
03_analyze_reads.py — on-target rates and genotype calls from mapped amplicon reads.

For each sample BAM (from 02_map_reads.sh):
  1. On-target rate: fraction of primary mapped reads overlapping the intended
     locus (+/- 200 bp padding), plus the contig with the most primary reads
     (detects off-target amplification).
  2. chr9 deletion-marker genotype (reaction A): read depth across the interior
     of the CC-1373 deletion (NC_057012.1:7,368,000-7,368,800). A full-length
     CC-124 allele covers the interior; the CC-1373 deletion allele leaves it
     empty. Reported: depth at the deletion midpoint, and mean/min depth across
     the interior window (samtools depth).
  3. cox1 genotype (reaction E): per-base counts at the two diagnostic SNPs,
     NC_001638.1:5,911 (C in CC-124, T in CC-1373) and 5,929 (C -> A).

Expects BAMs named <order>_<num>_<reaction>_<prep>_<sample>.bam
(e.g. AI-8632_1_A_HF-7_cc124_A3.bam).

Requires: pysam, pandas (pip); samtools on PATH (tested v1.22.1).

Usage:
    python 03_analyze_reads.py --bam-dir bams --out results
"""

import argparse
import glob
import os
import subprocess
from collections import Counter

import pandas as pd
import pysam

# Intended locus per reaction (RefSeq accession, start, end, label).
# Coordinates are the designed amplicon intervals on the v5.5 reference.
# The validated panel is A (chr9) + D (rbcL) + E (cox1). Reaction C (chr13
# candidate) is kept here because the sequenced run contains those samples;
# the on-target table is what revealed its off-target amplification.
INTENDED = {
    "A": ("NC_057012.1", 7367489, 7369276, "chr9 deletion marker"),
    "C": ("NC_057016.1", 4685483, 4687445, "chr13 deletion marker (failed candidate)"),
    "D": ("NC_005353.1", 123488, 123876, "chloroplast rbcL"),
    "E": ("NC_001638.1", 5729, 6183, "mitochondrion cox1"),
}
LOCUS_PAD = 200  # bp of padding around the intended locus for on-target calls

# chr9 deletion interior: sequence present in CC-124, absent in CC-1373.
CHR9_ACC = "NC_057012.1"
CHR9_DELETION = (7367489, 7369276)      # ~1.8 kb CC-1373 deletion (approx. breakpoints)
CHR9_INTERIOR = (7368000, 7368800)      # interior window used for genotyping
CHR9_DEL_MEAN_THRESHOLD = 10            # mean interior depth below this => deletion allele

# Diagnostic cox1 SNPs (1-based on NC_001638.1): CC-124 base -> CC-1373 base.
COX1_ACC = "NC_001638.1"
COX1_SNPS = {5911: ("C", "T"), 5929: ("C", "A")}

# v5.5 accession -> chromosome label (NC_057004.1 = chr1 ... NC_057020.1 = chr17)
ACC2CHR = {f"NC_0570{n:02d}.1": f"chr{n - 3}" for n in range(4, 21)}
ACC2CHR.update({"NC_005353.1": "chloroplast", "NC_001638.1": "mitochondrion"})


def parse_sample(bam_path):
    """<order>_<num>_<reaction>_<prep>_<sample>.bam -> (num, reaction, sample)."""
    base = os.path.basename(bam_path).replace(".bam", "")
    parts = base.split("_")
    return int(parts[1]), parts[2], "_".join(parts[4:])


def primary_count(bf):
    return bf.count(read_callback=lambda r: not (
        r.is_secondary or r.is_supplementary or r.is_unmapped))


def on_target_table(bam_dir):
    rows = []
    for bam in sorted(glob.glob(os.path.join(bam_dir, "*.bam"))):
        num, rxn, sample = parse_sample(bam)
        acc, s, e, locus = INTENDED[rxn]
        bf = pysam.AlignmentFile(bam, "rb")
        total = primary_count(bf)
        on = sum(1 for r in bf.fetch(acc, max(0, s - LOCUS_PAD), e + LOCUS_PAD)
                 if not (r.is_secondary or r.is_supplementary or r.is_unmapped))
        bf.close()
        # contig tally uses a fresh handle: a regional fetch() advances the file
        # position, so reusing the same handle for until_eof would miss earlier contigs
        bf2 = pysam.AlignmentFile(bam, "rb")
        counts = Counter()
        for r in bf2.fetch(until_eof=True):
            if not (r.is_secondary or r.is_supplementary or r.is_unmapped):
                counts[r.reference_name] += 1
        top = counts.most_common(1)[0] if counts else ("none", 0)
        bf2.close()
        rows.append(dict(num=num, rxn=rxn, sample=sample, locus=locus,
                         total_reads=total, on_target=on,
                         pct_on=round(100 * on / total, 1) if total else 0.0,
                         top_contig=ACC2CHR.get(top[0], top[0]), top_n=top[1]))
    return pd.DataFrame(rows).sort_values(["rxn", "num"]).reset_index(drop=True)


def depth_series(bam, acc, start, end):
    """Per-base depth over acc:start-end (1-based, inclusive) via samtools depth -a."""
    out = subprocess.run(
        ["samtools", "depth", "-a", "-r", f"{acc}:{start}-{end}", bam],
        capture_output=True, text=True, check=True).stdout
    depths = {}
    for line in out.splitlines():
        _, pos, d = line.split("\t")
        depths[int(pos)] = int(d)
    return [depths.get(p, 0) for p in range(start, end + 1)]


def chr9_genotypes(bam_dir):
    rows = []
    for bam in sorted(glob.glob(os.path.join(bam_dir, "*.bam"))):
        num, rxn, sample = parse_sample(bam)
        if rxn != "A":
            continue
        s, e = CHR9_INTERIOR
        interior = depth_series(bam, CHR9_ACC, s, e)
        mean_d = sum(interior) / len(interior)
        # deletion midpoint, NC_057012.1:7,368,383 (1-based, as in the methods text)
        mid = (CHR9_DELETION[0] + CHR9_DELETION[1] + 1) // 2
        mid_d = depth_series(bam, CHR9_ACC, mid, mid)[0]
        call = "CC-1373 (deletion)" if mean_d < CHR9_DEL_MEAN_THRESHOLD else "CC-124 (full-length)"
        rows.append(dict(num=num, sample=sample, midpoint_depth=mid_d,
                         interior_mean_depth=round(mean_d, 1),
                         interior_min_depth=min(interior),
                         interior_breadth=round(sum(d > 0 for d in interior) / len(interior), 4),
                         genotype=call))
    return pd.DataFrame(rows).sort_values("num").reset_index(drop=True)


def base_counts(bam, acc, pos1):
    """Primary-read base counts at a 1-based reference position."""
    counts = Counter()
    bf = pysam.AlignmentFile(bam, "rb")
    for col in bf.pileup(acc, pos1 - 1, pos1, truncate=True, min_base_quality=0):
        if col.reference_pos != pos1 - 1:
            continue
        for pr in col.pileups:
            r = pr.alignment
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            counts[r.query_sequence[pr.query_position]] += 1
    bf.close()
    return counts


def cox1_genotypes(bam_dir):
    rows = []
    for bam in sorted(glob.glob(os.path.join(bam_dir, "*.bam"))):
        num, rxn, sample = parse_sample(bam)
        if rxn != "E":
            continue
        calls = []
        observed = {}
        for pos, (ref124, alt1373) in sorted(COX1_SNPS.items()):
            c = base_counts(bam, COX1_ACC, pos)
            total = sum(c.values())
            base, n = c.most_common(1)[0] if total else ("N", 0)
            conc = round(100 * n / total, 1) if total else 0.0
            calls.append(f"{pos}:{base} ({conc}%)")
            observed[pos] = base
            rows.append(dict(num=num, sample=sample, snp_pos=pos,
                             cc124_base=ref124, cc1373_base=alt1373,
                             observed=base, pct_concordant=conc,
                             counts=dict(c)))
        genotype = ("CC-124" if all(observed[p] == v[0] for p, v in COX1_SNPS.items())
                    else "CC-1373" if all(observed[p] == v[1] for p, v in COX1_SNPS.items())
                    else "ambiguous")
        print(f"  {sample}: {'; '.join(calls)} -> {genotype}")
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam-dir", required=True, help="directory of BAMs from 02_map_reads.sh")
    ap.add_argument("--out", default="results", help="output directory for CSV tables")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    print("=== On-target rates (all samples) ===")
    ont = on_target_table(args.bam_dir)
    print(ont.to_string(index=False))
    ont.to_csv(os.path.join(args.out, "on_target_table.csv"), index=False)

    print("\n=== chr9 deletion-marker genotypes (reaction A) ===")
    chr9 = chr9_genotypes(args.bam_dir)
    print(chr9.to_string(index=False))
    chr9.to_csv(os.path.join(args.out, "chr9_genotypes.csv"), index=False)

    print("\n=== cox1 SNP genotypes (reaction E) ===")
    cox1 = cox1_genotypes(args.bam_dir)
    if len(cox1):
        cox1.to_csv(os.path.join(args.out, "cox1_genotypes.csv"), index=False)

    print(f"\nTables written to {args.out}/")


if __name__ == "__main__":
    main()
