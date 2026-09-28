#!/usr/bin/env python3
"""
design/02_find_strain_deletions.py — find large strain-specific deletions from WGS coverage.

Two-step screen used to discover the nuclear deletion-marker candidates (the
validated chr9 marker in the final panel came from this screen):

  1. Coverage-gap analysis over mosdepth 100 bp windows. Nuclear screen (flat
     thresholds): a candidate CC-1373 deletion is a run of windows with CC-124
     coverage >= 20x and CC-1373 coverage <= 5x, with sharp breakpoints (CC-1373
     flanking coverage >= 10x on both sides). An organellar/normalized mode
     (>= 20% of the strain's median vs <= 10% of its median) is provided for
     genomes with uneven or low coverage.
  2. Bridging read-pair confirmation: count proper pairs in the CC-1373 BAM with
     one mate on each side of the candidate gap (these span the deletion
     junction). Genuine deletions have many such pairs in CC-1373 and zero in
     CC-124.

Requires: pysam, pandas, numpy (pip); mosdepth output from 01_download_and_align_wgs.sh.

Usage:
    python design/02_find_strain_deletions.py \
        --ref-windows workdir/mosdepth/CC124.regions.bed.gz \
        --query-windows workdir/mosdepth/CC1373.regions.bed.gz \
        --ref-bam workdir/bam/CC124.bam --query-bam workdir/bam/CC1373.bam \
        --mode nuclear --min-gap-bp 500 --out deletions_CC1373.csv
"""

import argparse
import gzip

import numpy as np
import pandas as pd
import pysam


def load_windows(path):
    """mosdepth --by 100 regions.bed.gz -> per-chrom dict of (starts, means)."""
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start", "end", "mean"],
                     compression="gzip" if path.endswith(".gz") else None)
    out = {}
    for chrom, sub in df.groupby("chrom"):
        out[chrom] = (sub["start"].to_numpy(), sub["mean"].to_numpy())
    return out


def find_gaps(starts, ref_mean, query_mean, ref_thresh, query_thresh,
              flank_thresh, min_gap_bp, window=100):
    """Runs of windows with ref >= ref_thresh and query <= query_thresh,
    requiring query flanking coverage >= flank_thresh on both sides (sharp
    breakpoints). Returns list of (start, end, ref_mean, query_mean)."""
    gaps = []
    in_gap = False
    gap_start = 0
    for i in range(len(starts)):
        is_gap = ref_mean[i] >= ref_thresh and query_mean[i] <= query_thresh
        if is_gap and not in_gap:
            gap_start, in_gap = i, True
        elif not is_gap and in_gap:
            in_gap = False
            lo, hi = gap_start, i
            if (hi - lo) * window < min_gap_bp:
                continue
            # sharp breakpoints: query coverage recovers on both flanks
            left_ok = lo > 0 and query_mean[lo - 1] >= flank_thresh
            right_ok = hi < len(query_mean) and query_mean[hi] >= flank_thresh
            if left_ok and right_ok:
                gaps.append((int(starts[lo]), int(starts[hi - 1] + window),
                             float(ref_mean[lo:hi].mean()),
                             float(query_mean[lo:hi].mean())))
    return gaps


def bridging_pairs(bam_path, chrom, gap_start, gap_end, pad=1000):
    """Proper pairs with one mate left of the gap and the other right of it
    (junction-spanning pairs in the deletion-carrying strain)."""
    n = 0
    bf = pysam.AlignmentFile(bam_path, "rb")
    try:
        for r in bf.fetch(chrom, max(0, gap_start - pad), gap_end + pad):
            if not (r.is_read1 and r.is_proper_pair) or r.is_unmapped or r.mate_is_unmapped:
                continue
            if r.is_secondary or r.is_supplementary:
                continue
            if r.reference_start < gap_start and r.next_reference_start >= gap_end:
                n += 1
    except ValueError:
        pass
    bf.close()
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref-windows", required=True)
    ap.add_argument("--query-windows", required=True)
    ap.add_argument("--ref-bam", required=True)
    ap.add_argument("--query-bam", required=True)
    ap.add_argument("--mode", choices=["nuclear", "normalized"], default="nuclear",
                    help="nuclear: flat >=20x/<=5x with >=10x flanks; "
                         "normalized: >=20%% / <=10%% of each strain's median")
    ap.add_argument("--min-gap-bp", type=int, default=500)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ref = load_windows(args.ref_windows)
    query = load_windows(args.query_windows)

    rows = []
    for chrom, (starts, ref_mean) in ref.items():
        if chrom not in query:
            continue
        q_starts, q_mean = query[chrom]
        if not np.array_equal(starts, q_starts):
            continue
        if args.mode == "nuclear":
            rt, qt, ft = 20.0, 5.0, 10.0
        else:
            rt = np.median(ref_mean[ref_mean > 0]) * 0.20
            qt = np.median(q_mean[q_mean > 0]) * 0.10
            ft = qt * 2  # flanks must clearly recover above the gap threshold
        for (s, e, rm, qm) in find_gaps(starts, ref_mean, q_mean, rt, qt, ft,
                                        args.min_gap_bp):
            rows.append(dict(chrom=chrom, start=s, end=e, length_bp=e - s,
                             ref_mean_cov=round(rm, 1), query_mean_cov=round(qm, 1)))

    df = pd.DataFrame(rows)
    if len(df):
        print(f"{len(df)} candidate deletions; confirming with bridging read pairs ...")
        df["bridging_pairs_query"] = df.apply(
            lambda r: bridging_pairs(args.query_bam, r.chrom, r.start, r.end), axis=1)
        df["bridging_pairs_ref"] = df.apply(
            lambda r: bridging_pairs(args.ref_bam, r.chrom, r.start, r.end), axis=1)
        df["confirmed"] = (df.bridging_pairs_query > 0) & (df.bridging_pairs_ref == 0)
    else:
        print("No candidate deletions found.")
    df.to_csv(args.out, index=False)
    print(df.to_string(index=False))
    print(f"Written to {args.out}")


if __name__ == "__main__":
    main()
