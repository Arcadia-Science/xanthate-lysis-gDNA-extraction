#!/usr/bin/env python3
"""
03_lower_marker_comparison.py — compare DNA amount at or above 10 kb between
the two high-throughput cleanups using the TapeStation lower marker as the
per-lane loading reference.

Qubit dsDNA can be inflated by RNA, so the ratio of each lane's signal at or
above 10 kb (10 kb–1 Mb) to its lower-marker signal (50–150 bp; the 100 bp
marker is spiked into every well at a fixed amount) gives an RNA-independent
comparison of DNA amount. Reads results/tapestation_metrics.csv from 01.

Supports "By TapeStation, the two cleanups carried similar amounts of DNA at or
above 10 kb relative to the lower marker" (High-throughput track) and the last
sentence of "TapeStation trace analysis".

Requires: numpy, pandas, scipy (see requirements.txt).

Usage:
    python 03_lower_marker_comparison.py --results results
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import qclib as q


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--results", default="results")
    ap.add_argument("--set", default="ht_track")
    args = ap.parse_args()
    res = Path(args.results)

    ts = pd.read_csv(res / "tapestation_metrics.csv").set_index("sample_id")
    sets = q.load_sets()
    groups = {}
    rows = []
    for g in q.members(sets, args.set)["group"].unique():
        ids = q.members(sets, args.set, g)["sample_id"]
        v = ts.loc[ids, "ge10kb_per_lower_marker"].to_numpy()
        groups[g] = v
        d = q.describe(v)
        rows.append(dict(set=args.set, group=g, **d, median=float(np.median(v))))
    summ = pd.DataFrame(rows)
    (a, va), (b, vb) = list(groups.items())[:2]
    p = mannwhitneyu(va, vb, alternative="two-sided").pvalue
    summ["ratio_of_means_vs_first"] = summ["mean"] / summ["mean"].iloc[0]
    summ.to_csv(res / "lower_marker_comparison.csv", index=False, float_format="%.4f")

    print(summ[["group", "n", "mean", "sd", "median", "min", "max"]].to_string(index=False))
    print(f"\n{b} / {a} ratio of means: {vb.mean() / va.mean():.2f}")
    print(f"Mann-Whitney U, two-sided (descriptive; different strains per group): p = {p:.2f}")


if __name__ == "__main__":
    main()
