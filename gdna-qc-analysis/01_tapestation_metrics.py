#!/usr/bin/env python3
"""
01_tapestation_metrics.py — size every Genomic DNA ScreenTape lane against its
own ladder and compute the trace metrics reported in the pub.

For each run listed in metadata/samples.csv, the ladder lane's 14 bands
(100 bp–48.5 kb) are located as peaks and log10(size) is fitted to acquisition
row by monotone piecewise cubic (PCHIP) interpolation, extended linearly in
log10(size) beyond the end bands. Per lane, using signal between 400 bp and
1 Mb (above the lower marker and below the baseline at the extrapolated top of
the size scale):

  avg_size_kb      fluorescence-weighted mean fragment size
  pct_ge_10kb      percentage of that signal at or above 10 kb
  pct_gt_ladder    percentage of that signal above the ladder's top band
                   (48.5 kb), i.e. signal sized only by extrapolation

and, for the lower-marker comparison (03):

  signal_ge_10kb   summed fluorescence from 10 kb to 1 Mb
  lower_marker     summed fluorescence from 50 to 150 bp (the 100 bp marker)

Feeds Figure 2 A–C, Figure 3 A–C, Figure 4 A–C and the matching text.

Requires: numpy, pandas, scipy (see requirements.txt).

Usage:
    python 01_tapestation_metrics.py --data-dir /path/to/gDNA_extractions --out results
"""

import argparse
from pathlib import Path

import pandas as pd

import qclib as q


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--data-dir", required=True, help="folder holding the <prep>_<date>/ folders")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    samples = q.load_samples()
    ts = samples[samples["tapestation_file"] != ""]

    rows, ladders = [], {}
    for sid, s in ts.iterrows():
        path = str(q.find_file(args.data_dir, s["prep"], s["tapestation_file"]))
        run = q.load_run(path, "gDNA")
        if s["tapestation_file"] not in ladders:
            ladders[s["tapestation_file"]] = dict(
                prep=s["prep"], run=s["tapestation_file"], threshold=run["threshold"],
                n_bands=len(run["band_rows"]), band_rows=" ".join(map(str, run["band_rows"])),
                row_60kb=round(q.bp_to_row(run["bp"], 60_000), 1),
                max_size_bp=round(float(run["bp"].max())))
        bp, fu = q.lane(args.data_dir, s["prep"], s["tapestation_file"], s["tapestation_well"],
                        s["tapestation_desc"])
        hmw = q.window_sum(bp, fu, q.HMW_CUTOFF, q.METRIC_HI)
        lm = q.window_sum(bp, fu, *q.LOWER_MARKER)
        rows.append(dict(sample_id=sid, prep=s["prep"], run=s["tapestation_file"],
                         well=s["tapestation_well"], description=s["tapestation_desc"],
                         avg_size_kb=q.avg_size_kb(bp, fu), pct_ge_10kb=q.pct_ge(bp, fu),
                         pct_gt_ladder=q.pct_ge(bp, fu, cutoff=q.GDNA_LADDER[-1]),
                         signal_ge_10kb=hmw, lower_marker=lm, ge10kb_per_lower_marker=hmw / lm))

    res = pd.DataFrame(rows)
    res.to_csv(out / "tapestation_metrics.csv", index=False, float_format="%.4f")
    lad = pd.DataFrame(ladders.values())
    lad.to_csv(out / "tapestation_ladders.csv", index=False)
    bad = lad[lad["n_bands"] != len(q.GDNA_LADDER)]
    if len(bad):
        raise SystemExit(f"ladder calibration failed for: {list(bad['run'])}")
    print(f"{len(res)} lanes from {len(lad)} runs -> {out}/tapestation_metrics.csv")
    print(lad[["run", "threshold", "row_60kb"]].to_string(index=False))


if __name__ == "__main__":
    main()
