#!/usr/bin/env python3
"""
02_hsd5000_low_mw.py — check the High Sensitivity D5000 traces for sample DNA
that the Genomic DNA ScreenTape's lower marker could hide.

The Genomic DNA ScreenTape's 100 bp lower marker saturates the detector across
roughly 35–400 bp. The HS D5000 ScreenTape places its lower marker at 15 bp, so
the same region is visible. Each run is sized against its own HS D5000 ladder
(11 bands, 15 bp–10 kb; PCHIP in log10(size), linear in log10(size) beyond the
end bands). Per lane, both markers are removed (everything below 60 bp, which
holds the 15 bp marker and its tail, and 9–11 kb around the 10 kb upper
marker), and the remaining signal is the denominator for:

  pct_70_150   share of lane signal at 70–150 bp (where the gDNA marker sits)
  pct_60_400   share of lane signal at 60–400 bp (the full masked span)

Feeds the High Sensitivity D5000 sentence in "High-throughput track".

Requires: numpy, pandas, scipy (see requirements.txt).

Usage:
    python 02_hsd5000_low_mw.py --data-dir /path/to/gDNA_extractions --out results
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import qclib as q


def shares(bp, fu):
    y = np.clip(fu, 0, None)
    lo_u, hi_u = q.HSD_UPPER_MARKER
    keep = (bp >= q.HSD_LOWER_MARKER_MAX) & ~((bp >= lo_u) & (bp <= hi_u))
    total = y[keep].sum()
    out = {}
    for name, (lo, hi) in q.HSD_WINDOWS.items():
        m = keep & (bp >= lo) & (bp <= hi)
        out[f"pct_{name}"] = 100 * y[m].sum() / total
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    samples = q.load_samples()
    rows, ladders = [], {}
    for sid, s in samples[samples["hsd5000_file"] != ""].iterrows():
        path = str(q.find_file(args.data_dir, s["prep"], s["hsd5000_file"]))
        run = q.load_run(path, "HSD5000")
        ladders[s["hsd5000_file"]] = dict(run=s["hsd5000_file"], n_bands=len(run["band_rows"]),
                                          band_rows=" ".join(map(str, run["band_rows"])))
        # HS D5000 lanes are labelled with the plate well, which is the sample name
        bp, fu = q.lane(args.data_dir, s["prep"], s["hsd5000_file"], s["hsd5000_well"],
                        s["sample"], assay="HSD5000")
        rows.append(dict(sample_id=sid, run=s["hsd5000_file"], well=s["hsd5000_well"], **shares(bp, fu)))

    res = pd.DataFrame(rows)
    res.to_csv(out / "hsd5000_low_mw.csv", index=False, float_format="%.4f")
    pd.DataFrame(ladders.values()).to_csv(out / "hsd5000_ladders.csv", index=False)
    print(f"{len(res)} lanes -> {out}/hsd5000_low_mw.csv")


if __name__ == "__main__":
    main()
