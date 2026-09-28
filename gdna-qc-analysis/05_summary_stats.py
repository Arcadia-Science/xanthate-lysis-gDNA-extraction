#!/usr/bin/env python3
"""
05_summary_stats.py — every mean ± SD, range and n quoted in "Observed
performance", "Optimization notes" and the figure legends, recomputed from the
per-sample tables written by 01, 02 and 04, and checked against the values as
written in the pub (metadata/pub_values.csv).

Means are arithmetic; SD is the sample standard deviation (n − 1). Qubit RNA
reads at or below the detection limit are counted, not averaged. Values are
formatted to the decimals the pub uses, and each row reports whether the
formatted value matches the pub.

Also writes group_summary.csv: n, mean, SD, min and max of every metric for
every group in metadata/sets.csv.

Requires: numpy, pandas (see requirements.txt).

Usage:
    python 05_summary_stats.py --results results
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import qclib as q

TABLES = {"tapestation": "tapestation_metrics.csv", "hsd5000": "hsd5000_low_mw.csv",
          "nanodrop_qubit": "nanodrop_qubit.csv"}


def values(tab, sets, set_, group, metric):
    ids = q.members(sets, set_, group or None)["sample_id"]
    return tab.loc[ids, metric].to_numpy(dtype=float)


def fmt(x, d):
    return f"{x:.{int(d)}f}"


def evaluate(c, tabs, sets):
    tab = tabs[c["table"]]
    v = values(tab, sets, c["set"], c["group"], c["metric"])
    d = int(c["decimals"])
    st = c["stat"]
    if st == "count_below":
        status = tab.loc[q.members(sets, c["set"], c["group"])["sample_id"], "qubit_rna_status"]
        return f"{int((status == 'below_range').sum())} of {len(status)}", ""
    if st == "values":
        ok = v[np.isfinite(v)]
        return ", ".join(fmt(x, d) for x in ok), ""
    if st in ("ratio_of_means", "group_ratio", "pct_drop_of_means"):
        if st == "ratio_of_means":
            r = np.nanmean(v) / np.nanmean(values(tab, sets, c["set"], c["group"], c["metric2"]))
            return fmt(r, d) + "-fold", ""
        v2 = values(tab, sets, c["set"], c["group2"], c["metric"])
        if st == "group_ratio":
            return fmt(np.mean(v) / np.mean(v2), d) + f" ({c['group']} / {c['group2']})", ""
        return fmt(100 * (1 - np.mean(v2) / np.mean(v)), d), ""
    s = q.describe(v)
    detail = f"n={s['n']} mean={s['mean']:.3f} sd={s['sd']:.3f} min={s['min']:.3f} max={s['max']:.3f}"
    if st == "mean_sd":
        return f"{fmt(s['mean'], d)} ± {fmt(s['sd'], d)}", detail
    if st == "range":
        return f"{fmt(s['min'], d)}–{fmt(s['max'], d)}", detail
    if st == "mean":
        return fmt(s["mean"], d), detail
    if st == "n":
        return str(s["n"]), detail
    raise ValueError(st)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--results", default="results")
    args = ap.parse_args()
    res = Path(args.results)
    tabs = {k: pd.read_csv(res / f).set_index("sample_id") for k, f in TABLES.items()}
    sets = q.load_sets()
    claims = pd.read_csv(q.META / "pub_values.csv", dtype=str, keep_default_na=False)

    rows = []
    for _, c in claims.iterrows():
        computed, detail = evaluate(c, tabs, sets)
        rep = c["reported"].replace(" ", " ").strip()
        if not rep or not any(ch.isdigit() for ch in rep):
            match = "check by eye"
        else:
            match = "yes" if computed == rep else "NO"
        rows.append(dict(claim_id=c["claim_id"], where=c["where"], group=c["group"],
                         metric=c["metric"], stat=c["stat"], reported=c["reported"],
                         computed=computed, match=match, detail=detail))
    out = pd.DataFrame(rows)
    out.to_csv(res / "summary_stats.csv", index=False)

    grp = []
    for (set_, group), t in sets.groupby(["set", "group"], sort=False):
        for table, tab in tabs.items():
            ids = [i for i in t["sample_id"] if i in tab.index]
            if not ids:
                continue
            for m in tab.select_dtypes("number").columns:
                v = tab.loc[ids, m].to_numpy(dtype=float)
                if np.isfinite(v).any():
                    grp.append(dict(set=set_, group=group, table=table, metric=m, **q.describe(v)))
    pd.DataFrame(grp).to_csv(res / "group_summary.csv", index=False, float_format="%.4f")

    pd.set_option("display.width", 200, "display.max_colwidth", 60)
    print(out[["claim_id", "reported", "computed", "match"]].to_string(index=False))
    bad = out[out["match"] == "NO"]
    print(f"\n{len(out)} values checked; {len(bad)} differ from the pub text")


if __name__ == "__main__":
    main()
