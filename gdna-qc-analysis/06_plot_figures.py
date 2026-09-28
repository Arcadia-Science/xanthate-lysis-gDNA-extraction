#!/usr/bin/env python3
"""
06_plot_figures.py — draw the data panels of Figures 2, 3 and 4 from the raw
traces and the tables written by 01–04.

These are the data layers only. The published figures were laid out and
labelled in Adobe Illustrator from panels like these, so spacing, titles and
annotations differ; the plotted values are the same.

  Figure 2 (high-fidelity track)   A traces, B % ≥10 kb, C average size,
                                   D NanoDrop spectra, E Qubit, F NanoDrop:Qubit
  Figure 3 (high-throughput track) A traces, B % ≥10 kb, C average size,
                                   D NanoDrop spectra, E Qubit dsDNA,
                                   F Qubit RNA, G NanoDrop:Qubit
  Figure 4 (optimization notes)    A handling, B lysis time course, C thaw,
                                   D lysis additives

Traces are normalized to each lane's highest point above 400 bp, resampled onto
one reference run's acquisition rows by calibrated size, and drawn against row
position with size labels at their calibrated positions (the calibration is
not log-linear). Spectra get a 9-point moving average for display only.
Qubit RNA reads at or below the detection limit are drawn at 0.

Requires: numpy, pandas, scipy, matplotlib; arcadia-pycolor for house style
(tested with matplotlib 3.10.9 and arcadia-pycolor 0.8.0).

Usage:
    python 06_plot_figures.py --data-dir /path/to/gDNA_extractions \
        --results results --out figures
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import qclib as q

try:
    import arcadia_pycolor as apc

    apc.mpl.setup()
    C = {k: getattr(apc, k).hex_code for k in
         ("aegean", "dragon", "aster", "amber", "seaweed", "bark", "charcoal", "crow")}
    style = apc.mpl.style_plot
except Exception:  # plain matplotlib fallback, same hex codes
    C = dict(aegean="#5088C5", dragon="#C85152", aster="#7A77AB", amber="#F28360",
             seaweed="#3B9886", bark="#8F8885", charcoal="#484B50", crow="#292928")

    def style(ax, **_):
        ax.spines[["top", "right"]].set_visible(False)

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42

GROUP_COLORS = {"C. reinhardtii": C["aegean"], "C. smithii": C["dragon"],
                "RNA-free": C["amber"], "Single-SPRI": C["seaweed"]}
TICKS_BP = [100, 1000, 10000]


class Ctx:
    def __init__(self, data_dir, results):
        self.data_dir = data_dir
        self.samples = q.load_samples()
        self.sets = q.load_sets()
        r = Path(results)
        self.ts = pd.read_csv(r / "tapestation_metrics.csv").set_index("sample_id")
        self.nq = pd.read_csv(r / "nanodrop_qubit.csv").set_index("sample_id")
        self.spec = pd.read_csv(r / "nanodrop_spectra_normalized.csv")

    def ids(self, set_, group=None):
        return list(q.members(self.sets, set_, group)["sample_id"])

    def lane(self, sid):
        s = self.samples.loc[sid]
        return q.lane(self.data_dir, s["prep"], s["tapestation_file"], s["tapestation_well"],
                      s["tapestation_desc"])


SPECIES = ("C. reinhardtii", "C. smithii")


def italic_ticks(ax):
    for t in ax.get_xticklabels():
        if t.get_text() in SPECIES:
            t.set_fontstyle("italic")


def italic_legend(leg):
    for t in leg.get_texts():
        if t.get_text().split(" (")[0] in SPECIES:
            t.set_fontstyle("italic")


def title(ax, text, letter):
    ax.set_title(text, fontsize=14, fontweight="semibold", color=C["crow"], loc="left")
    ax.text(-0.25, 1.18, letter, transform=ax.transAxes, fontsize=28, color=C["bark"], va="top")


def points(ax, groups, ylabel, colors=None, whisker=None):
    """Individual points, mean crossbar and ± SD whisker per group."""
    labels = []
    for i, (lab, v, col) in enumerate(groups):
        v = np.asarray(v, float)
        cols = col if isinstance(col, list) else [col] * len(v)
        off = np.linspace(-0.13, 0.13, len(v)) if len(v) > 1 else np.zeros(1)
        ax.scatter(i + off, v, c=cols, alpha=0.55, s=30, lw=0, zorder=2, clip_on=False)
        f = v[np.isfinite(v)]
        m, sd = f.mean(), f.std(ddof=1)
        wc = whisker or (col if isinstance(col, str) else C["charcoal"])
        ax.errorbar([i], [m], yerr=[sd], fmt="none", ecolor=wc, elinewidth=1.6, capsize=4, zorder=3)
        ax.plot([i - 0.19, i + 0.19], [m, m], color=C["crow"], lw=2.4, zorder=4, solid_capstyle="butt")
        labels.append(lab)
    ax.set_xlim(-0.5, len(groups) - 0.5)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=0)
    style(ax)


def traces(ctx, ax, groups, ref_bp, ymax=None, ylabel="Normalized fluorescence"):
    """Mean ± SD of normalized traces per group on the reference run's rows."""
    rows = np.arange(len(ref_bp))
    keep = (ref_bp > 12) & (ref_bp < 1e6)
    top = 0
    for lab, ids, col in groups:
        tr = []
        for sid in ids:
            bp, fu = ctx.lane(sid)
            y = q.on_reference_grid(ref_bp, bp, fu)
            tr.append(y / y[ref_bp > q.METRIC_LO].max())
        tr = np.array(tr)[:, keep]
        m, sd = tr.mean(0), tr.std(0, ddof=1)
        ax.fill_between(rows[keep], np.clip(m - sd, 0, None), m + sd, color=col, alpha=0.28, lw=0)
        ax.plot(rows[keep], m, color=col, lw=1.5, label=f"{lab} (n = {len(ids)})")
        top = max(top, float((m + sd)[ref_bp[keep] > q.MARKER_SHADE[1]].max()))
    ymax = ymax or top * 1.3
    ax.axvspan(q.bp_to_row(ref_bp, q.MARKER_SHADE[0]), q.bp_to_row(ref_bp, q.MARKER_SHADE[1]),
               color=C["bark"], alpha=0.14, lw=0, zorder=0)
    ax.axvline(q.bp_to_row(ref_bp, 60000), color=C["bark"], lw=0.9, ls=":", zorder=0)
    ax.text(q.bp_to_row(ref_bp, 60000), ymax * 0.97, "60 kb ", ha="right", va="top", fontsize=10,
            style="italic", color=C["charcoal"])
    ax.text(q.bp_to_row(ref_bp, q.MARKER_SHADE[0]), ymax * 0.97, " Lower marker", ha="left",
            va="top", fontsize=10, style="italic", color=C["charcoal"])
    ax.set_xlim(rows[keep].min(), rows[keep].max())
    ax.set_ylim(0, ymax)
    ax.set_xticks([q.bp_to_row(ref_bp, b) for b in TICKS_BP])
    ax.set_xticklabels([f"{b:,}" for b in TICKS_BP])
    ax.set_xlabel("Size (bp)")
    ax.set_ylabel(ylabel)
    style(ax)


def spectra(ctx, ax, groups):
    for lab, ids, col in groups:
        sub = ctx.spec[ctx.spec["sample_id"].isin(ids)]
        w = np.sort(sub["wavelength_nm"].unique())
        a = np.array([q.smooth_edge(sub[sub["sample_id"] == i].sort_values("wavelength_nm")
                                    ["abs_norm_a260"].to_numpy()) for i in ids])
        m, sd = a.mean(0), a.std(0, ddof=1)
        ax.fill_between(w, np.clip(m - sd, 0, None), m + sd, color=col, alpha=0.28, lw=0)
        ax.plot(w, m, color=col, lw=1.5, label=lab)
    for x in (230, 260, 280):
        ax.axvline(x, color=C["bark"], lw=0.6, alpha=0.5, zorder=0)
    ax.set_xlim(220, 320)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Absorbance (norm. A$_{260}$)")
    style(ax)


def track_figure(ctx, set_, groups, fname, heading, rna=False):
    ids = {g: ctx.ids(set_, g) for g in groups}
    col = [GROUP_COLORS[g] for g in groups]
    n = 4
    fig = plt.figure(figsize=(4.2 * n, 9))
    gs = fig.add_gridspec(2, n, hspace=0.55, wspace=0.45)
    ref_bp = ctx.lane(ids[groups[0]][0])[0]
    axA = fig.add_subplot(gs[0, :n - 2])
    traces(ctx, axA, [(g, ids[g], c) for g, c in zip(groups, col)], ref_bp)
    italic_legend(axA.legend(frameon=False, fontsize=10, loc="upper center"))
    title(axA, "TapeStation gDNA electropherogram", "A")
    for k, (m, lab, t) in enumerate((("pct_ge_10kb", "DNA ≥10 kb (%)", "DNA ≥10 kb"),
                                     ("avg_size_kb", "Average size (kb)", "Average DNA size"))):
        ax = fig.add_subplot(gs[0, n - 2 + k])
        points(ax, [(g, ctx.ts.loc[ids[g], m], c) for g, c in zip(groups, col)], lab)
        italic_ticks(ax)
        title(ax, t, "BC"[k])
    axD = fig.add_subplot(gs[1, 0])
    spectra(ctx, axD, [(g, ids[g], c) for g, c in zip(groups, col)])
    title(axD, "NanoDrop spectra", "D")
    panels = [("qubit_dsdna", "Qubit dsDNA (ng/µL)", "Qubit dsDNA")]
    if rna:
        panels.append(("qubit_rna", "Qubit RNA (ng/µL)", "Qubit RNA"))
    panels.append(("nd_qubit_ratio", "NanoDrop:Qubit", "NanoDrop:Qubit ratio"))
    for k, (m, lab, t) in enumerate(panels):
        ax = fig.add_subplot(gs[1, 1 + k])
        vals = [ctx.nq.loc[ids[g], m].fillna(0) if m == "qubit_rna" else ctx.nq.loc[ids[g], m]
                for g in groups]
        points(ax, [(g, v, c) for g, v, c in zip(groups, vals, col)], lab)
        italic_ticks(ax)
        title(ax, t, "EFG"[k])
    fig.suptitle(heading, x=0.02, ha="left", fontsize=16, color=C["crow"])
    return fig


def figure4(ctx):
    fig = plt.figure(figsize=(17, 10))
    gs = fig.add_gridspec(2, 6, hspace=0.6, wspace=0.7)
    # row 1: A, B, C (four sub-panels); row 2: D (three sub-panels)
    sp = lambda sid: C["aegean"] if ctx.samples.loc[sid, "species"] == "C. reinhardtii" else C["dragon"]

    ax = fig.add_subplot(gs[0, 0])
    g = [(lab, ctx.ts.loc[ctx.ids("opt_handling", lab), "avg_size_kb"],
          [sp(i) for i in ctx.ids("opt_handling", lab)]) for lab in ("Before", "After")]
    points(ax, g, "Average size (kb)", whisker=C["charcoal"])
    title(ax, "Gentle handling (HF)", "A")

    ax = fig.add_subplot(gs[0, 1])
    t = q.members(ctx.sets, "opt_timecourse")
    times = sorted(t["group"].astype(int).unique())
    means = []
    for tm in times:
        v = ctx.ts.loc[ctx.ids("opt_timecourse", str(tm)), "avg_size_kb"].to_numpy()
        ax.scatter(tm + np.linspace(-0.6, 0.6, len(v)), v, color=C["aegean"], alpha=0.55, s=30, lw=0)
        means.append(v.mean())
    ax.plot(times, means, color=C["aegean"], lw=1.5)
    ax.set_xticks([0, 10, 20, 30, 40])
    ax.set_xlabel("Lysis time (min)")
    ax.set_ylabel("Average size (kb)")
    ax.set_ylim(0, 8)
    style(ax)
    title(ax, "Hot dry thaw (HT)", "B")

    for k, (set_, m, lab, tt) in enumerate((
            ("opt_thaw_ht", "pct_ge_10kb", "DNA ≥10 kb (%)", "HT, 10 min"),
            ("opt_thaw_ht", "avg_size_kb", "Average size (kb)", "HT, 10 min"),
            ("opt_thaw_hf", "pct_ge_10kb", "DNA ≥10 kb (%)", "HF"),
            ("opt_thaw_hf", "avg_size_kb", "Average size (kb)", "HF"))):
        ax = fig.add_subplot(gs[0, 2 + k])
        points(ax, [(g.replace(" ", "\n"), ctx.ts.loc[ctx.ids(set_, g), m], C["aegean"])
                    for g in ("Hot dry", "Cold wet")], lab)
        title(ax, tt, "C" if k == 0 else "")

    t = q.members(ctx.sets, "opt_additives")
    conds = list(dict.fromkeys(t["group"]))
    for k, (m, lab) in enumerate((("a260_a230", "A$_{260}$/A$_{230}$"),
                                  ("nd_qubit_ratio", "NanoDrop:Qubit"),
                                  ("qubit_dsdna", "Qubit (ng/µL)"))):
        ax = fig.add_subplot(gs[1, k])
        for pair, pt in t.groupby("pair", sort=False):
            pt = pt.set_index("group").loc[conds]
            v = ctx.nq.loc[pt["sample_id"], m].to_numpy()
            col = sp(pt["sample_id"].iloc[0])
            ax.plot(range(3), v, color=col, lw=1, alpha=0.4)
            ax.scatter(range(3), v, color=col, s=30, alpha=0.55, lw=0)
        for i, c in enumerate(conds):
            mu = ctx.nq.loc[ctx.ids("opt_additives", c), m].mean()
            ax.plot([i - 0.25, i + 0.25], [mu, mu], color=C["crow"], lw=2.4)
        ax.set_xticks(range(3))
        ax.set_xticklabels(["Lysis\nbuffer", "+PK", "+PK\n+RNase"])
        ax.set_ylabel(lab)
        ax.set_ylim(bottom=0)
        style(ax)
        title(ax, "Lysis additives (HT)" if k == 0 else "", "D" if k == 0 else "")
    fig.suptitle("Figure 4 data panels: handling, thaw and lysis-additive tests", x=0.02,
                 ha="left", fontsize=16, color=C["crow"])
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default="figures")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    ctx = Ctx(args.data_dir, args.results)
    figs = {
        "figure2_high_fidelity": track_figure(ctx, "hf_track", ["C. reinhardtii", "C. smithii"],
                                              "figure2", "Figure 2 data panels: high-fidelity track"),
        "figure3_high_throughput": track_figure(ctx, "ht_track", ["RNA-free", "Single-SPRI"],
                                                "figure3", "Figure 3 data panels: high-throughput track",
                                                rna=True),
        "figure4_optimization": figure4(ctx),
    }
    for name, fig in figs.items():
        for ext in ("svg", "pdf", "png"):
            fig.savefig(out / f"{name}.{ext}", bbox_inches="tight", facecolor="white",
                        **({"dpi": 200} if ext == "png" else {}))
        plt.close(fig)
        print(f"wrote {out}/{name}.svg/.pdf/.png")


if __name__ == "__main__":
    main()
