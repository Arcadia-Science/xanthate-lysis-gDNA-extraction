"""
qclib.py — shared readers, calibration, and metrics for the gDNA QC analysis.

Every number comes from raw instrument exports; no vendor-computed metric
(peak tables, DIN, calibrated concentrations) is read anywhere:

  TapeStation  <run>_Electropherogram.csv   raw fluorescence per lane, by row
  NanoDrop     dsDNA <date> <time>.tsv       long-format spectra export
  Qubit        QubitData_<date>_<time>.csv   instrument export
               <prep>_Qubit_from-Notion.csv  values transcribed from the lab
                                             notebook where no export exists

Data files are located by file name under data_dir/<prep>/ (any depth), so the
scripts work on the Zenodo archive or on a copy of the Drive folder regardless
of sub-folder naming.

Requires: numpy, pandas, scipy (tested with numpy 2.4.4, pandas 3.0.2, scipy 1.17.1).
"""

import csv
import re
import unicodedata
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.signal import find_peaks

REPO = Path(__file__).resolve().parent
META = REPO / "metadata"

# Ladder band sizes (bp)
GDNA_LADDER = [100, 250, 400, 600, 900, 1200, 1500, 2000, 2500, 3000, 4000, 7000, 15000, 48500]
HSD5000_LADDER = [15, 100, 250, 400, 600, 1000, 1500, 2500, 3500, 5000, 10000]

# Genomic DNA ScreenTape trace metrics
METRIC_LO, METRIC_HI = 400, 1e6       # bp window for average size and % >= 10 kb
HMW_CUTOFF = 10_000                   # bp
LOWER_MARKER = (50, 150)              # bp window integrated as the 100 bp lower marker
MARKER_SHADE = (35, 400)              # bp span the saturated lower marker occupies (plots only)

# HS D5000 ScreenTape: marker exclusion and low-molecular-weight windows
HSD_LOWER_MARKER_MAX = 60             # bp; the 15 bp marker's tail is excluded below this
HSD_UPPER_MARKER = (9000, 11000)      # bp; the 10 kb upper marker
HSD_WINDOWS = {"70_150": (70, 150), "60_400": (60, 400)}


# --------------------------------------------------------------------------
# File lookup and metadata
# --------------------------------------------------------------------------
def _nfc(s):
    return unicodedata.normalize("NFC", s)


def find_file(data_dir, prep, filename):
    """Return the single file called `filename` anywhere under data_dir/prep."""
    root = Path(data_dir) / prep
    hits = [p for p in root.rglob("*") if p.is_file() and _nfc(p.name) == _nfc(filename)]
    if len(hits) != 1:
        raise FileNotFoundError(f"expected one '{filename}' under {root}, found {len(hits)}")
    return hits[0]


def load_samples(meta_dir=META):
    s = pd.read_csv(Path(meta_dir) / "samples.csv", dtype=str, keep_default_na=False)
    if s["sample_id"].duplicated().any():
        raise ValueError("duplicate sample_id in samples.csv")
    return s.set_index("sample_id", drop=False)


def load_sets(meta_dir=META):
    t = pd.read_csv(Path(meta_dir) / "sets.csv", dtype=str, keep_default_na=False)
    t["order"] = t["order"].astype(int)
    return t


# --------------------------------------------------------------------------
# TapeStation
# --------------------------------------------------------------------------
def read_electropherogram(path):
    """Parse an Electropherogram.csv export.

    Header cells are '<well>: <sample description>'; each following row is one
    acquisition position (no size column). Returns (wells, descriptions, data)
    with data shaped (rows, lanes).
    """
    lines = [l for l in open(path, encoding="latin-1").read().splitlines() if l.strip()]
    hdr = next(csv.reader([lines[0]]))
    wells = [h.split(":", 1)[0].strip() for h in hdr]
    descs = [h.split(":", 1)[1].strip() for h in hdr]
    data = np.array([[float(x) if x else 0.0 for x in next(csv.reader([l]))] for l in lines[1:]])
    return wells, descs, data


def calibrate_ladder(trace, bands, thresholds=(0.25, 0.20, 0.15, 0.10, 0.05), distance=5):
    """Map acquisition row -> fragment size (bp) from the ladder lane.

    Ladder bands are found as peaks, stepping the height threshold (fraction of
    the lane maximum) down until at least len(bands) peaks are found; if more
    are found, the tallest len(bands) are kept. log10(size) is fitted to row by
    PCHIP through the bands and extended linearly in log10(size) beyond the
    first and last bands (PCHIP extrapolation is not monotone).

    Returns (bp_per_row, band_rows, threshold_used).
    """
    n = len(bands)
    for h in thresholds:
        pk, _ = find_peaks(trace, height=trace.max() * h, distance=distance)
        if len(pk) >= n:
            break
    else:
        raise ValueError(f"ladder: found {len(pk)} bands, expected {n}")
    if len(pk) > n:
        pk = np.sort(pk[np.argsort(trace[pk])[-n:]])
    ly = np.log10(bands)
    core = PchipInterpolator(pk, ly, extrapolate=False)
    rows = np.arange(len(trace), dtype=float)
    v = core(rows)
    lo, hi = rows < pk[0], rows > pk[-1]
    v[lo] = ly[0] + (rows[lo] - pk[0]) * (ly[1] - ly[0]) / (pk[1] - pk[0])
    v[hi] = ly[-1] + (rows[hi] - pk[-1]) * (ly[-1] - ly[-2]) / (pk[-1] - pk[-2])
    return 10 ** v, pk, h


@lru_cache(maxsize=None)
def load_run(path, assay="gDNA"):
    """Calibrate one TapeStation run against its own ladder lane.

    Returns dict with 'bp' (size per row), 'lanes' {well: (description, trace)},
    'band_rows' and 'threshold'.
    """
    wells, descs, data = read_electropherogram(path)
    if "Ladder" not in descs:
        raise ValueError(f"no ladder lane in {path}")
    bands = GDNA_LADDER if assay == "gDNA" else HSD5000_LADDER
    thresholds = (0.25, 0.20, 0.15, 0.10) if assay == "gDNA" else (0.05,)
    bp, band_rows, h = calibrate_ladder(data[:, descs.index("Ladder")], bands, thresholds)
    lanes = {w: (d, data[:, j]) for j, (w, d) in enumerate(zip(wells, descs)) if d != "Ladder"}
    return {"bp": bp, "lanes": lanes, "band_rows": band_rows, "threshold": h, "bands": bands}


def lane(data_dir, prep, filename, well, desc=None, assay="gDNA"):
    """(bp, fluorescence) for one lane; checks the lane description if given."""
    run = load_run(str(find_file(data_dir, prep, filename)), assay)
    d, y = run["lanes"][well]
    if desc is not None and d != desc:
        raise ValueError(f"{filename} {well}: description '{d}' != sample sheet '{desc}'")
    return run["bp"], y


def window_sum(bp, fu, lo, hi):
    """Summed fluorescence of rows with lo <= size < hi (negative values clipped)."""
    m = (bp >= lo) & (bp < hi)
    return float(np.clip(fu[m], 0, None).sum())


def avg_size_kb(bp, fu, lo=METRIC_LO, hi=METRIC_HI):
    """Fluorescence-weighted mean fragment size (kb) over lo < size < hi."""
    m = (bp > lo) & (bp < hi)
    return float(np.average(bp[m], weights=np.clip(fu[m], 0, None))) / 1000


def pct_ge(bp, fu, cutoff=HMW_CUTOFF, lo=METRIC_LO, hi=METRIC_HI):
    """Percentage of lo < size < hi signal at or above `cutoff`."""
    m = (bp > lo) & (bp < hi)
    w = np.clip(fu[m], 0, None)
    return 100 * float(w[bp[m] >= cutoff].sum() / w.sum())


def normalize_trace(bp, fu, above=METRIC_LO):
    """Trace scaled to its highest point above `above` bp."""
    return fu / fu[bp > above].max()


def bp_to_row(bp_per_row, size):
    """Row position of a size on a calibrated run (inverse of the calibration)."""
    return float(np.interp(np.log10(size), np.log10(bp_per_row), np.arange(len(bp_per_row))))


def on_reference_grid(bp_ref, bp, fu):
    """Resample a lane onto a reference run's rows by matching calibrated size."""
    rows_ref = np.arange(len(bp_ref))
    pos = np.interp(np.log10(bp), np.log10(bp_ref), rows_ref)
    return np.interp(rows_ref, pos, fu)


# --------------------------------------------------------------------------
# NanoDrop
# --------------------------------------------------------------------------
@lru_cache(maxsize=None)
def read_nanodrop(path):
    """Parse the long-format NanoDrop export ('dsDNA <date> <time>.tsv').

    One block per measurement: a timestamp line, key/value pairs (A260/A280,
    A260/A230, A260, A280) and wavelength/absorbance pairs. The wide
    '_table.tsv' export is column-misaligned and is not used.
    Returns {timestamp: record}.
    """
    recs, cur = {}, None
    for ln in open(path, encoding="utf-8", errors="replace").read().splitlines():
        if re.match(r"^\d{1,2}/\d{1,2}/\d{4} ", ln):
            cur = {"wl": [], "ab": []}
            recs[ln.strip()] = cur
            continue
        if cur is None:
            continue
        p = ln.split("\t")
        if len(p) != 2:
            continue
        k, v = p[0].strip(), p[1].strip()
        try:
            if k in ("A260", "A280", "A260/A280", "A260/A230"):
                cur[k] = float(v)
            else:
                x = float(k)
                if 200 <= x <= 900:
                    cur["wl"].append(x)
                    cur["ab"].append(float(v))
        except ValueError:
            pass
    for r in recs.values():
        r["wl"], r["ab"] = np.array(r["wl"]), np.array(r["ab"])
    return recs


def spectrum_normalized(rec, lo=220, hi=320):
    """Unsmoothed spectrum between lo and hi nm divided by its value at 260 nm."""
    w, a = rec["wl"], rec["ab"]
    m = (w >= lo) & (w <= hi)
    w, a = w[m], a[m]
    return w, a / a[np.argmin(np.abs(w - 260))]


def smooth_edge(a, k=9):
    """k-point moving average with edge-replicate padding (display only)."""
    p = k // 2
    return np.convolve(np.pad(a, p, mode="edge"), np.ones(k) / k, mode="valid")


# --------------------------------------------------------------------------
# Qubit
# --------------------------------------------------------------------------
@lru_cache(maxsize=None)
def read_qubit_export(path):
    """Parse a Qubit 4 'QubitData_*.csv' export into {Test Date: reading}.

    Exports carry no sample names; the sample sheet maps each well to its
    reading by acquisition time. 'Out of range' readings are classified by
    comparing the sample RFU with the standards: below the top standard means
    below the assay's quantifiable range (reported as at or below the
    detection limit), above it means over range.
    """
    out = {}
    for r in csv.DictReader(open(path, encoding="utf-8-sig")):
        rfu_col = "Far Red RFU" if r["Emission"] == "Far Red" else "Green RFU"
        rfu = float(r[rfu_col])
        stds = [float(r[c]) for c in ("Std 1 RFU", "Std 2 RFU", "Std 3 RFU") if r[c]]
        val = r["Original sample conc."].strip()
        if val == "Out of range":
            status = "below_range" if rfu < max(stds) else "above_range"
            conc = np.nan
        else:
            status, conc = "ok", float(val)
        out[r["Test Date"].strip()] = dict(assay=r["Assay Name"], conc=conc, status=status, rfu=rfu,
                                           std1_rfu=stds[0], units=r["Original sample conc. units"])
    return out


@lru_cache(maxsize=None)
def read_qubit_notion(path):
    """Parse a '<prep>_Qubit_from-Notion.csv' transcription into {Sample: reading}."""
    out = {}
    for r in csv.DictReader(open(path, encoding="utf-8-sig")):
        out[r["Sample"].strip()] = dict(assay=r["Assay"], conc=float(r["Qubit dsDNA (ng/uL)"]),
                                        status="ok", rfu=np.nan, std1_rfu=np.nan, units="ng/uL")
    return out


def qubit(data_dir, prep, filename, key):
    path = str(find_file(data_dir, prep, filename))
    table = read_qubit_notion(path) if "from-Notion" in filename else read_qubit_export(path)
    return table[key]


# --------------------------------------------------------------------------
# Summaries
# --------------------------------------------------------------------------
def describe(values):
    """n, mean, sample SD, min and max of the finite values."""
    v = np.asarray([x for x in values if np.isfinite(x)], float)
    sd = float(v.std(ddof=1)) if len(v) > 1 else np.nan
    return dict(n=len(v), mean=float(v.mean()), sd=sd, min=float(v.min()), max=float(v.max()))


def members(sets, name, group=None):
    t = sets[sets["set"] == name]
    if group is not None:
        t = t[t["group"] == group]
    return t.sort_values("order")
