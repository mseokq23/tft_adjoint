"""Digitise the measured transfer curves of arXiv:2609.12184 from the embedded bitmaps.

The figures are raster images. Measured data are hollow markers drawn with alpha
~0.75 (series colour blended with white), TCAD curves are opaque lines of the same
hue on top. Marker centres are recovered column by column from the vertical pixel
runs of the marker colour (digitize_raster.run_centers): top/bottom edge pairs,
full side edges, and side edges cut by an occluder completed from their free end.
Centres are binned on a 0.05 V grid (median).

Panels used
  Fig. 5(b) log/linear : 2 %-W, L = 1 um, VD = 1 V (unique maroon squares; cleanest)
  Fig. 6(a) log/linear : 2 %-W, L = 1 um, VD = 0.05/0.1/0.2/0.5/1 V (one colour per VD)
  Fig. 6(b) log/linear : 2 %-W, L = 1 and 2 um, VD = 0.05 and 1 V (two VDs per colour)
  Fig. 7(a) log/linear, 7(c) log : 2 %-W and 4 %-W, L = 1 um (two VDs per colour in 7a-log)

Where two VD curves share a colour they are separated with a reference: the
VD = 1 V curve of the same device read from a single-curve panel (5b for 2 %-W,
7c for 4 %-W); points within REF_TOL of it are VD = 1 V, points clearly below are
VD = 0.05 V. For L = 2 um (no reference) the upper/lower cluster in each VG bin is used.

Outputs: data/meas_panels.csv (per panel), data/meas.csv (merged, consensus median of
panels; linear panels only above 2e-5 A), data/digitize/consistency.md, overlays.

Usage:  python scripts/digitize_paper.py
"""
from __future__ import annotations

import argparse
import collections
import csv
import os

import numpy as np

from digitize_raster import bin_curve, color_mask, find_frame, load, marker_points, overlay, run_centers

ALPHA = 0.75
TAB = {"blue": (31, 119, 180), "green": (44, 160, 44), "orange": (255, 127, 14), "purple": (148, 103, 189),
       "red": (214, 39, 40), "cblue": (0, 114, 178)}
LIN_MIN_A = 2e-5            # linear panels are used only above this current (1 px < ~4 %)
REF_TOL = 0.12              # dec: a same-colour point within this of the VD = 1 V reference belongs to it
GRID = np.round(np.arange(-1.0, 2.0001, 0.05), 3)


def light(rgb, a=ALPHA):
    return tuple(a * np.asarray(rgb, float) + (1 - a) * 255.0)


def S(set_, L, vd, col, shape, size, thick):
    return dict(set=set_, L=L, vd=vd, col=col, shape=shape, size=size, thick=thick)


PANELS = {
    "5b_log": dict(img="p5_img0.png", box=(1340, 30, 2080, 780), gx=([1386, 1601.5, 1817.5, 2033], [-1, 0, 1, 2]),
                   gy=([207.5, 341.5, 475.5, 608.5], [-5, -7, -9, -11]), ykind="log",
                   series=[S("2W", 1.0, 1.0, "maroon", "square", 25, 4.5)]),
    "5b_lin": dict(img="p5_img0.png", box=(2090, 40, 2900, 780), gx=([2206.5, 2420.5, 2634.5, 2848.5], [-1, 0, 1, 2]),
                   gy=([107.5, 184, 260.5, 336.5, 413, 489.5, 565.5, 642], [0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05]),
                   ykind="mA", series=[S("2W", 1.0, 1.0, "maroon", "square", 25, 4.5)]),
    "6a_log": dict(img="p5_img1.png", box=(200, 200, 1120, 880), gx=([264.5, 534, 804, 1073.5], [-1, 0, 1, 2]),
                   gy=([311, 423, 535, 647, 759], [-4, -6, -8, -10, -12]), ykind="log",
                   series=[S("2W", 1.0, v, c, "square", 17, 3.5) for v, c in
                           ((0.05, "blue"), (0.1, "green"), (0.2, "orange"), (0.5, "purple"), (1.0, "red"))]),
    "6a_lin": dict(img="p5_img1.png", box=(1250, 250, 2120, 880), gx=([1309.5, 1568.5, 1827.5, 2087.5], [-1, 0, 1, 2]),
                   gy=([310.5, 437, 564.5, 691.5], [0.4, 0.3, 0.2, 0.1]), ykind="mA",
                   series=[S("2W", 1.0, v, c, "circle", 16, 3) for v, c in
                           ((0.05, "blue"), (0.1, "green"), (0.2, "orange"), (0.5, "purple"), (1.0, "red"))]),
    "6b_log": dict(img="p5_img1.png", box=(2330, 220, 3230, 900), gx=([2415.7, 2671.5, 2925.5, 3180.5], [-1, 0, 1, 2]),
                   gy=([326, 441.5, 556.5, 672.5, 787.5], [-4, -6, -8, -10, -12]), ykind="log",
                   series=[S("2W", 1.0, "split", "cblue", "square", 17, 3.5), S("2W", 2.0, "split", "red", "circle", 16, 3)]),
    "6b_lin": dict(img="p5_img1.png", box=(3330, 230, 4207, 900), gx=([3403.5, 3656.5, 3910.5, 4163.5], [-1, 0, 1, 2]),
                   gy=([295, 430.5, 566.5, 701.5], [0.4, 0.3, 0.2, 0.1]), ykind="mA",
                   series=[S("2W", 1.0, "split", "cblue", "square", 16, 3.5), S("2W", 2.0, "split", "red", "circle", 15, 3)]),
    "7a_log": dict(img="p6_img0.png", box=(240, 60, 1160, 720), gx=([292.5, 565.5, 839.5, 1113], [-1, 0, 1, 2]),
                   gy=([157.5, 271, 385, 499, 613], [-4, -6, -8, -10, -12]), ykind="log",
                   series=[S("2W", 1.0, "split", "blue", "square", 22, 4), S("4W", 1.0, "split", "red", "diamond", 30, 5)]),
    "7a_lin": dict(img="p6_img0.png", box=(1300, 60, 2210, 700), gx=([1358.5, 1626.5, 1893.5, 2161.5], [-1, 0, 1, 2]),
                   gy=([122.5, 254, 384.5, 515.5], [0.4, 0.3, 0.2, 0.1]), ykind="mA",
                   series=[S("2W", 1.0, "split", "blue", "square", 22, 4), S("4W", 1.0, "split", "red", "diamond", 30, 5)]),
    "7c_log": dict(img="p6_img0.png", box=(1290, 960, 2240, 1640), gx=([1343.5, 1624, 1905.5, 2186], [-1, 0, 1, 2]),
                   gy=([1061.5, 1179, 1295.5, 1412.5, 1529], [-4, -6, -8, -10, -12]), ykind="log",
                   series=[S("2W", 1.0, 1.0, "blue", "square", 22, 4), S("4W", 1.0, 1.0, "red", "diamond", 30, 5)]),
}

# which panels show which curve; the merged curve is the median of the available panels
SOURCES = {
    ("2W", 1.0, 1.0): ["5b_log", "5b_lin", "6a_log", "6a_lin", "6b_log", "6b_lin", "7a_log", "7a_lin", "7c_log"],
    ("2W", 1.0, 0.05): ["6a_log", "6b_log", "7a_log", "6a_lin", "6b_lin"],
    ("2W", 1.0, 0.1): ["6a_log", "6a_lin"],
    ("2W", 1.0, 0.2): ["6a_log", "6a_lin"],
    ("2W", 1.0, 0.5): ["6a_log", "6a_lin"],
    ("2W", 2.0, 1.0): ["6b_log", "6b_lin"],
    ("2W", 2.0, 0.05): ["6b_log", "6b_lin"],
    ("4W", 1.0, 1.0): ["7c_log", "7a_log", "7a_lin"],
    ("4W", 1.0, 0.05): ["7a_log"],
}
REFERENCE = {("2W", 1.0): ["5b_log", "5b_lin"], ("4W", 1.0): ["7c_log"]}    # single-curve panels for splitting


def fit_axis(px, val):
    a, b = np.polyfit(np.asarray(px, float), np.asarray(val, float), 1)
    return (lambda p: a * np.asarray(p, float) + b), float(np.max(np.abs(np.asarray(val) - (a * np.asarray(px) + b))))


def colour(name):
    if name == "maroon":
        return (130.0, 6.0, 6.0), 50
    return light(TAB[name]), 34


def occluders(img, own_rgb):
    sat = (img.max(axis=2) - img.min(axis=2)) > 40
    darkp = img.max(axis=2) < 120
    own = np.sqrt(((img - np.asarray(own_rgb)) ** 2).sum(axis=2)) < 40
    return (sat | darkp) & ~own


def to_log(y, kind):
    if kind == "log":
        return np.asarray(y, float)
    v = np.asarray(y, float) * 1e-3
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(v > 0, np.log10(np.where(v > 0, v, 1)), np.nan)


def extract_all(figdir, od):
    """Column marker centres of every series in every panel, then grid-anchored readout.
    Returns {((set, L, vd), panel): {vg: log10 id}}."""
    from scipy import ndimage
    from digitize_raster import grid_points_x
    cache, items = {}, []
    for key, P in PANELS.items():
        if P["img"] not in cache:
            cache[P["img"]] = load(os.path.join(figdir, P["img"]))
        img = cache[P["img"]]
        fr = find_frame(img, P["box"])
        a_x, b_x = np.polyfit(P["gx"][1], P["gx"][0], 1)                  # VG -> pixel
        a_y, b_y = np.polyfit(P["gy"][1], P["gy"][0], 1)                  # value -> pixel
        for s in P["series"]:
            rgb, tol = colour(s["col"])
            m = color_mask(img, fr, rgb, tol=tol)
            m = ndimage.binary_opening(m, structure=np.ones((2, 2)))       # drop 1-px anti-aliasing seams
            c = run_centers(m, occluders(img, rgb), s["size"], s["thick"], s["shape"])
            items.append(dict(panel=key, P=P, s=s, c=c, ax=(a_x, b_x), ay=(a_y, b_y), img=P["img"], frame=fr))
    per_panel, draw = {}, collections.defaultdict(list)
    cols = {}

    def readout(it, cen, vd):
        P, s = it["P"], it["s"]
        cols[((s["set"], s["L"], float(vd)), it["panel"])] = (cen, it["ax"], it["ay"], P["ykind"], s["size"])
        a_x, b_x = it["ax"]
        a_y, b_y = it["ay"]
        pts = grid_points_x(cen, lambda v: a_x * v + b_x, GRID, s["size"])
        if not pts:
            return
        vg = np.array(sorted(pts))
        ypx = np.array([pts[v] for v in vg])
        val = (ypx - b_y) / a_y
        ly = to_log(val, P["ykind"])
        ok = np.isfinite(ly)
        if P["ykind"] == "mA":
            ok &= ly >= np.log10(LIN_MIN_A)
        if ok.sum():
            per_panel[((s["set"], s["L"], float(vd)), it["panel"])] = dict(zip(np.round(vg[ok], 3), ly[ok]))
            dot = (0, 0, 0) if vd == 1.0 else (255, 0, 255)
            draw[(it["img"], it["panel"])].append((dot, list(zip(a_x * vg[ok] + b_x, ypx[ok]))))

    for it in items:                                                      # single-curve series first
        if it["s"]["vd"] != "split" and it["c"].size:
            readout(it, it["c"], it["s"]["vd"])
    refs = {}
    for (set_, L), panels in REFERENCE.items():
        vals = collections.defaultdict(list)
        for p in panels:
            for v, y in per_panel.get(((set_, L, 1.0), p), {}).items():
                vals[v].append(y)
        refs[(set_, L)] = {v: float(np.median(y)) for v, y in vals.items()}
    for it in items:                                                      # two curves of one colour
        if it["s"]["vd"] != "split" or not it["c"].size:
            continue
        c = it["c"]
        a_x, b_x = it["ax"]
        a_y, b_y = it["ay"]
        P, s = it["P"], it["s"]
        vg_c = (c[:, 0] - b_x) / a_x
        val = (c[:, 1] - b_y) / a_y
        lc = to_log(val, P["ykind"])
        ref = refs.get((s["set"], s["L"]))
        if ref:
            rv = np.array(sorted(ref))
            r = np.interp(vg_c, rv, [ref[v] for v in rv], left=np.nan, right=np.nan)
            hi = np.isfinite(r) & (np.abs(lc - r) <= REF_TOL)
            lo = np.isfinite(r) & (lc < r - REF_TOL)
        else:                                                             # split at the largest gap per 0.05 V bin
            hi = np.zeros(len(c), bool)
            lo = np.zeros(len(c), bool)
            idx = np.round((vg_c + 1.0) / 0.05).astype(int)
            for k in np.unique(idx):
                sel = np.where((idx == k) & np.isfinite(lc))[0]
                if sel.size < 2:
                    continue
                ys = np.sort(lc[sel])
                g = np.diff(ys)
                j = int(np.argmax(g))
                if g[j] > 0.25:
                    hi[sel[lc[sel] > ys[j]]] = True
                    lo[sel[lc[sel] <= ys[j]]] = True
        readout(it, c[hi], 1.0)
        readout(it, c[lo], 0.05)
    if os.environ.get("DIGI_REFINE"):
        per_panel = refine(per_panel, cols)
    draw = collections.defaultdict(list)
    for (key, panel), d in per_panel.items():
        cen, (a_x, b_x), (a_y, b_y), kind, size = cols[(key, panel)]
        img = next(it["img"] for it in items if it["panel"] == panel)
        vv = np.array(sorted(d))
        val = np.array([d[v] for v in vv]) if kind == "log" else 10 ** np.array([d[v] for v in vv]) * 1e3
        draw[(img, panel)].append(((0, 0, 0) if key[2] == 1.0 else (255, 0, 255), list(zip(a_x * vv + b_x, a_y * val + b_y))))
    for (imgname, panel), d in draw.items():
        fr = next(it["frame"] for it in items if it["panel"] == panel)
        overlay(cache[imgname], d, os.path.join(od, f"overlay_{panel}.png"),
                box=(int(fr[0]) - 25, int(fr[2]) - 25, int(fr[1]) + 25, int(fr[3]) + 25))
    for (key, panel), d in sorted(per_panel.items()):
        print(f"[{panel}] {key}: {len(d)} grid points")
    return per_panel


def consensus(per_panel, key):
    vals = collections.defaultdict(list)
    for (k, p), d in per_panel.items():
        if k == key:
            for v, y in d.items():
                vals[v].append(y)
    if not vals:
        return None
    vv = np.array(sorted(vals))
    yy = np.array([np.median(vals[v]) for v in vv])
    return vv, yy


def refine(per_panel, cols, iters=3, half_frac=0.45, gap=2.5, accept=0.3):
    """Predictor-guided readout: at each grid VG, take the cluster of column centres that is
    closest to the multi-panel consensus (within `accept` dec). Resolves which of the
    overlapping neighbour markers is the one at this VG, and fills gaps where the plain
    readout was ambiguous. Iterated a few times."""
    for _ in range(iters):
        new = {}
        for (key, panel), (cen, (a_x, b_x), (a_y, b_y), kind, size) in cols.items():
            cs = consensus(per_panel, key)
            if cs is None or cen.size == 0:
                continue
            vv, yy = cs
            spacing = abs(a_x) * 0.05
            half = half_frac * spacing
            xs, ys = cen[:, 0], cen[:, 1]
            d = {}
            for v in GRID:
                if v < vv[0] - 0.051 or v > vv[-1] + 0.051:
                    continue
                pred = float(np.interp(v, vv, yy))
                sel = np.abs(xs - (a_x * v + b_x)) <= half
                if sel.sum() < 2:
                    continue
                py = np.sort(ys[sel])
                clusters = np.split(py, np.where(np.diff(py) > gap)[0] + 1)
                best, bestd = None, 1e9
                for cl in clusters:
                    if len(cl) < 2:
                        continue
                    yv = float(np.median(cl))
                    val = (yv - b_y) / a_y
                    lv = float(to_log(val, kind)) if kind != "log" else val
                    if not np.isfinite(lv):
                        continue
                    if kind == "mA" and lv < np.log10(LIN_MIN_A):
                        continue
                    if abs(lv - pred) < bestd:
                        best, bestd = lv, abs(lv - pred)
                if best is not None and bestd <= accept:
                    d[round(float(v), 3)] = best
            if d:
                new[(key, panel)] = d
        per_panel = new
    return per_panel


def monotone_clean(cur, tol=0.05, floor=-10.3):
    """Drop points that break the monotonic increase of log10 ID with VG by more than `tol`
    (worst offender first). Points below `floor` (measurement noise floor) are not tested."""
    cur = sorted(cur)
    dropped = []
    while True:
        worst, wi = 0.0, None
        for i in range(len(cur)):
            v, y = cur[i][0], cur[i][1]
            if y < floor:
                continue
            viol = 0.0
            if i + 1 < len(cur) and cur[i + 1][1] >= floor:
                viol = max(viol, y - cur[i + 1][1] - tol)          # higher than the next point
            if i > 0 and cur[i - 1][1] >= floor:
                viol = max(viol, cur[i - 1][1] - y - tol)          # lower than the previous point
            if viol > worst:
                worst, wi = viol, i
        if wi is None:
            return cur, dropped
        # remove the point of the violating pair that deviates most from its other neighbours
        cand = [wi] + [j for j in (wi - 1, wi + 1) if 0 <= j < len(cur)]

        def dev(j):
            nb = [cur[k][1] for k in (j - 1, j + 1) if 0 <= k < len(cur)]
            return abs(cur[j][1] - np.mean(nb)) if nb else 0.0
        j = max(cand, key=dev)
        dropped.append(cur[j][0])
        cur.pop(j)


def merge(per_panel, od, out):
    lines = ["# Digitisation consistency: the same measured curve read from different panels", "",
             "Differences above the 10 pA floor (log10 ID > -11); linear panels only above 2e-5 A.", "",
             "| curve | panel A | panel B | common bins | RMS Δ (dec) | max |Δ| (dec) |", "|---|---|---|---|---|---|"]
    rows, allspread, removed = [], [], {}
    for key, panels in SOURCES.items():
        have = [p for p in panels if (key, p) in per_panel and per_panel[(key, p)]]
        for i in range(len(have)):
            for j in range(i + 1, len(have)):
                A, B = per_panel[(key, have[i])], per_panel[(key, have[j])]
                com = [v for v in sorted(set(A) & set(B)) if A[v] > -11 and B[v] > -11]
                if len(com) < 5:
                    continue
                d = np.array([A[v] - B[v] for v in com])
                allspread.append(d)
                lines.append(f"| {key[0]} L={key[1]:g} VD={key[2]:g} | {have[i]} | {have[j]} | {len(com)} | "
                             f"{np.sqrt(np.mean(d ** 2)):.3f} | {np.max(np.abs(d)):.3f} |")
        cur = []
        for v in GRID:
            vals = [per_panel[(key, p)][v] for p in have if v in per_panel[(key, p)]]
            if vals:
                med = float(np.median(vals))
                spread = float(np.max(vals) - np.min(vals)) if len(vals) > 1 else float("nan")
                cur.append([float(v), med, len(vals), spread])
        cur, dropped = monotone_clean(cur)
        removed[key] = dropped
        rows += [(key[0], key[1], key[2], v, 10 ** y, n, sp) for v, y, n, sp in cur]
    a = np.concatenate(allspread) if allspread else np.array([np.nan])
    rms = float(np.sqrt(np.nanmean(a ** 2)))
    lines += ["", f"All pairs: RMS {rms:.3f} dec → single-reading σ ≈ {rms / np.sqrt(2):.3f} dec", "",
              "Points removed by the monotonicity check (transfer curves must increase with VG above the 10 pA floor):", ""]
    lines += [f"- {k[0]} L={k[1]:g} VD={k[2]:g}: " + (", ".join(f"VG {v:+.2f}" for v in d) if d else "none")
              for k, d in removed.items()]
    with open(os.path.join(od, "consistency.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(os.path.join(out, "meas.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["set", "L_um", "vd", "vg", "id", "n_panels", "spread_dec"])
        w.writerows(rows)
    print("\n".join(lines))
    print(f"wrote {len(rows)} points to {out}/meas.csv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figdir", default="data/figs")
    ap.add_argument("--out", default="data")
    args = ap.parse_args()
    od = os.path.join(args.out, "digitize")
    os.makedirs(od, exist_ok=True)
    per_panel = extract_all(args.figdir, od)
    with open(os.path.join(args.out, "meas_panels.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["panel", "set", "L_um", "vd", "vg", "log10_id"])
        for ((set_, L, vd), panel), d in sorted(per_panel.items()):
            w.writerows([(panel, set_, L, vd, v, y) for v, y in sorted(d.items())])
    merge(per_panel, od, args.out)


if __name__ == "__main__":
    main()
