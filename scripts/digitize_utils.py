"""Helpers to turn the measured transfer curves of arXiv:2609.12184 into data/meas.csv.

Two routes:
1. Vector extraction (preferred): if the figure is vector graphics in the PDF,
   PyMuPDF returns the plotted polylines exactly.
       python scripts/digitize_utils.py list  paper.pdf --page 5
       python scripts/digitize_utils.py dump  paper.pdf --page 5 --out data/page5_paths.json
   Then calibrate the axes with two known ticks per axis (PDF coordinates ->
   data values; the current axis is logarithmic) and convert the chosen path:
       python scripts/digitize_utils.py convert data/page5_paths.json --path 12 \
           --x0 PX0 VG0 --x1 PX1 VG1 --y0 PY0 I0 --y1 PY1 I1 --logy \
           --set step5 --L 1 --vd 1.0 --append data/meas.csv
2. WebPlotDigitizer (raster figures): export "x,y" CSV per curve and run
       python scripts/digitize_utils.py wpd curve.csv --set step5 --L 1 --vd 1.0 --append data/meas.csv

Repeat the digitization 3 times and keep the spread: it is used as the
measurement noise (sigma of log10 ID) in the calibration/identifiability steps.

Requires: pip install pymupdf  (route 1 only)
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os


def cmd_list(args):
    import fitz  # PyMuPDF
    doc = fitz.open(args.pdf)
    page = doc[args.page - 1]
    drawings = page.get_drawings()
    print(f"page {args.page}: {len(drawings)} vector drawings")
    for k, d in enumerate(drawings):
        n = sum(len(it) for it in d["items"])
        r = d["rect"]
        print(f"[{k:4d}] items={len(d['items']):4d} pts~{n:5d} rect=({r.x0:.1f},{r.y0:.1f},{r.x1:.1f},{r.y1:.1f}) "
              f"color={d.get('color')} width={d.get('width')}")


def cmd_dump(args):
    import fitz
    doc = fitz.open(args.pdf)
    page = doc[args.page - 1]
    paths = []
    for d in page.get_drawings():
        pts = []
        for it in d["items"]:
            if it[0] == "l":
                pts += [(it[1].x, it[1].y), (it[2].x, it[2].y)]
            elif it[0] == "c":
                pts += [(it[1].x, it[1].y), (it[4].x, it[4].y)]
        paths.append({"rect": list(d["rect"]), "color": d.get("color"), "width": d.get("width"), "points": pts})
    with open(args.out, "w") as f:
        json.dump(paths, f)
    print(f"wrote {len(paths)} paths to {args.out}")


def _axis_map(p0, v0, p1, v1, log):
    if log:
        a = (math.log10(v1) - math.log10(v0)) / (p1 - p0)
        return lambda p: 10 ** (math.log10(v0) + a * (p - p0))
    a = (v1 - v0) / (p1 - p0)
    return lambda p: v0 + a * (p - p0)


def _append(rows, path):
    new = not os.path.exists(path)
    with open(path, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["set", "L_um", "vd", "vg", "id"])
        w.writerows(rows)
    print(f"appended {len(rows)} points to {path}")


def cmd_convert(args):
    paths = json.load(open(args.paths))
    pts = paths[args.path]["points"]
    fx = _axis_map(args.x0[0], args.x0[1], args.x1[0], args.x1[1], False)
    fy = _axis_map(args.y0[0], args.y0[1], args.y1[0], args.y1[1], args.logy)
    seen, rows = set(), []
    for px, py in pts:
        key = (round(px, 2), round(py, 2))
        if key in seen:
            continue
        seen.add(key)
        rows.append([args.set, args.L, args.vd, round(fx(px), 6), fy(py)])
    rows.sort(key=lambda r: r[3])
    _append(rows, args.append)


def cmd_wpd(args):
    rows = []
    with open(args.csv) as f:
        for line in f:
            parts = line.strip().split(",")
            if len(parts) < 2:
                continue
            try:
                x, y = float(parts[0]), float(parts[1])
            except ValueError:
                continue
            rows.append([args.set, args.L, args.vd, x, y])
    rows.sort(key=lambda r: r[3])
    _append(rows, args.append)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("list"); a.add_argument("pdf"); a.add_argument("--page", type=int, required=True)
    a = sub.add_parser("dump"); a.add_argument("pdf"); a.add_argument("--page", type=int, required=True); a.add_argument("--out", required=True)
    a = sub.add_parser("convert"); a.add_argument("paths"); a.add_argument("--path", type=int, required=True)
    for k in ("x0", "x1", "y0", "y1"):
        a.add_argument(f"--{k}", type=float, nargs=2, required=True, metavar=("PDF_COORD", "VALUE"))
    a.add_argument("--logy", action="store_true")
    for k in ("set", "L", "vd", "append"):
        a.add_argument(f"--{k}", required=True)
    a = sub.add_parser("wpd"); a.add_argument("csv")
    for k in ("set", "L", "vd", "append"):
        a.add_argument(f"--{k}", required=True)
    args = ap.parse_args()
    {"list": cmd_list, "dump": cmd_dump, "convert": cmd_convert, "wpd": cmd_wpd}[args.cmd](args)


if __name__ == "__main__":
    main()
