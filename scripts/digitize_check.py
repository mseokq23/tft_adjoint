"""Draw the merged measured curves (data/meas.csv) back onto the paper's figure panels
(data/digitize/check_<panel>.png) for a visual check of the digitisation."""
import collections
import csv
import os

import numpy as np
from PIL import Image, ImageDraw

from digitize_paper import PANELS, SOURCES

COL = {("2W", 1.0, 1.0): (0, 0, 0), ("2W", 1.0, 0.05): (255, 0, 255), ("2W", 1.0, 0.1): (0, 160, 0),
       ("2W", 1.0, 0.2): (160, 90, 0), ("2W", 1.0, 0.5): (90, 0, 160), ("2W", 2.0, 1.0): (0, 0, 0),
       ("2W", 2.0, 0.05): (255, 0, 255), ("4W", 1.0, 1.0): (0, 120, 120), ("4W", 1.0, 0.05): (255, 140, 0)}


def main(meas="data/meas.csv", figdir="data/figs", od="data/digitize"):
    mm = collections.defaultdict(dict)
    for r in csv.DictReader(open(meas)):
        mm[(r["set"], float(r["L_um"]), float(r["vd"]))][float(r["vg"])] = float(r["id"])
    imgs = {}
    for panel, P in PANELS.items():
        if P["img"] not in imgs:
            imgs[P["img"]] = Image.open(os.path.join(figdir, P["img"])).convert("RGB")
        im = imgs[P["img"]].copy()
        dr = ImageDraw.Draw(im)
        ax = np.polyfit(P["gx"][1], P["gx"][0], 1)
        ay = np.polyfit(P["gy"][1], P["gy"][0], 1)
        for key, panels in SOURCES.items():
            if panel not in panels or key not in mm:
                continue
            for v, i in mm[key].items():
                y = np.log10(i) if P["ykind"] == "log" else i * 1e3
                x, yy = np.polyval(ax, v), np.polyval(ay, y)
                dr.rectangle((x - 3, yy - 3, x + 3, yy + 3), fill=COL[key])
        x0, x1 = np.polyval(ax, -1.1), np.polyval(ax, 2.1)
        ylo, yhi = (np.polyval(ay, -13), np.polyval(ay, -3)) if P["ykind"] == "log" else (np.polyval(ay, -0.02), np.polyval(ay, 0.42))
        im.crop((int(x0), int(min(ylo, yhi)), int(x1), int(max(ylo, yhi)))).save(os.path.join(od, f"check_{panel}.png"))


if __name__ == "__main__":
    main()
