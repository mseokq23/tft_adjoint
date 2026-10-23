"""Raster digitiser for the measured transfer curves in the figures of arXiv:2609.12184.

The figures in the PDF are bitmaps (no vector paths), so the measured curves are
recovered from pixel colours:

* measured data are hollow markers drawn with alpha ~0.75, so their colour is the
  series colour blended with white (e.g. tab:blue (31,119,180) -> (87,153,199)),
  which separates them from the opaque TCAD lines of the same hue;
* the axes frame and the light-grey grid lines give the pixel <-> data calibration;
* for every pixel column, marker pixels within +-s/2 (s = marker size) are split
  into clusters along y (gap > s); the mean y of a cluster of ring pixels is the
  local marker-centre height. Columns are then binned in VG and the median taken.

Every extraction writes an overlay PNG so the result can be checked by eye, and
curves that appear in several panels are compared to estimate the digitisation
error (used as sigma in the calibration).
"""
from __future__ import annotations

import numpy as np
from PIL import Image


def load(path):
    return np.asarray(Image.open(path).convert("RGB")).astype(float)


def dark(img, thr=80):
    return img.max(axis=2) < thr


def find_frame(img, box, min_frac=0.6):
    """Axes frame inside an approximate box (x0, y0, x1, y1): the outermost dark
    vertical/horizontal lines that span most of the box."""
    x0, y0, x1, y1 = box
    d = dark(img[y0:y1, x0:x1])
    colfrac = d.mean(axis=0)
    rowfrac = d.mean(axis=1)
    cols = np.where(colfrac > min_frac)[0]
    rows = np.where(rowfrac > min_frac)[0]
    if cols.size < 2 or rows.size < 2:
        raise RuntimeError(f"frame not found in {box}: cols {cols.size}, rows {rows.size}")
    # left/right: the first and last runs of dark columns; take their inner edges
    def runs(idx):
        out, start = [], idx[0]
        for a, b in zip(idx[:-1], idx[1:]):
            if b != a + 1:
                out.append((start, a))
                start = b
        out.append((start, idx[-1]))
        return out
    cr, rr = runs(cols), runs(rows)
    left, right = cr[0][1] + 1, cr[-1][0] - 1
    top, bottom = rr[0][1] + 1, rr[-1][0] - 1
    return x0 + left, x0 + right, y0 + top, y0 + bottom


def grid_positions(img, frame, axis, lo=200, hi=248, min_frac=0.5):
    """Light-grey grid lines inside the frame (positions in image pixels)."""
    xl, xr, yt, yb = frame
    sub = img[yt + 3:yb - 3, xl + 3:xr - 3]
    grey = (sub.max(axis=2) - sub.min(axis=2) < 12) & (sub.mean(axis=2) > lo) & (sub.mean(axis=2) < hi)
    frac = grey.mean(axis=1) if axis == "y" else grey.mean(axis=0)
    idx = np.where(frac > min_frac)[0]
    if idx.size == 0:
        return np.array([])
    groups, cur = [], [idx[0]]
    for a in idx[1:]:
        if a - cur[-1] <= 2:
            cur.append(a)
        else:
            groups.append(cur)
            cur = [a]
    groups.append(cur)
    off = (yt + 3) if axis == "y" else (xl + 3)
    return np.array([off + np.mean(g) for g in groups])


class Axis:
    """Linear map pixel -> value (value may be log10 of the plotted quantity)."""

    def __init__(self, p0, v0, p1, v1):
        self.a = (v1 - v0) / (p1 - p0)
        self.b = v0 - self.a * p0

    def __call__(self, p):
        return self.a * np.asarray(p, float) + self.b

    def inv(self, v):
        return (np.asarray(v, float) - self.b) / self.a


def color_mask(img, frame, rgb, tol=28, inset=8, exclude=None, exclude_tol=22):
    """Pixels within `tol` (Euclidean RGB) of `rgb` inside the frame; optionally drop
    pixels closer to any colour in `exclude`."""
    xl, xr, yt, yb = frame
    sub = img[yt + inset:yb - inset, xl + inset:xr - inset]
    dist = np.sqrt(((sub - np.asarray(rgb, float)) ** 2).sum(axis=2))
    m = dist < tol
    if exclude:
        for c in exclude:
            m &= np.sqrt(((sub - np.asarray(c, float)) ** 2).sum(axis=2)) > exclude_tol
    full = np.zeros(img.shape[:2], bool)
    full[yt + inset:yb - inset, xl + inset:xr - inset] = m
    return full


def column_clusters(mask, size, min_px=4):
    """For every column x: clusters of marker pixels within the window x +- size/2,
    split along y where the gap exceeds `size`. Returns rows (x, y_mean, n_px)."""
    H, W = mask.shape
    ys_by_col = [np.where(mask[:, x])[0] for x in range(W)]
    half = size // 2
    out = []
    for x in range(W):
        lo, hi = max(0, x - half), min(W, x + half + 1)
        ys = np.concatenate([ys_by_col[k] for k in range(lo, hi)])
        if ys.size < min_px:
            continue
        ys.sort()
        cuts = np.where(np.diff(ys) > size)[0]
        for seg in np.split(ys, cuts + 1):
            if seg.size >= min_px:
                out.append((x, float(seg.mean()), int(seg.size), float(seg.min()), float(seg.max())))
    return np.array(out) if out else np.zeros((0, 5))


def bin_curve(xv, yv, edges, min_n=3):
    """Median of yv in bins of xv (returns bin centres with data)."""
    xc, yc, n = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (xv >= a) & (xv < b)
        if m.sum() >= min_n:
            xc.append(0.5 * (a + b))
            yc.append(float(np.median(yv[m])))
            n.append(int(m.sum()))
    return np.array(xc), np.array(yc), np.array(n)


def overlay(img, points_by_color, out_png, box=None):
    """Draw extracted points (pixel coords) on a copy of the image for visual checks."""
    from PIL import ImageDraw
    im = Image.fromarray(img.astype(np.uint8))
    dr = ImageDraw.Draw(im)
    for col, pts in points_by_color:
        for x, y in pts:
            dr.ellipse((x - 3, y - 3, x + 3, y + 3), fill=tuple(int(c) for c in col))
    if box:
        im = im.crop(box)
    im.save(out_png)


def ring_kernel(shape, size, thick):
    """Binary ring template (hollow square or circle) of outer size `size` px."""
    n = int(np.ceil(size)) + 2
    n += (n + 1) % 2                      # odd
    c = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n] - c
    if shape == "square":
        r = np.maximum(np.abs(xx), np.abs(yy))
    elif shape == "diamond":
        r = np.abs(xx) + np.abs(yy)
    else:
        r = np.hypot(xx, yy)
    R = size / 2.0
    return ((r <= R) & (r >= R - thick)).astype(float)


def detect_markers(mask, kernel, thr=0.45, nms=None):
    """Matched filter for hollow markers: peaks of the fraction of the template ring
    covered by marker pixels. Robust to partial occlusion (the peak stays at the
    centre). Returns arrays x, y, score (image pixel coordinates)."""
    from scipy import ndimage
    from scipy.signal import fftconvolve
    A = fftconvolve(mask.astype(float), kernel[::-1, ::-1], mode="same") / kernel.sum()
    size = nms or max(5, int(round(kernel.shape[0] * 0.55)))
    mx = ndimage.maximum_filter(A, size=size)
    peaks = (A >= mx - 1e-9) & (A >= thr)
    lab, n = ndimage.label(peaks)
    if n == 0:
        return np.zeros(0), np.zeros(0), np.zeros(0)
    cy, cx = np.array(ndimage.center_of_mass(peaks, lab, range(1, n + 1))).T
    sc = np.array(ndimage.maximum(A, lab, range(1, n + 1)))
    return cx, cy, sc


def _runs(col):
    """Start/end (inclusive) of True runs in a 1-D boolean array."""
    d = np.diff(np.concatenate([[0], col.astype(np.int8), [0]]))
    starts, ends = np.where(d == 1)[0], np.where(d == -1)[0] - 1
    return list(zip(starts, ends))


def run_centers(mask, occ, size, thick, shape="square", tol=3):
    """Marker centres from vertical pixel runs, column by column, with occlusion handling.

    In a column through a hollow marker of outer size s and edge thickness t one sees
    either two short runs (top/bottom edge, centres s - t apart) or, through a side
    edge, one run of length ~s. Pairs and full side-runs give the centre directly.
    A run cut short by an occluder (line or another series drawn on top) is completed
    from its free end (squares only). Unexplained fragments are ignored.
    Returns an array of (x, y_centre, kind) with kind 0 = pair, 1 = side, 2 = completed side.
    """
    H, W = mask.shape
    out = []
    gap = size - thick
    for x in range(W):
        runs = _runs(mask[:, x])
        if not runs:
            continue
        used = set()
        mids = [(a + b) / 2.0 for a, b in runs]
        lens = [b - a + 1 for a, b in runs]
        # full side edges
        for k, (a, b) in enumerate(runs):
            if abs(lens[k] - size) <= tol + 1:
                out.append((x, mids[k], 1))
                used.add(k)
        # top/bottom pairs (circles: any chord is symmetric about the centre row)
        gmin = gap - tol if shape == "square" else 0.3 * size
        for k in range(len(runs)):
            if k in used or lens[k] > thick + tol + (0 if shape == "square" else 4):
                continue
            for j in range(k + 1, len(runs)):
                if j in used or lens[j] > thick + tol + (0 if shape == "square" else 4):
                    continue
                if gmin <= (mids[j] - mids[k]) <= gap + tol:
                    out.append((x, 0.5 * (mids[j] + mids[k]), 0))
                    used.update((k, j))
                    break
        if shape != "square":
            continue
        # side edges truncated by an occluder at one end
        for k, (a, b) in enumerate(runs):
            if k in used or lens[k] < 0.4 * size or lens[k] > size + tol:
                continue
            below = occ[b + 1:min(H, b + 3), x].any()
            above = occ[max(0, a - 2):a, x].any()
            if below and not above:
                out.append((x, a + (size - 1) / 2.0, 2))
            elif above and not below:
                out.append((x, b - (size - 1) / 2.0, 2))
    return np.array(out) if out else np.zeros((0, 3))


def group_markers(centers, size, dy=1.6, dx=2):
    """Group column centres (x, y, kind) into individual markers: centres in nearby columns
    at the same height belong to one marker. Fragments of one marker split by an occluder
    are re-joined. Returns array (x_centre, y_centre, width_px, n_columns)."""
    if centers.size == 0:
        return np.zeros((0, 4))
    order = np.lexsort((centers[:, 1], centers[:, 0]))
    c = centers[order]
    n = len(c)
    parent = np.arange(n)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    xs = c[:, 0].astype(int)
    by_col = {}
    for i, x in enumerate(xs):
        by_col.setdefault(x, []).append(i)
    for i, x in enumerate(xs):
        for xx in range(x - dx, x):
            for j in by_col.get(xx, []):
                if abs(c[j, 1] - c[i, 1]) <= dy:
                    ri, rj = find(i), find(j)
                    if ri != rj:
                        parent[ri] = rj
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    g = []
    for idx in groups.values():
        x = c[idx, 0]
        y = c[idx, 1]
        g.append([float(np.median(y)), float(x.min()), float(x.max()), len(idx)])
    g.sort(key=lambda r: (r[1], r[0]))
    # re-join fragments of the same marker (same height, close in x, combined width <= size + 3)
    merged = []
    for yc, x0, x1, k in g:
        for m in merged:
            if abs(m[0] - yc) <= 2 and x0 - m[2] <= size * 0.6 and max(x1, m[2]) - min(x0, m[1]) <= size + 3:
                m[0] = (m[0] * m[3] + yc * k) / (m[3] + k)
                m[1], m[2], m[3] = min(m[1], x0), max(m[2], x1), m[3] + k
                break
        else:
            merged.append([yc, x0, x1, k])
    out = [((x0 + x1) / 2.0, yc, x1 - x0 + 1, k) for yc, x0, x1, k in merged]
    return np.array(out)


def marker_points(centers, size, x_of_v, v_grid, snap_px=None):
    """Marker data points from column centres: isolated markers give one point at the
    centre of their x-range; runs of overlapping markers at nearly equal height (flat
    parts of a curve) are fitted with a line and sampled at the grid columns inside them.
    Points are snapped to the measurement grid (`v_grid`, mapped to pixels by `x_of_v`).
    Returns {v: y_pixel}."""
    g = group_markers(centers, size)
    if g.size == 0:
        return {}
    xg = np.asarray([x_of_v(v) for v in v_grid])
    snap = snap_px or 0.3 * np.median(np.diff(xg))
    pts = {}
    order = np.lexsort((centers[:, 1], centers[:, 0]))
    c = centers[order]
    for xc, yc, w, k in g:
        if w < 0.5 * size:
            continue
        if w <= size + 3:
            j = int(np.argmin(np.abs(xg - xc)))
            if abs(xg[j] - xc) <= snap:
                pts.setdefault(v_grid[j], []).append(yc)
        else:
            x0, x1 = xc - w / 2, xc + w / 2
            sel = (c[:, 0] >= x0) & (c[:, 0] <= x1) & (np.abs(c[:, 1] - yc) <= 6)
            if sel.sum() < 4:
                continue
            a, b = np.polyfit(c[sel, 0], c[sel, 1], 1)
            for j, xv in enumerate(xg):
                if x0 + size / 2 - 2 <= xv <= x1 - size / 2 + 2:
                    pts.setdefault(v_grid[j], []).append(a * xv + b)
    return {v: float(np.median(ys)) for v, ys in pts.items()}


def grid_points(centers, x_of_v, v_grid, half=3.0, gap=2.5, min_cols=2):
    """Grid-anchored readout. The measurement VG grid is known, so the centre column of
    every marker is known; only its height is unknown. Column centres within +-`half`
    px of the marker column are clustered in y: marker v itself contributes most of
    them; overlapping neighbours v-1, v+1 add clusters below/above (steep parts), so
    with three clusters the middle one is taken, otherwise the largest. Returns {v: y_px}."""
    if centers.size == 0:
        return {}
    xs, ys = centers[:, 0], centers[:, 1]
    out = {}
    for v in v_grid:
        xv = x_of_v(v)
        sel = np.abs(xs - xv) <= half
        if sel.sum() < min_cols:
            continue
        yy = np.sort(ys[sel])
        clusters = np.split(yy, np.where(np.diff(yy) > gap)[0] + 1)
        pick = clusters[1] if len(clusters) == 3 else max(clusters, key=len)
        if len(pick) >= min_cols:
            out[float(v)] = float(np.median(pick))
    return out


def grid_points_x(centers, x_of_v, v_grid, size, gap=2.5, min_cols=3, pred=None, pred_tol_px=None):
    """Grid-anchored readout using the markers' x-extent: around the expected marker
    column x_v, column centres within +-size are clustered by height; for each height
    cluster the mean x of its columns is the x-centre of that marker (or of a flat run
    of markers). The cluster whose x-centre is closest to x_v (within 0.45 spacing) is
    the marker at this VG. `pred` (optional {v: y_px}) breaks near-ties."""
    if centers.size == 0:
        return {}
    xs, ys = centers[:, 0], centers[:, 1]
    xg = np.array([x_of_v(v) for v in v_grid])
    spacing = float(np.median(np.abs(np.diff(xg))))
    out = {}
    for v, xv in zip(v_grid, xg):
        sel = np.abs(xs - xv) <= size
        if sel.sum() < min_cols:
            continue
        sx, sy = xs[sel], ys[sel]
        o = np.argsort(sy)
        sx, sy = sx[o], sy[o]
        cuts = np.where(np.diff(sy) > gap)[0] + 1
        best = None
        for gx, gy in zip(np.split(sx, cuts), np.split(sy, cuts)):
            if gx.size < min_cols:
                continue
            # columns of this height cluster must cover x_v (a marker is ~size wide)
            if gx.min() > xv + 0.35 * size or gx.max() < xv - 0.35 * size:
                continue
            xm = 0.5 * (gx.min() + gx.max())
            dist = abs(xm - xv)
            if pred is not None and v in pred and pred_tol_px:
                dist += 0.5 * spacing * min(1.0, abs(np.median(gy) - pred[v]) / pred_tol_px)
            if best is None or dist < best[0]:
                best = (dist, float(np.median(gy)))
        if best is not None and best[0] <= 0.45 * spacing + (0.5 * spacing if pred is not None else 0):
            out[float(v)] = best[1]
    return out
