#!/usr/bin/env python3
"""
phenosight_aph.py — per-pot growth read-out from PhenoSight / PhenoVation APH images.

Replaces the PhenoSight software's ROI + mask step for small-plant trays
(e.g. Arabidopsis, 5 x 4 pots per tray). For every APH_<exp>_<traycode>_<round>.DAT
it:

  1. reads the two raw frames of the CropReporter "alpha" (absorptance) module:
     frame 0 = red reflectance (leaves dark), frame 1 = far-red reflectance (leaves
     bright); headerless uint16, stored column-major as (2 frames, ImageCols=2048,
     ImageRows=1500). Read row-major that is a 2048-tall x 1500-wide portrait image,
     i.e. the transpose of the landscape view PhenoVation / PlantCV show; all pot labels
     below are defined in this portrait frame,
  2. finds the tray grid in THIS image (template match + sub-pixel refinement of
     origin, x/y pitch and rotation) — so a changed camera height is handled
     automatically and every image carries its own scale (`pitch_px`),
  3. segments plant tissue with a background-corrected NIR/red index
         q = (NIR - k*R) / (NIR + k*R + off)
     where k is re-estimated per image as the modal NIR/R ratio of the tray
     (soil + plastic), so soil sits at q ~ 0 regardless of exposure, and plant
     pixels are taken by hysteresis (q > lo, connected to q > hi),
     v0.3: each tray is processed in time order and every pot carries a memory of its
     plant: pots with algae on the soil switch to a texture filter and a core threshold
     that rises with the algae level; objects must touch the plant's previous footprint
     or the pot centre; pots empty in the first images are called empty and later
     signal in them (algae, weeds) is zeroed,
  4. assigns plant pixels to pots: a connected object lying in one cell goes
     to that cell whole (leaves overhanging a wall stay with their plant); objects
     that bridge cells (touching rosettes, bolting stems) are split by a
     marker-controlled watershed seeded at the cell centres,
  5. writes one row per pot per image with area in pixels AND normalised to the
     cell size (area / pitch^2), which is what makes images taken at different
     camera heights comparable, and in mm^2 (pot pitch = 63.5 mm by default,
     the top width of the 63.5 mm pots used at BTI), plus shape and index
     traits; and a QC overlay.

Imaging time is read from the HDR_<exp>_<traycode>_<round>.INF file that the
CropReporter writes next to every image set (TimeStamp, Delphi date). The .DAT files
themselves carry no time. INF files found in RAW_DIR are used automatically.

Usage
-----
    python phenosight_aph.py RAW_DIR OUT_DIR [--workers 4] [--pitch-mm 63.5]
                             [--first-code 61633 --code-step 20]
                             [--lo 0.22 --hi 0.32] [--overlay-every 1]

Outputs (in OUT_DIR):
    pots.csv     one row per image x pot  (the table to analyse)
    images.csv   one row per image        (grid fit, scale, QC flags)
    overlays/    downscaled QC images: pot outlines + per-pot masks

Re-running is incremental: images already in images.csv are skipped, so you can
point it at a folder that is still syncing and re-run later.

Requires: numpy, scipy, scikit-image, opencv-python-headless, pandas.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import cv2
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.optimize import minimize
from skimage.filters import apply_hysteresis_threshold
from skimage.measure import regionprops
from skimage.segmentation import watershed

__version__ = "0.4.1"
SEG_VERSION = "0.4.0"       # bump only when per-image segmentation changes (invalidates checkpoints); 0.4.1 changed overlays only

# ----------------------------------------------------------------------------- constants
NFRAME, H, W = 2, 2048, 1500           # APH layout (validated on EXP30273 / EXP30281)
NROW, NCOL = 5, 4                      # pots per tray as seen in the APH image
FNAME_RE = re.compile(r"APH_(\d+)_(\d+)_(\d+)\.DAT$", re.I)

# APH image cell (row a = 0..4 top->bottom, col b = 0..3 left->right)  ->  layout label.
# Validated against the PhenoSight PS2 export of EXP30273 (per-pot area r = 0.98 for this
# mapping, <= 0.51 for any other orientation) together with the lab convention
# "A1 = ROI 20". The APH image is the PS2 image transposed, and 180 deg from the layout
# reading order: the letter runs bottom->top, the number right->left.
def cell_label(a: int, b: int) -> tuple[str, int]:
    area = "ABCDE"[NROW - 1 - a] + str(NCOL - b)
    roi = b * NROW + a + 1                     # = PhenoSight PS2 "Roi No"
    return area, roi


# ----------------------------------------------------------------------------- I/O
def read_aph(path: str) -> tuple[np.ndarray, np.ndarray]:
    raw = np.fromfile(path, dtype="<u2")
    if raw.size != NFRAME * H * W:
        raise ValueError(f"unexpected size {raw.size*2} bytes (expected {NFRAME*H*W*2})")
    a = (raw >> 4).reshape(NFRAME, H, W).astype(np.float32)
    return a[0], a[1]          # red, far-red (named 'nir' below: the plant-bright band)


def parse_name(path: str):
    m = FNAME_RE.search(os.path.basename(path))
    if not m:
        return None
    return dict(exp=int(m.group(1)), tray_code=int(m.group(2)), round=int(m.group(3)))


# ----------------------------------------------------------------------------- grid
def _local_norm(img, s=25):
    m = cv2.GaussianBlur(img, (0, 0), s)
    v = cv2.GaussianBlur(img * img, (0, 0), s) - m * m
    return (img - m) / np.sqrt(np.maximum(v, 1.0))


def _template(p, wf=0.12):
    Wt, Ht, w = int(round(NCOL * p)), int(round(NROW * p)), max(2, int(round(wf * p)))
    T = -np.ones((Ht + w, Wt + w), np.float32)
    for i in range(NCOL + 1):
        x = int(round(i * p)); T[:, x:x + w] = 1
    for j in range(NROW + 1):
        y = int(round(j * p)); T[y:y + w, :] = 1
    T -= T.mean(); T /= np.linalg.norm(T)
    return T


def _wall_score(z, x0, y0, px, py, th, ns=200):
    c, s = np.cos(th), np.sin(th)
    t = np.linspace(0, 1, ns); tot = 0.0
    for i in range(NCOL + 1):
        u = np.full(ns, i * px); v = t * NROW * py
        tot += ndi.map_coordinates(z, [y0 + s * u + c * v, x0 + c * u - s * v], order=1).mean()
    for j in range(NROW + 1):
        u = t * NCOL * px; v = np.full(ns, j * py)
        tot += ndi.map_coordinates(z, [y0 + s * u + c * v, x0 + c * u - s * v], order=1).mean()
    return tot / (NROW + NCOL + 2)


def find_grid(red, prange=(50, 140), ds=2):
    """Locate the 5x4 pot grid. Returns wall-centre-line origin, pitches, angle, score."""
    z = np.clip(_local_norm(cv2.GaussianBlur(red, (0, 0), 1.5)), -3, 3)
    zs = np.ascontiguousarray(z[::ds, ::ds])
    best = (-np.inf, None, None)
    for p in np.arange(prange[0], prange[1], 1.0):
        res = cv2.matchTemplate(zs, _template(p / ds), cv2.TM_CCORR)
        _, mx, _, loc = cv2.minMaxLoc(res)
        if mx > best[0]:
            best = (mx, p, loc)
    _, p, loc = best
    w = 0.12 * p
    x = np.array([loc[0] * ds + w / 2, loc[1] * ds + w / 2, p, p, 0.0])
    zr = np.clip(_local_norm(cv2.GaussianBlur(red, (0, 0), 2.5)), -3, 3)
    f = lambda v: -_wall_score(zr, *v)
    simplex = x + np.vstack([np.zeros(5), np.diag([4, 4, 2, 2, 0.01])])
    r = minimize(f, x, method="Nelder-Mead",
                 options=dict(xatol=0.2, fatol=1e-4, initial_simplex=simplex, maxiter=600))
    x0, y0, px, py, th = r.x
    return dict(x0=x0, y0=y0, px=px, py=py, theta=th, wall_score=-r.fun)


def grid_coords(g, shape=(H, W)):
    """Per-pixel tray coordinates (u along columns, v along rows), in cell units."""
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]].astype(np.float32)
    c, s = np.cos(g["theta"]), np.sin(g["theta"])
    dx, dy = xx - g["x0"], yy - g["y0"]
    u = (c * dx + s * dy) / g["px"]
    v = (-s * dx + c * dy) / g["py"]
    return u, v


def cell_center(g, a, b):
    c, s = np.cos(g["theta"]), np.sin(g["theta"])
    u, v = (b + 0.5) * g["px"], (a + 0.5) * g["py"]
    return g["x0"] + c * u - s * v, g["y0"] + s * u + c * v


# ----------------------------------------------------------------------------- segmentation
def _hmode(x, lo, hi, nb=400, sm=5):
    h, e = np.histogram(x, np.linspace(lo, hi, nb + 1))
    h = np.convolve(h, np.ones(sm) / sm, "same"); i = int(h.argmax())
    return 0.5 * (e[i] + e[i + 1])


def plant_index(red, nir, tray_mask, u, v, blur=0.8, off=30.0, max_offset=0.6, min_soil=0.08):
    """Background-corrected NIR/red index. The soil NIR/R ratio k is estimated PER POT
    (modal log-ratio of that pot's non-plant pixels), because soil reflectance differs
    between pots (moisture, algae, substrate) and drifts over time; the tray-wide value is
    the fallback when a pot is nearly covered by leaves."""
    R = cv2.GaussianBlur(red, (0, 0), blur); N = cv2.GaussianBlur(nir, (0, 0), blur)
    lr = np.log((N + 1.0) / (R + 1.0))
    bright = R > np.percentile(R[tray_mask], 20)          # ignore darkest, noisiest pixels
    sel = tray_mask & bright
    lk_tray = _hmode(lr[sel], -3, 1)
    ca = np.clip(np.floor(v).astype(int), 0, NROW - 1); cb = np.clip(np.floor(u).astype(int), 0, NCOL - 1)
    du = u - (np.floor(u) + 0.5); dv = v - (np.floor(v) + 0.5)
    interior = tray_mask & bright & (np.abs(du) < 0.38) & (np.abs(dv) < 0.38)
    lk = np.full((NROW, NCOL), lk_tray, np.float32)
    for a in range(NROW):
        for b in range(NCOL):
            m = interior & (ca == a) & (cb == b)
            x = lr[m]
            x = x[x < lk_tray + max_offset]
            if x.size >= min_soil * m.sum() and x.size > 50:
                lk[a, b] = _hmode(x, lk_tray - 1.0, lk_tray + max_offset, nb=160)
    k = np.exp(lk[ca, cb])
    q = (N - k * R) / (N + k * R + off)
    ndvi = (N - R) / (N + R + 1.0)
    return q, ndvi, float(np.exp(lk_tray)), np.exp(lk), R, N


def segment(q, roi, lo, hi, min_px):
    """v0.2 single-image segmentation (kept for reference / --no-tracking)."""
    m = apply_hysteresis_threshold(q, lo, hi) & roi
    m = ndi.binary_opening(m, iterations=1) | (m & (q > hi))   # drop 1-px speckle, keep cores
    return _drop_small(m, min_px)


def _drop_small(m, min_px):
    lab, n = ndi.label(m)
    if n:
        sizes = ndi.sum(m, lab, np.arange(1, n + 1))
        keep = np.zeros(n + 1, bool); keep[1:] = sizes >= min_px
        m = keep[lab]
    return m


def _disk(r):
    r = max(1, int(round(r)))
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return x * x + y * y <= r * r


# --- per-pot footprint memory, stored in cell coordinates so it survives camera moves ----
FP_PAD, FP_N = 0.35, 96          # footprint raster covers the cell +/- 0.35 cell, 96 x 96


def _fp_index(u, v, a, b):
    """Map tray coordinates to footprint-raster indices of pot (a, b); -1 if outside."""
    iu = np.floor((u - (b - FP_PAD)) / (1 + 2 * FP_PAD) * FP_N).astype(int)
    iv = np.floor((v - (a - FP_PAD)) / (1 + 2 * FP_PAD) * FP_N).astype(int)
    ok = (iu >= 0) & (iu < FP_N) & (iv >= 0) & (iv < FP_N)
    return np.where(ok, iv, -1), np.where(ok, iu, -1)


def _fp_sample(F, iv, iu):
    out = np.zeros(iv.shape, bool)
    ok = iv >= 0
    out[ok] = F[iv[ok], iu[ok]]
    return out


def segment_tracked(q, u, v, pitch, state, dt_days, p):
    """Segmentation with texture filtering, per-pot adaptive thresholds and a temporal gate.

    * Texture: algae / moss on the soil is a fine-grained mesh of 1-4 px granules at
      intermediate q, whereas leaf blades are smooth and several pixels wide. Cores are
      therefore found on a grey-level OPENING of q (disk of open_frac * pitch), which
      flattens the granular algae signal but leaves leaf blades intact.
    * Adaptive threshold: per pot, the 99th percentile of opened q on that pot's soil
      (interior minus the plant cores) measures how much algae there is; the core
      threshold is raised to that level + margin (between --hi and --hi-max), so it climbs
      as algae develops, and only in the pots that have it. The edge threshold follows the
      95th percentile of raw q on the soil. hi_used / lo_used / algae_frac are reported.
    * Leaf margins (mixed pixels) are added back only within edge_frac * pitch of a core.
    * Temporal gate (see gate_and_update): objects must touch the plant's previous
      footprint, grown by a margin that scales with time since the last image, or the
      centre of the pot.
    """
    ca = np.floor(v).astype(int); cb = np.floor(u).astype(int)
    du = u - (np.floor(u) + 0.5); dv = v - (np.floor(v) + 0.5)
    roi = (u >= -FP_PAD) & (u <= NCOL + FP_PAD) & (v >= -FP_PAD) & (v <= NROW + FP_PAD)
    inside = (ca >= 0) & (ca < NROW) & (cb >= 0) & (cb < NCOL)
    interior = inside & (np.abs(du) < 0.38) & (np.abs(dv) < 0.38)
    min_px = max(4, int(round(p["min_cellfrac"] * pitch * pitch)))

    ys0, xs0 = np.nonzero(roi)
    y0, y1, x0, x1 = ys0.min(), ys0.max() + 1, xs0.min(), xs0.max() + 1
    qo = np.full(q.shape, -1.0, np.float32)
    qo[y0:y1, x0:x1] = ndi.grey_opening(q[y0:y1, x0:x1], footprint=_disk(p["open_frac"] * pitch))

    # Algae indicator per pot: fraction of the pot's outer soil ring (0.22-0.38 cell from
    # the centre, away from where the seedling sits, minus anything near plant tissue)
    # with q above the edge threshold. Clean soil gives 0.000; algae/moss 0.01-0.06.
    cac = np.clip(ca, 0, NROW - 1); cbc = np.clip(cb, 0, NCOL - 1)
    s = cv2.GaussianBlur(q, (0, 0), p["smooth"])
    near_plant = ndi.binary_dilation(_drop_small((s > p["hi"]) & roi, min_px), structure=_disk(0.08 * pitch))
    cheb = np.maximum(np.abs(du), np.abs(dv))
    ring = inside & (cheb > 0.22) & (cheb < 0.38) & ~near_plant
    algae = np.full((NROW, NCOL), np.nan)
    for a in range(NROW):
        for b in range(NCOL):
            sel = ring & (ca == a) & (cb == b)
            if sel.sum() > 100:
                algae[a, b] = float(np.mean(q[sel] > p["lo"]))
    # Texture filter (grey opening) per pot: on once the rosette is large enough that its
    # leaves are wider than the filter (>= seedling_cellfrac), or earlier (>= algae_min_cellfrac)
    # in a pot that has shown algae. Below that, seedlings are about as small as algae
    # granules and the filter would erase them, so pots keep plain hysteresis on q (v0.2
    # behaviour, best agreement with the PS2 masks); the temporal gate keeps algae away from
    # the seedling there.
    use_open = np.zeros((NROW, NCOL), bool)
    for a in range(NROW):
        for b in range(NCOL):
            st = state.setdefault((a, b), {})
            if np.isfinite(algae[a, b]) and algae[a, b] >= p["algae_switch"]:
                st["algae_seen"] = True
            la = st.get("last_area", 0.0)
            use_open[a, b] = la >= p["seedling_cellfrac"] or (st.get("algae_seen", False) and la >= p["algae_min_cellfrac"])
    sig = np.where(use_open[cac, cbc], qo, s)
    prelim = _drop_small((sig > p["hi"]) & roi, min_px)
    excl = ndi.binary_dilation(prelim, structure=_disk(0.10 * pitch))

    hi_map = np.full((NROW, NCOL), p["hi"], np.float32); lo_map = np.full((NROW, NCOL), p["lo"], np.float32)
    diag = {}
    for a in range(NROW):
        for b in range(NCOL):
            soil = interior & (ca == a) & (cb == b) & ~excl
            if soil.sum() > 200:
                o99 = float(np.percentile(sig[soil], 98))
                r95 = float(np.percentile(q[soil], 95))
                frac = float(np.mean(q[soil] > p["lo"]))
            else:                                   # pot covered by leaves: keep last level
                o99, r95, frac = state.get((a, b), {}).get("soil", (np.nan, np.nan, np.nan))
            base_hi = p["hi"] if use_open[a, b] else p["seed_hi"]
            base_lo = p["lo"] if use_open[a, b] else p["seed_lo"]
            hi_map[a, b], lo_map[a, b] = base_hi, base_lo
            if np.isfinite(o99):
                hi_map[a, b] = np.clip(o99 + p["hi_margin"], base_hi, p["hi_max"])
                lo_map[a, b] = np.clip(r95 + p["lo_margin"], base_lo, hi_map[a, b] - 0.04)
            state.setdefault((a, b), {})["soil"] = (o99, r95, frac)
            diag[(a, b)] = dict(algae_ring_frac=algae[a, b], soil_sig_p98=o99, soil_q_p95=r95,
                                hi_used=float(hi_map[a, b]), lo_used=float(lo_map[a, b]))
    cac = np.clip(ca, 0, NROW - 1); cbc = np.clip(cb, 0, NCOL - 1)
    core = _drop_small((sig > hi_map[cac, cbc]) & roi, min_px)
    for (a, b), d in diag.items():
        d["texture_filter"] = bool(use_open[a, b])
    edge = ndi.binary_dilation(core, structure=_disk(max(1.5, p["edge_frac"] * pitch))) & (q > lo_map[cac, cbc])
    big = core | (edge & roi)
    # seedling pots: plain hysteresis on q (v0.2 behaviour; best agreement with the PS2
    # masks for 1-3-leaf seedlings)
    low = (q > lo_map[cac, cbc]) & roi
    lab_low, _ = ndi.label(low)
    hit = np.unique(lab_low[(q > hi_map[cac, cbc]) & low]); hit = hit[hit > 0]
    hyst = np.isin(lab_low, hit) & low
    hyst = ndi.binary_opening(hyst, iterations=1) | (hyst & (q > hi_map[cac, cbc]))
    seedling_px = ~use_open[cac, cbc]
    mask = _drop_small(np.where(seedling_px, hyst, big), min_px)

    # previous footprints, grown by a time-dependent margin, for the temporal gate
    grow_cell = min(p["grow_max"], p["grow_base"] + p["grow_per_day"] * dt_days)
    idx = {}
    for a in range(NROW):
        for b in range(NCOL):
            sel = roi & (u >= b - FP_PAD) & (u < b + 1 + FP_PAD) & (v >= a - FP_PAD) & (v < a + 1 + FP_PAD)
            ys, xs = np.nonzero(sel)
            iv, iu = _fp_index(u[ys, xs], v[ys, xs], a, b)
            idx[(a, b)] = (ys, xs, iv, iu)
            F = state.get((a, b), {}).get("F")
            if F is not None:
                state[(a, b)]["Fg"] = ndi.binary_dilation(F, structure=_disk(grow_cell / (1 + 2 * FP_PAD) * FP_N))
    # weak "bridge" signal (q above the pot's edge threshold, granular or not): never
    # counted as plant, only used to decide which leaf blades belong to a tracked rosette,
    # so blades joined to it by thin petioles are not dropped by the temporal gate
    bridge = (q > lo_map[cac, cbc]) & roi
    return mask, bridge, idx, diag


def gate_and_update(lab, u, v, idx, state, p, bridge=None, pitch=100.0):
    """Drop pot objects that touch neither the grown footprint nor the pot centre;
    then store each pot's new footprint (in cell coordinates)."""
    du = u - (np.floor(u) + 0.5); dv = v - (np.floor(v) + 0.5)
    rejected = {}
    for a in range(NROW):
        for b in range(NCOL):
            pid = a * NCOL + b + 1
            ys, xs, iv, iu = idx[(a, b)]
            st = state.setdefault((a, b), {})
            pm = lab == pid
            if not pm.any():
                rejected[(a, b)] = 0
                st["last_area"] = 0.0
                continue
            comp, n = ndi.label(pm if bridge is None else (pm | bridge))
            r_c = p["centre_r_tracked"] if st.get("F") is not None else p["centre_r_new"]
            centre = (np.floor(u) == b) & (np.floor(v) == a) & (du * du + dv * dv <= r_c ** 2)
            gate = centre.copy()
            if st.get("Fg") is not None:
                gate[ys, xs] |= _fp_sample(st["Fg"], iv, iu)
            keep_ids = np.unique(comp[gate & pm]); keep_ids = keep_ids[keep_ids > 0]
            keep = np.isin(comp, keep_ids) & pm
            rejected[(a, b)] = int(pm.sum() - keep.sum())
            lab[pm & ~keep] = 0
            st["last_area"] = keep.sum() / pitch ** 2
            if keep.any():
                F = np.zeros((FP_N, FP_N), bool)
                inpot = keep[ys, xs] & (iv >= 0)
                F[iv[inpot], iu[inpot]] = True
                F = ndi.binary_closing(F, structure=_disk(2))
                st["F"] = F
            st.pop("Fg", None)
    return lab, rejected


# ----------------------------------------------------------------------------- pot assignment
def assign_pots(mask, g, u, v, core_r=0.30, own_frac=0.85):
    """Return a label image: 0 = background, 1..20 = pot index (a*NCOL + b + 1)."""
    cell_a = np.floor(v).astype(int); cell_b = np.floor(u).astype(int)
    inside = (cell_a >= 0) & (cell_a < NROW) & (cell_b >= 0) & (cell_b < NCOL)
    cell_id = np.where(inside, cell_a * NCOL + cell_b + 1, 0)
    # central disc of each cell = seed region for splitting bridged objects
    du = u - (np.floor(u) + 0.5); dv = v - (np.floor(v) + 0.5)
    core = inside & (du * du + dv * dv <= core_r ** 2)

    out = np.zeros(mask.shape, np.int16)
    lab, n = ndi.label(mask)
    if n == 0:
        return out, 0
    objs = ndi.find_objects(lab)
    n_split = 0
    for i, sl in enumerate(objs, start=1):
        obj = lab[sl] == i
        cid = cell_id[sl][obj]
        cid_in = cid[cid > 0]
        if cid_in.size == 0:
            continue                                     # entirely outside the tray
        counts = np.bincount(cid_in, minlength=NROW * NCOL + 1)
        top = int(counts.argmax())
        # seeds: cells whose central disc this object occupies
        seeds = np.unique(cell_id[sl][obj & core[sl]])
        seeds = seeds[seeds > 0]
        if counts[top] >= own_frac * cid.size or len(seeds) <= 1:
            owner = top if len(seeds) == 0 else int(seeds[0]) if len(seeds) == 1 else top
            out[sl][obj] = owner
            continue
        # bridged object: watershed on distance transform from seeded cores
        n_split += 1
        markers = np.zeros(obj.shape, np.int32)
        cs = core[sl] & obj
        markers[cs] = cell_id[sl][cs]
        dist = ndi.distance_transform_edt(obj)
        ws = watershed(-dist, markers, mask=obj)
        tgt = out[sl]
        tgt[obj] = ws[obj]
    return out, n_split


# ----------------------------------------------------------------------------- PS2 / spectral traits
# Sibling files written by the CropReporter for the same image set (same stem, same pixel grid):
#   PMT_* (9 frames): 0 Fdark, 1 F0, 2 Fm, 3 Fdark+sat, 4 Flight, 5 Fs', 6 Fm', 7 Flight+sat, 8 F0'
#   SPC_* (3 frames): Anth (~540 nm), FarRed / red edge (~710 nm), NIR (~770 nm)
#   CLR_* (3 frames): Red, Green, Blue
# Reflectance = DN / 65535 (bands are pre-calibrated by their shutters). The PMT frame order was
# checked against the PhenoVation PS2 export of EXP30273: F0/Fm/Fs'/Fm'/F0' frame means track the
# export's values at a constant gain (~0.72) and the derived traits agree pot by pot (README).
PMT_FRAMES, SPC_FRAMES, CLR_FRAMES = 9, 3, 3
FLUO_RANGES = dict(Fv_Fm=(0, 1), Fq_Fm=(0, 1), NPQ=(0, 10), qP=(0, 1.5), qN=(0, 1.5), qL=(0, 1.5),
                   phiNO=(0, 1), phiNPQ=(-0.5, 1), npq_t=(0, 15))
FLUO_MIN_FM = 500        # pixels must have Fm above this (DN) to enter the fluorescence averages


def _read_crop(path, n, ys, xs):
    """Window of an n-frame CropReporter file; None if the file is missing, truncated or not
    readable (e.g. a cloud placeholder that has not been downloaded)."""
    if not path:
        return None
    try:
        if os.path.getsize(path) != n * H * W * 2:
            return None
        out = np.empty((n, ys.stop - ys.start, xs.stop - xs.start), np.float32)
        with open(path, "rb") as fh:          # plain reads (an mmap of a cloud placeholder can crash)
            for i in range(n):
                fh.seek((i * H + ys.start) * W * 2)
                buf = np.frombuffer(fh.read((ys.stop - ys.start) * W * 2), dtype="<u2")
                out[i] = buf.reshape(-1, W)[:, xs]
        return out
    except (OSError, ValueError):
        return None


def sibling(path, prefix):
    d, b = os.path.split(path)
    p = os.path.join(d, prefix + b[3:])
    return p if os.path.exists(p) else None


def pixel_traits(pmt=None, spc=None, clr=None):
    T = {}
    with np.errstate(all="ignore"):
        if pmt is not None:
            F0, Fm, Fs, Fmp, F0p = [np.maximum(pmt[i], 1.0) for i in (1, 2, 5, 6, 8)]
            fqfm = (Fmp - Fs) / Fmp; npq = Fm / Fmp - 1
            qp = (Fmp - Fs) / (Fmp - F0p); ql = qp * F0p / Fs
            phino = 1 / (npq + 1 + ql * (Fm / F0 - 1))
            T.update(Fv_Fm=(Fm - F0) / Fm, Fq_Fm=fqfm, NPQ=npq, qP=qp, qN=1 - (Fmp - F0p) / (Fm - F0),
                     qL=ql, phiNO=phino, phiNPQ=1 - fqfm - phino, npq_t=4.88 / (Fmp / F0p - 1) - 1)
        eps = 1e-6
        if spc is not None:
            Ra, Rfr, Rn = [np.maximum(spc[i] / 65535.0, eps) for i in range(3)]
            T.update(ChlIdx=(1 / Rfr - 1 / Rn) * Rn, AriIdx=(1 / Ra - 1 / Rfr) * Rn, refl_anth=Ra,
                     refl_farred=Rfr, refl_nir=Rn)
            if clr is not None:
                Rr = np.maximum(clr[0] / 65535.0, eps)
                T["NDVI"] = (Rn - Rr) / (Rn + Rr)
        if clr is not None:
            Rr, Rg, Rb = [np.maximum(clr[i] / 65535.0, eps) for i in range(3)]
            mx, mn = np.maximum(np.maximum(Rr, Rg), Rb), np.minimum(np.minimum(Rr, Rg), Rb)
            hue = np.where(mx == Rr, (Rg - Rb) / (mx - mn + eps) % 6,
                  np.where(mx == Rg, (Rb - Rr) / (mx - mn + eps) + 2, (Rr - Rg) / (mx - mn + eps) + 4)) * 60
            T.update(refl_red=Rr, refl_green=Rg, refl_blue=Rb, Hue=hue, GreenIdx=(2 * Rg - Rr - Rb) / (Rr + Rg + Rb))
    return T


def pot_trait_values(m, T, Fm=None):
    """Mean of each pixel trait over pot mask m. Fluorescence traits use only pixels with
    Fm > FLUO_MIN_FM and values inside FLUO_RANGES (vendor-style gating)."""
    out = {}
    mf = m & (Fm > FLUO_MIN_FM) if Fm is not None else None
    if mf is not None:
        out["n_px_fluo"] = int(mf.sum())
    for k, v in T.items():
        if k in FLUO_RANGES:
            x = v[mf]; lo, hi = FLUO_RANGES[k]; x = x[(x >= lo) & (x <= hi)]
            out[k] = float(x.mean()) if x.size >= 5 else np.nan
        else:
            x = v[m]
            out[k] = float(x.mean()) if x.size else np.nan
    return out


# ----------------------------------------------------------------------------- per-image
DEFAULTS = dict(lo=0.22, hi=0.32, seed_lo=0.22, seed_hi=0.32, hi_max=0.45, hi_margin=0.03, lo_margin=0.02, smooth=1.0,
                edge_frac=0.025, open_frac=0.022, seedling_cellfrac=0.06, algae_switch=0.005, algae_min_cellfrac=0.02, min_cellfrac=2e-4,
                grow_base=0.06, grow_per_day=0.05, grow_max=0.35,     # footprint margin, cell units
                centre_r_new=0.40, centre_r_tracked=0.25,
                establish_days=7.0, establish_n=14, establish_k=3, establish_min=0.002,     # empty-pot call
                tracking=True, traits=True)


def process_image(path, state=None, dt_days=0.5, params=None, overlay_dir=None):
    p = dict(DEFAULTS, **(params or {}))
    state = {} if state is None else state
    meta = parse_name(path)
    t0 = time.time()
    red, nir = read_aph(path)
    # grid: search near the previous image's pitch first (fast); fall back to the full
    # search when the camera has moved (fit clearly worse than the previous image's)
    g = None
    last = state.get("_grid")
    if last is not None:
        g = find_grid(red, prange=(last["px"] - 8, last["px"] + 8))
        if g["wall_score"] < 0.85 * last["wall_score"]:
            g = None
    if g is None:
        g = find_grid(red)
    state["_grid"] = dict(px=g["px"], wall_score=g["wall_score"])
    # crop to the tray (+ overhang margin); all later steps work on this window
    c, s_ = np.cos(g["theta"]), np.sin(g["theta"])
    xs_, ys_ = [], []
    for uu, vv in [(-0.7, -0.7), (NCOL + 0.7, -0.7), (-0.7, NROW + 0.7), (NCOL + 0.7, NROW + 0.7)]:
        xs_.append(g["x0"] + c * uu * g["px"] - s_ * vv * g["py"]); ys_.append(g["y0"] + s_ * uu * g["px"] + c * vv * g["py"])
    cx0, cx1 = max(0, int(min(xs_))), min(W, int(max(xs_)) + 1)
    cy0, cy1 = max(0, int(min(ys_))), min(H, int(max(ys_)) + 1)
    red, nir = red[cy0:cy1, cx0:cx1], nir[cy0:cy1, cx0:cx1]
    win = (slice(cy0, cy1), slice(cx0, cx1))
    pmt = spc = clr = None
    if p.get("traits", True):
        pmt = _read_crop(sibling(path, "PMT"), PMT_FRAMES, *win)
        spc = _read_crop(sibling(path, "SPC"), SPC_FRAMES, *win)
        clr = _read_crop(sibling(path, "CLR"), CLR_FRAMES, *win)
    g_full = dict(g)
    g = dict(g, x0=g["x0"] - cx0, y0=g["y0"] - cy0)
    u, v = grid_coords(g, shape=red.shape)
    pitch = 0.5 * (g["px"] + g["py"])
    tray = (u >= 0) & (u <= NCOL) & (v >= 0) & (v <= NROW)
    q, ndvi, k, kcell, R, N = plant_index(red, nir, tray, u, v)
    if p["tracking"]:
        mask, bridge, idx, diag = segment_tracked(q, u, v, pitch, state, dt_days, p)
        lab, n_split = assign_pots(mask, g, u, v)
        lab, rejected = gate_and_update(lab, u, v, idx, state, p, bridge, pitch)
    else:
        roi = (u >= -0.35) & (u <= NCOL + 0.35) & (v >= -0.35) & (v <= NROW + 0.35)
        mask = segment(q, roi, p["lo"], p["hi"], max(4, int(round(p["min_cellfrac"] * pitch ** 2))))
        lab, n_split = assign_pots(mask, g, u, v)
        diag, rejected = {}, {}

    rows = []
    PT = pixel_traits(pmt, spc, clr) if any(x is not None for x in (pmt, spc, clr)) else {}
    props = {pp.label: pp for pp in regionprops(lab.astype(np.int32))}
    for a in range(NROW):
        for b in range(NCOL):
            pid = a * NCOL + b + 1
            area_lbl, roi_no = cell_label(a, b)
            cx, cy = cell_center(g, a, b)
            r = dict(**meta, file=os.path.basename(path), cell_row=a, cell_col=b,
                     Area=area_lbl, Roi_No=roi_no, pitch_px=pitch,
                     soil_fr_red=float(kcell[a, b]), rejected_cellfrac=rejected.get((a, b), 0) / pitch ** 2,
                     **diag.get((a, b), {}))
            pp = props.get(pid)
            if pp is None:
                r.update(area_px=0, area_cellfrac=0.0)
            else:
                m = lab[pp.slice] == pid
                qv = q[pp.slice][m]
                hull = pp.area_convex
                r.update(
                    area_px=int(pp.area),
                    area_cellfrac=pp.area / pitch ** 2,
                    hull_cellfrac=hull / pitch ** 2,
                    solidity=pp.area / hull if hull else np.nan,
                    perimeter_cell=pp.perimeter / pitch,
                    diameter_cell=pp.feret_diameter_max / pitch,
                    n_objects=int(ndi.label(m)[1]),
                    centroid_offset_cell=float(np.hypot(pp.centroid[1] - cx, pp.centroid[0] - cy) / pitch),
                    q_mean=float(qv.mean()), q_median=float(np.median(qv)),
                    abs_mean=float(np.mean(1.0 - R[pp.slice][m] / np.maximum(N[pp.slice][m], 1.0))),
                    farred_mean=float(N[pp.slice][m].mean()), red_mean=float(R[pp.slice][m].mean()),
                )
            if PT and pp is not None:
                r.update(pot_trait_values(lab[pp.slice] == pid, {k: v[pp.slice] for k, v in PT.items()},
                                          None if pmt is None else pmt[2][pp.slice]))
            rows.append(r)

    img = dict(**meta, file=os.path.basename(path), x0=g_full["x0"], y0=g_full["y0"], px=g["px"], py=g["py"],
               theta_deg=np.degrees(g["theta"]), pitch_px=pitch, wall_score=g["wall_score"],
               k_bg=k, k_cell_sd=float(np.std(np.log(kcell))), red_mean=float(R[tray].mean()),
               farred_mean=float(N[tray].mean()), plant_frac=float((lab > 0)[tray].mean()), n_split=n_split,
               pots_with_plant=sum(1 for r in rows if r["area_px"] > 0), dt_days=dt_days,
               has_pmt=pmt is not None, has_spc=spc is not None, has_clr=clr is not None,
               seconds=round(time.time() - t0, 2), version=__version__)
    if overlay_dir:
        write_overlay(os.path.join(overlay_dir, os.path.basename(path)[:-4] + ".jpg"), nir, lab, g)
    return img, rows


def process_tray(paths, times=None, params=None, overlay_dir=None, ckpt=None, deadline=None):
    """Process one tray's images in time order, carrying each pot's footprint forward.
    `times` maps file basename -> Timestamp (missing: rounds are assumed 12 h apart).
    With `ckpt` (a file path) progress is saved after every image and a later call resumes
    where this one stopped; `deadline` (time.time() value) stops cleanly before it.
    Returns (imgs, rows, complete)."""
    import pickle
    p = dict(DEFAULTS, **(params or {}))
    times = times or {}
    def key(f):
        t = times.get(os.path.basename(f))
        return (0, t.value, 0) if t is not None and pd.notna(t) else (1, 0, parse_name(f)["round"])
    paths = sorted(paths, key=key)
    names = [os.path.basename(f) for f in paths]
    C = dict(done=[], state={}, imgs=[], rows=[], prev_t=None, prev_r=None)
    if ckpt and os.path.exists(ckpt):
        with open(ckpt, "rb") as fh:
            old = pickle.load(fh)
        # resume only if the same version/params processed a prefix of the current order
        segp = {k: v for k, v in p.items() if not k.startswith("establish")}
        oldp = {k: v for k, v in (old.get("params") or {}).items() if not k.startswith("establish")}
        if (old.get("version") in (SEG_VERSION,) and oldp == segp
                and old["done"] == names[:len(old["done"])]):
            C = old
    for f in paths[len(C["done"]):]:
        if deadline is not None and time.time() > deadline:
            break
        t = times.get(os.path.basename(f)); rnd = parse_name(f)["round"]
        if C["prev_t"] is not None and t is not None and pd.notna(t):
            dt = max(0.0, (t - C["prev_t"]).total_seconds() / 86400)
        elif C["prev_r"] is not None:
            dt = 0.5 * max(1, rnd - C["prev_r"])
        else:
            dt = 0.5
        img, r = process_image(f, C["state"], dt, p, overlay_dir)
        C["imgs"].append(img); C["rows"].extend(r); C["done"].append(os.path.basename(f))
        C["prev_t"], C["prev_r"] = t, rnd
        if ckpt:
            C["version"], C["params"] = SEG_VERSION, p
            tmp = ckpt + ".tmp"
            with open(tmp, "wb") as fh:
                pickle.dump(C, fh)
            os.replace(tmp, ckpt)
    return C["imgs"], C["rows"], len(C["done"]) == len(paths)


def call_pot_status(P, p):
    """Empty-pot call per tray x pot. A pot holds a plant if it is detected (area >=
    establish_min of a cell, ~8 mm2) in at least establish_k images during the tray's first
    establish_days days (or its first establish_n images when imaging times are unknown).
    Seedlings that come up a few days late are therefore kept, while signal that first
    appears after that window in a pot with no plant (algae, a weed, a neighbour's leaf) is
    zeroed. The raw values stay in area_cellfrac_raw; late_signal_cellfrac holds the largest
    late signal for review. Trays imaged for less than the window are 'undetermined'."""
    P = P.copy()
    if "area_cellfrac_raw" not in P:
        P["area_cellfrac_raw"] = P["area_cellfrac"]
    col = "area_cellfrac_raw"
    timed = "datetime" in P and P["datetime"].notna().all()
    P = P.sort_values(["exp", "tray_code", "datetime" if timed else "round"])
    if timed:
        t0 = P.groupby(["exp", "tray_code"])["datetime"].transform("min")
        P["_day"] = (pd.to_datetime(P["datetime"]) - pd.to_datetime(t0)).dt.total_seconds() / 86400
    else:
        P["_day"] = np.nan
    def one(d):
        if timed:
            first = d[d["_day"] <= p["establish_days"]]
            covered = d["_day"].max() >= p["establish_days"]
        else:
            first = d.head(p["establish_n"])
            covered = len(d) >= p["establish_n"]
        hits = int((first[col] >= p["establish_min"]).sum())
        if hits >= p["establish_k"]:
            st = "plant"
        elif not covered:
            st = "undetermined"            # not imaged long enough yet to call the pot empty
        else:
            st = "empty"
        return pd.Series(dict(pot_status=st,
                              late_signal_cellfrac=float(d[col].max()) if st == "empty" else 0.0))
    S = P.groupby(["exp", "tray_code", "Area"]).apply(one).reset_index()
    P = P.drop(columns=["pot_status", "late_signal_cellfrac", "_day"], errors="ignore").merge(
        S, on=["exp", "tray_code", "Area"])
    P["area_cellfrac"] = np.where(P["pot_status"] == "empty", 0.0, P["area_cellfrac_raw"])
    return P


# ----------------------------------------------------------------------------- QC overlay
_PALETTE = np.array([[0, 0, 0]] + [[int(255 * c) for c in col] for col in
                    [(0.12, .47, .71), (1, .5, .05), (.17, .63, .17), (.84, .15, .16), (.58, .4, .74),
                     (.55, .34, .29), (.89, .47, .76), (.5, .5, .5), (.74, .74, .13), (.09, .75, .81)] * 2],
                    np.uint8)


def write_overlay(path, nir, lab, g, pad=0.6, out_w=700, door_view=True):
    """QC overlay. door_view=True (default from v0.4.1) draws the tray as seen standing at the
    chamber door: the camera's own view (the transpose of the array read here) turned 180 deg,
    white column at the top, A1 upper left, E4 lower right. door_view=False keeps the raw
    array orientation used up to v0.4.0 (a mirror image of the camera view)."""
    c, s = np.cos(g["theta"]), np.sin(g["theta"])
    corners = []
    for uu, vv in [(-pad, -pad), (NCOL + pad, -pad), (-pad, NROW + pad), (NCOL + pad, NROW + pad)]:
        corners.append((g["x0"] + c * uu * g["px"] - s * vv * g["py"], g["y0"] + s * uu * g["px"] + c * vv * g["py"]))
    xs, ys = zip(*corners)
    x0, x1 = max(0, int(min(xs))), min(W, int(max(xs)))
    y0, y1 = max(0, int(min(ys))), min(H, int(max(ys)))
    base = nir[y0:y1, x0:x1]
    base = np.clip(base / max(np.percentile(base, 99.7), 1) * 255, 0, 255).astype(np.uint8)
    rgb = np.dstack([base] * 3)
    L = lab[y0:y1, x0:x1]
    col = _PALETTE[L.clip(0, 20)]
    edge = (L > 0) & (ndi.grey_erosion(L, size=3) != L)
    fill = L > 0
    rgb[fill] = (0.65 * rgb[fill] + 0.35 * col[fill]).astype(np.uint8)
    rgb[edge] = col[edge]
    for i in range(NCOL + 1):                                    # grid lines
        p1 = (g["x0"] + c * i * g["px"] - x0, g["y0"] + s * i * g["px"] - y0)
        p2 = (p1[0] - s * NROW * g["py"], p1[1] + c * NROW * g["py"])
        cv2.line(rgb, tuple(map(int, p1)), tuple(map(int, p2)), (255, 255, 255), 1)
    for j in range(NROW + 1):
        p1 = (g["x0"] - s * j * g["py"] - x0, g["y0"] + c * j * g["py"] - y0)
        p2 = (p1[0] + c * NCOL * g["px"], p1[1] + s * NCOL * g["px"])
        cv2.line(rgb, tuple(map(int, p1)), tuple(map(int, p2)), (255, 255, 255), 1)
    if not door_view:
        for a in range(NROW):
            for b in range(NCOL):
                cx, cy = cell_center(g, a, b)
                cv2.putText(rgb, cell_label(a, b)[0], (int(cx - x0 - 0.42 * g["px"]), int(cy - y0 - 0.3 * g["py"])),
                            cv2.FONT_HERSHEY_SIMPLEX, 0.0045 * g["px"], (255, 255, 0), 1, cv2.LINE_AA)
    sc = out_w / rgb.shape[1]
    rgb = cv2.resize(rgb, (out_w, int(rgb.shape[0] * sc)), interpolation=cv2.INTER_AREA)
    if door_view:
        rgb = _door_view(rgb, g, (x0, y0), sc)
    cv2.imwrite(path, cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, 85])


def _door_view(rgb, g, origin, sc):
    """Transpose + 180 deg turn of an overlay (new[i, j] = old[h-1-j, w-1-i]), cell labels in the
    upper-left corner of each pot, and a 'white column' / 'DOOR' strip at the top / bottom."""
    h, w = rgb.shape[:2]
    c, s = np.cos(g["theta"]), np.sin(g["theta"])
    def to_new(u, v):
        x = (g["x0"] + c * u * g["px"] - s * v * g["py"] - origin[0]) * sc
        y = (g["y0"] + s * u * g["px"] + c * v * g["py"] - origin[1]) * sc
        return h - 1 - y, w - 1 - x
    new = np.ascontiguousarray(rgb.transpose(1, 0, 2)[::-1, ::-1])
    nw = new.shape[1]
    fs = 0.0045 * g["px"] * sc; th = max(1, int(round(fs * 2)))
    for a in range(NROW):
        for b in range(NCOL):
            cs = [to_new(b + du, a + dv) for du in (0, 1) for dv in (0, 1)]
            xl = min(p[0] for p in cs); yt = min(p[1] for p in cs); cell = max(p[0] for p in cs) - xl
            org = (int(xl + 0.08 * cell), int(yt + 0.22 * cell))
            cv2.putText(new, cell_label(a, b)[0], org, cv2.FONT_HERSHEY_SIMPLEX, fs, (0, 0, 0), th + 2, cv2.LINE_AA)
            cv2.putText(new, cell_label(a, b)[0], org, cv2.FONT_HERSHEY_SIMPLEX, fs, (255, 255, 0), th, cv2.LINE_AA)
    band = max(18, nw // 28); fsb = band / 30
    top = np.full((band, nw, 3), 235, np.uint8); bot = np.full((band, nw, 3), 60, np.uint8)
    for im, txt, colr in [(top, "white column (far side)", (40, 40, 40)), (bot, "DOOR (viewer stands here)", (255, 255, 255))]:
        (tw, tht), _ = cv2.getTextSize(txt, cv2.FONT_HERSHEY_SIMPLEX, fsb, 1)
        cv2.putText(im, txt, ((nw - tw) // 2, (band + tht) // 2), cv2.FONT_HERSHEY_SIMPLEX, fsb, colr, 1, cv2.LINE_AA)
    return np.vstack([top, new, bot])



# ----------------------------------------------------------------------------- imaging times
def read_inf_time(path):
    """TimeStamp from a CropReporter HDR_*.INF (Delphi serial date, decimal comma;
    instrument local time)."""
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("TimeStamp="):
                days = float(line.split("=", 1)[1].strip().replace(",", "."))
                return (pd.Timestamp("1899-12-30") + pd.to_timedelta(days, unit="D")).round("s")
    return pd.NaT


INF_RE = re.compile(r"HDR_(\d+)_(\d+)_(\d+)\.INF$", re.I)


def times_from_inf_dir(folder):
    rows = []
    for f in os.listdir(folder):
        m = INF_RE.search(f)
        if m:
            rows.append(dict(tray_code=int(m.group(2)), round=int(m.group(3)),
                             datetime=read_inf_time(os.path.join(folder, f))))
    return pd.DataFrame(rows, columns=["tray_code", "round", "datetime"])


def read_times(path):
    """Accepts: a folder of HDR_*.INF files, a single .INF, a PS2 export, or a CSV with
    tray_code, round, datetime.
    """
    if os.path.isdir(path):
        return times_from_inf_dir(path)
    if path.upper().endswith(".INF"):
        m = INF_RE.search(os.path.basename(path))
        return pd.DataFrame([dict(tray_code=int(m.group(2)), round=int(m.group(3)), datetime=read_inf_time(path))])
    # PS2 export (.TXT with instrument header, or header-stripped .csv: File, Date, Time)
    # or a plain CSV with columns tray_code, round, datetime
    with open(path, encoding="utf-8-sig", errors="replace") as fh:
        lines = fh.readlines(4_000_000)
    hdr = next((i for i, l in enumerate(lines) if l.lstrip("\ufeff").startswith(("File", "tray_code"))), 0)
    sep = "\t" if "\t" in lines[hdr] else ","
    if lines[hdr].startswith("tray_code"):
        T = pd.read_csv(path, sep=sep, skiprows=hdr)
        T["datetime"] = pd.to_datetime(T["datetime"])
        return T[["tray_code", "round", "datetime"]]
    T = pd.read_csv(path, sep=sep, skiprows=hdr, usecols=[0, 1, 2], dtype=str, encoding="utf-8-sig")
    T.columns = ["File", "Date", "Time"]
    T = T[T["File"].str.match(r"^[A-Z]+_\d+_\d+_\d+", na=False)].drop_duplicates("File")
    parts = T["File"].str.replace(r"\.\w+$", "", regex=True).str.split("_")
    T["tray_code"] = parts.str[2].astype(int); T["round"] = parts.str[3].astype(int)
    T["datetime"] = pd.to_datetime(T["Date"] + " " + T["Time"], format="%Y%m%d %H:%M")
    return T.groupby(["tray_code", "round"], as_index=False)["datetime"].min()

# ----------------------------------------------------------------------------- driver
def _job(args):
    tray_key, paths, times, params, ov, ckpt, deadline = args
    try:
        return tray_key, process_tray(paths, times, params, ov, ckpt, deadline), None
    except Exception as e:                              # keep going; report at the end
        import traceback
        return tray_key, None, f"{type(e).__name__}: {e}\n{traceback.format_exc()}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("raw_dir"); ap.add_argument("out_dir")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--lo", type=float, default=DEFAULTS["lo"], help="base low (edge) threshold on q")
    ap.add_argument("--hi", type=float, default=DEFAULTS["hi"], help="base core threshold on q (raised per pot for algae)")
    ap.add_argument("--hi-max", type=float, default=DEFAULTS["hi_max"], help="ceiling for the per-pot core threshold")
    ap.add_argument("--establish-days", type=float, default=DEFAULTS["establish_days"],
                    help="window (days from the tray's first image) used to call empty pots")
    ap.add_argument("--establish-n", type=int, default=DEFAULTS["establish_n"],
                    help="same window in images, used only when imaging times are unknown")
    ap.add_argument("--establish-k", type=int, default=DEFAULTS["establish_k"],
                    help="a pot holds a plant if detected in >= K images within that window")
    ap.add_argument("--no-tracking", action="store_true",
                    help="v0.2 behaviour: independent images, fixed thresholds, no empty-pot call")
    ap.add_argument("--pitch-mm", type=float, default=63.5,
                    help="pot pitch in mm (default 63.5 = top width of the standard tray pots); area_mm2 = area_cellfrac * pitch_mm^2")
    ap.add_argument("--first-code", type=int, default=None,
                    help="tray code of Tray01 (default: smallest code seen)")
    ap.add_argument("--code-step", type=int, default=20, help="increment between consecutive tray codes")
    ap.add_argument("--times", nargs="*", default=[],
                    help="extra time sources: folder(s) of HDR_*.INF, PS2 export(s), or a "
                         "tray_code,round,datetime CSV. HDR_*.INF files in RAW_DIR are always used.")
    ap.add_argument("--no-overlays", action="store_true")
    ap.add_argument("--no-traits", action="store_true",
                    help="skip the PS2 / spectral / colour traits from the PMT_, SPC_ and CLR_ files")
    ap.add_argument("--trays", nargs="*", type=int, default=None, help="only these tray codes (testing)")
    ap.add_argument("--max-minutes", type=float, default=None,
                    help="stop starting new images after this long; rerun to continue (progress is "
                         "checkpointed per image)")
    a = ap.parse_args(argv)
    deadline = time.time() + 60 * a.max_minutes if a.max_minutes else None
    params = dict(lo=a.lo, hi=a.hi, hi_max=a.hi_max, establish_n=a.establish_n,
                  establish_days=a.establish_days,
                  establish_k=a.establish_k, tracking=not a.no_tracking, traits=not a.no_traits)

    os.makedirs(a.out_dir, exist_ok=True)
    ov = None if a.no_overlays else os.path.join(a.out_dir, "overlays")
    if ov: os.makedirs(ov, exist_ok=True)
    f_img, f_pot = os.path.join(a.out_dir, "images.csv"), os.path.join(a.out_dir, "pots.csv")
    ckdir = os.path.join(a.out_dir, "checkpoints"); os.makedirs(ckdir, exist_ok=True)

    # imaging times first: tracking needs the true order and spacing of images
    T = pd.concat([read_times(t) for t in [a.raw_dir] + list(a.times)])
    T = T.dropna().groupby(["tray_code", "round"], as_index=False)["datetime"].min()
    tmap = {(r.tray_code, r.round): r.datetime for r in T.itertuples()}

    trays = {}
    for f in sorted(os.listdir(a.raw_dir)):
        m = parse_name(f)
        if m and (a.trays is None or m["tray_code"] in a.trays):
            trays.setdefault((m["exp"], m["tray_code"]), []).append(os.path.join(a.raw_dir, f))
    print(f"{sum(len(v) for v in trays.values())} APH files in {len(trays)} trays", flush=True)

    # Every tray is handed to process_tray with its checkpoint: finished trays return at once,
    # partly done trays resume, and a tray whose image list changed (new images earlier in
    # time, or a new pipeline version / parameters) is redone from the start.
    imgs, pots, errs, incomplete = [], [], [], []
    jobs = []
    for (exp, code), paths in trays.items():
        times = {os.path.basename(f): tmap.get((code, parse_name(f)["round"])) for f in paths}
        jobs.append(((exp, code), paths, times, params, ov, os.path.join(ckdir, f"{exp}_{code}.pkl"), deadline))
    with ProcessPoolExecutor(a.workers) as ex:
        futs = [ex.submit(_job, j) for j in jobs]
        for i, fu in enumerate(as_completed(futs), 1):
            key, res, err = fu.result()
            if err:
                errs.append((key, err)); continue
            imgs.extend(res[0]); pots.extend(res[1])
            if not res[2]:
                incomplete.append(key[1])
            print(f"  tray {key[1]}: {len(res[0])}/{len(trays[key])} images ({i}/{len(jobs)})", flush=True)

    if imgs:
        I = pd.DataFrame(imgs); P = pd.DataFrame(pots)
        first = a.first_code if a.first_code is not None else int(I["tray_code"].min())
        for D in (I, P):     # tray code -> Tray01.. by arithmetic, not by rank
            n = (D["tray_code"] - first) / a.code_step + 1
            D["TrayID"] = ["Tray%02d" % int(x) if float(x).is_integer() and x >= 1 else "Tray??" for x in n]
        I = I.drop(columns=["datetime", "TOE_h"], errors="ignore").merge(T, on=["tray_code", "round"], how="left")
        I["datetime"] = pd.to_datetime(I["datetime"])
        I["TOE_h"] = (I["datetime"] - I["datetime"].min()).dt.total_seconds() / 3600
        P = P.drop(columns=["datetime", "TOE_h"], errors="ignore").merge(
            I[["file", "datetime", "TOE_h"]], on="file", how="left")
        miss = int(I["datetime"].isna().sum())
        if miss:
            print(f"WARNING: {miss} of {len(I)} images have no imaging time "
                  f"(no matching HDR_*.INF in RAW_DIR or --times); rounds assumed 12 h apart", file=sys.stderr)
        if params["tracking"]:
            P = call_pot_status(P, dict(DEFAULTS, **params))
        P = P.drop(columns=["area_mm2", "area_mm2_raw"], errors="ignore")
        if a.pitch_mm:
            P["area_mm2"] = P["area_cellfrac"] * a.pitch_mm ** 2
            if "area_cellfrac_raw" in P:
                P["area_mm2_raw"] = P["area_cellfrac_raw"] * a.pitch_mm ** 2
        med = I.groupby(["exp", "tray_code"])["wall_score"].transform("median")
        I["qc_grid_weak"] = I["wall_score"] < 0.6 * med
        I = I.sort_values(["exp", "tray_code", "round"]); P = P.sort_values(["exp", "tray_code", "round", "Roi_No"])
        I.to_csv(f_img, index=False); P.to_csv(f_pot, index=False)
        if "pot_status" in P:
            E = P[P.pot_status == "empty"].groupby(["TrayID", "Area"])["late_signal_cellfrac"].first()
            print(f"empty pots: {len(E)}" + ("".join(f"\n  {t} {ar}  (largest later signal {x * a.pitch_mm**2:.0f} mm2)"
                                                     for (t, ar), x in E.items()) if len(E) else ""))
    for k, e in errs:
        print("ERROR tray", k, e, file=sys.stderr)
    if incomplete:
        print(f"NOT FINISHED: {len(incomplete)} trays stopped at the time limit; run the same command "
              f"again to continue. Tables written so far cover the images done.")
    print(f"{len(imgs)} images in tables, {len(errs)} tray errors -> {a.out_dir}")


if __name__ == "__main__":
    main()
