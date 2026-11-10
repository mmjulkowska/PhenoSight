#!/usr/bin/env python3
"""
phenosight_raw.py — open PhenoVation CropReporter / PhenoSight raw exports
(HDR_*.INF + APH/CHL/CLR/PMT/SPC_*.DAT) without the vendor software, build a
plant mask from chlorophyll fluorescence, split it into user-defined ROIs
(one per pot), and extract the PhenoVation trait set per ROI.

PhenoSight Plant Phenotyping Facility, Boyce Thompson Institute.  Python >= 3.9; numpy, scipy, scikit-image, pandas,
pillow, matplotlib, tifffile.

------------------------------------------------------------------------
Raw format (verified on CropReporter 7.3.1 files, Sep 2026):
  * every .DAT is headerless little-endian uint16
  * frames are stored as (n_frames, ImageCols, ImageRows), i.e. column-major;
    transpose to (n_frames, ImageRows, ImageCols) = 1500 x 2048 landscape
    (same convention as PlantCV's read_cropreporter)
  * reflectance bands are pre-calibrated by the shutter settings:
    reflectance = DN / 65535  (leaf NIR ~0.50, red ~0.06, NDVI ~0.77)
  * HDR TimeStamp is a Delphi/Excel serial date with a decimal comma
    (days since 1899-12-30, local instrument time)

Frame maps (SaveAllFrames = 0).  PMT = the "TmPam" protocol:
  0 Fdark   1 F0   2 Fm   3 Fdarksat   4 Flight   5 Fs'   6 Fm'   7 Flightsat   8 F0'
  Evidence: background levels (raw ambient frames 0/3/4/7 carry the 240 DN
  offset + light leak: dark 240, sat-only 288, actinic 720, actinic+sat 752;
  fluorescence frames 1/2/5/6/8 are ambient-subtracted, bg ~20-50) and an
  independent NPQ check: (Flightsat-Flight)/(Fdarksat-Fdark) = Fm'/Fm
  (0.68) agrees with 1/(1+NPQ) from frames 2 and 6 (0.65).
  CLR = Red, Green, Blue      SPC = Anth(~540), RedEdge(~710), NIR(~770)
  APH = Red, FarRed (alpha/absorptance)      CHL = Fdark, Chl
------------------------------------------------------------------------
Usage
  python3 phenosight_raw.py <data_folder> [--rois rois.csv] [--out DIR]
                            [--tiff] [--mask-level N] [--n-pots 5]
"""
from __future__ import annotations

import argparse
import datetime as dt
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import measure, morphology, segmentation
from skimage.draw import disk as draw_disk, polygon as draw_polygon

DN_MAX = 65535.0

FRAME_MAPS = {
    "PMT": ["Fdark", "F0", "Fm", "Fdarksat", "Flight", "Fs_p", "Fm_p", "Flightsat", "F0_p"],
    "CLR": ["Red", "Green", "Blue"],
    "SPC": ["Anth", "RedEdge", "NIR"],
    "APH": ["AphRed", "AphFarRed"],
    "CHL": ["ChlDark", "Chl"],
    # vendor protocols not used at PhenoSight so far, known from PlantCV:
    "PMD": ["Fdark", "F0", "Fm", "Fdarksat"],
    "PML": ["Flight", "Fs_p", "Fm_p", "Flightsat"],
}

TRAITS_PS2 = ["Fv.Fm", "Fq.Fm", "NPQ", "qP", "qN", "qL", "phiNO", "phiNPQ", "npq.t"]
TRAITS_SPEC = ["ChlIdx", "AriIdx", "NDVI", "Abs"]
TRAITS_RGB = ["R", "G", "B", "Hue"]


# ---------------------------------------------------------------- reading
def read_inf(path: str) -> dict:
    meta = {}
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if "=" in line:
                k, v = line.rstrip("\r\n").split("=", 1)
                meta[k.strip()] = v.strip()
    return meta


def delphi_time(s: str) -> dt.datetime:
    days = float(s.replace(",", "."))
    return dt.datetime(1899, 12, 30) + dt.timedelta(days=days)


def read_dat(path: str, rows: int, cols: int) -> np.ndarray:
    raw = np.fromfile(path, dtype="<u2")
    n = raw.size // (rows * cols)
    if n * rows * cols != raw.size:
        raise ValueError(f"{path}: size {raw.size} not a multiple of {rows}x{cols}")
    return raw.reshape(n, cols, rows).transpose(0, 2, 1)  # -> (n, rows, cols)


def load_measurement(inf_path: str) -> dict:
    """Return dict with meta, timestamp and named float32 frames."""
    meta = read_inf(inf_path)
    rows, cols = int(meta["ImageRows"]), int(meta["ImageCols"])
    folder, name = os.path.split(inf_path)
    frames, kinds = {}, {}
    for kind, labels in FRAME_MAPS.items():
        p = os.path.join(folder, name.replace("HDR_", f"{kind}_", 1).rsplit(".", 1)[0] + ".DAT")
        if not os.path.exists(p):
            continue
        cube = read_dat(p, rows, cols)
        if len(labels) != cube.shape[0]:
            print(f"  ! {os.path.basename(p)} has {cube.shape[0]} frames, expected "
                  f"{len(labels)} — labelled generically", file=sys.stderr)
            labels = [f"{kind}{i}" for i in range(cube.shape[0])]
        kinds[kind] = labels
        for lab, fr in zip(labels, cube):
            frames[lab] = fr.astype(np.float32)
    m = re.search(r"E(\d+)P(\d+)N(\d+)", name)
    return dict(
        inf=inf_path, stem=name.rsplit(".", 1)[0].replace("HDR_", "", 1), meta=meta,
        time=delphi_time(meta["TimeStamp"]) if "TimeStamp" in meta else None,
        experiment=int(m.group(1)) if m else None,
        plant_no=int(m.group(2)) if m else None,
        image_no=int(m.group(3)) if m else None,
        frames=frames, kinds=kinds, shape=(rows, cols),
    )


# ---------------------------------------------------------------- mask
def fluorescence_mask(fr: dict, meta: dict, level: float | None = None,
                      min_object_px: int = 30) -> np.ndarray:
    """Plant mask = Fm (dark-adapted maximal fluorescence) above the mask level
    stored in the INF (TmPamMaskLevel; PhenoVation's own threshold), cleaned."""
    if level is None:
        level = float(meta.get("TmPamMaskLevel", meta.get("DkPamMaskLevel", 2000)))
    base = fr["Fm"] if "Fm" in fr else fr["Fm_p"]
    m = base > level
    m = morphology.remove_small_objects(m, min_object_px)
    m = morphology.remove_small_holes(m, 200)   # fill specks, keep gaps between leaflets
    return m


# ---------------------------------------------------------------- ROIs
def read_rois(path: str) -> pd.DataFrame:
    """ROI file (written by roi_editor.html, or by hand):
         image, roi, label, shape, x, y, r
       image  = HDR stem (e.g. E0001P0007N0001_GCB22060646_20260923), a PlantNo
                (e.g. P0007) or 'default' (used for any image without its own ROIs)
       roi    = 1..n (pot number)
       label  = free text (treatment / genotype); joined to traits as-is
       shape  = square | rect | polygon : one row per vertex, in drawing order
                circle                  : ONE row, x/y = centre, r = radius
                (column optional; missing = polygon)
       x, y   = pixel coordinates in the 2048 x 1500 landscape image"""
    r = pd.read_csv(path)
    r.columns = [c.strip().lower() for c in r.columns]
    need = {"image", "roi", "x", "y"}
    if not need.issubset(r.columns):
        raise ValueError(f"ROI file needs columns {need}, has {list(r.columns)}")
    if "label" not in r.columns:
        r["label"] = ""
    r["label"] = r["label"].fillna("").astype(str)
    if "shape" not in r.columns:
        r["shape"] = "polygon"
    r["shape"] = r["shape"].fillna("polygon").astype(str).str.strip().str.lower().replace("", "polygon")
    if "r" not in r.columns:
        r["r"] = np.nan
    bad = r[(r["shape"] == "circle") & ~(pd.to_numeric(r["r"], errors="coerce") > 0)]
    if len(bad):
        raise ValueError(f"circle ROI(s) without a positive radius r: {bad[['image','roi']].values.tolist()}")
    r["image"] = r["image"].astype(str)
    return r


def rois_for_image(rois: pd.DataFrame, m: dict) -> pd.DataFrame | None:
    keys = [m["stem"], f"P{m['plant_no']:04d}" if m["plant_no"] is not None else None, "default"]
    for k in keys:
        if k is None:
            continue
        sub = rois[rois["image"] == k]
        if len(sub):
            return sub
    return None


def labels_from_rois(mask: np.ndarray, roi_df: pd.DataFrame,
                     assign: str = "inside") -> tuple[np.ndarray, dict]:
    """Rasterise polygons -> label image.  assign='inside': only mask pixels
    inside a polygon count.  assign='grow': polygons act as seeds and every
    mask pixel connected to them is assigned to the geodesically nearest ROI
    (useful when leaves sprawl outside the drawn box)."""
    shape = mask.shape
    roi_img = np.zeros(shape, np.int32)
    info = {}
    for rid, g in roi_df.groupby("roi", sort=True):
        if g["shape"].iloc[0] == "circle":
            cx, cy, rad = float(g["x"].iloc[0]), float(g["y"].iloc[0]), float(g["r"].iloc[0])
            rr, cc = draw_disk((cy, cx), rad, shape=shape)
        else:
            rr, cc = draw_polygon(g["y"].to_numpy(float), g["x"].to_numpy(float), shape)
        free = roi_img[rr, cc] == 0          # first ROI wins on overlap
        roi_img[rr[free], cc[free]] = int(rid)
        info[int(rid)] = str(g["label"].iloc[0]) if "label" in g else ""
    if assign == "grow":
        reach = morphology.binary_dilation(mask, morphology.disk(3)) | (roi_img > 0)
        lab = segmentation.watershed(np.zeros(shape, np.uint8), markers=roi_img, mask=reach)
        return lab * mask, info
    return roi_img * mask, info


def auto_pot_rois(fr: dict, mask: np.ndarray, n_pots: int = 5,
                  min_sep: int = 250) -> tuple[np.ndarray, dict, list]:
    """Fallback when no ROI file: find pots from the brown root-tube/soil
    colour in the RGB frames, then grow plants from them geodesically."""
    R, G, B = fr["Red"], fr["Green"], fr["Blue"]
    brown = (R > 1.6 * B) & (R > 1.15 * G) & (R > 3000) & ~mask
    h, w = brown.shape   # ignore the cabinet edges (orange cables / hardware)
    brown[: int(0.12 * h)] = brown[int(0.88 * h):] = False
    brown[:, : int(0.05 * w)] = brown[:, int(0.95 * w):] = False
    dens = ndi.uniform_filter(brown.astype(np.float32), 121)
    pk = []
    d = dens.copy()
    for _ in range(n_pots):
        y, x = np.unravel_index(np.argmax(d), d.shape)
        if d[y, x] <= 0.01:
            break
        pk.append((y, x))
        yy, xx = np.ogrid[:d.shape[0], :d.shape[1]]
        d[(yy - y) ** 2 + (xx - x) ** 2 < min_sep ** 2] = 0
    # number pots in reading order: top row left->right, then bottom row
    pk.sort(key=lambda p: p[0])
    row, rows_ = 0, []
    for i, p in enumerate(pk):
        if i and p[0] - pk[i - 1][0] > 200:
            row += 1
        rows_.append(row)
    pk = [p for _, p in sorted(zip(rows_, pk), key=lambda t: (t[0], t[1][1]))]
    markers = np.zeros(mask.shape, np.int32)
    for i, (y, x) in enumerate(pk, 1):
        yy, xx = np.ogrid[:mask.shape[0], :mask.shape[1]]
        markers[(yy - y) ** 2 + (xx - x) ** 2 < 70 ** 2] = i
    reach = morphology.binary_dilation(mask, morphology.disk(3)) | (markers > 0)
    lab = segmentation.watershed(np.zeros(mask.shape, np.uint8), markers=markers, mask=reach)
    lab = lab * mask
    # pieces not connected to any pot: nearest pot by Euclidean distance
    orphan = mask & (lab == 0)
    if orphan.any() and pk:
        _, (iy, ix) = ndi.distance_transform_edt(markers == 0, return_indices=True)
        lab[orphan] = markers[iy[orphan], ix[orphan]]
    return lab, {i: f"pot{i}" for i in range(1, len(pk) + 1)}, pk


# ---------------------------------------------------------------- tissues
def classify_tissue(mask: np.ndarray, stem_radius: int = 4, min_stem_px: int = 40,
                    min_elongation: float = 3.0, rim_px: int = 2):
    """Split the plant mask into leaf blade, leaf rim and stem/petiole.
    Stems are found by SHAPE, not by NDVI: a morphological opening with a disk of
    radius `stem_radius` keeps the leaf blades and drops every structure narrower
    than ~2*stem_radius px; of the dropped pieces, those that are big enough and
    elongated (major/minor axis >= min_elongation) are stems/petioles, the rest
    (leaf tips, serrations) go back to the leaf.  NDVI alone cannot do this,
    because the 1-2 px rim of every leaf is a mixed leaf/background pixel with
    the same low NDVI as a stem.
    Returns boolean images (leaf, stem, rim); rim is a subset of leaf."""
    blade = morphology.binary_opening(mask, morphology.disk(stem_radius))
    thin = mask & ~blade
    lab = measure.label(thin)
    stem = np.zeros_like(mask)
    for p in measure.regionprops(lab):
        if p.area >= min_stem_px and p.major_axis_length / max(p.minor_axis_length, 1) >= min_elongation:
            stem[lab == p.label] = True
    leaf = mask & ~stem
    rim = leaf & ~morphology.binary_erosion(blade, morphology.disk(rim_px))
    return leaf, stem, rim


# ---------------------------------------------------------------- traits
def _div(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        out = a / b
    out[~np.isfinite(out)] = np.nan
    return out


def trait_maps(fr: dict) -> dict:
    """Per-pixel trait images, PhenoVation formulas
    (https://www.phenovation.com/formulas).  Frames that are absent give NaN."""
    t = {}
    if {"F0", "Fm"} <= fr.keys():
        F0, Fm = fr["F0"], fr["Fm"]
        t["Fv.Fm"] = _div(Fm - F0, Fm)
    if {"Fs_p", "Fm_p"} <= fr.keys():
        Fs, Fmp = fr["Fs_p"], fr["Fm_p"]
        t["Fq.Fm"] = _div(Fmp - Fs, Fmp)
        if "Fm" in fr:
            t["NPQ"] = _div(fr["Fm"] - Fmp, Fmp)
        if "F0_p" in fr:
            F0p = fr["F0_p"]
            t["qP"] = _div(Fmp - Fs, Fmp - F0p)
            t["qL"] = t["qP"] * _div(F0p, Fs)
            t["npq.t"] = _div(np.full_like(Fmp, 4.88), _div(Fmp, F0p) - 1) - 1
            if {"F0", "Fm"} <= fr.keys():
                t["qN"] = 1 - _div(Fmp - F0p, fr["Fm"] - fr["F0"])
        if {"NPQ", "qL"} <= t.keys():
            t["phiNO"] = _div(np.ones_like(Fmp), t["NPQ"] + 1 + t["qL"] * _div(fr["Fm"], fr["F0"]))
            t["phiNPQ"] = 1 - t["Fq.Fm"] - t["phiNO"]
    # qE / qI need Fm'' (second dark-adapted saturation) — not saved when
    # SaveAllFrames=0, so they are reported as NaN rather than guessed.
    if {"Anth", "RedEdge", "NIR"} <= fr.keys():
        Ra, Rre, Rn = fr["Anth"] / DN_MAX, fr["RedEdge"] / DN_MAX, fr["NIR"] / DN_MAX
        t["ChlIdx"] = (_div(np.ones_like(Rre), Rre) - _div(np.ones_like(Rn), Rn)) * Rn
        t["AriIdx"] = (_div(np.ones_like(Ra), Ra) - _div(np.ones_like(Rre), Rre)) * Rn
        if "Red" in fr:
            Rr = fr["Red"] / DN_MAX
            t["NDVI"] = _div(Rn - Rr, Rn + Rr)
    if {"AphRed", "AphFarRed"} <= fr.keys():
        # absorptance as PhenoVation defines it; APH bands use their own
        # shutters, so treat Abs as relative until checked against an export
        t["Abs"] = 1 - _div(fr["AphRed"], fr["AphFarRed"])
    if {"Red", "Green", "Blue"} <= fr.keys():
        r, g, b = (fr[k] / DN_MAX for k in ("Red", "Green", "Blue"))
        t["R"], t["G"], t["B"] = r, g, b
        mx, mn = np.maximum(np.maximum(r, g), b), np.minimum(np.minimum(r, g), b)
        c = mx - mn
        with np.errstate(divide="ignore", invalid="ignore"):
            h = np.where(mx == r, ((g - b) / c) % 6, np.where(mx == g, (b - r) / c + 2, (r - g) / c + 4))
        t["Hue"] = np.where(c > 0, h * 60.0, np.nan).astype(np.float32)
    return t


# plausibility gates applied per pixel before averaging (keeps stems/edge
# pixels with near-zero signal from dragging the means)
GATES = {"Fv.Fm": (0, 1), "Fq.Fm": (0, 1), "NPQ": (0, 10), "qP": (0, 1.2), "qN": (0, 1.2),
         "qL": (0, 1.5), "phiNO": (0, 1), "phiNPQ": (-0.2, 1), "npq.t": (-1, 20),
         "ChlIdx": (-1, 20), "AriIdx": (-5, 20), "NDVI": (-1, 1), "Abs": (-1, 1)}


def summarise(labels: np.ndarray, roi_info: dict, maps: dict, m: dict,
              px_per_mm: float | None = None, include_all: bool = True) -> pd.DataFrame:
    rows = []
    ids = list(roi_info.keys())
    if include_all:
        ids = ["All"] + ids
    for rid in ids:
        sel = labels > 0 if rid == "All" else labels == rid
        n = int(sel.sum())
        rec = dict(File=os.path.basename(m["inf"]), Stem=m["stem"],
                   Date=m["time"].strftime("%Y-%m-%d") if m["time"] else None,
                   Time=m["time"].strftime("%H:%M:%S") if m["time"] else None,
                   Experiment=m["experiment"], PlantNo=m["plant_no"], ImageNo=m["image_no"],
                   Roi=rid, Label="All" if rid == "All" else roi_info[rid],
                   **{"Area.px": n})
        if px_per_mm:
            rec["Area.mm2"] = n / px_per_mm ** 2
        if n and rid != "All":
            props = measure.regionprops(sel.astype(np.uint8))[0]
            rec.update({"Perimeter.px": props.perimeter, "ConvexArea.px": props.area_convex,
                        "Solidity": props.solidity,
                        "Height.px": props.bbox[2] - props.bbox[0],
                        "Width.px": props.bbox[3] - props.bbox[1],
                        "Centroid.x": props.centroid[1], "Centroid.y": props.centroid[0],
                        "nPieces": int(measure.label(sel).max())})
        for k, img in maps.items():
            v = img[sel]
            if k in GATES:
                lo, hi = GATES[k]
                v = v[(v >= lo) & (v <= hi)]
            v = v[np.isfinite(v)]
            rec[k] = float(v.mean()) if v.size else np.nan
            if k in ("Fv.Fm", "Fq.Fm", "NPQ"):
                rec[k + ".median"] = float(np.median(v)) if v.size else np.nan
        rec["qE"] = np.nan
        rec["qI"] = np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- outputs
def rgb_preview(fr: dict, gamma: float = 0.6) -> np.ndarray:
    rgb = np.stack([fr["Red"], fr["Green"], fr["Blue"]], -1) / DN_MAX
    rgb = np.clip(rgb / np.percentile(rgb, 99.5), 0, 1) ** gamma
    return (rgb * 255).astype(np.uint8)


def save_qc(m: dict, mask, labels, roi_info, maps, out_png, roi_df=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    rgb = rgb_preview(m["frames"])
    fig, ax = plt.subplots(2, 3, figsize=(18, 9.5), dpi=90)
    ax = ax.ravel()
    ax[0].imshow(rgb); ax[0].set_title("RGB + ROIs")
    if roi_df is not None:
        for rid, g in roi_df.groupby("roi"):
            if g["shape"].iloc[0] == "circle":
                ax[0].add_patch(plt.Circle((g["x"].iloc[0], g["y"].iloc[0]), g["r"].iloc[0],
                                           fill=False, lw=1.5, color="yellow"))
                continue
            xs, ys = list(g["x"]) + [g["x"].iloc[0]], list(g["y"]) + [g["y"].iloc[0]]
            ax[0].plot(xs, ys, "-", lw=1.5, color="yellow")
    for rid, lab in roi_info.items():
        yy, xx = np.nonzero(labels == rid)
        if yy.size:
            ax[0].text(xx.mean(), yy.mean(), f"{rid}\n{lab}", color="white", ha="center",
                       fontsize=11, weight="bold",
                       bbox=dict(fc="black", alpha=0.5, lw=0))
    cols = plt.get_cmap("tab10")(np.arange(10))
    lab_rgb = np.zeros(labels.shape + (3,))
    for rid in roi_info:
        lab_rgb[labels == rid] = cols[(rid - 1) % 10, :3]
    unassigned = mask & (labels == 0)
    lab_rgb[unassigned] = (0.5, 0.5, 0.5)
    ax[1].imshow(lab_rgb); ax[1].set_title(f"Fluorescence mask by ROI (grey = outside ROIs: {unassigned.sum()} px)")
    for a, (k, lo, hi, cm) in zip(ax[2:], [("Fv.Fm", 0.6, 0.85, "viridis"), ("Fq.Fm", 0, 0.3, "viridis"),
                                           ("NPQ", 0, 1.5, "magma"), ("NDVI", 0.5, 0.9, "RdYlGn")]):
        if k in maps:
            img = np.where(mask, maps[k], np.nan)
            im = a.imshow(img, vmin=lo, vmax=hi, cmap=cm)
            a.set_facecolor("black"); fig.colorbar(im, ax=a, fraction=0.035)
        a.set_title(k)
    for a in ax:
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle(f"{m['stem']}   {m['time']}", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_png); plt.close(fig)


def export_tiffs(m: dict, out_dir: str):
    import tifffile
    os.makedirs(out_dir, exist_ok=True)
    for kind, labels in m["kinds"].items():   # only files that exist on disk
        stack = [m["frames"][l] for l in labels]
        if stack:
            tifffile.imwrite(os.path.join(out_dir, f"{kind}_{m['stem']}.tif"),
                             np.stack(stack).astype(np.uint16), imagej=True,
                             metadata={"Labels": labels, "axes": "ZYX"})


# fixed display ranges for the colour trait images, so images from different days
# and trays are directly comparable (edit here if a trait runs off-scale)
COLOR_RANGES = {
    "Fv.Fm": (0.60, 0.85, "viridis"), "Fq.Fm": (0.0, 0.30, "viridis"),
    "NPQ": (0.0, 1.5, "magma"), "qP": (0.0, 0.4, "viridis"), "qN": (0.0, 0.8, "magma"),
    "qL": (0.0, 0.2, "viridis"), "phiNO": (0.4, 0.9, "cividis"), "phiNPQ": (0.0, 0.5, "magma"),
    "npq.t": (0.0, 3.0, "magma"), "ChlIdx": (0.0, 4.0, "YlGn"), "AriIdx": (0.0, 2.0, "RdPu"),
    "NDVI": (0.5, 0.9, "RdYlGn"), "Abs": (0.3, 0.8, "viridis"),
}


def export_trait_tiffs(m: dict, mask, labels, maps: dict, out_dir: str, tissue=None):
    """Per image: <out_dir>/<stem>/
         <trait>.png        colour image on the fixed scale in COLOR_RANGES, black background
                            (opens in Preview / PowerPoint / Illustrator)
         mask.png           white = plant (Fm above mask level), black = background
         roi_labels.png     plant pixels coloured by pot, pot numbers printed on top
         tissue.png         leaf blade (green), leaf rim (dark green), stem/petiole (orange)
         color_scales.png   the colour bar and range used for every trait
         values_for_Fiji/   the numbers behind the pictures, for measuring:
             <trait>.tif      32-bit float trait value per plant pixel, NaN outside the mask
             mask.tif         8-bit, 255 = plant
             roi_labels.tif   8-bit, pixel value = pot number (0 = outside all ROIs)
             tissue.tif       8-bit, 1 = leaf core, 2 = leaf rim, 3 = stem/petiole
    All 2048 x 1500, same pixel grid as the RGB preview and the ROI coordinates."""
    import tifffile
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image, ImageDraw, ImageFont
    d = os.path.join(out_dir, m["stem"])
    v = os.path.join(d, "values_for_Fiji")
    os.makedirs(v, exist_ok=True)
    try:                       # LZW opens in Fiji and macOS Preview; needs imagecodecs
        import imagecodecs     # noqa: F401
        kw = dict(compression="lzw")
    except ImportError:
        kw = {}
    tifffile.imwrite(os.path.join(v, "mask.tif"), (mask * 255).astype(np.uint8), **kw)
    tifffile.imwrite(os.path.join(v, "roi_labels.tif"), labels.astype(np.uint8), **kw)
    Image.fromarray((mask * 255).astype(np.uint8)).save(os.path.join(d, "mask.png"))
    if tissue is not None:
        leaf, stem, rim = tissue
        code = np.zeros(mask.shape, np.uint8)
        code[leaf] = 1; code[rim] = 2; code[stem] = 3
        tifffile.imwrite(os.path.join(v, "tissue.tif"), code, **kw)
        pal = np.array([(0, 0, 0), (90, 180, 90), (35, 100, 35), (255, 140, 0)], np.uint8)
        Image.fromarray(pal[code]).save(os.path.join(d, "tissue.png"))

    # coloured ROI view with pot numbers
    cols = (plt.get_cmap("tab10")(np.arange(10))[:, :3] * 255).astype(np.uint8)
    lab_rgb = np.zeros(labels.shape + (3,), np.uint8)
    lab_rgb[mask & (labels == 0)] = (110, 110, 110)
    ids = [int(i) for i in np.unique(labels) if i > 0]
    for i in ids:
        lab_rgb[labels == i] = cols[(i - 1) % 10]
    im = Image.fromarray(lab_rgb)
    dr = ImageDraw.Draw(im)
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", 56)
    except OSError:
        try:
            font = ImageFont.truetype("Arial Bold.ttf", 56)
        except OSError:
            font = ImageFont.load_default()
    for i in ids:
        yy, xx = np.nonzero(labels == i)
        dr.text((xx.mean(), yy.mean()), str(i), fill=(255, 255, 255), font=font, anchor="mm",
                stroke_width=4, stroke_fill=(0, 0, 0))
    im.save(os.path.join(d, "roi_labels.png"))

    shown = []
    for k, img in maps.items():
        if k in ("R", "G", "B", "Hue"):
            continue
        fname = k.replace(".", "_")
        val = np.where(mask, img, np.nan).astype(np.float32)
        tifffile.imwrite(os.path.join(v, f"{fname}.tif"), val, **kw)
        if k in COLOR_RANGES:
            lo, hi, cm = COLOR_RANGES[k]
            norm = np.clip((val - lo) / (hi - lo), 0, 1)
            rgb = (matplotlib.colormaps[cm](np.nan_to_num(norm))[..., :3] * 255).astype(np.uint8)
            rgb[~mask | ~np.isfinite(val)] = 0
            Image.fromarray(rgb).save(os.path.join(d, f"{fname}.png"))
            shown.append((k, lo, hi, cm))

    # one sheet with every colour scale
    fig, axes = plt.subplots(len(shown), 1, figsize=(6, 0.75 * len(shown) + 0.4), dpi=150)
    for ax, (k, lo, hi, cm) in zip(np.atleast_1d(axes), shown):
        ax.imshow(np.linspace(0, 1, 256)[None, :], aspect="auto", cmap=cm, extent=(lo, hi, 0, 1))
        ax.set_yticks([]); ax.set_ylabel(k, rotation=0, ha="right", va="center", fontsize=9)
        ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(d, "color_scales.png")); plt.close(fig)


# ---------------------------------------------------------------- driver
def process_folder(folder, out=None, rois_path=None, mask_level=None, n_pots=5,
                   assign="inside", tiff=False, px_per_mm=None, previews_only=False,
                   trait_tiffs=True, tissues=True, stem_radius=4):
    out = out or os.path.join(folder, "open_pipeline")
    os.makedirs(os.path.join(out, "qc"), exist_ok=True)
    os.makedirs(os.path.join(out, "rgb_for_roi_picking"), exist_ok=True)
    rois = read_rois(rois_path) if rois_path else None
    infs = sorted(glob.glob(os.path.join(folder, "HDR_*.INF")))
    if not infs:
        sys.exit(f"no HDR_*.INF files in {folder}")
    tables = []
    for inf in infs:
        m = load_measurement(inf)
        print(f"{m['stem']}: {m['time']}  frames={list(m['frames'])}")
        from PIL import Image
        Image.fromarray(rgb_preview(m["frames"])).save(
            os.path.join(out, "rgb_for_roi_picking", f"{m['stem']}_RGB.png"))
        if previews_only:
            continue
        mask = fluorescence_mask(m["frames"], m["meta"], mask_level)
        roi_df = rois_for_image(rois, m) if rois is not None else None
        if roi_df is not None:
            labels, info = labels_from_rois(mask, roi_df, assign)
            mode = f"ROI file ({assign})"
        else:
            labels, info, _ = auto_pot_rois(m["frames"], mask, n_pots)
            mode = "auto pot detection"
        maps = trait_maps(m["frames"])
        if tissues:
            leaf, stem, rim = classify_tissue(mask, stem_radius)
            parts = [("whole", mask), ("leaf", leaf), ("leaf_core", leaf & ~rim), ("stem", stem)]
        else:
            leaf = stem = rim = None
            parts = [("whole", mask)]
        dfs = []
        for tis, sel in parts:
            d_ = summarise(labels * sel, info, maps, m, px_per_mm)
            d_.insert(d_.columns.get_loc("Label") + 1, "Tissue", tis)
            dfs.append(d_)
        df = pd.concat(dfs, ignore_index=True)
        df["RoiMode"] = mode
        df["MaskPx.total"] = int(mask.sum())
        df["MaskPx.outsideROIs"] = int((mask & (labels == 0)).sum())
        tables.append(df)
        save_qc(m, mask, labels, info, maps, os.path.join(out, "qc", f"{m['stem']}_QC.png"), roi_df)
        np.save(os.path.join(out, "qc", f"{m['stem']}_labels.npy"), labels.astype(np.uint8))
        if tiff:
            export_tiffs(m, os.path.join(out, "tiff"))
        if trait_tiffs:
            export_trait_tiffs(m, mask, labels, maps, os.path.join(out, "trait_maps"),
                               tissue=(leaf, stem, rim) if tissues else None)
        w = df[df["Tissue"] == "whole"]
        print(f"   mask {mask.sum()} px; {mode}; ROI areas: "
              + ", ".join(f"{a}={b}" for a, b in zip(w['Roi'], w['Area.px']))
              + (f"; stem {int(stem.sum())} px ({stem.sum() / mask.sum():.1%})" if tissues else ""))
    if tables:
        res = pd.concat(tables, ignore_index=True)
        t0 = pd.to_datetime(res["Date"] + " " + res["Time"]).min()
        res["TOE.h"] = (pd.to_datetime(res["Date"] + " " + res["Time"]) - t0).dt.total_seconds() / 3600
        res.to_csv(os.path.join(out, "PS2_traits_per_ROI.csv"), index=False)
        print(f"wrote {os.path.join(out, 'PS2_traits_per_ROI.csv')}  ({len(res)} rows)")
        return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder")
    ap.add_argument("--rois", help="ROI CSV (image, roi, label, shape, x, y, r) from roi_editor.html")
    ap.add_argument("--assign", choices=["inside", "grow"], default="inside",
                    help="inside: only plant pixels inside each polygon; grow: polygons seed a "
                         "geodesic growth along the plant mask (for sprawling plants)")
    ap.add_argument("--out")
    ap.add_argument("--mask-level", type=float, help="Fm threshold (default: TmPamMaskLevel from INF)")
    ap.add_argument("--n-pots", type=int, default=5, help="pots per image for auto mode")
    ap.add_argument("--px-per-mm", type=float, help="optional scale to add Area.mm2")
    ap.add_argument("--tiff", action="store_true", help="also export every DAT as ImageJ TIFF stack")
    ap.add_argument("--no-trait-tiffs", action="store_true",
                    help="skip trait_maps/ (per-pixel Fv/Fm, NPQ ... images + mask + ROI labels)")
    ap.add_argument("--no-tissues", action="store_true",
                    help="skip the leaf / leaf_core / stem split (whole plant only)")
    ap.add_argument("--stem-width", type=int, default=8,
                    help="structures narrower than this (px) count as stem/petiole (default 8)")
    ap.add_argument("--previews-only", action="store_true",
                    help="only write RGB PNGs for drawing ROIs, no analysis")
    a = ap.parse_args()
    process_folder(a.folder, a.out, a.rois, a.mask_level, a.n_pots, a.assign, a.tiff,
                   a.px_per_mm, a.previews_only, not a.no_trait_tiffs,
                   not a.no_tissues, max(1, a.stem_width // 2))


if __name__ == "__main__":
    main()
