#!/usr/bin/env python3
"""Redraw existing APH QC overlays in the 'standing at the door' orientation.

The pipeline reads the raw .DAT as a 2048 x 1500 portrait image, which is the transpose (a mirror
image) of the camera's own view. This script turns each overlay into the camera view rotated 180 deg,
i.e. as seen from the door of the PhenoSight chamber: white column at the top, door at the bottom,
A1 upper left, E4 lower right. Plant masks are not recomputed; old labels are erased and redrawn.

usage: python reorient_overlays.py <aph_output_dir> [--out overlays_door] [--first-code N]
"""
import os, sys, argparse
import numpy as np, pandas as pd, cv2

NROW, NCOL, W, H = 5, 4, 1500, 2048
PAD, OUT_W = 0.6, 700


def cell_label(a, b):
    return "ABCDE"[NROW - 1 - a] + str(NCOL - b)


def crop_box(g):
    c, s = np.cos(g["theta"]), np.sin(g["theta"])
    pts = [(g["x0"] + c * u * g["px"] - s * v * g["py"], g["y0"] + s * u * g["px"] + c * v * g["py"])
           for u, v in [(-PAD, -PAD), (NCOL + PAD, -PAD), (-PAD, NROW + PAD), (NCOL + PAD, NROW + PAD)]]
    xs, ys = zip(*pts)
    return max(0, int(min(xs))), max(0, int(min(ys))), min(W, int(max(xs))), min(H, int(max(ys)))


def to_overlay(g, box, sc, u, v):
    """grid coords (u along columns, v along rows, in cells) -> pixel in the old overlay."""
    c, s = np.cos(g["theta"]), np.sin(g["theta"])
    x = g["x0"] + c * u * g["px"] - s * v * g["py"] - box[0]
    y = g["y0"] + s * u * g["px"] + c * v * g["py"] - box[1]
    return x * sc, y * sc


def reorient(src, dst, g):
    img = cv2.imread(src)
    if img is None:
        return False
    box = crop_box(g); sc = OUT_W / (box[2] - box[0])
    h, w = img.shape[:2]
    # 1. erase the old yellow labels (drawn in the upper-left of each cell)
    mask = np.zeros((h, w), np.uint8)
    B, G, R = [img[..., i].astype(int) for i in range(3)]
    yellow = (R > 80) & (G > 80) & (B < 0.7 * np.minimum(R, G)) & (np.abs(R - G) < 50)
    for a in range(NROW):
        for b in range(NCOL):
            x1, y1 = to_overlay(g, box, sc, b + 0.0, a + 0.0)
            x2, y2 = to_overlay(g, box, sc, b + 0.5, a + 0.32)
            xa, xb = int(max(0, x1)), int(min(w, x2)); ya, yb = int(max(0, y1)), int(min(h, y2))
            mask[ya:yb, xa:xb] = yellow[ya:yb, xa:xb]
    mask = cv2.dilate(mask, np.ones((5, 5), np.uint8))
    img = cv2.inpaint(img, mask, 3, cv2.INPAINT_TELEA)
    # 2. camera view (transpose) turned 180 deg  ->  new[i, j] = old[h-1-j, w-1-i]
    new = np.ascontiguousarray(img.transpose(1, 0, 2)[::-1, ::-1])
    nh, nw = new.shape[:2]                                  # nh = w, nw = h
    T = lambda x, y: (h - 1 - y, w - 1 - x)                 # old (x, y) -> new (x, y)
    # 3. labels in the upper-left corner of each cell, in the new orientation
    fs = 0.0045 * g["px"] * sc * 1.0
    th = max(1, int(round(fs * 2)))
    for a in range(NROW):
        for b in range(NCOL):
            cs = [T(*to_overlay(g, box, sc, b + du, a + dv)) for du in (0, 1) for dv in (0, 1)]
            xl = min(p[0] for p in cs); yt = min(p[1] for p in cs)
            cell = (max(p[0] for p in cs) - xl)
            org = (int(xl + 0.08 * cell), int(yt + 0.22 * cell))
            cv2.putText(new, cell_label(a, b), org, cv2.FONT_HERSHEY_SIMPLEX, fs, (0, 0, 0), th + 2, cv2.LINE_AA)
            cv2.putText(new, cell_label(a, b), org, cv2.FONT_HERSHEY_SIMPLEX, fs, (0, 255, 255), th, cv2.LINE_AA)
    # 4. orientation strips
    band = max(18, nw // 28)
    top = np.full((band, nw, 3), 235, np.uint8); bot = np.full((band, nw, 3), 60, np.uint8)
    fsb = band / 30
    for im, txt, col in [(top, "white column (far side)", (40, 40, 40)), (bot, "DOOR (viewer stands here)", (255, 255, 255))]:
        (tw, tht), _ = cv2.getTextSize(txt, cv2.FONT_HERSHEY_SIMPLEX, fsb, 1)
        cv2.putText(im, txt, ((nw - tw) // 2, (band + tht) // 2), cv2.FONT_HERSHEY_SIMPLEX, fsb, col, 1, cv2.LINE_AA)
    new = np.vstack([top, new, bot])
    cv2.imwrite(dst, new, [cv2.IMWRITE_JPEG_QUALITY, 85])
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", help="pipeline output folder (holds images.csv and overlays/)")
    ap.add_argument("--out", default="overlays_door")
    ap.add_argument("--only", nargs="*", default=None, help="process only these overlay basenames")
    a = ap.parse_args()
    imgs = pd.read_csv(os.path.join(a.out_dir, "images.csv"))
    od = os.path.join(a.out_dir, a.out); os.makedirs(od, exist_ok=True)
    n = 0
    for _, r in imgs.iterrows():
        base = r["file"][:-4] + ".jpg"
        if a.only and base not in a.only:
            continue
        dst = os.path.join(od, base)
        if os.path.exists(dst) and not a.only:
            continue
        g = dict(x0=r.x0, y0=r.y0, px=r.px, py=r.py, theta=np.deg2rad(r.theta_deg))
        n += reorient(os.path.join(a.out_dir, "overlays", base), dst, g)
    print(f"{n} overlays written to {od}")


if __name__ == "__main__":
    main()
