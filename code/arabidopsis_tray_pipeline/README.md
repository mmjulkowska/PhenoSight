# Arabidopsis tray pipeline (`phenosight_aph.py`)

Image analysis for small plants grown in 5 × 4 trays (one plant per pot) on the PhenoSight
PhenoVation CropReporter. The pipeline reads the raw `.DAT` files directly, finds the pots in
every image, segments each rosette, follows every pot through time, and writes one row per
plant per imaging round with growth, chlorophyll-fluorescence, spectral and color traits.
It replaces the ROI and mask step of the vendor software for this type of experiment.

Current version: **0.4.1**. Step-by-step instructions are in [SOP.md](SOP.md).

## Files

| File | Purpose |
| --- | --- |
| `phenosight_aph.py` | The pipeline (command-line tool) |
| `reorient_overlays.py` | Redraws QC overlays written by versions before 0.4.1 in the door view |
| `analysis/treatment_vs_control.Rmd` | Statistics and figures: each treatment level against its control |
| `templates/layout_template.csv`, `templates/coding_template.csv` | Formats of the design files the notebook reads |
| `requirements.txt` | Python packages |
| `SOP.md` | Standard operating procedure |

## Why it exists

- **The camera moves.** The camera is raised as the plants grow. In the vendor software every
  move needs new ROIs drawn by hand; with the wrong ROIs, per-pot areas agree with the plants at
  r ≈ 0.5.
- **Pixel areas change with camera height.** One pot can span 74 px early in an experiment and
  over 100 px later, so pixel counts mix plant growth with camera position.

The pipeline fits the pot grid again in **every** image and uses the pot as the ruler:
`area_cellfrac` = pixels / pot pitch², and `area_mm2` = `area_cellfrac` × 63.5² (63.5 mm is the
top width of the standard tray pots). Both compare across camera heights.

## Raw files

Every tray × imaging round produces files sharing the stem `<prefix>_<experiment>_<tray code>_<round>`:

- `APH_*.DAT` — 2 frames, red and far-red reflectance (12-bit values shifted left by 4). Used for
  the grid, the masks and the growth traits.
- `PMT_*.DAT` — 9 chlorophyll-fluorescence frames: 0 Fdark, 1 F0, 2 Fm, 3 Fdarksat, 4 Flight,
  5 Fs', 6 Fm', 7 Flightsat, 8 F0'.
- `SPC_*.DAT` — 550 nm, far-red and NIR reflectance. `CLR_*.DAT` — red, green, blue.
- `HDR_*.INF` — text header; its `TimeStamp=` line (a Delphi serial date, decimal comma or dot)
  gives the imaging time. INF files in the raw folder are read automatically.

All `.DAT` files are headerless little-endian uint16 arrays of shape (frames, 2048, 1500) on one
shared pixel grid. The pipeline reads them as 2048-tall portrait images, which is the transpose
(a mirror image) of the camera's own landscape view. Tray codes increase by 20 per tray; the
smallest code is Tray01. Imaging rounds are not always synchronized across trays, so times are
joined per tray and round.

## Pot positions

Positions are named as seen **standing at the chamber door**, with the white column on the far
side: **A1 upper left, E4 lower right**; letters run left to right and numbers from the far side
toward the door. In the portrait array the pipeline works in, this is A1 bottom-right and E4
top-left. `Roi_No` gives the matching vendor ROI number (A1 = 20, E4 = 1). QC overlays are
drawn in the door view from version 0.4.1.

The labels were checked two ways: they are the only orientation that matches the vendor ROI
export pot by pot (r = 0.98, against ≤ 0.51 for the other three), and in a 28-tray treatment
experiment only this orientation recovers the treatment effects (with the other three, each
label mixes 12–16 treatments and the effects disappear). Trays must always be loaded the same
way round: the tray is symmetric, so a tray loaded rotated by 180° cannot be detected from the
images.

## How it works

Each tray's images are processed in time order; every pot keeps a memory of its plant.

1. **Grid.** A 5 × 4 wall template is matched at pot pitches from 50 to 140 px and refined to
   sub-pixel accuracy for origin, x/y pitch and rotation.
2. **Plant index.** `q = (FR − k·R) / (FR + k·R + 30)`, with `k` the soil's far-red/red ratio
   estimated **per pot** from that pot's non-plant pixels, so soil sits at q ≈ 0 whatever its
   moisture, exposure or algae.
3. **Mask.**
   - *Seedling mode* (clean soil, rosette under 6% of a pot): hysteresis, q > 0.22 connected to
     q > 0.32.
   - *Texture mode* (once a pot shows algae on its soil — `algae_ring_frac` ≥ 0.005 with a plant
     of ≥ 2% of a pot — or once the rosette passes 6% of a pot): plant cores are found on a
     grey-level opening of q (disk of 2.2% of the pitch), which removes the fine-grained algae
     texture but keeps leaf blades; leaf margins are added back within 2.5% of the pitch of a
     core.
   - The core threshold rises with algae in each pot and image: max(0.32, 98th percentile of
     the filtered soil signal + 0.03), capped at 0.45 (`hi_used`, `lo_used` are recorded).
4. **Assignment to pots.** An object inside one pot belongs to it; objects that bridge pots
   (touching rosettes, bolting stems) are split by watershed from each pot's center. An object
   is kept only if it touches the plant's footprint from the previous image (grown by 6% of a
   pot plus 5% per day) or the pot center; rejected area goes to `rejected_cellfrac`.
5. **Empty pots.** A pot is `empty` if a plant (≥ 0.2% of a pot, ≈ 8 mm²) is found in fewer
   than 3 images during the tray's first 7 days. Later signal in an empty pot is set to zero in
   `area_mm2`; the raw value stays in `area_mm2_raw`, and the largest late signal in
   `late_signal_cellfrac`. Trays imaged for less than 7 days are `undetermined`.
6. **Growth traits.** Area, convex hull, solidity, perimeter, diameter, number of objects,
   centroid offset (large values flag objects that may belong to a neighbor), band means.
7. **Fluorescence, spectral and color traits.** The matching `PMT_`, `SPC_` and `CLR_` files are
   read on the same grid; each trait is computed per pixel (PhenoVation definitions) and
   averaged over the pot's mask, using only pixels with Fm > 500 counts and values inside the
   valid range (`n_px_fluo` = pixels used):
   - Fv/Fm = (Fm − F0)/Fm; Fq'/Fm' = (Fm' − Fs')/Fm'; NPQ = Fm/Fm' − 1
   - qP = (Fm' − Fs')/(Fm' − F0'); qN = 1 − (Fm' − F0')/(Fm − F0); qL = qP · F0'/Fs'
   - ϕNO = 1/(NPQ + 1 + qL(Fm/F0 − 1)); ϕNPQ = 1 − Fq'/Fm' − ϕNO; npq(t) = 4.88/(Fm'/F0' − 1) − 1
   - ChlIdx = (1/R_FR − 1/R_NIR)·R_NIR; AriIdx = (1/R_550 − 1/R_FR)·R_NIR; NDVI = (R_NIR − R_red)/(R_NIR + R_red)
   - Hue from mean RGB; green index = (2G − R − B)/(R + G + B)

Absolute mm² values are a few percent low because the rosette lies below the pot rim, from which
the scale is taken; comparisons between treatments imaged together are unaffected.

## Validation

On a 28-tray Arabidopsis experiment (1032 images, 20,640 pot-images), compared pot by pot with
the vendor software, using for every image the ROI set that matched its camera height:

| Measurement | Pooled r |
| --- | --- |
| Rosette area (median per image 0.976) | 0.97 |
| Fv/Fm | 1.00 |
| Anthocyanin index | 0.99 |
| Chlorophyll index | 0.98 |
| NPQ | 0.97 |
| NDVI | 0.96 |
| qL | 0.91 |
| ϕNO | 0.90 |
| Fq'/Fm', qP, npq(t) | 0.88 |
| ϕNPQ | 0.71 |
| qN | 0.70 |

Fq'/Fm' and qP read about 0.05 higher than the vendor values because the plant masks differ.
rETR (Fq'/Fm' × a fixed light level) and qE/qI (frames not saved) are not computed.

## Running it

```bash
pip install -r requirements.txt
python phenosight_aph.py <experiment>/raw <experiment>/aph_output --workers 4
```

- Expect 2–4 s per image; trays run in parallel, images within a tray in sequence.
- The run is checkpointed per image (`aph_output/checkpoints/`). `--max-minutes M` stops cleanly
  after M minutes; the same command continues where it stopped.
- The run is incremental per tray: a tray is reprocessed when any of its images is new. Times,
  the empty-pot call and mm² values are recomputed on every run.
- Output from an older segmentation version is reprocessed automatically. To change thresholds,
  move the old output folder aside and run again.

| Option | Default | Use |
| --- | --- | --- |
| `--workers N` | CPUs − 1 | Parallel trays |
| `--max-minutes M` | off | Time-limited sessions; rerun to continue |
| `--trays CODE …` | all | Test or rerun selected trays |
| `--pitch-mm` | 63.5 | Pot top width in mm |
| `--first-code`, `--code-step` | smallest code, 20 | Tray numbering |
| `--times FILE …` | INF files | Extra time sources: INF folders, vendor exports or a `tray_code,round,datetime` CSV |
| `--establish-days`, `--establish-k` | 7, 3 | Empty-pot rule |
| `--hi`, `--hi-max`, `--lo` | 0.32, 0.45, 0.22 | Mask thresholds |
| `--no-tracking` | off | Independent images, fixed thresholds, no empty-pot call |
| `--no-traits` | off | Growth traits only |
| `--no-overlays` | off | Skip QC overlays |

If the raw folder is in Dropbox or Box, make it available offline first: online-only
placeholders download slowly or fail to read.

## Outputs

- `pots.csv` — one row per pot per image: `TrayID`, `Area`, `Roi_No`, `round`, `datetime`,
  `TOE_h`, growth traits, fluorescence/spectral/color traits, `pot_status` and mask QC columns.
- `images.csv` — one row per image: grid fit (`x0`, `y0`, `px`, `py`, `theta_deg`,
  `wall_score`), `pitch_px`, exposure, `pots_with_plant`, `qc_grid_weak`, `has_pmt`,
  `has_spc`, `has_clr`, time.
- `overlays/*.jpg` — one QC image per tray image: grid, labels and masks, as seen from the door.

## Known limitations

- Overlapping rosettes late in a trial are split geometrically; thin bolting stems are partly
  counted.
- Where algae is present, the raised threshold trims a little off leaf margins. Algae follows
  watering, so check `hi_used` and `algae_ring_frac` by treatment in watering experiments.
- 1–2-leaf seedlings are near the resolution limit, so the first days are noisier.
- `q` and the band ratios are uncalibrated: compare within an experiment and camera setting.
- Built for 5 × 4 trays; other formats need `NROW`, `NCOL` and the pot size changed.

## Version history

- **0.4.1** — QC overlays drawn as seen from the chamber door.
- **0.4** — Fluorescence, spectral and color traits from the `PMT_`, `SPC_` and `CLR_` files.
- **0.3** — Each tray followed through time: empty-pot call, algae filter, rising threshold,
  footprint gate.
- **0.2** — Per-image grid fit and hysteresis masks.

Contact: Magdalena M. Julkowska, Director, PhenoSight Facility, Boyce Thompson Institute
(mmj55@cornell.edu).
