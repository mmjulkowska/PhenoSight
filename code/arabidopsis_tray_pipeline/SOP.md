# SOP: Arabidopsis tray image analysis with the PhenoSight APH pipeline

**BOYCE THOMPSON INSTITUTE · PHENOSIGHT PLANT PHENOTYPING FACILITY**

Standard operating procedure · pipeline `phenosight_aph.py` v0.4.1 · prepared by Magdalena M. Julkowska, Director, PhenoSight Facility (mmj55@cornell.edu)

## 1. Purpose and scope

This SOP turns the raw CropReporter files of a PhenoSight experiment on Arabidopsis into one table of measurements per plant per imaging round, and from there into statistics and figures. It replaces the ROI and mask step of the PhenoVation software for small plants grown in trays.

**Use it for:** Arabidopsis (or similar rosette plants) in trays of 20 pots arranged 5 × 4, one plant per pot, imaged repeatedly over days to weeks, including experiments in which the camera height changes.

**What it gives you, per pot per imaging round:**

- **Growth:** rosette area in mm², convex hull area, solidity (compactness).
- **Chlorophyll fluorescence:** Fv/Fm, Fq'/Fm', NPQ, qP, qN, qL, ϕNO, ϕNPQ, npq(t).
- **Pigment and color indices:** chlorophyll index (ChlIdx), anthocyanin index (AriIdx), NDVI, hue, green index.
- **Time:** read from the instrument's own `HDR_*.INF` files.
- **QC:** an overlay image per tray image showing the pot grid, labels and plant masks.

**Why it exists:** the camera is raised as the plants grow. In the vendor software every camera move needs new ROIs drawn by hand, and pixel areas change with camera height. This pipeline finds the pot grid again in every image, uses the pot itself as the ruler (63.5 mm top width), follows each pot through time so empty pots stay empty, and filters algae on the soil.

**Not covered:** single large plants per pot, whole-tray monocultures, and side-view (RGB) imaging. Those use the facility's PS2/RGB export workflows.

## 2. Before you start

### 2.1 Raw files from the CropReporter

For every tray and imaging round the instrument writes a set of files that share one name stem, `<prefix>_<experiment>_<tray code>_<round>`:

| File | Contents | Used for |
| --- | --- | --- |
| `APH_*.DAT` | Red and far-red reflectance (2 frames) | Pot grid, plant masks, growth traits |
| `PMT_*.DAT` | Chlorophyll fluorescence (9 frames) | Fv/Fm, Fq'/Fm', NPQ and other quenching traits |
| `SPC_*.DAT` | Spectral reflectance: 550 nm, far-red, NIR | ChlIdx, AriIdx, NDVI |
| `CLR_*.DAT` | RGB color | NDVI (red band), hue, green index |
| `HDR_*.INF` | Text header with `TimeStamp=` | Imaging time |
| `CHL_*.DAT` | Not used | — |

Keep all of them in one folder (usually `raw/`). Tray codes increase by 20 per tray; the smallest code is Tray01.

### 2.2 Design files

- **Tray layout** (Excel or CSV, see `templates/layout_template.csv`): one row per pot with `TrayID` (Tray01…), `Area` (A1…E4) and the treatment code (`Treatment`, or `PlantName` as in PhenoVation layout files).
- **Treatment coding** (CSV, see `templates/coding_template.csv`): treatment code → label, series, level and role (Section 5).
- **Balance check:** number of treatments × replicates must equal the number of pots (for example 28 × 20 = 560). Run it before any analysis.

### 2.3 Pot positions

Positions are named as seen **standing at the chamber door**, with the white column on the far side: **A1 is the upper-left pot and E4 the lower-right pot**. Letters A–E run left to right and numbers 1–4 run from the far side toward the door. The layout file must use the same convention.

### 2.4 Software

- Python ≥ 3.9 with `numpy`, `scipy`, `scikit-image`, `opencv-python-headless`, `pandas`: `pip install numpy scipy scikit-image opencv-python-headless pandas`
- R ≥ 4.2 with `data.table`, `ggplot2`, `cowplot`, `scales`, `readxl`, `rmarkdown` for the statistics.
- The scripts in this folder: `phenosight_aph.py` (image analysis), `reorient_overlays.py` (redraws overlays from versions before 0.4.1), and `analysis/treatment_vs_control.Rmd` (statistics and figures). Install the Python packages with `pip install -r requirements.txt`.

### 2.5 Cloud-synced folders

If the raw folder lives in Dropbox or Box, mark it **Available offline** before the run. Online-only placeholders either download slowly during the run or fail to read; images whose PMT, SPC or CLR files fail to read get growth traits only.

## 3. Procedure

1. **Set up the experiment folder.** Put the raw files in `<experiment>/raw/` and the layout and coding files in `<experiment>/`. The pipeline writes to `<experiment>/aph_output/`.
2. **Check the raw files are complete.** Count the files of each type; APH, PMT, SPC, CLR and HDR should have the same count (trays × rounds): `ls raw | cut -d_ -f1 | sort | uniq -c` A shortfall usually means the upload or sync has not finished. The pipeline still runs, but missing rounds or trait files will be missing from the output.
3. **Test on one tray first.** Run the smallest tray code alone and look through its overlays (Section 4) before committing to the whole experiment: `python phenosight_aph.py <experiment>/raw <experiment>/aph_output --trays 61633 --workers 1`
4. **Run the full experiment:** `python phenosight_aph.py <experiment>/raw <experiment>/aph_output --workers 4` Expect 2–4 seconds per image; trays run in parallel, images within a tray in time order. A 28-tray, 35-round experiment (980 images) takes roughly 10–20 minutes on 4 workers.
5. **If the session has a time limit,** add `--max-minutes 2.5` (or any limit). The run stops cleanly, and the same command continues where it stopped: every image is checkpointed in `aph_output/checkpoints/`. Repeat until it reports that all trays are done.
6. **Read the end-of-run summary.** It lists the pots called empty and warns if any image has no time stamp. Keep this output with the run.
7. **When more files arrive later,** run the same command again. Only trays with new images are reprocessed; times, the empty-pot call and mm² values are recomputed for all.

**Options you may need**

| Option | Default | When to change it |
| --- | --- | --- |
| `--workers N` | CPUs − 1 | Lower it on a shared or small computer |
| `--max-minutes M` | off | Sessions with a time limit; rerun to continue |
| `--pitch-mm` | 63.5 | Pots other than the 63.5 mm square pots |
| `--first-code` | smallest code | Tray01 missing from the folder |
| `--times FILE` | INF files | No `HDR_*.INF` files: give PS2 exports or a `tray_code,round,datetime` CSV |
| `--establish-days`, `--establish-k` | 7 days, 3 images | Change the empty-pot rule (Section 7) |
| `--hi`, `--hi-max`, `--lo` | 0.32, 0.45, 0.22 | Mask thresholds; change only after inspecting overlays |
| `--no-traits` | off | Growth traits only (faster, no PMT/SPC/CLR needed) |
| `--trays CODE …` | all | Test or rerun selected trays |

To rerun from scratch with new thresholds, move the old `aph_output/` aside and run again.

## 4. Quality control

Do every check below before any statistics. A decoding or mask problem produces numbers that look plausible.

- [ ] **Overlays.** Open `aph_output/overlays/` and page through at least one tray from start to end, plus a few images from the last rounds of every tray. Each overlay shows the tray as seen from the door (white column strip at the top, DOOR strip at the bottom), the pot grid in white, the labels A1 (upper left) to E4 (lower right), and each pot's mask in its own color. Check that:
  - the grid sits on the pot walls in every image, including after camera moves;
  - masks cover the leaves and not algae, moss or the pot rim;
  - overlapping rosettes late in the trial are split sensibly between pots;
  - empty pots carry no mask.
- [ ] **Counts.** After joining to the layout, the number of plants must equal the number of pots in the layout minus the empty pots, and every treatment must have its expected number of plants. Any shortfall is a decoding error until proven otherwise; look at *which* positions are missing.
- [ ] **Position check.** Treatment effects should be visible with the layout as decoded. If a strong known contrast (for example two watering levels) shows nothing, check the pot-position convention before anything else.
- [ ] **Empty pots.** Compare the printed list of empty pots with the experiment notes. A late-germinating plant can be reinstated from `area_mm2_raw`.
- [ ] **Grid flags.** Images with `qc_grid_weak = TRUE` in `images.csv` (for example a tray not under the camera) are excluded from the analysis.
- [ ] **Missing trait files.** In `images.csv`, `has_pmt`, `has_spc` and `has_clr` show which images had their fluorescence, spectral and color files. Count the FALSE values and report them; download the files and rerun the affected trays if needed.
- [ ] **Time stamps.** Every image should have a `datetime`. Rounds are not synchronized across trays in every experiment, so always join times per tray and round.

## 5. Statistical analysis

The facility's standard analysis compares each treatment level on its own with the control grown under the same conditions, at every imaging round and at the end of the trial. It is implemented in `analysis/treatment_vs_control.Rmd`.

1. **Copy the notebook** next to the experiment and set its `params`: an experiment ID, the pipeline output folder, the layout file, the coding file, the minimum number of fluorescence pixels (default 20) and the end-of-trial window (default the final 3 days).
2. **Define the comparison groups** in the coding file (`templates/coding_template.csv`): for each treatment code a `Label`, a `Series` (for example the product, which shares a figure), a `Level` (for example the rate, which orders the panels) and a `Role` — `treatment`, `control` (exactly one) or `reference` (for example a well-watered mock).
3. **Knit the notebook.** It:
   - drops images with a weak grid fit, sets fluorescence traits to missing where fewer than 20 valid pixels exist, and drops every trait once a pot has no plant;
   - balance-checks the layout and prints plants per group;
   - tests each treatment level against its control with Welch's t-test at every imaging round, with Cohen's d (pooled SD) and Benjamini–Hochberg correction across the treatment levels within each trait and round;
   - repeats the test on each plant's mean over the final days;
   - writes the tables and one figure per series per trait (levels side by side, each with the control, a −log10(p\_BH) strip and a Cohen's d strip), plus an effect-size heatmap.
4. **Report two numbers for every effect:** the difference at the end of the trial and the number of imaging rounds at which it was significant. One alone can mislead: a difference present during growth may close by the end.
5. **Write a README into the output folder** with the experiment, the decode, the exclusions, the file list and any known issues, so others can reuse the tables without the notebook.

**Notes on the design**

- Use the tests on the measured rounds when all trays are imaged together in each round. When rounds are not synchronized across trays, smooth each plant's time course onto a common time grid first (smoothing spline, 12 h grid).
- Fluorescence and spectral traits change between pre-dawn and evening rounds. Compare treatments only within a round, and keep the same rounds for every plant in the end-of-trial window.
- If the design has blocks (groups of trays each with their own control), compare each treatment with the control in its own block.

## 6. Outputs

The pipeline writes four things to `aph_output/`:

| File | One row per | Contents |
| --- | --- | --- |
| `pots.csv` | pot × image (20 per image) | All traits, time, tray and position; the table to analyze |
| `images.csv` | image | Grid fit, pot size in pixels, exposure, time, `qc_grid_weak`, `has_pmt`/`has_spc`/`has_clr` |
| `overlays/*.jpg` | image | QC picture: grid, labels, masks, as seen from the door |
| `checkpoints/` | tray | Progress files that let a run resume; safe to delete after the run |

**Main columns of `pots.csv`**

| Column | Meaning |
| --- | --- |
| `TrayID`, `Area`, `Roi_No` | Tray (Tray01…), position (A1…E4), matching PhenoVation ROI number |
| `round`, `datetime`, `TOE_h` | Imaging round, time from the INF file, hours since the first image |
| `area_mm2` | Rosette area in mm²; zero in pots called empty |
| `area_mm2_raw` | The same before the empty-pot rule |
| `area_cellfrac`, `hull_cellfrac`, `solidity` | Area and convex hull as a fraction of one pot; compactness |
| `Fv_Fm`, `Fq_Fm`, `NPQ`, `qP`, `qN`, `qL`, `phiNO`, `phiNPQ`, `npq_t` | Chlorophyll fluorescence, mean over the plant's pixels |
| `n_px_fluo` | Number of pixels behind the fluorescence values |
| `ChlIdx`, `AriIdx`, `NDVI`, `Hue`, `GreenIdx` | Pigment and color indices |
| `pot_status` | `plant`, `empty` or `undetermined` |
| `algae_ring_frac`, `hi_used`, `rejected_cellfrac` | Mask QC: algae on the soil, the threshold used, area rejected as not plant |

The analysis notebook adds, in its own output folder: the per-round and end-of-trial result tables, each plant's end-of-trial values, the time-course figures and the effect-size heatmap.

## 7. How it works

Each tray's images are processed in time order, and each pot keeps a memory of its plant from one image to the next.

1. **Pot grid.** A 5 × 4 template of pot walls is matched to the red image and refined for position, pot spacing and rotation. This is redone in every image, so camera moves need no manual step.
2. **Plant signal.** Leaves reflect far-red and absorb red. The index `q = (FR − k·R) / (FR + k·R + 30)` is computed with `k` estimated from each pot's own soil, so bare soil sits at q ≈ 0 whatever its moisture or exposure.
3. **Mask.** Young seedlings on clean soil are segmented with two thresholds (q > 0.22 connected to q > 0.32). Once a pot shows algae on its soil, or the rosette passes 6% of the pot, the mask is built on a smoothed signal that removes the fine-grained algae texture but keeps leaf blades. The threshold rises in pots, and at times, where algae is present (up to q = 0.45).
4. **Assignment to pots.** Objects within one pot belong to it; objects bridging pots (touching rosettes, bolting stems) are split between them. An object is kept only if it touches the plant's footprint from the previous image (allowing for growth) or the pot center, so algae patches and debris away from the rosette are rejected.
5. **Empty pots.** A pot is called empty if a plant is found in fewer than 3 images during the tray's first 7 days. Anything later detected in an empty pot is set to zero in `area_mm2`; the raw value is kept.
6. **Scale.** Area in pot fractions is converted to mm² with the pot's 63.5 mm top width, which makes it comparable across camera heights. Absolute values are a few percent low because the rosette lies below the rim; comparisons between treatments imaged together are unaffected.
7. **Fluorescence, spectral and color traits.** For each image the matching PMT, SPC and CLR files are read on the same pixel grid. Each trait is calculated per pixel with the PhenoVation formulas and averaged over the plant's mask, using only pixels with Fm above 500 counts and values within the valid range.

**Validation against the PhenoVation software** (a 28-tray Arabidopsis experiment, same pots, correlation across pot-images)

| Measurement | r |
| --- | --- |
| Rosette area (rounds with correctly placed vendor ROIs) | 0.98–0.99 |
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

Fq'/Fm' and qP read about 0.05 higher than the vendor values, because the plant masks differ; comparisons within an experiment are unaffected.

## 8. Troubleshooting and known limitations

| Symptom | Likely cause | What to do |
| --- | --- | --- |
| Run fails with "Resource deadlock avoided" or files read as empty | Cloud-sync placeholders not downloaded | Mark the folder Available offline, wait for the sync, rerun |
| Some images have no fluorescence or color values | PMT/SPC/CLR file missing or not yet synced | Check `has_pmt`/`has_spc`/`has_clr`; rerun the tray once the files are there |
| `TOE_h` empty for some images | No `HDR_*.INF` file | Copy the INF files into `raw/` or give `--times` |
| Grid misplaced on a few images | Tray not under the camera, or an unusual view | Those images carry `qc_grid_weak`; exclude them |
| Masks spread over the soil | Heavy algae or moss beyond the filter | Inspect `hi_used` and `algae_ring_frac`; raise `--hi-max` and rerun that tray |
| Seedlings missed in the first days | Plants near the resolution limit (1–2 leaves) | Expected; early areas are noisier. Do not base conclusions on the first rounds alone |
| Too many or too few empty pots | Germination later than the 7-day window | Check the list against the notes; adjust `--establish-days` or reinstate from `area_mm2_raw` |
| Treatment effects absent where expected | Layout or position convention mismatch | Re-check Section 2.3 and the counts in Section 4 before anything else |
| Run stopped partway | Session or time limit | Rerun the same command; it resumes from the checkpoints |

**Known limitations**

- Overlapping rosettes late in a trial are split between pots geometrically, which is an approximation. Thin bolting stems are partly counted.
- In pots with algae, the raised threshold trims a little off leaf margins. Algae follows watering, so in watering experiments check `hi_used` and `algae_ring_frac` by treatment.
- `q` and the band ratios are uncalibrated; they compare within an experiment and camera setting, not across instruments.
- ϕNPQ and qN agree least with the vendor software (r ≈ 0.7).
- The pipeline is built for 5 × 4 trays of 63.5 mm pots. Other tray formats need the grid dimensions and pot size changed in the code.

Questions and bug reports: mmj55@cornell.edu.
