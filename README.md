# PhenoSight

Code, protocols and hardware for the **PhenoSight Plant Phenotyping Facility** at the
Boyce Thompson Institute (BTI), Ithaca, NY.

PhenoSight runs a PhenoVation CropReporter imaging system for high-throughput phenotyping of
plants in trays and pots: chlorophyll fluorescence, multispectral reflectance, RGB color and
growth, imaged repeatedly over the course of an experiment. This repository collects what the
facility uses around the instrument: 3D-printed accessories, data-handling scripts and analysis
pipelines.

Facility director: Magdalena M. Julkowska (mmj55@cornell.edu).

## Contents

| Folder | What it holds |
| --- | --- |
| [`code/arabidopsis_tray_pipeline/`](code/arabidopsis_tray_pipeline/) | Image analysis of small plants in 5 × 4 trays (e.g. Arabidopsis) from the raw CropReporter files: pot grid, plant masks that follow each pot through time, growth, chlorophyll-fluorescence, spectral and color traits per plant per imaging round, QC overlays, and an R Markdown template that compares each treatment level with its control. |
| [`code/`](code/) | Scripts for organizing raw CropReporter and RGB camera files into one folder per experiment. |
| [`3D printed accessories/`](3D%20printed%20accessories/) | STL files for adapter rings that hold pots of different sizes (inner diameter 89, 113 and 135 mm) in the imaging trays. |

## Arabidopsis tray pipeline — quick start

```bash
pip install -r code/arabidopsis_tray_pipeline/requirements.txt

# 1. raw files -> one table per pot per imaging round (+ QC overlays)
python code/arabidopsis_tray_pipeline/phenosight_aph.py <experiment>/raw <experiment>/aph_output --workers 4

# 2. statistics and figures: copy the template next to the experiment, set its params, knit
Rscript -e 'rmarkdown::render("treatment_vs_control.Rmd")'
```

The full procedure, quality checks and outputs are described in the
[pipeline README](code/arabidopsis_tray_pipeline/README.md) and the
[standard operating procedure](code/arabidopsis_tray_pipeline/SOP.md).

## Conventions used across the facility

- **Pot positions** in 5 × 4 trays are named as seen standing at the imaging chamber door,
  with the white column on the far side: **A1 is the upper-left pot and E4 the lower-right
  pot**; letters A–E run left to right, numbers 1–4 from the far side toward the door.
- **Units** are metric throughout. Rosette areas are reported in mm², scaled by the pot's top
  width (63.5 mm for the standard trays).
- **Tray codes** in CropReporter file names increase by 20 per tray; the smallest code is Tray01.

## Citing and reuse

The contents of this repository are released under CC0 1.0 (see [LICENSE](LICENSE)). If you
use the analysis pipeline in published work, please acknowledge the PhenoSight Plant
Phenotyping Facility at the Boyce Thompson Institute.
