The codes in this folder are as follow:

RGB_Data_Exp??????.py - the code to move the files made by RGB/CF camera from individual folders containing rgb and jpeg images per imaging round + automatically calculated and masked files into one folder per experiment

move_files_Exo??????.py - the code to move CropReporter files into one experimental folder for analysis through PlantCV and/or PhenoVation DataAnalysis reporter


arabidopsis_tray_pipeline/ - image analysis of small plants (e.g. Arabidopsis) in 5 x 4 trays from the raw CropReporter files: per-plant growth, chlorophyll-fluorescence, spectral and color traits per imaging round, QC overlays, and an R Markdown template for treatment-vs-control statistics. See arabidopsis_tray_pipeline/README.md and SOP.md.

raw_roi_pipeline/ - opens raw CropReporter files of any layout without the vendor software; regions drawn in roi_editor.html; fluorescence, spectral and color traits per region, plus per-pixel trait maps. See raw_roi_pipeline/README.md.
