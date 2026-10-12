# Machine-learning emulation of RUSLE-derived erosion classes, Sultanpur district (India)

R code and sample data for the machine-learning part of *A GIS-based approach for modelling soil erosion
susceptibility using integrated empirical and machine learning models* (Environment, Development and
Sustainability). C5.0 decision tree, SVM and MLPNN are trained on RUSLE-derived erosion classes and evaluated
by nested, spatially blocked cross-validation.

## Contents
| Path | Content |
|---|---|
| `data/sample_points_15000.shp` | 15,000 sample points (WGS 84 / UTM 44N): `id`, RUSLE soil loss `E` (t ha⁻¹ yr⁻¹), `stratum` (0: E < 1, 1: 1 ≤ E < 5, 2: E ≥ 5), `label` (E ≥ 5), `weight` (cells represented), GV01–GV17 and the R, K, LS and C values |
| `data/spatial_folds_5km.csv` | outer fold, 5 km block and buffer flags of each point |
| `R/00_setup.R` | settings, paths, packages and shared functions |
| `R/01_screening.R` | association matrix, GVIF, overlap with the RUSLE factors |
| `R/02_spatial_folds.R` | 5 km blocks, five folds, 1.2 km buffer |
| `R/03_nested_cv.R` | nested spatial CV, both designs, both variable sets |
| `R/04_evaluation.R` | performance, CIs, paired comparisons, confusion matrices |
| `R/05_importance.R` | permutation importance |
| `R/06_final_maps.R` | final Set B models and susceptibility maps – needs the Zenodo rasters |
| `R/07_autocorrelation.R` | autocorrelation ranges used for the block size – needs the Zenodo rasters |
| `run_all.R` | runs all steps in order |

The points were drawn from the baseline RUSLE soil-loss map in three strata (5,000 per stratum, seed 2026,
no two points in adjacent cells), and the variable values were extracted at the cell centres of the 30 m
rasters archived on Zenodo: https://doi.org/10.5281/zenodo.23297006

## Requirements
R 4.3.3 with C50 0.2.0, kernlab 0.9-33, nnet 7.3-19, blockCV 4.0-1, sf, terra, FNN and corrplot.

## How to run
Open R in the repository root and run `source("run_all.R")`, or run the scripts in `R/` one at a time.
Steps 01–05 need only `data/`. For steps 06–07, download the rasters from Zenodo into `rasters/`.
Results are written to `results/`. Finished cross-validation tasks are kept in `results/models/`, so an
interrupted run of `03_nested_cv.R` continues where it stopped.

## Licence
Code: MIT. Data: CC BY 4.0.

## Citation

Majid, S. I. (2026). Integrating RUSLE and Machine Learning for Soil Erosion Susceptibility Mapping in the Central Indo-Gangetic Plains: R Code and Sample Data (Version v1.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23314482

