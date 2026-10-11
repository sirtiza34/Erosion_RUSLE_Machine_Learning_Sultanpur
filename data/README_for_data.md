## Data Availability

The datasets required to run the scripts are archived on Zenodo, accessible at
https://doi.org/10.5281/zenodo.XXXXXXX. The archive contains the 30 m rasters (WGS 84 / UTM zone 44N) of the
seventeen geo-environmental variables, the RUSLE soil-loss maps (baseline and two sensitivity scenarios) and
the susceptibility probability maps of the final C5.0 DT, SVM and MLPNN models. The 15,000 sample points in
`data/` were extracted from these rasters at the cell centres, so steps 01–05 run without them. For steps
06–07, all rasters should be downloaded into the `rasters/` folder of the repository and then used.

Data Citation: Majid, S. I. (2026). Rasters for RUSLE soil loss and machine-learning erosion susceptibility,
Sultanpur district, India [Data set]. Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX
