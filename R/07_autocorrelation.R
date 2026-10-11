## =============================================================================================
##  07_autocorrelation.R  --  range of spatial autocorrelation used to choose the block size
##  Variograms fitted by blockCV (automap) to 5,000 random cells of each raster.
##  Needs the Zenodo rasters in rasters/. Output: autocorrelation_ranges.csv (Table S15)
## =============================================================================================
source("R/00_setup.R")

cont <- c("GV01_Elevation", "GV02_Slope", "GV03_PlanCurvature", "GV04_ProfileCurvature",
          "GV05_DrainageDensity", "GV06_DrainageDistance", "GV07_TWI", "GV08_SPI", "GV09_STI",
          "GV13_LineamentDensity", "GV15_NDVI", "GV16_SBI")
r <- rast(file.path(RASTERS, paste0(cont, ".tif"))); names(r) <- cont
logE <- log(rast(E_RASTER)); names(logE) <- "log_soil_loss"

set.seed(SEED); ac_gv <- cv_spatial_autocor(r = r, num_sample = 5000, plot = FALSE, progress = FALSE)
set.seed(SEED); ac_e  <- cv_spatial_autocor(r = logE, num_sample = 5000, plot = FALSE, progress = FALSE)
print(ac_gv$range_table); print(ac_e$range_table)
cat("\nMedian range of the continuous variables:", round(ac_gv$range), "m\n")
cat("Range of log soil loss:", round(ac_e$range), "m\n")
write.csv(rbind(ac_gv$range_table, ac_e$range_table), file.path(OUT, "autocorrelation_ranges.csv"), row.names = FALSE)
