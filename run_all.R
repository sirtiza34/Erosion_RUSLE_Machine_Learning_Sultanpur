## =============================================================================================
##  run_all.R  --  runs the complete machine-learning analysis from the repository root
##  Steps 01-05 need only data/; steps 06-07 also need the Zenodo rasters in rasters/.
## =============================================================================================
steps <- c("01_screening", "02_spatial_folds", "03_nested_cv", "04_evaluation", "05_importance")
if (file.exists("rasters/E_P1_masked.tif")) steps <- c(steps, "06_final_maps", "07_autocorrelation")
for (s in steps) {
  cat("\n####################", s, format(Sys.time()), "####################\n")
  source(file.path("R", paste0(s, ".R")), local = new.env())
}
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
