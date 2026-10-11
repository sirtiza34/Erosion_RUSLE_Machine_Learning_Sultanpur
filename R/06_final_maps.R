## =============================================================================================
##  06_final_maps.R  --  final Set B models and susceptibility maps
##  The Set B models are refitted to all non-erosion and erosion points with the setting selected
##  most often across the outer folds (ties: higher inner ROC-AUC; SVM: most frequent C and median
##  sigma; SVM Platt coefficients averaged over the folds) and applied to every analysed cell.
##  Needs: results of 03_nested_cv.R and the Zenodo rasters in rasters/ (RUSLE soil loss and the
##  nine Set B variables, named as in SET_B, e.g. rasters/GV01_Elevation.tif).
##  Outputs: maps/SUS_setB_{C50,SVM,MLP}.tif, M4_final_params.csv,
##           M4_class_tables.csv (Table S24), M4_map_agreement.csv
## =============================================================================================
source("R/00_setup.R")
dir.create(file.path(OUT, "maps"), showWarnings = FALSE)

d <- prep(read_samples())
d$y <- factor(ifelse(d$label == 1, "erosion", "none"), levels = c("none", "erosion"))
tun <- read.csv(file.path(OUT, "M1_tuning_strict.csv"))
parse_p <- function(s) { kv <- strsplit(strsplit(s, "; ")[[1]], "=")
  setNames(lapply(kv, function(z) as.numeric(z[2])), sapply(kv, `[`, 1)) }
final_params <- function(st, al) {
  x <- tun[tun$set == st & tun$algo == al, ]
  if (al == "SVM") { p <- lapply(x$params, parse_p); Cs <- sapply(p, `[[`, "C")
    return(list(sigma = median(sapply(p, `[[`, "sigma")), C = as.numeric(names(which.max(table(Cs)))))) }
  tab <- table(x$params); best <- names(tab)[tab == max(tab)]
  if (length(best) > 1) { y <- x[x$params %in% best, ]; best <- y$params[which.max(y$inner_wAUC)] }
  p <- parse_p(best); if (al == "C50") p$winnow <- as.logical(p$winnow); p
}

## ---- analysed cells ----
Em  <- rast(E_RASTER)
gvs <- rast(file.path(RASTERS, paste0(SET_B, ".tif"))); names(gvs) <- SET_B
compareGeom(Em, gvs)
ev_all <- values(Em, mat = FALSE); idx <- which(!is.na(ev_all)); Ev <- ev_all[idx]; rm(ev_all)
vals <- prep(as.data.frame(extract(gvs, idx)))
cat("Analysed cells:", length(idx), "| incomplete rows:", sum(!complete.cases(vals)), "\n")

## ---- final models and probability maps ----
trn <- d[d$stratum != 1, ]
probs <- list(); used <- list()
for (al in c("C50", "MLP", "SVM")) {
  p  <- final_params("B", al)
  used[[al]] <- data.frame(algo = al, params = paste(names(p), unlist(p), sep = "=", collapse = "; "))
  fm <- fit_model(al, p, trn, SET_B)
  if (al == "SVM") cf <- rowMeans(sapply(1:K, function(k) coef(readRDS(model_file("strict", "B", "SVM", k))$platt)))
  f <- switch(al,
              C50 = function(ii) predict(fm$model, vals[ii, SET_B], type = "prob")[, "erosion"],
              MLP = function(ii) as.numeric(predict(fm$model, apply_sc(vals[ii, ], fm$scaler))),
              SVM = function(ii) plogis(cf[1] + cf[2] * fm$sgn *
                                          as.numeric(predict(fm$model, apply_sc(vals[ii, ], fm$scaler), type = "decision"))))
  size <- if (al == "SVM") 10000 else 250000
  pr <- numeric(length(idx))
  for (s in seq(1, length(idx), by = size)) { ii <- s:min(length(idx), s + size - 1); pr[ii] <- f(ii) }
  probs[[al]] <- pr
  r <- rast(Em); values(r) <- NA; r[idx] <- pr
  writeRaster(r, file.path(OUT, "maps", sprintf("SUS_setB_%s.tif", al)), overwrite = TRUE, datatype = "FLT4S")
  cat("map saved:", al, "\n")
}
used <- do.call(rbind, used)
write.csv(used, file.path(OUT, "M4_final_params.csv"), row.names = FALSE)
cat("\nSelected final settings:\n"); print(used, row.names = FALSE)

## ---- susceptibility classes: area, mean soil loss and capture of E >= 5 ----
lv <- c("Very low", "Low", "Moderate", "High", "Very high")
cls_tab <- function(pr, breaks) {
  cl <- cut(pr, breaks, include.lowest = TRUE, labels = lv)
  n  <- tabulate(cl, 5); e5 <- tapply(Ev >= 5, cl, sum); e5[is.na(e5)] <- 0
  data.frame(class = lv, area_km2 = round(n * CELL_KM2, 1), share_pct = round(100 * n / length(pr), 2),
             mean_E = round(as.numeric(tapply(Ev, cl, mean)), 2),
             pct_of_class_with_E5 = round(100 * as.numeric(e5) / pmax(n, 1), 2),
             pct_of_all_E5_captured = round(100 * as.numeric(e5) / sum(Ev >= 5), 1))
}
tabs <- do.call(rbind, lapply(names(probs), function(al) rbind(
  cbind(algo = al, scheme = "equal interval 0.2", cls_tab(probs[[al]], seq(0, 1, 0.2))),
  cbind(algo = al, scheme = "quintiles", cls_tab(probs[[al]], quantile(probs[[al]], seq(0, 1, 0.2)))))))
write.csv(tabs, file.path(OUT, "M4_class_tables.csv"), row.names = FALSE)
cat("\nSusceptibility classes:\n"); print(tabs, row.names = FALSE)

## ---- agreement between the three maps ----
set.seed(SEED); sub <- sample(length(idx), min(100000, length(idx)))
rho <- round(cor(sapply(probs, function(p) p[sub]), method = "spearman"), 3)
hi  <- sapply(probs, function(p) p >= 0.6)
agr <- data.frame(measure = c("p >= 0.6 in at least one model (km2)", "p >= 0.6 in all three models (km2)"),
                  value = c(sum(rowSums(hi) >= 1), sum(rowSums(hi) == 3)) * CELL_KM2)
write.csv(rbind(agr, data.frame(measure = paste("Spearman", rownames(rho)[row(rho)], colnames(rho)[col(rho)])[upper.tri(rho)],
                                value = rho[upper.tri(rho)])),
          file.path(OUT, "M4_map_agreement.csv"), row.names = FALSE)
cat("\nSpearman correlation between the probability maps (100,000 random cells):\n"); print(rho)
print(agr, row.names = FALSE)
