## =============================================================================================
##  04_evaluation.R  --  performance of the pooled out-of-fold predictions
##  Weighted ROC-AUC (all points and extremes), PR-AUC, Brier skill score (probabilities rescaled
##  to the district prevalence) and threshold measures at p = 0.5, with 95% CIs from 1,000
##  bootstrap resamples of the spatial blocks, paired model comparisons, confusion matrices in km2
##  and the mean probability by RUSLE soil-loss class.
##  Outputs: M2_metrics.csv (Table 4, Tables S19-S20), M2_comparisons.csv (Table S21),
##           M2_confusion_km2.csv (Table S22), M2_prob_by_soil_loss.csv,
##           M2_roc.png (Fig. 5), M2_prob_by_soil_loss.png (Fig. 6), M2_calibration.png (Fig. 7)
## =============================================================================================
source("R/00_setup.R")
DESIGNS <- c("strict", "nogap")[file.exists(file.path(OUT, sprintf("M1_oof_predictions_%s.csv", c("strict", "nogap"))))]
cat("Designs found:", DESIGNS, "\n")

## ---- data, blocks, erosion share in the training data of each fold ----
d <- merge(read_samples(), read.csv(FOLDS), by = "id")
if (is.null(d$block)) {                                     # block of each point from the polygons
  blk <- st_read(BLOCKS, quiet = TRUE)
  pts <- st_as_sf(d, coords = c("x", "y"), crs = st_crs(blk))
  d$block <- sapply(st_intersects(pts, blk), function(i) i[1])
}
cat("Points without a block:", sum(is.na(d$block)), "\n")
pi_s <- function(design, k) {
  tr <- d$fold != k & d[[paste0("buf_excl_f", k)]] == 0
  if (design == "strict") tr <- tr & d$stratum != 1
  mean(d$label[tr])
}
P <- do.call(rbind, lapply(DESIGNS, function(ds) {
  p <- read.csv(file.path(OUT, sprintf("M1_oof_predictions_%s.csv", ds))); p$design <- ds
  p$pi_s <- sapply(1:K, function(k) pi_s(ds, k))[p$fold]; p }))
P <- merge(P, d[, c("id", "block")], by = "id")
a <- P$prob * PI / P$pi_s; b <- (1 - P$prob) * (1 - PI) / (1 - P$pi_s)
P$p_adj <- a / (a + b)                                      # probability rescaled to the district prevalence
P$key <- paste(P$design, P$set, P$algo, sep = "|")
keys <- unique(P$key)
L <- lapply(setNames(keys, keys), function(k) { x <- P[P$key == k, ]; x[order(x$id), ] })

metrics <- function(x) {
  ex <- x$stratum != 1
  br <- sum(x$weight * (x$p_adj - x$label)^2) / sum(x$weight)
  c(ROC_AUC = wauc(x$prob, x$label, x$weight),
    ROC_AUC_extremes = wauc(x$prob[ex], x$label[ex], x$weight[ex]),
    PR_AUC = wap(x$prob, x$label, x$weight),
    Brier = br, BSS = 1 - br / (PI * (1 - PI)),
    cls(x$prob, x$label, x$weight)[1:6])
}

## ---- estimates, fold statistics and block bootstrap ----
est <- t(sapply(keys, function(k) metrics(L[[k]])))
mn  <- colnames(est)
fs  <- do.call(rbind, lapply(keys, function(k) {
  m <- t(sapply(1:K, function(f) metrics(L[[k]][L[[k]]$fold == f, ])))
  data.frame(key = k, metric = mn, fold_mean = colMeans(m), fold_sd = apply(m, 2, sd),
             fold_min = apply(m, 2, min), fold_max = apply(m, 2, max))
}))
set.seed(SEED)
rows_by_block <- split(seq_len(nrow(L[[1]])), L[[1]]$block); nb <- length(rows_by_block)
boot <- array(NA, c(N_BOOT, length(keys), length(mn)), dimnames = list(NULL, keys, mn))
for (bb in 1:N_BOOT) {
  ii <- unlist(rows_by_block[sample.int(nb, nb, replace = TRUE)], use.names = FALSE)
  for (k in keys) boot[bb, k, ] <- metrics(L[[k]][ii, ])
  if (bb %% 100 == 0) cat("bootstrap", bb, "of", N_BOOT, "\n")
}
ci  <- apply(boot, c(2, 3), quantile, probs = c(0.025, 0.975), na.rm = TRUE)
lab <- do.call(rbind, strsplit(keys, "\\|"))
long <- do.call(rbind, lapply(seq_along(keys), function(i) data.frame(
  key = keys[i], design = lab[i, 1], set = lab[i, 2], algo = lab[i, 3], metric = mn,
  estimate = est[keys[i], mn], ci_low = ci[1, keys[i], mn], ci_high = ci[2, keys[i], mn])))
long <- merge(long, fs, by = c("key", "metric"))
write.csv(long[, -1], file.path(OUT, "M2_metrics.csv"), row.names = FALSE)
fmt <- function(m) sprintf("%.3f (%.3f-%.3f)", est[keys, m], ci[1, keys, m], ci[2, keys, m])
w1 <- data.frame(design = lab[, 1], set = lab[, 2], algo = lab[, 3]); w2 <- w1
for (m in c("ROC_AUC", "ROC_AUC_extremes", "PR_AUC", "BSS")) w1[[m]] <- fmt(m)
for (m in c("sensitivity", "specificity", "balanced_accuracy", "precision", "F1", "MCC")) w2[[m]] <- fmt(m)
cat("\nThreshold-free measures, pooled out-of-fold, 95% block-bootstrap CI:\n"); print(w1, row.names = FALSE)
cat("\nMeasures at probability 0.5:\n"); print(w2, row.names = FALSE)

## ---- paired comparisons on the same bootstrap resamples ----
cmp <- list()
add_cmp <- function(lbl, k1, k2) for (m in c("ROC_AUC", "PR_AUC", "MCC")) {
  dif <- boot[, k1, m] - boot[, k2, m]
  cmp[[length(cmp) + 1]] <<- data.frame(comparison = lbl, metric = m, difference = est[k1, m] - est[k2, m],
                                        ci_low = quantile(dif, 0.025, na.rm = TRUE),
                                        ci_high = quantile(dif, 0.975, na.rm = TRUE))
}
for (ds in DESIGNS) for (st in names(SETS)) for (pr in list(c("C50", "SVM"), c("C50", "MLP"), c("SVM", "MLP")))
  add_cmp(sprintf("%s, set %s: %s - %s", ds, st, pr[1], pr[2]),
          paste(ds, st, pr[1], sep = "|"), paste(ds, st, pr[2], sep = "|"))
if (length(DESIGNS) == 2) for (st in names(SETS)) for (al in ALGOS)
  add_cmp(sprintf("set %s, %s: strict - nogap", st, al),
          paste("strict", st, al, sep = "|"), paste("nogap", st, al, sep = "|"))
cmp <- do.call(rbind, cmp); rownames(cmp) <- NULL
cmp$significant <- cmp$ci_low > 0 | cmp$ci_high < 0
write.csv(cmp, file.path(OUT, "M2_comparisons.csv"), row.names = FALSE)
cat("\nPaired differences (95% block-bootstrap CI):\n")
print(transform(cmp, difference = round(difference, 4), ci_low = round(ci_low, 4), ci_high = round(ci_high, 4)),
      row.names = FALSE)

## ---- confusion matrices in km2 and mean probability by soil loss ----
cm <- do.call(rbind, lapply(keys, function(k) {
  x <- L[[k]]; v <- cls(x$prob, x$label, x$weight); p1 <- x$prob >= 0.5
  data.frame(key = k, TP_km2 = v[["TP"]], FN_km2 = v[["FN"]], FP_km2 = v[["FP"]], TN_km2 = v[["TN"]],
             FP_from_E_below1_km2 = sum(x$weight[p1 & x$stratum == 0]),
             FP_from_E_1to5_km2 = sum(x$weight[p1 & x$stratum == 1]))
}))
cm[, -1] <- round(cm[, -1] * CELL_KM2, 1)
write.csv(cm, file.path(OUT, "M2_confusion_km2.csv"), row.names = FALSE)
cat("\nConfusion matrices at probability 0.5 (km2 of the analysed area):\n"); print(cm, row.names = FALSE)

brk    <- c(0, 0.25, 0.5, 1, 2, 3, 4, 5, 10, Inf)
binlab <- c("<0.25", "0.25-0.5", "0.5-1", "1-2", "2-3", "3-4", "4-5", "5-10", ">=10")
pb <- do.call(rbind, lapply(keys, function(k) {
  x <- L[[k]]; bin <- cut(x$E, brk, right = FALSE, labels = binlab)
  data.frame(key = k, bin = binlab, mean_prob = round(as.numeric(tapply(x$prob, bin, mean)), 3))
}))
write.csv(pb, file.path(OUT, "M2_prob_by_soil_loss.csv"), row.names = FALSE)
cat("\nMean predicted probability by RUSLE soil-loss class:\n")
print(reshape(pb, idvar = "key", timevar = "bin", direction = "wide"), row.names = FALSE)

## ---- figures ----
cols <- c(C50 = "#1b9e77", SVM = "#d95f02", MLP = "#7570b3")
wroc <- function(s, y, w) { o <- order(-s); y <- y[o]; w <- w[o]
  list(fpr = c(0, cumsum(w * (1 - y)) / sum(w * (1 - y))), tpr = c(0, cumsum(w * y) / sum(w * y))) }
panels <- function(file, fun) {
  png(file.path(OUT, file), 2400, 1150 * length(DESIGNS), res = 300)
  par(mfrow = c(length(DESIGNS), 2), mar = c(4.5, 4.5, 2.5, 1))
  for (ds in DESIGNS) for (st in names(SETS)) fun(ds, st)
  invisible(dev.off())
}
panels("M2_roc.png", function(ds, st) {
  plot(0:1, 0:1, type = "n", xlab = "False positive rate", ylab = "True positive rate",
       main = sprintf("%s training, set %s", ds, st)); abline(0, 1, lty = 2, col = "grey60")
  for (al in names(cols)) { r <- with(L[[paste(ds, st, al, sep = "|")]], wroc(prob, label, weight))
    lines(r$fpr, r$tpr, col = cols[al], lwd = 1.5) }
  legend("bottomright", sprintf("%s  AUC %.3f", names(cols), est[paste(ds, st, names(cols), sep = "|"), "ROC_AUC"]),
         col = cols, lwd = 1.5, bty = "n", cex = 0.85)
})
panels("M2_calibration.png", function(ds, st) {
  plot(0:1, 0:1, type = "n", xlab = "Predicted probability (district prevalence)",
       ylab = "Observed share with E >= 5", main = sprintf("%s training, set %s", ds, st))
  abline(0, 1, lty = 2, col = "grey60")
  for (al in names(cols)) { x <- L[[paste(ds, st, al, sep = "|")]]
    bin <- cut(x$p_adj, unique(quantile(x$p_adj, seq(0, 1, 0.1))), include.lowest = TRUE)
    lines(tapply(x$p_adj * x$weight, bin, sum) / tapply(x$weight, bin, sum),
          tapply(x$label * x$weight, bin, sum) / tapply(x$weight, bin, sum), type = "b", pch = 16, col = cols[al]) }
  legend("topleft", names(cols), col = cols, lty = 1, pch = 16, bty = "n", cex = 0.85)
})
panels("M2_prob_by_soil_loss.png", function(ds, st) {
  plot(1:9, rep(0, 9), type = "n", ylim = c(0, 1), xaxt = "n", xlab = "", ylab = "Mean predicted probability",
       main = sprintf("%s training, set %s", ds, st))
  axis(1, 1:9, binlab, las = 2, cex.axis = 0.75); mtext("RUSLE soil loss (t/ha/yr)", 1, 3.3, cex = 0.8)
  abline(v = c(3.5, 7.5), lty = 3, col = "grey50")
  for (al in names(cols)) lines(1:9, pb$mean_prob[pb$key == paste(ds, st, al, sep = "|")], type = "b", pch = 16, col = cols[al])
  legend("topleft", names(cols), col = cols, lty = 1, pch = 16, bty = "n", cex = 0.85)
})
cat("\nEvaluation results saved in", OUT, "\n")
