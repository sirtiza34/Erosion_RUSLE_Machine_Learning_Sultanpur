## =============================================================================================
##  05_importance.R  --  permutation importance on the held-out spatial folds (strict design)
##  Each variable is permuted three times within the held-out points of each outer fold, and the
##  decrease in weighted ROC-AUC is averaged over the repetitions and then over the folds.
##  Needs the fold models written by 03_nested_cv.R.
##  Outputs: M3_importance_raw.csv, M3_importance.csv (Table S23), M3_importance_set{A,B}.png (Fig. 8)
## =============================================================================================
source("R/00_setup.R")

d <- prep(merge(read_samples(), read.csv(FOLDS), by = "id"))
imp <- list()
for (st in c("B", "A")) for (al in ALGOS) for (k in 1:K) {
  o  <- readRDS(model_file("strict", st, al, k))
  te <- d[d$fold == k, ]
  base <- wauc(score(o, al, te), te$label, te$weight)
  for (v in o$vars) for (r in 1:R_PERM) {
    set.seed(SEED + 100 * k + r); tp <- te; tp[[v]] <- tp[[v]][sample(nrow(tp))]
    imp[[length(imp) + 1]] <- data.frame(set = st, algo = al, fold = k, variable = v, rep = r,
                                         drop = base - wauc(score(o, al, tp), tp$label, tp$weight))
  }
  cat("importance: set", st, al, "fold", k, "| baseline ROC-AUC", round(base, 3), "\n")
}
imp <- do.call(rbind, imp)
write.csv(imp, file.path(OUT, "M3_importance_raw.csv"), row.names = FALSE)
f1 <- aggregate(drop ~ set + algo + variable + fold, imp, mean)                    # mean over repetitions
si <- do.call(data.frame, aggregate(drop ~ set + algo + variable, f1, function(x) c(mean = mean(x), sd = sd(x))))
si <- si[order(si$set, si$algo, -si$drop.mean), ]
write.csv(si, file.path(OUT, "M3_importance.csv"), row.names = FALSE)
cat("\nPermutation importance, set B (mean decrease in weighted ROC-AUC over the folds):\n")
print(si[si$set == "B", ], row.names = FALSE)
cat("\nPermutation importance, set A:\n"); print(si[si$set == "A", ], row.names = FALSE)

cols <- c(C50 = "#1b9e77", SVM = "#d95f02", MLP = "#7570b3")
for (st in names(SETS)) {
  x <- si[si$set == st, ]
  m <- tapply(x$drop.mean, list(sub("^GV\\d+_", "", x$variable), x$algo), mean)[, names(cols)]
  m <- m[order(rowMeans(m)), , drop = FALSE]
  png(file.path(OUT, sprintf("M3_importance_set%s.png", st)), 2000, 500 + 120 * nrow(m), res = 300)
  par(mar = c(4, 10, 1, 1))
  barplot(t(m), beside = TRUE, horiz = TRUE, las = 1, cex.names = 0.65, col = cols,
          xlab = "Mean decrease in weighted ROC-AUC", legend.text = names(cols),
          args.legend = list(x = "bottomright", bty = "n", cex = 0.7))
  invisible(dev.off())
}
