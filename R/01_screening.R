## =============================================================================================
##  01_screening.R  --  redundancy among the 17 variables and their overlap with the RUSLE factors
##  Outputs: screen_association_matrix.csv/.png (Fig. S5), screen_gvif.csv and
##           screen_rusle_overlap.csv (Table S13)
## =============================================================================================
source("R/00_setup.R")

d <- read_samples()
fac <- c("R", "K", "LS", "C")
if (!all(fac %in% names(d)) && file.exists("data/rusle_factors_at_points.csv"))
  d <- merge(d, read.csv("data/rusle_factors_at_points.csv")[, c("id", fac)], by = "id")
con_gv <- setdiff(GV, CAT_GV)
for (v in CAT_GV) d[[v]] <- factor(d[[v]])
print(table(stratum = d$stratum))
for (v in CAT_GV) { cat("\n", v, "(rows = class, columns = stratum)\n"); print(table(d[[v]], d$stratum)) }

## district-representative subsample (erosion back to its district share)
set.seed(SEED)
s <- d[sample(nrow(d), 10000, replace = TRUE, prob = d$weight), ]
cat("\nErosion share in the representative subsample:", round(100 * mean(s$label), 2),
    "% (district", round(100 * PI, 2), "%)\n")

## ---- association measures, all on a 0-1 scale ----
spear  <- function(a, b) abs(cor(a, b, method = "spearman"))           # continuous - continuous
eta    <- function(y, f) sqrt(summary(lm(rank(y) ~ f))$r.squared)       # continuous - class
cramer <- function(f, g) {                                               # class - class
  tb  <- table(droplevels(f), droplevels(g))
  chi <- suppressWarnings(chisq.test(tb, correct = FALSE)$statistic)
  as.numeric(sqrt(chi / (sum(tb) * (min(dim(tb)) - 1))))
}
assoc <- function(a, b, ca, cb) {
  if (!ca && !cb) spear(a, b) else if (ca && cb) cramer(a, b)
  else if (ca) eta(b, a) else eta(a, b)
}
## generalized VIF (Fox & Monette 1992)
ld   <- function(M) as.numeric(determinant(M, logarithm = TRUE)$modulus)
gvif <- function(df, vars) {
  mm  <- model.matrix(as.formula(paste("~", paste(vars, collapse = " + "))), df)
  asg <- attr(mm, "assign")[-1]; X <- mm[, -1, drop = FALSE]
  keep <- apply(X, 2, sd) > 0; X <- X[, keep, drop = FALSE]; asg <- asg[keep]
  R <- cor(X)
  out <- t(sapply(seq_along(vars), function(k) {
    i <- which(asg == k)
    g <- exp(ld(R[i, i, drop = FALSE]) + ld(R[-i, -i, drop = FALSE]) - ld(R))
    c(GVIF = g, Df = length(i), GVIF_adj = g^(1 / (2 * length(i))))
  }))
  rownames(out) <- vars; out
}

## ---- variables against each other ----
n <- length(GV)
A <- matrix(1, n, n, dimnames = list(GV, GV))
for (i in 1:(n - 1)) for (j in (i + 1):n)
  A[i, j] <- A[j, i] <- assoc(s[[GV[i]]], s[[GV[j]]], GV[i] %in% CAT_GV, GV[j] %in% CAT_GV)
hi <- which(A >= 0.7 & upper.tri(A), arr.ind = TRUE)
cat("\nPairs with association >= 0.7:\n")
print(data.frame(var1 = GV[hi[, 1]], var2 = GV[hi[, 2]], assoc = round(A[hi], 3))[order(-A[hi]), ], row.names = FALSE)

vf <- gvif(s, GV)
cat("\nGVIF (adjusted GVIF = GVIF^(1/(2 Df)); 2.24 corresponds to VIF = 5):\n")
print(round(vf, 2))

## ---- overlap with the RUSLE factors (needs the R, K, LS and C values at the points) ----
if (all(fac %in% names(s))) {
  O <- sapply(fac, function(f) sapply(GV, function(v)
    if (v %in% CAT_GV) eta(s[[f]], s[[v]]) else spear(s[[v]], s[[f]])))
  O  <- round(O, 3)
  ov <- data.frame(variable = GV, O, max = apply(O, 1, max), strongest = fac[apply(O, 1, which.max)],
                   set_B = GV %in% SET_B)
  cat("\nAssociation of each variable with the RUSLE factors:\n")
  print(ov[order(-ov$max), ], row.names = FALSE)
  write.csv(ov, file.path(OUT, "screen_rusle_overlap.csv"), row.names = FALSE)
} else cat("\nR, K, LS and C not found at the points: overlap with the RUSLE factors skipped\n")

## ---- save ----
write.csv(round(A, 3), file.path(OUT, "screen_association_matrix.csv"))
write.csv(round(vf, 3), file.path(OUT, "screen_gvif.csv"))
png(file.path(OUT, "screen_association_matrix.png"), 2400, 2400, res = 300)
corrplot(A, method = "color", type = "upper", diag = FALSE, is.corr = FALSE, col.lim = c(0, 1),
         col = colorRampPalette(c("white", "orange", "red"))(100),
         tl.col = "black", tl.cex = 0.6, addCoef.col = "black", number.cex = 0.45)
invisible(dev.off())
cat("\nScreening results saved in", OUT, "\n")
