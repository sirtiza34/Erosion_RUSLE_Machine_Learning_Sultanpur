## =============================================================================================
##  00_setup.R  --  settings, packages, data loading and helper functions
##  Sourced by every other script. Run all scripts from the repository root, e.g.
##  source("R/03_nested_cv.R"), or everything with source("run_all.R").
## =============================================================================================

## ---- settings ----
SEED    <- 2026      # random seed (plus the fold number for fold-specific steps)
K       <- 5         # outer spatial folds
K_IN    <- 3         # inner spatial folds for hyperparameter tuning
BLOCK   <- 5000      # side of the square spatial blocks (m)
BUF     <- 1200      # training points closer than this to a held-out point are removed (m)
N_BOOT  <- 1000      # block-bootstrap resamples
R_PERM  <- 3         # permutations per variable and fold (permutation importance)
N_CORES <- max(1, parallel::detectCores() - 1)
CRS_UTM <- 32644     # WGS 84 / UTM zone 44N

## ---- paths (relative to the repository root) ----
SHP      <- "data/sample_points_15000.shp"   # the 15,000 sample points
FOLDS    <- "data/spatial_folds_5km.csv"     # outer folds, buffer flags and blocks (02_spatial_folds.R)
BLOCKS   <- "data/spatial_blocks_5km.gpkg"
RASTERS  <- "rasters"                        # Zenodo rasters, needed only by 06 and 07
E_RASTER <- file.path(RASTERS, "E_P1_masked.tif")   # baseline RUSLE soil loss (t/ha/yr)
OUT      <- "results"
MODELS   <- file.path(OUT, "models")
for (p in c(OUT, MODELS)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

## ---- district totals of the analysed area (30 m cells) ----
CELLS    <- c("0" = 1253732, "1" = 1126174, "2" = 69886)   # E < 1, 1 <= E < 5, E >= 5 t/ha/yr
PI       <- CELLS[["2"]] / sum(CELLS)                        # district share of E >= 5 (0.0285)
CELL_KM2 <- 9e-4

## ---- variables ----
GV <- c("GV01_Elevation", "GV02_Slope", "GV03_PlanCurvature", "GV04_ProfileCurvature",
        "GV05_DrainageDensity", "GV06_DrainageDistance", "GV07_TWI", "GV08_SPI", "GV09_STI",
        "GV10_HSG", "GV11_Lithology", "GV12_Geomorphology", "GV13_LineamentDensity",
        "GV14_SoilTexture", "GV15_NDVI", "GV16_SBI", "GV17_LULC")
CAT_GV <- c("GV10_HSG", "GV11_Lithology", "GV12_Geomorphology", "GV14_SoilTexture", "GV17_LULC")
SET_A  <- GV                                                  # all seventeen variables
SET_B  <- c("GV01_Elevation", "GV03_PlanCurvature", "GV04_ProfileCurvature", "GV05_DrainageDensity",
            "GV06_DrainageDistance", "GV10_HSG", "GV11_Lithology", "GV12_Geomorphology",
            "GV13_LineamentDensity")                          # nine variables independent of RUSLE
SETS   <- list(A = SET_A, B = SET_B)
ALGOS  <- c("C50", "SVM", "MLP")

## ---- packages ----
pk   <- c("sf", "terra", "blockCV", "FNN", "C50", "kernlab", "nnet", "corrplot")
miss <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install the missing packages first: install.packages(c(\"",
                       paste(miss, collapse = "\", \""), "\"))")
suppressPackageStartupMessages(invisible(lapply(pk, library, character.only = TRUE)))

## ---- sample points ----
## Reads the point shapefile. Shapefile field names are cut to 10 characters, so the full
## variable names are restored from their GVxx prefix. stratum, label and weight are derived
## from E if they are not stored in the file.
read_samples <- function(path = SHP) {
  p  <- st_read(path, quiet = TRUE)
  xy <- st_coordinates(p)
  d  <- st_drop_geometry(p)
  nm <- names(d)
  for (v in GV) { i <- which(toupper(substr(nm, 1, 4)) == substr(v, 1, 4)); if (length(i) == 1) nm[i] <- v }
  for (v in c("id", "stratum", "label", "weight")) { i <- which(tolower(nm) == v); if (length(i) == 1) nm[i] <- v }
  for (v in c("E", "R", "K", "LS", "C")) { i <- which(toupper(nm) == v); if (length(i) == 1) nm[i] <- v }
  names(d) <- nm
  miss <- setdiff(c("id", "E", GV), nm)
  if (length(miss)) stop("Fields not found in ", path, ": ", paste(miss, collapse = ", "))
  d$x <- xy[, 1]; d$y <- xy[, 2]
  if (is.null(d$stratum)) d$stratum <- findInterval(d$E, c(1, 5))   # 0: E < 1, 1: 1 <= E < 5, 2: E >= 5
  if (is.null(d$label))   d$label   <- as.integer(d$E >= 5)          # erosion class
  if (is.null(d$weight)) {                                           # cells represented by each point
    n <- table(d$stratum)
    d$weight <- as.numeric(CELLS[as.character(d$stratum)] / n[as.character(d$stratum)])
  }
  d <- d[order(d$id), ]; rownames(d) <- NULL
  stopifnot(nrow(d) == 15000, all(table(d$stratum) == 5000), !anyNA(d[, GV]))
  d
}

## Coding used for modelling: geomorphology 9 units -> 4 classes (Table S12), HSG ordinal,
## lithology, soil texture and LU/LC as classes, STI log-transformed.
GEO <- c("1" = "Active floodplain", "5" = "Active floodplain", "6" = "Active floodplain",
         "2" = "Older alluvial plain", "3" = "Gullied/ravinous area",
         "4" = "Abandoned channels and depressions", "7" = "Abandoned channels and depressions",
         "8" = "Abandoned channels and depressions", "9" = "Abandoned channels and depressions")
prep <- function(x) {
  x$GV12_Geomorphology <- factor(GEO[as.character(x$GV12_Geomorphology)],
                                 levels = c("Older alluvial plain", "Active floodplain",
                                            "Abandoned channels and depressions", "Gullied/ravinous area"))
  x$GV11_Lithology <- factor(x$GV11_Lithology, levels = 1:3)
  x$GV10_HSG <- as.numeric(x$GV10_HSG)
  if (!is.null(x$GV14_SoilTexture)) x$GV14_SoilTexture <- factor(x$GV14_SoilTexture, levels = 2:7)
  if (!is.null(x$GV17_LULC))        x$GV17_LULC        <- factor(x$GV17_LULC, levels = 3:9)
  if (!is.null(x$GV09_STI))         x$GV09_STI         <- log1p(x$GV09_STI)
  x
}

## ---- weighted performance measures (weights = cells represented by each point) ----
wauc <- function(s, y, w) {                     # ROC-AUC (weighted Mann-Whitney statistic)
  g  <- factor(s)
  wp <- tapply(w * (y == 1), g, sum); wn <- tapply(w * (y == 0), g, sum)
  wp[is.na(wp)] <- 0; wn[is.na(wn)] <- 0
  sum(wp * (cumsum(wn) - wn + 0.5 * wn)) / (sum(wp) * sum(wn))
}
wap <- function(s, y, w) {                      # PR-AUC (weighted average precision)
  o <- order(-s); s <- s[o]; y <- y[o]; w <- w[o]
  g <- cumsum(c(TRUE, diff(s) != 0))
  tp <- tapply(w * y, g, sum); fp <- tapply(w * (1 - y), g, sum)
  ctp <- cumsum(tp); cfp <- cumsum(fp)
  sum((tp / sum(tp)) * ctp / (ctp + cfp))
}
cls <- function(s, y, w, thr = 0.5) {           # measures at a probability threshold
  p1 <- s >= thr
  TP <- sum(w[p1 & y == 1]); FP <- sum(w[p1 & y == 0]); FN <- sum(w[!p1 & y == 1]); TN <- sum(w[!p1 & y == 0])
  sens <- TP / (TP + FN); spec <- TN / (TN + FP); prec <- TP / (TP + FP)
  c(sensitivity = sens, specificity = spec, balanced_accuracy = (sens + spec) / 2, precision = prec,
    F1 = 2 * prec * sens / (prec + sens),
    MCC = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
    TP = TP, FP = FP, FN = FN, TN = TN)
}

## ---- models ----
## SVM and MLPNN: indicator coding of classes, centring and scaling with the training data
scaler <- function(df, vars) {
  f <- as.formula(paste("~", paste(vars, collapse = " + ")))
  X <- model.matrix(f, df)[, -1, drop = FALSE]
  s <- apply(X, 2, sd); s[!is.finite(s) | s == 0] <- 1
  list(f = f, mu = colMeans(X), s = s)
}
apply_sc <- function(df, sc) {
  X <- model.matrix(sc$f, df)[, -1, drop = FALSE]
  sweep(sweep(X, 2, sc$mu, "-"), 2, sc$s, "/")
}
## fit one model; p = list of hyperparameters; trn needs the factor y ("none", "erosion")
fit_model <- function(al, p, trn, vars) {
  if (al == "C50")
    return(list(model = C5.0(x = trn[vars], y = trn$y, trials = p$trials,
                             control = C5.0Control(winnow = as.logical(p$winnow), minCases = p$minCases,
                                                   seed = SEED)), vars = vars))
  sc <- scaler(trn, vars); X <- apply_sc(trn, sc)
  if (al == "SVM") {
    m  <- ksvm(X, trn$y, kernel = "rbfdot", kpar = list(sigma = p$sigma), C = p$C, scaled = FALSE)
    dv <- as.numeric(predict(m, X, type = "decision"))
    return(list(model = m, scaler = sc, vars = vars,       # sgn orients the decision values to "erosion"
                sgn = if (mean(dv[trn$y == "erosion"]) >= mean(dv[trn$y == "none"])) 1 else -1))
  }
  set.seed(SEED)
  list(model = nnet(X, as.numeric(trn$y == "erosion"), size = p$size, decay = p$decay, maxit = 500,
                    entropy = TRUE, MaxNWts = 10000, trace = FALSE), scaler = sc, vars = vars)
}
## ranking score (C5.0 and MLPNN: probability; SVM: oriented decision value)
score <- function(o, al, nd) {
  if (al == "C50") return(predict(o$model, nd[o$vars], type = "prob")[, "erosion"])
  X <- apply_sc(nd, o$scaler)
  if (al == "SVM") return(o$sgn * as.numeric(predict(o$model, X, type = "decision")))
  as.numeric(predict(o$model, X))
}
## probability of the erosion class (SVM: Platt scaling of the decision values)
prob <- function(o, al, nd) {
  s <- score(o, al, nd)
  if (al == "SVM") plogis(coef(o$platt)[1] + coef(o$platt)[2] * s) else s
}
model_file <- function(design, set, algo, k) file.path(MODELS, sprintf("M1_%s_set%s_%s_fold%d.rds", design, set, algo, k))
