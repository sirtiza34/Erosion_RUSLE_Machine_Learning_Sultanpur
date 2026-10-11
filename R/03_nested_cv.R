## =============================================================================================
##  03_nested_cv.R  --  nested, spatially blocked cross-validation of C5.0 DT, SVM and MLPNN
##  Designs: strict (trained on E < 1 and E >= 5) and no-gap (trained on all points);
##  both are tested on all held-out points. Sets: A (17 variables) and B (9 variables).
##  Each outer fold: training points within 1.2 km of the held-out fold removed; hyperparameters
##  tuned by three-fold spatial CV inside the training data (mean weighted ROC-AUC); the selected
##  setting refitted on the whole training data; SVM probabilities by Platt scaling of the inner
##  out-of-fold decision values.
##  Outputs: M1_oof_predictions_<design>.csv, M1_tuning_<design>.csv (Table S18), fold models.
##  Runtime: several hours (about 1 h strict, 3 h no-gap on 3 cores). Finished tasks are kept in
##  results/models/ and skipped when the script is run again.
## =============================================================================================
source("R/00_setup.R")
DESIGNS <- c("strict", "nogap")

d <- prep(merge(read_samples(), read.csv(FOLDS), by = "id"))
d$y <- factor(ifelse(d$label == 1, "erosion", "none"), levels = c("none", "erosion"))

## ---- hyperparameter search (Table S17) ----
grid_c50 <- expand.grid(trials = c(10, 30, 60, 100), winnow = c(0, 1), minCases = c(2, 20, 50))
grid_mlp <- expand.grid(size = c(5, 10, 15, 20), decay = c(0.01, 0.1, 0.3, 1))
C_svm    <- c(0.25, 1, 4, 16, 64)
SIG_MULT <- c(0.25, 0.5, 1, 2, 4)          # multiples of the median of sigest() on the training data

## inner spatial folds inside the training data of outer fold k
inner_folds <- function(trn, k) {
  pts <- st_as_sf(trn, coords = c("x", "y"), crs = CRS_UTM)
  cv_spatial(x = pts, column = "label", k = K_IN, hexagon = FALSE, size = BLOCK, selection = "random",
             iteration = 100, balance = TRUE, seed = SEED + k, plot = FALSE, report = FALSE,
             progress = FALSE)$folds_ids
}

## tune one algorithm on the inner folds, refit with the selected setting
run_algo <- function(al, trn, vars, k) {
  grid <- switch(al, C50 = grid_c50, MLP = grid_mlp, SVM = {
    set.seed(SEED + k)
    sg <- sigest(apply_sc(trn, scaler(trn, vars)), scaled = FALSE)[2]
    expand.grid(sigma = as.numeric(sg) * SIG_MULT, C = C_svm)
  })
  inn <- trn$inner
  oof <- matrix(NA_real_, nrow(trn), nrow(grid))
  for (g in seq_len(nrow(grid))) for (j in 1:K_IN) {
    o <- fit_model(al, as.list(grid[g, ]), trn[inn != j, ], vars)
    oof[inn == j, g] <- score(o, al, trn[inn == j, ])
  }
  auc <- apply(oof, 2, function(s) mean(sapply(1:K_IN, function(j)
    wauc(s[inn == j], trn$label[inn == j], trn$weight[inn == j]))))
  b <- which.max(auc); p <- as.list(grid[b, ])
  o <- fit_model(al, p, trn, vars)
  if (al == "SVM") o$platt <- glm(y ~ dv, family = binomial, data = data.frame(y = trn$label, dv = oof[, b]))
  list(model = o, params = paste(names(p), unlist(p), sep = "=", collapse = "; "), auc = auc[b])
}

task_file <- function(ds, st, k) file.path(MODELS, sprintf("M1_task_%s_set%s_fold%d.rds", ds, st, k))

## one task = one design, one variable set and one outer fold (all three algorithms)
run_task <- function(k, st, ds) {
  vars <- SETS[[st]]
  tr <- d$fold != k & d[[paste0("buf_excl_f", k)]] == 0
  if (ds == "strict") tr <- tr & d$stratum != 1
  trn <- d[tr, ]; tst <- d[d$fold == k, ]
  trn$inner <- inner_folds(trn, k)
  pred <- list(); tune <- list(); errs <- character()
  for (al in ALGOS) {
    r <- tryCatch(run_algo(al, trn, vars, k), error = function(e) e)
    if (inherits(r, "error")) {
      errs <- c(errs, sprintf("%s set %s fold %d %s: %s", ds, st, k, al, conditionMessage(r))); next
    }
    saveRDS(r$model, model_file(ds, st, al, k))
    pred[[al]] <- data.frame(id = tst$id, set = st, algo = al, fold = k, prob = prob(r$model, al, tst))
    tune[[al]] <- data.frame(fold = k, set = st, algo = al, n_train = nrow(trn), params = r$params,
                             inner_wAUC = round(r$auc, 4))
  }
  out <- list(pred = do.call(rbind, pred), tune = do.call(rbind, tune), errors = errs)
  if (!length(errs)) saveRDS(out, task_file(ds, st, k))
  out
}

## ---- run the tasks in parallel ----
tasks <- expand.grid(k = 1:K, set = names(SETS), design = DESIGNS, stringsAsFactors = FALSE)
tfile <- task_file(tasks$design, tasks$set, tasks$k)
todo  <- tasks[!file.exists(tfile), ]
cat("Start:", format(Sys.time()), "| tasks to run:", nrow(todo), "of", nrow(tasks), "\n")
if (nrow(todo)) {
  cl <- parallel::makeCluster(min(N_CORES, nrow(todo)))
  parallel::clusterCall(cl, setwd, getwd())
  invisible(parallel::clusterEvalQ(cl, source("R/00_setup.R")))
  parallel::clusterExport(cl, c("d", "grid_c50", "grid_mlp", "C_svm", "SIG_MULT", "inner_folds",
                                "run_algo", "task_file", "run_task"), envir = environment())
  res <- parallel::parLapply(cl, split(todo, seq_len(nrow(todo))), function(t) run_task(t$k, t$set, t$design))
  parallel::stopCluster(cl)
  errs <- unlist(lapply(res, `[[`, "errors"))
  if (length(errs)) { cat("Errors (these tasks will be run again next time):\n"); print(errs) }
}
cat("End:", format(Sys.time()), "\n")

## ---- combine, save and summarise each design ----
for (ds in DESIGNS) {
  f <- tfile[tasks$design == ds]; f <- f[file.exists(f)]
  if (!length(f)) next
  r <- lapply(f, readRDS)
  pred <- do.call(rbind, lapply(r, `[[`, "pred")); tune <- do.call(rbind, lapply(r, `[[`, "tune"))
  cat("\n==========", ds, "design:", length(f), "of", K * length(SETS), "tasks,", nrow(pred), "predictions\n")
  pred <- merge(pred, d[, c("id", "E", "stratum", "label", "weight")], by = "id")
  write.csv(pred, file.path(OUT, sprintf("M1_oof_predictions_%s.csv", ds)), row.names = FALSE)
  write.csv(tune, file.path(OUT, sprintf("M1_tuning_%s.csv", ds)), row.names = FALSE)
  q <- do.call(rbind, lapply(split(pred, list(pred$set, pred$algo, pred$fold), drop = TRUE), function(p)
    data.frame(set = p$set[1], algo = p$algo[1], fold = p$fold[1],
               AUC_all = wauc(p$prob, p$label, p$weight),
               AUC_extremes = with(p[p$stratum != 1, ], wauc(prob, label, weight)))))
  cat("Weighted ROC-AUC on the held-out spatial folds (mean, min, max over the folds):\n")
  print(aggregate(cbind(AUC_all, AUC_extremes) ~ set + algo, q,
                  function(v) round(c(mean = mean(v), min = min(v), max = max(v)), 3)))
  cat("Selected settings per fold:\n")
  print(tune[order(tune$set, tune$algo, tune$fold), ], row.names = FALSE)
}
