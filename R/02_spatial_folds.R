## =============================================================================================
##  02_spatial_folds.R  --  5 km spatial blocks in five folds and the 1.2 km buffer
##  Writes data/spatial_folds_5km.csv (id, fold, block, buf_excl_f1..f5) and the block polygons.
##  If the deposited folds file already exists, it is kept and the new folds are compared with it.
##  Outputs: Table S16 (printed), spatial_folds_5km.png (Fig. S6)
## =============================================================================================
source("R/00_setup.R")

d   <- read_samples()
pts <- st_as_sf(d, coords = c("x", "y"), crs = CRS_UTM, remove = FALSE)
sb  <- cv_spatial(x = pts, column = "label", k = K, hexagon = FALSE, size = BLOCK,
                  selection = "random", iteration = 100, balance = TRUE, seed = SEED,
                  plot = FALSE, report = TRUE, progress = FALSE)
f <- data.frame(id = d$id, fold = sb$folds_ids,
                block = sapply(st_intersects(pts, sb$blocks), function(i) i[1]))
cat("\nBlocks:", nrow(sb$blocks), "\n\nPoints per fold and stratum:\n")
print(table(fold = f$fold, stratum = d$stratum))

## buffer: training points within BUF of any held-out point are removed
xy <- as.matrix(d[, c("x", "y")])
for (k in 1:K) {
  te <- which(f$fold == k); tr <- which(f$fold != k)
  nn <- get.knnx(xy[te, ], xy[tr, ], k = 1)$nn.dist[, 1]
  ex <- integer(nrow(d)); ex[tr[nn < BUF]] <- 1L
  f[[paste0("buf_excl_f", k)]] <- ex
  keep <- f$fold != k & ex == 0
  cat(sprintf("fold %d: held-out %d | removed by buffer %d (%.1f%%) | strict training: non-erosion %d, erosion %d | no-gap training %d\n",
              k, length(te), sum(ex), 100 * sum(ex) / length(tr),
              sum(keep & d$stratum == 0), sum(keep & d$stratum == 2), sum(keep)))
}

## save (the deposited folds are the reference: random block assignment can differ between versions)
if (file.exists(FOLDS)) {
  old  <- read.csv(FOLDS)
  same <- identical(as.integer(old$fold[match(f$id, old$id)]), as.integer(f$fold))
  cat("\nDeposited folds", if (same) "reproduced exactly" else "NOT reproduced; the deposited folds are used", "\n")
  write.csv(f, file.path(OUT, "spatial_folds_5km_regenerated.csv"), row.names = FALSE)
  if (same && is.null(old$block)) {             # add the block of each point to the deposited file
    write.csv(f, FOLDS, row.names = FALSE)
    st_write(sb$blocks, BLOCKS, delete_dsn = TRUE, quiet = TRUE)
  }
} else {
  write.csv(f, FOLDS, row.names = FALSE)
  st_write(sb$blocks, BLOCKS, delete_dsn = TRUE, quiet = TRUE)
  cat("\nSaved", FOLDS, "and", BLOCKS, "\n")
}

png(file.path(OUT, "spatial_folds_5km.png"), 2400, 2000, res = 300)
plot(st_geometry(sb$blocks), border = "grey75")
cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
points(d$x, d$y, pch = 16, cex = 0.25, col = cols[f$fold])
legend("bottomleft", legend = paste("Fold", 1:K), col = cols, pch = 16, bty = "n", cex = 0.8)
invisible(dev.off())
