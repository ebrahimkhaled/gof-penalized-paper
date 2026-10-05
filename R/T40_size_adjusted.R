# Size-adjusted power for T40 (the convention of SI Table S5): for each method and each departure cell, the
# threshold is the 5th percentile (type 1) of that method's p-values in the matched null cell (same panel,
# p, rho, signal), and power is the share of departure p-values at or below it. Every test is then compared
# at exactly 5% size. Only complete cells, and only departures whose null cell is complete.
#   Rscript T40_size_adjusted.R  -> data/T40_size_adjusted.csv
src <- readLines("T40_highdim_grid.R")
eval(parse(text = src[1:(which(startsWith(src, "if (MODE"))[1] - 1)]))
m <- read.csv(file.path("..", "data", "T40_main_pvalues.csv"))
meth <- c(HL = "p_hl", EDGE_u = "p_edge_u", SC.HL = "p_schl", SC.EDGE = "p_scedge",
          GRP5 = "p_grp5", GRP1 = "p_grp1", PLS = "p_pls")
complete <- function(c) sum(m$cell == c) >= cells$reps[cells$cell == c]
out <- NULL
for (c in cells$cell[cells$dep != "null"]) {
  ce <- cells[cells$cell == c, ]
  c0 <- cells$cell[cells$panel == ce$panel & cells$p == ce$p & cells$rho == ce$rho &
                   cells$signal == ce$signal & cells$dep == "null"]
  if (!complete(c) || !complete(c0)) next
  z <- m[m$cell == c, ]; z0 <- m[m$cell == c0, ]
  row <- data.frame(ce[, c("cell", "panel", "p", "rho", "signal", "dep", "a")], kappa = ce$p / N)
  for (nm in names(meth)) {
    v <- z[[meth[[nm]]]]; v0 <- z0[[meth[[nm]]]]; v <- v[!is.na(v)]; v0 <- v0[!is.na(v0)]
    thr <- unname(quantile(v0, 0.05, type = 1))
    pw <- mean(v <= thr)
    row[[nm]] <- pw; row[[paste0("se_", nm)]] <- sqrt(pw * (1 - pw) / length(v))
    row[[paste0("nom_", nm)]] <- mean(v < 0.05); row[[paste0("thr_", nm)]] <- thr
  }
  out <- rbind(out, row)
}
write.csv(out, file.path("..", "data", "T40_size_adjusted.csv"), row.names = FALSE)
show <- out[, c("panel", "p", "kappa", "dep", "a", "SC.HL", "SC.EDGE", "GRP5", "GRP1", "PLS", "HL", "EDGE_u")]
print(show, row.names = FALSE, digits = 2)
