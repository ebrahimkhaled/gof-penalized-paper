# Summary of T40 (with T41 and BAGofT beside it): one row per COMPLETE cell, rejection rate at p < 0.05 for
# every method, with its Monte Carlo SE. Incomplete cells are skipped, never reported early.
#   Rscript T40_summarise.R            writes data/T40_summary.csv and prints the size and power tables
src <- readLines("T40_highdim_grid.R")
eval(parse(text = src[1:(which(startsWith(src, "if (MODE"))[1] - 1)]))
m <- read.csv(file.path("..", "data", "T40_main_pvalues.csv"))
stopifnot(!anyDuplicated(paste(m$cell, m$rep)))
mle <- if (file.exists(f <- file.path("..", "data", "T41_mle_pvalues.csv"))) read.csv(f) else NULL
bag <- if (file.exists(f <- file.path("..", "data", "T40_bagoft_pvalues.csv"))) read.csv(f) else NULL
if (!is.null(mle)) {                                                   # pairing check with T41
  k <- merge(m[, c("cell", "rep", "fp_y", "fp_x")], mle[, c("cell", "rep", "fp_y", "fp_x")], by = c("cell", "rep"))
  stopifnot(all(k$fp_y.x == k$fp_y.y), all(abs(k$fp_x.x - k$fp_x.y) < 1e-8))
}
rate <- function(v) { v <- v[!is.na(v)]; if (!length(v)) return(c(NA, NA, 0))
  r <- mean(v < 0.05); c(r, sqrt(r * (1 - r) / length(v)), length(v)) }
meth <- c(HL = "p_hl", EDGE_u = "p_edge_u", SC.HL = "p_schl", SC.EDGE = "p_scedge",
          GRP5 = "p_grp5", GRP1 = "p_grp1", PLS = "p_pls")
out <- NULL
for (c in sort(unique(m$cell))) {
  ce <- cells[cells$cell == c, ]; z <- m[m$cell == c, ]
  if (nrow(z) < ce$reps) next
  row <- data.frame(ce[, c("cell", "panel", "p", "rho", "signal", "dep", "a")], kappa = ce$p / N, n = nrow(z))
  for (nm in names(meth)) { r <- rate(z[[meth[[nm]]]]); row[[nm]] <- r[1]; row[[paste0("se_", nm)]] <- r[2] }
  row$grp_NA <- sum(is.na(z$p_grp5)); row$pls_NA <- sum(is.na(z$p_pls))
  if (!is.null(mle)) { e <- mle[mle$cell == c & mle$mle_exists, ]
    row$mle_exists <- mean(mle$mle_exists[mle$cell == c]); row$HL_mle <- rate(e$p_hl_mle)[1]
    row$EDGE_mle <- rate(e$p_edge_mle)[1] }
  if (!is.null(bag) && any(bag$cell == c)) { b <- bag[bag$cell == c, ]
    row$BAG <- rate(b$p_bag)[1]; row$se_BAG <- rate(b$p_bag)[2]; row$BAG_n <- sum(!is.na(b$p_bag)) }
  row$t_ours <- mean(z$t_ours); row$t_grp5 <- mean(z$t_grp5); row$t_pls <- mean(z$t_pls)
  out <- rbind(out, row)
}
write.csv(out, file.path("..", "data", "T40_summary.csv"), row.names = FALSE)
cols <- intersect(c("panel", "rho", "signal", "p", "dep", "a", "HL", "EDGE_u", "SC.HL", "se_SC.HL", "SC.EDGE",
                    "se_SC.EDGE", "GRP5", "GRP1", "PLS", "mle_exists", "HL_mle", "BAG"), names(out))
cat("\n== SIZE (null cells complete so far) ==\n")
print(out[out$dep == "null", cols], row.names = FALSE, digits = 3)
cat("\n== POWER (departure cells complete so far) ==\n")
print(out[out$dep != "null", cols], row.names = FALSE, digits = 3)
cat(sprintf("\n%d of %d cells complete; %d of %d replicates on disk\n", nrow(out), nrow(cells), nrow(m), sum(cells$reps)))
