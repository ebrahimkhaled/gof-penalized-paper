# Why did the pilot's BAGofT workers die? One BAGofT test per p (20, 100, 180), serially, in one plain R
# process: time, peak memory, and any R error. A segfault kills the process; the log then shows the cell
# that was running. Not a result file: output goes to data/T40_bag_diagnose.log only.
src <- readLines("T40_highdim_grid.R")
eval(parse(text = src[1:(which(startsWith(src, "if (MODE"))[1] - 1)]))
suppressPackageStartupMessages(library(BAGofT))
for (c in c(1L, 15L, 29L)) {
  ce <- cells[cells$cell == c, ]
  cat(sprintf("%s  start cell %d (p = %d)\n", format(Sys.time(), "%H:%M:%S"), c, ce$p)); flush.console()
  gc(reset = TRUE)
  t0 <- proc.time()[[3]]
  res <- tryCatch(bag_rep(ce, 1L), error = function(e) { cat("R ERROR:", conditionMessage(e), "\n"); NULL })
  g <- gc()
  cat(sprintf("%s  end   cell %d: %.1f min, peak R memory %.0f MB, p.value %s\n", format(Sys.time(), "%H:%M:%S"),
              c, (proc.time()[[3]] - t0) / 60, sum(g[, 6]), if (is.null(res)) "NA" else format(res$p_bag)))
  flush.console()
}
