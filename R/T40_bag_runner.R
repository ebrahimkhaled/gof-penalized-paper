# Runs the T40 BAGofT job list with one R process per test (Amendment 3), NW at a time.
# A test whose process ends without writing its row is retried once; a second failure is recorded as a
# failure (NA p-value) and counted in the report. When every job is resolved the rows are merged into
# data/T40_bagoft_pvalues.csv. Restartable: finished rows on disk are never re-run.
#   Rscript T40_bag_runner.R        env T40_NW = concurrent processes (default 22)
suppressPackageStartupMessages(library(processx))
src <- readLines("T40_highdim_grid.R")
eval(parse(text = src[1:(which(startsWith(src, "if (MODE"))[1] - 1)]))
NW  <- as.integer(Sys.getenv("T40_NW", "22"))
dir <- file.path("..", "data", "T40_bag_rows"); dir.create(dir, showWarnings = FALSE)
row_file <- function(c, r) file.path(dir, sprintf("c%d_r%d.csv", c, r))

jobs <- do.call(rbind, lapply(cells$cell[cells$bag_reps > 0], function(c)
  data.frame(cell = c, rep = seq_len(cells$bag_reps[cells$cell == c]))))
jobs$tries <- 0L
todo <- which(!file.exists(row_file(jobs$cell, jobs$rep)))
cat(sprintf("%s  %d BAGofT tests, %d to run\n", format(Sys.time()), nrow(jobs), length(todo))); flush.console()

running <- list(); failed <- integer(0); T0 <- Sys.time()
while (length(todo) || length(running)) {
  while (length(todo) && length(running) < NW) {
    j <- todo[1]; todo <- todo[-1]; jobs$tries[j] <- jobs$tries[j] + 1L
    running[[as.character(j)]] <- process$new("Rscript", c("T40_bag_one.R", jobs$cell[j], jobs$rep[j]),
                                               wd = getwd(), stdout = NULL, stderr = NULL)
  }
  Sys.sleep(15)
  for (k in names(running)) if (!running[[k]]$is_alive()) {
    j <- as.integer(k); running[[k]] <- NULL
    if (!file.exists(row_file(jobs$cell[j], jobs$rep[j]))) {
      if (jobs$tries[j] < 2L) todo <- c(todo, j) else failed <- c(failed, j)
      cat(sprintf("%s  cell %d rep %d ended without a row (try %d)\n", format(Sys.time()),
                  jobs$cell[j], jobs$rep[j], jobs$tries[j])); flush.console()
    }
  }
  done <- sum(file.exists(row_file(jobs$cell, jobs$rep)))
  if (done %% 10 == 0) { cat(sprintf("[%7.1f min] %d / %d done, %d running, %d failed\n",
      as.numeric(difftime(Sys.time(), T0, units = "mins")), done, nrow(jobs), length(running), length(failed)))
    flush.console() }
}

rows <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(j) {
  f <- row_file(jobs$cell[j], jobs$rep[j])
  if (file.exists(f)) read.csv(f) else
    data.frame(cell = jobs$cell[j], rep = jobs$rep[j], fp_y = NA, fp_x = NA, p_bag = NA, p_bag_min = NA,
               stat_bag_mean = NA, t_bag = NA)
}))
rows$failed <- is.na(rows$fp_y)
write.csv(rows, file.path("..", "data", "T40_bagoft_pvalues.csv"), row.names = FALSE)
cat(sprintf("%s  merged %d rows (%d failed twice)\n", format(Sys.time()), nrow(rows), sum(rows$failed)))
