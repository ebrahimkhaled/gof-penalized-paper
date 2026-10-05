# One BAGofT test of T40 in its own R process (Amendment 3): Rscript T40_bag_one.R <cell> <rep>
# Writes data/T40_bag_rows/c<cell>_r<rep>.csv through a temporary file and a rename, so a crash leaves no row
# (the runner then retries) and a finished row is never half-written.
src <- readLines("T40_highdim_grid.R")
eval(parse(text = src[1:(which(startsWith(src, "if (MODE"))[1] - 1)]))
a <- as.integer(commandArgs(TRUE)); cell <- a[1]; rep <- a[2]
dir <- file.path("..", "data", "T40_bag_rows"); dir.create(dir, showWarnings = FALSE)
out <- file.path(dir, sprintf("c%d_r%d.csv", cell, rep))
if (!file.exists(out)) {
  row <- bag_rep(cells[cells$cell == cell, ], rep)
  tmp <- paste0(out, ".tmp")
  write.csv(row, tmp, row.names = FALSE)
  file.rename(tmp, out)
}
