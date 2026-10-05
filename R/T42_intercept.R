# T42 -- T40's main panel with a non-zero intercept (PREDECLARATION_T40_highdim.md, Amendment 4).
# The T40 truth has no intercept, which happens to match PLStests' intercept-free refit. Here the truth has
# intercept alpha = -1.530973 (event fraction exactly 0.25 under the null, since eta0 ~ N(0, 1.5^2)), so every
# method meets a model with an intercept and a lower prevalence. Tests, functions and settings are T40's and
# T41's own (read from their scripts); only the data generator gains the intercept. Cells 201-212.
#   Rscript T42_intercept.R smoke | main        env T42_NW = workers (default 22)
src <- readLines("T40_highdim_grid.R")
eval(parse(text = src[1:(which(startsWith(src, "if (MODE"))[1] - 1)]))
s41 <- readLines("T41_mle_baseline.R")
eval(parse(text = s41[grep("^mle_rep <- function", s41):(which(startsWith(s41, "if (MODE"))[1] - 1)]))
MODE <- commandArgs(TRUE)[1]; if (is.na(MODE)) MODE <- "smoke"
NW <- as.integer(Sys.getenv("T42_NW", "22"))
ALPHA <- -1.530973
OUT42 <- file.path("..", "data", "T42_intercept_pvalues.csv")

cells <- do.call(rbind, lapply(seq_along(c(20L, 100L, 180L)), function(i) {
  p <- c(20L, 100L, 180L)[i]
  data.frame(cell = 200L + 4L * (i - 1L) + 1:4, panel = "I", p = p, rho = 0.7, signal = "dense",
             dep = c("null", "index", "coord", "inter"), a = c(0, 1, 1, 1),
             reps = c(1000L, 500L, 500L, 500L), bag_reps = 0L)
}))

make_data <- function(ce, r) {                  # T40's generator with the intercept added to the truth
  set.seed(ce$cell * 1e6 + r)
  p <- ce$p; S <- ce$rho^abs(outer(seq_len(p), seq_len(p), "-"))
  b <- if (ce$signal == "dense") rep(c(0.35, -0.30, 0.25, -0.20, 0.15), length.out = p)
       else c(rep(1, 5), rep(0, p - 5))
  b <- b * GAMMA / sqrt(drop(t(b) %*% S %*% b))
  X <- matrix(rnorm(N * p), N, p) %*% chol(S); colnames(X) <- paste0("x", seq_len(p))
  e0 <- drop(X %*% b)
  dep <- switch(ce$dep,
    null  = 0,
    index = ((e0 / GAMMA)^2 - 1) / sqrt(2),
    coord = (X[, 1]^2 - 1) / sqrt(2),
    inter = X[, 1] * X[, 2] / sqrt(1 + ce$rho^2))
  y <- rbinom(N, 1, plogis(ALPHA + e0 + ce$a * dep))
  list(X = X, y = y, fp_y = sum(y), fp_x = round(sum(X[1:3, 1:3]), 10))
}

# T40's run_jobs ships only T40's names to the workers; T42 also needs ALPHA, mle_rep and both_rep there.
.run_jobs_t40 <- run_jobs
run_jobs <- function(jobs, fun, path, nw) {
  done <- done_pairs(path)
  jobs <- jobs[!paste(jobs$cell, jobs$rep) %in% done, , drop = FALSE]
  cat(sprintf("%s: %d jobs to run (%d already on disk)\n", basename(path), nrow(jobs), length(done)))
  if (!nrow(jobs)) return(invisible())
  cl <- makeCluster(nw); on.exit(stopCluster(cl))
  clusterExport(cl, c("cells", "make_data", "ridge_fit", "pieces", "mc_edge_p", "main_rep", "mle_rep",
                      "both_rep", "N", "G", "NB", "GAMMA", "ALPHA"), envir = globalenv())
  T0 <- Sys.time(); chunk <- 4L * nw
  for (s in seq(1L, nrow(jobs), by = chunk)) {
    J <- jobs[s:min(s + chunk - 1L, nrow(jobs)), ]
    R <- parLapplyLB(cl, seq_len(nrow(J)), function(i, J, job_fun) {
      ce <- cells[cells$cell == J$cell[i], ]; job_fun(ce, J$rep[i]) }, J = J, job_fun = fun)
    append_rows(do.call(rbind, R), path)
    cat(sprintf("[%7.1f min] %s %d / %d\n", as.numeric(difftime(Sys.time(), T0, units = "mins")),
                basename(path), min(s + chunk - 1L, nrow(jobs)), nrow(jobs))); flush.console()
  }
}

both_rep <- function(ce, r) {                    # T40's tests and T41's MLE baseline on one dataset
  a <- main_rep(ce, r); b <- mle_rep(ce, r)
  stopifnot(a$fp_y == b$fp_y)
  cbind(a, b[, setdiff(names(b), c("cell", "rep", "fp_y", "fp_x"))])
}

if (MODE == "smoke") {
  for (c in c(201L, 205L, 209L)) { t0 <- Sys.time(); x <- both_rep(cells[cells$cell == c, ], 1L)
    cat(sprintf("cell %d p=%d | events %d | HL %.3f SC.HL %.3f SC.EDGE %.3f GRP5 %.3f PLS %.3f | MLE %s HL_mle %.3f | %.0fs\n",
        c, cells$p[cells$cell == c], x$fp_y, x$p_hl, x$p_schl, x$p_scedge, x$p_grp5, x$p_pls, x$mle_exists,
        x$p_hl_mle, as.numeric(difftime(Sys.time(), t0, units = "secs")))) }
  ev <- mean(sapply(1:200, function(r) make_data(cells[1, ], r)$fp_y)) / N
  cat(sprintf("event fraction over 200 null datasets: %.3f (target 0.25)\n", ev))
  tmp <- tempfile(fileext = ".csv")
  run_jobs(data.frame(cell = c(201L, 202L), rep = c(1L, 1L)), both_rep, tmp, 2L)
  run_jobs(data.frame(cell = c(201L, 202L), rep = c(1L, 1L)), both_rep, tmp, 2L)
  stopifnot(nrow(read.csv(tmp)) == 2)
  cat("smoke OK\n")
} else if (MODE == "main") {
  jobs <- do.call(rbind, lapply(cells$cell, function(c) data.frame(cell = c, rep = seq_len(cells$reps[cells$cell == c]))))
  run_jobs(jobs, both_rep, OUT42, NW)
} else if (MODE == "recheck") {
  # Amendment 6: independent second block of the p/n = 0.25 null cell (205), replicates 1001..2000,
  # seeds the main run never used
  out <- Sys.getenv("T42_OUT", file.path("..", "data", "T42_recheck_pvalues.csv"))
  nrep <- as.integer(Sys.getenv("T42_REPS", "1000"))
  run_jobs(data.frame(cell = 205L, rep = 1000L + seq_len(nrep)), both_rep, out, NW)
}
