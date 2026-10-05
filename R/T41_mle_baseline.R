# T41 -- the maximum-likelihood baseline on the T40 datasets (PREDECLARATION_T40_highdim.md, Amendment 2).
# "Why not skip the penalty and use the classical test?" Same cells, same seeds, same data as T40: the
# definitions (cells, make_data, pieces, mc_edge_p, run_jobs ...) are read from T40_highdim_grid.R itself,
# up to its driver, so the two studies cannot drift apart. Paired with T40 through (cell, rep) and the
# data fingerprint.
#
# Usage (from R/):  Rscript T41_mle_baseline.R smoke | main      env T41_NW = workers (default 2)
src <- readLines("T40_highdim_grid.R")
eval(parse(text = src[1:(which(startsWith(src, "if (MODE"))[1] - 1)]))
MODE <- commandArgs(TRUE)[1]; if (is.na(MODE)) MODE <- "smoke"
NW <- as.integer(Sys.getenv("T41_NW", "2"))
OUT_MLE <- file.path("..", "data", "T41_mle_pvalues.csv")

mle_rep <- function(ce, r) {
  d <- make_data(ce, r); X1 <- cbind(1, d$X); y <- d$y
  fit <- suppressWarnings(tryCatch(glm.fit(X1, y, family = binomial(), control = list(maxit = 100)),
                                   error = function(e) NULL))
  na_row <- data.frame(cell = ce$cell, rep = r, fp_y = d$fp_y, fp_x = d$fp_x, mle_exists = FALSE,
                       converged = NA, min_pi = NA_real_, max_pi = NA_real_, max_coef = NA_real_,
                       p_hl_mle = NA_real_, p_edge_mle = NA_real_)
  if (is.null(fit)) return(na_row)
  pi <- fit$fitted.values; mc <- max(abs(fit$coefficients), na.rm = TRUE)
  exists <- isTRUE(fit$converged) && all(pi > 1e-8 & pi < 1 - 1e-8) && mc < 50
  # lambda = 0 makes pieces() the classical (Chernoff-Lehmann / Moore-Spruill) quantities on the MLE
  st <- pieces(X1, y, list(beta = fit$coefficients, pi = pi), 0, diag(ncol(X1)))
  hl <- sum((st$o - st$e)^2 / (st$e * (1 - st$e / st$n_g)))
  data.frame(cell = ce$cell, rep = r, fp_y = d$fp_y, fp_x = d$fp_x, mle_exists = exists,
             converged = fit$converged, min_pi = min(pi), max_pi = max(pi), max_coef = mc,
             p_hl_mle = pchisq(hl, G - 2, lower.tail = FALSE),
             p_edge_mle = tryCatch(mc_edge_p(st$edge_u, st$Om, st$Z), error = function(e) NA_real_))
}

if (MODE == "smoke") {
  for (c in c(1L, 15L, 29L, 36L, 76L)) {
    x <- mle_rep(cells[cells$cell == c, ], 1L)
    cat(sprintf("cell %2d p=%3d | exists %s conv %s max|b| %.1f | HL_mle %.3f EDGE_mle %.3f\n", c,
                cells$p[cells$cell == c], x$mle_exists, x$converged, x$max_coef, x$p_hl_mle, x$p_edge_mle))
  }
  # pairing with T40: the fingerprint of replicate 1 of cell 1 must equal T40's
  t40 <- file.path("..", "data", "T40_main_pvalues.csv")
  if (file.exists(t40)) {
    a <- read.csv(t40); a <- a[a$cell == 1 & a$rep == 1, ]
    b <- mle_rep(cells[cells$cell == 1, ], 1L)
    stopifnot(nrow(a) == 1, a$fp_y == b$fp_y, abs(a$fp_x - b$fp_x) < 1e-8)
    cat("paired with T40: fingerprints match\n")
  }
  tmp <- tempfile(fileext = ".csv")
  run_jobs(data.frame(cell = c(1L, 29L), rep = c(1L, 1L)), mle_rep, tmp, 2L)
  run_jobs(data.frame(cell = c(1L, 29L), rep = c(1L, 1L)), mle_rep, tmp, 2L)
  stopifnot(nrow(read.csv(tmp)) == 2)
  cat("smoke OK\n")
} else if (MODE == "main") {
  jobs <- do.call(rbind, lapply(cells$cell, function(c) data.frame(cell = c, rep = seq_len(cells$reps[cells$cell == c]))))
  run_jobs(jobs, mle_rep, OUT_MLE, NW)
}
