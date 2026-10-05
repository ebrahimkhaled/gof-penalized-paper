# T40 -- the corrected test across p/n, against packaged high-dimensional rivals.
# Protocol: PREDECLARATION_T40_highdim.md (written before any data). Paired design: every method sees the
# same dataset in each replicate. Replicate r of cell c is generated after set.seed(c * 1e6 + r), so results do
# not depend on the number of workers, and BAGofT (run as its own job list because it costs ~50 min per test)
# regenerates exactly the same data; both output files carry a data fingerprint to prove it.
# Every replicate is appended to disk as it finishes; a restart skips what is already there.
#
# Usage (from R/):  Rscript T40_highdim_grid.R smoke | pilot | main | bag
#   env T40_NW = workers (default 22); T40_ONLY = comma list of cell ids to restrict a run.
suppressPackageStartupMessages({library(glmnet); library(parallel)})

MODE <- commandArgs(TRUE)[1]; if (is.na(MODE)) MODE <- "smoke"
NW   <- as.integer(Sys.getenv("T40_NW", "22"))
# T40_OUT redirects a run to a scratch file, so a replicator can recompute a few cells (T40_ONLY) without
# touching the archived results, which a run would otherwise skip as already done
OUT_MAIN <- Sys.getenv("T40_OUT", file.path("..", "data", "T40_main_pvalues.csv"))
OUT_BAG  <- file.path("..", "data", "T40_bagoft_pvalues.csv")
N <- 400L; G <- 10L; NB <- 499L; GAMMA <- 1.5

# ---- the cells of the pre-declaration (Section 2) ---------------------------------------------------
cells <- local({
  ps <- c(20L, 40L, 100L, 140L, 180L); out <- list(); id <- 0L
  add <- function(panel, p, rho, signal, dep, a, reps, bag) {
    id <<- id + 1L
    # Amendment 1 (2026-10-02): BAGofT runs last, on p/n = 0.05, 0.25, 0.45 only, for the null and the
    # index- and coordinate-aligned departures at a = 1.0, 100 replicates each, at its published defaults
    bag_cell <- bag && p %in% c(20L, 100L, 180L) && (dep == "null" || (dep %in% c("index", "coord") && a == 1.0))
    out[[id]] <<- data.frame(cell = id, panel = panel, p = p, rho = rho, signal = signal, dep = dep, a = a,
                             reps = reps, bag_reps = if (bag_cell) 100L else 0L)
  }
  for (p in ps) {                                                     # panel M
    add("M", p, 0.7, "dense", "null", 0, 1000L, TRUE)
    for (d in c("index", "coord", "inter")) for (a in c(0.5, 1.0)) add("M", p, 0.7, "dense", d, a, 500L, TRUE)
  }
  for (rho in c(0.4, 0.8)) for (p in ps) {                              # panel R
    add("R", p, rho, "dense", "null", 0, 1000L, FALSE)
    for (d in c("index", "coord", "inter")) add("R", p, rho, "dense", d, 1.0, 500L, FALSE)
  }
  for (p in ps) {                                                     # panel S
    add("S", p, 0.7, "sparse", "null", 0, 1000L, FALSE)
    for (d in c("index", "coord", "inter")) add("S", p, 0.7, "sparse", d, 1.0, 500L, FALSE)
  }
  do.call(rbind, out)
})

# ---- data -------------------------------------------------------------------------------------------
make_data <- function(ce, r) {
  set.seed(ce$cell * 1e6 + r)
  p <- ce$p; S <- ce$rho^abs(outer(seq_len(p), seq_len(p), "-"))
  b <- if (ce$signal == "dense") rep(c(0.35, -0.30, 0.25, -0.20, 0.15), length.out = p)
       else c(rep(1, 5), rep(0, p - 5))
  b <- b * GAMMA / sqrt(drop(t(b) %*% S %*% b))                       # sd(eta0) = 1.5 exactly
  X <- matrix(rnorm(N * p), N, p) %*% chol(S); colnames(X) <- paste0("x", seq_len(p))
  e0 <- drop(X %*% b)
  dep <- switch(ce$dep,
    null  = 0,
    index = ((e0 / GAMMA)^2 - 1) / sqrt(2),
    coord = (X[, 1]^2 - 1) / sqrt(2),
    inter = X[, 1] * X[, 2] / sqrt(1 + ce$rho^2))
  y <- rbinom(N, 1, plogis(e0 + ce$a * dep))
  list(X = X, y = y, fp_y = sum(y), fp_x = round(sum(X[1:3, 1:3]), 10))
}

# ---- the paper's tests (functions as in R/T24_rivals.R and R/T31_cv_lambda_size.R) -------------------
ridge_fit <- function(X1, y, lam, D, tol = 1e-9, maxit = 100) {
  beta <- rep(0, ncol(X1))
  for (it in seq_len(maxit)) {
    eta <- drop(X1 %*% beta); pi <- 1 / (1 + exp(-eta))
    w <- pmax(pi * (1 - pi), 1e-8); z <- eta + (y - pi) / w
    bn <- drop(solve(crossprod(X1, w * X1) + lam * D, crossprod(X1, w * z)))
    if (max(abs(bn - beta)) < tol) { beta <- bn; break }
    beta <- bn
  }
  list(beta = beta, pi = 1 / (1 + exp(-drop(X1 %*% beta))))
}
pieces <- function(X1, y, fit, lam, D) {
  pi <- pmin(pmax(fit$pi, 1e-6), 1 - 1e-6); w <- pmax(pi * (1 - pi), 1e-8)
  g  <- pmin(ceiling(rank(pi, ties.method = "first") / (length(y) / G)), G)
  idx <- split(seq_along(y), g)
  Vg <- vapply(idx, function(I) sum(w[I]), 0.0)
  r  <- (vapply(idx, function(I) sum(y[I]), 0.0) - vapply(idx, function(I) sum(pi[I]), 0.0)) / sqrt(Vg)
  U  <- t(vapply(idx, function(I) colSums(w[I] * X1[I, , drop = FALSE]), numeric(ncol(X1)))) / sqrt(Vg)
  Fm <- crossprod(X1, w * X1); K <- lam * D
  bt <- fit$beta + drop(solve(Fm, K %*% fit$beta))
  mu <- drop(U %*% solve(Fm + K, K %*% bt))
  rc <- r - mu
  pbar <- vapply(idx, function(I) mean(pi[I]), 0.0)
  Z <- tryCatch(as.matrix(stats::poly(pbar, 3)), error = function(e) NULL)
  qf <- function(v) if (is.null(Z)) NA_real_ else {
    zr <- crossprod(Z, v); as.numeric(t(zr) %*% solve(crossprod(Z)) %*% zr) }
  list(dec_c = sum(rc^2), edge_c = qf(rc), dec_u = sum(r^2), edge_u = qf(r), r = r, Z = Z, bt = bt,
       Om = diag(G) - U %*% solve(Fm, t(U)), o = vapply(idx, function(I) sum(y[I]), 0.0),
       e = vapply(idx, function(I) sum(pi[I]), 0.0), n_g = lengths(idx))
}
mc_edge_p <- function(S, Om, Z, nd = 20000) {
  if (is.na(S) || is.null(Z)) return(NA_real_)
  R <- chol(Om + diag(1e-9, nrow(Om))); W <- matrix(rnorm(nd * nrow(Om)), nd) %*% R
  A <- solve(crossprod(Z)); Q <- W %*% Z
  mean(rowSums((Q %*% A) * Q) >= S)
}

main_rep <- function(ce, r) {
  suppressPackageStartupMessages({library(glmnet); library(GRPtests); library(PLStests)})
  d <- make_data(ce, r); X <- d$X; y <- d$y; p <- ncol(X)
  X1 <- cbind(1, X); D <- diag(c(0, rep(1, p)))
  tm <- c(); tic <- function() proc.time()[[3]]
  t0 <- tic()
  cv <- cv.glmnet(X, y, family = "binomial", alpha = 0, nfolds = 10)
  lam <- cv$lambda.1se * N
  fit <- ridge_fit(X1, y, lam, D); st <- pieces(X1, y, fit, lam, D)
  # textbook HL: sum (o - e)^2 / (e (1 - e/n_g)) on the ridge deciles, chi-square G - 2
  hl <- sum((st$o - st$e)^2 / (st$e * (1 - st$e / st$n_g)))
  p_hl <- pchisq(hl, G - 2, lower.tail = FALSE)
  p_edge_u <- mc_edge_p(st$edge_u, st$Om, st$Z)
  gen <- pmin(pmax(as.numeric(1 / (1 + exp(-X1 %*% st$bt))), 1e-6), 1 - 1e-6)
  Sd <- Se <- numeric(NB)
  for (bb in seq_len(NB)) {
    ys <- rbinom(N, 1, gen); s2 <- pieces(X1, ys, ridge_fit(X1, ys, lam, D), lam, D)
    Sd[bb] <- s2$dec_c; Se[bb] <- s2$edge_c
  }
  p_schl <- (1 + sum(Sd >= st$dec_c)) / (NB + 1)
  p_scedge <- (1 + sum(Se >= st$edge_c, na.rm = TRUE)) / (NB + 1)
  tm["ours"] <- tic() - t0
  t0 <- tic(); g5 <- tryCatch(as.numeric(GRPtest(X, y, fam = "binomial", nsplits = 5L)), error = function(e) NA_real_)
  tm["grp5"] <- tic() - t0
  t0 <- tic(); g1 <- tryCatch(as.numeric(GRPtest(X, y, fam = "binomial", nsplits = 1L)), error = function(e) NA_real_)
  tm["grp1"] <- tic() - t0
  t0 <- tic(); pl <- tryCatch(suppressMessages(PLStests(y = matrix(y), x = X, family = "binomial")),
                              error = function(e) list(T_cauchy = NA, T_alpha = NA, T_hmp = NA, T_beta = NA))
  tm["pls"] <- tic() - t0
  gm <- suppressWarnings(tryCatch(glm.fit(X1, y, family = binomial()), error = function(e) NULL))
  data.frame(cell = ce$cell, rep = r, fp_y = d$fp_y, fp_x = d$fp_x, lambda = lam,
             p_hl = p_hl, p_edge_u = p_edge_u, p_schl = p_schl, p_scedge = p_scedge,
             p_grp5 = g5, p_grp1 = g1,
             p_pls = as.numeric(pl$T_cauchy), p_pls_alpha = as.numeric(pl$T_alpha),
             p_pls_hmp = as.numeric(pl$T_hmp), p_pls_beta = as.numeric(pl$T_beta),
             glm_conv = if (is.null(gm)) NA else gm$converged,
             glm_maxcoef = if (is.null(gm)) NA_real_ else max(abs(gm$coefficients), na.rm = TRUE),
             gen_min = min(gen), gen_max = max(gen),
             t_ours = tm[["ours"]], t_grp5 = tm[["grp5"]], t_grp1 = tm[["grp1"]], t_pls = tm[["pls"]])
}

bag_rep <- function(ce, r) {
  suppressPackageStartupMessages(library(BAGofT))
  d <- make_data(ce, r)
  df <- data.frame(y = d$y, d$X)
  t0 <- proc.time()[[3]]
  res <- tryCatch(suppressMessages(BAGofT(testModel = testGlmnet(formula = y ~ ., alpha = 0), data = df)),
                  error = function(e) NULL)
  data.frame(cell = ce$cell, rep = r, fp_y = d$fp_y, fp_x = d$fp_x,
             p_bag = if (is.null(res)) NA_real_ else res$p.value,
             p_bag_min = if (is.null(res)) NA_real_ else res$p.value3,
             stat_bag_mean = if (is.null(res)) NA_real_ else res$pmean,
             t_bag = proc.time()[[3]] - t0)
}

# ---- driver -------------------------------------------------------------------------------------------
append_rows <- function(df, path) {
  write.table(df, path, sep = ",", row.names = FALSE, col.names = !file.exists(path), append = file.exists(path))
}
done_pairs <- function(path) {
  if (!file.exists(path)) return(character(0))
  x <- read.csv(path); paste(x$cell, x$rep)
}
run_jobs <- function(jobs, fun, path, nw) {
  done <- done_pairs(path)
  jobs <- jobs[!paste(jobs$cell, jobs$rep) %in% done, , drop = FALSE]
  cat(sprintf("%s: %d jobs to run (%d already on disk)\n", basename(path), nrow(jobs), length(done)))
  if (!nrow(jobs)) return(invisible())
  cl <- makeCluster(nw); on.exit(stopCluster(cl))
  clusterExport(cl, c("cells", "make_data", "ridge_fit", "pieces", "mc_edge_p", "main_rep", "bag_rep",
                      "N", "G", "NB", "GAMMA"), envir = globalenv())
  T0 <- Sys.time(); chunk <- 4L * nw
  for (s in seq(1L, nrow(jobs), by = chunk)) {
    J <- jobs[s:min(s + chunk - 1L, nrow(jobs)), ]
    # the job function travels as 'job_fun': 'fun' is parLapplyLB's own argument name
    R <- parLapplyLB(cl, seq_len(nrow(J)), function(i, J, job_fun) {
      ce <- cells[cells$cell == J$cell[i], ]; job_fun(ce, J$rep[i]) }, J = J, job_fun = fun)
    append_rows(do.call(rbind, R), path)
    cat(sprintf("[%7.1f min] %s %d / %d\n", as.numeric(difftime(Sys.time(), T0, units = "mins")),
                basename(path), min(s + chunk - 1L, nrow(jobs)), nrow(jobs))); flush.console()
  }
}

if (MODE == "smoke") {
  # every kind of cell, one replicate, serially: catches a broken branch before hours are spent
  pick <- do.call(rbind, lapply(split(cells, paste(cells$panel, cells$dep, cells$signal)), function(z) z[1, ]))
  pick <- rbind(pick, cells[cells$p == 180 & cells$panel == "M" & cells$dep == "null", ])
  for (i in seq_len(nrow(pick))) {
    t0 <- Sys.time(); x <- main_rep(pick[i, ], 1L)
    cat(sprintf("cell %3d %s p=%3d rho=%.1f %-6s %-5s a=%.1f | HL %.3f EDGEu %.3f SC.HL %.3f SC.EDGE %.3f GRP5 %.3f GRP1 %.3f PLS %.3f | %.1fs\n",
        x$cell, pick$panel[i], pick$p[i], pick$rho[i], pick$signal[i], pick$dep[i], pick$a[i], x$p_hl, x$p_edge_u,
        x$p_schl, x$p_scedge, x$p_grp5, x$p_grp1, x$p_pls, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  # data identity between the two job lists
  ce <- cells[cells$cell == 1, ]; a <- make_data(ce, 7L); b <- make_data(ce, 7L)
  stopifnot(identical(a$X, b$X), identical(a$y, b$y))
  # the parallel driver too (its own scratch file), and resume: a second call must find nothing to do
  tmp <- tempfile(fileext = ".csv")
  run_jobs(data.frame(cell = c(1L, 2L), rep = c(1L, 1L)), main_rep, tmp, 2L)
  run_jobs(data.frame(cell = c(1L, 2L), rep = c(1L, 1L)), main_rep, tmp, 2L)
  chk <- read.csv(tmp); stopifnot(nrow(chk) == 2, chk$fp_y[chk$cell == 1] == make_data(cells[1, ], 1L)$fp_y)
  cat("smoke OK:", nrow(cells), "cells,", sum(cells$reps), "main replicates,", sum(cells$bag_reps), "BAGofT replicates\n")
} else if (MODE == "pilot") {
  # timing on the most expensive cells, written to their own file so the main run does not inherit them
  OUT_MAIN <- file.path("..", "data", "T40_pilot_main.csv"); OUT_BAG <- file.path("..", "data", "T40_pilot_bag.csv")
  big <- cells[cells$panel == "M" & cells$dep == "null", ]
  run_jobs(expand.grid(cell = big$cell, rep = 1:4), main_rep, OUT_MAIN, NW)
  run_jobs(expand.grid(cell = big$cell[big$p %in% c(20, 100, 180)], rep = 1:2), bag_rep, OUT_BAG, NW)
} else if (MODE == "main") {
  only <- Sys.getenv("T40_ONLY", ""); cs <- if (nzchar(only)) as.integer(strsplit(only, ",")[[1]]) else cells$cell
  jobs <- do.call(rbind, lapply(cs, function(c) data.frame(cell = c, rep = seq_len(cells$reps[cells$cell == c]))))
  run_jobs(jobs, main_rep, OUT_MAIN, NW)
} else if (MODE == "recheck") {
  # Amendment 5: an independent second block for every null cell at p/n >= 0.25. Replicates 1001..2000 use
  # seeds cell * 1e6 + 1001.. , which the main run never touched, so the two blocks are independent.
  out <- Sys.getenv("T40_OUT", file.path("..", "data", "T40_recheck_pvalues.csv"))
  nrep <- as.integer(Sys.getenv("T40_REPS", "1000"))
  cs <- as.integer(strsplit(Sys.getenv("T40_ONLY", "15,22,29,44,48,52,64,68,72,84,88,92"), ",")[[1]])
  stopifnot(all(cells$dep[match(cs, cells$cell)] == "null"))
  jobs <- do.call(rbind, lapply(cs, function(c) data.frame(cell = c, rep = 1000L + seq_len(nrep))))
  run_jobs(jobs, main_rep, out, NW)
} else if (MODE == "bag") {
  only <- Sys.getenv("T40_ONLY", ""); cs <- if (nzchar(only)) as.integer(strsplit(only, ",")[[1]]) else cells$cell
  cs <- cs[cells$bag_reps[match(cs, cells$cell)] > 0]
  jobs <- do.call(rbind, lapply(cs, function(c) data.frame(cell = c, rep = seq_len(cells$bag_reps[cells$cell == c]))))
  run_jobs(jobs, bag_rep, OUT_BAG, NW)
}
