# T44 -- PILOT (not confirmatory) of a rescaled bootstrap generator for the corrected test at large p/n.
#
# Why: T43 showed that the debiased estimate beta_tilde over-states the signal as p/n grows (sd of its linear
# predictor 0.96x the truth at p/n = 0.05, 1.10-1.21x at 0.45, worst for a sparse truth) -- the Sur-Candes
# inflation of the maximum likelihood estimate, which beta_tilde equals to first order. The bootstrap then draws
# its null world from a model that is too sharp, and the corrected test becomes liberal in exactly those cells.
#
# Fix tested here: estimate the inflation at the generator itself and undo it (one step of bootstrap bias
# correction). Draw K0 = 50 responses from pi(beta_tilde), refit at the same lambda, debias, and take
# kappa = median sd(X beta_tilde*) / sd(X beta_tilde) over the slopes. Shrink the slopes by 1/kappa and re-solve
# the intercept so the generator's expected events equal the observed events. Everything else is the T40 test.
#
# Paired with T40/T42: same cells, same seeds, replicates 1..R, so each new p-value sits beside the stored one.
# Pre-stated reading (pilot): the fix is worth a confirmatory run if every failing cell falls to <= 0.08 and the
# two cells that were valid stay >= 0.03 (100 replicates per cell: s.e. about 0.02, so this is a screen only).
#   Rscript T44_rescaled_generator_pilot.R        env T44_NW (default 2), T44_R (default 100)
suppressPackageStartupMessages({library(glmnet); library(parallel)})
src <- readLines("T40_highdim_grid.R")
eval(parse(text = src[seq_len(grep("^# ---- driver", src) - 1)]))   # cells, ridge_fit, pieces, N, G, NB, GAMMA
NW <- as.integer(Sys.getenv("T44_NW", "2")); RR <- as.integer(Sys.getenv("T44_R", "100")); K0 <- 50L
OUT <- file.path("..", "data", "T44_rescaled_pilot.csv")

# cells: T40 ids keep their own spec; T42 ids (205, 209) are panel M with intercept ALPHA
ALPHA42 <- -1.530973
spec <- rbind(cells[cells$cell %in% c(15L, 72L, 88L, 92L), c("cell", "p", "rho", "signal", "dep", "a")],
              data.frame(cell = c(205L, 209L), p = c(100L, 180L), rho = 0.7, signal = "dense", dep = "null", a = 0))
spec$alpha <- ifelse(spec$cell > 200, ALPHA42, 0)

make_data_a <- function(ce, r) {                 # T40/T42 generator, identical seeds and draws
  set.seed(ce$cell * 1e6 + r)
  p <- ce$p; S <- ce$rho^abs(outer(seq_len(p), seq_len(p), "-"))
  b <- if (ce$signal == "dense") rep(c(0.35, -0.30, 0.25, -0.20, 0.15), length.out = p)
       else c(rep(1, 5), rep(0, p - 5))
  b <- b * GAMMA / sqrt(drop(t(b) %*% S %*% b))
  X <- matrix(rnorm(N * p), N, p) %*% chol(S)
  e0 <- drop(X %*% b)
  y <- rbinom(N, 1, plogis(ce$alpha + e0))
  list(X = X, y = y, fp_y = sum(y))
}

rescaled_rep <- function(ce, r) {
  d <- make_data_a(ce, r); X <- d$X; y <- d$y; p <- ncol(X)
  X1 <- cbind(1, X); D <- diag(c(0, rep(1, p)))
  cv <- cv.glmnet(X, y, family = "binomial", alpha = 0, nfolds = 10)
  lam <- cv$lambda.1se * N
  fit <- ridge_fit(X1, y, lam, D); st <- pieces(X1, y, fit, lam, D)
  t0 <- proc.time()[[3]]
  # stage 1: inflation of the debiased slopes, measured at the unadjusted generator
  gen0 <- pmin(pmax(plogis(drop(X1 %*% st$bt)), 1e-6), 1 - 1e-6)
  sd0 <- sd(drop(X %*% st$bt[-1]))
  ratios <- vapply(seq_len(K0), function(k) {
    ys <- rbinom(N, 1, gen0); s2 <- pieces(X1, ys, ridge_fit(X1, ys, lam, D), lam, D)
    sd(drop(X %*% s2$bt[-1])) / sd0 }, 0.0)
  kappa <- median(ratios)
  # stage 2: shrink the slopes, then match the observed number of events with the intercept
  eta_s <- drop(X %*% st$bt[-1]) / kappa
  a0 <- uniroot(function(a) sum(plogis(a + eta_s)) - sum(y), c(-20, 20))$root
  gen <- pmin(pmax(plogis(a0 + eta_s), 1e-6), 1 - 1e-6)
  Sd <- Se <- numeric(NB)
  for (bb in seq_len(NB)) {
    ys <- rbinom(N, 1, gen); s2 <- pieces(X1, ys, ridge_fit(X1, ys, lam, D), lam, D)
    Sd[bb] <- s2$dec_c; Se[bb] <- s2$edge_c
  }
  data.frame(cell = ce$cell, rep = r, fp_y = d$fp_y, kappa = kappa,
             p_schl_rs = (1 + sum(Sd >= st$dec_c)) / (NB + 1),
             p_scedge_rs = (1 + sum(Se >= st$edge_c, na.rm = TRUE)) / (NB + 1),
             t_rs = proc.time()[[3]] - t0)
}

jobs <- expand.grid(rep = seq_len(RR), cell = spec$cell)[, c("cell", "rep")]
# interleave cells so an early stop still gives every cell some replicates
jobs <- jobs[order(jobs$rep, jobs$cell), ]
if (file.exists(OUT)) { x <- read.csv(OUT); jobs <- jobs[!paste(jobs$cell, jobs$rep) %in% paste(x$cell, x$rep), ] }
cat(nrow(jobs), "jobs to run\n")
if (nrow(jobs)) {
  cl <- makeCluster(NW); on.exit(stopCluster(cl))
  clusterExport(cl, c("spec", "make_data_a", "rescaled_rep", "ridge_fit", "pieces", "N", "G", "NB", "GAMMA", "K0"))
  T0 <- Sys.time()
  for (s in seq(1L, nrow(jobs), by = 2L * NW)) {
    J <- jobs[s:min(s + 2L * NW - 1L, nrow(jobs)), ]
    R <- parLapplyLB(cl, seq_len(nrow(J)), function(i, J) {
      suppressPackageStartupMessages(library(glmnet))
      rescaled_rep(spec[spec$cell == J$cell[i], ], J$rep[i]) }, J = J)
    R <- do.call(rbind, R)
    write.table(R, OUT, sep = ",", row.names = FALSE, col.names = !file.exists(OUT), append = file.exists(OUT))
    cat(sprintf("[%6.1f min] %d / %d\n", as.numeric(difftime(Sys.time(), T0, units = "mins")),
                min(s + 2L * NW - 1L, nrow(jobs)), nrow(jobs))); flush.console()
  }
}
