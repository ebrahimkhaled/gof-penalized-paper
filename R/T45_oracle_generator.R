# T45 -- ORACLE DIAGNOSTIC (not a usable test): what about the bootstrap generator makes the corrected test
# liberal at large p/n? T44 showed the bootstrap cannot see the excess spread of the debiased fit (its own
# estimate of the signal strength is as inflated as the naive one), so before building an estimator we ask
# which property of the generator matters, using the true signal strength sd(eta0) = 1.5 that a user never has.
#   V_dir_tilde: debiased direction, slopes rescaled so sd(X b) = 1.5   -> tests "strength is the problem"
#   V_dir_ridge: ridge direction,    slopes rescaled so sd(X b) = 1.5   -> tests "noisy direction is the problem"
# Both re-solve the intercept to match the observed events. Paired with T40/T42 on replicates 1..R.
# NB = 199 (size resolution is ample for a diagnostic). Rscript T45_oracle_generator.R  env T45_NW, T45_R
suppressPackageStartupMessages({library(glmnet); library(parallel)})
src <- readLines("T44_rescaled_generator_pilot.R")
eval(parse(text = src[seq_len(grep("^rescaled_rep <- function", src) - 1)]))
NW <- as.integer(Sys.getenv("T45_NW", "2")); RR <- as.integer(Sys.getenv("T45_R", "60")); NB45 <- 199L
OUT <- file.path("..", "data", "T45_oracle_generator.csv")

boot_p <- function(gen, X1, lam, D, obs) {
  Sd <- numeric(NB45)
  for (bb in seq_len(NB45)) {
    ys <- rbinom(N, 1, gen); Sd[bb] <- pieces(X1, ys, ridge_fit(X1, ys, lam, D), lam, D)$dec_c }
  (1 + sum(Sd >= obs)) / (NB45 + 1)
}
oracle_rep <- function(ce, r) {
  d <- make_data_a(ce, r); X <- d$X; y <- d$y; p <- ncol(X)
  X1 <- cbind(1, X); D <- diag(c(0, rep(1, p)))
  lam <- cv.glmnet(X, y, family = "binomial", alpha = 0, nfolds = 10)$lambda.1se * N
  fit <- ridge_fit(X1, y, lam, D); st <- pieces(X1, y, fit, lam, D)
  gen_at <- function(slopes) {
    eta <- drop(X %*% slopes); eta <- eta * GAMMA / sd(eta)
    a0 <- uniroot(function(a) sum(plogis(a + eta)) - sum(y), c(-20, 20))$root
    pmin(pmax(plogis(a0 + eta), 1e-6), 1 - 1e-6)
  }
  data.frame(cell = ce$cell, rep = r, fp_y = d$fp_y,
             p_tilde_oracle = boot_p(gen_at(st$bt[-1]), X1, lam, D, st$dec_c),
             p_ridge_oracle = boot_p(gen_at(fit$beta[-1]), X1, lam, D, st$dec_c))
}

jobs <- expand.grid(rep = seq_len(RR), cell = c(15L, 92L, 209L))[, c("cell", "rep")]
jobs <- jobs[order(jobs$rep, jobs$cell), ]
if (file.exists(OUT)) { x <- read.csv(OUT); jobs <- jobs[!paste(jobs$cell, jobs$rep) %in% paste(x$cell, x$rep), ] }
cat(nrow(jobs), "jobs to run\n")
if (nrow(jobs)) {
  cl <- makeCluster(NW); on.exit(stopCluster(cl))
  clusterExport(cl, c("spec", "make_data_a", "oracle_rep", "boot_p", "ridge_fit", "pieces", "N", "G", "GAMMA", "NB45"))
  T0 <- Sys.time()
  for (s in seq(1L, nrow(jobs), by = 2L * NW)) {
    J <- jobs[s:min(s + 2L * NW - 1L, nrow(jobs)), ]
    R <- do.call(rbind, parLapplyLB(cl, seq_len(nrow(J)), function(i, J) {
      suppressPackageStartupMessages(library(glmnet))
      oracle_rep(spec[spec$cell == J$cell[i], ], J$rep[i]) }, J = J))
    write.table(R, OUT, sep = ",", row.names = FALSE, col.names = !file.exists(OUT), append = file.exists(OUT))
    cat(sprintf("[%6.1f min] %d / %d\n", as.numeric(difftime(Sys.time(), T0, units = "mins")),
                min(s + 2L * NW - 1L, nrow(jobs)), nrow(jobs))); flush.console()
  }
}
