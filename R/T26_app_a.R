# T2.6(a) -- RE-MEASURE THE THREE AMBER CONSTANTS OF MANIFEST 7.3
#
# Claim under test: "Remove mu.hat, keep prepivoting and refits identical:
#   size 0.020, power 0.000 against every departure, because ||r||^2 is then
#   dominated by ||mu_K||^2 ~= 57.5.  The corrected statistic removes 87% of the
#   shrinkage displacement, cutting ||E(r|H0)|| from 7.867 -> 1.005."
# All three constants are AMBER (single agent run). Re-measure them.
#
# TWO EXPERIMENTS, both at DGP B, lambda = 416.
#   (A1) cheap, fixed design: ||mu_K||^2, ||E(r|H0)||, ||E(r-mu.hat|H0)||, % removed
#   (A2) expensive, random design: prepivoted test with and WITHOUT mu.hat
#        subtracted, PAIRED on identical fits and identical bootstrap refits,
#        under H0 and under the quadratic and cubic alternatives.
#
# Functions ridge_fit / pieces / prepivot are reused VERBATIM from
# edge_basis_study.R and oracle_knownnull.R (only the return list is widened to
# carry the uncorrected statistics as well, which `pieces` already computed).
suppressPackageStartupMessages(library(parallel))

G <- 10; NB <- 149; LAM <- 416
B_CELL <- as.integer(Sys.getenv("B_CELL", "600"))
NW     <- as.integer(Sys.getenv("NW", "10"))

## ---- DGP B (frozen) --------------------------------------------------------
nB <- 400; pB <- 100
SigB <- 0.7^abs(outer(1:pB, 1:pB, "-")); ChB <- chol(SigB)
bB <- rep(c(.35, -.30, .25, -.20, .15), length.out = pB) * 0.8887
DB <- diag(c(0, rep(1, pB)))
vB <- as.numeric(t(bB) %*% SigB %*% bB)          # var(eta0)

## ---- verbatim workhorses ---------------------------------------------------
ridge_fit <- function(X1, y, lam, D, tol = 1e-9, maxit = 60) {
  beta <- rep(0, ncol(X1))
  for (it in 1:maxit) { eta <- drop(X1 %*% beta); pi <- 1/(1+exp(-eta))
    w <- pmax(pi*(1-pi), 1e-8); z <- eta + (y-pi)/w
    bn <- drop(solve(crossprod(X1, w*X1) + lam*D, crossprod(X1, w*z)))
    if (max(abs(bn-beta)) < tol) { beta <- bn; break }; beta <- bn }
  list(beta = beta, pi = 1/(1+exp(-drop(X1 %*% beta))))
}

# identical to edge_basis_study.R::pieces, but also returns r, mu, and lets the
# grouping vector be supplied (used by A1 for the fixed-grouping variant)
pieces <- function(X1, y, fit, lam, D, g = NULL) {
  pi <- pmin(pmax(fit$pi, 1e-6), 1-1e-6); w <- pmax(pi*(1-pi), 1e-8)
  if (is.null(g)) g <- pmin(ceiling(rank(pi, ties.method="first")/(length(y)/G)), G)
  idx  <- split(seq_along(y), g)
  Vg   <- vapply(idx, function(I) sum(w[I]), 0.0)
  og   <- vapply(idx, function(I) sum(y[I]), 0.0)
  eg   <- vapply(idx, function(I) sum(pi[I]), 0.0)
  pbar <- vapply(idx, function(I) mean(pi[I]), 0.0)
  r    <- (og - eg)/sqrt(Vg)
  U    <- t(vapply(idx, function(I) colSums(w[I]*X1[I,,drop=FALSE]), numeric(ncol(X1))))/sqrt(Vg)
  Fm   <- crossprod(X1, w*X1); K <- lam*D
  bt   <- fit$beta + drop(solve(Fm, K %*% fit$beta))
  mu   <- drop(U %*% solve(Fm + K, K %*% bt))
  rc   <- r - mu
  Z  <- tryCatch(as.matrix(stats::poly(pbar, 3)), error = function(e) NULL)
  qf <- function(v, Z) { Zr <- crossprod(Z, v); as.numeric(t(Zr) %*% solve(crossprod(Z)) %*% Zr) }
  list(dec_c = sum(rc^2), edge_c = if (is.null(Z)) NA_real_ else qf(rc, Z),
       dec_p = sum(r^2),  edge_p = if (is.null(Z)) NA_real_ else qf(r, Z),
       r = r, mu = mu, bt = bt, Z = Z)
}

# prepivoted p-values for FOUR statistics from ONE set of refits:
#   corrected decile / corrected EDGE / UNcorrected decile / UNcorrected EDGE
p_prepivot4 <- function(X1, y, lam, D) {
  ft <- ridge_fit(X1, y, lam, D); st <- pieces(X1, y, ft, lam, D)
  gen <- pmin(pmax(as.numeric(1/(1+exp(-X1 %*% st$bt))), 1e-6), 1-1e-6)
  S <- matrix(NA_real_, NB, 4)
  for (bb in 1:NB) {
    ys <- rbinom(nrow(X1), 1, gen)
    s2 <- pieces(X1, ys, ridge_fit(X1, ys, lam, D), lam, D)
    S[bb, ] <- c(s2$dec_c, s2$edge_c, s2$dec_p, s2$edge_p)
  }
  obs <- c(st$dec_c, st$edge_c, st$dec_p, st$edge_p)
  p <- (1 + colSums(sweep(S, 2, obs, ">="), na.rm = TRUE))/(NB + 1)
  names(p) <- c("corr_dec", "corr_edge", "unc_dec", "unc_edge")
  p
}

## =================== A1: the three constants (fixed design) =================
cat("=== A1  fixed-design Monte Carlo of E(r | H0), DGP B, lambda = 416 ===\n")
set.seed(2026)
Xf  <- matrix(rnorm(nB*pB), nB, pB) %*% ChB
X1f <- cbind(1, Xf)
eta0 <- drop(Xf %*% bB); pi0 <- 1/(1+exp(-eta0))
gfix <- pmin(ceiling(rank(pi0, ties.method="first")/(nB/G)), G)   # deciles of TRUE pi

# theoretical mu_K at the truth: mu_K = U (F+K)^{-1} K beta0, everything at pi0,
# grouped by the true-probability deciles.
w0  <- pi0*(1-pi0); idx0 <- split(seq_len(nB), gfix)
Vg0 <- vapply(idx0, function(I) sum(w0[I]), 0.0)
U0  <- t(vapply(idx0, function(I) colSums(w0[I]*X1f[I,,drop=FALSE]), numeric(pB+1)))/sqrt(Vg0)
F0  <- crossprod(X1f, w0*X1f); K0 <- LAM*DB
b0f <- c(0, bB)
muK <- drop(U0 %*% solve(F0 + K0, K0 %*% b0f))
cat(sprintf("  ||mu_K||^2 = %.3f   (||mu_K|| = %.3f)\n", sum(muK^2), sqrt(sum(muK^2))))

B0 <- 2000
Rf <- Rc <- matrix(NA_real_, B0, G)   # fixed grouping (deciles of true pi)
Rf2 <- Rc2 <- matrix(NA_real_, B0, G) # data-driven grouping (deciles of pi.hat)
for (b in 1:B0) {
  ys <- rbinom(nB, 1, pi0)
  ft <- ridge_fit(X1f, ys, LAM, DB)
  s1 <- pieces(X1f, ys, ft, LAM, DB, g = gfix); Rf[b,] <- s1$r; Rc[b,] <- s1$r - s1$mu
  s2 <- pieces(X1f, ys, ft, LAM, DB);           Rf2[b,] <- s2$r; Rc2[b,] <- s2$r - s2$mu
}
nrm <- function(M) sqrt(sum(colMeans(M)^2))
se_nrm <- function(M) sqrt(sum(apply(M,2,var))/nrow(M))   # rough SE of the norm
a1 <- data.frame(
  grouping = c("true-pi deciles (fixed)", "pi.hat deciles (data-driven)"),
  norm_Er_before = c(nrm(Rf), nrm(Rf2)),
  norm_Er_after  = c(nrm(Rc), nrm(Rc2)),
  pct_removed    = c(100*(1 - nrm(Rc)/nrm(Rf)), 100*(1 - nrm(Rc2)/nrm(Rf2))),
  mc_noise_floor = c(se_nrm(Rf), se_nrm(Rf2)))
print(a1, row.names = FALSE)
cat(sprintf("\n  ||mu_K||^2 = %.3f ; B0 = %d\n\n", sum(muK^2), B0))

## =================== A2: the test with and without mu.hat ==================
cat("=== A2  prepivoted test, mu.hat subtracted vs NOT (paired), B =", B_CELL, "===\n")
cells <- list(list(k="null",  g=0.0), list(k="quad", g=1.0),
              list(k="cubic", g=1.0), list(k="cubic", g=1.5))

mkB <- function(g, k) {
  X <- matrix(rnorm(nB*pB), nB, pB) %*% ChB; e0 <- drop(X %*% bB)
  eta <- if (k=="cubic") e0 + g*(e0^3 - 3*vB*e0)/sqrt(6*vB^3)
         else if (k=="quad") e0 + g*(X[,1]^2 - 1) else e0
  list(X1 = cbind(1, X), y = rbinom(nB, 1, 1/(1+exp(-eta))))
}

cl <- makeCluster(NW)
clusterExport(cl, c("ridge_fit","pieces","p_prepivot4","mkB","G","NB","LAM",
                    "nB","pB","ChB","bB","DB","vB"))
clusterSetRNGStream(cl, 2026)
t0 <- Sys.time(); out <- NULL
for (ce in cells) {
  tt <- Sys.time()
  P <- parLapply(cl, seq_len(B_CELL), function(i, k, g) {
         d <- mkB(g, k); p_prepivot4(d$X1, d$y, LAM, DB) }, k = ce$k, g = ce$g)
  P <- do.call(rbind, P)
  rej <- colMeans(P < 0.05)
  se  <- sqrt(rej*(1-rej)/B_CELL)
  out <- rbind(out, data.frame(alt = ce$k, gamma = ce$g, B = B_CELL,
    corr_dec = rej[1], corr_edge = rej[2], unc_dec = rej[3], unc_edge = rej[4],
    se_corr_dec = se[1], se_corr_edge = se[2], se_unc_dec = se[3], se_unc_edge = se[4],
    mins = as.numeric(difftime(Sys.time(), tt, units="mins"))))
  cat(sprintf("  %-5s g=%.1f | CORRECTED dec=%.3f edge=%.3f | UNCORRECTED dec=%.3f edge=%.3f  (%.1f min)\n",
      ce$k, ce$g, rej[1], rej[2], rej[3], rej[4],
      as.numeric(difftime(Sys.time(), tt, units="mins"))))
  saveRDS(list(a1=a1, muK2=sum(muK^2), a2=out), "T26_app_a_results.rds")
}
stopCluster(cl)
cat(sprintf("\ntotal A2 wall clock %.1f min (NW=%d workers)\n",
    as.numeric(difftime(Sys.time(), t0, units="mins")), NW))
print(out, row.names = FALSE)
saveRDS(list(a1=a1, muK2=sum(muK^2), a2=out), "T26_app_a_results.rds")
