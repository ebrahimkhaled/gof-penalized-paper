# E13 DIAGNOSTIC — why is EDGE conservative in design B at the LIGHT penalty (lambda=50)?
#
# Manifest 8f(b) asks for "the eigenvalue spectrum of Omega_K in both designs at
# lambda in {50,416,1000}". That is computed here (part 1). But a PREPIVOTED test cannot
# be conservative merely because a reference distribution has the wrong shape -- the
# reference is generated, not assumed. So part 2 tests the mechanism that CAN do it:
# the bootstrap generator pi(beta.tilde) with beta.tilde = beta.hat + F^{-1} K beta.hat.
# When lambda is light and p/n large, F is nearly singular, so F^{-1}K beta.hat can
# OVERSHOOT. A generator further from the origin than the truth makes the bootstrap
# world more variable than the real world, so S* > S on average and p-values inflate.
#
# Part 2 is the decisive one: it measures E[S*] - E[S] directly.
suppressPackageStartupMessages({library(parallel)})
G <- 10; NB <- 99; REPS <- 60

design <- function(panel) {
  if (panel == "A") {
    n <- 500; p <- 5; b0 <- c(0.5,-0.4,0.3,-0.3,0.2)
    X <- matrix(rnorm(n*p), n, p)
  } else {
    n <- 400; p <- 100; S <- 0.7^abs(outer(1:p, 1:p, "-")); Ch <- chol(S)
    b0 <- rep(c(.35,-.30,.25,-.20,.15), length.out = p) * 0.8887
    X <- matrix(rnorm(n*p), n, p) %*% Ch
  }
  list(X1 = cbind(1, X), y = rbinom(n, 1, 1/(1+exp(-drop(X %*% b0)))),
       D = diag(c(0, rep(1, p))), b0 = c(0, b0), n = n, p = p)
}

ridge_fit <- function(X1, y, lam, D, tol = 1e-9, maxit = 60) {
  beta <- rep(0, ncol(X1))
  for (it in 1:maxit) {
    eta <- drop(X1 %*% beta); pi <- 1/(1+exp(-eta))
    w <- pmax(pi*(1-pi), 1e-8); z <- eta + (y-pi)/w
    bn <- drop(solve(crossprod(X1, w*X1) + lam*D, crossprod(X1, w*z)))
    if (max(abs(bn-beta)) < tol) { beta <- bn; break }; beta <- bn
  }
  list(beta = beta, pi = 1/(1+exp(-drop(X1 %*% beta))))
}

# returns the statistics AND the geometry
pieces <- function(X1, y, fit, lam, D, geom = FALSE) {
  pi <- pmin(pmax(fit$pi, 1e-6), 1-1e-6); w <- pmax(pi*(1-pi), 1e-8)
  g <- pmin(ceiling(rank(pi, ties.method="first")/(length(y)/G)), G)
  idx <- split(seq_along(y), g); Vg <- vapply(idx, function(I) sum(w[I]), 0.0)
  r <- (vapply(idx, function(I) sum(y[I]), 0.0) -
        vapply(idx, function(I) sum(pi[I]), 0.0))/sqrt(Vg)
  U <- t(vapply(idx, function(I) colSums(w[I]*X1[I,,drop=FALSE]),
                numeric(ncol(X1))))/sqrt(Vg)
  Fm <- crossprod(X1, w*X1); K <- lam*D; M <- Fm + K
  bt <- fit$beta + drop(solve(Fm, K %*% fit$beta))
  mu <- drop(U %*% solve(M, K %*% bt))
  pbar <- vapply(idx, function(I) mean(pi[I]), 0.0); v <- r - mu
  Z <- tryCatch(as.matrix(stats::poly(pbar, 3)), error = function(e) NULL)
  ed <- if (is.null(Z)) NA_real_ else {
    zr <- crossprod(Z, v); as.numeric(t(zr) %*% solve(crossprod(Z)) %*% zr) }
  out <- list(dec = sum(v^2), edge = ed, bt = bt)
  if (geom) {
    Om_K   <- diag(G) - U %*% solve(M, (Fm + 2*K) %*% solve(M, t(U)))  # Prop 1
    Om_MLE <- diag(G) - U %*% solve(Fm, t(U))                          # Prop 2 target
    out$ev_K   <- sort(eigen(Om_K,   symmetric = TRUE, only.values = TRUE)$values)
    out$ev_MLE <- sort(eigen(Om_MLE, symmetric = TRUE, only.values = TRUE)$values)
    out$df_K   <- sum(diag(Om_K)); out$df_MLE <- sum(diag(Om_MLE))
    if (!is.null(Z)) {                # effective df the 3 EDGE directions actually get
      Zn <- Z %*% solve(chol(crossprod(Z)))
      out$df_edge <- sum(diag(crossprod(Zn, Om_MLE %*% Zn)))
    }
  }
  out
}

worker <- function(seed, panel, lam) {
  set.seed(seed)
  d  <- design(panel)
  ft <- ridge_fit(d$X1, d$y, lam, d$D)
  st <- pieces(d$X1, d$y, ft, lam, d$D, geom = TRUE)
  gen <- pmin(pmax(as.numeric(1/(1+exp(-d$X1 %*% st$bt))), 1e-6), 1-1e-6)
  Sd <- Se <- numeric(NB)
  for (b in 1:NB) {
    ys <- rbinom(d$n, 1, gen); s2 <- pieces(d$X1, ys, ridge_fit(d$X1, ys, lam, d$D), lam, d$D)
    Sd[b] <- s2$dec; Se[b] <- s2$edge
  }
  c(dec = st$dec, edge = st$edge,
    boot_dec = mean(Sd), boot_edge = mean(Se, na.rm = TRUE),
    df_K = st$df_K, df_MLE = st$df_MLE, df_edge = st$df_edge,
    evmin_MLE = st$ev_MLE[1], evmax_MLE = st$ev_MLE[G],
    nb0 = sqrt(sum(d$b0^2)), nbh = sqrt(sum(ft$beta^2)), nbt = sqrt(sum(st$bt^2)),
    p_dec = (1+sum(Sd >= st$dec))/(NB+1),
    p_edge = (1+sum(Se >= st$edge, na.rm = TRUE))/(NB+1))
}

nc <- max(1, detectCores() - 2)
cl <- makeCluster(nc); on.exit(stopCluster(cl))
clusterExport(cl, c("G","NB","design","ridge_fit","pieces"))
clusterSetRNGStream(cl, 20260808)
cat("E13 diagnostic on", nc, "workers; REPS =", REPS, "NB =", NB, "\n\n")

cells <- list(list(p="A",lam=100), list(p="A",lam=137), list(p="A",lam=200),
              list(p="B",lam=50),  list(p="B",lam=416), list(p="B",lam=1000))
res <- list()
for (ce in cells) {
  t0 <- Sys.time()
  M <- simplify2array(parLapply(cl, seq_len(REPS), worker, panel = ce$p, lam = ce$lam))
  m <- rowMeans(M, na.rm = TRUE)
  # rejection RATE, not the mean p-value
  rej_d <- mean(M["p_dec", ]  <= 0.05, na.rm = TRUE)
  rej_e <- mean(M["p_edge", ] <= 0.05, na.rm = TRUE)
  res[[paste0(ce$p, ce$lam)]] <- c(panel = ce$p, lam = ce$lam, m,
                                   rej_dec = rej_d, rej_edge = rej_e)
  cat(sprintf(
    "%s lam=%-5g | df: Om_K %.2f  Om_MLE %.2f  EDGE-3 %.3f | ev_MLE [%.3f, %.3f]\n",
    ce$p, ce$lam, m["df_K"], m["df_MLE"], m["df_edge"], m["evmin_MLE"], m["evmax_MLE"]))
  cat(sprintf("            | E[S]  dec %.2f  edge %.3f | E[S*] dec %.2f  edge %.3f",
              m["dec"], m["edge"], m["boot_dec"], m["boot_edge"]))
  cat(sprintf("  <- edge ratio %.3f\n", m["boot_edge"]/m["edge"]))
  cat(sprintf("            | ||b0|| %.2f  ||bhat|| %.2f  ||btilde|| %.2f  overshoot %.2fx",
              m["nb0"], m["nbh"], m["nbt"], m["nbt"]/m["nb0"]))
  cat(sprintf(" | size dec %.3f edge %.3f  (%.0fs)\n\n", rej_d, rej_e,
              as.numeric(difftime(Sys.time(), t0, units="secs"))))
}
saveRDS(res, "omega_spectrum_results.rds")
cat("saved omega_spectrum_results.rds\n")
