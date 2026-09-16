# =====================================================================
# shrink.gof() -- the shrinkage-corrected Hosmer-Lemeshow test
#
#   SC.HL    the correction on the decile (Hosmer-Lemeshow) basis
#   SC.EDGE  the correction on the EDGE basis
#
# This is the reference implementation for
#
#   Ebrahim, E.K. and Arashi, M. (2026). Shrinkage invalidates the
#   Hosmer-Lemeshow test: goodness of fit for ridge logistic regression,
#   with an application to glaucoma diagnosis.
#
# It is SELF-CONTAINED: it depends only on base R and stats, so reproducing
# the paper does not depend on any package version. The same test is
# exported as shrink.gof() by the 'ebrahim.gof' package on CRAN from version
# 2.6.0; this copy adds the eigenvalue floor and the diagnostics that the
# revised paper reports, which the package does not yet have.
#
# The procedure (Sections 3.1 and 3.2 of the paper):
#   1. fit the ridge model by penalized IRLS;
#   2. form the grouped standardized residuals r;
#   3. subtract the estimated shrinkage non-centrality mu-hat (Proposition 3.1);
#   4. refer the result to a bootstrap built from the debiased generator
#      pi(beta-tilde) -- Beran prepivoting.
#
# lambda is on the THEORY scale: lambda = n * lambda_glmnet.
# =====================================================================

# ---- penalized IRLS ---------------------------------------------------
.sg_ridge <- function(X1, y, lambda, D, tol = 1e-10, maxit = 100) {
  beta <- rep(0, ncol(X1))
  for (it in seq_len(maxit)) {
    eta <- drop(X1 %*% beta)
    pi  <- 1 / (1 + exp(-eta))
    w   <- pmax(pi * (1 - pi), 1e-8)
    z   <- eta + (y - pi) / w
    bn  <- drop(solve(crossprod(X1, w * X1) + lambda * D,
                      crossprod(X1, w * z)))
    if (max(abs(bn - beta)) < tol) { beta <- bn; break }
    beta <- bn
  }
  list(beta = beta, pi = 1 / (1 + exp(-drop(X1 %*% beta))))
}

# ---- F^{-1} b, with the eigenvalues of F optionally floored -----------
# F is inverted in exactly one place, the debiasing step below, and on a
# strongly collinear design that inversion is the fragile part of the whole
# procedure. tau = 0 is a plain solve(); a positive tau raises every
# eigenvalue below tau * (largest one) to that value before inverting.
.sg_finv <- function(Fm, b, tau = 0) {
  if (tau <= 0) return(drop(solve(Fm, b)))
  e <- eigen(Fm, symmetric = TRUE)
  d <- pmax(e$values, tau * e$values[1])
  drop(e$vectors %*% (crossprod(e$vectors, b) / d))
}

# ---- grouped residuals, mu-hat, and both statistics -------------------
.sg_pieces <- function(X1, y, fit, lambda, D, G, finv.floor = 0) {
  pi <- pmin(pmax(fit$pi, 1e-6), 1 - 1e-6)
  w  <- pmax(pi * (1 - pi), 1e-8)
  g  <- pmin(ceiling(rank(pi, ties.method = "first") / (length(y) / G)), G)
  idx <- split(seq_along(y), g)
  Vg <- vapply(idx, function(I) sum(w[I]), 0.0)
  r  <- (vapply(idx, function(I) sum(y[I]),  0.0) -
         vapply(idx, function(I) sum(pi[I]), 0.0)) / sqrt(Vg)
  U  <- t(vapply(idx, function(I) colSums(w[I] * X1[I, , drop = FALSE]),
                 numeric(ncol(X1)))) / sqrt(Vg)
  Fm <- crossprod(X1, w * X1)
  K  <- lambda * D
  # one-step debias, then the shrinkage non-centrality it implies
  bt <- fit$beta + .sg_finv(Fm, drop(K %*% fit$beta), finv.floor)
  mu <- drop(U %*% solve(Fm + K, K %*% bt))
  v  <- r - mu                                   # the corrected residual
  pbar <- vapply(idx, function(I) mean(pi[I]), 0.0)
  Z <- tryCatch(as.matrix(stats::poly(pbar, 3)), error = function(e) NULL)
  qf <- function(u, Z) {
    zr <- crossprod(Z, u)
    as.numeric(t(zr) %*% solve(crossprod(Z)) %*% zr)
  }
  list(S_dec_unc = sum(r^2),
       S_dec     = sum(v^2),
       S_edge_unc = if (is.null(Z)) NA_real_ else qf(r, Z),
       S_edge     = if (is.null(Z)) NA_real_ else qf(v, Z),
       Om = diag(G) - U %*% solve(Fm, t(U)),     # the MLE covariance
       Z = Z, beta_tilde = bt, Fm = Fm)
}

# ---- the uncorrected reference, by direct Monte Carlo -----------------
.sg_mc_p <- function(S, Om, Z = NULL, nd = 20000) {
  if (is.na(S)) return(NA_real_)
  R <- chol(Om + diag(1e-9, nrow(Om)))
  W <- matrix(rnorm(nd * nrow(Om)), nd) %*% R
  if (is.null(Z)) return(mean(rowSums(W^2) >= S))
  A <- solve(crossprod(Z))
  mean(apply(W, 1, function(u) {
    zr <- crossprod(Z, u); as.numeric(t(zr) %*% A %*% zr)
  }) >= S)
}

#' The shrinkage-corrected Hosmer-Lemeshow test
#'
#' @param X   numeric matrix of covariates, WITHOUT an intercept column.
#' @param y   0/1 response.
#' @param lambda ridge penalty on the theory scale, lambda = n * lambda_glmnet.
#' @param G   number of groups (default 10).
#' @param basis "edge" (SC.EDGE, the default and the more powerful),
#'   "decile" (SC.HL), or both.
#' @param B   bootstrap replicates for the prepivoted reference.
#' @param seed optional integer; set it for a reproducible p-value.
#' @param penalize logical vector of length ncol(X): which columns are
#'   penalized. Defaults to all of them; the intercept is never penalized.
#' @param uncorrected if TRUE, also return the uncorrected p-value, which is
#'   the quantity the paper shows to be invalid. For comparison only.
#' @param finv.floor eigenvalue floor for the inversion of F in the debiasing
#'   step, as a fraction of the largest eigenvalue of F; 0, the default, is a
#'   plain solve(). On a badly conditioned design, repeat the test at 1e-9 and
#'   1e-7 and report the sensitivity: a floor changes the generator, so it is
#'   not a neutral numerical safeguard.
#'
#' @return a list with the statistic and p-value for each basis requested, and
#'   a $diagnostics element to report alongside them.
shrink.gof <- function(X, y, lambda, G = 10,
                       basis = c("edge", "decile"), B = 999,
                       seed = NULL, penalize = NULL, uncorrected = FALSE,
                       finv.floor = 0) {
  basis <- match.arg(basis, c("edge", "decile"), several.ok = TRUE)
  X <- as.matrix(X); y <- as.numeric(y)
  stopifnot(all(y %in% c(0, 1)), nrow(X) == length(y), lambda >= 0, G >= 3)
  if (!is.null(seed)) set.seed(seed)
  n  <- length(y)
  X1 <- cbind(1, X)
  if (is.null(penalize)) penalize <- rep(TRUE, ncol(X))
  D  <- diag(c(0, as.numeric(penalize)))       # intercept never penalized

  fit <- .sg_ridge(X1, y, lambda, D)
  st  <- .sg_pieces(X1, y, fit, lambda, D, G, finv.floor)
  if ("edge" %in% basis && is.null(st$Z))
    stop("the EDGE basis needs at least 4 distinct group means; raise G")

  # prepivot: resample from the DEBIASED generator pi(beta-tilde), not pi(beta-hat)
  gen <- pmin(pmax(as.numeric(1 / (1 + exp(-X1 %*% st$beta_tilde))), 1e-6), 1 - 1e-6)
  Sd <- Se <- numeric(B)
  for (b in seq_len(B)) {
    ys <- rbinom(n, 1, gen)
    s2 <- .sg_pieces(X1, ys, .sg_ridge(X1, ys, lambda, D), lambda, D, G,
                     finv.floor)
    Sd[b] <- s2$S_dec; Se[b] <- s2$S_edge
  }

  out <- list(lambda = lambda, G = G, B = B, n = n, p = ncol(X))
  if ("decile" %in% basis) {
    out$SC.HL <- list(statistic = st$S_dec,
                      p.value = (1 + sum(Sd >= st$S_dec)) / (B + 1))
    if (uncorrected)
      out$SC.HL$p.uncorrected <- .sg_mc_p(st$S_dec_unc, st$Om)
  }
  if ("edge" %in% basis) {
    out$SC.EDGE <- list(statistic = st$S_edge,
                        p.value = (1 + sum(Se >= st$S_edge, na.rm = TRUE)) / (B + 1))
    if (uncorrected)
      out$SC.EDGE$p.uncorrected <- .sg_mc_p(st$S_edge_unc, st$Om, st$Z)
  }
  # ---- diagnostics, to be reported with the p-value -------------------
  # ||beta.tilde|| is NOT the quantity to watch: it is huge whenever X has a
  # near-null space, and that on its own costs the test nothing. beta.tilde
  # enters only through the generator, so the generator is what to check. If
  # it pins observations at 0 or 1, or implies a number of events far from
  # the observed one, the bootstrap is not sampling the null world you meant.
  out$diagnostics <- list(
    cond.F     = kappa(st$Fm, exact = TRUE),
    gen.range  = range(gen),
    gen.sd     = sd(gen),
    pinned     = sum(gen <= 1e-6 | gen >= 1 - 1e-6),
    events.gen = sum(gen),
    events.obs = sum(y),
    norm.ratio = sqrt(sum(st$beta_tilde[-1]^2)) / sqrt(sum(fit$beta[-1]^2)),
    finv.floor = finv.floor)
  if (out$diagnostics$cond.F > 1e10 || out$diagnostics$pinned > 0)
    warning("F is ill conditioned (cond = ",
            signif(out$diagnostics$cond.F, 3), "; ",
            out$diagnostics$pinned,
            " generator probabilities pinned at 0 or 1). Repeat the test ",
            "with finv.floor = 1e-9 and 1e-7 and report the sensitivity.",
            call. = FALSE)

  class(out) <- "shrink.gof"
  out
}

print.shrink.gof <- function(x, ...) {
  cat("\nShrinkage-corrected Hosmer-Lemeshow test\n")
  cat(sprintf("n = %d, p = %d, lambda = %.4g, G = %d, B = %d\n\n",
              x$n, x$p, x$lambda, x$G, x$B))
  for (nm in c("SC.HL", "SC.EDGE")) if (!is.null(x[[nm]])) {
    cat(sprintf("  %-8s statistic = %8.4f   p = %.4f", nm,
                x[[nm]]$statistic, x[[nm]]$p.value))
    if (!is.null(x[[nm]]$p.uncorrected))
      cat(sprintf("   (uncorrected p = %.5f)", x[[nm]]$p.uncorrected))
    cat("\n")
  }
  if (!is.null(x$diagnostics)) {
    d <- x$diagnostics
    cat(sprintf("\n  cond(F) = %.3g, F inverse floored at %g\n",
                d$cond.F, d$finv.floor))
    cat(sprintf("  generator in [%.4f, %.4f], %d pinned, %.1f events against %d observed\n",
                d$gen.range[1], d$gen.range[2], d$pinned,
                d$events.gen, d$events.obs))
  }
  cat("\n")
  invisible(x)
}
